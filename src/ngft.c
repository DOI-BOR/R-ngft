/*
 *  ngft.c
 *  The New GFT Framework
 *
 *	Written by Chris Wood on 2/1/2017
 *	This software is a reimplementation of the GFT framework
 *	library by Robert Brown et al (2008). The interface for
 *	this new implementation is incompatible with the original GFT.
 *
 *  --- Original copyright notice from the GFT Framework is below
 *  >	Created by Robert Brown on 30/05/08.
 *	>	This software is copyright © 2010 UTI Limited Partnership.
 *	>	The original authors are Robert A. Brown, M. Louis Lauzon
 *	>	and Richard Frayne.  This software is licensed in the terms
 *	>	set forth in the “FST License Notice.txt” file, which is
 *	>	included in the LICENSE directory of this distribution.
 *
 */

#include "ngft.h"

#ifdef DllImport
#  undef DllImport
#  define DllImport DllExport
# endif
#include "ngft_proto.h"


// For odd width partitions, the center point is unambiguous. For even width, the center point
// can be at the end of the left half (left_bias=T), or the start of the right half (left_bias=F).
BOOL left_bias = TRUE;


// fftw wrapper from original GFT library
static void fft(int N, DCMPLX *in, int stride) {
	fftw_plan p;
	p = fftw_plan_many_dft(1, &N, 1, (fftw_complex *)in, NULL, stride, 0, (fftw_complex *)in,
													NULL, stride, 0, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}


// fftw wrapper from original GFT library
static void ifft(int N, DCMPLX *in, int stride) {
	int ii;
	fftw_plan p;
	p = fftw_plan_many_dft(1, &N, 1, (fftw_complex *)in, NULL, stride, 0, (fftw_complex *)in,
													NULL, stride, 0, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	for ( ii = 0 ; ii < N ; ii++ ) {
		in[ii*stride].r /= N;
		in[ii*stride].i /= N;
	}
}


// return the modulus of a DCMPLX value
static double modulus( DCMPLX *x ) {
	return sqrt( x->r*x->r + x->i*x->i );
}


// multiply 2 DCMPLX numbers, replacing first DCMPLX with the result
static void cmul(DCMPLX *x, DCMPLX *y) {
	double r = x->r * y->r - x->i * y->i;
	double i = x->r * y->i + x->i * y->r;
	x->r = r;
	x->i = i;
}


// divide 2 DCMPLX numbers, replacing first DCMPLX with the result
static void cdiv(DCMPLX *x, DCMPLX *y) {
	double r = x->r * y->r + x->i * y->i;
	double i = -x->r * y->i + x->i * y->r;
	double norm = y->r * y->r + y->i * y->i;
	x->r = r / norm;
	x->i = i / norm;
}


// multiply a DCMPLX number by a double, replacing the number with the result
static void cmulByReal(DCMPLX *x, double multiplier) {
	x->r *= multiplier;
	x->i *= multiplier;
}


// Circular shift an array of DCMPLX values
// Shift is to the right if amount > 0, and to the left if amount < 0
// Shift amount is modulo N
static void shift(DCMPLX *sig, int N, int amount) {
	int len;
	DCMPLX *tmp;
	int s1, s2, s3, s4;
	int tsiz = sizeof( *tmp );

	if ( sig == NULL || N == 0 || amount == 0 )
		return; // silently do nothing

	amount %= N;	// circular shifts are modulo N
	if ( amount < 0 )
		amount = N - ABS(amount);	// circular left shift is the same as a right shift of N - |amount|

	// approach is to shift by using a single memmove of the bigger piece, and two memcpy of the smaller piece.

	// shift right by |amount|
	if ( amount < N / 2 ) {
		len = amount; // length of smaller piece is |amount|
		// bigger piece begins at offset s3 = 0, smaller piece begins at offset s1 = N - |amount|
		// bigger piece destination is offset s2 = |amount|, smaller piece destination is offset s4 = 0
		s1 = N - len; s3 = 0;
		s2 = len; s4 = 0;
	} else {
		len = N - amount; // length of smaller (or equivalent) piece is N - |amount|
		// bigger piece begins at offset s3 = |amount|, smaller piece begins at offset s1 = 0
		// bigger piece destination is offset s2 = 0, smaller piece destination is offset s4 = N - |amount|
		s1 = 0; s3 = len;
		s2 = 0; s4 = N - len;
	}
	tmp = malloc( len * tsiz );
	memcpy( tmp, sig + s1, len * tsiz );
	memmove( sig + s2, sig + s3, (N - len) * tsiz );
	memcpy( sig + s4, tmp, len * tsiz );
	free( tmp );
}


// Return the FT of a time-domain Gaussian(t,sig-t) windowing function that is centered at t=0, has
// sigma-t = 1 / freq, and is normalized. The FT will be centered on frequency 0, and have a specified length N.
// Note: The FT of a Gaussian is also a Gaussian:
//	FT(Gaussian(t,sig-t)) = sqrt(2 * pi) * sigma-f * Gaussian(f,sig-f)
// where sigma-f = 1/(2*pi*sigma-t).  Therefore one could compute the frequency-domain window directly
// (paying careful attention to the phase shifts). Instead, this function will compute the DFT of the
// time-domain Gaussian. Using t -> i * dt, and freq -> j / (N * dt), we have
//	DFT(w(t,sig-t)) = [|freq| / sqrt(2*pi)] * DFT(exp(-(i/N)^2 * j^2 / 2)) = [.] * DFT(exp(x^2 * j^2 / 2))
// where x = [i * 2 * x_max / (N - 1) - x_max] is in the range [-x_max, +x_max], and x_max = 0.5.
static DCMPLX *gaussian_dft(int N, int freq) {
	int ii;
	// Should we avoid exp(-x) underflow for x > 708 on IEEE hardware? On underflow, exp() returns 0, and
	//		raises an FE_UNDERFLOW exception. To avoid underflow, we need x_max * freq_max < sqrt(2 * 708) = 38.
	//		The largest frequency is N/2. For dyadic scaling, freq_max ~ 3 * N / 8, which gives x_max < 100 / N.
	//		For freq_max = N/2, then x_max < 75 / N to avoid underflow. However, when x_max < 0.5, the frequency
	//		resolution is reduced by an amount gamma = x_max / 0.5. The best choice seems to be to allow underflow
	//		(resulting terms are 0), and avoid the decrease in resolution.
	double x, x_max = 0.5;  // MIN( 0.5, 100. / N );
	double dx = 2 * x_max / (N - 1);
	DCMPLX *win = calloc(N, sizeof( *win ));
	double win_sum;

	// get a Gaussian in time domain from -x_max to +x_max, centered at N/2
	double scale = ABS( (double)freq / N ) / sqrt( 2 * M_PI );
	for ( win_sum = 0 , ii = 0 ; ii < N ; ii++ ) {
		double g;
		x = (ii * dx - x_max) * freq;
		g = scale * exp(-0.5 * x * x);
		win[ii].r = g;
		win[ii].i = 0;
		win_sum += win[ii].r;
	}

	// Normalize time window so that the sum = 1, since the Gaussian may be significantly
	// truncated for small values of |freq|
	for ( ii = 0; ii < N ; ii++ )
		win[ii].r /= win_sum;

#ifdef NGFT_VERBOSE_DEBUG
	fprintf( stderr, "Gaussian Window: N=%2d freq=%2d, gamma=%.3f, x_max=%.3f, win_sum=%.3f\n", N, freq, x_max/0.5, x_max, win_sum );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	// left-shift the time window so that 0 time is at index 0, and negative times start at
	// index N/2+1. Note: left-shifting by N/2 handles even and odd N; right-shifting
	// by N/2 is correct only for even N
	shift(win, N, -N/2);

#ifdef NGFT_VERBOSE_DEBUG
	fprintf( stderr, "Shifted %d:\n", N/2 );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	// take the DFT of the 0-centered, time-domain Gaussian
	// Note: FT of a real, even function is also a real, even function. But, see notes below
	fft(N,win,1);

	// Note: If N is odd, then index 0 will be exactly at time 0, and the imaginary part
	// of the DFT output will be about 1e-14 less than the real, as expected.
	// If N is even, however, then time 0 is shifted by 1/2 of a sample period from index 0,
	// which introduces a phase shift, so the imaginary parts of the DFT aren't zero.

#ifdef NGFT_VERBOSE_DEBUG
	for ( win_sum = 0, ii = 0 ; ii < N ; ii++ )
		win_sum += modulus(win + ii);
	fprintf( stderr, "FFT: freq=%2d, win_sum=%.3f\n", freq, win_sum );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	return(win); // length is N
}


// return the FT of a time-domain sinc windowing function of
// width N. 
static DCMPLX *box_dft(int N, int freq) {
	int ii, lshift, abs_freq = ABS(freq);	// window width is abs_freq if epsilon=0; otherwise, it's wider
	DCMPLX *win = calloc(N, sizeof( *win ));

	// rather than re-creating the window, which will be wider than
	// abs_freq if epsilon > 0, just set all elements to 1 and
	// let the calling function select the elements it wants.
	for ( ii = 0 ; ii < N ; ii++ )
		win[ii].r = 1; // elements outside window are actually 0, but we can set everything to 1 here

	// this isn't needed if we're not setting the width of the boxcar here
	lshift = abs_freq / 2 + (IS_EVEN(abs_freq) && left_bias ? -1 : 0);
	shift(win, N, -lshift);

	if ( IS_EVEN(N) ) {
		// Note: If N is odd, then index 0 of the time-domain sinc function will be exactly
		// at time 0. If N is even, however, then time 0 is shifted by 1/2 of a sample period
		// from index 0, which introduces a phase shift, so the imaginary parts of the DFT
		// aren't zero. That phase shift is applied here.
		double fscale = M_PI / N;
		for ( ii = 1 ; ii <= N / 2 ; ii++ ) {
			DCMPLX ps = {cos(ii * fscale), (left_bias ? -1 : 1) * sin(ii * fscale)};
			if ( ii < N / 2 || left_bias )
				cmul(win + ii, &ps);
			if ( ii < N / 2 || ! left_bias ) {
				ps.i = -ps.i;	// need complex conjugate for corresponding negative frequency
				cmul(win + N - ii, &ps);
			}
		}
	}

#ifdef NGFT_VERBOSE_DEBUG
	{
		double win_sum = 0;
		for ( ii = 0 ; ii < N ; ii++ )
			win_sum += modulus(win + ii);
		fprintf(stderr, "FFT: freq=%2d, win_sum=%.3f\n", freq, win_sum);
		for ( ii = 0 ; ii < N ; ii++ )
			fprintf(stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus(win + ii));
		fflush(stderr);
	}
#endif

	return(win);
}


// conversion between indicies defining the start, center, and end of windows. The mapping
// between centers and endpoints is 1-1 for non-overlapping windows, and cases where the
// width is odd. For even width, and overlapping windows, there are ambiguities at the 0 and
// midpoint frequency.

static int shift_left(int idx, BOOL left_bias, BOOL is_f, int N) {
	// for frequency indicies, left_bias for positive frequencies implies right_bias for negative frequencies,
	// if the windows are to be symmetric between positive and negative frequencies.
	return is_f ? ((left_bias && idx <= N / 2) || (!left_bias && idx > N / 2)) : left_bias;
}
static int start_2_center(int start, int width, BOOL left_bias, BOOL is_f, int N) {
	int center = start + width / 2;

	// Note: there is ambiguity in which center to return for the cases where width is even,
	// and the true center is between the most positive and negative frequencies. However,
	// that case should not occur in this code. But, bomb if it does
	if ( is_f && IS_EVEN(width) && center == N / 2 + 1 )
		oops("start_2_center", "application exception: indeterminate center");

	if ( IS_EVEN(width) && shift_left(center, left_bias, is_f, N) )
		center -= 1;
	center %= N;
	if ( center < 0 )
		center += N;
	return center;
}
static int center_2_start(int center, int width, BOOL left_bias, BOOL is_f, int N) {
	int start = center - width / 2;
	// Note: for even width, there is ambiguity for the case of true center = +/- 1/2 (center = 0),
	// but taking the +1/2 case allows glossing over it for this code
	if ( IS_EVEN(width) && shift_left(center, left_bias, is_f, N) )
		start += 1;
	start %= N;
	if ( start < 0 )
		start += N;
	return start;
}
static int center_2_end(int center, int width, BOOL left_bias, BOOL is_f, int N) {
	int end = center + width / 2;
	// Note: for even width, there is ambiguity for the case of true center = +/- 1/2 (center = 0),
	// but taking the +1/2 case allow glossing over it for this code
	if ( IS_EVEN(width) && ! shift_left(center, left_bias, is_f, N) )
		end -= 1;
	end %= N;
	if ( end < 0 )
		end += N;
	return end;
}


static void make_tdset(int N, FPART *fpart, TDSET *tdset) {
	int ii, tc;
	double tc_exact;
	int fwin_width = fpart->win_len;
	double exact_tp_width = (double)N / fwin_width;

	if ( exact_tp_width < 1 )
		oops("make_tdset", "Invalid argument: bad N or bad partition window-length");

	// partitions for this time-decimation set are assumed to be null; make it so
	tdset->partitions = NULL;
	tdset->pcount = 0;

	// time-axis decimation (scaling) factor = N / frequency-window width
	tdset->decimation = exact_tp_width;

	// add time partitions for this decimation factor, rounding as necessary to fill N
	for ( tc = 0, tc_exact = 0, ii = 0 ; ii < fwin_width ; ii++, tc_exact += exact_tp_width ) {
		TPART* partition;
		tdset->partitions = realloc(tdset->partitions, ++(tdset->pcount) * sizeof(*(tdset->partitions)));
		partition = tdset->partitions + tdset->pcount - 1;
		partition->start = tc;
		partition->width = ROUND( tc_exact + exact_tp_width - tc );
		partition->end = partition->start + partition->width - 1;
		partition->center = start_2_center( partition->start, partition->width, left_bias, FALSE, N );

		tc += partition->width;
	}

	// sanity check
	if ( tdset->pcount != fwin_width )
		oops("make_tdset", "Application exception: tdset->pcount != fwin_width");
}


// create the time-partition collection corresponding to a specified frequency partition collection
static TPCOL *timePartitions(FPCOL *pars) {
	int ii, N;
	TPCOL *tpcol;

	if ( pars == NULL )
		oops("timePartitions", "Invalid argument: pars is NULL");
	N = pars->N;

	tpcol = calloc( 1, sizeof( *tpcol ) ); // note: everything in tpcol initialized to NULL or 0
	tpcol->N = N;

	// loop over the non-negative (index 0) and the negative (index 1)
	// frequency partition sets. Ordering doesn't matter
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the frequency partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			BOOL found = FALSE;
			if ( tpcol->tdcount > 0 ) {
				// search existing time-partition sets for a time-axis decimation factor
				// corresponding to this frequency-partition width
				int kk;
				for ( kk = 0 ; kk < tpcol->tdcount ; kk++ ) {
					if ( (int)ROUND(fpart->win_len * tpcol->tdsets[kk].decimation) == N ) {
						found = TRUE;
						break;
					}
				}
			}
			if ( ! found ) {
				// no time-partition sets with the correct decimation found, so create one and fill,
				// using decimation implied by this frequency-partition
				tpcol->tdsets = realloc(tpcol->tdsets, ++(tpcol->tdcount) * sizeof(*(tpcol->tdsets)));
				make_tdset(N, fpart, tpcol->tdsets + tpcol->tdcount - 1);
			}
		}
	}
	return(tpcol);
}


// free the space
static void freeTimePartitions(TPCOL *tpcol) {
	int ii;
	if ( tpcol == NULL )
		return;

	for ( ii = 0 ; ii < tpcol->tdcount ; ii++ ) {
		if ( tpcol->tdsets == NULL )
			continue;
		TDSET *tdset = tpcol->tdsets + ii;
		if ( tdset->pcount > 0 && tdset->partitions != NULL )
			free( tdset->partitions );
	}
	free( tpcol->tdsets );
	free( tpcol );
}


// set the limits of the window corresponding to a frequency partition. The window may exactly coincide
// with the partition (epsilon = 0), or it may be larger - up to N
static void setWindowLimits(FPART *partition, int N, double epsilon) {
	int win_start, win_end, win_len;

	if ( epsilon > 0 ) {
		BOOL is_f = TRUE;
		win_len = MIN(ROUND((1 + epsilon) * partition->width), N);	// window can't exceed length
		win_start = center_2_start( partition->center, win_len, left_bias, is_f, N );
		win_end = center_2_end( partition->center, win_len, left_bias, is_f, N );
	} else {
		// use partition start and end (no overlap)
		win_start = partition->start;
		win_end = partition->end;
		win_len = partition->width;
	}

	partition->win_start = win_start;
	partition->win_end = win_end;
	partition->win_len = win_len;
}


// for specified starting and ending frequencies fs and fe (both non-negative or negative)
// of a partitioning, define the start, center, and end indices into the N frequencies,
// and allocate a new partition structure to store them
static void add_freq_partition(FPART **ppart, int *pcount, int fs, int fe, int N, double epsilon) {
	FPART *p;
	int start, center, end, width;

	if ( fs * fe <= 0 && fs != fe )
		oops("add_partition", "argument exception: non-zero fs and fe must have the same sign, or both be zero");

	// add space for another partition
	*ppart = realloc(*ppart, ++(*pcount) * sizeof(**ppart));
	p = *ppart + *pcount - 1;

	if ( fs >= 0 ) {
		// handle non-negative frequencies
		// check that fs < fe, and silently fix if not
		if ( fs > fe ) {
			int tmp = fe;
			fe = fs;
			fs = tmp;
		}
		start = NONNEG_F_IND(fs, N);
		end = NONNEG_F_IND(fe, N);
	} else {
		// handle negative frequencies
		// check that fs > fe, and silently fix if not
		if ( fs < fe ) {
			int tmp = fe;
			fe = fs;
			fs = tmp;
		}
		start = NEG_F_IND(fe, N);	// reverse start and end indices for negative
		end = NEG_F_IND(fs, N);	// frequencies so that start < end for array ops
	}
	width = end - start + 1;
	// for odd width, center is unambiguous. For even width, center can be at end
	// of the left half (left_bias=T), or the start of the right half (left_bias=F).
	center = start_2_center( start, width, left_bias, TRUE, N );

	// set partition indices [0, N-1] into the N-length frequency array, using standard DFT ordering
	p->start = start;
	p->center = center;
	p->end = end;
	p->width = width;

	// set window indices [0, N-1] into the N-length frequency array, using standard DFT ordering
	setWindowLimits(p, N, epsilon);

	// initialize pointers to NULL
	p->tdset = NULL;
	p->gft = NULL;

	return;
}


// frequency-partition comparator for sorting on increasing start frequency.
static int part_cmp(const void* vp1, const void* vp2) {
	const FPART *p1 = vp1, *p2 = vp2;
	return p1->start - p2->start;
}


// Dyadic partitioning of a DFT frequency array of length N, where N can be odd or even,
//	and does not need to be a power of 2. The partitioning will span the array, with
//	no overlap in partitions.
static FPCOL *dyadicPartitions(int N, double epsilon) {
	int fs, pcount, n_pcount;
	FPART *partitions, *n_partitions;
	FPCOL *partition_collection;

	if ( N <= 0 )
		oops( "dyadicPartitions", "Invalid argument: N <= 0" );

	// initialize
	partitions = NULL, pcount = 0;	// non-negative frequency partitions
	n_partitions = NULL, n_pcount = 0; // negative frequency partitions

	// handle 0-frequency as a special case
	add_freq_partition( &partitions, &pcount, 0, 0, N, epsilon );

	// loop over powers of 2 in (positive) frequency, from 1 to Nyquist,
	// to make dyadic partitions for positive and negative frequencies
	for ( fs = 1 ; fs <= N / 2 ; fs *= 2 ) {
		// add positive frequencies
		add_freq_partition( &partitions, &pcount, fs, 2 * fs - 1, N, epsilon );
		// add negative frequencies (except for Nyquist, unless N is odd)
		if ( fs < N/2 || IS_ODD(N) )
			add_freq_partition( &n_partitions, &n_pcount, -fs, -2 * fs + 1, N, epsilon );
	}

	// sort the negative partitions by decreasing frequency to be
	// consistent with normal frequency ordering
	qsort(n_partitions, n_pcount, sizeof(*n_partitions), part_cmp);

	// put results into a partition collection, and return a pointer to it
	partition_collection = calloc(1, sizeof(*partition_collection));
	partition_collection->N = N;
	partition_collection->fpset[0].partitions = partitions;
	partition_collection->fpset[0].pcount = pcount;
	partition_collection->fpset[1].partitions = n_partitions;
	partition_collection->fpset[1].pcount = n_pcount;
	partition_collection->partition_type = FP_DYADIC;

	return partition_collection;
}


// Equal division of the octave (EDO) partitioning: Define a reference octave frequency,
// and subdivide the frequency range into octaves such that one of them is equal to the
// the reference octave. Next, subdivide the octaves into T equal log-intervals. Also known
// as equal temperament. Returns the discrete frequencies best approximating the computed
// continuous frequencies. Each of the N discrete frequencies will belong to exactly one partition.
// If epsilon > 0, then partition windows also will be defined, where the amount of window
// overlap is specified by epsilon.  This is defined by:
//		window_width / partition_width = 1 + epsilon
// The underlying partitions remain unchanged, and never overlap, regardless of epsilon.
// Note: f_ref and T can be set <= 0 to use default values
static FPCOL *edoPartitions(int N, double epsilon, int f_ref, int T) {
	int fs_prev, fe_prev, pcount, n_pcount;
	double f_ratio, ff_prev;
	FPART *partitions, *n_partitions;
	FPCOL *partition_collection;

	// complain about invalid arguments we can't fix
	if ( N < 1  )
		oops( "edoPartitions", "Invalid argument: must have N > 0" );

	if ( f_ref <= 0 ) // set reference frequency to default value, if not provided
		f_ref = N / 2;
	if ( f_ref > N / 2 )	// sanity check
		f_ref = N / 2;

	if ( T <= 0 ) // set T to default value, if not provided
		T = 5;
	if ( T > N / 4 )	// sanity check
		T = N / 4;

	f_ratio = pow(2, 1 / (double)T);	// frequency ratio between successive continuous frequencies

	// initialize
	partitions = NULL, pcount = 0;	// non-negative frequency partitions
	n_partitions = NULL, n_pcount = 0; // negative frequency partitions

	// get partitions with frequencies less than or equal to the reference
	ff_prev = f_ref;
	fs_prev = f_ref + 1;
	do {
		int fe = fs_prev - 1;
		double ff = ff_prev / f_ratio;
		int fs = ROUND(ff);
		ff_prev = ff;
		if ( fs > fe )
			continue;
		fs_prev = fs;
		// add positive frequencies
		add_freq_partition(&partitions, &pcount, fs, fe, N, epsilon);
		// add negative frequencies (except for Nyquist - unless N is odd - or 0)
		if ( fs > 0 && ( fs < N / 2  || IS_ODD(N) ) )
			add_freq_partition(&n_partitions, &n_pcount, -fs, -fe, N, epsilon);
	} while ( fs_prev > 0 );

	// get partitions with frequencies greater than the reference
	ff_prev = f_ref;
	fe_prev = f_ref;
	while ( fe_prev < N / 2 ) {
		int fs = fe_prev + 1;
		double ff = ff_prev * f_ratio;
		int fe = ROUND(ff);
		ff_prev = ff;
		if ( fe < fs )
			continue;
		fe_prev = fe;
		// add positive frequencies
		add_freq_partition(&partitions, &pcount, fs, fe, N, epsilon);
		// add negative frequencies (except for Nyquist, unless N is odd)
		if ( fs < N / 2 || IS_ODD(N) )
			add_freq_partition(&n_partitions, &n_pcount, -fs, -fe, N, epsilon);
	}

	// sort the partitions to be consistent with normal frequency ordering
	qsort(partitions, pcount, sizeof(*partitions), part_cmp);
	qsort(n_partitions, n_pcount, sizeof(*n_partitions), part_cmp);

	// put results into a partition collection, and return a pointer to it
	partition_collection = calloc(1, sizeof(*partition_collection));
	partition_collection->N = N;
	partition_collection->fpset[0].partitions = partitions;
	partition_collection->fpset[0].pcount = pcount;
	partition_collection->fpset[1].partitions = n_partitions;
	partition_collection->fpset[1].pcount = n_pcount;
	partition_collection->partition_type = FP_EDO;

	return partition_collection;
}


// Fixed width partitioning
// The underlying partitions remain unchanged, and never overlap, regardless of epsilon.
// Note: f_ref and T can be set <= 0 to use default values
static FPCOL *fwPartitions(int N, double epsilon, int W) {
	int fe_prev, pcount, n_pcount;
	FPART *partitions, *n_partitions;
	FPCOL *partition_collection;

	// complain about invalid arguments we can't fix
	if ( N < 1  )
		oops( "fwPartitions", "Invalid argument: must have N > 0" );

	if ( W <= 0 ) // set W to default value, if not provided
		W = 10;
	if ( W > N / 4 )	// sanity check
		W = N / 4;

	// initialize
	partitions = NULL, pcount = 0;	// non-negative frequency partitions
	n_partitions = NULL, n_pcount = 0; // negative frequency partitions

	// handle 0-frequency as a special case
	add_freq_partition( &partitions, &pcount, 0, 0, N, epsilon );

	// get fixed-width partitions starting from frequency = 1
	fe_prev = 0;
	do {
		int fs = fe_prev + 1;
		int fe = fs + W - 1;
		if ( fs > N / 2 )
			fe = N / 2;
		fe_prev = fe;
		// add positive frequencies
		add_freq_partition(&partitions, &pcount, fs, fe, N, epsilon);
		// add negative frequencies (except for Nyquist - unless N is odd)
		if ( fs < N / 2  || IS_ODD(N) )
			add_freq_partition(&n_partitions, &n_pcount, -fs, -fe, N, epsilon);
	} while ( fe_prev < N / 2 );

	// sort the partitions to be consistent with normal frequency ordering
	qsort(n_partitions, n_pcount, sizeof(*n_partitions), part_cmp);

	// put results into a partition collection, and return a pointer to it
	partition_collection = calloc(1, sizeof(*partition_collection));
	partition_collection->N = N;
	partition_collection->fpset[0].partitions = partitions;
	partition_collection->fpset[0].pcount = pcount;
	partition_collection->fpset[1].partitions = n_partitions;
	partition_collection->fpset[1].pcount = n_pcount;
	partition_collection->partition_type = FP_FW;

	return partition_collection;
}


// create frequency partitions
DllExport FPCOL *ngft_FrequencyPartitions(int N, double epsilon, FreqPartitionType ptype, FreqWindowType wtype,
																					int W, int f_ref, int T) {
	int ii;
	TPCOL *tpcol;
	FPCOL *pars;

	if ( N <= 0 )
		oops( "ngft_FrequencyPartitions", "Invalid argument: N <= 0" );

	// require fractional amount of overlap to be between 0 and N - 1
	epsilon = MIN( MAX( epsilon, 0 ), N - 1 );

	// make the selected type of frequency partition
	if ( ptype == FP_DYADIC )
		pars = dyadicPartitions(N, epsilon);
	else if ( ptype == FP_EDO )
		pars = edoPartitions(N, epsilon, f_ref, T);
	else if ( ptype == FP_FW )
		pars = fwPartitions(N, epsilon, W);
	else
		oops("ngft_FrequencyPartitions", "Invalid argument: unknown partition type");

	// initialize window type
	pars->window_type = wtype;

	// make a set of corresponding time partitions
	tpcol = timePartitions( pars );

	// assign a time partition to each corresponding frequency partition
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			BOOL found;
			FPART *partition = fpset->partitions + jj;

			// find the time-partition corresponding to this frequency partition
			for ( found = FALSE, kk = 0 ; kk < tpcol->tdcount ; kk++ ) {
				if ( (int)ROUND(tpcol->tdsets[kk].decimation * partition->win_len) == N ) {
					found = TRUE;
					partition->tdset = tpcol->tdsets + kk;
					break;
				}
			}
			if ( ! found ) {
				oops( "ngft_FrequencyPartitions", "Application exception: can't find corresponding time partition" );
			}
		}
	}

	// stash pointer to the collection of time partitions so we can find it later
	pars->tpcol = tpcol;

	return pars;
}


// free frequency-partition space
DllExport void ngft_FreeFreqPartitions(FPCOL *pars) {
	int ii;
	if ( pars == NULL )
		return;

	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET* fpset = pars->fpset + ii;
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			freeDClist(fpart->gft);
		}
		free( fpset->partitions );
	}
	if ( pars->tpcol != NULL )
		freeTimePartitions( pars->tpcol );
	free( pars );
}


// get shifted DFT of N-length time-domain window function. Shift is such that
// f = 0 of the original window is shifted to fcenter (which may be positive or negative)
static DCMPLX *getWindowDFT(int fcenter, int width, int N, windowFunction *window_fn) {
	DCMPLX *win;

	// get the N-length FT of the time window, centered about 0
	if ( width < 1 ) { // error: width must be >= 1
		oops( "getWindowDFT", "Application Error: width < 1" );
	} else if ( width == 1 ) { // length-1 windows are always the same
		win = calloc( N, sizeof( *win ) );
		win[0].r = 1;
	} else { // width > 1
		if ( fcenter == 0 ) {
			// sigma-t is infinite in this case, so the Gaussian is ill-defined. Need a boxcar from
			// start to end, with a height of 1 (TODO: should this be normalized to 1/width?)
			int ii;
			win = calloc( N, sizeof( *win ) );
			// rather than re-creating the window, or passing the window arguments, just set all N elements
			// to 1, and let the calling function grab 'width' of them from wherever it wants.  
			for ( ii = 0 ; ii < N ; ii++ )
				win[ii].r = 1.;
		} else {
			// get DFT of time-domain window of length N, centered on freq=0
			win = window_fn(N, fcenter);
		}
	}

	// right-shift DFT of window so that original 0-frequency index is moved to fcenter (positive or negative)
	shift(win, N, FREQ_2_INDEX(fcenter,N));

	return win; // length is N
}


// take fast discrete transform of a complex double 1D signal
DllExport FPCOL *ngft_1dComplex64(DCLIST *sig, double epsilon, FreqPartitionType ptype, FreqWindowType wtype,
																	int W, int f_ref, int T) {
	int ii, N, stride = 1;
	windowFunction *window_fn = NULL;
	FPCOL *fpcol;
	DCMPLX *signal;

	if ( sig == NULL || sig->count <= 0 )
		oops( "ngft_1dComplex64", "Invalid argument: signal is null or empty" );
	signal = sig->values;
	N = sig->count;

	fpcol = ngft_FrequencyPartitions( N, epsilon, ptype, wtype, W, f_ref, T );

	if ( fpcol->window_type == FWT_GAUSSIAN )
		window_fn = gaussian_dft;
	else if ( fpcol->window_type == FWT_BOX )
		window_fn = box_dft;
	else
		oops("ngft_1dComplex64", "Invalid argument: unknown frequency-window type");

	// get the spectrum of the signal (warning: overwrites input array). Result is length N, centered on freq=0
	fft(N, signal, stride);

	// loop over the non-negative (index 0) and the negative (index 1)
	// frequency partition sets. Ordering doesn't matter
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			DCMPLX *win, *dst;
			DCLIST* gft;
			FPART *partition = fpset->partitions + jj;
			int fcenter = INDEX_2_FREQ(partition->center, N);	// need actual frequency, not index
			int win_start = partition->win_start;
			int win_len = partition->win_len;

			// get DFT of N-length time-domain window function for this center frequency
			win = getWindowDFT(fcenter, win_len, N, window_fn); // length is N, original f=0 is centered on fcenter

			// make DST and apply the window to the spectrum
			dst = calloc( win_len, sizeof( *dst ) );
			for ( kk = 0 ; kk < win_len ; kk++ ) {
				int f_idx = (win_start + kk) % N; // get frequency index, handling wrap-around
				dst[kk] = signal[f_idx];	// copy and shift the signal FFT (0-index is at the window start)
				cmul( dst + kk, win + f_idx );	// apply the window
			}
			free( win );	// done with window

			// inverse FFT the windowed and transformed data to get this piece of S-space
			ifft( win_len, dst, stride );

			// save this row of the S-transform.
			gft = calloc(1, sizeof(*gft));
			gft->count = win_len;
			gft->values = dst;
			partition->gft = gft;

			// TODO: handle strides > 1
		}
	}

	return fpcol;	// this structure holds the partitions and the gft
}


// make a linear array of all the GFT elements
DllExport DCLIST *ngft_makeGftArray(FPCOL *fpcol) {
	int ii, dst_len = 0;
	DCMPLX *dst = NULL;
	DCLIST *gft;

	if ( fpcol == NULL )
		oops("ngft_makeGftArray", "Argument exception: fpcol is null");

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the frequency partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			int p_len = fpart->gft->count;
			dst_len += p_len;
			dst = realloc(dst, dst_len * sizeof(*dst));
			memcpy(dst + dst_len - p_len, fpart->gft->values, p_len * sizeof(*dst));
		}
	}
	gft = calloc(1, sizeof(*gft));
	gft->count = dst_len;
	gft->values = dst;
	return gft;
}


// unpack a linear array of all the GFT elements
DllExport void ngft_unpackGftArray(DCLIST *gft, FPCOL *fpcol) {
	int ii, dst_current = 0, dst_len;
	DCMPLX *dst;

	if ( gft == NULL || fpcol == NULL )
		oops("ngft_unpackGftArray", "Argument exception: gft or fpcol is null");
	dst_len = gft->count;
	dst = gft->values;

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the frequency partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			int p_len = fpart->tdset->pcount;
			if ( dst_current + p_len > dst_len )
				oops("ngft_unpackGftArray", "Application exception: gft length is too short for partitions");
			fpart->gft = calloc(1, sizeof(*(fpart->gft)));
			fpart->gft->values = calloc(p_len, sizeof(*dst));
			memcpy(fpart->gft->values, dst + dst_current, p_len * sizeof(*dst));
			fpart->gft->count = p_len;
			dst_current += p_len;
		}
	}

	return;
}


// take inverse fast discrete transform of a complex DST
DllExport DCLIST *ngft_1dComplex64Inv( FPCOL *fpcol ) {
	int ii, N, stride = 1;
	windowFunction *window_fn = NULL;
	DCMPLX *signal;
	DCLIST *sig;

	if ( fpcol == NULL )
		oops( "ngft_1dComplex64Inv", "Invalid argument: fpcol is null" );
	N = fpcol->N;

	if ( fpcol->window_type == FWT_GAUSSIAN )
		window_fn = gaussian_dft;
	else if ( fpcol->window_type == FWT_BOX )
		window_fn = box_dft;
	else
		oops("ngft_1dComplex64Inv", "Invalid argument: unknown frequency-window type");

	// create initialized space for signal
	signal = calloc( N, sizeof( *signal ) );

	// loop over the non-negative (index 0) and the negative (index 1)
	// frequency partition sets. Ordering doesn't matter
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			DCMPLX *dst_copy, *win;
			FPART *partition = fpset->partitions + jj;
			int fcenter = INDEX_2_FREQ(partition->center, N);	// need actual frequency, not index
			int win_start = partition->win_start;
			int win_len = partition->win_len;
			int start = partition->start;
			int end = partition->end;

			// make a copy of this partition's part of the S-transform, so as not to modify original
			dst_copy = calloc( win_len, sizeof( *dst_copy ) );
			memcpy(dst_copy, partition->gft->values, partition->gft->count * sizeof(*dst_copy));

			// FFT this partition's part of the S-transform
			fft(win_len, dst_copy, stride);

			// get DFT of N-length time-domain window function for this center frequency
			win = getWindowDFT(fcenter, win_len, N, window_fn); // length is N, original f=0 is centered on fcenter

			// remove partition window from the transformed dst
			for ( kk = 0 ; kk < win_len ; kk++ ) {
				int f_idx = (win_start + kk) % N; // get frequency index, handling wrap-around
				cdiv( dst_copy + kk, win + f_idx );
				// subset out those frequencies contained within the partition (drop overlapped points)
				if ( start <= f_idx && f_idx <= end )
					signal[f_idx] = dst_copy[kk];
			}
			free( dst_copy );
			free( win );
		}
	}

	// inverse FFT to recover the signal
	ifft(N, signal, stride);

	// encapsulate signal into structure and return
	sig = calloc(1, sizeof(*sig));
	sig->values = signal;
	sig->count = N;

	return sig;
}


#ifdef IS_PORTED
// 2D transform, mostly following the original GFT code.
// Need to write an inverse function to be useful for filtering applications
// Need an explicit mapping between S-transform indices and image partitions(wavenumber and pixel)
// This function compiles, but has not been tested by me (CKW)
DllExport void ngft_2dComplex64(DCMPLX *image, int N, int M, windowFunction *window_fn) {
	int row, col;
	double epsilon = -1;
	double f_ref = -1, T = -1, W = -1;
	FPCOL *pars;

	if ( image == NULL || N <= 0 || M <= 0 )
		oops( "ngft_2dComplex64", "Invalid argument: image is NULL or N <= 0 or M <= 0" );
	if ( window_fn == NULL )
		window_fn = gaussian;

	pars = ngft_FrequencyPartitions(N, epsilon);

	for ( row = 0 ; row < N ; row++ )
		ngft_1dComplex64(image + row * N, N, pars, window_fn, 1);

	if ( M != N ) {
		free( pars );
		pars = ngft_FrequencyPartitions(M, epsilon);
	}

	for ( col = 0 ; col < M ; col++ )
		ngft_1dComplex64(image + col, M, pars, window_fn, N);

	free( pars );
}
#endif


// free an ILIST structure
DllExport void freeIlist( ILIST *ilist ) {
	if ( ilist != NULL ) {
		if ( ilist->values != NULL && ilist->count > 0 )
			free( ilist->values );
		free( ilist );
	}
}


// free a DCLIST structure
DllExport void freeDClist( DCLIST *dclist ) {
	if ( dclist != NULL ) {
		if ( dclist->values != NULL && dclist->count > 0 )
			free( dclist->values );
		free( dclist );
	}
}


// get the count of frequency partitions
static int getFreqDim( FPCOL *fpcol, BOOL all_freqs ) {
	int ii, imax, f_dim;
	// use non-negative (index 0), and optionally, negative (index 1) frequency partition sets
	imax = all_freqs ? 2 : 1;
	for ( f_dim = 0, ii = 0 ; ii < imax ; ii++ ) {
		FPSET *fpset = fpcol->fpset + ii;
		f_dim += fpset->pcount;
	}
	return f_dim;
}


// get the center-values of the frequency partitions, in order of the partitions
DllExport ILIST *freqPartitionCenters(FPCOL *fpcol, BOOL all_freqs) {
	int ii, imax, f_count, f_dim, N;
	ILIST *f_centers;
	int *values;

	N = fpcol->N;
	f_dim = getFreqDim(fpcol, all_freqs);
	values = calloc(f_dim, sizeof(*values));

	// use non-negative (index 0), and optionally, negative (index 1) frequency partition sets
	imax = all_freqs ? 2 : 1;
	for ( f_count = 0, ii = 0 ; ii < imax ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *part = fpset->partitions + jj;
			values[f_count++] = INDEX_2_FREQ(part->center, N); // need actual frequency, not array index
		}
	}

	f_centers = calloc( 1, sizeof( *f_centers ) );
	f_centers->values = values;
	f_centers->count = f_count;

	return f_centers;
}


// get the center-values of the time-partition with the most partitions, in order of the partitions
DllExport ILIST *timePartitionCenters( FPCOL *fpcol ) {
	int ii, t_dim;
	TDSET *tdset_max;
	ILIST *t_centers;
	int *values;

	TPCOL *tpcol = fpcol->tpcol;

	// find the partition set with max number of partitions
	for ( t_dim = 0, tdset_max = NULL, ii = 0 ; ii < tpcol->tdcount ; ii++ ) {
		TDSET *tdset = tpcol->tdsets + ii;
		if ( tdset->pcount > t_dim ) {
			t_dim = tdset->pcount;
			tdset_max = tdset;
		}
	}

	// grab the centers from these partitions
	values = calloc( t_dim, sizeof( *values ) );
	for ( ii = 0 ; ii < tdset_max->pcount ; ii++ ) {
		TPART *part = tdset_max->partitions + ii;
		values[ii] = part->center;
	}
	t_centers = calloc( 1, sizeof( *t_centers ) );
	t_centers->values = values;
	t_centers->count = t_dim;

	return t_centers;
}


// get the maximum partition count over the sets of partitions with common decimation factors
static int getTimeDim(TPCOL *tpcol) {
	int ii, t_dim;
	for ( t_dim = 0, ii = 0 ; ii < tpcol->tdcount ; ii++ ) {
		TDSET *tdset = tpcol->tdsets + ii;
		if ( tdset->pcount > t_dim )
			t_dim = tdset->pcount;
	}
	return t_dim;
}


// find DST element corresponding to frequency ff and time tt.
static DCMPLX *find_element(FPCOL *fpcol, int ff, int tt, BOOL all_freqs) {
	int ii, imax;
	int N = fpcol->N;
	int f_target = FREQ_2_INDEX(ff, N);
	int t_target = tt;

	// loop over the non-negative (index 0) and (optionally) the negative (index 1) frequency partition sets
	// Note: This function requires that the frequency partition sets are in normal frequency
	// order, i,e, positive sets: (0, 1, ... , max_f, [Nyquist]) negative sets: (-max_f, ... -1)
	imax = all_freqs ? 2 : 1;
	for ( ii = 0 ; ii < imax ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the frequency partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			if ( fpart->start <= f_target && f_target <= fpart->end ) {
				// frequency is in this partition; now check time partitions
				int kk;
				TDSET *tdset = fpart->tdset;
				if ( tdset == NULL )
					oops( "find_index", "Application exception: tdset for this partition is NULL" );
				for ( kk = 0 ; kk < tdset->pcount ; kk++ ) {
					// search for the partition with matching time
					TPART *tpart = tdset->partitions + kk;
					if ( tpart->start <= t_target && t_target <= tpart->end )
						return fpart->gft->values + kk;	// cool
				}
			}
		}
	}
	return NULL;	// bummer
}


// Make a 2D time-frequency image using nearest-neighbor sampling of the fast S transform.
// Sampling is based either on: (1) the centers of the time and frequency partitions, by setting
// by_part = TRUE, or, (2) the complete (redundant) set of times and frequencies implied by the
// original signal sampling and FT (by_part = FALSE). If the original signal was N points long,
// then using by_part = TRUE results in an image dimension of roughly log2N x log2N; otherwise the
// image will be NxN. Negative frequencies will be included if all_freqs is TRUE. A down-sampled
// image can be obtained by specifying M < N (log2N, if by_part = TRUE), resulting in an image of
// approximately MxM. If M <= 0, then down-sampling will not be used. If make_ind_map = TRUE, then
// the image will consist of the Fast S-Transform nearest-neighbor indicies at the specified sampling.
// The returned image is stored by row, with time on the horizontal axis, and frequency on the vertical.
// ILISTs of the center times and frequencies are also provided, in units of the original NxN time-freq space
DllExport DIMAGE *ngft_1d_InterpolateNN(FPCOL *fpcol, int M, BOOL by_part, BOOL all_freqs) {
	int ii, f_dim, nf_dim, t_dim, nt_dim, img_len, N;
	ILIST *f_centers = NULL, *t_centers = NULL;
	DCMPLX *image_arr;
	DIMAGE *image;
	TPCOL *tpcol;
	double t_factor, f_factor;
	BOOL down_sample;

	if ( fpcol == NULL )
		oops( "ngft_1d_InterpolateNN", "Invalid argument: fpcol is NULL" );
	N = fpcol->N;

	if ( fpcol->tpcol == NULL )
		oops( "ngft_1d_InterpolateNN", "Application exception: tpcol not set in fpcol" );
	tpcol = fpcol->tpcol;

	t_dim = by_part ? getTimeDim(tpcol) : N;
	if ( M <= 0 )
		M = t_dim;
	M = MIN( M, t_dim );	// don't allow augmentation
	down_sample = (M < t_dim);
	t_factor = (double)M / t_dim;
	nt_dim = M;

	f_dim = by_part ? getFreqDim(fpcol, all_freqs) : all_freqs ? N : N / 2;
	f_factor = down_sample ? MIN( (all_freqs ? 1 : 2) * t_factor, 1 ) : 1;
	nf_dim = ROUND(f_factor * f_dim);

	img_len = nf_dim * nt_dim;
	image_arr = calloc(img_len, sizeof(*image_arr));

	// get centers of frequency and time partitions
	if ( by_part ) {
		// TODO: this implementation is a kludge if resampling - clean it up by having separate lists
		f_centers = freqPartitionCenters( fpcol, all_freqs );
		t_centers = timePartitionCenters( fpcol );
	} else {
		f_centers = calloc( 1, sizeof( *f_centers ) );
		f_centers->values = calloc( f_dim, sizeof( *(f_centers->values) ) );
		f_centers->count = f_dim;
		t_centers = calloc( 1, sizeof( *t_centers ) );
		t_centers->values = calloc( t_dim, sizeof( *(t_centers->values) ) );
		t_centers->count = t_dim;
	}

	// get image. TODO: if downsampling, need to take averages, not just use a single sample
	for ( ii = 0 ; ii < nf_dim ; ii++ ) {
		int jj;
		int is = down_sample ? ROUND(ii / f_factor) : ii;
		int ff = by_part ? f_centers->values[is] : INDEX_2_FREQ( is, N );
		if ( ! by_part || f_dim != nf_dim )
			f_centers->values[ii] = ff;	// set f-center (or possibly overwrite, if by_part). TODO: clean up this kludge
		for ( jj = 0 ; jj < nt_dim ; jj++ ) {
			DCMPLX *gft;
			int img_ind;
			int js = down_sample ? ROUND(jj / t_factor) : jj;
			int tt = by_part ? t_centers->values[js] : js;
			if ( (gft = find_element(fpcol, ff, tt, all_freqs)) == NULL )
				oops( "ngft_1d_interpolateNN", "Application exception: can't find gft element" );
			if ( ii == nf_dim - 1 && ( ! by_part || t_dim != nt_dim ) )
				t_centers->values[jj] = tt;	// set t-center (or possibly overwrite, if by_part). TODO: clean up this kludge
			img_ind = ii * nt_dim + jj; // image stored by row, with time on the horizontal axis, and frequency on the vertical
			memcpy( image_arr + img_ind, gft, sizeof( *image_arr ) );
		}
	}

	if ( by_part ) {
		// re-size f and t centers if downsampling. TODO: clean up this kludge
		if ( t_dim != nt_dim ) {
			t_centers->values = realloc(t_centers->values, nt_dim * sizeof(*(t_centers->values)));
			t_centers->count = nt_dim;
		}
		if ( f_dim != nf_dim ) {
			f_centers->values = realloc(f_centers->values, nf_dim * sizeof(*(f_centers->values)));
			f_centers->count = nf_dim;
		}
	}

	image = calloc( 1, sizeof( *image ) );
	image->img = calloc(1, sizeof(*(image->img)));
	image->img->values = image_arr;
	image->img->count = img_len;
	image->x_centers = t_centers;
	image->y_centers = f_centers;

	return image;
}


// free space for an allocated DIMAGE pointer
DllExport void freeDImage( DIMAGE *image ) {
	if ( image != NULL ) {
		if ( image->img != NULL ) {
			if ( image->img->values != NULL && image->img->count > 0 )
				free( image->img->values );
			free( image->img );
		}
		if ( image->x_centers != NULL )
			freeIlist( image->x_centers );
		if ( image->y_centers != NULL )
			freeIlist( image->y_centers );
		free( image );
	}
}
