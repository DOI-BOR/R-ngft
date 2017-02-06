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


// For odd width partitions, center is unambiguous. For even width, center can be at
// the end of the left half (left_bias=T), or the start of the right half (left_bias=F).
BOOL left_bias = FALSE;	// FALSE seems to give better contrast for resolving impulses, and better inverse


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


// Create a Gaussian, with sigma-f = freq, in the frequency domain by taking
// the DFT of a Gaussian in the time domain with sigma-t = 1/(2*pi*freq).
// For this case, g(f) = FT[exp(-4 * pi^2 * t^2 * freq^2 / 2)]. See, e.g., Brigham (1974).
// For the discrete case, t -> i * dt, freq -> j / (N * dt), so we have
// Gaussian(f) = DFT(exp(-4 * pi^2 * (i/N)^2 * j^2 / 2)] = DFT[exp(x^2 * j^2 / 2)]
// where x = [i * 2 * x_max / (N - 1) - x_max] is in the range [-x_max, +x_max], and x_max = pi.
DllExport DCMPLX *gaussian(int N, int freq) {
	int ii;
	// Avoid exp(-x) underflow for x > 708 on IEEE hardware? On underflow, exp() returns 0,
	//		and raises an FE_UNDERFLOW exception. Need x_max * freq_max < sqrt(2 * 708) = 38,
	//		to avoid underflow. The largest frequency is N/2. For dyadic scaling, freq_max ~ 3 * N / 8),
	//		which gives x_max < 100 / N. For freq_max = N/2, then x_max < 75 / N. Original GFT code used
	//		x_max = 0.5. When x_max < pi, frequency resolution is reduced by an amount gamma = x_max / pi.
	//		Best choice seems to be to allow underflow (terms are 0), or else resolution is decreased.
	double x, x_max = M_PI;  // MIN( M_PI, 100. / N );
	double dx = 2 * x_max / (N - 1);
	DCMPLX *win = calloc(N, sizeof( *win ));
	double win_sum;
#ifdef NGFT_VERBOSE_DEBUG
	BOOL N_is_odd = (N % 2 == 0 ? FALSE : TRUE);
#endif

	// TODO: seems more obvious to construct windows directly in the frequency domain, which
	// is easy to do analytically since the FT of a Gaussian is also a Gaussian. Would
	// need to include phase shift (for even N) corresponding to time shift of 1/2 delta-t at t=0.
	// Other difference between FT and DFT will be sinc(x) convolution due to windowing over N.

	// get Gaussian in time domain from -x_max to +x_max, centered at N/2
	double scale = ABS( freq ) / sqrt( 2 * M_PI );
	for ( win_sum = 0 , ii = 0 ; ii < N ; ii++ ) {
		double g;
		x = (ii * dx - x_max) * freq;
		g = scale * exp(-0.5 * x * x);
		win[ii].r = g;
		win[ii].i = 0;
		win_sum += win[ii].r;
	}

	// Normalize window.
	for ( ii = 0; ii < N ; ii++ )
		win[ii].r /= win_sum;

#ifdef NGFT_VERBOSE_DEBUG
	fprintf( stderr, "Gaussian Window: N=%2d freq=%2d, gamma=%.3f, x_max=%.3f, win_sum=%.3f\n", N, freq, x_max/M_PI, x_max, win_sum );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	// left-shift time window so that 0 time is at index 0, and negative times start at
	// index N/2+1. Note: left-shifting by N/2 handles even and odd N; right-shifting
	// by N/2 is correct only for even N
	shift(win, N, -N/2);

#ifdef NGFT_VERBOSE_DEBUG
	fprintf( stderr, "Shifted %d:\n", N/2 );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	// take the DFT
	// Note: FT of a real, even function is also a real, even function. But, see below
	fft(N,win,1);

#ifdef NGFT_VERBOSE_DEBUG
	fprintf( stderr, "FFT: freq=%2d\n", freq );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	// If N is odd, then index 0 will be exactly at time 0, and the imaginary part of
	// the DFT output indeed seems to be about 1e-14 less than the real, as expected.
	// If N is even, however, then time 0 is shifted by 1/2 of a sample period from index 0,
	// which introduces a phase shift, so the imaginary parts won't be zero.

#ifdef NGFT_TEST_FORCE_REAL
	// For even N, should I therefore replace the real part with the modulus, and set imaginary
	// part to 0? Window won't be exactly symmetric, but it will be real.
	// Tests indicate that NOT doing this performs better for the fwd -> inv test. Also, going to
	// an odd value of N didn't improve the fwd -> inv test results
	if ( ! N_is_odd ) {
		for ( ii = 0 ; ii < N ; ii++ ) {
			win[ii].r = modulus( win + ii );
			win[ii].i = 0;
		}
	}
#endif

	// right-shift frequency window so that 0-index is centered at freq (positive or negative)
	shift(win, N, FREQ_2_INDEX(freq,N));

#ifdef NGFT_VERBOSE_DEBUG
	fprintf( stderr, "Final FFT after shift and setting to real: freq=%2d\n", freq );
	for ( ii = 0 ; ii < N ; ii++ )
		fprintf( stderr, "\t%4d:  %10.3e %10.3e\t%9.3e\n", ii, win[ii].r, win[ii].i, modulus( win + ii ) );
	fflush( stderr );
#endif

	return(win);
}


DllExport DCMPLX *box(int N, int freq) {
	int ii;
	DCMPLX *win = calloc(N, sizeof( *win ));
	for ( ii = 0 ; ii < N ; ii++ ) {
		win[ii].r = 1.;
		win[ii].i = 0.0;
	}

	return(win);
}


static void add_dyadic_partition(FPART **ppart, int *pcount, int fs, int N) {
	FPART *p, *partitions = *ppart;
	int start, center, end, width;

	partitions = realloc(partitions, ++(*pcount) * sizeof(*partitions));
	p = partitions + *pcount - 1;

	if ( fs >= 0 ) {
		// handle non-negative frequencies
		start = NONNEG_F_IND(fs, N);
		end = NONNEG_F_IND(2 * fs - 1, N);
		width = end - start + 1;
		// for odd width, center is unambiguous. For even width, center can be at end
		// of the left half (left_bias=T), or the start of the right half (left_bias=F).
		center = left_bias ? end - width / 2 : start + width / 2;
	} else {
		// handle negative frequencies
		start = NEG_F_IND(2 * fs + 1, N);	// reverse start and end indices for negative
		end = NEG_F_IND(fs, N);	// frequencies so that start < end for array ops
		width = end - start + 1;
		center = left_bias ? start + width / 2 : end - width / 2;
	}
	p->start = start;
	p->center = center;
	p->end = end;
	p->width = width;
	p->window = NULL;
	p->win_len = 0;

	*ppart = partitions;

	return;
}


// Dyadic partitioning of a DFT frequency array of length N, where N can be odd or even,
//	and does not need to be a power of 2. The partitioning will span the array, with
//	no overlap in partitions.
DllExport FPCOL *ngft_DyadicPartitions(int N) {
	int fs, pcount, n_pcount;
	FPART *partitions, *n_partitions, *tmp;
	FPCOL *partition_set;
	BOOL N_is_odd = (N % 2 == 0 ? FALSE : TRUE);

	if ( N <= 0 )
		oops( "ngft_DyadicPartitions", "Invalid argument: N <= 0" );

	partitions = NULL, pcount = 0;	// non-negative frequency partitions
	n_partitions = NULL, n_pcount = 0; // negative frequency partitions
	// handle 0-frequency as a special case
	fs = 0, add_dyadic_partition( &partitions, &pcount, fs, N );
	// loop over powers of 2 in (positive) frequency, from 1 to Nyquist,
	// to make dyadic partitions
	for ( fs = 1 ; fs <= N / 2 ; fs *= 2 ) {
		// add positive frequencies
		add_dyadic_partition( &partitions, &pcount, fs, N );
		// add negative frequencies (except for Nyquist, unless N is odd)
		if ( fs < N/2 || N_is_odd )
			add_dyadic_partition( &n_partitions, &n_pcount, -fs, N );
	}

	// reverse the order of the negative partitions to be consistent with normal frequency ordering
	tmp = calloc( n_pcount, sizeof( *tmp ) );
	for ( fs = n_pcount ; fs > 0 ; fs-- )
		memcpy(tmp + n_pcount - fs, n_partitions + fs - 1, sizeof( *n_partitions ));
	memcpy( n_partitions, tmp, n_pcount * sizeof( *tmp ) );
	free( tmp );

	partition_set = calloc(1, sizeof(*partition_set));
	partition_set->N = N;
	partition_set->fpset[0].partitions = partitions;
	partition_set->fpset[0].pcount = pcount;
	partition_set->fpset[1].partitions = n_partitions;
	partition_set->fpset[1].pcount = n_pcount;

	return partition_set;
}


// Partitioning along the lines of a 12-tone scale (12 log divisions per doubling),
// mostly following the original GFT code.
//	cents * 100 is the desired logFrequency increment, in units of semitones
//	samplerate is the sampling rate, in Hz. For audio, this should be above 2000 Hz, and
//		if it is, then the lower-frequency reference will be about 2 octaves below 440 Hz.
// This function compiles, but has not been tested by me (CWood)
DllExport FPCOL *ngft_1dMusicPartitions(int N, double samplerate, int cents) {
	int ii, ii_max, pcount, n_pcount;
	FPART *p, *partitions, *np, *n_partitions, *tmp;
	FPCOL *partition_set;

	double fSpacing;
	double reference = 440. / (2 * 2); // make reference 2 octaves below 440 Hz
	double logcent = 1. / (12 * 100);	// assume a 12-tone scale
	double logdelta, minlog2f, maxlog2f;

	if ( N <= 0 )
		oops( "ngft_1dMusicPartitions", "Invalid argument: N <= 0" );

	logdelta = logcent * cents;
	if ( samplerate < 2000 ) {
		// adjust reference for rates that appear to be well below the human audio band
		reference = samplerate / 100;
	}
	minlog2f = log2(reference) - logdelta * floor(log2(reference) / logdelta);
	maxlog2f = log2(0.5 * samplerate);	// Nyquist
	fSpacing = samplerate / N;
	while ( pow(2, minlog2f + logdelta) - pow(2, minlog2f) < fSpacing )
		minlog2f += logdelta;

	ii_max = (int)(floor((maxlog2f - minlog2f) / logdelta)) + 1;
	for ( partitions = NULL, n_partitions = NULL, pcount = 0, n_pcount = 0, ii = 0 ; ii < ii_max ; ii++ ) {
		double log2f;
		int fs, fe;
		// handle non-negative frequencies
		partitions = realloc(partitions, ++pcount * sizeof(*partitions));
		p = partitions + pcount - 1;
		log2f = minlog2f - 0.5 * logdelta + ii * logdelta;
		fs = ROUND( pow( 2, log2f ) / fSpacing );
		log2f = minlog2f - 0.5 * logdelta + (ii + 1) * logdelta;
		fe = ROUND( pow( 2, log2f ) / fSpacing ) - 1;
		p->start = NONNEG_F_IND(fs, N);
		p->end = NONNEG_F_IND(fe, N);
		p->width = p->end - p->start + 1;
		p->center = NONNEG_F_IND(p->start + p->width / 2, N);
		p->window = NULL;
		p->win_len = 0;
		if ( fs > 0 ) {
			// handle negative frequencies
			n_partitions = realloc(n_partitions, ++n_pcount * sizeof(*n_partitions));
			np = n_partitions + n_pcount - 1;
			np->start = NEG_F_IND(-p->start, N);
			np->end = NEG_F_IND(-p->end, N);
			np->width = np->start - np->end + 1;	// np->width <= p->width
			np->center = NEG_F_IND(p->start + np->width / 2, N);
			np->window = NULL;
			np->win_len = 0;
		}
	}

	// reverse the order of the negative partitions to be consistent with normal frequency ordering
	tmp = calloc( n_pcount, sizeof( *tmp ) );
	for ( ii = n_pcount ; ii > 0 ; ii-- )
		memcpy(tmp + n_pcount - ii, n_partitions + ii - 1, sizeof( *n_partitions ));
	memcpy( n_partitions, tmp, n_pcount * sizeof( *tmp ) );
	free( tmp );

	partition_set = calloc(1, sizeof(*partition_set));
	partition_set->N = N;
	partition_set->fpset[0].partitions = partitions;
	partition_set->fpset[0].pcount = pcount;
	partition_set->fpset[1].partitions = n_partitions;
	partition_set->fpset[1].pcount = n_pcount;

	return partition_set;
}


// free the space
DllExport void ngft_FreeFreqPartitions(FPCOL *pars) {
	int ii;
	if ( pars == NULL )
		return;

	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET* fpset = pars->fpset + ii;
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART* partition = fpset->partitions + jj;
			if ( partition == NULL )
				continue;
			if ( partition->window != NULL && partition->win_len > 0 )
				free( partition->window );
		}
		free( fpset->partitions );
	}
	free( pars );
}


static void fill_tdset(int N, FPART *fpart, TDSET *tdset) {
	int ii, tc;
	double tc_exact, exact_width;
	int fpwidth = fpart->width;
	exact_width = (double)N / fpwidth;

	// partitions for this time-decimation set are assumed to be null; make it so
	tdset->partitions = NULL;
	tdset->pcount = 0;

	// time-axis decimation factor = N / frequency-partition width
	tdset->decimation = (int)floor(exact_width);

	// add time partitions for this decimation factor, rounding as necessary to fill N
	for ( tc = 0, tc_exact = 0, ii = 0 ; ii < fpwidth ; ii++, tc_exact += exact_width ) {
		TPART* partition;
		tdset->partitions = realloc(tdset->partitions, ++(tdset->pcount) * sizeof(*(tdset->partitions)));
		partition = tdset->partitions + tdset->pcount - 1;
		partition->start = tc;
		partition->width = ROUND( tc_exact + exact_width - tc );
		partition->end = partition->start + partition->width - 1;
		partition->center = left_bias ? partition->end - partition->width / 2 :
				partition->start + partition->width / 2;

		tc += partition->width;
	}
}


// create the time-partition collection corresponding to a specified frequency partition collection
DllExport TPCOL *ngft_TimePartitions(FPCOL *pars) {
	int ii, N;
	TPCOL *tpcol;

	if ( pars == NULL )
		oops("ngft_TimePartitions", "Invalid argument: pars is NULL");
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
					if ( tpcol->tdsets[kk].decimation == N / fpart->width ) {
						found = TRUE;
						break;
					}
				}
			}
			if ( ! found ) {
				// no time-partition sets with the correct decimation found, so create one and fill,
				// using decimation implied by this frequency-partition
				tpcol->tdsets = realloc(tpcol->tdsets, ++(tpcol->tdcount) * sizeof(*(tpcol->tdsets)));
				fill_tdset(N, fpart, tpcol->tdsets + tpcol->tdcount - 1);
			}
		}
	}
	return(tpcol);
}


// free the space
DllExport void ngft_FreeTimePartitions(TPCOL *tpcol) {
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


// create partitions and corresponding windows from a specified N and window function
DllExport FPCOL *ngft_MakePartsAndWindows(int N, windowFunction *window_fn) {
	if ( N <= 0 )
		oops( "ngft_MakePartsAndWindows", "Invalid argument: N <= 0" );
	if ( window_fn == NULL )
		window_fn = gaussian;

	FPCOL *partition_set = ngft_DyadicPartitions( N );
	ngft_AddWindowsToParts(partition_set, window_fn);
	return partition_set;
}


// Create windows and associate with specified partitions
// This should be the only place in the ngft library where windows are created
DllExport void ngft_AddWindowsToParts(FPCOL *pars, windowFunction *window_fn) {
	int ii;
	int N;

	if ( pars == NULL )
		oops("ngft_AddWindowsToParts", "Invalid argument: pars is NULL");
	if ( window_fn == NULL )
		window_fn = gaussian;

	N = pars->N;

	// loop over the non-negative (index 0) and the negative (index 1)
	// frequency partition sets. Ordering doesn't matter
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			DCMPLX *win;
			FPART *partition = fpset->partitions + jj;
			int win_len = partition->width;
			int fcenter = INDEX_2_FREQ(partition->center,N);	// need actual frequency, not array index
			// create a window of length N in the frequency domain, centered on fcenter
			if ( win_len > 1 ) {
				int kk;
				DCMPLX *full_win = window_fn( N, fcenter );

				// Subset out the part of the full window overlapping this partition. The partition
				// values are frequencies, expressed as non-negative indices into the N-length
				// spectrum of the full window. Full window is centered on fcenter, which is positioned
				// at array index N/2 - 1 + is_odd
				win = calloc( win_len, sizeof( *win ) );
				memcpy( win, full_win + partition->start, win_len * sizeof( *win ) );
				free( full_win );
			} else {
				if ( win_len < 1 )
					oops( "ngft_AddWindowsToParts", "Application Error: win_len < 1" );	// win_len must be >= 1
				win = calloc( win_len, sizeof( *win ) );	// win_len is 1
				win[0].r = 1; win[0].i = 0;
			}

			partition->window = win;
			partition->win_len = win_len;
		}
	}
}


// take fast discrete transform of a complex double 1D signal
DllExport void ngft_1dComplex64(DCMPLX *signal, int N, FPCOL *pars, int stride) {
	int ii;
	BOOL free_pars = FALSE;

	if ( signal == NULL )
		oops( "ngft_1dComplex64", "Invalid argument: Signal is null" );
	if ( N <= 0 && pars == NULL )
		oops( "ngft_1dComplex64", "Invalid argument: Must specify either non-NULL pars, or N > 0" );

	if ( pars == NULL ) {
		pars = ngft_MakePartsAndWindows( N, gaussian );	// N is used only if pars is NULL
		free_pars = TRUE;
	} else if ( N > 0 && pars->N != N )
		smsg( "ngft_1dComplex64", "Warning: Invalid argument: N in pars conflicts with argument N (ignored)" );
	N = pars->N;

	stride = MAX(1, stride);	// if stride <= 0, set to default value of 1

	// Do the initial FFT of the signal
	fft(N, signal, stride);

	// Apply the windows to the FFT, and take the
	//	inverse FFT to get the S transform

	// loop over the non-negative (index 0) and the negative (index 1)
	// frequency partition sets. Ordering doesn't matter
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			FPART *partition = fpset->partitions + jj;
			int fstart = partition->start;
			// apply partition window to the transformed data
			for ( kk = 0 ; kk < partition->win_len ; kk++ ) {
				int ff = fstart + kk;
				cmul(signal + ff * stride, partition->window + kk);
			}
			// inverse FFT the windowed and transformed data to get this piece of S-space
			ifft(partition->width, signal + fstart * stride, stride);
		}
	}

	if ( free_pars )
		free( pars );
}


// take inverse fast discrete transform of a complex DST
DllExport void ngft_1dComplex64Inv( DCMPLX *dst, FPCOL *pars, int stride ) {
	int ii, N;

	if ( dst == NULL )
		oops( "ngft_1dComplex64Inv", "Invalid argument: dst is NULL" );
	if ( pars == NULL )
		oops( "ngft_1dComplex64Inv", "Invalid argument: Pars is NULL" );

	stride = MAX( 1, stride );	// if stride <= 0, set to default value of 1

	N = pars->N;

	// loop over the non-negative (index 0) and the negative (index 1)
	// frequency partition sets. Ordering doesn't matter
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			FPART *partition = fpset->partitions + jj;
			int fstart = partition->start;
			// FFT the S-transform over this partition
			fft(partition->width, dst + fstart * stride, stride);
			// remove partition window from the transformed data
			for ( kk = 0 ; kk < partition->win_len ; kk++ ) {
				int ff = fstart + kk;
				cdiv(dst + ff * stride, partition->window + kk);
			}
		}
	}

	// inverse FFT to recover the signal
	ifft(N, dst, stride);

}


// 2D transform, mostly following the original GFT code.
// Need to write an inverse function to be useful for filtering applications
// Need an explicit mapping between S-transform indices and image partitions(wavenumber and pixel)
// This function compiles, but has not been tested by me (CKW)
DllExport void ngft_2dComplex64(DCMPLX *image, int N, int M, windowFunction *window_fn) {
	int row, col;
	FPCOL *pars;

	if ( image == NULL || N <= 0 || M <= 0 )
		oops( "ngft_2dComplex64", "Invalid argument: image is NULL or N <= 0 or M <= 0" );
	if ( window_fn == NULL )
		window_fn = gaussian;

	pars = ngft_MakePartsAndWindows(N, window_fn);

	for ( row = 0 ; row < N ; row++ )
		ngft_1dComplex64(image + row * N, N, pars, 1);

	if ( M != N ) {
		free( pars );
		pars = ngft_MakePartsAndWindows(M, window_fn);
	}

	for ( col = 0 ; col < M ; col++ )
		ngft_1dComplex64(image + col, M, pars, N);

	free( pars );
}


// free an ILIST structure
DllExport void freeIlist( ILIST *ilist ) {
	if ( ilist != NULL ) {
		if ( ilist->values != NULL && ilist->count > 0 )
			free( ilist->values );
		free( ilist );
	}
}


// get the count of frequency partitions
static int getFreqDim( FPCOL *pars, BOOL all_freqs ) {
	int ii, imax, f_dim;
	// use non-negative (index 0), and optionally, negative (index 1) frequency partition sets
	imax = all_freqs ? 2 : 1;
	for ( f_dim = 0, ii = 0 ; ii < imax ; ii++ ) {
		FPSET *fpset = pars->fpset + ii;
		f_dim += fpset->pcount;
	}
	return f_dim;
}


// get the center-values of the frequency partitions, in order of the partitions
DllExport ILIST *freqPartitionCenters(FPCOL *pars, BOOL all_freqs) {
	int ii, imax, f_count, f_dim, N;
	ILIST *f_centers;
	int *values;

	N = pars->N;
	f_dim = getFreqDim(pars, all_freqs);
	values = calloc(f_dim, sizeof(*values));

	// use non-negative (index 0), and optionally, negative (index 1) frequency partition sets
	imax = all_freqs ? 2 : 1;
	for ( f_count = 0, ii = 0 ; ii < imax ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
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


// get the center-values of the time-partition with the most partitions, in order of the partitions
DllExport ILIST *timePartitionCenters( TPCOL *tpcol ) {
	int ii, t_dim;
	TDSET *tdset_max;
	ILIST *t_centers;
	int *values;

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


// find the index in the DST corresponding to frequency ff and time tt.
static int find_index(FPCOL *pars, TPCOL *tpcol, int ff, int tt, BOOL all_freqs) {
	int ii, imax, ind;
	int N = pars->N;
	int f_target = FREQ_2_INDEX(ff, N);
	int t_target = tt;

	// loop over the non-negative (index 0) and (optionally) the negative (index 1) frequency partition sets
	// Note: This function requires that the frequency partition sets are in normal frequency
	// order, i,e, positive sets: (0, 1, ... , max_f, [Nyquist]) negative sets: (-max_f, ... -1)
	imax = all_freqs ? 2 : 1;
	for ( ind = 0, ii = 0 ; ii < imax ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the frequency partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			if ( fpart->start <= f_target && f_target <= fpart->end ) {
				// frequency is in this partition; now check time partitions
				int kk;
				// first need to find the time-partition set having a decimation that
				// corresponds to this frequency-partition width
				for ( kk = 0 ; kk < tpcol->tdcount ; kk++ ) {
					TDSET *tdset = tpcol->tdsets + kk;
					// check that time-partition decimation and frequency-partition width correspond
					if ( tdset->decimation == N / fpart->width ) {
						// have a match, so now check all time partitions with this decimation
						int ll;
						for ( ll = 0 ; ll < tdset->pcount ; ll++, ind++ ) {
							// search for the partition with matching time
							TPART *tpart = tdset->partitions + ll;
							if ( tpart->start <= t_target && t_target <= tpart->end )
								return ind;	// cool
						}
					}
				}
			} else {
				ind += fpart->width;	// fwidth = N / decimation = number of time partitions for this frequency partition
			}
		}
	}
	return -1;	// bummer
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
DllExport DIMAGE *ngft_1d_InterpolateNN(DCMPLX *signal, FPCOL *pars, TPCOL *tpcol, int M,
												BOOL by_part, BOOL all_freqs, BOOL make_ind_map) {
	int ii, f_dim, nf_dim, t_dim, nt_dim, img_len, N;
	ILIST *f_centers = NULL, *t_centers = NULL;
	DCMPLX *image_arr;
	DIMAGE *image;
	double t_factor, f_factor;
	BOOL down_sample, free_tpcol = FALSE;

	if ( (signal == NULL && make_ind_map == FALSE) || pars == NULL )
		oops( "ngft_1d_InterpolateNN", "Invalid argument: Signal and/or pars is NULL" );

	if ( tpcol == NULL ) {
		tpcol = ngft_TimePartitions( pars );
		free_tpcol = TRUE;
	}

	N = pars->N;

	t_dim = by_part ? getTimeDim(tpcol) : N;
	if ( M <= 0 )
		M = t_dim;
	M = MIN( M, t_dim );	// don't allow augmentation
	down_sample = (M < t_dim);
	t_factor = (double)M / t_dim;
	nt_dim = M;

	f_dim = by_part ? getFreqDim(pars, all_freqs) : all_freqs ? N : N / 2;
	f_factor = down_sample ? MIN( (all_freqs ? 1 : 2) * t_factor, 1 ) : 1;
	nf_dim = ROUND(f_factor * f_dim);

	img_len = nf_dim * nt_dim;
	image_arr = calloc(img_len, sizeof(*image_arr));

	if ( by_part ) {
		f_centers = freqPartitionCenters( pars, all_freqs );
		t_centers = timePartitionCenters( tpcol );
	} else {
		f_centers = calloc( 1, sizeof( *f_centers ) );
		f_centers->values = calloc( f_dim, sizeof( *(f_centers->values) ) );
		f_centers->count = f_dim;
		t_centers = calloc( 1, sizeof( *t_centers ) );
		t_centers->values = calloc( t_dim, sizeof( *(t_centers->values) ) );
		t_centers->count = t_dim;
	}

	// get image
	for ( ii = 0 ; ii < nf_dim ; ii++ ) {
		int jj;
		int is = down_sample ? ROUND(ii / f_factor) : ii;
		int ff = by_part ? f_centers->values[is] : INDEX_2_FREQ( is, N );
		if ( ! by_part )
			f_centers->values[ii] = ff;
		for ( jj = 0 ; jj < nt_dim ; jj++ ) {
			int img_ind;
			int js = down_sample ? ROUND(jj / t_factor) : jj;
			int tt = by_part ? t_centers->values[js] : js;
			if ( ! by_part )
				t_centers->values[jj] = tt;
			int kk = find_index(pars, tpcol, ff, tt, all_freqs);
			if ( kk < 0 || kk >= N )
				oops( "ngft_1d_interpolateNN", "Application exception: can't find image index" );
			img_ind = ii * nt_dim + jj; // image stored by row, with time on the horizontal axis, and frequency on the vertical
			if ( make_ind_map ) {
				image_arr[img_ind].r = kk;	// image value is just the DST index
				image_arr[img_ind].i = 0;
			} else {
				memcpy( image_arr + img_ind, signal + kk, sizeof( *image_arr ) );
			}
		}
	}

	if ( free_tpcol )
		ngft_FreeTimePartitions(tpcol);

	image = calloc( 1, sizeof( *image ) );
	image->img = calloc(1, sizeof(*(image->img)));
	image->img->values = image_arr;
	image->img->count = img_len;
	image->x_centers = t_centers;
	image->y_centers = f_centers;

	return image;
}

