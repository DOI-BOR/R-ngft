/*
 *  gft.c
 *  GFT Framework
 *
 *  Created by Robert Brown on 30/05/08.
 *	This software is copyright © 2010 UTI Limited Partnership.  
 *	The original authors are Robert A. Brown, M. Louis Lauzon 
 *	and Richard Frayne.  This software is licensed in the terms 
 *	set forth in the “FST License Notice.txt” file, which is 
 *	included in the LICENSE directory of this distribution.
 *
 *	Substantially modified by Chris Wood on 1/26/2017
 *
 */

#include "gft.h"

#ifdef DllImport
#  undef DllImport
#  define DllImport DllExport
# endif
#include "gft_proto.h"

#include "fftw3.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>


static void fft(int N, DCMPLX *in, int stride) {
	fftw_plan p;
	p = fftw_plan_many_dft(1, &N, 1, (fftw_complex *)in, NULL, stride, 0, (fftw_complex *)in, NULL, stride, 0, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);	
}


static void ifft(int N, DCMPLX *in, int stride) {
	int ii;
	fftw_plan p;
	p = fftw_plan_many_dft(1, &N, 1, (fftw_complex *)in, NULL, stride, 0, (fftw_complex *)in, NULL, stride, 0, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	for ( ii = 0 ; ii < N ; ii++ ) {
		in[ii*stride].r /= N;
		in[ii*stride].i /= N;
	}
}	


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


// efficient circular shift of an array of complex doubles
// Shift is to the right if amount > 0, and to the left if amount < 0
static void shift(DCMPLX *sig, int N, int amount) {
	int len;
	DCMPLX *tmp;
	int s1, s2, s3, s4;
	unsigned tsiz = sizeof( *tmp );

	if ( sig == NULL || N == 0 || amount == 0 )
		return; // do nothing

	amount %= N;	// circular shifts are modulo N

	if ( amount < 0 )
		amount = N + amount;	// circular left shift is the same as a right shift of N - |amount|

	// shift by using a single memmove of the big piece, and two memcpy of the small piece.
	// shift right by |amount|
	if ( amount < N / 2 ) {
		len = amount; // length of small piece is |amount|
		// big piece begins at offset s3 = 0, small piece begins at offset s1 = N - |amount|
		// big piece destination is offset s2 = |amount|, small piece destination is offset s4 = 0
		s1 = N - len; s3 = 0;
		s2 = len; s4 = 0;	
	} else {
		len = N - amount; // length of small piece is N - |amount|
		// big piece begins at offset s3 = |amount|, small piece begins at offset s1 = 0
		// big piece destination is offset s2 = 0, small piece destination is offset s4 = N - |amount|
		s1 = 0; s3 = len;
		s2 = 0; s4 = N - len;
	}
	tmp = malloc( len * tsiz );
	memcpy( tmp, sig + s1 * tsiz, len * tsiz );
	memmove( sig + s2 * tsiz, sig + s3 * tsiz, (N - len) * tsiz );
	memcpy( sig + s4 * tsiz, tmp, len * tsiz );
	free( tmp );
}


DllExport DCMPLX *gaussian(int N, int freq) {
	int ii;
	// Original code set x_max = 0.5, but exp(-x) underflows for x > 708 on IEEE hardware. Thus we
	//		need x_max * freq < sqrt(2 * 708) = 38, to avoid underflow. The largest center frequency
	//		in a dyadic scaling is about freq_max = 3 * N / 8, so x_max must be scaled such that
	//		x_max * freq_max < sqrt(2 * 708), or x_max < 100 / N. Select x_max < 10 / N.
	double x, x_max = MIN( 0.5, 10. / N );
	BOOL N_is_odd = (N % 2 == 0 ? FALSE : TRUE);
	double dx = 2 * x_max / (N_is_odd ? N - 1 : N);
	int i0 = (N_is_odd ? 0 : 1);
	DCMPLX *win = calloc(N, sizeof( *win ));
	double win_sum;

	for ( win_sum = 0 , ii = i0 ; ii < i0 + N ; ii++ ) {
		x = (ii * dx - x_max) * freq;
		win[ii].r = exp(-x * x / 2) * ABS(freq) / sqrt(2*M_PI);
		win[ii].i = 0;
		win_sum += win[ii].r;
	}

	// Make sure the window area is 1.0
	for ( ii = 0; ii < N ; ii++ )
		win[ii].r /= win_sum;

	// shift window so that 0 time is at index 0, and negative times start at index N/2+1 (standard DFT convention)
	shift(win, N, N/2 + 1);

	// take the DFT of an even Real function
	fft(N,win,1);

	return(win);
}


DllExport DCMPLX *box(int N, int freq) {
	int ii;
	DCMPLX *win = calloc(N, sizeof( *win ));
	for ( ii = 0 ; ii < N ; ii++ ) {
		win[ii].r = 1.0;
		win[ii].i = 0.0;
	}

	return(win);
}


// Dyadic partitioning of a DFT frequency array of length N, where N can be odd or even,
//	and does not need to be a power of 2. The partitioning will span the array, with
//	no overlap in partitions.
DllExport FPCOL *ngft_DyadicPartitions(int N) {
	int fs, pcount, n_pcount;
	FPART *p, *partitions, *np, *n_partitions;
	FPCOL *partition_set;

	if ( N <= 0 )
		oops( "ngft_DyadicPartitions", "Invalid argument: N <= 0" );

	for ( partitions = NULL, n_partitions = NULL, pcount = 0, n_pcount = 0, fs = 0 ; fs < N / 2 ; fs *= 2 ) {
		// handle non-negative frequencies
		partitions = realloc(partitions, ++pcount * sizeof(*partitions));
		p = partitions + pcount - 1;
		p->start = NONNEG_F_IND(fs, N);
		p->end = NONNEG_F_IND(2 * fs - 1, N);
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
// This code has not been tested much (or even thought about much) by me (CWood)
DllExport FPCOL *ngft_1dMusicPartitions(int N, double samplerate, int cents) {
	int ii, ii_max, pcount, n_pcount;
	FPART *p, *partitions, *np, *n_partitions;
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

	partition_set = calloc(1, sizeof(*partition_set));
	partition_set->N = N;
	partition_set->fpset[0].partitions = partitions;
	partition_set->fpset[0].pcount = pcount;
	partition_set->fpset[1].partitions = n_partitions;
	partition_set->fpset[1].pcount = n_pcount;
}


// free the space
DllExport void ngft_FreeFreqPartitions(FPCOL *pars) {
	int ii;
	if ( pars == NULL )
		return;

	for ( ii = 0 ; ii < 2 ; ii++ ) {
		if ( pars->fpset[ii].partitions != NULL && pars->fpset[ii].pcount > 0 ) {
			if ( pars->fpset[ii].partitions->window != NULL && pars->fpset[ii].partitions->win_len > 0 )
				free( pars->fpset[ii].partitions->window );
			free( pars->fpset[ii].partitions );
		}
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

	// set values implied by this frequency partition
	tdset->decimation = fpwidth;	// critical point: time-axis decimation factor = frequency-partition width
	tdset->partitions = NULL;
	tdset->pcount = 0;

	// add time partitions for this decimation factor, rounding as necessary to fill N
	for ( tc = 0, tc_exact = 0, ii = 0 ; ii < fpwidth ; ii++, tc_exact += exact_width ) {
		TPART* partition;
		tdset->partitions = realloc(tdset->partitions, ++(tdset->pcount) * sizeof(*(tdset->partitions)));
		partition = tdset->partitions + ii;
		partition->start = tc;
		partition->width = ROUND( tc_exact + exact_width - tc );
		partition->center = partition->start + partition->width / 2;
		partition->end = partition->start + partition->width - 1;
		tc += partition->width;
	}
}


DllExport TPCOL *ngft_TimePartitions(FPCOL *pars) {
	int ii, N;
	TPCOL *tpcol;

	if ( pars == NULL )
		oops("ngft_TimePartitions", "Invalid argument: pars is NULL");
	N = pars->N;

	tpcol = calloc( 1, sizeof( *tpcol ) ); // note: everything in tpcol initialized to NULL or 0

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the frequency partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			BOOL found = FALSE;
			if ( tpcol->tdcount > 0 ) {
				// search existing time-partition sets for a time-axis decimation factor equal to
				// this frequency-partition width
				int kk;
				for ( kk = 0 ; kk < tpcol->tdcount ; kk++ ) {
					if ( tpcol->tdsets[kk].decimation == fpart->width ) {
						found = TRUE;
						break;
					}
				}
			}
			if ( ! found ) {
				// no time-partition sets with the correct decimation found, so create one and fill,
				// using decimation implied by this frequency-partition
				tpcol->tdsets = realloc(tpcol->tdsets, ++(tpcol->tdcount) * sizeof(*(tpcol->tdsets)));
				fill_tdset(N, fpart, tpcol->tdsets);
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
		free( tdset );
	}
	free( tpcol );
}


DllExport FPCOL *ngft_MakePartsAndWindows(int N, windowFunction *window_fn) {
	if ( N <= 0 )
		oops( "ngft_MakePartsAndWindows", "Invalid argument: N <= 0" );
	if ( window_fn == NULL )
		window_fn = gaussian;

	FPCOL *partition_set = ngft_DyadicPartitions( N );
	ngft_AddWindowsToParts(partition_set, window_fn);
	return partition_set;
}


// Create windows and associate with corresponding partitions
// This should be the one and only place in the ngft library where windows are created
DllExport void ngft_AddWindowsToParts(FPCOL *pars, windowFunction *window_fn) {
	int ii;
	int N;

	if ( pars == NULL )
		oops("ngft_AddWindowsToParts", "Invalid argument: pars is NULL");
	if ( window_fn == NULL )
		window_fn = gaussian;

	N = pars->N;

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			double sum, norm;
			FPART *partition = fpset->partitions + jj;
			// create a centered window of length N in the frequency domain, possibly using sigma ~ 1/fcenter
			int fcenter = partition->center;
			DCMPLX *full_win = window_fn(N, fcenter);

			// subset out the part of the window that overlaps this partition. Note: this code doesn't support
			//	overlapping windows, but it would be straightforward to implement
			int win_len = partition->width;
			DCMPLX *win = calloc( win_len, sizeof( *win ) );
			memcpy( win, full_win + partition->start, win_len * sizeof( *win ) );

			// need to renormalize window so that sum is 1
			for ( sum = 0, kk = 0; kk < win_len ; kk++ )
				sum += modulus(win + kk);
			norm = 1 / sum;
			for ( sum = 0, kk = 0; kk < win_len ; kk++ )
				cmulByReal(win+kk, norm);

			partition->window = win;
			partition->win_len = win_len;
			free( full_win );
		}
	}
}


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

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
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


DllExport void ngft_1dComplex64Inv( DCMPLX *signal, FPCOL *pars, int stride ) {
	int ii, N;

	if ( signal == NULL )
		oops( "ngft_1dComplex64Inv", "Invalid argument: Signal is NULL" );
	if ( pars == NULL )
		oops( "ngft_1dComplex64Inv", "Invalid argument: Pars is NULL" );

	stride = MAX( 1, stride );	// if stride <= 0, set to default value of 1

	N = pars->N;

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			FPART *partition = fpset->partitions + jj;
			int fstart = partition->start;
			// FFT the S-transform over this partition
			fft(partition->width, signal + fstart * stride, stride);
			// remove partition window from the transformed data 
			for ( kk = 0 ; kk < partition->win_len ; kk++ ) {
				int ff = fstart + kk;
				cdiv(signal + ff * stride, partition->window + kk);
			}
		}
	}

	// inverse FFT to recover the signal
	ifft(N, signal, stride);

}


// 2D transform, mostly following the original GFT code.
// This code has not been tested (or even thought about much) by me (CWood)
// Need to write an inverse function to be useful for filtering applications
// Need an explicit mapping between S-transform indices and image partitions(wavenumber and pixel)
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


static int find_index(FPCOL *pars, TPCOL *tpcol, int ff, int tt, BOOL down_sample, double factor) {
	int ii, ind;
	int f_target = (down_sample ? ROUND(factor * ff) : ff);
	int t_target = (down_sample ? ROUND(factor * tt) : tt);

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ind = 0, ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		// loop over the partitions in this set
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *fpart = fpset->partitions + jj;
			if ( fpart->start <= f_target && f_target <= fpart->end ) {
				int kk;
				// find the time-partitions corresponding to this frequency-partition width
				for ( kk = 0 ; kk < tpcol->tdcount ; kk++ ) {
					// search for the partition with matching time
					TDSET *tdset = tpcol->tdsets + kk;
					if ( tdset->decimation == fpart->width ) {
						int ll;
						for ( ll = 0 ; ll < tdset->pcount ; ll++, ind++ ) {
							TPART *tpart = tdset->partitions + ll;
							if ( tpart->start <= t_target && t_target <= tpart->end )
								return ind;	// Eureka!
						}
					}
				}
			} else {
				ind += fpart->width;	// fwidth is also the number of time partitions for this frequency partition
			}
		}
	}
	return -1;
}


// Make a 2D time-frequency image using the original linear spacing in time and frequency
// using a nearest neighbor sampling of the fast S transform.
// If the original signal was N points long, the image will be NxN. Optionally, if 0 < M <=N
// is specified, then the image will be down-sampled to MxM (can set M<=0 to use N)
// Since the transform only has N points, the NxN image will be highly undersampled, and appear
// pixelated.
DllExport DIMAGE *ngft_1d_interpolateNN(DCMPLX *signal, FPCOL *pars, TPCOL *tpcol, int M) {
	int N, ff;
	DCMPLX *image_arr;
	DIMAGE *image;
	BOOL free_tpcol = FALSE;
	double factor;
	BOOL down_sample;
	
	if ( signal == NULL || pars == NULL )
		oops( "ngft_1d_interpolateNN", "Invalid argument: Signal and/or pars is NULL" );

	if ( tpcol == NULL ) {
		tpcol = ngft_TimePartitions( pars );
		free_tpcol = TRUE;
	}

	N = pars->N;
	if ( M <= 0 )
		M = N;
	M = MIN( M, N );	// don't allow augmentation
	down_sample = (M < N);
	factor = (double)N / M;

	image_arr = calloc(sizeof(*image_arr), M * M);

	// get image
	for ( ff = 0 ; ff < M ; ff++ ) {
		int tt;
		for ( tt = 0 ; tt < M ; tt++ ) {
			int kk = find_index(pars, tpcol, ff, tt, down_sample, factor);
			if ( kk < 0 || kk >= N )
				oops( "ngft_1d_interpolateNN", "Application error: can't find image index" );
			memcpy(image_arr + (M - 1 - ff) * M + tt, signal + kk, sizeof(*image));
		}
	}

	if ( free_tpcol )
		ngft_FreeTimePartitions(tpcol);

	image = calloc( 1, sizeof( *image ) );
	image->img = calloc(1, sizeof(*(image->img)));
	image->img->values = calloc(1, sizeof(*(image->img->values)));
	image->img->count = M * M;
	image->ht = M;
	image->wd = M;

	return image;
}


// get the count of frequency partitions
static int getFreqDim( FPCOL *pars ) {
	int ii, f_dim;
	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( f_dim = 0, ii = 0 ; ii < 2 ; ii++ ) {
		FPSET *fpset = pars->fpset + ii;
		f_dim += fpset->pcount;
	}
	return f_dim;
}


// get the center-values of the frequency partitions
DllExport ILIST *getFreqCenters(FPCOL *pars) {
	int ii, f_count, f_dim;
	ILIST *f_centers;
	int *values;

	f_dim = getFreqDim(pars);
	values = calloc(f_dim, sizeof(*values));

	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( f_count = 0, ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = pars->fpset + ii;
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			FPART *part = fpset->partitions + jj;
			values[f_count++] = part->center;
		}
	}

	f_centers = calloc( 1, sizeof( *f_centers ) );
	f_centers->values = values;
	f_centers->count = f_count;

	return f_centers;
}


DllExport void freeIlist( ILIST *ilist ) {
	if ( ilist != NULL ) {
		if ( ilist->values != NULL && ilist->count > 0 )
			free( ilist->values );
		free( ilist );
	}
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


// get the center-values of the maximal time-partition
DllExport ILIST *getTimeCenters( TPCOL *tpcol ) {
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


// Make a 2D time-frequency image using a linear spacing in time and log spacing
// in frequency (non-negative frequencies only) using a nearest neighbor sampling
// of the fast S transform. The frequencies will be the centers of the frequency
// partitions, and the times will be the centers of the time partition with the greatest
// number of partitions. If the original signal was N points long, the image will
// be log2N x log2N. Optionally, if 0 < M <= log2N is specified, then the image will be
// down-sampled to MxM (can set M<=0 to use default)
DllExport DIMAGE *ngft_1d_logfInterpolateNN(DCMPLX *signal, FPCOL *pars, TPCOL *tpcol, int M) {
	int dim, ii, f_dim, t_dim, N;
	ILIST *f_centers;
	ILIST *t_centers;
	DCMPLX *image_arr;
	DIMAGE *image;
	BOOL free_tpcol = FALSE;
	double factor;
	BOOL down_sample;

	if ( signal == NULL || pars == NULL )
		oops( "ngft_1d_logfInterpolateNN", "Invalid argument: Signal and/or pars is NULL" );

	if ( tpcol == NULL ) {
		tpcol = ngft_TimePartitions( pars );
		free_tpcol = TRUE;
	}

	f_dim = getFreqDim(pars);
	t_dim = getTimeDim(tpcol);
	if ( f_dim != 2 * t_dim )
		oops( "ngft_1d_logfInterpolateNN", "Application error: f_dim != t_dim" );
	dim = t_dim;
	if ( M <= 0 )
		M = dim;
	M = MIN( M, dim );	// don't allow augmentation
	down_sample = (M < dim);
	factor = (double)dim / M;

	N = pars->N;

	image_arr = calloc(sizeof(*image_arr), M * M);

	f_centers = getFreqCenters( pars );
	t_centers = getTimeCenters( tpcol );

	// get image
	for ( ii = 0 ; ii < dim ; ii++ ) {
		int jj;
		for ( jj = 0 ; jj < dim ; jj++ ) {
			int ff = f_centers->values[ii];
			int tt = t_centers->values[jj];
			int kk = find_index(pars, tpcol, ff, tt, down_sample, factor);
			if ( kk < 0 || kk >= N )
				oops( "ngft_1d_logfInterpolateNN", "Application error: can't find image index" );
			memcpy(image_arr + (M - 1 - ff) * M + tt, signal + kk, sizeof(*image));
		}
	}

	if ( free_tpcol )
		ngft_FreeTimePartitions(tpcol);

	image = calloc( 1, sizeof( *image ) );
	image->img = calloc(1, sizeof(*(image->img)));
	image->img->values = calloc(1, sizeof(*(image->img->values)));
	image->img->count = M * M;
	image->ht = M;
	image->wd = M;

	return image;
}


DllExport void freeDImage( DIMAGE *image ) {
	if ( image != NULL ) {
		if ( image->img != NULL ) {
			if ( image->img->values != NULL && image->img->count > 0 )
				free( image->img );
			free( image->img );
		}
		free( image );
	}
}


