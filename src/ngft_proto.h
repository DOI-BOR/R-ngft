#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );

// Windows
DllExport DCMPLX *gaussian(int N, int freq);
DllExport DCMPLX *box(int N, int freq);

// GFT partition and window constructors and destructors
DllExport FPCOL *ngft_DyadicPartitions(int N);
DllExport FPCOL *ngft_1dMusicPartitions(int N, double samplerate, int cents);
DllExport FPCOL *ngft_MakePartsAndWindows(int N, windowFunction *window_fn);
DllExport void ngft_AddWindowsToParts(FPCOL *pars, windowFunction *window_fn);
DllExport TPCOL *ngft_TimePartitions(FPCOL *pars);
DllExport void ngft_FreeFreqPartitions(FPCOL *pars);
DllExport void ngft_FreeTimePartitions(TPCOL *tpcol);
DllExport void freeIlist(ILIST *ilist);
DllExport void freeDImage(DIMAGE *image);

// 1D GFT Functions
DllExport void ngft_1dComplex64(DCMPLX *signal, int N, FPCOL *pars, int stride);
DllExport void ngft_1dComplex64Inv(DCMPLX *dst, FPCOL *pars, int stride);

// 2D GFT Functions
DllExport void ngft_2dComplex64(DCMPLX *image, int N, int M, windowFunction *window_fn);

// Interpolation functions
DllExport DIMAGE *ngft_1d_InterpolateNN(DCMPLX *signal, FPCOL *pars, TPCOL *tpcol, int M,
														BOOL by_part, BOOL all_freqs, BOOL make_ind_map);
DllExport ILIST *freqPartitionCenters(FPCOL *pars, BOOL all_freqs);
DllExport ILIST *timePartitionCenters(TPCOL *tpcol);
