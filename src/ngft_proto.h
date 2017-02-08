#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );

// Windows
DllExport DCMPLX *gaussian(int N, int freq);
DllExport DCMPLX *box(int N, int freq);

// GFT partition and window constructors and destructors
DllExport FPCOL *ngft_FrequencyPartitions(int N, double epsilon);
DllExport void ngft_FreeFreqPartitions(FPCOL *pars);
DllExport void freeIlist(ILIST *ilist);
DllExport void freeDImage(DIMAGE *image);

// 1D GFT Functions
DllExport void ngft_1dComplex64(DCMPLX **signal, int *N, FPCOL **fpcol, windowFunction *window_fn, int stride);
DllExport void ngft_1dComplex64Inv( DCMPLX **dst, FPCOL **fpcol, windowFunction *window_fn, int stride );

#ifdef IS_PORTED
// 2D GFT Functions
DllExport void ngft_2dComplex64(DCMPLX *image, int N, int M, windowFunction *window_fn);
#endif

// Interpolation functions
DllExport DIMAGE *ngft_1d_InterpolateNN(DCMPLX *signal, FPCOL *pars, int M,
														BOOL by_part, BOOL all_freqs, BOOL make_ind_map);
DllExport ILIST *freqPartitionCenters(FPCOL *pars, BOOL all_freqs);
DllExport ILIST *timePartitionCenters(FPCOL *pars);
