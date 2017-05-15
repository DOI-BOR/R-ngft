#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );

// GFT partition constructors and destructors
DllExport FPCOL *ngft_FrequencyPartitions(int N, double epsilon, FreqPartitionType ptype, FreqWindowType wtype,
																					int W, int f_ref, int T);
DllExport void ngft_FreeFreqPartitions(FPCOL *fpcol);
DllExport void freeIlist(ILIST *ilist);
DllExport void freeDClist(DCLIST *ilist);
DllExport void freeDImage(DIMAGE *image);

// 1D GFT Functions
DllExport FPCOL *ngft_1dComplex64(DCLIST *sig, double epsilon, FreqPartitionType ptype, FreqWindowType wtype,
																	int W, int f_ref, int T);
DllExport DCLIST *ngft_makeGftArray(FPCOL *fpcol);
DllExport DCLIST *ngft_1dComplex64Inv( FPCOL *fpcol );

#ifdef IS_PORTED
// 2D GFT Functions
DllExport void ngft_2dComplex64(DCMPLX *image, int N, int M, windowFunction *window_fn);
#endif

// Interpolation functions
DllExport DIMAGE *ngft_1d_InterpolateNN(FPCOL *fpcol, int M, BOOL by_part, BOOL all_freqs);
DllExport ILIST *freqPartitionCenters(FPCOL *fpcol, BOOL all_freqs);
DllExport ILIST *timePartitionCenters(FPCOL *fpcol);

// support functions
DllExport DCLIST *ngft_makeGftArray(FPCOL *fpcol);
DllExport void ngft_unpackGftArray(DCLIST *gft, FPCOL *fpcol);
