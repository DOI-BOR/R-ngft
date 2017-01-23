#pragma once

extern void oops ( char *, char * );
extern void smsg ( char *, char * );

// Windows
DllExport void gaussian(double *win, int N, int freq);
DllExport void box(double *win, int N, int freq);

// GFT partition and window screen generators
DllExport int gft_1dSizeOfPartitions(unsigned int N);
DllExport int *gft_1dPartitions(unsigned int N);
DllExport int *gft_1dMusicPartitions(unsigned int N, float samplerate, int cents);
DllExport double *windows(int N, windowFunction *window);
DllExport double *windowsFromPars(int N, windowFunction *window, int *pars);

// 1D GFT Functions
DllExport void gft_1dComplex64(double *signal, unsigned int N, double *win, int *pars, int stride);

// 2D GFT Functions
DllExport void gft_2dComplex64(double *signal, unsigned int N, unsigned int M, windowFunction *window);

// Utility Functions
DllExport void gft_1d_shift(double *signal, unsigned int N, unsigned int shiftBy);

// Interpolation functions
DllExport double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M);
