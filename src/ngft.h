#pragma once

#include "common.h"
#include "fftw3.h"

/*
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

// collection of double values
typedef struct {
	double *values;
	int count;
} DLIST;

// collection of integer values
typedef struct {
	int *values;
	int count;
} ILIST;

typedef struct {
	double r;	/* real part */
	double i;	/* imaginary part */
} DCMPLX;

// collection of DCMPLX values
typedef struct {
	DCMPLX *values;
	int count;
} DCLIST;

typedef DCMPLX *(windowFunction)(int,int);

// In the standard DFT array layout of length N, frequencies are ordered as follows (parentheses refer
// to negative frequencies):
//		0, 1, 2, ..., N/2, N/2+1 (-N/2 for N odd, or -(N/2-1) for N even), ..., N-2 (-2), N-1 (-1)
//	If N is odd, this notation assumes the value of N/2 is truncated per the integer division rules of C.
//	For this layout:
//		index(f >= 0) = f, which spans [0, N/2]    (even or odd N)
//		index(f < 0) = N - |f|, which spans [N/2 + 1, N - 1]    (even or odd N)

// Get the bounded index for a non-negative frequency f in an array of length N, where N can be odd or even
#define NONNEG_F_IND(f,N)	MIN(MAX(0, (f)), (N)/2)

// Get the bounded index for a negative frequency -|f| in an array of length N, where N can be odd or even
#define NEG_F_IND(f,N)	MIN(MAX((N) - ABS(f), (N)/2 + 1), (N) - 1)

// convert an index to a frequency; assumes 0 <= index < N, and -(N+1)/2 < f <= N/2
#define INDEX_2_FREQ(i,N)	((i) <= (N)/2 ? (i) : (i) - (N))
#define FREQ_2_INDEX(f,N) ((f) < 0 ? (N) - ABS(f) : (f))

// Macros to test if an integer is even or odd
#define IS_EVEN(i)	((i) % 2 == 0)
#define IS_ODD(i)	((i) % 2 == 1)

// structure defining a time partition
typedef struct {
	int start;	// index into N of partition starting point
	int center;	// index into N of partition center point (odd width), or next point after center-value (even width)
	int end;	// index into N of partition ending point
	int width;	// number of points in partition
} TPART;

// structure defining the set of time partitions obtained for a particular decimation factor
typedef struct {
	TPART *partitions;	// pointer to partitions
	int pcount;	// number of partitions
	double decimation;	// the decimation (scaling) factor used to obtain these partitions
} TDSET;

// structure defining a collection of time partitions
typedef struct {
	int N;	// the length spanned by the partitions
	TDSET *tdsets;	// pointer to sets of time partitions with common decimation factors
	int tdcount;	// number of tdsets
} TPCOL;

// structure defining a frequency partition
// Note: a frequency partition must not overlap any other partition, however, the windows may overlap
typedef struct {
	int start;	// index into N of partition starting point
	int center;	// index into N of partition center point (odd width), or next point after center-value (even width)
	int end;	// index into N of partition ending point
	int width;	// number of points in partition
	int win_start;	// index into N of frequency window starting point
	int win_end;	// index into N of frequency window ending point
	int win_len;	// length of frequency window to be applied: for non-overlapping windows (fast GFT), win_len = width
	TDSET *tdset;	// pointer to the set of GFT time partitions for this frequency partition
	DCLIST *gft;	// pointer to GFT(t,f) for f = this center frequency. Length is win_len (= tdset->pcount)
} FPART;

// set of frequency partitions
// Note: the set of partitions must defined to completely cover the set of discrete frequencies without
// overlap. However, a window is defined for each partition that can overlap other partitions.
typedef struct {
	FPART *partitions;	// pointer to partitions
	int pcount;	// number of partitions
} FPSET;

// enum defining partitioning types
typedef enum { FP_UNKNOWN = 0, FP_DYADIC, FP_EDO, FP_FW } FreqPartitionType;
// enum defining window types
typedef enum { FWT_UNKNOWN = 0, FWT_GAUSSIAN, FWT_BOX } FreqWindowType;

// structure defining a collection of frequency partitions
typedef struct {
	// This structure stores frequency partitions, allowing all pertinent info to be passed or returned using a single pointer.
	// The non-negative and negative partitions are managed separately because they won't necessarily have the symmetry that the
	// frequency array does when N is not a power of 2. For example, consider an array of length 6, where the non-neg sets are
	// {0}, {1}, and {2,3}, but the negative sets are {-1} and {-2} (rather {-2,-3}, because there is no f = -3). Even when N is
	// a power of 2, separate partitions allows for simpler code, and there's no performance penalty.
	int N;	// the data length spanned by the partitions
	double epsilon; // fractional amount of overlap for frequency partitions; 0 if none
	FreqPartitionType partition_type;	// type of frequency partition (e.g., dyadic, equal-temperament, etc.)
	FreqWindowType window_type;	// type of frequency-domain window (e.g., Gaussian, box, etc.)
	FPSET fpset[2];	// array of partitions for non-negative (index 0) and negative (index 1) frequencies
	TPCOL *tpcol;	// pointer to a time-partition collection used by this frequency partitioning
} FPCOL;

// Image
typedef struct {
	ILIST *y_centers;
	ILIST *x_centers;
	DCLIST *img;
} DIMAGE;
