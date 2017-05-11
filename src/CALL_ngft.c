#include "ngft.h"
#include "ngft_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

DllExport SEXP CALLngft_1dComplex64(SEXP ts_d, SEXP dt_d, SEXP eps_d, SEXP ptype_s, SEXP wtype_s,
																		SEXP image_dim_i, SEXP by_part_l, SEXP all_freqs_l)
{
	int ii, jj, pcnt = 0, nf, nt, *ip;
	DCMPLX *cts;
	DCLIST *signal, *gft;
	DIMAGE *image;
	Rcomplex *ctp;
	FPCOL *partitions;
	SEXP ret_l, names_s;
	SEXP deltat_d, epsilon_d;
	SEXP part_name_s, win_name_s;
	SEXP gft_c, cts_img_c;
	SEXP img_wd_i, img_ht_i;
	SEXP freq_centers_i, time_centers_i;
	SEXP ts_len_i;

	double *ts = REAL(ts_d);
	int n_samples = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	double epsilon = ! R_FINITE(asReal(eps_d)) ? 0 : REAL(eps_d)[0];
	FreqPartitionType ptype = asChar(ptype_s) == NA_STRING ? FP_DYADIC :
			strncasecmp(CHAR(STRING_ELT(ptype_s,0)), "d", 1) == 0 ? FP_DYADIC : /* Dyadic */
			FP_DYADIC; /* silently ignore anything else, and use Dyadic */
	FreqWindowType wtype = asChar(wtype_s) == NA_STRING ? FWT_GAUSSIAN :
			strncasecmp(CHAR(STRING_ELT(wtype_s,0)), "g", 1) == 0 ? FWT_GAUSSIAN : /* Gaussian */
			FWT_GAUSSIAN; /* silently ignore anything else, and use Gaussian */
	int image_dim = asInteger(image_dim_i) == NA_INTEGER ? -1 : INTEGER(image_dim_i)[0];
	BOOL by_part = asLogical( by_part_l ) == NA_LOGICAL ? TRUE : LOGICAL( by_part_l )[0];
	BOOL all_freqs = asLogical( all_freqs_l ) == NA_LOGICAL ? FALSE : LOGICAL( all_freqs_l )[0];

	if ( n_samples < 3 )
		oops("CALLngft_1dComplex64", "time series must have at least 3 points");
	if ( dt <= 0 )
		oops("CALLngft_1dComplex64", "dt must be positive");

	// copy the data
	if ( (cts = calloc( n_samples, sizeof( *cts ) )) == NULL )
		oops( "CALLngft_1dComplex64","can't get space" );
	for ( ii = 0 ; ii < n_samples ; ii++ )
		cts[ii].r = ts[ii];

	// Call 1D GFT Function
	signal = calloc(1, sizeof(*signal));
	signal->count = n_samples;
	signal->values = cts;
	partitions = ngft_1dComplex64(signal, epsilon, ptype, wtype);
	freeDClist(signal); // also frees space pointed to by cts

	// get complex image
	image = ngft_1d_InterpolateNN(partitions, image_dim, by_part, all_freqs);
	nf = image->y_centers->count;
	nt = image->x_centers->count;

	// get the gft as a linear array - only used to pass gft to inverse function)
	gft = ngft_makeGftArray(partitions);
	ngft_FreeFreqPartitions(partitions);

	/* allocate space for R structures and copy output */

	// gft array (generally not useful by itself, except to pass to inverse function)
	gft_c = PROTECT(allocVector(CPLXSXP, gft->count)); pcnt++;
	ctp = COMPLEX(gft_c);
	for ( ii = 0 ; ii < gft->count ; ii++ ) {
		ctp[ii].r = gft->values[ii].r;
		ctp[ii].i = gft->values[ii].i;
	}
	freeDClist(gft);

	// Image, with frequencies on the vertical axis and times on the horizontal axis
	cts_img_c = PROTECT(allocMatrix(CPLXSXP, nf, nt)); pcnt++;
	ctp = COMPLEX(cts_img_c);
	for ( ii = 0 ; ii < nf ; ii++ ) {
		for ( jj = 0 ; jj < nt ; jj++ ) {
			int kk = ii * nt + jj;	// NGFT stores images by row
			int ll = jj * nf + ii;	// R stores images by column
			ctp[ll].r = image->img->values[kk].r;
			ctp[ll].i = image->img->values[kk].i;
		}
	}

	// image width and height
	img_wd_i = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(img_wd_i)[0] = nt; // time
	img_ht_i = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(img_ht_i)[0] = nf; // frequency

	// frequency centers
	freq_centers_i = PROTECT(allocVector(INTSXP, nf)); pcnt++;
	ip = INTEGER( freq_centers_i );
	for ( ii = 0 ; ii < nf ; ii++ )
		ip[ii] = image->y_centers->values[ii];

	// time centers
	time_centers_i = PROTECT(allocVector(INTSXP, nt)); pcnt++;
	ip = INTEGER( time_centers_i );
	for ( ii = 0 ; ii < nt ; ii++ )
		ip[ii] = image->x_centers->values[ii];

	freeDImage(image);

	// other returned items
	ts_len_i = PROTECT(allocVector(INTSXP, 1)); pcnt++; INTEGER(ts_len_i)[0] = n_samples;
	deltat_d = PROTECT(allocVector(REALSXP, 1)); pcnt++; REAL(deltat_d)[0] = dt;
	epsilon_d = PROTECT(allocVector(REALSXP, 1)); pcnt++; REAL(epsilon_d)[0] = epsilon;
	part_name_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(part_name_s, 0, mkChar(ptype == FP_DYADIC ? "Dyadic" : "Unknown"));
	win_name_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(win_name_s, 0, mkChar(wtype == FWT_GAUSSIAN ? "Gaussian" : "Unknown"));

	/* put the return values into a list */
	ret_l = PROTECT(allocVector(VECSXP, 13)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, 13)); pcnt++;
	SET_VECTOR_ELT(ret_l, 0, ts_len_i); SET_VECTOR_ELT(names_s, 0, mkChar("ts.len"));
	SET_VECTOR_ELT(ret_l, 1, deltat_d); SET_VECTOR_ELT(names_s, 1, mkChar("dt"));
	SET_VECTOR_ELT(ret_l, 2, epsilon_d); SET_VECTOR_ELT(names_s, 2, mkChar("eps"));
	SET_VECTOR_ELT(ret_l, 3, part_name_s); SET_VECTOR_ELT(names_s, 3, mkChar("part.type"));
	SET_VECTOR_ELT(ret_l, 4, win_name_s); SET_VECTOR_ELT(names_s, 4, mkChar("win.type"));
	SET_VECTOR_ELT(ret_l, 5, gft_c); SET_VECTOR_ELT(names_s, 5, mkChar("gft"));
	SET_VECTOR_ELT(ret_l, 6, cts_img_c); SET_VECTOR_ELT(names_s, 6, mkChar("image"));
	SET_VECTOR_ELT(ret_l, 7, img_wd_i); SET_VECTOR_ELT(names_s, 7, mkChar("wd"));
	SET_VECTOR_ELT(ret_l, 8, img_ht_i); SET_VECTOR_ELT(names_s, 8, mkChar("ht"));
	SET_VECTOR_ELT(ret_l, 9, freq_centers_i); SET_VECTOR_ELT(names_s, 9, mkChar("f.centers"));
	SET_VECTOR_ELT(ret_l, 10, time_centers_i); SET_VECTOR_ELT(names_s, 10, mkChar("t.centers"));
	SET_VECTOR_ELT(ret_l, 11, by_part_l); SET_VECTOR_ELT(names_s, 11, mkChar("by.partition"));
	SET_VECTOR_ELT(ret_l, 12, all_freqs_l); SET_VECTOR_ELT(names_s, 12, mkChar("all.freqs"));
	setAttrib(ret_l, R_NamesSymbol, names_s);

	UNPROTECT(pcnt);

	return ret_l;
}


DllExport SEXP CALLngft_1dComplex64Inv(SEXP gft_c, SEXP ts_len_i, SEXP dt_d,
																			 SEXP eps_d, SEXP ptype_s, SEXP wtype_s)
{
	int ii, pcnt = 0;
	DCMPLX *dst;
	DCLIST *gft, *cts;
	FPCOL *partitions;
	Rcomplex *ctp_r;
	SEXP cts_c;
	SEXP ret_l, names_s;
	SEXP epsilon_d;
	SEXP part_name_s, win_name_s;

	Rcomplex *gft_r = COMPLEX(gft_c);
	int dst_len = length(gft_c);
	int ts_len = asInteger(ts_len_i) == NA_INTEGER ? -1 : INTEGER(ts_len_i)[0];
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	double epsilon = ! R_FINITE(asReal(eps_d)) ? 0 : REAL(dt_d)[0];
	FreqPartitionType ptype = asChar(ptype_s) == NA_STRING ? FP_DYADIC :
		strncasecmp(CHAR(STRING_ELT(ptype_s,0)), "d", 1) == 0 ? FP_DYADIC : /* Dyadic */
		FP_DYADIC; /* silently ignore anything else, and use Dyadic */
	FreqWindowType wtype = asChar(wtype_s) == NA_STRING ? FWT_GAUSSIAN :
		strncasecmp(CHAR(STRING_ELT(wtype_s,0)), "g", 1) == 0 ? FWT_GAUSSIAN : /* Gaussian */
		FWT_GAUSSIAN; /* silently ignore anything else, and use Gaussian */

	if ( dst_len < 3 || ts_len < 3 )
		error("S-transform and original time series must have at least 3 points");
	if ( dt <= 0 )
		error("dt must be positive");

	/* initialize */

	// create partitions
	partitions = ngft_FrequencyPartitions(ts_len, epsilon, ptype, wtype);

	// copy the input dst, and unpack into partitions
	if ( (dst = calloc( dst_len, sizeof( *dst ) )) == NULL )
		oops( "CALLngft_1dComplex64Inv","can't get space" );
	for ( ii = 0 ; ii < dst_len ; ii++ ) {
		dst[ii].r = gft_r[ii].r;
		dst[ii].i = gft_r[ii].i;
	}
	gft = calloc(1, sizeof(*gft));
	gft->count = dst_len;
	gft->values = dst;
	ngft_unpackGftArray(gft, partitions);
	freeDClist(gft); // also frees space pointed to by dst

	// Call 1D GFT Inverse Function
	cts = ngft_1dComplex64Inv(partitions);
	ngft_FreeFreqPartitions(partitions);

	/* allocate space for R structures for complex time series, and copy inv_gft output */
	cts_c = PROTECT(allocVector(CPLXSXP, cts->count)); pcnt++;
	ctp_r = COMPLEX(cts_c);
	for ( ii = 0 ; ii < cts->count ; ii++ ) {
		ctp_r[ii].r = cts->values[ii].r;
		ctp_r[ii].i = cts->values[ii].i;
	}
	freeDClist(cts);

	// other returned items
	epsilon_d = PROTECT(allocVector(REALSXP, 1)); pcnt++; REAL(epsilon_d)[0] = epsilon;
	part_name_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(part_name_s, 0, mkChar(ptype == FP_DYADIC ? "Dyadic" : "Unknown"));
	win_name_s = PROTECT(allocVector(STRSXP, 1)); pcnt++;
	SET_STRING_ELT(win_name_s, 0, mkChar(wtype == FWT_GAUSSIAN ? "Gaussian" : "Unknown"));

	/* put the return values into a list */
	ret_l = PROTECT(allocVector(VECSXP, 4)); pcnt++;
	names_s = PROTECT(allocVector(VECSXP, 4)); pcnt++;
	SET_VECTOR_ELT(ret_l, 0, cts_c); SET_VECTOR_ELT(names_s, 0, mkChar("TS"));
	SET_VECTOR_ELT(ret_l, 1, epsilon_d); SET_VECTOR_ELT(names_s, 1, mkChar("eps"));
	SET_VECTOR_ELT(ret_l, 2, part_name_s); SET_VECTOR_ELT(names_s, 2, mkChar("part.type"));
	SET_VECTOR_ELT(ret_l, 3, win_name_s); SET_VECTOR_ELT(names_s, 3, mkChar("win.type"));
	setAttrib(ret_l, R_NamesSymbol, names_s);

	UNPROTECT(pcnt);

	return ret_l;
}
