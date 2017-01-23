#include "gft.h"
#include "gft_proto.h"

/* be sure to include these last so they can override some windows defs */
#include <R.h>
#include <Rinternals.h>

DllExport SEXP CALLgft_1dComplex64(SEXP ts_d, SEXP dt_d, SEXP nd_i, SEXP order_i)
{
	int ii, pcnt = 0;
	double *dts = NULL, *dtp;
	SEXP dts_d;

	double *ts = REAL(ts_d);
	int len = length(ts_d);
	double dt = ! R_FINITE(asReal(dt_d)) ? .005 : REAL(dt_d)[0];
	int nd = asInteger(nd_i) == NA_INTEGER ? 1 : INTEGER(nd_i)[0];
	FdOrder order = asInteger(order_i) == NA_INTEGER ? FD_ORDER_8 : INTEGER(order_i)[0];

	if ( len < 3 )
		error("time series must have at least 3 points");
	/* if ( nd < 1 || nd > 2 )
	error("only first and second derivatives supported"); */
	if ( dt <= 0 )
		error("dt must be positive");

	dts = fd_deriv(ts, len, dt, nd, order);

	/* allocate space for R structures for derivative, and copy fd_deriv output */
	dts_d = PROTECT(allocVector(REALSXP, len)); pcnt++;
	dtp = REAL(dts_d);
	for ( ii = 0 ; ii < len ; ii++ )
		dtp[ii] = dts[ii];

	free(dts);

	UNPROTECT(pcnt);

	return dts_d;
}
