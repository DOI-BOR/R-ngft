#include "gft.h"
#include "gft_proto.h"

static double modulus( DCMPLX *x ) {
	return sqrt( x->r*x->r + x->i*x->i );
}
static void print_freq_partitions(FPCOL *fpcol) {
	int ii;
	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the partitions in this set
		fprintf( stderr, "%s frequencies:\n", ii == 0 ? "non-negative" : "negative" );
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			FPART *partition = fpset->partitions + jj;
			int fstart = partition->start;
			fprintf( stderr, "\tpartition %2d: start=%3d center=%3d end=%3d width=%3d\n",
							 jj,fstart, partition->center, partition->end, partition->width);
			// apply partition window to the transformed data 
			for ( kk = 0 ; kk < partition->win_len ; kk++ ) {
				int ff = fstart + kk;
				DCMPLX *val = partition->window + kk;
				fprintf( stderr, "\t\t%3d %10.3e %10.3e\t%9.3e\n", ff, val->r, val->i, modulus(val));
			}
		}
	}
	fflush( stderr );
}

static void print_time_partitions( TPCOL* tpcol ) {
	int ii;
	fprintf( stderr, "Time partitions: N=%4d, tdcount=%3d\n", tpcol->N, tpcol->tdcount);
	for ( ii = 0 ; ii < tpcol->tdcount ; ii++ ) {
		int jj;
		if ( tpcol->tdsets == NULL )
			continue;
		TDSET *tdset = tpcol->tdsets + ii;
		fprintf( stderr, "\tSet %3d: decimation=%3d, pcount=%2d\n", ii, tdset->decimation, tdset->pcount );
		for ( jj = 0 ; jj < tdset->pcount ; jj++ ) {
			TPART* tpart = tdset->partitions + jj;
			fprintf( stderr, "\t\tPartition %3d: start=%3d, center=%3d, end=%3d, width=%d\n",
							 jj, tpart->start, tpart->center, tpart->end, tpart->width );
		}
	}

}

int main( int argc, char **argv ) {
	char *progname, *cp;
	char line[2048];
	int ii, line_no, dcount;
	double dt = 0.01;
	DCMPLX *cts;
	DIMAGE *image;
	FILE *infile = stdin, *outfile = stdout;
	BOOL do_inverse = FALSE;
	BOOL do_image = FALSE;
	BOOL gaussian_window = TRUE;
	int stride = 1, image_dim = -1;
	windowFunction *window_fn = NULL;
	FPCOL *partitions;
	TPCOL *tpcol;
	ILIST *f_centers;
	ILIST *t_centers;

	progname = (cp = strrchr(*argv,PATH_SEP)) != (char *)NULL ? cp+1 : *argv;
	while ( --argc ) {
		if ( **++argv == '-' ) {
			switch ( *(*argv + 1) ) {
				case 'h':
				case 'H':
					fprintf( stdout, "USAGE: %s [-h|-H] [-i] [-I] [-o outfile] [infile]\n", progname );
					break;
				case 'i': // do inverse
					do_inverse = TRUE;
					do_image = FALSE;
					break;
				case 'I':	// output image
					if ( ! do_inverse )
						do_image = TRUE;
					break;
				case 'o':
					outfile = fopen( *++argv, "w" );
					argc--;
					break;
			}
		} else {
			infile = fopen( *argv, "r" );
		}
	}
	if ( infile == NULL )
		return 1;

	line_no = 0;
	dcount = 0;
	cts = NULL;
	while ( fgets( line, sizeof( line ), infile ) != NULL ) {
		line_no++;
		/* skip comment lines */
		if ( *line == '#' || *line == '*' || strlen( line ) < 1 )
			continue;
		/* strip the newline */
		line[strlen( line ) - 1] = '\0';
		cts = realloc( cts, ++dcount * sizeof( *cts ) );
		if ( do_inverse ) {
			sscanf( line, "%lf %lf", &cts[dcount - 1].r, &cts[dcount - 1].i );
		} else {
			sscanf( line, "%lf", &cts[dcount - 1].r);
			cts[dcount - 1].i = 0;
		}
	}
	fprintf( stderr, "read %d lines, %d data points\n", line_no, dcount );

	if ( dcount < 4 )
		return 1;

	if ( sizeof( DCMPLX ) != sizeof( fftw_complex ) ) {
		fprintf( stderr, "Warning: sizeof(DCMPLX) = % is not equal to sizeof(fftw_complex) = %d",
						 sizeof( DCMPLX ), sizeof( fftw_complex ) );
		fflush( stderr );
	}

	/* initialize */
	window_fn = gaussian_window ? gaussian : box;
	partitions = ngft_DyadicPartitions(dcount);
	ngft_AddWindowsToParts(partitions, window_fn);
	print_freq_partitions(partitions);

	// Call 1D GFT Function
	if ( do_inverse ) {
		ngft_1dComplex64Inv( cts, partitions, stride );
		fprintf( outfile, "# inv-FST\n" );
	} else {
		ngft_1dComplex64( cts, dcount, partitions, stride );
		fprintf( outfile, "# FST\n" );
	}
	for ( ii = 0 ; ii < dcount ; ii++ )
		fprintf( outfile, "%10.3e %10.3e\t %9.3e\n", cts[ii].r, cts[ii].i, modulus(cts+ii) );

	if ( do_image ) {
		BOOL do_log = TRUE, make_ind_map = FALSE;
		// get complex image
		tpcol = ngft_TimePartitions( partitions );
		print_time_partitions(tpcol);

		image_dim = -1;
		if ( do_log )
			image = ngft_1d_logfInterpolateNN(cts, partitions, tpcol, image_dim, make_ind_map);
		else
			image = ngft_1d_interpolateNN(cts, partitions, tpcol, image_dim, make_ind_map);

		// get time and frequency centers
		f_centers = getFreqCenters( partitions );
		t_centers = getTimeCenters( tpcol );

		fprintf( stderr, "\n%s:\n", make_ind_map ? "Index map" : "Image" );
		for ( ii = 0 ; ii < image->ht ; ii++ ) {
			int jj;
			fprintf( outfile, "%3d: ", do_log ? f_centers->values[ii] : image_dim > 0 ? ROUND(ii * dcount / (double)image_dim) : ii);
			for ( jj = 0 ; jj < image->wd ; jj++ )
				fprintf( outfile, " %8.1e", modulus( image->img->values + ii*image->wd + jj ) );
			fprintf( outfile, "\n" );
		}
		fprintf( outfile, "     " );
		for ( ii = 0 ; ii < image->wd ; ii++ )
			fprintf( outfile, " %8d", do_log ? t_centers->values[ii] : image_dim > 0 ? ROUND(ii * dcount / (double)image_dim) : ii);
		fprintf( outfile, "\n" );
	}

	free(cts);
	ngft_FreeFreqPartitions( partitions );
	if ( do_image ) {
		freeDImage( image );
		ngft_FreeTimePartitions( tpcol );
		freeIlist( f_centers );
		freeIlist( t_centers );
	}

	return 0;
}