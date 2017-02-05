#ifdef TEST_NGFT

#include "ngft.h"
#include "ngft_proto.h"

static double modulus( DCMPLX *x ) {
	return sqrt( x->r*x->r + x->i*x->i );
}
static void print_freq_partitions(FPCOL *fpcol, FILE *ofile) {
	int ii;
	fprintf( ofile, "\nFrequency Partitions: N=%d\n", fpcol->N );
	// loop over the non-negative (index 0) and the negative (index 1) frequency partition sets
	for ( ii = 0 ; ii < 2 ; ii++ ) {
		int jj;
		FPSET *fpset = fpcol->fpset + ii;
		// loop over the partitions in this set
		fprintf( ofile, "    %s frequencies:\n", ii == 0 ? "non-negative" : "negative" );
		for ( jj = 0 ; jj < fpset->pcount ; jj++ ) {
			int kk;
			FPART *partition = fpset->partitions + jj;
			int fstart = partition->start;
			fprintf( ofile, "\tpartition %2d: start=%3d center=%3d end=%3d width=%3d\n",
							 jj,fstart, partition->center, partition->end, partition->width);
			// apply partition window to the transformed data
			for ( kk = 0 ; kk < partition->win_len ; kk++ ) {
				int ff = fstart + kk;
				DCMPLX *val = partition->window + kk;
				fprintf( ofile, "\t\t%3d %10.3e %10.3e\t%9.3e\n", ff, val->r, val->i, modulus(val));
			}
		}
	}
	fflush( ofile );
}

static void print_time_partitions( TPCOL* tpcol, FILE *ofile ) {
	int ii;
	fprintf( ofile, "\nTime partitions: N=%d, tdcount=%3d\n", tpcol->N, tpcol->tdcount);
	for ( ii = 0 ; ii < tpcol->tdcount ; ii++ ) {
		int jj;
		if ( tpcol->tdsets == NULL )
			continue;
		TDSET *tdset = tpcol->tdsets + ii;
		fprintf( ofile, "\tSet %3d: decimation=%3d, pcount=%2d\n", ii, tdset->decimation, tdset->pcount );
		for ( jj = 0 ; jj < tdset->pcount ; jj++ ) {
			TPART* tpart = tdset->partitions + jj;
			fprintf( ofile, "\t\tPartition %3d: start=%3d, center=%3d, end=%3d, width=%d\n",
							 jj, tpart->start, tpart->center, tpart->end, tpart->width );
		}
	}
	fflush( ofile );
}

int main( int argc, char **argv ) {
	char *progname, *cp;
	char line[2048];
	int ii, line_no, dcount;
	DCMPLX *cts;
	DIMAGE *image;
	FILE *infile = stdin, *outfile = stdout;
	BOOL do_inverse = FALSE, do_image = FALSE, gaussian_window = TRUE;
	BOOL by_part = TRUE, all_freqs = FALSE, ind_map = FALSE;
	int stride = 1, image_dim = -1;
	windowFunction *window_fn = NULL;
	FPCOL *partitions;
	TPCOL *tpcol;

	progname = (cp = strrchr(*argv,PATH_SEP)) != (char *)NULL ? cp+1 : *argv;
	while ( --argc ) {
		if ( **++argv == '-' ) {
			switch ( *(*argv + 1) ) {
				case 'h':
				case 'H':
					fprintf( stdout, "USAGE: %s [-h|-H] [-i] [-I] [-a] [-A] [-m] [-o outfile] [infile]\n", progname );
					break;
				case 'i': // do inverse
					do_inverse = TRUE;
					do_image = FALSE;
					break;
				case 'I':	// output image
					if ( ! do_inverse )
						do_image = TRUE;
					break;
				case 'a': // output image for all times and freqs, instead of by partition
					by_part = FALSE;
					if ( !do_inverse )
						do_image = TRUE;
					break;
				case 'A': // output image results for all frequencies
					all_freqs = TRUE;
					if ( !do_inverse )
						do_image = TRUE;
					break;
				case 'm':	// output image as index map
					ind_map = TRUE;
					if ( !do_inverse )
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
		fprintf( stderr, "Warning: sizeof(DCMPLX) = %d is not equal to sizeof(fftw_complex) = %d",
						 sizeof( DCMPLX ), sizeof( fftw_complex ) );
		fflush( stderr );
	}

	/* initialize */
	window_fn = gaussian_window ? gaussian : box;
	partitions = ngft_DyadicPartitions(dcount);
	ngft_AddWindowsToParts(partitions, window_fn);
	
	// Call 1D GFT Function
	if ( do_inverse ) {
		ngft_1dComplex64Inv( cts, partitions, stride );
		fprintf( outfile, "# INV: Re(z)   Im(z)        |z|\n" );
	} else {
		ngft_1dComplex64( cts, dcount, partitions, stride );
		fprintf( outfile, "# FST: Re(z)   Im(z)        |z|\n" );
	}
	for ( ii = 0 ; ii < dcount ; ii++ )
		fprintf( outfile, "%10.3e %10.3e\t %9.3e\n", cts[ii].r, cts[ii].i, modulus(cts+ii) );

	if ( do_image ) {
		// get time partitions
		tpcol = ngft_TimePartitions( partitions );

		// get complex image
		image = ngft_1d_InterpolateNN(cts, partitions, tpcol, image_dim, by_part, all_freqs, ind_map);

		// print the image or index map
		fprintf( outfile, "\n%s%s:\n", ind_map ? "Index map" : "Image",  all_freqs ? " (all frequencies)" : " (non-negative frequencies)");
		fprintf( outfile, "  f\\t" );
		for ( ii = 0 ; ii < image->x_centers->count ; ii++ )
			fprintf( outfile, " %8d", image->x_centers->values[ii]);
		fprintf( outfile, "\n---- " );
		for ( ii = 0 ; ii < image->x_centers->count ; ii++ )
			fprintf( outfile, " --------");
		fprintf( outfile, "\n" );
		for ( ii = 0 ; ii < image->y_centers->count ; ii++ ) {
			int jj;
			char *fmt = ind_map ? " %8.0f" : " %8.1e";
			fprintf( outfile, "%3d: ", image->y_centers->values[ii]);
			for ( jj = 0 ; jj < image->x_centers->count ; jj++ )
				fprintf( outfile, fmt, modulus( image->img->values + ii * image->x_centers->count + jj ) );
			fprintf( outfile, "\n" );
		}

		// print the frequency and time partitions
		print_freq_partitions(partitions, outfile);
		print_time_partitions(tpcol, outfile);
	}

	free(cts);
	ngft_FreeFreqPartitions( partitions );
	if ( do_image ) {
		freeDImage( image );
		ngft_FreeTimePartitions( tpcol );
	}

	return 0;
}
#endif