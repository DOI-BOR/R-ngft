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
			FPART *partition = fpset->partitions + jj;
			int fstart = partition->start;
			fprintf( ofile, "\tpart %2d: start=%3d center=%3d end=%3d width=%3d   ",
							 jj,fstart, partition->center, partition->end, partition->width);
			fprintf( ofile, "win_start=%3d win_end=%3d win_len=%3d\n",
							partition->win_start, partition->win_end, partition->win_len);
		}
	}
	fflush( ofile );
}

static void print_time_partitions( FPCOL* fpcol, FILE *ofile ) {
	int ii;
	TPCOL* tpcol = fpcol->tpcol;

	fprintf( ofile, "\nTime partitions: N=%d, tdcount=%3d\n", tpcol->N, tpcol->tdcount);
	for ( ii = 0 ; ii < tpcol->tdcount ; ii++ ) {
		int jj;
		if ( tpcol->tdsets == NULL )
			continue;
		TDSET *tdset = tpcol->tdsets + ii;
		fprintf( ofile, "\tSet %3d: decimation=%.3e, pcount=%2d\n", ii, tdset->decimation, tdset->pcount );
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
	int ii, line_no, dcount, ts_len = -1;
	DCMPLX *data_in;
	FILE *infile = stdin, *outfile = stdout;
	BOOL do_inverse = FALSE, do_image = FALSE, gaussian_window = TRUE;
	BOOL by_part = TRUE, all_freqs = FALSE;
	int stride = 1, image_dim = -1, fw_width = -1, edo_fref = -1, edo_nd = -1;
	FPCOL *partitions;
	FreqPartitionType ptype = FP_DYADIC;
	FreqWindowType wtype = FWT_GAUSSIAN;
	double epsilon = -1;

	progname = (cp = strrchr(*argv,PATH_SEP)) != (char *)NULL ? cp+1 : *argv;
	while ( --argc ) {
		if ( **++argv == '-' ) {
			switch ( *(*argv + 1) ) {
				case 'h':
				case 'H':
					fprintf( stdout, "USAGE: %s [-h|-H] [-i] [-I] [-a] [-A] [-M dim] [-e eps]\n\t\t[-p ptype] [-w wtype] [-W fw_width] [-f edo_fref] [-n edo_nd]\n\t\t[-o outfile] [infile]\n", progname );
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
					break;
				case 'A': // output image results for all frequencies
					all_freqs = TRUE;
					break;
				case 'e': // epsilon
					epsilon = atof(*++argv);
					argc--;
					break;
				case 'M': // image dimension
					image_dim = atoi(*++argv);
					argc--;
					break;
				case 'p':	// partition type: 1 - dyadic, 2 - EDO, 3 - fixed width
					ptype = atoi(*++argv);
					argc--;
					break;
				case 'w':	// window type: 1 - Gaussian, 2 - Box
					wtype = atoi(*++argv);
					argc--;
					break;
				case 'W': // width for fixed width partitions
					fw_width = atoi(*++argv);
					argc--;
					break;
				case 'f': // EDO reference frequency
					edo_fref = atoi(*++argv);
					argc--;
					break;
				case 'n': // EDO number of divisions per octave
					edo_nd = atoi(*++argv);
					argc--;
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
	data_in = NULL;
	while ( fgets( line, sizeof( line ), infile ) != NULL ) {
		line_no++;
		/* strip the newline and any other whitespace */
		while ( strlen( line ) > 0 && isspace( line[strlen( line ) - 1] ) )
			line[strlen( line ) - 1] = '\0';
		/* skip comment lines */
		if ( *line == '#' || *line == '*' || strlen( line ) < 1 )
			continue;
		if ( do_inverse ) {
			if ( ts_len <= 0 ) {
				// first line is ts_len
				sscanf(line, "%d", &ts_len);
			} else {
				// subsequent lines are complex dst data
				data_in = realloc( data_in, ++dcount * sizeof( *data_in ) );
				sscanf(line, "%lf %lf", &data_in[dcount - 1].r, &data_in[dcount - 1].i);
			}
		} else {
			data_in = realloc( data_in, ++dcount * sizeof( *data_in ) );
			sscanf( line, "%lf", &data_in[dcount - 1].r);
			data_in[dcount - 1].i = 0;
		}
	}
	fprintf( stderr, "read %d lines, %d data points\n", line_no, dcount );

	if ( dcount < 4 )
		return 1;

	if ( sizeof( DCMPLX ) != sizeof( fftw_complex ) ) {
		fprintf( stderr, "Warning: sizeof(DCMPLX) = %zd is not equal to sizeof(fftw_complex) = %zd",
						 sizeof( DCMPLX ), sizeof( fftw_complex ) );
		fflush( stderr );
	}

	/* initialize */

	// Call 1D GFT Function
	if ( ! do_inverse ) {
		DCLIST *gft;
		DCLIST *signal = calloc(1, sizeof(*signal));
		signal->count = dcount;
		signal->values = data_in;
		partitions = ngft_1dComplex64(signal, epsilon, ptype, wtype, fw_width, edo_fref, edo_nd);
		freeDClist(signal); // also frees space pointed to by data_in
		gft = ngft_makeGftArray(partitions);
		fprintf( outfile, "%d\n", partitions->N );
		fprintf( outfile, "# FST: Re(z)         Im(z)            |z|\n" );
		for ( ii = 0 ; ii < gft->count ; ii++ )
			fprintf( outfile, "%15.8e %15.8e\t %14.8e\n",
							gft->values[ii].r, gft->values[ii].i, modulus(gft->values+ii) );
	} else {
		DCLIST *gft, *cts;
		// create partitions
		partitions = ngft_FrequencyPartitions(ts_len, epsilon, ptype, wtype, fw_width, edo_fref, edo_nd);
		gft = calloc(1, sizeof(*gft));
		gft->count = dcount;
		gft->values = data_in;
		ngft_unpackGftArray(gft, partitions);
		freeDClist(gft); // also frees space pointed to by data_in
		cts = ngft_1dComplex64Inv(partitions);
		fprintf( outfile, "# INV: Re(z)         Im(z)            |z|\n" );
		for ( ii = 0 ; ii < cts->count ; ii++ )
			fprintf( outfile, "%15.8e %15.8e\t %14.8e\n",
							cts->values[ii].r, cts->values[ii].i, modulus(cts->values+ii) );
		freeDClist(cts);
	}

	if ( do_image ) {
		// get complex image
		DIMAGE *image = ngft_1d_InterpolateNN(partitions, image_dim, by_part, all_freqs);

		// print the image
		fprintf( outfile, "\n%s%s:\n", "Image",  all_freqs ? "   (all frequencies)" : " (non-negative frequencies)");
		fprintf( outfile, "  f\\t" );
		for ( ii = 0 ; ii < image->x_centers->count ; ii++ )
			fprintf( outfile, " %8d", image->x_centers->values[ii]);
		fprintf( outfile, "\n------ " );
		for ( ii = 0 ; ii < image->x_centers->count ; ii++ )
			fprintf( outfile, " --------");
		fprintf( outfile, "\n" );
		for ( ii = 0 ; ii < image->y_centers->count ; ii++ ) {
			int jj;
			char *fmt = " %8.1e";
			fprintf( outfile, "%5d: ", image->y_centers->values[ii]);
			for ( jj = 0 ; jj < image->x_centers->count ; jj++ )
				fprintf( outfile, fmt, modulus( image->img->values + ii * image->x_centers->count + jj ) );
			fprintf( outfile, "\n" );
		}

		// print the frequency and time partitions
		print_freq_partitions(partitions, outfile);
		print_time_partitions(partitions, outfile);

		freeDImage( image );
	}

	ngft_FreeFreqPartitions( partitions );

	return 0;
}
#endif