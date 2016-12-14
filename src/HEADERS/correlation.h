#ifndef CORRELATION_H
#define CORRELATION_H

void
correlations( double **correlation , 
	      const struct resampled *data ,
	      const corrtype CORRFIT ,
	      const size_t NDATA ) ;

double **
correlations_inv( const struct resampled *data ,
		  const corrtype CORRFIT ,
		  const size_t NDATA ) ;

void
write_corrmatrix( double **correlation ,
		  const size_t NCUT ) ;

void
write_corrmatrix_mathematica( double **correlation ,
			      const size_t NCUT ) ;

void
write_corrmatrix_to_file( FILE *outfile , 
			  const double **correlation ,
			  const size_t NCUT ) ;

#endif
