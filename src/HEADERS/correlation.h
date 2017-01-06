#ifndef CORRELATION_H
#define CORRELATION_H

int
inverse_correlation( struct data_info *Data ,
		     const struct fit_info Fit ) ;

void
correlations( double **correlation , 
	      const struct resampled *data ,
	      const corrtype CORRFIT ,
	      const size_t Ndata ,
	      const bool Divided ) ;

void
write_corrmatrix( const double **correlation ,
		  const size_t NCUT ) ;

void
write_corrmatrix_mathematica( const double **correlation ,
			      const size_t NCUT ) ;

void
write_corrmatrix_to_file( FILE *outfile , 
			  const double **correlation ,
			  const size_t NCUT ) ;

#endif
