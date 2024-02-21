#ifndef CORRELATION_H
#define CORRELATION_H

void
compute_upper_correlation( double **correlation , 
			   const struct resampled *data ,
			   const size_t NDATA ,
			   const corrtype CORRFIT ) ;

void
fill_lower_triangular( double **correlation ,
		       const size_t NDATA ) ;

void
modified_covariance( double **correlation ,
		     const size_t Ndata ) ;

int
inverse_correlation( struct data_info *Data ,
		     const struct fit_info Fit ) ;

void
write_corrmatrix( const double **correlation ,
		  const size_t NCUT ,
		  const corrtype Corrfit ) ;

#endif
