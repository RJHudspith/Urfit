#ifndef CORRELATION_H
#define CORRELATION_H

int
inverse_correlation( struct data_info *Data ,
		     const struct fit_info Fit ) ;

void
write_corrmatrix( const double **correlation ,
		  const size_t NCUT ,
		  const corrtype Corrfit ) ;

#endif
