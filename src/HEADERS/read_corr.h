#ifndef READ_CORR_H
#define READ_CORR_H

int
get_correlator( double complex *C ,
		const char *str ,
		const size_t snk ,
		const size_t src ,
		const double *mompoint ,
		const size_t Nlt ) ;

int
read_corr( struct input_params *Input ) ;

#endif
