#ifndef READ_CORR_H
#define READ_CORR_H

int
get_correlator( double complex *C ,
		const char *str ,
		const size_t snk ,
		const size_t src ,
		const uint32_t *mompoint ,
		const size_t Nlt ) ;

int
read_corr( struct input_params *Input ,
	   const fold_type fold ,
	   const size_t src ,
	   const size_t snk ) ;

#endif
