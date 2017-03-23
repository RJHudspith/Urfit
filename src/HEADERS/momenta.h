#ifndef MOMENTA_H
#define MOMENTA_H

double
lattmom( const size_t *Dimensions ,
	 const size_t Nd ,
	 const int *N ,
	 const int Selection ) ;

int
average_equivalent( struct input_params *Input ) ;

#endif
