#ifndef LM_H
#define LM_H

int
ml_iter( struct fit_descriptor *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL ) ;

#endif
