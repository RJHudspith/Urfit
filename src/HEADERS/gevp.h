#ifndef GEVP_H
#define GEVP_H

int
solve_gevp( double *re_evalues , 
	    const double *A ,
	    const double *B , 
	    const size_t n ,
	    const bool before ,
	    const bool write_evalues ) ;

#endif
