#ifndef GEVP_H
#define GEVP_H

struct resampled *
solve_GEVP( const struct resampled *y ,
	    const size_t Ndata ,
	    const size_t N ,
	    const size_t M ,
	    const size_t t0 ,
	    const size_t td ) ;

#endif
