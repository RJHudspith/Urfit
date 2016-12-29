#ifndef PMAP_H
#define PMAP_H

struct pmap *
parammap( const size_t Nparams ,
	  const size_t Nsims ,
	  const size_t *Ndata ,
	  const bool *sim_params ) ;

void
free_pmap( struct pmap *map ,
	   const size_t n ) ;

#endif
