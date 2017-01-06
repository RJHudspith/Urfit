#ifndef PMAP_H
#define PMAP_H

struct pmap *
parammap( const struct data_info Data ,
	  const struct fit_info Fit ) ;

void
free_pmap( struct pmap *map ,
	   const size_t n ) ;

#endif
