#ifndef BLACKBOX_H
#define BLACKBOX_H

void
blackbox( const double *data ,
	  const size_t NDATA ,
	  const size_t NSTATES ,
	  double masses[ NSTATES ][ NDATA ] ) ;

#endif
