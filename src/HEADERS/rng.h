#ifndef RNG_H
#define RNG_H

void
free_rng( void ) ;

void
init_rng( size_t Seed ) ;

double
rng_double( void ) ;

double
rng_gaussian( const double sigma ) ;

size_t
rng_int( const size_t max_idx ) ;

void
rng_reseed( void ) ;

#endif
