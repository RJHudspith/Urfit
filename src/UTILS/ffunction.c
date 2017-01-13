#include "gens.h"

#include <string.h>

// allocate the fit function
struct ffunction
allocate_ffunction( const size_t NPARAMS ,
		    const size_t NDATA )
{
  struct ffunction f ;
  f.NPARAMS = NPARAMS ;
  f.N = NDATA ;
  f.fparams = malloc( NPARAMS * sizeof( double ) ) ;
  f.f = malloc( NDATA * sizeof( double ) ) ;
  f.df = malloc( NPARAMS * sizeof( double* ) ) ;
  f.d2f = malloc( NPARAMS * NPARAMS * sizeof( double* ) ) ;
  size_t i ;
  for( i = 0 ; i < NPARAMS ; i++ ) {
    f.df[i] = malloc( NDATA * sizeof( double ) ) ;
  }
  for( i = 0 ; i < NPARAMS*NPARAMS ; i++ ) {
    f.d2f[i] = malloc( NDATA * sizeof( double ) ) ;
  }
  // initialise priors
  for( i = 0 ; i < NPARAMS ; i++ ) {
    f.fparams[i] = UNINIT_FLAG ;
  }
  return f ;
}

#define mcpy

// copy the ffunction
void
copy_ffunction( struct ffunction *f1 ,
		const struct ffunction f )
{
  size_t i , j ;
  f1 -> N = f.N ;
  f1 -> NPARAMS = f.NPARAMS ;
  f1 -> CORRFIT = f.CORRFIT ;
  f1 -> Prior = f.Prior ;
  // copy the data
  for( i = 0 ; i < f.N ; i++ ) {
    f1 -> f[i] = f.f[i] ;
  }
  // first deriv
  for( i = 0 ; i < f.NPARAMS ; i++ ) {
    for( j = 0 ; j < f.N ; j++ ) {
      f1 -> df[i][j] = f.df[i][j] ;
    }
    f1 -> fparams[i]   = f.fparams[i] ;
  }
  // second deriv
  for( i = 0 ; i < f.NPARAMS*f.NPARAMS ; i++ ) {
    for( j = 0 ; j < f.N ; j++ ) {
      f1 -> d2f[i][j] = f.d2f[i][j] ;
    }
  }
  return ;
}

// free the fit function
void
free_ffunction( struct ffunction *f , 
		const size_t NPARAMS )
{
  size_t i ;
  for( i = 0 ; i < NPARAMS ; i++ ) {
    free( f -> df[i] ) ;
  }
  for( i = 0 ; i < NPARAMS*NPARAMS ; i++ ) {
    free( f -> d2f[i] ) ;
  }
  free( f -> fparams ) ;
  free( f -> f ) ;
  free( f -> df ) ;
  free( f -> d2f ) ;
  return ;
}
