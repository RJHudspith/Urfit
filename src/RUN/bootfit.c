/**
   @file bootfit.c
   @brief perform a bootstrap fit
 */
#include "gens.h"

#include "ffunction.h"
#include "fit_chooser.h"
#include "pmap.h"
#include "resampled_ops.h"
#include "stats.h"

// fit types
#include "GA.h"
#include "CG.h"
#include "SD.h"
#include "LM.h"

// perform a fit over bootstraps
struct resampled *
perform_bootfit( size_t *Npars ,
		 const struct resampled *x ,
		 const struct resampled *y ,
		 const double **W ,
		 const size_t *Ndata ,
		 const size_t Nsims ,
		 const bool *sims , // simultaneous fit indices
		 const size_t LT ,
		 const fittype fit ,
		 const corrtype CORRFIT )
{
  if( x[0].NSAMPLES != y[0].NSAMPLES ) {
    printf( "[BOOTFIT] number of x  and y samples different" ) ;
    return NULL ;
  }

  // compute the total amount of data there is
  size_t i , N = 0 ;
  for( i = 0 ; i < Nsims ; i++ ) {
    N += Ndata[i] ;
  }

  // initialise the fit
  struct fit_descriptor fdesc = init_fit( fit , N , CORRFIT ,
					  Nsims , sims ) ;
  *Npars = fdesc.Nlogic ;

  // initialise the type of fit
  int (*f) ( struct fit_descriptor *fdesc ,
	     const void *data ,
	     const double **W ,
	     const double TOL ) ;
  f = ml_iter ;

  // set up the param map
  struct pmap *map = parammap( fdesc.Nparam , Nsims , Ndata , sims ) ;
  
  // allocate the fitparams
  struct resampled *fitparams = malloc( fdesc.Nlogic * sizeof( struct resampled ) ) ; 
  for( i = 0 ; i < fdesc.Nlogic ; i++ ) {
    fitparams[i] = init_dist( NULL , y[0].NSAMPLES , y[0].restype ) ;
  }

  // allocate the chisq
  struct resampled chisq = init_dist( NULL , y[0].NSAMPLES , y[0].restype ) ;
  
  // loop boots
  for( i = 0 ; i < chisq.NSAMPLES ; i++ ) { 

    // initialise the data we will fit
    double *yloc = malloc( N * sizeof( double ) ) ;
    double *xloc = malloc( N * sizeof( double ) ) ;
    size_t j ;
    for( j = 0 ; j < N ; j++ ) {
      xloc[j] = x[j].resampled[i] ;
      yloc[j] = y[j].resampled[i] ;
    }
    struct data d = { N , xloc , yloc , LT , fdesc.Nparam , map } ;

    // do the fit, compute the chisq
    f( &fdesc , &d , (const double**)W , 1E-10 ) ;
    chisq.resampled[i] = fdesc.f.chisq ;

    for( j = 0 ; j < fdesc.Nlogic ; j++ ) {
      fitparams[j].resampled[i] = fdesc.f.fparams[j] ;
    }
    free( yloc ) ;
    free( xloc ) ;
  }

  // multiply by the range
  mult_constant( &chisq , 1.0 / ( N - fdesc.Nlogic ) ) ;

  // compute and free the chisq
  compute_err( &chisq ) ;

  printf( "[CHISQ] %f %f \n" , chisq.avg , chisq.err ) ;

  free( chisq.resampled ) ;

  // tell us what we have computed
  for( i = 0 ; i < fdesc.Nlogic ; i++ ) {
    compute_err( &fitparams[i] ) ;
    if( sims[i] == true ) {
      printf( "-> SIMUL " ) ;
    }
    printf( "PARAM_%zu %f %f \n" , i , fitparams[i].avg , fitparams[i].err ) ;
  }

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  // free the parameter map
  free_pmap( map , N ) ;

  return fitparams ;
}
