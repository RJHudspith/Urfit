/**
   @file bootfit.c
   @brief perform a bootstrap fit
 */
#include "gens.h"

#include "fit_chooser.h"
#include "resampled_ops.h"

// fit types
#include "GA.h"
#include "CG.h"
#include "SD.h"
#include "LM.h"

// perform a fit over bootstraps
struct resampled *
perform_bootfit( const struct resampled *x ,
		 const struct resampled *y ,
		 const double **W ,
		 const size_t Ndata ,
		 const size_t LT ,
		 const fittype fit ,
		 const corrtype CORRFIT )
{
  if( x[0].NSAMPLES != y[0].NSAMPLES ) {
    printf( "[BOOTFIT] number of x  and y samples different" ) ;
    return NULL ;
  }

  // initialise the fit
  struct fit_descriptor fdesc = init_fit( fit , Ndata , CORRFIT ) ;

  int (*f) ( struct fit_descriptor *fdesc ,
	     const void *data ,
	     const double **W ,
	     const double TOL ) ;
  f = ml_iter ;

  // allocate the fitparams
  struct resampled *fitparams = malloc( fdesc.NPARAMS * sizeof( struct resampled ) ) ;
  size_t i ;
  for( i = 0 ; i < fdesc.NPARAMS ; i++ ) {
    fitparams[i] = init_dist( NULL , y[0].NSAMPLES , y[0].restype ) ;
  }

  // allocate the chisq
  struct resampled chisq = init_dist( NULL , y[0].NSAMPLES , y[0].restype ) ;

  // loop boots
  for( i = 0 ; i < chisq.NSAMPLES ; i++ ) { 

    // initialise the data we will fit
    double yloc[ Ndata ] , xloc[ Ndata ] ;
    size_t j ;
    for( j = 0 ; j < Ndata ; j++ ) {
      xloc[j] = x[j].resampled[i] ;
      yloc[j] = y[j].resampled[i] ;
    }
    struct data d = { Ndata , xloc , yloc , LT , fdesc.NPARAMS } ;

    // do the fit, compute the chisq
    f( &fdesc , &d , (const double**)W , 1E-10 ) ;
    chisq.resampled[i] = fdesc.f.chisq ;

    for( j = 0 ; j < fdesc.NPARAMS ; j++ ) {
      fitparams[j].resampled[i] = fdesc.f.fparams[j] ;
    }
  }

  // multiply by the range
  mult_constant( &chisq , 1.0 / ( Ndata - fdesc.NPARAMS ) ) ;

  // compute and free the chisq
  compute_err( &chisq ) ;

  printf( "[CHISQ] %f %f \n" , chisq.avg , chisq.err ) ;

  free( chisq.resampled ) ;

  for( i = 0 ; i < fdesc.NPARAMS ; i++ ) {
    compute_err( &fitparams[i] ) ;
    printf( "PARAM_%zu %f %f \n" , i , fitparams[i].avg , fitparams[i].err ) ;
  }

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.NPARAMS ) ;

  return fitparams ;
}
