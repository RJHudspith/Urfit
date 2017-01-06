/**
   @file bootfit.c
   @brief perform a bootstrap fit
 */
#include "gens.h"

#include "ffunction.h"
#include "fit_chooser.h"
#include "resampled_ops.h"
#include "stats.h"

// perform a fit over bootstraps
struct resampled *
perform_bootfit( const struct data_info Data ,
		 const struct fit_info Fit )
{
  if( Data.x[0].NSAMPLES != Data.y[0].NSAMPLES ) {
    printf( "[BOOTFIT] number of x  and y samples different" ) ;
    return NULL ;
  }

  // bootstrap counter
  size_t i ;

  // initialise the fit
  struct fit_descriptor fdesc = init_fit( Data , Fit ) ;
  
  // allocate the fitparams
  struct resampled *fitparams = malloc( fdesc.Nlogic * sizeof( struct resampled ) ) ; 
  for( i = 0 ; i < fdesc.Nlogic ; i++ ) {
    fitparams[i] = init_dist( NULL , Data.y[0].NSAMPLES , Data.y[0].restype ) ;
  }

  // allocate the chisq
  struct resampled chisq = init_dist( NULL , Data.y[0].NSAMPLES , Data.y[0].restype ) ;
  
  // loop boots
  for( i = 0 ; i < chisq.NSAMPLES ; i++ ) { 

    // initialise the data we will fit
    double *yloc = malloc( Data.Ntot * sizeof( double ) ) ;
    double *xloc = malloc( Data.Ntot * sizeof( double ) ) ;
    size_t j ;
    for( j = 0 ; j < Data.Ntot ; j++ ) {
      xloc[j] = Data.x[j].resampled[i] ;
      yloc[j] = Data.y[j].resampled[i] ;
    }
    struct data d = { Data.Ntot , xloc , yloc , Data.LT ,
		      fdesc.Nparam , Fit.map } ;

    // do the fit, compute the chisq
    Fit.Minimize( &fdesc , &d , (const double**)Data.Cov.W , Fit.Tol ) ;
    chisq.resampled[i] = fdesc.f.chisq ;

    for( j = 0 ; j < fdesc.Nlogic ; j++ ) {
      fitparams[j].resampled[i] = fdesc.f.fparams[j] ;
    }
    free( yloc ) ;
    free( xloc ) ;
  }

  // divide out the number of degrees of freedom
  mult_constant( &chisq , 1.0 / ( Data.Ntot - fdesc.Nlogic ) ) ;

  // compute and free the chisq
  compute_err( &chisq ) ;

  printf( "[CHISQ (d.o.f)] %f %f \n" , chisq.avg , chisq.err ) ;

  free( chisq.resampled ) ;

  // tell us what we have computed
  for( i = 0 ; i < fdesc.Nlogic ; i++ ) {
    compute_err( &fitparams[i] ) ;
    if( Fit.Sims[i] == true ) {
      printf( "-> SIMUL " ) ;
    }
    printf( "PARAM_%zu %f %f \n" , i , fitparams[i].avg , fitparams[i].err ) ;
  }

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  return fitparams ;
}
