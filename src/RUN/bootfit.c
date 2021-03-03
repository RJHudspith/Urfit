/**
   @file bootfit.c
   @brief perform a bootstrap fit

   Uses the fit result of the average values as a starting guess and 
   does the remaining bootstrap samples in parallel, threaded with 
   openmp
 */
#include "gens.h"

#include "ffunction.h"
#include "fit_chooser.h"
#include "resampled_ops.h"
#include "stats.h"

#include <gsl/gsl_cdf.h> // pvalue

// perform a single bootstrap fit to our data
static int
single_fit( struct resampled *fitparams ,
	    struct resampled *chisq ,
	    struct fit_descriptor fdesc ,
	    const struct data_info Data ,
	    const struct fit_info Fit ,
	    const size_t sample_idx ,
	    const bool is_average )
{
  if( Data.Ntot == 0 || fdesc.Nlogic == 0 ) return FAILURE ;
  
  // initialise the data we will fit
  double *yloc = malloc( Data.Ntot * sizeof( double ) ) ;
  double *xloc = malloc( Data.Ntot * sizeof( double ) ) ;
  int Flag = SUCCESS ;
  
  // initialise the data we are fitting
  size_t j ;
  for( j = 0 ; j < Data.Ntot ; j++ ) {
    if( is_average == true ) {
      xloc[j] = Data.x[j].avg ;
      yloc[j] = Data.y[j].avg ;
    } else {
      xloc[j] = Data.x[j].resampled[sample_idx] ;
      yloc[j] = Data.y[j].resampled[sample_idx] ;
    }
  }
  
  struct data d = { Data.Ntot , xloc , yloc , Data.LT ,
		    fdesc.Nparam , Fit.map , Fit.N , Fit.M } ;

  // set the data to the fit params average for a guess
  // guesses are either generated in the fit function or by
  // the user in the input file
  if( is_average == true ) {
    if( Fit.Guesses_Initialised == false ) {
      fdesc.guesses( fdesc.f.fparams , Data , Fit ) ;
    } else {
      for( j = 0 ; j < fdesc.Nlogic ; j++ ) {
	fdesc.f.fparams[j] = Fit.Guess[j] ;
      }
    }
  } else {
    for( j = 0 ; j < fdesc.Nlogic ; j++ ) {
      fdesc.f.fparams[j] = fitparams[j].avg ;
    }
  }
  
  // do the fit, compute the chisq
  if( Fit.Minimize( &fdesc , &d , (const double**)Data.Cov.W ,
		    Fit.Tol ) == FAILURE ) {
    Flag = FAILURE ;
  }

  // set the chisq
  if( is_average == true ) {
    chisq -> avg = fdesc.f.chisq ;
  } else {
    chisq -> resampled[sample_idx] = fdesc.f.chisq ;
  }

  // set the fit parameters
  for( j = 0 ; j < fdesc.Nlogic ; j++ ) {
    if( is_average == true ) {
      fitparams[j].avg = fdesc.f.fparams[j] ;
    } else {
      fitparams[j].resampled[sample_idx] = fdesc.f.fparams[j] ;
    }
  }

  free( yloc ) ;
  free( xloc ) ;
  
  return Flag ;
}	    

// perform a fit over bootstraps
struct resampled *
perform_bootfit( const struct data_info Data ,
		 const struct fit_info Fit ,
		 double *Chi )
{
  // TODO :: sanity check all indices?
  if( Data.x[0].NSAMPLES != Data.y[0].NSAMPLES ) {
    fprintf( stderr , "[BOOTFIT] number of x  and y samples different" ) ;
    return NULL ;
  }

  // bootstrap counter
  size_t i ;

  // initialise the fit
  struct fit_descriptor fdesc = init_fit( Data , Fit ) ;
  fdesc.Prior = (const struct prior*)Fit.Prior ;

  fprintf( stdout , "[FIT] check Nlogic %zu\n" , fdesc.Nlogic ) ;
  
  // allocate the fitparams
  struct resampled *fitparams = malloc( fdesc.Nlogic * sizeof( struct resampled ) ) ; 
  for( i = 0 ; i < fdesc.Nlogic ; i++ ) {
    fitparams[i] = init_dist( NULL , Data.y[0].NSAMPLES ,
			      Data.y[0].restype ) ;
    fitparams[i].avg = UNINIT_FLAG ;
  }

  // allocate the chisq
  struct resampled chisq = init_dist( NULL , Data.y[0].NSAMPLES ,
				      Data.y[0].restype ) ;

  if( fitparams == NULL ) goto memfree ;

  fprintf( stdout , "[FIT] single fit for the average\n" ) ;
  
  // do the average first
  single_fit( fitparams , &chisq , fdesc , Data , Fit , 0 , true ) ;

  // do the other boots in parallel
  #pragma omp parallel
  {
    // allocate another fit descriptor for the loop over boots
    struct fit_descriptor fdesc_boot = init_fit( Data , Fit ) ;
    fdesc_boot.Prior = Fit.Prior ;
    
    // loop boots
    #pragma omp for private(i)
    for( i = 0 ; i < chisq.NSAMPLES ; i++ ) { 
      single_fit( fitparams , &chisq , fdesc_boot ,
		  Data , Fit , i , false ) ;
    }

    // free the fitfunction
    free_ffunction( &fdesc_boot.f , fdesc.Nlogic ) ;
  }
  
  // divide out the number of degrees of freedom
  const size_t Dof = ( Data.Ntot - fdesc.Nlogic + Fit.Nprior ) ;
  if( Dof != 0 ) {
    divide_constant( &chisq , Dof ) ;
    fprintf( stdout , "\n[ CHISQ/dof ] %e %e (Ndof) %zu\n" ,
	     chisq.avg , chisq.err , Dof ) ;
  }
  
  // set the chi value
  *Chi = chisq.avg ;
  
  // tell us what we have computed
  for( i = 0 ; i < fdesc.Nlogic ; i++ ) {
    compute_err( &fitparams[i] ) ;
    if( i < fdesc.Nparam ) {
      if( Fit.Sims[i] == true ) {
	fprintf( stdout , "-> SIMUL " ) ;
      }
    }
    fprintf( stdout , "[FIT] PARAM_%zu %e %e \n" ,
	     i , fitparams[i].avg , fitparams[i].err ) ;
  }

  // following we could compute a p-value
  // http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node4.html  
  fprintf( stdout , "[FIT] pvalue %f\n" , 1-gsl_cdf_chisq_P( Dof*(*Chi) , Dof ) ) ;
  
 memfree :

  // free the chisq
  if( chisq.resampled != NULL ) {
    free( chisq.resampled ) ;
  }
  
  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  return fitparams ;
}
