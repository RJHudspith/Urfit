/**
   @file stats.c
   @brief statistical sampling
 */
#include "gens.h"

#include "bootstrap.h"
#include "jacknife.h"
#include "raw.h"

#include "rng.h"

// a wrapper for the resampling
void
compute_err( struct resampled *replicas )
{
  switch( replicas -> restype ) {
  case Raw :
    raw_err( replicas ) ;
    break ;
  case JackKnife :
    jackknife_error( replicas ) ;
    break ;
  case BootStrap :
    bootstrap_error( replicas ) ;
    break ;
  }
  return ;
}

// bootstrap or jackknife or whatever
int
resample_data( struct input_params *Input )
{
  size_t i ;
  bool must_resample = false ;
  
  // compute the error and check if we must resample
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    compute_err( &Input -> Data.x[i] ) ;
    compute_err( &Input -> Data.y[i] ) ;
    if( Input -> Data.x[i].restype != Input -> Data.Restype ||
	Input -> Data.x[i].restype != Input -> Data.Restype ) {
      must_resample = true ;
    }
  }
  
  // bootstrap it
  if( must_resample == true ) {
    switch( Input -> Data.Restype ) {
    case BootStrap :
      #ifdef VERBOSE
      fprintf( stdout , "[STATS] Bootstrapping \n" ) ;
      #endif
      init_rng( 123456 ) ;
      bootstrap_full( Input ) ;
      free_rng() ;
      break ;
    case JackKnife :
      #ifdef VERBOSE
      fprintf( stdout , "[STATS] JackKnifing \n" ) ;
      #endif
      jackknife_full( Input ) ;
      break ;
    case Raw :
      #ifdef VERBOSE
      fprintf( stdout , "[STATS] Raw data\n" ) ;
      #endif
      break ;
    }
  }

  #ifdef VERBOSE
  fprintf( stdout , "[STATS] Resampling finished\n" ) ;
  #endif

  return SUCCESS ;
}
