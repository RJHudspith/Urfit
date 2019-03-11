/**
   @file io_wrapper.c
   @brief wraps the IO
 */
#include "gens.h"

#include "distribution.h"
#include "fake.h"
#include "init.h"
#include "read_corr.h"
#include "read_flat.h"
#include "read_GLU.h"
#include "read_GLU_tcorr.h"
#include "read_GLU_Qmoment.h"

// wrapper function for the IO
int
io_wrap( struct input_params *Input )
{
  // IO switch
  switch( Input -> FileType ) {
  case Corr_File :
    if( read_corr( Input ) == FAILURE ) {
      return FAILURE ;
    }
    // set Lt
    if( init_LT( &Input -> Data , Input -> Traj ) == FAILURE ) {
      return FAILURE ;
    }
    return SUCCESS ;
  case Distribution_File :
    if( read_distribution_old( Input ) == FAILURE ) {
      fprintf( stderr , "[IO] Dist reading failed\n" ) ;
      return FAILURE ;
    }
    return SUCCESS ;
  case Fake_File :
    if( generate_fake_data( &Input -> Data , Input -> Fit ,
			    Input -> Traj , 0.0 , 0.005 ) == FAILURE ) {
      return FAILURE ;
    }
    return SUCCESS ;
  case Flat_File :
    if( read_flat( Input ) == FAILURE ) {
      return FAILURE ;
    }
    // set Lt
    if( init_LT( &Input -> Data , Input -> Traj ) == FAILURE ) {
      return FAILURE ;
    }
    return SUCCESS ;
  case GLU_File :
    if( read_GLU( Input ) == FAILURE ) {
      return FAILURE ;
    }
    // set Lt
    if( init_LT( &Input -> Data , Input -> Traj ) == FAILURE ) {
      return FAILURE ;
    }
    return SUCCESS ;
  case GLU_Tcorr_File :
    if( read_GLU_tcorr( Input ) == FAILURE ) {
      return FAILURE ;
    }
    // set Lt
    if( init_LT( &Input -> Data , Input -> Traj ) == FAILURE ) {
      return FAILURE ;
    }
    return SUCCESS ;
  case GLU_Qmoment_File :
    if( read_GLU_Qmoment( Input ) == FAILURE ) {
      return FAILURE ;
    }
    // set Lt
    if( init_LT( &Input -> Data , Input -> Traj ) == FAILURE ) {
      return FAILURE ;
    }
    return SUCCESS ;
  }
  return SUCCESS ;
}
