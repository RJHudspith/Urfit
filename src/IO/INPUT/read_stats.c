/**
   @file read_stats.c
   @brief read the input file for the statistics information

   Expected inputs
   Resample = {BootStrap,JackKnife,Raw}
   Nboots = %d
   
   // correlation matrix stuff
   CovDiv = {true,false}
   CovBal = {true,false}
   CovEva = %f
 */
#include "gens.h"

#include "read_inputs.h"

// get the resampling type
static int
get_resample( struct input_params *Input ,
	      const struct flat_file *Flat ,
	      const size_t Ntags )
{
  size_t tag = 0 ;
  if( ( tag = tag_search( Flat , "Resample" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] Resample not found in input file!\n" ) ;
    return FAILURE ;
  }
  if( are_equal( Flat[tag].Value , "BootStrap" ) ) {
    Input -> Data.Restype = BootStrap ;
  } else if( are_equal( Flat[tag].Value , "JackKnife" ) ) {
    Input -> Data.Restype = JackKnife ;
  } else if( are_equal( Flat[tag].Value , "Raw" ) ) {
    Input -> Data.Restype = Raw ;
  } else {
    return FAILURE ;
  }
  return SUCCESS ;
}

int
get_stats( struct input_params *Input ,
	   const struct flat_file *Flat ,
	   const size_t Ntags )
{
  size_t tag = 0 ;
  // set the resample type
  if( get_resample( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }
  if( ( tag = tag_search( Flat , "Nboots" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] Nboots not found in input file!\n" ) ;
    return FAILURE ;
  }
  char *endptr ;
  Input -> Data.Nboots = strtod( Flat[tag].Value , &endptr ) ;

  fprintf( stdout , "\n[INPUTS] summary for statistics\n" ) ;
  switch( Input -> Data.Restype ) {
  case Raw : fprintf( stdout , "[INPUTS] Raw resampling\n" ) ; break ;
  case JackKnife :
    fprintf( stdout , "[INPUTS] JackKnife resampling\n" ) ; break ;
  case BootStrap :
    fprintf( stdout , "[INPUTS] BootStrap resampling\n" ) ;
    fprintf( stdout , "[INPUTS] Using %zu bootstraps\n" ,
	     Input -> Data.Nboots ) ;
    break ;
  }

  // if we are doing a correlated fit then we care about the correlation
  // matrix stuff, otherwise not so much
  Input -> Data.Cov.Divided_Covariance = false ;
  Input -> Data.Cov.Column_Balanced = false ;
  Input -> Data.Cov.Eigenvalue_Tol = 1E-8 ;
  
  if( Input -> Fit.Corrfit == CORRELATED ) {
    // are we performing divided covariance?
    if( ( tag = tag_search( Flat , "CovDiv" , 0 , Ntags ) ) == Ntags ) {
      fprintf( stderr , "[INPUTS] CovDiv not found in input file!\n" ) ;
      return FAILURE ;
    }
    Input -> Data.Cov.Divided_Covariance = \
      are_equal( Flat[tag].Value , "true" ) ;
    // are we balancing the columns of the covariance matrix?
    if( ( tag = tag_search( Flat , "CovBal" , 0 , Ntags ) ) == Ntags ) {
      fprintf( stderr , "[INPUTS] CovBal not found in input file!\n" ) ;
      return FAILURE ;
    }
    Input -> Data.Cov.Column_Balanced = are_equal( Flat[tag].Value , "true" ) ;
    // svd eigenvalue tolerance
    if( ( tag = tag_search( Flat , "CovEva" , 0 , Ntags ) ) == Ntags ) {
      fprintf( stderr , "[INPUTS] CovEva not found in input file!\n" ) ;
      return FAILURE ;
    }
    Input -> Data.Cov.Eigenvalue_Tol = strtod( Flat[tag].Value , &endptr ) ;

    if( Input -> Data.Cov.Divided_Covariance == true ) {
      fprintf( stdout , "[INPUTS] Using Divided correlation matrix ala Michaels\n" ) ;
    }
    if( Input -> Data.Cov.Column_Balanced == true ) {
      fprintf( stdout , "[INPUTS] Using Column-Balanced SVD for inverse correlation matrix\n" ) ;
    }
    fprintf( stdout , "[INPUTS] Filtering out SVD 'Eigenvalues' that are less than %e \n" , Input -> Data.Cov.Eigenvalue_Tol ) ;
  }
  
  return SUCCESS ;
}
