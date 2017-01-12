/**
   @file read_fit.c
   @brief read the fit struct from the input file

   Fit info expected to exist somewhere in the input file:

   FitDef = function we are fitting to
   Fit_NM = n,m
   FitCorr = {CORRELATED,UNCORRELATED,UNWEIGHTED}
   FitSims = , , , , , -> Up to Nlogic of these if there is only NULL at the start we have 0
   Prior = index,val,err -- can have loads of these
   FitTol = tolerance we minimize to
   FitMin = minimizer we use {CG,GA,LM,SD}

   @warning code needs to be called after the trajectories have been figured out
 */
#include "gens.h"

#include <string.h>

// minimizers
#include "CG.h"
#include "GA.h"
#include "LM.h"
#include "SD.h"

#include "fit_chooser.h"
#include "read_inputs.h"

// linked list for the dimensions
struct node {
  size_t sim ;
  struct node *next ;
} ;

// get the function we are fitting our data to
static int
get_fitDef(  struct input_params *Input ,
	     const struct flat_file *Flat ,
	     const size_t Ntags )
{
  size_t tag = 0 ;
  if( ( tag = tag_search( Flat , "FitDef" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] FitDef not found in input file!\n" ) ;
    return FAILURE ;
  }
  if( are_equal( Flat[tag].Value , "EXP" ) ) {
    Input -> Fit.Fitdef = EXP ;
  } else if( are_equal( Flat[tag].Value , "EXP_PLUSC" ) ) {
    Input -> Fit.Fitdef = EXP_PLUSC ;
  } else if( are_equal( Flat[tag].Value , "POLY" ) ) {
    Input -> Fit.Fitdef = POLY ;
  } else if( are_equal( Flat[tag].Value , "PADE" ) ) {
    Input -> Fit.Fitdef = PADE ;
  } else if( are_equal( Flat[tag].Value , "COSH" ) ) {
    Input -> Fit.Fitdef = COSH ;
  } else if( are_equal( Flat[tag].Value , "SINH" ) ) {
    Input -> Fit.Fitdef = SINH ;
  } else {
    return FAILURE ;
  }
  return SUCCESS ;
}

// get the function we are fitting our data to
static int
get_fitCorr( struct input_params *Input ,
	     const struct flat_file *Flat ,
	     const size_t Ntags )
{
  size_t tag = 0 ;
  if( ( tag = tag_search( Flat , "FitCorr" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] FitCorr not found in input file!\n" ) ;
    return FAILURE ;
  }
  if( are_equal( Flat[tag].Value , "CORRELATED" ) ) {
    Input -> Fit.Corrfit = CORRELATED ;
  } else if( are_equal( Flat[tag].Value , "UNCORRELATED" ) ) {
    Input -> Fit.Corrfit = UNCORRELATED ;
  } else if( are_equal( Flat[tag].Value , "UNWEIGHTED" ) ) {
    Input -> Fit.Corrfit = UNWEIGHTED ;
  } else {
    return FAILURE ;
  }
  return SUCCESS ;
}

// get the minimiser for the fit
static int
get_fitMin( struct input_params *Input ,
	    const struct flat_file *Flat ,
	    const size_t Ntags )
{
  size_t tag = 0 ;
  if( ( tag = tag_search( Flat , "FitMin" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] FitMin not found in input file!\n" ) ;
    return FAILURE ;
  }
  if( are_equal( Flat[tag].Value , "CG" ) ) {
    Input -> Fit.Minimize = cg_iter ;
  } else if( are_equal( Flat[tag].Value , "GA" ) ) {
    Input -> Fit.Minimize = ga_iter ;
  } else if( are_equal( Flat[tag].Value , "LM" ) ) {
    Input -> Fit.Minimize = lm_iter ;
  } else if( are_equal( Flat[tag].Value , "SD" ) ) {
    Input -> Fit.Minimize = sd_iter ;
  } else {
    return FAILURE ;
  }
  return SUCCESS ;
}

// get the number of sim params and put them in a list
static struct node *
get_Nsims( size_t *Nsims , const char *Value )
{
  struct node *head = NULL , *curr ;
  char *tok = strtok( (char*)Value , "," ) , *endptr ;
  if( tok == NULL ) return NULL ;
  *Nsims = 1 ;
  curr = (struct node*)malloc( sizeof( struct node ) ) ;
  curr -> sim = strtol( tok , &endptr , 10 ) ;
  curr -> next = head ;
  head = curr ;
  while( ( tok = strtok( NULL , "," ) ) != NULL ) {
    curr = (struct node*)malloc( sizeof( struct node ) ) ;
    curr -> sim = strtol( tok , &endptr , 10 ) ;
    curr -> next = head ;
    head = curr ;
    *Nsims = *Nsims + 1 ;
  }
  return head ;
}

// get the priors from the file
static int
get_Priors( struct input_params *Input ,
	    const struct flat_file *Flat ,
	    const size_t Ntags )
{
  size_t i , Nprior = 0 , block_idx = 0 , Ncommon = 0 ;
  // compute the nqumber of logical parameters we can have
  for( i = 0 ; i < Input -> Fit.Nparam ; i++ ) {
    if( Input -> Fit.Sims[i] == true ) {
      Ncommon++ ;
    }
  }
  Input-> Fit.Nlogic = Input -> Ntraj * ( Input -> Fit.Nparam - Ncommon )
    + Ncommon ;
  
  // allocate and initialise to no priors
  Input -> Fit.Prior = malloc( Input -> Fit.Nlogic * sizeof( struct prior ) ) ;
  for( i = 0 ; i < Input -> Fit.Nlogic ; i++ ) {
    Input -> Fit.Prior[i].Initialised = false ;
    Input -> Fit.Prior[i].Val = UNINIT_FLAG ;
    Input -> Fit.Prior[i].Err = UNINIT_FLAG ;
  }
  
  // count how many priors there are
  while( ( block_idx = tag_search( Flat , "Prior" , block_idx , Ntags ) )
	 != Ntags ) {

    // tokenize the priors
    char *tok = strtok( Flat[block_idx].Value , "," ) , *endptr ;
    const size_t idx = strtol( tok , &endptr , 10 ) ;

    // poke into the prior struct
    Input -> Fit.Prior[idx].Initialised = true ;
    tok = strtok( NULL , "," ) ;
    Input -> Fit.Prior[idx].Val = strtod( tok , &endptr ) ;
    tok = strtok( NULL , "," ) ;
    Input -> Fit.Prior[idx].Err = strtod( tok , &endptr ) ;
    tok = strtok( NULL , "," ) ;

    if( Nprior >= Input -> Fit.Nlogic ) {
      fprintf( stderr , "[INPUTS] too many priors in input file\n" ) ;
      return FAILURE ;
    }
    Nprior++ ;
    block_idx ++ ;
  }
  return SUCCESS ;
}

// fill out the fit_info struct of Input
int
get_fit( struct input_params *Input ,
	 const struct flat_file *Flat ,
	 const size_t Ntags )
{
  // get the function we are fitting to
  if( get_fitDef( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }
  // get what type of fit we are doing, correlated, uncorrelated, unweighted
  if( get_fitCorr( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }
  // get the minimizer
  if( get_fitMin( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }
  // directly read the tolerance
  char *endptr = NULL ;
  size_t tag = 0 ;
  if( ( tag = tag_search( Flat , "FitTol" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] FitTol not found in input file!\n" ) ;
    return FAILURE ;
  }
  Input -> Fit.Tol = strtod( Flat[tag].Value , &endptr ) ;
  
  // directly read fit n and m
  if( ( tag = tag_search( Flat , "Fit_NM" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] Fit_NM not found in input file!\n" ) ;
    return FAILURE ;
  }
  char *tok = strtok( Flat[ tag ].Value , "," ) ;
  Input -> Fit.N = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;
  Input -> Fit.M = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;

  // set Nparam
  Input -> Fit.Nparam = get_Nparam( Input -> Fit ) ;

  // read and allocate the number of simparams
  if( ( tag = tag_search( Flat , "FitSims" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] FitSims not found in input file!\n" ) ;
    return FAILURE ;
  }
  size_t Nsims = 0 ;
  struct node *Node = NULL ;
  if( !are_equal( Flat[tag].Value , "NULL" ) ) {
    Node = get_Nsims( &Nsims , Flat[tag].Value ) ;
  }
  Input -> Fit.Sims = malloc( Input -> Fit.Nparam * sizeof( bool ) ) ;

  // set them all to false
  size_t i ;
  for( i = 0 ; i < Input -> Fit.Nparam ; i++ ) {
    Input -> Fit.Sims[i] = false ;
  }
  // we have simultaneous parameters and we poke them in
  if( Node != NULL ) {
    for( i = 0 ; i < Nsims ; i++ ) {
      if( Node -> sim < Input -> Fit.Nparam ) {
	Input -> Fit.Sims[ Node -> sim ] = true ;
	free( Node ) ;
	Node = Node -> next ;
      } else {
	fprintf( stderr , "[INPUTS] specified simultaneous fit parameter %zu"
		 "is greater than the number of fit parameters %zu \n" ,
		 Node -> sim , Input -> Fit.Nparam ) ;
	return FAILURE ;
      }
    }
  }

  // get the priors
  if( get_Priors( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }

  // print out a summary of the priors and simultaneous parameters
  fprintf( stdout , "\n[INPUTS] Fit parameter summary\n" ) ;
  fprintf( stdout , "\n[INPUTS] Fit tolerance %e \n" , Input -> Fit.Tol ) ;
  for( i = 0 ; i < Input -> Fit.Nlogic ; i++ ) {
    fprintf( stdout , "[INPUTS] PRIOR %zu %f %f \n" , i ,
	     Input -> Fit.Prior[i].Initialised ,
	     Input -> Fit.Prior[i].Val ,
	     Input -> Fit.Prior[i].Err ) ;
  }
  for( i = 0 ; i < Input -> Fit.Nparam ; i++ ) {
    if( Input -> Fit.Sims[i] == true ) {
      fprintf( stdout , "[INPUTS] Parameter %zu IS simultaneous \n" , i ) ;
    } else {
      fprintf( stdout , "[INPUTS] Parameter %zu IS NOT simultaneous \n" , i ) ;
    }
  }
  fprintf( stdout , "\n" ) ;
  
  return SUCCESS ;
}
