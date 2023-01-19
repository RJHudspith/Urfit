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

   Guess_0 = val
   Guess_1 = val
   ....
   Must have NLOGIC of these, need to be in order

   @warning code needs to be called after the trajectories have been figured out
 */
#include "gens.h"

#include <string.h>

// minimizers
#include "CG.h"
#include "GA.h"
#include "GLS.h"
#include "GLS_pade.h"
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
  if( are_equal( Flat[tag].Value , "ALPHA_D0" ) ) {
    Input -> Fit.Fitdef = ALPHA_D0 ;
  } else if( are_equal( Flat[tag].Value , "ALPHA_D0_MULTI" ) ) {
    Input -> Fit.Fitdef = ALPHA_D0_MULTI ;
  } else if( are_equal( Flat[tag].Value , "ADLERALPHA_D0" ) ) {
    Input -> Fit.Fitdef = ADLERALPHA_D0 ;
  } else if( are_equal( Flat[tag].Value , "ADLERALPHA_D0_MULTI" ) ) {
    Input -> Fit.Fitdef = ADLERALPHA_D0_MULTI ;
  } else if( are_equal( Flat[tag].Value , "C4C7" ) ) {
    Input -> Fit.Fitdef = C4C7 ;
  } else if( are_equal( Flat[tag].Value , "CORNELL" ) ) {
    Input -> Fit.Fitdef = CORNELL ;
  } else if( are_equal( Flat[tag].Value , "CORNELL_V2" ) ) {
    Input -> Fit.Fitdef = CORNELL_V2 ;
  } else if( are_equal( Flat[tag].Value , "COSH" ) ) {
    Input -> Fit.Fitdef = COSH ;
  } else if( are_equal( Flat[tag].Value , "COSH_PLUSC" ) ) {
    Input -> Fit.Fitdef = COSH_PLUSC ;
  } else if( are_equal( Flat[tag].Value , "EXP" ) ) {
    Input -> Fit.Fitdef = EXP ;
  } else if( are_equal( Flat[tag].Value , "EXP_XINV" ) ) {
    Input -> Fit.Fitdef = EXP_XINV ;
  } else if( are_equal( Flat[tag].Value , "EXP_PLUSC" ) ) {
    Input -> Fit.Fitdef = EXP_PLUSC ;
  } else if( are_equal( Flat[tag].Value , "FVOLCC" ) ) {
    Input -> Fit.Fitdef = FVOLCC ;
  } else if( are_equal( Flat[tag].Value , "FVOL1" ) ) {
    Input -> Fit.Fitdef = FVOL1 ;
  } else if( are_equal( Flat[tag].Value , "FVOL2" ) ) {
    Input -> Fit.Fitdef = FVOL2 ;
  } else if( are_equal( Flat[tag].Value , "FVOL3" ) ) {
    Input -> Fit.Fitdef = FVOL3 ;
  } else if( are_equal( Flat[tag].Value , "HALEXP" ) ) {
    Input -> Fit.Fitdef = HALEXP ;
  } else if( are_equal( Flat[tag].Value , "HLBL_CONT" ) ) {
    Input -> Fit.Fitdef = HLBL_CONT ;
  } else if( are_equal( Flat[tag].Value , "LARGENB" ) ) {
    Input -> Fit.Fitdef = LARGENB ;
  } else if( are_equal( Flat[tag].Value , "NRQCD_EXP" ) ) {
    Input -> Fit.Fitdef = NRQCD_EXP ;
  } else if( are_equal( Flat[tag].Value , "NRQCD_EXP2" ) ) {
    Input -> Fit.Fitdef = NRQCD_EXP2 ;
  } else if( are_equal( Flat[tag].Value , "NOFIT" ) ) {
    Input -> Fit.Fitdef = NOFIT ;
  } else if( are_equal( Flat[tag].Value , "PADE" ) ) {
    Input -> Fit.Fitdef = PADE ;
  } else if( are_equal( Flat[tag].Value , "PEXP" ) ) {
    Input -> Fit.Fitdef = PEXP ;
  } else if( are_equal( Flat[tag].Value , "POLY" ) ) {
    Input -> Fit.Fitdef = POLY ;
  } else if( are_equal( Flat[tag].Value , "POLES" ) ) {
    Input -> Fit.Fitdef = POLES ; 
  } else if( are_equal( Flat[tag].Value , "PP_AA" ) ) {
    Input -> Fit.Fitdef = PP_AA ;
  } else if( are_equal( Flat[tag].Value , "PP_AA_EXP" ) ) {
    Input -> Fit.Fitdef = PP_AA_EXP ;
  } else if( are_equal( Flat[tag].Value , "PP_AA_WW" ) ) {
    Input -> Fit.Fitdef = PP_AA_WW ;
  } else if( are_equal( Flat[tag].Value , "PP_AA_WW_R2" ) ) {
    Input -> Fit.Fitdef = PP_AA_WW_R2 ;
  } else if( are_equal( Flat[tag].Value , "PPAA" ) ) {
    Input -> Fit.Fitdef = PPAA ;
  } else if( are_equal( Flat[tag].Value , "SINH" ) ) {
    Input -> Fit.Fitdef = SINH ;
  } else if( are_equal( Flat[tag].Value , "TANH" ) ) {
    Input -> Fit.Fitdef = TANH ;
  } else if( are_equal( Flat[tag].Value , "SOL" ) ) {
    Input -> Fit.Fitdef = SOL ;
  } else if( are_equal( Flat[tag].Value , "SU2_SHITFIT" ) ) {
    Input -> Fit.Fitdef = SU2_SHITFIT ;
  } else if( are_equal( Flat[tag].Value , "SUN_CONT" ) ) {
    Input -> Fit.Fitdef = SUN_CONT ;
  } else if( are_equal( Flat[tag].Value , "QCORR_BESSEL" ) ) {
    Input -> Fit.Fitdef = QCORR_BESSEL ;
  } else if( are_equal( Flat[tag].Value , "QSLAB" ) ) {
    Input -> Fit.Fitdef = QSLAB ;
  } else if( are_equal( Flat[tag].Value , "QSLAB_FIXED" ) ) {
    Input -> Fit.Fitdef = QSLAB_FIXED ;
  } else if( are_equal( Flat[tag].Value , "QSUSC_SU2" ) ) {
    Input -> Fit.Fitdef = QSUSC_SU2 ;
  } else if( are_equal( Flat[tag].Value , "UDCB_HEAVY" ) ) {
    Input -> Fit.Fitdef = UDCB_HEAVY ;
  } else if( are_equal( Flat[tag].Value , "ZV_EXP" ) ) {
    Input -> Fit.Fitdef = ZV_EXP ;
  } else {
    fprintf( stderr , "[INPUT] Fit %s not recognised\n" , Flat[tag].Value ) ;
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
    fprintf( stderr , "[INPUT] FitCorr %s not recognised\n" , Flat[tag].Value ) ;
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
  } else if( are_equal( Flat[tag].Value , "GLS" ) ) {
    if( Input -> Fit.Fitdef == POLY ||
	Input -> Fit.Fitdef == NOFIT ||
	Input -> Fit.Fitdef == POLES ||
	Input -> Fit.Fitdef == LARGENB ) {
      Input -> Fit.Minimize = gls_iter ;
    } else {
      fprintf( stderr , "[INPUTS] GLS only supports POLY/POLES type fit\n" ) ;
      return FAILURE ;
    }
  } else if( are_equal( Flat[tag].Value , "GLS_pade" ) ) {
    if( Input -> Fit.Fitdef == PADE || Input -> Fit.Fitdef == NOFIT ) {
      Input -> Fit.Minimize = gls_pade_iter ;
    } else {
      fprintf( stderr , "[INPUTS] GLS only supports PADE type fit\n" ) ;
      return FAILURE ;
    }
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
  Input -> Fit.Nlogic = Input -> Data.Nsim * ( Input -> Fit.Nparam - Ncommon )
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

    if( idx >= Input -> Fit.Nlogic ) {
      fprintf( stderr ,
	       "[INPUTS] prior index %zu is greater than Nlogic %zu \n" ,
	       idx , Input -> Fit.Nlogic ) ;
      return FAILURE ;
    }

    // poke into the prior struct
    Input -> Fit.Prior[idx].Initialised = true ;
    tok = strtok( NULL , "," ) ;
    Input -> Fit.Prior[idx].Val = strtod( tok , &endptr ) ;
    tok = strtok( NULL , "," ) ;
    Input -> Fit.Prior[idx].Err = strtod( tok , &endptr ) ;

    if( Nprior >= Input -> Fit.Nlogic ) {
      fprintf( stderr , "[INPUTS] too many priors in input file\n" ) ;
      return FAILURE ;
    }
    Nprior++ ;
    block_idx ++ ;
  }
  Input -> Fit.Nprior = Nprior ;
  
  return SUCCESS ;
}

// get the guesses from the input file, called after get_Priors
static int
get_Guesses( struct input_params *Input ,
	     const struct flat_file *Flat ,
	     const size_t Ntags )
{
  // allocate and initialise to no guesses
  Input -> Fit.Guesses_Initialised = false ;
  Input -> Fit.Guess = malloc( Input -> Fit.Nlogic * sizeof( double ) ) ;

  // as the guesses should be in order this becomes a finite-state machine
  size_t block_idx = 0 , i ;
  for( i = 0 ; i < Input -> Fit.Nlogic ; i++ ) {
    Input -> Fit.Guess[i] = UNINIT_FLAG ;
    char guess_str[ 32 ] , *endptr ;
    sprintf( guess_str , "Guess_%zu" , i ) ;
    
    block_idx = tag_search( Flat , guess_str , block_idx , Ntags ) ;

    // if we can't fint the specific tag we leave not initialising any
    // of the guesses
    if( block_idx == Ntags ) {
      printf( "%s not found\n" , guess_str ) ;
      goto end ;
    }
    
    // set the value
    Input -> Fit.Guess[i] = strtod( Flat[block_idx].Value , &endptr ) ;
  }

  // we have initialised all the guesses so we set this flag
  Input -> Fit.Guesses_Initialised = true ;

 end :
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
  Input -> Fit.M = strtol( tok , &endptr , 10 ) ;

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
	//free( Node ) ;
	Node = Node -> next ;
      } else {
	fprintf( stderr , "[INPUTS] specified simultaneous fit parameter %zu"
		 "is greater than the number of fit parameters %zu \n" ,
		 Node -> sim , Input -> Fit.Nparam ) ;
	while( Node != NULL ) {
	  free( Node ) ;
	  Node = Node -> next ;
	}
	return FAILURE ;
      }
    }
  }

  // get the priors
  if( get_Priors( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }

  // get the guesses
  if( get_Guesses( Input , Flat , Ntags ) == FAILURE ) {
    return FAILURE ;
  }

  // print out a summary of the priors and simultaneous parameters
  fprintf( stdout , "\n[INPUTS] Summary for fit parameters\n" ) ;
  fprintf( stdout , "[INPUTS] Fit tolerance %e \n" , Input -> Fit.Tol ) ;
  for( i = 0 ; i < Input -> Fit.Nlogic ; i++ ) {
    if( Input -> Fit.Prior[i].Initialised == true ) {
      fprintf( stdout , "[INPUTS] PRIOR %zu %f %f \n" , i ,
	       Input -> Fit.Prior[i].Val ,
	       Input -> Fit.Prior[i].Err ) ;
    }
  }
  if( Input -> Fit.Guesses_Initialised == true ) {
    for( i = 0 ; i < Input -> Fit.Nlogic ; i++ ) {
      fprintf( stdout , "[INPUTS] Guess %zu %e \n" , i ,
	       Input -> Fit.Guess[i] ) ;
    }
  }
  for( i = 0 ; i < Input -> Fit.Nparam ; i++ ) {
    if( Input -> Fit.Sims[i] == true ) {
      fprintf( stdout , "[INPUTS] Parameter %zu IS simultaneous \n" , i ) ;
    } else {
      fprintf( stdout , "[INPUTS] Parameter %zu IS NOT simultaneous \n" , i ) ;
    }
  }
  switch( Input -> Fit.Corrfit ) {
  case CORRELATED :
    fprintf( stdout , "[INPUTS] performing a CORRELATED fit \n" ) ;
    break ;
  case UNCORRELATED :
    fprintf( stdout , "[INPUTS] performing an UNCORRELATED fit \n" ) ;
    break ;
  case UNWEIGHTED :
    fprintf( stdout , "[INPUTS] performing an UNWEIGHTED fit \n" ) ;
    break ;
  }
  
  return SUCCESS ;
}
