/**
   genetic algorithm fit

   The idea is to reduce the chisq function

   We start with an initial population randomly assigned near the 
   guesses for the fit.

   We sort this w.r.t the chisq and choose this as our breeding 
   population as they are the best fit.

   Our breeding step is defined by an elitism model 

   Keep the top NBREED - these are the parents
   
   For the breeding Pop ( Ngen - NBREED - NMUTANTS )
   inheritance is an average of the mother and fathers genes

   For the mutants we add random gaussian noise to the NBREED
   population, pulled from the width of the breeding population

   The best fit parameter is to be found in the first gene

   We stop our fit when all the members of NBREED have the same chisq
   meaning we are just adding noise to the mutants

   Addendum : I don't know how genetics works but this seems to give 
   the right answer, please excuse use of genes when chromosomes or 
   whatever should be used.

   5/01/2017 - Compute the noise based on the s.d of the kept population, changed to a max and min of the population (decent idea but could do with tweaking)

   6/01/2017 - Made the mutation population smaller and tuned the ratio of NBREED to NCHILD for a standard use case, troubling how dependent on the initial population this thing is.
 */
#include "gens.h"

#include <string.h> // memcpy

#include "chisq.h"
#include "ffunction.h"

//#define VERBOSE

#define NGEN (2048) // number of generations in gene pool

// percentage that persist from the previous generation
#define NBREED ((int)(4)) 

// number of parents of a child
#define NPARENTS (int)(4)

#define NTOURNAMENT (int)(15)

// probability of mutation of the whole pop
#define PMUTANT (0.15)

// It turns out this one is very important!! WHY???
#define NOISE (0.1) // guesses * gaussian of width NOISE to start our run

// defaults to selection sort but can use insertion sort here too
#define INSERTION_SORT

// struct for holding our gene pool
struct genes {
  double *g ;
  double chisq ;
  size_t Nlogic ;
} ;

// struct for storing the standard deviation of each parameter
struct sigma {
  double *s ;
} ;

// compute the chi^2 of our fit functions
static double
compute_chi( struct ffunction *f2 , 
	     const struct ffunction f1 ,
	     const struct fit_descriptor *fdesc ,
	     const double *fparam ,
	     const void *data ,
	     const double **W )
{
  // copy f1 to f2 our temporary fit function
  copy_ffunction( f2 , f1 ) ;
  size_t i ;
  for( i = 0 ; i < fdesc -> Nlogic ; i++ ) {
    f2 -> fparams[i] = fparam[i] ;
  }
  fdesc -> F( f2 -> f , data , f2 -> fparams ) ;
  return compute_chisq( *f2 , W , f2 -> CORRFIT ) ;
}

// simple insertion sort is enough for small population which is
// likely to be roughly sorted already
static void
insertion_sort_GA( struct genes *G )
{
  size_t i , j , Nlogic = G -> Nlogic ;
  double gtemp[ Nlogic ] , x ;
  int hole ;
  for( i = 1 ; i < NGEN ; i++ ) {
    // set up the temporaries
    x = G[i].chisq ;
    for( j = 0 ; j < Nlogic ; j++ ) {
      gtemp[j] = G[i].g[j] ;
    }
    // set the hole we will insert into
    hole = (int)i - 1 ;
    while( hole >= 0 && G[hole].chisq > x ) {
      for( j = 0 ; j < Nlogic ; j++ ) {
	G[ hole + 1 ].g[ j ] = G[ hole ].g[ j ] ;
      }
      G[ hole + 1 ].chisq = G[ hole ].chisq ;
      
      hole-- ;
    }
    // swap back
    for( j = 0 ; j < Nlogic ; j++ ) {
      G[hole+1].g[j] = gtemp[j] ;
    }
    G[hole+1].chisq = x ;
  }
  return ;
}

#ifdef verbose
// print out the whole population
static void
print_population( struct genes *G )
{
  size_t i , j ;
  for( i = 0 ; i < NGEN ; i++ ) {
    for( j = 0 ; j < G[i].Nlogic ; j++ ) {
      printf( " %f " , G[i].g[j] ) ;
    }
    printf( " :: %e \n" , G[i].chisq ) ;
  }
  return ;
}
#endif

// computes the width of the distribution
static void
compute_sigma( struct sigma *Sig ,
	       const struct genes *G ,
	       const size_t Nlogic )
{
  size_t i , j ;
  for( i = 0 ; i < Nlogic ; i++ ) {
    double max = G[0].g[i] , min = G[0].g[i] ;
    for( j = 1 ; j < NGEN ; j++ ) {
      if( G[j].g[i] < min ) {
	min = G[j].g[i] ;
      }
      if( G[j].g[i] > max ) {
	max = G[j].g[i] ;
      }
    }
    Sig -> s[i] = ( max - min ) ;
  }
  return ;
}

// perform simple selection
static size_t
tournament_selection( const struct genes *G ,
		      gsl_rng *r )
{
  // select a bunch of NTOURNAMENT elements at random
  size_t i , idx_best = gsl_rng_uniform_int( r , NGEN ) ;
  double tbest = G[ idx_best ].chisq ;
  for( i = 0 ; i < NTOURNAMENT ; i++ ) {
    size_t idx = gsl_rng_uniform_int( r , NGEN ) ;
    if( G[ idx ].chisq < tbest ) {
      tbest = G[idx].chisq ;
      idx_best = idx ;
    }
  }
  return idx_best ;
}

#ifdef NORMALISE_POP
static void
normalise_pop( struct genes *G )
{
  size_t i ;
  register double Nsum = 0.0 ;
  for( i = 0 ; i < NGEN ; i++ ) {
    Nsum += G[i].chisq ;
  }
  for( i = 0 ; i < NGEN ; i++ ) {
    G[i].chisq /= Nsum ;
  }
  for( i = 1 ; i < NGEN ; i++ ) {
    G[i].chisq += G[i-1].chisq ;
  }
  return ;
}
#endif

// perform a minimisation using a genetic algorithm, parameters are 
// at the top to be played around with. Does not use any derivative 
// information!
int
ga_iter( void *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // point to the fit descriptor struct
  struct fit_descriptor *Fit = (struct fit_descriptor*)fdesc ;
    
  if( Fit -> Nlogic == 0 ) return SUCCESS ; // do nothing
  
  // counters and max iterations GAMAX
  size_t iters = 0 , i , j ;
  const size_t GAMAX = 5000 ;

  // allocate the fitfunction
  struct ffunction f2 = allocate_ffunction( Fit -> Nlogic , 
					    Fit -> f.N ) ;

  // get priors
  Fit -> f.Prior = Fit -> Prior ;

  // evaluate the function, its first and second derivatives
  Fit -> F( Fit -> f.f , data , Fit -> f.fparams ) ;
  Fit -> f.chisq = compute_chisq( Fit -> f , W , Fit -> f.CORRFIT ) ;

#ifdef VERBOSE
  printf( "[GA] Using a population of %d \n" , NGEN ) ;
  printf( "[GA] Keeping %d elites \n" , NBREED ) ;
  printf( "[GA] Breeding elites into %d Children\n" , NCHILD ) ;
  printf( "[GA] Mutating the above into a pop of %d \n" , NMUTANTS ) ;
#endif
  
  // gene pool
  struct genes *G = NULL ;
  struct sigma Sig ;
  gsl_rng *r = NULL ;

  G = malloc( NGEN * sizeof( struct genes ) ) ;
  for( i = 0 ; i < NGEN ; i++ ) {
    G[i].g = NULL ;
    G[i].g = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  }

  Sig.s = malloc( Fit -> Nlogic * sizeof( double ) ) ;
  
  // get a seed from urandom
  size_t Seed ;
  FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
  if( fread( &Seed , sizeof( Seed ) , 1 , urandom ) != 1 ) {
    fprintf( stderr , "[GA] urandom read failure! \n" ) ;
    goto memfree ;
  }
  fclose( urandom ) ;

  // initialise gsl rng is the mersenne twister I believe
  gsl_rng_env_setup( ) ;
  r = gsl_rng_alloc( gsl_rng_default ) ;
  gsl_rng_set( r , Seed ) ;

  // initialise the population as gaussian noise around initial
  // guesses that are gaussian distibuted with sigma of NOISE
  // defined at the top of the file
  for( i = 0 ; i < NGEN ; i++ ) {
    G[i].Nlogic = Fit -> Nlogic ;
    for( j = 0 ; j < Fit -> Nlogic ; j++ ) {
      // if we have priors we use their errors like noise vectors
      if( Fit -> Prior[j].Initialised == true ) {
	G[i].g[j] = Fit -> f.fparams[j] * 
	  ( 1 + gsl_ran_gaussian( r , fabs( Fit -> f.fparams[j] ) *
				  Fit -> Prior[j].Err ) ) ;
      } else {
	G[i].g[j] = Fit -> f.fparams[j] * 
	  ( 1 + gsl_ran_gaussian( r , fabs( Fit -> f.fparams[j] ) * NOISE ) ) ;
      }
    }
    G[i].chisq = compute_chi( &f2 , Fit -> f , Fit , G[i].g , data , W ) ;
  }

  // normalise it
  #ifdef NORMALISE_POP
  normalise_pop( G ) ;
  #endif

  // sort by chisq
  insertion_sort_GA( G ) ;
  
  #ifdef verbose
  print_population( G ) ;
  #endif
  
  // iterate the algorithm
  double chisq_diff = 10 ;
  while( chisq_diff > TOL && iters < GAMAX ) {


    for( i = NBREED ; i < NGEN ; i++ ) {
      size_t parent[ NPARENTS ] , k ;
      for( k = 0 ; k < NPARENTS ; k++ ) {
	parent[k] = tournament_selection( G , r ) ;
      }
      // breed into the new population
      for( j = 0 ; j < Fit -> Nlogic ; j++ ) {
	G[i].g[j] = 0.0 ;
	for( k = 0 ; k < NPARENTS ; k++ ) {
	  G[i].g[j] += G[ parent[k] ].g[j] ;
	}
        G[i].g[j] /= NPARENTS ;
      }
      G[i].chisq = compute_chi( &f2 , Fit -> f , Fit ,
				G[i].g , data , W ) ;
    }


        // sort for the next breeding population
    #ifdef INSERTION_SORT
    insertion_sort_GA( G ) ;
    #else
    selection_sort_GA( G ) ;
    #endif
    
    // look for the elites to be static as an 
    // opportunity to turn off the fit
    //printf( "%e %e \n" , G[0].chisq , G[NGEN-1].chisq  ) ;
    chisq_diff = fabs( G[0].chisq - G[NGEN-1].chisq ) ;

    // compute the width of the distribution as-is and use
    // that as a guide for how to add noise
    compute_sigma( &Sig , G , Fit -> Nlogic ) ;

    // perform a mutation on the entire pop with mutation
    // chance PMUTANT
    for( i = 0 ; i < NGEN ; i++ ) {
      if( gsl_rng_uniform( r ) < PMUTANT ) {
	for( j = 0 ; j < Fit -> Nlogic ; j++ ) {
	  // if we have prior information we use it
	  if( Fit -> Prior[j].Initialised == true ) {
	    G[i].g[j] = G[ i ].g[j] * 
	      ( 1 + gsl_ran_gaussian( r , Fit -> Prior[j].Err ) ) ;
	  } else {
	    G[i].g[j] = G[ i ].g[j] * 
	      ( 1 + gsl_ran_gaussian( r , Sig.s[j] ) ) ;
	  }
	}
	G[i].chisq = compute_chi( &f2 , Fit -> f , Fit ,
				  G[i].g , data , W ) ;
      }
    }    

    // sort for the next breeding population
    insertion_sort_GA( G ) ;
    
    iters++ ;
  }

  // set the fit parameters as the best gene and set the chisq
  for( i = 0 ; i < Fit -> Nlogic ; i++ ) {
    Fit -> f.fparams[i] = G[0].g[i] ;
    #ifdef VERBOSE
    printf( "FPARAMS_%zu :: %f \n" , i , Fit -> f.fparams[i] ) ;
    #endif
  }
  Fit -> f.chisq = G[0].chisq ;
  
  // best fit parameter will be in population 1
  if( iters == GAMAX ) {
    printf( "\n[GA] chisq :: %e || %zu iterations (GAMAX) || %e \n" ,
	    G[0].chisq , iters , chisq_diff ) ;
  } else {
    printf( "\n[GA] chisq :: %e || %zu iterations \n" ,
	    G[0].chisq , iters ) ;
  }

 memfree :
  
  // cleanse the gene pool
  if( G != NULL ) {
    for( i = 0 ; i < NGEN ; i++ ) {
      if( G[i].g != NULL ) {
	free( G[i].g ) ;
      }
    }
    free( G ) ;
  }

  // free the fitfunction
  free_ffunction( &f2 , Fit -> Nlogic ) ;

  // free the rng
  if( r != NULL ) {
    gsl_rng_free( r ) ;
  }
  
  // free the variance
  if( Sig.s != NULL ) {
    free( Sig.s ) ;
  }

  return iters ;
}
