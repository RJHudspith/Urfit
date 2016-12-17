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
     P(1/3) inheritance comes from the mother
     P(1/3) inheritance comes from the father
     P(1/3) inheritance is an average of the mother and fathers genes

   For the mutants we add random gaussian noise to the NBREED
   population

   The best fit parameter is to be found in the first gene

   We stop our fit when all the members of NBREED have the same chisq
   meaning we are just adding noise to the mutants

   Addendum : I don't know how genetics works but this seems to give 
   the right answer, please excuse use of genes when chromosomes or 
   whatever should be used.
 */
#include "gens.h"
#include <string.h> // memcpy

#include "chisq.h"
#include "ffunction.h"

// verbose output
//#define VERBOSE

#define Ngen (400) // number of generations in gene pool

#define NBREED (40) // number that persist from the previous generation

#define NCHILD (120) // number of children

#define NMUTANTS (Ngen - NBREED - NCHILD) // number of mutants

#define NOISE (0.2) // guesses * gaussian of width NOISE to start our run

// struct for holding our gene pool
struct genes {
  double *g ;
  double chisq ;
  size_t NPARAMS ;
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
  size_t j ;
  for( j = 0 ; j < f1.NPARAMS ; j++ ) {
    f2 -> fparams[ j ] = fparam[ j ] ;
  }
  fdesc -> F( f2 -> f , data , f2 -> fparams ) ;
  return compute_chisq( *f2 , W , f2 -> CORRFIT ) ;
}

// swap two genes in the gene pool for the insertion sort
static void
swap_genes( struct genes *G1 ,
	    struct genes *G2 )
{
  double tmp[ G1->NPARAMS ] , tchi = G1->chisq ;
  // make a temporary
  memcpy( tmp , G1->g , G1->NPARAMS * sizeof( double ) ) ;
  // swap genes
  memcpy( G1->g , G2->g , G1->NPARAMS * sizeof( double ) ) ;
  G1->chisq = G2->chisq ;
  // 
  memcpy( G2->g , tmp , G1->NPARAMS * sizeof( double ) ) ;
  G2->chisq = tchi ;
  return ;
}

// simple insertion sort is enough for small population which is
// likely to be roughly sorted already
static void
insertion_sort_GA( struct genes *G )
{
  size_t i , hole ;
  for( i = 1 ; i < Ngen ; i++ ) {
    hole = i ;
    while( hole > 0 && G[hole].chisq < G[hole-1].chisq ) {
      swap_genes( G + hole - 1 , G + hole ) ;
      hole-- ;
    }
  }
  return ;
}

#ifdef VERBOSE
// print out the whole population
static void
print_population( struct genes *G )
{
  size_t i , j ;
  for( i = 0 ; i < Ngen ; i++ ) {
    for( j = 0 ; j < G[i].NPARAMS ; j++ ) {
      printf( " %f " , G[i].g[j] ) ;
    }
    printf( " :: %e \n" , G[i].chisq ) ;
  }
  return ;
}
#endif

// perform a minimisation using a genetic algorithm, parameters are 
// at the top to be played around with. Does not use any derivative 
// information!
int
ga_iter( struct fit_descriptor *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL )
{
  // counters and max iterations GAMAX
  size_t iters = 0 , i ;
  const size_t GAMAX = 10000 ;

  // allocate the fitfunction
  struct ffunction f2 = allocate_ffunction( fdesc -> NPARAMS , 
					    fdesc -> f.N ) ;

  // some guesses
  fdesc -> guesses( fdesc -> f.fparams ) ;

  // get priors
  fdesc -> set_priors( fdesc -> f.prior , fdesc -> f.err_prior ) ;

  // evaluate the function, its first and second derivatives
  fdesc -> F( fdesc -> f.f , data , fdesc -> f.fparams ) ;
  fdesc -> f.chisq = compute_chisq( fdesc -> f , W , fdesc -> f.CORRFIT ) ;

  // gene pool
  struct genes *G = malloc( Ngen * sizeof( struct genes ) ) ;

  gsl_rng_env_setup( ) ;
  gsl_rng *r = gsl_rng_alloc( gsl_rng_default ) ;

  size_t Seed ;
  FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
  fread( &Seed , sizeof( Seed ) , 1 , urandom ) ;
  fclose( urandom ) ;
  gsl_rng_set( r , Seed ) ;

  // initialise the population as gaussian noise around initial
  // guesses
  for( i = 0 ; i < Ngen ; i++ ) {
    G[i].NPARAMS = fdesc -> NPARAMS ;
    G[i].g = malloc( fdesc -> NPARAMS * sizeof( double ) ) ;
    size_t j ;
    for( j = 0 ; j < fdesc -> NPARAMS ; j++ ) {
      G[i].g[j] = fdesc -> f.fparams[j] * 
	( 1 + gsl_ran_gaussian( r , NOISE ) ) ;
    }
    G[i].chisq = compute_chi( &f2 , fdesc -> f , fdesc ,
			       G[i].g , data , W ) ;
  }

  // sort by chisq
  insertion_sort_GA( G ) ;

#ifdef VERBOSE
  print_population( G ) ;
#endif

  // iterate the algorithm
  double chisq_diff = 10 ;
  while( chisq_diff > TOL && iters < GAMAX ) {

    // keep the top NBREED breed them to create children
    for( i = NBREED ; i < Ngen - NMUTANTS ; i++ ) {
      size_t j ;
      const size_t father = gsl_rng_uniform_int( r , NBREED ) ;
      const size_t mother = gsl_rng_uniform_int( r , NBREED ) ;
      // if it is above we average we combine otherwise we pick 
      // randomly from the father or mother
      for( j = 0 ; j < fdesc -> NPARAMS ; j++ ) {
	//const double diceroll = gsl_rng_uniform( r ) ;

	// roll the dice to decide where the genes go!
	G[i].g[j] = 0.5 * ( G[ father ].g[j] + G[ mother ].g[j] ) ;
      }
      G[i].chisq = compute_chi( &f2 , fdesc -> f , fdesc ,
				G[i].g , data , W ) ;
    }

    insertion_sort_GA( G ) ;
 
    // recompute sigma, estimating from the breeding population
#if 0
    double sigma[ fdesc -> NPARAMS ] ;
    for( i = 0 ; i < fdesc -> NPARAMS ; i++ ) {
      sigma[i] = fabs( G[0].g[i] - G[NBREED+NCHILD-1].g[i] ) ;
    }
#endif

    // bottom 20 can be randomly mutated from the top 20
    for( i = Ngen - NMUTANTS ; i < Ngen ; i++ ) {
      //const size_t mutate = gsl_rng_uniform_int( r , NBREED ) ;
      const size_t mutate = gsl_rng_uniform_int( r , NBREED ) ;
      size_t j ;
      for( j = 0 ; j < fdesc -> NPARAMS ; j++ ) {
	// perform some noise attenuation to help vary around 
	// most popular solution
	G[i].g[j] = G[ mutate ].g[j] * 
	  ( 1 + gsl_ran_gaussian( r , NOISE ) ) ;
      }
      G[i].chisq = compute_chi( &f2 , fdesc -> f , fdesc ,
				G[i].g , data , W ) ;
    }

    insertion_sort_GA( G ) ;

    // look for a static solution in the elites as an 
    // opportunity to turn off the fit
    chisq_diff = fabs( G[0].chisq - G[NBREED+NCHILD-1].chisq ) ;

    iters++ ;

    #ifdef VERBOSE
    print_population( G ) ;
    printf( "\nCHISQ_DIFF :: %e || %e %e \n\n" , chisq_diff , 
	    G[0].chisq , G[NBREED-1].chisq ) ;
    #endif
  }

  // best fit parameter will be in population 1
  if( iters == GAMAX ) {
    printf( "[GA] finished in GAMAX %zu iterations \n" , iters ) ;
  }


  //#ifdef VERBOSE
  printf( "[GA] finished in %zu iterations \n" , iters ) ;
  for( i = 0 ; i < fdesc -> NPARAMS ; i++ ) {
    fdesc -> f.fparams[i] = G[0].g[i] ;
    printf( "FPARAMS_%zu :: %f \n" , i , fdesc -> f.fparams[i] ) ; 
  }
  fdesc -> f.chisq = G[0].chisq ;
  printf( "CHISQ :: %f \n" , G[0].chisq ) ;
  //#endif

  // cleanse the gene pool
  for( i = 0 ; i < Ngen ; i++ ) {
    free( G[i].g ) ;
  }
  free( G ) ;

  // free the fitfunction
  free_ffunction( &f2 , fdesc -> NPARAMS ) ;

  // free the rng
  gsl_rng_free( r ) ;

  return iters ;
}
