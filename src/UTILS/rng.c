/**
   @file rng.c
   @brief random number generator
 */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

static gsl_rng *r ;
static size_t SEED = 0 ;

// free the gsl rng
void
free_rng( void )
{
  gsl_rng_free( r ) ;
  return ;
}

// initialise the rng from urandom
void
init_rng( size_t Seed )
{
  gsl_rng_env_setup( ) ;
  
  const gsl_rng_type *type = gsl_rng_default ;
  r = gsl_rng_alloc( type ) ;

  if( Seed == 0 ) {
    FILE *urandom = fopen( "/dev/urandom" , "r" ) ;
    if( urandom == NULL ) exit( 1 ) ;
    if( fread( &Seed , sizeof( Seed ) , 1 , urandom ) != 1 ) exit(1) ;
    fclose( urandom ) ;
  }

  SEED = Seed ;
  printf( "USING GSL seed %lu \n" , SEED ) ;

  gsl_rng_set( r , Seed ) ;

  return ;
}

// returns a double
double
rng_double( void ) { return gsl_rng_uniform( r ) ; }

// returns a gaussian
double
rng_gaussian( const double sigma ) 
{ 
  return gsl_ran_gaussian( r , sigma ) ; 
}

// return an int
size_t
rng_int( const size_t max_idx ) 
{ 
  return gsl_rng_uniform_int( r , max_idx ) ;
}

// reseed the SEED to what it was before
void
rng_reseed( void ) 
{
  if( SEED != 0 ) {
    gsl_rng_set( r , SEED ) ;
  } else {
    printf( "GSL rng not seeded properly!\n" );
    exit(1) ;
  }
  return ;
}

