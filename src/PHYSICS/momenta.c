/**
   @file momenta.c
   @brief computes lattice momenta
 */
#include "gens.h"

#include "resampled_ops.h"
#include "stats.h"
#include "sort.h"

//#define VERBOSE

//#define INSERTION

// compute a lattice momentum from fourier modes n
double
lattmom( const size_t *Dimensions ,
	 const size_t Nd ,
	 const int *N ,
	 const int Selection )
{
  size_t mu ;
  double psq = 0.0 , tmp ;
  switch( Selection ) {
  case 0 : // PSQMOM
    for( mu = 0 ; mu < Nd ; mu++ ) {
      tmp = 2 * M_PI * N[mu] / (double)Dimensions[mu] ;
      psq += tmp * tmp ;
    }
    return psq ;
  case 1 : // 2sin mom
    for( mu = 0 ; mu < Nd ; mu++ ) {
      tmp = 2 * sin( M_PI * N[mu] / (double)Dimensions[mu] ) ;
      psq += tmp * tmp ;
    }
    return psq ;
  case 2 : // sinmom
    for( mu = 0 ; mu < Nd ; mu++ ) {
      tmp = sin( 2 * M_PI * N[mu] / (double)Dimensions[mu] ) ;
      psq += tmp * tmp ;
    }
    return psq ;
  default :
    fprintf( stderr , "Mom selection %d not known\n" , Selection ) ;
    return sqrt(-1) ;
  }
}

// function for averaging equivalent x-values
int
average_equivalent( struct input_params *Input )
{
  size_t i , j ;
  
  // compute error if it hasn't been done already
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    compute_err( &Input -> Data.x[i] ) ;
    compute_err( &Input -> Data.y[i] ) ;
    #ifdef VERBOSE
    fprintf( stdout , "IN %e %e | %e %e \n" ,
	     Input -> Data.x[i].avg , Input -> Data.x[i].err ,
	     Input -> Data.y[i].avg , Input -> Data.y[i].err ) ;
    #endif
  }

  printf( "SOrting\n" ) ;

#ifdef INSERTION
  // sort the data
  if( insertion_sort_data( Input ) == FAILURE ) {
    return FAILURE ;
  }
#else
  // sort the data
  if( quick_sort_data( Input ) == FAILURE ) {
    return FAILURE ;
  }
#endif
  
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    #ifdef VERBOSE
    fprintf( stdout , "OUT %e %e \n" ,
	     Input -> Data.x[i].avg ,
	     Input -> Data.y[i].avg ) ;
    #endif
  }

  // count number of distinct momenta
  size_t NewNdata[ Input -> Data.Nsim ] , Ntot = 0 , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    NewNdata[i] = 0 ;
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      size_t k = j+1 ;
      while( ( fabs( Input -> Data.x[k].avg - Input -> Data.x[j].avg ) < 1E-12 )
	     &&
	     ( k < ( shift + Input -> Data.Ndata[i] ) ) ) {
	#ifdef VERBOSE
	printf( "%f degenerate with %f (%zu %zu) \n" ,
		Input -> Data.x[k].avg ,
		Input -> Data.x[j].avg ,
		j , k ) ;
	#endif
	k++ ;
      }
      NewNdata[i]++ ;
      j = k-1 ;
    }
    #ifdef VERBOSE
    printf( "New Ndata %zu \n" , NewNdata[i] ) ;
    #endif
    shift = j ;
    Ntot += NewNdata[i] ;
  }
  
  // momentum average into temporary distributions
  struct resampled *tmpx = malloc( Ntot * sizeof( struct resampled ) ) ;
  struct resampled *tmpy = malloc( Ntot * sizeof( struct resampled ) ) ;
  size_t idx = 0 ;
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      // initialise x and y      
      tmpx[idx] = init_dist( &Input -> Data.x[j] ,
			     Input -> Data.x[j].NSAMPLES ,
			     Input -> Data.x[j].restype ) ;
      tmpy[idx] = init_dist( &Input -> Data.y[j] ,
			     Input -> Data.y[j].NSAMPLES ,
			     Input -> Data.y[j].restype ) ;

      // loop degenerate momenta
      size_t k = 1 ;
      while( ( fabs( Input -> Data.x[j+k].avg - Input -> Data.x[j].avg ) < 1E-12 )
	     &&
	     ( k < ( shift + Input -> Data.Ndata[i] ) ) ) {
	add( &tmpx[idx] , Input -> Data.x[j+k] ) ;
	add( &tmpy[idx] , Input -> Data.y[j+k] ) ;
	k++ ;
      }
      fprintf( stdout , "Average on x = %e has %zu self-averages\n" ,
	       Input -> Data.x[j].avg , k ) ;

      mult_constant( &tmpx[idx] , 1/(double)k ) ;
      mult_constant( &tmpy[idx] , 1/(double)k ) ;
      
      idx++ ;
      
      j += (k-1) ;
    }
    shift = j ;
  }

  // free all of x and y 
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( Input -> Data.x[i].resampled ) ;
    free( Input -> Data.y[i].resampled ) ;
  }
  free( Input -> Data.x ) ;
  free( Input -> Data.y ) ;

  // allocate new versions
  Input -> Data.x = malloc( Ntot * sizeof( struct resampled ) ) ;
  Input -> Data.y = malloc( Ntot * sizeof( struct resampled ) ) ;
  
  // reallocate Data
  for( i = 0 ; i < Ntot ; i++ ) {

    Input -> Data.x[i] = init_dist( &tmpx[i] , tmpx[i].NSAMPLES , tmpx[i].restype ) ;
    Input -> Data.y[i] = init_dist( &tmpy[i] , tmpx[i].NSAMPLES , tmpx[i].restype ) ;

    #ifdef VERBOSE
    printf( "Averaged %e %e | %e %e \n" ,
	    Input -> Data.x[i].avg , Input -> Data.x[i].err ,
	    Input -> Data.y[i].avg , Input -> Data.y[i].err ) ;
    #endif
  }

  // finally set the data struct
  Input -> Data.Ntot = Ntot ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    Input -> Data.Ndata[i] = NewNdata[i] ;
  }

  // free tmpx and tmpy
  for( i = 0 ; i < Ntot ; i++ ) {
    free( tmpx[i].resampled ) ;
    free( tmpy[i].resampled ) ;
  }
  free( tmpx ) ;
  free( tmpy ) ;
  
  return SUCCESS ;
}
