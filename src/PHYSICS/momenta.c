/**
   @file momenta.c
   @brief computes lattice momenta
 */
#include "gens.h"

#include "resampled_ops.h"
#include "stats.h"

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

// sort the data by x.avg
int
sort_data( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
    
  // sort in ascending order the data by x-avg using an insertion sort
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    struct resampled tempx = init_dist( &Input -> Data.x[ shift ] ,
					Input -> Data.x[ shift ].NSAMPLES ,
					Input -> Data.x[ shift ].restype ) ;
    struct resampled tempy = init_dist( &Input -> Data.y[ shift ] ,
					Input -> Data.y[ shift ].NSAMPLES ,
					Input -> Data.y[ shift ].restype ) ;
    
    // insertion sort
    for( j = shift + 1 ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      double x = Input -> Data.x[ j ].avg ;
      equate( &tempx , Input -> Data.x[j] ) ;
      equate( &tempy , Input -> Data.y[j] ) ;
      
      int hole = (int)j - 1 ;
      while( hole >= shift && Input -> Data.x[hole].avg > x ) {
	// copy data
	equate( &Input -> Data.x[hole+1] , Input -> Data.x[hole] ) ;
	equate( &Input -> Data.y[hole+1] , Input -> Data.y[hole] ) ;
	hole-- ;
      }
      equate( &Input -> Data.x[hole+1] , tempx ) ;
      equate( &Input -> Data.y[hole+1] , tempy ) ;
      #ifdef VERBOSE
      fprintf( stdout , "Sorting :: %zu %zu \n" , j , Input -> Data.Ndata[i] ) ;
      #endif
    }
    shift = j ;
    
    free( tempx.resampled ) ;
    free( tempy.resampled ) ;
  }

#ifdef VERBOSE
  // test
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      fprintf( stdout , "TEST -> %f %f %f %f \n" ,
	       Input -> Data.x[j].avg ,
	       Input -> Data.x[j].err_lo ,
	       Input -> Data.x[j].err_hi ,
	       Input -> Data.y[j].avg ) ;
    }
    shift = j ;
  }
#endif
  
  return SUCCESS ;
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
    fprintf( stdout , "IN %f %f \n" ,
	     Input -> Data.x[i].avg ,
	     Input -> Data.y[i].avg ) ;
    #endif
  }

  // sort the data
  if( sort_data( Input ) == FAILURE ) {
    return FAILURE ;
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
		Input -> Data.x[k].err_lo ,
		Input -> Data.x[j].err_hi ,
		j , k ) ;
	#endif
	k++ ;
      }
      NewNdata[i]++ ;
      j = k ;
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
      mult_constant( &tmpx[idx] , 1/(double)k ) ;
      mult_constant( &tmpy[idx] , 1/(double)k ) ;
      
      idx++ ;
      
      j += k ;
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
    printf( "Averaged %f %f | %f %f \n" ,
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
