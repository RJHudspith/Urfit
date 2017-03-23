/**
   @file distribution.c
   @brief read a distribution
 */
#include "gens.h"

#include "GLU_bswap.h"
#include "momenta.h"
#include "resampled_ops.h"
#include "stats.h"

enum{ DO_NOT_ADD , ADD_TO_LIST } list_creation ;

////////// Cylinder cutting procedurals //////////
// gets the body diagonal vectors for our lattice
static inline void
get_diagonal( ND , n , i , DIMS )
     const int ND ;
     double n[ ND ] ;
     const int i , DIMS ;
{
  int mu , subvol = 1 ;
  for( mu = 0 ; mu < ND ; mu++ ) {
    if( mu < DIMS ) {
      n[ mu ] = ( ( i - i % subvol ) / subvol ) % 2 ;
      if( n[ mu ] == 0 ) {
	n[ mu ] -- ;
      }
      subvol *= 2 ;
    } else {// set it to 0?
      n[ mu ] = 0 ;
    }
  }
  return ;
}

// generic cylinder calc
static int 
cylinder_DF( const int ND ,
	     const double q[ ND ] ,
	     const int DIMS ,
	     const double cyl_width )
{
  // test that this satisfies the correct cylinder cutting procedure
  const double norm = 1. / ( 2. * sqrt( DIMS ) ) ;
  const int diagonals = 2 << ( ND - 1 ) ;
 
  // generix for loop over the diagonals
  int mu ;
  double x[ ND ] ;
  for( mu = 0 ; mu < diagonals ; mu ++ ) {
    // inline for getting a diagonal for lexi order
    get_diagonal( ND , x , mu , DIMS ) ;

    double scalar_prod = 0. ; 
    int nu ;
    for( nu = 0 ; nu < ND ; nu ++ ) {
      scalar_prod += q[ nu ] * x[ nu ] ;
    }
    scalar_prod *= norm ;

    double mod = 0.0 ; 
    for( nu = 0 ; nu < ND ; nu ++ ) {
      register const double temp = q[ nu ] - scalar_prod * x[ nu ] ;
      mod += temp * temp ;
    }

    if( sqrt( mod ) <= cyl_width ) {
      return ADD_TO_LIST ;
    }
  }
  return DO_NOT_ADD ;
}

//
bool **
init_dists( struct input_params *Input )
{
  bool **inlist = malloc( Input -> Data.Nsim * sizeof( bool* ) ) ;
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  Input -> Data.Ntot = 0 ;
  
  size_t i , j , mu , Nsamples[ Input -> Data.Nsim ] ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // open the file
    FILE *file = NULL ;
    file = fopen( Input -> Traj[i].FileY , "rb" ) ;
    if( file == NULL ) {
      return NULL ;
    }
    
    // read the first few numbers
    int tmp[1] ;
    if( fread( tmp , (sizeof( int )) , 1 , file ) != 1 ) return NULL ;
    #if WORDS_BIGENDIAN
    bswap_32( 1 , tmp ) ;
    #endif
    
    uint32_t size[1] ;
    if( fread( size , (sizeof( uint32_t )) , 1 , file ) != 1 ) return NULL ;
    #if WORDS_BIGENDIAN
    bswap_32( 1 , size ) ;
    #endif

    size_t sum = 0 ;
    
    // read and skip the momentum list
    Input -> Data.Ndata[i] = size[0] ;
    inlist[i] = malloc( size[0] * sizeof( bool ) ) ;

    // loop Ndata reading in the momenta
    for( j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
      // read in a momentum
      if( fread( size , (sizeof( int )) , 1 , file ) != 1 ) {
	return NULL ;
      }
      #if WORDS_BIGENDIAN
      bswap_32( 1 , size ) ;
      #endif

      int mom[ size[0] ] ;
      if( fread( mom , (sizeof( int )) , size[0] , file ) != size[0] ) {
	return NULL ;
      }
      #if WORDS_BIGENDIAN
      bswap_32( size[0] , mom ) ;
      #endif

      // does it pass the filter?
      double q[ size[0] ] ;
      for( mu = 0 ; mu < size[0] ; mu++ ) {
	q[mu] = 2 * M_PI * mom[mu] / Input -> Traj[i].Dimensions[mu] ;
      }

      if( cylinder_DF( size[0] , q , size[0] , 0.24 ) == ADD_TO_LIST ) {
	sum++ ;
	inlist[i][j] = true ;
      } else {
	inlist[i][j] = false ;
      }
    }

    Input -> Data.Ndata[i] = sum ;
    
    // read the number of samples
    if( fread( size , (sizeof( uint32_t )) , 1 , file ) != 1 ) {
      return NULL ;
    }
    #if WORDS_BIGENDIAN
    bswap_32( 1 , size ) ;
    #endif

    Nsamples[i] = size[0] ;

    Input -> Data.Ntot += Input -> Data.Ndata[i] ;

    fclose( file ) ;
  }

  // allocate all of the resampled stuff
  Input -> Data.x = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;
  Input -> Data.y = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;
  size_t shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      Input -> Data.x[j].resampled = malloc( Nsamples[i] * sizeof( double ) ) ;
      Input -> Data.x[j].NSAMPLES  = Nsamples[i] ;
      Input -> Data.x[j].restype   = Raw ;
      
      Input -> Data.y[j].resampled = malloc( Nsamples[i] * sizeof( double ) ) ;
      Input -> Data.y[j].NSAMPLES  = Nsamples[i] ;
      Input -> Data.y[j].restype   = Raw ;
    }
    shift += Input -> Data.Ndata[i] ;
  }

  return inlist ;
}

// read a distribution
int
read_distribution_old( struct input_params *Input )
{
  // sanity check filenames
  size_t i , j , k , idx = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    if( Input -> Traj[i].FileX != NULL ) {
      fprintf( stderr , "[IO] read_distribution_old only reads Y data \n" ) ;
      return FAILURE ;
    }
    if( Input -> Traj[i].FileY == NULL ) {
      fprintf( stderr , "[IO] traj_%zu entry has NULL Y value \n" , i ) ;
      return FAILURE ;
    }
  }

  // initialise the Input struct
  bool **inlist = NULL ;
  if( ( inlist = init_dists( Input ) ) == NULL ) {
    fprintf( stderr , "Failed to init Input\n" ) ;
    return FAILURE ;
  }

  //const double ainv[ 3 ] = { 3.148 , 2.3833 , 1.7848 } ;
  //const double ainv[ 3 ] = { 1.7848 , 2.3833 , 1.7848 } ;
  //const double ainv[ 3 ] = { 2.3833 , 2.3833 , 2.3833 } ;
  //const double ainv[ 3 ] = { 1.7848 , 1.7848 , 1.7848 } ;

  const double ainv[ 3*12 ] = { 3.148 , 3.148 , 3.148 , 3.148 ,
				3.148 , 3.148 , 3.148 , 3.148 ,
				3.148 , 3.148 , 3.148 , 3.148 ,
				2.3833 , 2.3833 , 2.3833 , 2.3833 ,
				2.3833 , 2.3833 , 2.3833 , 2.3833 ,
				2.3833 , 2.3833 , 2.3833 , 2.3833 ,
				1.7848 , 1.7848 , 1.7848 , 1.7848 ,
				1.7848 , 1.7848 , 1.7848 , 1.7848 ,
				1.7848 , 1.7848 , 1.7848 , 1.7848 } ;
  
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // open the file
    FILE *file = NULL ;
    file = fopen( Input -> Traj[i].FileY , "rb" ) ;
    if( file == NULL ) {
      return FAILURE ;
    }

    // read the first few numbers
    int tmp[1] ;
    if( fread( tmp , (sizeof( int )) , 1 , file ) != 1 ) return FAILURE ;
    #if WORDS_BIGENDIAN
    bswap_32( 1 , tmp ) ;
    #endif
    
    uint32_t size[1] ;
    if( fread( size , (sizeof( uint32_t )) , 1 , file ) != 1 ) return FAILURE ;
    #if WORDS_BIGENDIAN
    bswap_32( 1 , size ) ;
    #endif

    const size_t Nmoms = size[0] , idx_set = idx ;

    // loop Ndata reading in the momenta
    for( j = 0 ; j < Nmoms ; j++ ) {
      // read in a momentum
      if( fread( size , (sizeof( int )) , 1 , file ) != 1 ) {
	return FAILURE ;
      }
      #if WORDS_BIGENDIAN
      bswap_32( 1 , size ) ;
      #endif

      int mom[ size[0] ] ;
      if( fread( mom , (sizeof( int )) , size[0] , file ) != size[0] ) {
	return FAILURE ;
      }
      #if WORDS_BIGENDIAN
      bswap_32( size[0] , mom ) ;
      #endif

      #if 0
      printf( "here %u \n" , size[0] ) ;
      printf( "%zu %zu %zu %zu \n" ,
	      Input -> Traj[i].Dimensions[0] ,
	      Input -> Traj[i].Dimensions[1] ,
	      Input -> Traj[i].Dimensions[2] ,
	      Input -> Traj[i].Dimensions[3] ) ;
      
      printf( "%d %d %d %d -> %f \n" , mom[0] , mom[1] , mom[2] , mom[3] , 
	      lattmom( Input -> Traj[i].Dimensions , size[0] , mom , 1 ) ) ;
      #endif

      if( inlist[i][j] == true ) {
	const double p2 = lattmom( Input -> Traj[i].Dimensions ,
				   size[0] , mom , 1 ) * ainv[i] * ainv[i] ;
      
	// compute the momenta
	equate_constant( &(Input -> Data.x[idx]) , p2 ,
			 Input -> Data.x[idx].NSAMPLES , Raw ) ;
	
	idx++ ;
      }
    }
    
    // fill in the distribution
    idx = idx_set ;
    
    for( j = 0 ; j < Nmoms ; j++ ) {

      // read the number of samples
      if( fread( size , (sizeof( uint32_t )) , 1 , file ) != 1 ) {
	return FAILURE ;
      }
      #if WORDS_BIGENDIAN
      bswap_32( 1 , size ) ;
      #endif
      
      // fill the resamples
      double tmp[ size[0] ] ;
      if( fread( tmp , sizeof( double ) , size[0] , file ) != size[0] ) {
	return FAILURE ;
      }
      #if WORDS_BIGENDIAN
      bswap_32( size[0] , tmp ) ;
      #endif
      // loop number of samples
      if( inlist[i][j] == true ) {
	for( k = 0 ; k < size[0] ; k++ ) {
	  Input -> Data.y[idx].resampled[k] = tmp[k] ;
	}
	compute_err( &Input -> Data.y[idx] ) ;
	idx++ ;
      }
    }
    
    fclose( file ) ;
  }

  // sort and average the rest
  average_equivalent( Input ) ;
    
  return SUCCESS ;
}
