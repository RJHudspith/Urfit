/**
   @file read_GLU.c
   @brief read a GLU propagator
 */
#include "gens.h"

#include "GLU_bswap.h"
#include "momenta.h"
#include "resampled_ops.h"
#include "sort.h"
#include "stats.h"

static bool
keep_r( const int32_t r[] ,
	const int Nd )
{
  // keep only those that are parallel to 2220
  const int Nr = 2 ;
  const int keep[][4] = { {1,1,1,3} ,
			  {2,2,2,0} } ;
  int i ;
  for( i = 0 ; i < Nr ; i++ ) {
    int mu , a = 0 , b = 0 , ab = 0 ;
    for( mu = 0 ; mu < Nd ; mu++ ) {
      const int ar = abs( (int)r[mu]) ;
      ab += keep[i][mu] * ar ;
      a += ar*ar ;
      b += keep[i][mu]*keep[i][mu] ;
    }
    #ifdef VERBOSE
    printf( "AB %d , A %d , B %d \n" , ab , a , b ) ;
    #endif
    if( fabs( ab/(sqrt(a)*sqrt(b)) - 1. ) < 1E-12 ) {
      #ifdef VERBOSE
      fprintf( stdout , "Kept (%d %d %d %d)\n" ,
	       r[0] , r[1] , r[2] , r[3] ) ;
      #endif
      return true ;
    }
  }
  return true ;
  //return false ;
}

static int
init_GLU( struct input_params *Input )
{
  // loop trajectories
  size_t i , j ;
  Input -> Data.Ntot = 0 ;
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    char str[ 256 ] ;
    sprintf( str , Input -> Traj[i].FileY , Input -> Traj[i].Begin ) ;

    FILE *Infile = fopen( str , "rb" ) ;
    if( Infile == NULL ) {
      fprintf( stderr , "[IO] file %s does not exist!\n" , str ) ;
      return FAILURE ;
    }

    uint32_t rlist[1] ;
    if( fread( rlist , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
      fprintf( stderr , "[IO] cannot read file's first entry\n" ) ;
      return FAILURE ;
    }
    #ifndef WORDS_BIGENDIAN
    bswap_32( 1 , rlist ) ;
    #endif
    #ifdef VERBOSE
    printf( "[IO] file length -> %u \n" , rlist[0] ) ;
    #endif
    Input -> Data.Ndata[i] = rlist[0] ;
    Input -> Data.Ntot += Input -> Data.Ndata[i] ;
    
    fclose( Infile ) ;
  }

  Input -> Data.x = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;
  Input -> Data.y = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;

  size_t shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // allocate Nsamples
    const size_t Nsamples = ( Input -> Traj[i].End -
			      Input -> Traj[i].Begin )
      / Input -> Traj[i].Increment ;
    
    for( j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
      Input -> Data.x[ shift ].resampled = malloc( Nsamples * sizeof( double ) ) ;
      Input -> Data.x[ shift ].NSAMPLES = Nsamples ;
      Input -> Data.x[ shift ].restype = Raw ;
      
      Input -> Data.y[ shift ].resampled = malloc( Nsamples * sizeof( double ) ) ;
      Input -> Data.y[ shift ].NSAMPLES = Nsamples ;
      Input -> Data.y[ shift ].restype = Raw ;

      shift++ ;
    }
  }
  return SUCCESS ;
}

// read the GLU file
int
read_GLU( struct input_params *Input )
{  
  if( init_GLU( Input ) == FAILURE ) {
    fprintf( stderr , "[IO] failed to init GLU_File\n" ) ;
    return FAILURE ;
  }

  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    char str[ 256 ] ;
    
    // loop files
    size_t idx = 0 ;
    for( j = Input -> Traj[i].Begin ;
	 j < Input -> Traj[i].End ;
	 j += Input -> Traj[i].Increment ) {

      // print in the trajectory index
      sprintf( str , Input -> Traj[i].FileY , j ) ;

      // open the file
      FILE *Infile = fopen( str , "rb" ) ;
      if( Infile == NULL ) {
	fprintf( Infile , "[IO] File %s does not exist\n" , str ) ;
	return FAILURE ;
      }

      // Lt is the length of the time correlator
      uint32_t rlist[ 1 ] ;
      if( fread( rlist , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
	fprintf( stderr , "[IO] rlist length reaad failure\n" ) ;
	return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_32( 1 , rlist ) ;
      #endif

      // start reading the r-list
      int32_t r[ rlist[0] ][ 4 ] ;
      for( k = 0 ; k < rlist[0] ; k++ ) {
	r[k][0] = r[k][1] = r[k][2] = r[k][3] = 0 ;
	uint32_t Nd[1] ;
	if( fread( Nd , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
	  printf( "[IO] Failed to read Nd @ %zu \n" , idx ) ;
	  return FAILURE ;
	}
        #ifndef WORDS_BIGENDIAN
	bswap_32( 1 , Nd ) ;
        #endif
	//printf( "WTF %u\n" , Nd[0] ) ;
	if( fread( r[k] , sizeof( int32_t ) , Nd[0] , Infile ) != Nd[0] ) {
	  printf( "[IO] Failed to read momentum list @ %zu \n" , idx ) ;
	  return FAILURE ;
	}
        #ifndef WORDS_BIGENDIAN
	bswap_32( Nd[0] , r[k] ) ;
        #endif
	Input -> Data.x[shift+k].resampled[idx] = 0.0 ;
	#ifdef VERBOSE
	printf( "MOM %d %d %d %d\n" , r[k][0] , r[k][1] , r[k][2] , r[k][3] ) ;
	#endif
      }

      uint32_t Newrlist[ 1 ] ;
      if( fread( Newrlist , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
	fprintf( stderr , "[IO] rlist length read failure\n" ) ;
	return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_32( 1 , Newrlist ) ;
      #endif

      if( Newrlist[0] != rlist[0] ) {
	fprintf( stderr , "[IO] Lt[0] misread \n" ) ;
	return FAILURE ;
      }

      double Cr[ rlist[0] ] ;
      if( fread( Cr , sizeof( double ) , rlist[0] , Infile ) != rlist[0] ) {
	fprintf( stderr , "[IO] cannot read GLU correlator\n" ) ;
	return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_64( rlist[0] , Cr ) ;
      #endif

      size_t kidx = 0 ;
      for( k = 0 ; k < rlist[0] ; k++ ) {
	if( keep_r( r[k] , 4 ) ) {
	  Input -> Data.x[shift+kidx].resampled[idx] = 0.0 ;
	  Input -> Data.y[shift+kidx].resampled[ idx ] = 0 ;
	  for( size_t l = 0 ; l < 4 ; l++ ) {	  
	    Input -> Data.x[shift+kidx].resampled[idx] += r[k][l]*r[k][l] ;
	  }
	  Input -> Data.y[ shift + kidx ].resampled[ idx ] = Cr[k] ;
	  #ifdef VERBOSE
	  printf( "[IO] data -> %e %e \n" ,
		  Input -> Data.x[ shift + kidx ].resampled[ idx ] ,
		  Input -> Data.y[ shift + kidx ].resampled[ idx ] ) ;
          #endif
	  kidx++ ;
	}
      }
      Input -> Data.Ndata[i] = kidx ;
      
      idx++ ;
      
      fclose( Infile ) ;
    }
    
    shift += Input -> Data.Ndata[i] ;
  }

  Input -> Data.Ntot = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    Input -> Data.Ntot += Input -> Data.Ndata[i] ;
  }

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

  printf( "Averaging equivalent\n" ) ;
  average_equivalent( Input ) ;
  
  return SUCCESS ;
}
