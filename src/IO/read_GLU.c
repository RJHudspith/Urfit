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

#define ALLR

//#define VERBOSE

static bool
is_parallel( const int r[4] )
{
#ifndef ALLR
  const int Ncut = 1 ;
  int v[][4] = { //{ 2, 2, 2, 0 } ,
		 //{ 1 ,1 ,1, 3 } ,
    { 1 ,1 ,1, 2 } ,
    //{ 1 ,1 ,1, 1 } 
		 };
  bool isin = false ;
  for( int i = 0 ; i < Ncut ; i++ ) {
    int mu , rdotv = 0 , vsq = 0 , rsq = 0 ;
    for( mu = 0 ; mu < 4 ; mu++ ) {
      vsq += v[i][mu]*v[i][mu] ;
      rsq += r[mu]*r[mu] ;
      rdotv += abs(r[mu])*v[i][mu] ;
    }
    if( fabs( rdotv/(sqrt(rsq)*sqrt(vsq)) - 1. ) < 1E-10 ) {
      isin = true ;
    }
  }
  return isin ;
#else
  return true ;
#endif
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

    // start reading the r-list
    size_t Nfilter = 0 , k ;
    for( k = 0 ; k < rlist[0] ; k++ ) {
      uint32_t Nd[1] ;
      if( fread( Nd , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
	return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_32( 1 , Nd ) ;
      #endif
      int32_t r[ Nd[0] ] ;
      if( fread( r , sizeof( int32_t ) , Nd[0] , Infile ) != Nd[0] ) {
	return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_32( Nd[0] , r ) ;
      #endif
      
      if( is_parallel( r ) ) {
	Nfilter++ ;
	#ifdef VERBOSE
	fprintf( stdout , "MOM %d %d %d %d\n" ,
		 r[0] , r[1] , r[2] , r[3] ) ;
	#endif
      }
    }
    Input -> Data.Ndata[i] = Nfilter ;
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

  fprintf( stdout , "GLU inited\n" ) ;
  
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

      fprintf( stdout , "[IO] reading file %s\n" , str ) ;

      // Lt is the length of the time correlator
      uint32_t rlist[ 1 ] ;
      if( fread( rlist , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
	fprintf( stderr , "[IO] rlist length reaad failure -> %zu\n" , j ) ;
	return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_32( 1 , rlist ) ;
      #endif

      // start reading the r-list
      size_t Nfilter = 0 ;
      bool *filter = malloc( rlist[0]*sizeof(bool) ) ;
      for( k = 0 ; k < rlist[0] ; k++ ) {
	filter[k] = false ;
	size_t l ;
	uint32_t Nd[1] ;
	if( fread( Nd , sizeof( uint32_t ) , 1 , Infile ) != 1 ) {
	  printf( "[IO] Failed to read Nd @ %zu \n" , idx ) ;
	  return FAILURE ;
	}
        #ifndef WORDS_BIGENDIAN
	bswap_32( 1 , Nd ) ;
        #endif
	int32_t r[4] ;
	if( fread( r , sizeof( int32_t ) , Nd[0] , Infile ) != Nd[0] ) {
	  printf( "[IO] Failed to read momentum list @ %zu \n" , idx ) ;
	  return FAILURE ;
	}
        #ifndef WORDS_BIGENDIAN
	bswap_32( Nd[0] , r ) ;
        #endif
	
	if( is_parallel( r ) ) {
	  Input -> Data.x[shift+Nfilter].resampled[idx] = 0.0 ;
	  for( l = 0 ; l < Nd[0] ; l++ ) {
	    Input -> Data.x[shift+Nfilter].resampled[idx] += r[l]*r[l] ;
	  }
	  Nfilter++ ;
	  filter[k] = true ;
	}
	#ifdef VERBOSE
	fprintf( stdout , "MOM %d %d %d %d\n" ,
		 r[0] , r[1] , r[2] , r[3] ) ;
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

      Nfilter = 0 ;
      for( k = 0 ; k < rlist[0] ; k++ ) {

	if( filter[k] == true ) {
	  Input -> Data.y[ shift + Nfilter ].resampled[ idx ] = Cr[k] ;
	  Nfilter++ ;
	}
	
	#ifdef VERBOSE
	fprintf( stdout  , "[IO] data -> %e %e \n" ,
		 Input -> Data.x[ shift + k ].resampled[ idx ] ,
		 Input -> Data.y[ shift + k ].resampled[ idx ] ) ;
	#endif
      }
      idx++ ;

      free( filter ) ;
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

  fprintf( stdout , "Averaging equivalent\n" ) ;
  average_equivalent( Input ) ;
  
  return SUCCESS ;
}
