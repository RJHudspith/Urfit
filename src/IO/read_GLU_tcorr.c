/**
   @file read_GLU_tcorr.c
   @brief read a temporal correlator produced by GLU
 */
#include "gens.h"

#include "GLU_bswap.h"

static int
init_GLU_tcorr( struct input_params *Input )
{
  // loop trajectories
  size_t i , j ;
  Input -> Data.Ntot = 0 ;
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    Input -> Data.Ndata[i] = Input -> Traj[i].Dimensions[ 3 ] ;
    Input -> Data.Ntot += Input -> Data.Ndata[i] ;
  }

  Input -> Data.x = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;
  Input -> Data.y = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;

  size_t shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // allocate Nsamples
    const size_t Nsamples = ( Input -> Traj[i].End -
			      Input -> Traj[i].Begin ) / Input -> Traj[i].Increment ;

    printf( "%zu %zu \n" , Nsamples , Input -> Data.Ntot ) ;
    
    for( j = 0 ; j < Input -> Traj[i].Dimensions[ 3 ] ; j++ ) {
      Input -> Data.x[ shift ].resampled = malloc( Nsamples * sizeof( double ) ) ;
      Input -> Data.x[ shift ].NSAMPLES = Nsamples ;
      Input -> Data.x[ shift ].restype = Raw ;
      
      Input -> Data.y[ shift ].resampled = malloc( Nsamples * sizeof( double ) ) ;
      Input -> Data.y[ shift ].NSAMPLES = Nsamples ;
      Input -> Data.y[ shift ].restype = Raw ;

      shift++ ;
    }
    printf( "SHIFT %zu \n" , shift ) ;
  }
  return SUCCESS ;
}

int
read_GLU_tcorr( struct input_params *Input )
{
  char str[ 256 ] ;
  
  if( init_GLU_tcorr( Input ) == FAILURE ) {
    return FAILURE ;
  }
  
  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // loop files
    size_t idx = 0 ;
    for( j = Input -> Traj[i].Begin ;
	 j < Input -> Traj[i].End ;
	 j += Input -> Traj[i].Increment ) {

      // print in the trajectory index
      sprintf( str , Input -> Traj[i].FileY , j ) ;

      printf( "%s \n" , str ) ;

      // open the file
      FILE *Infile = fopen( str , "rb" ) ;
      if( Infile == NULL ) {
	return FAILURE ;
      }

      // Lt is the length of the time correlator
      uint32_t Lt[ 1 ] ;
      fread( Lt , sizeof( uint32_t ) , 1 , Infile ) ;
      #ifndef WORDS_BIGENDIAN
      bswap_32( 1 , Lt ) ;
      #endif

      // simple check
      if( Lt[0] != Input -> Traj[i].Dimensions[ 3 ] ) {
	fprintf( stderr , "[IO] lt and input file differ %u %zu \n" ,
		 Lt[0] , Input -> Traj[i].Dimensions[ 3 ] ) ;
	return FAILURE ;
      }

      // t-list
      uint32_t tlist[ Lt[0] ] ;
      if( fread( tlist , sizeof(uint32_t) , Lt[0] , Infile ) != Lt[0] ) {
	fprintf( stderr , "[IO] Tlist read failure \n" ) ;
        return FAILURE ;
      }
      #ifndef WORDS_BIGENDIAN
      bswap_32( Lt[0] , tlist ) ;
      #endif

      uint32_t NewLt[ 1 ] ;
      fread( NewLt , sizeof( uint32_t ) , 1 , Infile ) ;
      #ifndef WORDS_BIGENDIAN
      bswap_32( 1 , NewLt ) ;
      #endif

      if( NewLt[0] != Lt[0] ) {
	fprintf( stderr , "[IO] Lt[0] misread \n" ) ;
	return FAILURE ;
      }

      // read in y data
      double Ct[ Lt[0] ] ;
      fread( Ct , sizeof( double ) , Lt[0] , Infile ) ;
      #ifndef WORDS_BIGENDIAN
      bswap_64( Lt[0] , Ct ) ;
      #endif
      for( k = 0 ; k < Lt[0] ; k++ ) {
	printf( "Corr %e \n" , creal( Ct[k] ) ) ;
      }

      // set x and y data
      for( k = 0 ; k < Lt[0] ; k++ ) {
      Input -> Data.x[ shift + k ].resampled[idx] = k ;
      Input -> Data.y[ shift + k ].resampled[idx] = Ct[k] ;
      printf( "Check :: %f %f \n" ,
	Input -> Data.x[ shift + k ].resampled[idx] ,
	Input -> Data.y[ shift + k ].resampled[idx] ) ;
      }

      idx++ ;

      printf( "Closing file \n" ) ;
      
      fclose( Infile ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }

  return SUCCESS ;
}