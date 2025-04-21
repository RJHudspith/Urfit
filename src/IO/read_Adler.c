/**
   @file read_GLU_tcorr.c
   @brief read a temporal correlator produced by GLU
 */
#include "gens.h"

#include "GLU_bswap.h"
#include "stats.h"
#include "crc32c.h"

static bool must_swap = false ;

// read 32 bytes
static int
FREAD32( void *d , const size_t size , const size_t N , FILE *file )
{
  if( fread( d , size , N , file ) != N ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_32( N , d ) ;
  return SUCCESS ;
}

// read 64 bytes
static int
FREAD64( void *d , const size_t size , const size_t N , FILE *file )
{
  if( fread( d , size , N , file ) != N ) {
    return FAILURE ;
  }
  if( must_swap ) bswap_64( N , d ) ;
  return SUCCESS ;
}

// read the gamma matrices
static int
read_magic_Nmom( FILE *file , 
		 uint32_t *nmom )
{
  uint32_t magic[ 1 ] = { 0 } , NMOM[ 1 ] ;
  int flag ;
  if( ( flag = fread( magic , sizeof( uint32_t ) , 1 , file ) ) != 1 ) {
    fprintf( stderr , "[IO] magic read failure %u %zu \n" , magic[0] , flag ) ;
    return FAILURE ;
  }
  // check the magic number, tells us the edianness
  if( magic[0] != 67678282 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 67678282 ) {
      fprintf( stderr , "[IO] Magic number read failure %u \n" , magic[0] ) ;
      return FAILURE ;
    }
    must_swap = true ;
  }

  // read the momentum list
  FREAD32( nmom , sizeof( uint32_t ) , 1 , file ) ;

  return SUCCESS ;
}

// read correlation files
static int
pre_allocate( struct input_params *Input )
{
  uint32_t mommatch = 0 ;
  size_t i ;
  
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  Input -> Data.Ntot = 0 ;
  // loop number of trajectories reading in the header information from the beginning ones
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // open a file
    char str[ strlen( Input -> Traj[i].FileY ) + 6 ] ;
    sprintf( str , Input -> Traj[i].FileY , Input -> Traj[i].Begin ) ;
    FILE *file = fopen( str , "rb" ) ;
    if( file == NULL ) {
      fprintf( stderr , "[IO] cannot open %s \n" , str ) ;
      return FAILURE ;
    }
    // read a file to figure out how long it is
    uint32_t Nmom;
    if( read_magic_Nmom( file , &Nmom ) == FAILURE ) {
      fprintf( stderr , "[IO] magic gammas failure\n" ) ;
      return FAILURE ;
    }

    Input -> Data.Ndata[i]  = Nmom ;
    Input -> Data.Ntot     += Input -> Data.Ndata[i] ;
    
    fclose(file) ;
  }

  // allocate x and y
  Input -> Data.x = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;
  Input -> Data.y = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;

  size_t shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t j , t = 0 ;
    const size_t Nmeas = ( Input -> Traj[i].End - Input -> Traj[i].Begin ) \
      / Input -> Traj[i].Increment ;
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      // allocate and set x and y
      Input -> Data.x[j].resampled = malloc( Nmeas * sizeof( double ) ) ;
      Input -> Data.x[j].NSAMPLES  = Nmeas ;
      Input -> Data.x[j].restype   = Raw ;
      
      Input -> Data.y[j].resampled = malloc( Nmeas * sizeof( double ) ) ;
      Input -> Data.y[j].NSAMPLES  = Nmeas ;
      Input -> Data.y[j].restype   = Raw ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  
  return SUCCESS ;
}

int
read_Adler( struct input_params *Input )
{
  // sanity check filenames
  size_t i , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    if( Input -> Traj[i].FileX != NULL ) {
      fprintf( stderr , "[IO] read_Adler only reads Y data \n" ) ;
      return FAILURE ;
    }
    if( Input -> Traj[i].FileY == NULL ) {
      fprintf( stderr , "[IO] traj_%zu entry has NULL Y value \n" , i ) ;
      return FAILURE ;
    }
  }

  if( pre_allocate( Input ) == FAILURE ) {
    return FAILURE ;
  }

  fprintf( stdout , "[IO] pre-allocated Adler data\n") ;

  // reread files and poke in the data
  int Flag = SUCCESS ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // linearised measurement index
    register size_t meas = 0 ;
#pragma omp parallel for private(k)
    for( k = Input -> Traj[i].Begin ;
	 k < Input -> Traj[i].End ;
	 k += Input -> Traj[i].Increment ) {

      // temporary string
      char str[ strlen( Input -> Traj[i].FileY ) + 6 ] ;

      const size_t meas = (k-Input -> Traj[i].Begin)/Input -> Traj[i].Increment ;
      
      // open up a file
      sprintf( str , Input -> Traj[i].FileY , k ) ;

      FILE *file = fopen( str , "rb" ) ;

      // read a file to figure out how long it is
      uint32_t Nmom;
      if( read_magic_Nmom( file , &Nmom ) == FAILURE ) {
	fprintf( stderr , "[IO] magic gammas failure\n" ) ;
	Flag = FAILURE ;
      }

      // Read X values
      double pval[ Nmom ] = {} ;
      FREAD64( pval , sizeof( double ) , Nmom , file ) ;
      for( int p = 0 ; p < Nmom ; p++ ) {
	Input -> Data.x[p+shift].resampled[meas] = pval[p] ;
      }

      // should be Nmom again
      uint32_t Nmom_tmp = 0 ;
      FREAD32( &Nmom_tmp , sizeof( uint32_t ) , 1 , file ) ;

      printf( "Here %u == %u\n" , Nmom , Nmom_tmp ) ;
      if( Nmom_tmp != Nmom ) {
	Flag = FAILURE ;
      }
      // these are complex in the file
      double Adl[ 2*Nmom ] = {} ;
      FREAD64( Adl , sizeof( double ) , 2*Nmom , file ) ;
      // just keep the real part
      for( int p = 0 ; p < Nmom ; p++ ) {

	printf( "Here %f %f\n" , pval[p] , Adl[2*p] ) ;
	
	Input -> Data.y[p+shift].resampled[meas] = Adl[2*p] ;
      }
      uint32_t cksuma = 0 , cksumb = 0 ;
      DML_checksum_accum_crc32c( &cksuma , &cksumb , 1 , Adl , 2*sizeof( double )*Nmom ) ;
      uint32_t csum[2] = {} ;
      FREAD32( csum , sizeof( uint32_t ) , 2 , file ) ;
      if( cksuma != csum[0] || cksumb != csum[1] ) {
	fprintf( stderr , "[IO] mismatched checksums %x %x | %x %x\n" ,
		 cksuma , cksumb , csum[0] , csum[1] ) ;
	Flag = FAILURE ;
      }

      
      fclose( file ) ;
    }    
    shift += Input -> Data.Ndata[i] ;
  }

  for( size_t i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    compute_err( &Input -> Data.x[i] ) ;
    compute_err( &Input -> Data.y[i] ) ;
    printf( "%zu %f %f\n" , i , Input -> Data.x[i].avg , Input -> Data.y[i].avg ) ;
  }
  
  return Flag ;
}
