/**
   @file read_corr.c
   @brief read a CORR output file
 */
#include "gens.h"

#include "GLU_bswap.h"
#include "resampled_ops.h"
#include "tfold.h"

// do we have to byte swap?
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
read_magic_gammas( FILE *file , 
		   uint32_t *NGSRC ,
		   uint32_t *NGSNK ,
		   uint32_t *LT ,
		   uint32_t *mommatch ,
		   const uint32_t *mompoint )
{
  uint32_t magic[ 1 ] = { } , NMOM[ 1 ] , n[ 4 ] ;
  size_t p ;
  
  if( fread( magic , sizeof( uint32_t ) , 1 , file ) != 1 ) {
    printf( "[IO] magic read failure \n" ) ;
    return FAILURE ;
  }
  // check the magic number, tells us the edianness
  if( magic[0] != 67798233 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 67798233 ) {
      printf( "Magic number read failure\n" ) ;
      return FAILURE ;
    }
    must_swap = true ;
  }

  // read the momentum list
  FREAD32( NMOM , sizeof( uint32_t ) , 1 , file ) ;

  mommatch = 0 ;
  for( p = 0 ; p < NMOM[0] ; p++ ) {
    FREAD32( n , sizeof( uint32_t ) , 4 , file ) ;
    if( n[ 0 ] != 4-1 ) {
      printf( "[MOMLIST] %d should be %d \n" , n[ 0 ] , 3 ) ;
      return FAILURE ;
    }
    size_t mu , matches = 0 ;
    for( mu = 0 ; mu < 3 ; mu++ ) {
      if( n[mu] == mompoint[mu] ) matches++ ;
      if( matches == 3 ) *mommatch = p ;
    }
  }

  FREAD32( NMOM , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( NGSRC , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( NGSNK , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( LT , sizeof( uint32_t ) , 1 , file ) ;

  return SUCCESS ;
}

// read correlation files
int
pre_allocate( struct input_params *Input ,
	      uint32_t *mommatch )
{
  size_t i ;
  uint32_t mompoint[ 3 ] = { 0 , 0 , 0 } ;
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  Input -> Data.Ntot = 0 ;
  // loop number of trajectories reading in the header information from the beginning ones
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // open a file
    char str[ Input -> Traj[i].Filename_Length + 6 ] ;
    sprintf( str , Input -> Traj[i].Filename , Input -> Traj[i].Begin ) ;
    FILE *file = fopen( str , "rb" ) ;
    if( file == NULL ) {
      fprintf( stderr , "[IO] cannot open %s \n" , str ) ;
      return FAILURE ;
    }

    // read a file to figure out how long it is
    uint32_t Ngsrc , Ngsnk , LT ;
    if( read_magic_gammas( file , &Ngsrc , &Ngsnk , &LT , mommatch , mompoint ) == FAILURE ) {
      return FAILURE ;
    }

    // sanity check LT
    if( (size_t)LT != Input -> Traj[i].Dimensions[ Input -> Traj[i].Nd - 1 ] ) {
      fprintf( stderr , "[IO] File LT and input file LT do not match!\n" ) ;
      return FAILURE ;
    }
    Input -> Data.Ndata[i] = LT ;
    Input -> Data.Ntot += LT ;

    fclose( file ) ;
  }

  // now set Ndata to be LT and allocate x and y
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
      Input -> Data.y[j].NSAMPLES = Nmeas ;
      Input -> Data.y[j].restype = Raw ;
      
      // set x to "t"
      equate_constant( &( Input -> Data.x[j] ) , t , Nmeas , Raw ) ;
      t++ ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  
  return SUCCESS ;
}

// read in the correlators into our data struct
int
read_corr( struct input_params *Input ,
	   const foldtype fold ,
	   const size_t src ,
	   const size_t snk )
{
  // go through the files and allocate x and y
  uint32_t mompoint[ 3 ] = { 0 , 0 , 0 } , mommatch = 0 ;
  if( pre_allocate( Input , &mommatch ) == FAILURE ) {
    return FAILURE ;
  }
  
  // reread files and poke in the data
  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // temporary string
    char str[ Input -> Traj[i].Filename_Length + 6 ] ;

    // linearised measurement index
    register size_t meas = 0 ;
    for( k = Input -> Traj[i].Begin ;
	 k < Input -> Traj[i].End ;
	 k += Input -> Traj[i].Increment ) {
      
      // open up a file
      sprintf( str , Input -> Traj[i].Filename , k ) ;
      FILE *file = fopen( str , "rb" ) ;
      if( file == NULL ) {
	fprintf( stderr , "[IO] cannot open %s \n" , str ) ;
	return FAILURE ;
      }

      // read in the correlator data and poke into "y"
      uint32_t NGSRC , NGSNK, LT ;
      if( read_magic_gammas( file , &NGSRC , &NGSNK , &LT ,
			     &mommatch , mompoint ) == FAILURE ) {
	return FAILURE ;
      }
	
      // mommatch is at the zero point for the moment
      const size_t Gseek = (size_t)( snk + NGSRC * src ) * mommatch ;
      // skip all of the correlators that are in the file
      const size_t Cseek = Gseek * (size_t)( LT * sizeof( double complex ) ) ;
      // skip all of the LT's that are in the file up to the one we want to check
      const size_t Lseek = (Gseek-1) * (size_t)( sizeof( uint32_t ) ) ;
      // skip the file to the correlator we need
      fseek( file , Lseek + Cseek , SEEK_CUR ) ;

      // read the initial LT and make sure that it is the same as Ndata
      if( FREAD32( &LT , sizeof( uint32_t ) , 1 , file ) == FAILURE ) {
	fprintf( stderr , "[IO] Fread failure (LT)\n" ) ;
	return FAILURE ;
      }
      if( LT != Input -> Data.Ndata[i] ) {
	fprintf( stderr , "[IO] LT mismatch (read %d) (Ndata %zu)\n" ,
		 LT , Input -> Data.Ndata[i] ) ;
	return FAILURE ;
      }
	
      double complex C[ LT ] ;
      if( FREAD64( C , sizeof( double complex ) , LT , file ) == FAILURE ) {
	fprintf( stderr , "[IO] Fread failure C(t) \n" ) ;
	return FAILURE ;
      }

      // otherwise poke C into y
      for( j = 0 ; j < LT ; j++ ) {
	Input -> Data.y[ shift + j ].resampled[ meas ] = -creal( C[j] ) ;
      }

      meas ++ ;
      fclose( file ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  
  return SUCCESS ;
}

