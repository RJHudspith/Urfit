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
  if( must_swap ) bswap32( N , d ) ;
  return SUCCESS ;
}

// read 64 bytes
static int
FREAD64( void *d , const size_t size , const size_t N , FILE *file )
{
  if( fread( d , size , N , file ) != N ) {
    return FAILURE ;
  }
  if( must_swap ) bswap64( N , d ) ;
  return SUCCESS ;
}

// read the gamma matrices
static int
read_magic_gammas( FILE *file , 
		   uint32_t *NGSRC ,
		   uint32_t *NGSNK ,
		   uint32_t *LT )
{
  uint32_t magic[ 1 ] = { } , NMOM[ 1 ] , n[ 4 ] ;
  size_t p , mu ;
  
  if( fread( magic , sizeof( uint32_t ) , 1 , file ) != 1 ) {
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

  for( p = 0 ; p < NMOM[0] ; p++ ) {
    FREAD( n , sizeof( uint32_t ) , 4 , file ) ;
    if( n[ 0 ] != 4-1 ) {
      printf( "[MOMLIST] %d should be %d \n" , n[ 0 ] , 3 ) ;
      return FAILURE ;
    }
  }

  FREAD32( NMOM , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( NGSRC , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( NGSNK , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( LT , sizeof( uint32_t ) , 1 , file ) ;

  return SUCCESS ;
}

// read a single file
static int
read_corrfile( struct resampled *sample ,
	       FILE *file ,
	       const size_t meas ,
	       const size_t src ,
	       const size_t snk ,
	       const foldtype fold )
{
  uint32_t NGSRC , NGSNK , LT ;
  if( read_magic_gammas( file , &NGSRC , &NGSNK , &LT ) == FAILURE ) {
    return FAILURE ;
  }
  const size_t gseek = (size_t)( snk + NGSNK * src ) ;
  const size_t cseek = gseek * ( LT * sizeof( double complex ) ) ;
  const size_t Lseek = ( gseek - 1 ) * ( LT * sizeof( uint32_t ) ) ;
  fseek( file , Lseek + cseek , SEEK_CUR ) ;

  FREAD32( &LT , sizeof( uint32_t ) , 1 , file ) ;

  double complex C[ LT ] ;
  FREAD64( C , sizeof( double complex ) , LT , file ) ;

  // perform some folding
  tfold( sample , C , LT , fold , meas ) ;
  
  return SUCCESS ;
}

// read correlation files
struct resampled *
read_rawcorr( struct input_params *INPARAMS ,
	      const char *filename ,
	      const size_t fileno ,
	      const foldtype fold ,
	      const size_t src ,
	      const size_t snk ,
	      const size_t nfile )
{
  // number of measurements
  const size_t Nmeas = ( INPARAMS -> traj_end[ fileno ] - 
			 INPARAMS -> traj_beg[ fileno ] ) / 
    INPARAMS -> traj_inc[fileno ] ;

  char filestr[ 256 ] ;
  sprintf( filestr , filename , INPARAMS -> traj_beg[fileno] ) ;

  // open the initial file
  FILE *file = fopen( filestr , "rb" ) ;
  if( file == NULL ) {
    fprintf( stdout , "[IO] read rawcorr -> cannot open %s\n" ,
	     filestr ) ;
    return NULL ;
  }

  // read in the magic gammas
  uint32_t NGSRC , NGSNK , LT ;
  if( read_magic_gammas( file , &NGSRC , &NGSNK , &LT ) == FAILURE ) {
    return NULL ;
  }

  // set how long the data is
  if( fold == NOFOLD ) {
    INPARAMS -> NDATA[ nfile ] = (size_t)LT ;
  } else {
    INPARAMS -> NDATA[ nfile ] = (size_t)LT/2 ;
  }
  
  fclose( file ) ;

  // read in all the 
  struct resampled *sample = malloc( INPARAMS -> NDATA[ nfile ] * 
				     sizeof( struct resampled ) ) ;
  size_t i ;
  for( i = 0 ; i < INPARAMS -> NDATA[ nfile ] ; i++ ) {
    sample[i].resampled = malloc( Nmeas * sizeof( double ) ) ;
    sample[i].restype   = RAWDATA ;
    sample[i].NSAMPLES  = Nmeas ;
  }

  
  // loop the files
  size_t meas = 0 ;
  for( i =  INPARAMS -> traj_beg[nfile] ; 
       i <  INPARAMS -> traj_end[nfile] ;
       i += INPARAMS -> traj_inc[nfile] ) {
    char loc_filestr[ 256 ] ;
    sprintf( loc_filestr , filename , i ) ;
    FILE *loc_file = fopen( loc_filestr , "rb" ) ;
    if( loc_file == NULL ) {
      fprintf( stdout , "[IO] read_rawcorr cannot open -> %s \n" ,
	       loc_filestr ) ;
      return NULL ;
    }

    if( read_corrfile( sample , loc_file , meas , src , snk , fold ) 
	== FAILURE ) {
      return NULL ;
    }

    meas++ ; // increment the measurement index
  }
  
  return sample ;
}

// read s single correlator
int
read_corr( struct resampled **x ,
	   struct resampled **y ,
	   struct input_params *INPARAMS ,
	   const char *filename ,
	   const foldtype fold ,
	   const size_t src ,
	   const size_t snk ,
	   const size_t fileno )
{
  *y = read_rawcorr( INPARAMS ,filename , fileno , fold , 
		     src , snk , fileno ) ;
  if( *y == NULL ) {
    return FAILURE ;
  }
  size_t i ;
  for( i = 0 ; i < INPARAMS -> NDATA[ fileno ] ; i++ ) {
    init_dist( *y+i , y[i] -> NSAMPLES , y[i] -> restype ) ;
  }
  return SUCCESS ;
}
	   
