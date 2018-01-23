/**
   @file read_corr.c
   @brief read a CORR output file
 */
#include "gens.h"

#include <stdint.h>

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
  size_t p , flag ;
  
  if( ( flag = fread( magic , sizeof( uint32_t ) , 1 , file ) ) != 1 ) {
    fprintf( stderr , "[IO] magic read failure %u %zu \n" , magic[0] , flag ) ;
    return FAILURE ;
  }
  // check the magic number, tells us the edianness
  if( magic[0] != 67798233 ) {
    bswap_32( 1 , magic ) ;
    if( magic[0] != 67798233 ) {
      fprintf( stderr , "[IO] Magic number read failure %u \n" , magic[0] ) ;
      return FAILURE ;
    }
    must_swap = true ;
  }

  // read the momentum list
  FREAD32( NMOM , sizeof( uint32_t ) , 1 , file ) ;

  *mommatch = 0 ;
  for( p = 0 ; p < NMOM[0] ; p++ ) {
    FREAD32( n , sizeof( uint32_t ) , 4 , file ) ;
    if( n[ 0 ] != 4-1 ) {
      fprintf( stderr , "[IO] momlist :: %d should be %d \n" , n[ 0 ] , 3 ) ;
      return FAILURE ;
    };
    size_t mu , matches = 0 ;
    for( mu = 1 ; mu < 4 ; mu++ ) {
      if( n[mu] == mompoint[mu-1] ) matches++ ;
      if( matches == 3 ) *mommatch = (uint32_t)p ;
    }
    #ifdef VERBOSE
    fprintf( stdout , "%u %u %u %u\n" , n[0] , n[1] , n[2] , n[3] ) ;
    #endif
  }

  FREAD32( NMOM , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( NGSRC , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( NGSNK , sizeof( uint32_t ) , 1 , file ) ;
  FREAD32( LT , sizeof( uint32_t ) , 1 , file ) ;

  #ifdef VERBOSE
  printf( "%u %u %u %u\n" , *NMOM , *NGSRC , *NGSNK , *LT ) ;
  #endif

  return SUCCESS ;
}

// read correlation files
static int
pre_allocate( struct input_params *Input ,
	      const uint32_t *mompoint )
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
    //printf( "file :: %s \n" , str ) ;

    // read a file to figure out how long it is
    uint32_t Ngsrc , Ngsnk , LT ;
    if( read_magic_gammas( file , &Ngsrc , &Ngsnk , &LT ,
			   &mommatch , mompoint ) == FAILURE ) {
      fprintf( stderr , "[IO] magic gammas failure\n" ) ;
      return FAILURE ;
    }

    //printf( "WHAT? %zu %zu %zu \n" , Ngsrc , Ngsnk , LT ) ;

    // sanity check LT
    if( (size_t)(2*LT) == Input -> Traj[i].Dimensions[ Input -> Traj[i].Nd - 1 ] ) {
      fprintf( stdout , "[IO] File LT is twice input time length\n" ) ;
    } else if( (size_t)LT != Input -> Traj[i].Dimensions[ Input -> Traj[i].Nd - 1 ] ) {
      fprintf( stderr , "[IO] File LT and input file LT do not match!\n" ) ;
      return FAILURE ;
    }

    // if we are folding the data this is halved
    if( Input -> Traj[i].Fold != NOFOLD &&
	Input -> Traj[i].Fold != NOFOLD_MINUS ) {
      Input -> Data.Ndata[i] = LT/2 ;
      Input -> Data.Ntot += LT/2 ;
    } else {
      Input -> Data.Ndata[i] = LT ;
      Input -> Data.Ntot += LT ;
    }
    
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

// get the correlation function
int
get_correlator( double complex *C ,
		const char *str ,
		const size_t snk ,
		const size_t src ,
		const uint32_t *mompoint ,
		const size_t Nlt )
{
  uint32_t mommatch = 0 ;
  
  // open the file
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
  if( snk >= (size_t)NGSNK || src >= (size_t)NGSRC ) {
    fprintf( stderr , "[IO] source and sink provided are out of bounds\n"
	     "[IO] ( %zu , %zu ) >= ( %du , %du ) \n" ,
	     src , snk , NGSRC , NGSNK ) ;
    return FAILURE ;
  }

  double complex *Ctmp = calloc( LT , sizeof( double complex ) ) ;
  size_t i ;

  // seeking code
  const size_t goffset = snk + src * (size_t)NGSNK ;
  const size_t OFFSET = goffset * ( Nlt * sizeof( double complex ) + sizeof( uint32_t ) ) ;
  
  if( fseek( file , OFFSET , SEEK_CUR ) != 0 ) {
    fprintf( stderr , "[IO] Fseek failed\n" ) ;
  }

  if( FREAD64( Ctmp , sizeof( double complex ) , LT , file ) == FAILURE ) {
    fprintf( stderr , "[IO] Fread failure C(t) \n" ) ;
    return FAILURE ;
  }
  
  // read the final LT and make sure that it is the same as Ndata
  if( src != (NGSRC-1) && snk != (NGSNK-1) ) {
    if( FREAD32( &LT , sizeof( uint32_t ) , 1 , file ) == FAILURE ) {
      fprintf( stderr , "[IO] Fread failure (LT) %u \n" , LT ) ;
      return FAILURE ;
    }
  }
  if( LT != Nlt ) {
    fprintf( stderr , "[IO] LT mismatch (read %u) (Ndata %zu)\n" ,
	     LT , Nlt ) ;
    return FAILURE ;
  }

  
 corr :

  // sum into C
  for( i = 0 ; i < LT ; i++ ) {
    C[ i ] += Ctmp[ i ] ;
  }
  fclose( file ) ;

  free( Ctmp ) ;
  
  return SUCCESS ;
}

// read in the correlators into our data struct
int
read_corr( struct input_params *Input ,
	   const fold_type fold ,
	   const size_t src ,
	   const size_t snk )
{
  // sanity check filenames
  size_t i , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    if( Input -> Traj[i].FileX != NULL ) {
      fprintf( stderr , "[IO] read_corr only reads Y data \n" ) ;
      return FAILURE ;
    }
    if( Input -> Traj[i].FileY == NULL ) {
      fprintf( stderr , "[IO] traj_%zu entry has NULL Y value \n" , i ) ;
      return FAILURE ;
    }
  }
  
  // go through the files and allocate x and y
  uint32_t *mompoint = malloc( 3 * sizeof( uint32_t ) ) ;
  mompoint[0] = 0 ; mompoint[1] = 0 ; mompoint[2] = 0 ;
  
  if( pre_allocate( Input , mompoint ) == FAILURE ) {
    free( mompoint ) ;
    return FAILURE ;
  }

  // reread files and poke in the data
  int Flag = SUCCESS ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // temporary string
    char str[ strlen( Input -> Traj[i].FileY ) + 6 ] ;

    // set the temporary correlator
    size_t Nlt = Input -> Data.Ndata[i] ;
    if( Input -> Traj[i].Fold != NOFOLD &&
	Input -> Traj[i].Fold != NOFOLD_MINUS ) {
      Nlt *= 2 ;
    }
    
    // linearised measurement index
    register size_t meas = 0 ;
    for( k = Input -> Traj[i].Begin ;
	 k < Input -> Traj[i].End ;
	 k += Input -> Traj[i].Increment ) {
      
      // open up a file
      sprintf( str , Input -> Traj[i].FileY , k ) ;

      double complex *C = NULL ;
      
      // apply the correct summation map from the correlator data
      if( ( C = map_correlator( Input -> Traj[i] , str ,
				mompoint , Nlt ) ) == NULL ) {
	Flag = FAILURE ;
	break ;
      }

      // poke C correctly into y, doing the folding if required
      if( time_fold( Input -> Data.y + shift , C , Nlt ,
		     Input -> Traj[i].Fold , meas ) == FAILURE ) {
	Flag = FAILURE ;
	break ;
      }

      fprintf( stdout , "[IO] measurement %zu done\n" , meas ) ;

      free( C ) ;
      meas ++ ;
    }    
    shift += Input -> Data.Ndata[i] ;
  }
  
  free( mompoint ) ;
  				     
  return Flag ;
}

