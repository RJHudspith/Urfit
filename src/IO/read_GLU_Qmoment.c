#include "gens.h"
#include "GLU_bswap.h"

// read 32 bytes
static int
FREAD32( void *d , const size_t size , const size_t N , FILE *file )
{
  if( fread( d , size , N , file ) != N ) {
    return FAILURE ;
  }
#ifndef WORDS_BIGENDIAN
  bswap_32( N , d ) ;
#endif
  return SUCCESS ;
}

// read 64 bytes
static int
FREAD64( void *d , const size_t size , const size_t N , FILE *file )
{
  if( fread( d , size , N , file ) != N ) {
    return FAILURE ;
  }
#ifndef WORDS_BIGENDIAN
  bswap_64( N , d ) ;
#endif
  return SUCCESS ;
}

static int
pre_allocate( struct input_params *Input )
{
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  Input -> Data.Ntot = 0 ;
  
  size_t i ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // open a file
    char str[ strlen( Input -> Traj[i].FileY ) + 6 ] ;
    sprintf( str , Input -> Traj[i].FileY , Input -> Traj[i].Begin ) ;
    FILE *file = fopen( str , "rb" ) ;
    if( file == NULL ) {
      fprintf( stderr , "[IO] cannot open %s \n" , str ) ;
      return FAILURE ;
    }
    
    // format has nmeas and then Ncorr
    uint32_t Nmeas[1] , Nmoments[1] ;
    FREAD32( Nmeas , sizeof( uint32_t ) , 1 , file ) ;
    
    if( Input -> Traj[i].Gs > Nmeas[0] ) {
      fprintf( stderr , "[IO] Input file GSRC %zu greater than"
	       " number of measurements %ud\n" ,
	       Input -> Traj[i].Gs , Nmeas[0] ) ;
      return FAILURE ;
    }
    
    FREAD32( Nmoments , sizeof( uint32_t ) , 1 , file ) ;
    
    Input -> Data.Ndata[i] = (size_t)Nmoments[0] ;
    Input -> Data.Ntot += Input -> Data.Ndata[i] ;

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

      equate_constant( &( Input -> Data.x[j] ) , 0.0 , Nmeas , Raw ) ;
      equate_constant( &( Input -> Data.y[j] ) , 0.0 , Nmeas , Raw ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  
  return SUCCESS ;
}

int
read_GLU_Qmoment( struct input_params *Input )
{
  // sanity check filenames
  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    if( Input -> Traj[i].FileX != NULL ) {
      fprintf( stderr , "[IO] read_GLU_qmoment only reads Y data \n" ) ;
      return FAILURE ;
    }
    if( Input -> Traj[i].FileY == NULL ) {
      fprintf( stderr , "[IO] traj_%zu entry has NULL Y value \n" , i ) ;
      return FAILURE ;
    }
  }

  // read all of the files and allocate data
  if( pre_allocate( Input ) == FAILURE ) {
    return FAILURE ;
  }

  FILE *file = NULL ;
  int Flag = SUCCESS ;

  // reread files and poke in the data
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // temporary string
    char str[ strlen( Input -> Traj[i].FileY ) + 6 ] ;

    double *Qmoments = malloc( Input -> Data.Ndata[i] * sizeof( double ) ) ;

    // linearised measurement index
    register size_t meas = 0 ;
    for( k = Input -> Traj[i].Begin ;
	 k < Input -> Traj[i].End ;
	 k += Input -> Traj[i].Increment ) {
      
      // open up a file
      sprintf( str , Input -> Traj[i].FileY , k ) ;

      if( ( file = fopen( str , "rb" ) ) == NULL ) {
	fprintf( stderr , "[IO] cannot open %s\n" , str ) ;
	Flag = FAILURE ;
	break ;
      }

      // read the first line
      uint32_t Nmeas[1] , Nmoments[1] ;

      // keep reading until we hit GSRC
      FREAD32( Nmeas , sizeof( uint32_t ) , 1 , file ) ;

      size_t Nmeas_counter = 0 ;
      for( Nmeas_counter = 0 ; Nmeas_counter < Nmeas[0] ; Nmeas_counter++ ) {

	FREAD32( Nmoments , sizeof( uint32_t ) , 1 , file ) ;
	
	FREAD64( Qmoments , sizeof( double ) , Nmoments[0] , file ) ;

	if( Nmeas_counter == Input -> Traj[i].Gs ) break ;
      }

      for( j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
	Input -> Data.x[j + shift].resampled[meas] = (double)j ;
	Input -> Data.y[j + shift].resampled[meas] = Qmoments[j] ; //- Qmoments[0] ;

	#ifdef verbose
	printf( "WHAT %f %f \n" ,
		Input -> Data.x[j + shift].resampled[meas] ,
		Input -> Data.y[j + shift].resampled[meas] ) ;
	#endif
      }

      fclose( file ) ;

      meas ++ ;
    }
    
    shift += Input -> Data.Ndata[i] ;
    
    free( Qmoments ) ;
  }

  return Flag ;
}
