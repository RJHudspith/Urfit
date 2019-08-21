/**
   @file read_flat.c
   @brief read (plain text!) data from a file

   Data order is :
   Sample Type - size_t
   Ndata - size_t
   Nsamples - size_t

   .. data

   AVG %lf %lf
 */
#include "gens.h"

#include "stats.h"

// reads flat single data 
struct resampled*
read_flat_single( struct input_params *Input )
{
  return NULL ;
}

// 
static int
read_initial( FILE *file ,
	      size_t *Restype ,
	      size_t *Ndata )
{
  fscanf( file , "%zu\n" , Restype ) ;
  if( *Restype > 2 ) {
    fprintf( stderr , "[IO] unknown sample index %zu \n" , *Restype ) ;
    return FAILURE ;
  }
  fscanf( file , "%zu\n" , Ndata ) ;
  if( *Ndata == 0 ) {
    fprintf( stderr , "[IO] Ndata is zero \n" ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// expects x y data layout with a single space between x and y
static int
read_flat_double( struct input_params *Input ,
		  const size_t Traj_idx ,
		  const size_t xy ,
		  const size_t shift )
{
  FILE *file = NULL ;
  if( xy == 0 ) {
    file = fopen( Input -> Traj[ Traj_idx ].FileX , "r" ) ;
  } else {
    file = fopen( Input -> Traj[ Traj_idx ].FileY , "r" ) ;
  }
  
  if( file == NULL ) {
    fprintf( stderr , "[IO] cannot open file %s %s \n" ,
	     Input -> Traj[ Traj_idx ].FileX , Input -> Traj[ Traj_idx ].FileY ) ;
    return FAILURE ;
  }

  // read the sample type
  size_t Restype , Ndata , Nsamples ;
  read_initial( file , &Restype , &Ndata ) ;

  // do some loops, allocate x and y
  Input -> Data.Ndata[ Traj_idx ] = Ndata ;
  
  size_t i , j ;
  for( i = 0 ; i < Ndata ; i++ ) {
    // should be another NSAMPLES here
    fscanf( file , "%zu" , &Nsamples ) ;

    printf( "Here NSAMPLES %zu\n" , Nsamples ) ;
    
    for( j = 0 ; j < Nsamples ; j++ ) {
      fscanf( file , "%lf %lf\n" ,
	      &Input -> Data.x[ shift + i ].resampled[j] ,
	      &Input -> Data.y[ shift + i ].resampled[j] ) ;
      #ifdef VERBOSE
      printf( "THIS :: %f %f \n" ,
	      Input -> Data.x[ shift + i ].resampled[j] ,
	      Input -> Data.y[ shift + i ].resampled[j] ) ;
      #endif
    }
    // do the average
    if( Input -> Data.x[ shift + i ].restype != Raw) {
      fscanf( file , "AVG %lf %lf\n" ,
	      &Input -> Data.x[ shift + i ].avg ,
	      &Input -> Data.y[ shift + i ].avg ) ;
    }
    compute_err( &Input -> Data.x[ shift + i ] ) ;
    compute_err( &Input -> Data.y[ shift + i ] ) ;

    printf( "AVERAGE read %e %e \n" , Input -> Data.x[ shift + i ].avg ,
	    Input -> Data.y[ shift + i ].avg ) ;
  }

  fclose(file) ;

  return SUCCESS ;
}

static int
init_data( struct input_params *Input )
{
  Input -> Data.Ndata = malloc( Input -> Data.Nsim * sizeof( size_t ) ) ;
  Input -> Data.Ntot = 0 ;

  // set these up
  size_t
    Restype[ Input -> Data.Nsim ] ,
    Ndata[ Input -> Data.Nsim ] ;
  
  size_t i ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    FILE *file = NULL ;
    file = fopen( Input -> Traj[ i ].FileY , "r" ) ;
  
    if( file == NULL ) {
      fprintf( stderr , "[IO] cannot open file %s %s \n" ,
	     Input -> Traj[ i ].FileX , Input -> Traj[ i ].FileY ) ;
      return FAILURE ;
    }

    // read the sample type
    if( read_initial( file , &Restype[i] , &Ndata[i] ) == FAILURE ) {
      fprintf( stderr , "[IO] failed to read initial sample type \n" ) ;
      return FAILURE ;
    }

    Input -> Data.Ndata[i] = Ndata[i] ;
    Input -> Data.Ntot += Ndata[i] ;
    
    fclose( file ) ;
  }

  // allocate the x and y data
  Input -> Data.x = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;
  Input -> Data.y = malloc( Input -> Data.Ntot * sizeof( struct resampled ) ) ;

  size_t shift = 0 , j ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    FILE *file = NULL ;
    file = fopen( Input -> Traj[ i ].FileY , "r" ) ;

    // read the sample type
    if( read_initial( file , &Restype[i] , &Ndata[i] ) == FAILURE ) {
      return FAILURE ;
    }
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      size_t Nsamples ;
      fscanf( file , "%zu" , &Nsamples ) ;
      printf( "NSAMPLE :: %zu \n" , Nsamples ) ;
      
      Input -> Data.x[ j ].resampled = malloc( Nsamples * sizeof( double ) ) ;
      Input -> Data.x[ j ].restype = Restype[i] ;
      Input -> Data.x[ j ].NSAMPLES = Nsamples ;
    
      Input -> Data.y[ j ].resampled = malloc( Nsamples * sizeof( double ) ) ;
      Input -> Data.y[ j ].restype = Restype[i] ;
      Input -> Data.y[ j ].NSAMPLES = Nsamples ;

      size_t k ;
      for( k = 0 ; k < Nsamples ; k++ ) {
	double a , b ;
	fscanf( file , "%lf %lf\n" , &a , &b ) ; 
      }
    }
    shift = j ;

    fclose( file ) ;
  }
  
  return SUCCESS ;
}

int
read_flat( struct input_params *Input )
{
  if( init_data( Input ) == FAILURE ) {
    fprintf( stderr , "[IO] data initialisation failure \n" ) ;
    return FAILURE ;
  }

  size_t i , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    // x and y data in different files

    // xy data in one file
    if( Input -> Traj[i].FileX == NULL && Input -> Traj[i].FileY != NULL ) {
      if( read_flat_double( Input , i , 1 , shift ) == FAILURE ) {
	return FAILURE ;
      }
    } else if( Input -> Traj[i].FileY == NULL && Input -> Traj[i].FileX != NULL ) {
      if( read_flat_double( Input , i , 0 , shift ) == FAILURE ) {
	return FAILURE ;
      }
    }
    shift += Input -> Data.Ndata[i] ;
  }
  Input -> Data.Ntot = shift ;

  #ifdef VERBOSE
  printf( "WHAT :: %zu %zu %zu \n" , Input -> Data.Nsim , Input -> Data.Ndata[0] , Input -> Data.Ntot ) ;

  // test
  shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    size_t j , k ;
    for( j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
      for( k = 0 ; k < Input -> Data.x[ shift + j ].NSAMPLES ; k++ ) {
	printf( "TEST :: %f %f \n" ,
		Input -> Data.x[ shift + j ].resampled[ k ] ,
		Input -> Data.y[ shift + j ].resampled[ k ] ) ;
      }
    }
    shift += Input -> Data.Ndata[i] ;
  }
  #endif
  
  return SUCCESS ;  
}
