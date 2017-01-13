/**
   @file read_inputs.c
   @brief read an input file
 */
#include "gens.h"

#include <string.h>

#include "init.h"      // free_fit && free_data
#include "read_fit.h"
#include "read_graph.h"
#include "read_traj.h"
#include "read_stats.h"

// set these to something reasonable
#define STR1_LENGTH (64)
#define STR2_LENGTH (256)

// get the file length
static size_t
get_file_length( FILE *Infile )
{
  size_t Ntags = 0 ;
  char str1[ STR1_LENGTH ] , str2[ STR2_LENGTH ] ;
  while( true ) {
    const int test = fscanf( Infile , "%s = %s" , str1 , str2 ) ;
    if( test == EOF ) break ; // if the file is finished we leave
    if( test != 2 ) continue ; // if there are comments and stuff
    Ntags ++ ;
  }
  rewind( Infile ) ;
  return Ntags ;
}

// allocate the input file struct
static struct flat_file *
pack_inputs( FILE *Infile ,
	     const size_t Ntags )
{
  struct flat_file *Flat = malloc( Ntags * sizeof( struct flat_file ) ) ;
  char str1[ STR1_LENGTH ] , str2[ STR2_LENGTH ] ;
  size_t i = 0 ;
  while( true ) {
    const int test = fscanf( Infile , "%s = %s" , str1 , str2 ) ;

    if( test == EOF ) break ; // if the file is finished we leave
    if( test != 2 ) continue ; // if there are comments and stuff
    
    Flat[i].Token_Length = strlen( str1 ) ;
    Flat[i].Token = malloc( Flat[i].Token_Length * sizeof( char ) ) ;
    sprintf( Flat[i].Token , "%s" , str1 ) ;
    
    Flat[i].Value_Length = strlen( str2 ) ;
    Flat[i].Value = malloc( Flat[i].Value_Length * sizeof( char ) ) ;
    sprintf( Flat[i].Value , "%s" , str2 ) ;
    
    i++ ;
  }
  return Flat ;
}

// frees up the struct
static void
unpack_inputs( struct flat_file *Flat ,
	       const size_t Ntags )
{
  size_t i ;
  for( i = 0 ; i < Ntags ; i++ ) {
    free( Flat[i].Token ) ;
    free( Flat[i].Value ) ;
  }
  free( Flat ) ;
  return ;
}

// test if two strings are equivalent
int
are_equal( const char *str_1 , const char *str_2 ) 
{
  return ( strcmp( str_1 , str_2 ) != 0 ) ? 0 : 1 ;
}

void
free_inputs( struct input_params *Input )
{ 
  size_t i ;
  if( Input -> Traj != NULL ) {
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      free( Input -> Traj[i].Dimensions ) ;
      free( Input -> Traj[i].Filename ) ;
    }
    free( Input -> Traj ) ;
  }
  free( Input -> Graph.Name ) ;
  free( Input -> Graph.Xaxis ) ;
  free( Input -> Graph.Yaxis ) ;
  return ;
}

// read the input file
int
read_inputs( struct input_params *Input ,
	     const char *filename )
{
  // open the file
  FILE *infile = fopen( filename , "r" ) ;
  int Flag = SUCCESS ;

  // initialise any arrays of Input to NULL
  Input -> Data.x = NULL ;
  Input -> Data.y = NULL ;
  Input -> Data.Ndata = NULL ;
  Input -> Data.Cov.W = NULL ;
  Input -> Data.Ntot = 0 ;
  Input -> Data.Nsim = 0 ;
  
  Input -> Fit.Sims = NULL ;
  Input -> Traj = NULL ;
  Input -> Fit.Prior = NULL ;
  Input -> Fit.map = NULL ;

  Input -> Graph.Name = NULL ;
  Input -> Graph.Xaxis = NULL ;
  Input -> Graph.Yaxis = NULL ;
  
  if( infile == NULL ) {
    printf( "[INPUT] Cannot find input file %s \n" , filename ) ;
    return Flag ;
  }

  // get the flat file length
  const size_t Ntags = get_file_length( infile ) ;

  printf( "[INPUTS] Ntags = %zu \n" , Ntags ) ;

  // pack the input file
  struct flat_file *Flat = pack_inputs( infile , Ntags ) ;

  if( get_traj( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }

  if( get_fit( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }

  if( get_stats( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }

  if( get_graph( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }
  printf( "\n" ) ;
  
  // close and unpack
  unpack_inputs( Flat , Ntags ) ;
  
  fclose( infile ) ;

  return Flag ;
}

// search for a specific tag
size_t
tag_search( const struct flat_file *Flat ,
	    const char *tag ,
	    const size_t Start ,
	    const size_t Ntags ) 
{
  size_t i ;
  for( i = Start ; i < Ntags ; i++ ) {
    if( are_equal( Flat[i].Token , tag ) ) {
      return i ;
    }
  }
  return i ;
}
