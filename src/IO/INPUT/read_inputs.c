/**
   @file read_inputs.c
   @brief read an input file
 */
#include "gens.h"

#include <string.h>

#include "read_fit.h"
#include "read_traj.h"

// set these to something reasonable
#define STR1_LENGTH (64)
#define STR2_LENGTH (256)

// get the file length
static size_t
get_file_length( FILE *Infile )
{
  size_t Ntags = 0 ;
  char str1[ STR1_LENGTH ] , str2[ STR2_LENGTH ] ;
  while( Ntags++ , fscanf( Infile , "%s = %s" , str1 , str2 )  != EOF ) { }
  rewind( Infile ) ;
  return Ntags-1 ;
}

// allocate the input file struct
static struct flat_file *
pack_inputs( FILE *Infile ,
	     const size_t Ntags )
{
  struct flat_file *Flat = malloc( Ntags * sizeof( struct flat_file ) ) ;
  char str1[ STR1_LENGTH ] , str2[ STR2_LENGTH ] ;
  size_t i ;
  for( i = 0 ; i < Ntags ; i++ ) {
    if( fscanf( Infile , "%s = %s" , str1 , str2 ) != 2 ) {
      continue ;
    }
    Flat[i].Token_Length = strlen( str1 ) ;
    Flat[i].Token = malloc( Flat[i].Token_Length * sizeof( char ) ) ;
    sprintf( Flat[i].Token , "%s" , str1 ) ;
    
    Flat[i].Value_Length = strlen( str2 ) ;
    Flat[i].Value = malloc( Flat[i].Value_Length * sizeof( char ) ) ;
    sprintf( Flat[i].Value , "%s" , str2 ) ;
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
  for( i = 0 ; i < Input -> Ntraj ; i++ ) {
    free( Input -> Traj[i].Dimensions ) ;
    free( Input -> Traj[i].Filename ) ;
  }
  free( Input -> Traj ) ;
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
    Flag == FAILURE ;
  }
  
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
