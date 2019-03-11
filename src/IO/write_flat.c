/**
   @file write_flat.c
   @brief write out a flat file
 */
#include "gens.h"

int
write_flat_file( const struct input_params Input ,
		 const char *name )
{
  // write out a flat file with the resamples?
  printf( "Writing a flat file to %s \n" , name ) ;
  
  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input.Data.Nsim ; i++ ) {
    char str[256] ;
    sprintf( str , "%s.%zu.flat" , name , i ) ;
    
    FILE *outfile = fopen( str , "w" ) ;
    fprintf( outfile , "%d\n" , Input.Data.y[shift].restype ) ;
    fprintf( outfile , "%zu\n" , Input.Data.Ndata[i] ) ;

    for( j = shift ; j < shift + Input.Data.Ndata[i] ; j++ ) {
      fprintf( outfile , "%zu\n" , Input.Data.y[j].NSAMPLES ) ;
      for( k = 0 ; k < Input.Data.y[j].NSAMPLES ; k++ ) {
	fprintf( outfile , "%1.15e %1.15e\n" ,
		 Input.Data.x[j].resampled[k] ,
		 Input.Data.y[j].resampled[k] ) ;
      }
      fprintf( outfile , "AVG %1.15e %1.15e\n" ,
	       Input.Data.x[j].avg ,
	       Input.Data.y[j].avg ) ;
    }
    shift += Input.Data.Ndata[i] ;
    fclose( outfile ) ;
  }

  return FAILURE ;
}

int
write_flat_dist( const struct resampled *y ,
		 const struct resampled *x ,
		 const size_t Ndata ,
		 const char *name )
{
  printf( "%s \n" , name ) ;
  // write out a flat file with the resamples?
  printf( "Writing a flat file to %s \n" , name ) ;
  
  size_t j , k ;

  FILE *outfile = fopen( name , "w" ) ;
  
  fprintf( outfile , "%zu\n" , (size_t)y[0].restype ) ;
  fprintf( outfile , "%zu\n" , Ndata ) ;

  for( j = 0 ; j < Ndata ; j++ ) {
    fprintf( outfile , "%zu\n" , y[j].NSAMPLES ) ;
    for( k = 0 ; k < y[j].NSAMPLES ; k++ ) {
      fprintf( outfile , "%1.15e %1.15e\n" ,
	       x[j].resampled[k] ,
	       y[j].resampled[k] ) ;
    }
    fprintf( outfile , "AVG %1.15e %1.15e\n" ,
	     x[j].avg , y[j].avg ) ;
  }
  fclose( outfile ) ;

  return FAILURE ;
}

int
write_flat_single( const struct resampled *y ,
		   const char *name )
{
  // write out a flat file with the resamples?
  fprintf( stdout , "Writing a flat file to %s \n" , name ) ;
  
  size_t k ;

  FILE *outfile = fopen( name , "w" ) ;
  
  fprintf( outfile , "%zu\n" , (size_t)y[0].restype ) ;
  fprintf( outfile , "1\n" ) ;

  fprintf( outfile , "%zu\n" , y[0].NSAMPLES ) ;
  for( k = 0 ; k < y[0].NSAMPLES ; k++ ) {
    fprintf( outfile , "%1.15e\n" , y[0].resampled[k] ) ;
  }
  fprintf( outfile , "AVG %1.15e\n" , y[0].avg ) ;

  fclose( outfile ) ;

  return FAILURE ;
}
