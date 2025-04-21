/**
   @file reweight.c
   @brief reweight the data we have measured

   expects a single flat RW file with

   cnfg_idx rw_factor

   in it
 **/
#include "gens.h"

#include "read_inputs.h"
#include "resampled_ops.h"
#include "stats.h"

#define CONTINUOUS

static double
get_rw( double *rwfac ,
	FILE *file ,
	const size_t cnfg )
{
  *rwfac = 0.0 ;
  size_t cidx ;

  //fprintf( stdout , "searching for cnfg %zu\n" , cnfg ) ;

#ifndef CONTINUOUS
  while( fscanf( file , "%zu %le\n" , &cidx , rwfac ) != EOF ) {
    printf( "%zu %1.15e\n" , cidx , *rwfac ) ;
    if( cidx == cnfg ) return SUCCESS ;
  }
  fprintf( stderr , "[RW] matching rw factor for cnfg %zu not found\n" ,
	   cnfg ) ;
  return FAILURE ;
#else
  fscanf( file , "%zu %le\n" , &cidx , rwfac ) ;
  printf( "%zu %1.15e\n" , cidx , *rwfac ) ;
  return SUCCESS ;
#endif
}

int
reweight_data( struct input_params *Input )
{  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    
    // NULL means no-reweighting
    if( Input -> Traj[i].RW == NULL ) {
      fprintf( stdout , "[RW] traj_%zu not being reweighted\n" , i ) ;
      continue ;
    }
    FILE *file = fopen( Input -> Traj[i].RW , "r" ) ;
    if( file == NULL ) {
      fprintf( stderr , "[IO] RW file %s is empty\n" , Input -> Traj[i].RW ) ;
      return FAILURE ;
    }
    
    size_t cnfg , k = 0 ;
    for( cnfg = Input -> Traj[i].Begin ;
	 cnfg < Input -> Traj[i].End ;
	 cnfg += Input -> Traj[i].Increment ) {
      
      double rwfac ;
      if( get_rw( &rwfac , file , cnfg ) == FAILURE ) {
	fclose( file ) ;
	return FAILURE ;
      }

      // and multiply all data with it
      for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
	Input -> Data.y[j].resampled[k] *= rwfac ;
      }
      k++ ;
    }
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      compute_err( &Input -> Data.y[j] ) ;
    }    
    shift += Input -> Data.Ndata[i] ;
    
    fclose( file ) ;
  }
  
  return SUCCESS ;
}
