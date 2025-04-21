/**
   @file general_ops.c
   @brief general operations on data
 */
#include "gens.h"

#include "resampled_ops.h"
#include "write_flat.h"

static int
do_op( const struct resampled A ,
       const struct resampled B ,
       void (*f)( struct resampled *a ,
		  const struct resampled b ) ,
       const char *s ,
       const double x ,
       const size_t j )
{
  struct resampled res = init_dist( &A , A.NSAMPLES , A.restype ) ;

  fprintf( stdout , "\n--------------------------------\n"
	   "[GEN] A %e %e :: B %e %e \n" ,
	   A.avg , A.err , B.avg , B.err ) ;

  f( &res , B ) ;

  fprintf( stdout , "[GEN] %s %e +/- %e \n\n" , s , res.avg , res.err ) ;

  // write out a flat file with the moniker
  char *str = malloc( 256 * sizeof( char ) ) ;
  sprintf( str , "%s_%zu.flat" , s , j ) ;

  printf( "Flat %s \n" , str ) ;
  //write_flat_single( res , str ) ;

  struct resampled mpi2 = init_dist( NULL , A.NSAMPLES , A.restype ) ;
  equate_constant( &mpi2 , x , A.NSAMPLES , A.restype ) ;
  write_flat_dist( &res , &mpi2 , 1 , str ) ;

  free( str ) ;
  free( mpi2.resampled ) ;
  free( res.resampled ) ;

  return SUCCESS ;
}

int
gen_ops( struct input_params *Input )
{
  if( Input -> Data.Nsim == 2 ) {
    if( Input -> Data.Ndata[0] != Input -> Data.Ndata[1] ) {
      fprintf( stderr , "[GEN] two distributions have different lengths\n" ) ;
      return FAILURE ;
    }

    printf( "WTF? %zu %zu %zu \n" , Input -> Data.Nsim , Input -> Data.Ndata[0] ,
	    Input -> Data.Ntot ) ;
    
    size_t j ;
    for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {

      const double mpi2 = Input -> Data.x[j].avg ; //Input -> Traj[0].Fit_Low ;
      
      // add
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     add , "Add" , mpi2 , j ) ;
      
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     subtract , "Sub" , mpi2 , j ) ;
      
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     mult , "Mult" , mpi2 , j ) ;

      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     divide , "Div" , mpi2 , j ) ;

      /*
      mult_constant( &Input -> Data.y[j+Input -> Data.Ndata[0]] , 6/2. ) ;
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     subtract , "SubNC" , mpi2 , j ) ;
      */

      /*
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     phi4_comb , "Phi4_comb" , mpi2 , j ) ;
      */
      
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     spin_average , "SpinAve" , mpi2 , j ) ;

      // computes x^1/2 / y
      raise( &Input -> Data.y[j] , 0.5 ) ;
      //mult_constant( &Input -> Data.y[j] , 12 ) ;
      //raise( &Input -> Data.y[j+Input -> Data.Ndata[0]] , 2 ) ;
      
      do_op( Input -> Data.y[j] ,
	     Input -> Data.y[j+Input -> Data.Ndata[0]] ,
	     divide , "x^2Dy^2" , mpi2 , j ) ;
      //#endif
    }
  }
  return SUCCESS ;
}
