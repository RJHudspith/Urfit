/**
   @file adler.c
   @brief adler function analysis
 */
#include "gens.h"

#include "adler_alpha_D0.h"
#include "adler_alpha_D0_multi.h"
#include "cruel_runnings.h"
#include "fit_and_plot.h"
#include "init.h"
#include "resampled_ops.h"
#include "rng.h"
#include "stats.h"
#include "write_flat.h"

static void
write( struct input_params *Input )
{
  write_flat_file( *Input , "Adler" ) ;
}

static void
analysis( struct input_params *Input )
{
  // perform a fit perchance? or do something better
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
}

static void
write_renorm( struct input_params *Input )
{
  const int Nmom = Input -> Data.Ndata[0] ;
  printf( "Nmom -> %zu \n" , Nmom ) ;
  raise( &Input -> Data.y[Nmom+0] , 2 ) ;
  printf( "WTF1 %e %e \n" , Input->Data.y[Nmom].avg , Input -> Data.y[Nmom].err ) ;
  raise( &Input -> Data.y[Nmom+1] , 2 ) ;
  printf( "WTF2 %e %e \n" , Input->Data.y[Nmom+1].avg , Input -> Data.y[Nmom+1].err ) ;
  for( int p = 0 ; p < Nmom ; p++ ) {
    mult( &Input -> Data.x[p] , Input -> Data.y[Nmom+1] ) ;
    mult( &Input -> Data.y[p] , Input -> Data.y[Nmom+0] ) ;
  }

  write_flat_dist( Input -> Data.y , Input -> Data.x , Input -> Data.Ndata[0] , "Adler_ren.flat") ;
}

// this is for renormalised
// will write a flat file of t_0/a^2 \Pi for each t0 Q^2
static void
meta_analysis( struct input_params *Input )
{
  const size_t NQ = Input -> Data.Ndata[0] ;
  const size_t NT = Input -> Data.Nsim/2 ;

  printf( "NQ %zu | NT %zu\n" , NQ , NT ) ;

  // half because we interleave t0
  struct resampled *Y = malloc( NT*sizeof( struct resampled ) ) ;
  struct resampled *Inv_t0 = malloc( NT*sizeof( struct resampled ) ) ;
  for( size_t j = 0 ; j < NT ; j++ ) {
    Y[j] = init_dist( NULL , Input -> Data.y[0].NSAMPLES , Input -> Data.y[0].restype) ;
    Inv_t0[j] = init_dist( NULL , Input -> Data.y[0].NSAMPLES , Input -> Data.y[0].restype) ;
  }
  
  for( size_t i = 0 ; i < NQ ; i++ ) {
    // t0 Q^2-indexed file name
    char str[64] ;
    sprintf( str , "Adep_%f.flat" , Input -> Data.x[i].avg ) ;

    for( size_t j = 0 ; j < NT ; j++ ) {
      equate( &Y[j] , Input -> Data.y[ j + j*NQ + i ] ) ;
      add_constant( &Y[j] , -1 ) ;
      equate( &Inv_t0[j] , Input -> Data.y[ (j+1)*NQ+j ] ) ;
      raise( &Inv_t0[j] , -2 ) ;

      printf( "WTF (%zu %zu) %f %f \n" , j + j*NQ + i , (j+1)*NQ+j , Inv_t0[j].avg , Y[j].avg ) ;
    }
    //exit(1) ;
    write_flat_dist( Y , Inv_t0 , NT , str ) ;    
  }
  return ;
}

static void
test_sub( struct input_params *Input )
{
  const int Nmom = Input -> Data.Ndata[0] ;
  
  //const double a = 1.635 , ZV = 9.736695e-01 ;
  //const double a = 1.97  , ZV = 9.752028e-01 ;
  //const double a = 2.46  , ZV = 9.790892e-01 ;
  //const double a = 3.140 , ZV = 9.816336e-01 ;//, ZV2 = 0.995 ;
  //const double a = 3.990 , ZV = 0.988 ;//, ZV2 = 0.995 ;

  
  //mult_constant( &Input -> Data.y[Nmom] , ZV2*ZV2 ) ;
  
  for( int i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    //mult_constant( &Input -> Data.y[p] , ZV*ZV ) ;
    //subtract( &Input -> Data.y[p] , Input -> Data.y[p+Nmom] ) ;
    add_constant( &Input -> Data.y[i] , -1.0 ) ;
    //mult_constant( &Input -> Data.y[p] , -1.0 ) ;
    mult_constant( &Input -> Data.x[i] , 1/0.5372839752835529 ) ;
    raise( &Input -> Data.x[i] , 0.5 ) ;
  }
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
}

int
adler_analysis( struct input_params *Input )
{
  //write( Input ) ;
  //analysis( Input ) ;
  //write_renorm( Input ) ;
  test_sub( Input ) ;

  //meta_analysis( Input ) ;
  return SUCCESS ;
}
