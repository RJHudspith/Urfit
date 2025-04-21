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

static void
test_sub( struct input_params *Input )
{
  const int Nmom = Input -> Data.Ndata[0] ;
  
  //const double a = 1.635 , ZV = 9.736695e-01 ;
  //const double a = 1.97  , ZV = 9.752028e-01 ;
  const double a = 2.46  , ZV = 9.790892e-01 ;
  //const double a = 3.140 , ZV = 9.816336e-01 ;//, ZV2 = 0.995 ;
  //const double a = 3.990 , ZV = 0.988 ;//, ZV2 = 0.995 ;

  
  //mult_constant( &Input -> Data.y[Nmom] , ZV2*ZV2 ) ;
  

  for( int p = 0 ; p < Nmom ; p++ ) {
    mult_constant( &Input -> Data.y[p] , ZV*ZV ) ;
    //subtract( &Input -> Data.y[p] , Input -> Data.y[p+Nmom] ) ;
    add_constant( &Input -> Data.y[p] , -1.0 ) ;
    //mult_constant( &Input -> Data.y[p] , -1.0 ) ;
    mult_constant( &Input -> Data.x[p] , a*a ) ;
  }
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
}

int
adler_analysis( struct input_params *Input )
{
  write( Input ) ;
  //analysis( Input ) ;
  //write_renorm( Input ) ;
  //test_sub( Input ) ;
  return SUCCESS ;
}
