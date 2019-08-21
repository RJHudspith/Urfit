/**
   @file HLBL.c
   @brief integrates the HLBL from Jamie's data
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "Nint.h"
#include "resampled_ops.h"
#include "write_flat.h"

// beta value of the ensemble of interest
#define B355

// power of local current renormalisation
#define PL (4)

#define IN_FERMI

static void
integrate( struct input_params *Input )
{
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    fprintf( stdout , "\nPerforming integrals\n" ) ;
    struct resampled Int[ Input -> Data.Ndata[i] ] ;
    Int[0] = init_dist( NULL ,
			Input -> Data.x[shift].NSAMPLES ,
			Input -> Data.x[shift].restype ) ;
    size_t k ;
    for( k = 1 ; k < Input -> Data.Ndata[i] ; k++ ) {
      // do the numberical integration
      Int[k] = Nint( Input -> Data.x + shift ,
		     Input -> Data.y + shift ,
		     k+1 , true ) ;
      #ifdef IN_FERMI
      fprintf( stdout , "Integral %zu %e %e %e \n" ,
	       i ,
	       Input -> Data.x[shift+k].avg ,
	       Int[k].avg , Int[k].err ) ;
      #else
      fprintf( stdout , "Integral %zu %e %e %e \n" ,
	       i , Input -> Data.x[shift+k].avg*a ,
	       Int[k].avg , Int[k].err ) ;
      #endif      
    }

    char str[256] ;
    sprintf( str , "Integral_%zu.flat" , i ) ;
    write_flat_dist( Int , Input->Data.x+shift ,
		     Input -> Data.Ndata[i] , str ) ;

    // write out a flat distribution
    for( k = 0 ; k < Input -> Data.Ndata[i] ; k++ ) {
      free( Int[k].resampled ) ;
    }

    shift += Input -> Data.Ndata[i] ;
  }
  return ;
}

int
HLBL_analysis( struct input_params *Input )
{
  const double NUM[2] = { 17 , 17 } ;
  const double DEN = 81. ;

  fprintf( stdout , "[HLBL] using %g/%g charge factor\n" , NUM[0] , DEN ) ;
  
  //const double Q = 17/81. ; // -> 25/81. for the disconnected!
  // if we do the symmetric point this factor is 18/81.
  // values from the g-2 paper and the ZV paper
  // https://arxiv.org/pdf/1904.03120.pdf
  // https://arxiv.org/pdf/1811.08209.pdf
#ifdef B340
  const double amu = 0.04624130625372472 ;
  const double a = 0.08636 ;
  const double ZV = pow(0.70912,PL) ;
#elif (defined B346)
  const double amu = 0.040876115324332385 ;
  const double a = 0.07634 ;
  const double ZV = pow(0.71998,PL) ;
#elif (defined B355)
  const double amu = 0.034407901110055004 ;
  const double a = 0.06426 ;
  const double ZV = pow(0.73453,PL) ;
#elif (defined B370)
  const double amu = 0.02667067466996327 ;
  const double a = 0.04981 ;
  const double ZV = pow(0.75413,PL) ;
#else
  fprintf( stderr , "I do not understand your selected beta\n" ) ;
  return FAILURE ;
#endif
  const double alpha_QED = 1.0/137.035999 ;
  const double prefac = ZV*amu*pow(4*M_PI*alpha_QED,3)/3.*2*M_PI*M_PI ;
  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    const double Q = NUM[i%2]/DEN ;
    
    // multiply by |x|^3*prefac
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      
      // x is |y| when comparing to note
      root( &Input -> Data.x[j] ) ;

      #ifdef IN_FERMI
      mult_constant( &Input -> Data.y[j] ,
		     Q*prefac*pow( Input -> Data.x[j].avg , 3 )/a ) ;
      mult_constant( &Input -> Data.x[j] , a ) ;
      #else
      mult_constant( &Input -> Data.y[j] ,
		     Q*prefac*pow( Input -> Data.x[j].avg , 3 ) ) ;
      #endif

      fprintf( stdout , "%e %e %e \n" , Input->Data.x[j].avg ,
	       Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
    }    
    shift += Input -> Data.Ndata[i] ;
  }

  // add the two
  /*
  if( Input -> Data.Nsim == 2 ) {
    for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
      add( &Input -> Data.y[j] ,
	   Input -> Data.y[j+Input->Data.Ndata[0]] ) ;
    }
  }
  */
  integrate( Input ) ;
  
  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot_and_Nint( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;  
  
  return SUCCESS ;
}
