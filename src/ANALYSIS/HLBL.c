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

// power of local current renormalisation
#define PL (4)

#define IN_FERMI
#define PUT_ZERO

// beta value of the ensemble of interest
#define B346
#define B450
//#define CONNECTED

#ifdef B340
  #define BETALABEL "Beta 3.4"
  #if (defined H101)
    #define SYMMETRIC
    #define CONFIGLABEL "H101"
  #elif (defined U103)
    #define SYMMETRIC
    #define CONFIGLABEL "U103"
  #elif (defined H102)
    #define CONFIGLABEL "H102"
  #elif (defined H105)
    #define CONFIGLABEL "H105"
  #elif (defined C101)
    #define CONFIGLABEL "C101"
  #else
    #error
  #endif
#elif (defined B346)
  #define SYMMETRIC
  #define BETALABEL "Beta 3.46"
  #define CONFIGLABEL "B450"
#elif (defined B355)
  #define BETALABEL "Beta 3.55"
  #if (defined H200)
    #define SYMMETRIC
    #define CONFIGLABEL "H200"
  #elif (defined N202)
    #define SYMMETRIC
    #define CONFIGLABEL "N202"
  #elif (defined N203)
    #define CONFIGLABEL "N203"
  #elif (defined N200)
    #define CONFIGLABEL "N200"
  #elif (defined D200)
    #define CONFIGLABEL "D200"
  #else
    #error
  #endif
#elif (defined B370)
  #define SYMMETRIC
  #define BETALABEL "Beta 3.7"
  #define CONFIGLABEL "N300"
#else
  #error
#endif

static void
integrate( struct input_params *Input ,
	   FILE *txtfile )
{
  const bool trap = true ;
  size_t i , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    if( trap == true ) {
      fprintf( txtfile , "\n\nIntegrating with the trapezoid rule\n" ) ;
    } else {
      fprintf( txtfile , "\n\nIntegrating with Simpson's rule\n" ) ;
    }
    fprintf( txtfile , "         |y|_max [fm]\t I(|y|_max)\t error\n" ) ;
    struct resampled Int[ Input -> Data.Ndata[i] ] ;

    // do a little trap for the implicit zero
    Int[0] = init_dist( &Input -> Data.y[shift] ,
			Input -> Data.x[shift].NSAMPLES ,
			Input -> Data.x[shift].restype ) ;
    #ifdef PUT_ZERO
    mult( &Int[0] , Input -> Data.x[shift] ) ;
    mult_constant( &Int[0] , 0.5 ) ;
    fprintf( txtfile , "Integral %e %e %e \n" ,
	       Input -> Data.x[shift].avg , Int[0].avg , Int[0].err ) ;
    #endif
    
    size_t k ;
    for( k = 1 ; k < Input -> Data.Ndata[i] ; k++ ) {
      // do the numberical integration
      Int[k] = Nint( Input -> Data.x + shift ,
		     Input -> Data.y + shift ,
		     k+1 , trap ) ;
      #ifdef IN_FERMI
      fprintf( txtfile , "Integral %e %e %e \n" ,
	       Input -> Data.x[shift+k].avg ,
	       Int[k].avg , Int[k].err ) ;
      #else
      fprintf( txtfile , "Integral %e %e %e \n" ,
	       Input -> Data.x[shift+k].avg*a ,
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
#ifdef SYMMETRIC
  #ifdef CONNECTED
  const double NUM[2] = { 18 , 18 } ;
  #else
  const double NUM[2] = { -36 , -36 } ;
  #endif
#else
  #ifdef CONNECTED
  const double NUM[2] = { 17 , 17 } ;
  #else
  const double NUM[2] = { 25 , 25 } ;
  #endif
#endif
  const double DEN = 81. ;

  FILE *txtfile = fopen( "HLBL_data.txt" , "w" ) ;
    
  //const double Q = 17/81. ; // -> 25/81. for the disconnected!
  // if we do the symmetric point this factor is 18/81. & 36/81.
  // values from the g-2 paper and the ZV paper
  // https://arxiv.org/pdf/1904.03120.pdf
  // https://arxiv.org/pdf/1811.08209.pdf
#ifdef B340
  const double amu = 0.04624130625372472 ;
  const double a = 0.08636 ;
  #if (defined H101)
  const double Z = 0.71562 ;
  #elif (defined U103)
  const double Z = 0.71562 ;
  #elif (defined H102)
  const double Z = 0.71226 ;
  #elif (defined H105)
  const double Z = 0.70908 ;
  #elif (defined C101)
  const double Z = 0.70717 ;
  #endif
#elif (defined B346)
  const double amu = 0.040876115324332385 ;
  const double a = 0.07634 ;
  const double Z = 0.71998 ;
#elif (defined B355)
  const double amu = 0.034407901110055004 ;
  const double a = 0.06426 ;
  #if (defined H200) 
  const double Z = 0.74028 ;
  #elif (defined N202)
  const double Z = 0.74028 ;
  #elif (defined N203)
  const double Z = 0.73792 ;
  #elif (defined N200)
  const double Z = 0.73614 ;
  #elif (defined D200)
  const double Z = 0.73429 ;
  #else
  fprintf( stdout , "Ensemble not recognised for Z_V\n") ;
  return FAILURE ;
  #endif
#elif (defined B370)
  const double amu = 0.02667067466996327 ;
  const double a = 0.04981 ;
  const double Z = 0.75413 ;
#else
  fprintf( stderr , "I do not understand your selected beta\n" ) ;
  return FAILURE ;
#endif
  const double ZV = pow(Z,PL) ;
  const double alpha_QED = 1.0/137.035999 ;
  const double prefac = ZV*amu*pow(4*M_PI*alpha_QED,3)/3.*2*M_PI*M_PI ;

  fprintf( txtfile , "Using %g/%g charge factor\n" , NUM[0] , DEN ) ;
  fprintf( txtfile , "Config is %s, %s\n" , CONFIGLABEL , BETALABEL ) ;
  fprintf( txtfile , "a = %g [fm], amu %g\n" , a , amu ) ;
  fprintf( txtfile , "Renormalisation (%g)^%d\n" , Z , PL ) ;
  fprintf( txtfile , "          |y| [fm]\tf(|y|)[fm^-1]\terror\n" ) ;
  
  size_t i , j , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    const double Q = NUM[i%2]/DEN ;

    #ifdef PUT_ZERO
    fprintf( txtfile , "Integrand %e %e %e\n" , 0. , 0. , 0. ) ;
    #endif
    
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

      fprintf( txtfile , "Integrand %e %e %e\n" ,
	       Input->Data.x[j].avg ,
	       Input -> Data.y[j].avg ,
	       Input -> Data.y[j].err ) ;
    }    
    shift += Input -> Data.Ndata[i] ;
  }

  /*
  // add the two
  if( Input -> Data.Nsim == 2 ) {
    for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
      add( &Input -> Data.y[j] ,
	   Input -> Data.y[j+Input->Data.Ndata[0]] ) ;
    }
  }
  */
  
  integrate( Input , txtfile ) ;

  fclose( txtfile ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot_and_Nint( *Input , &Chi ) ;

  free_fitparams( Fit , Input -> Fit.Nlogic ) ;  
  
  return SUCCESS ;
}
