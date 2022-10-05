/**
   @file HLBL.c
   @brief integrates the HLBL from Jamie's data
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "Nint.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"
#include "write_flat.h"

// power of local current renormalisation
#define PL (4)

//#define NINT_HLBL
#define PUT_ZERO
#define LIGHT_ONLY

#ifndef NINT_HLBL
  #define IN_FERMI
  #define MULT_X3
#endif

// beta value of the ensemble of interest
#define B385
#define J500
#define CONNECTED

//#define STRANGE
#define CHARM
#define CHARMZV

//#define LS
//#define LC
//#define CC

static const double lerp_pt = 1 ;

#ifdef B334
  #define BETALABEL "Beta 3.34"
  #if (defined A654)
    #define CONFIGLABEL "A654"
  #elif (defined A653)
    #define SYMMETRIC
    #define CONFIGLABEL "A653"
  #else
    #error
  #endif
#elif (defined B340)
  #define BETALABEL "Beta 3.4"
  #if (defined H101)
    #define SYMMETRIC
    #define CONFIGLABEL "H101"
  #elif (defined U103)
    #define SYMMETRIC
    #define CONFIGLABEL "U103"
  #elif (defined U102)
    #define CONFIGLABEL "U102"
  #elif (defined H102)
    #define CONFIGLABEL "H102"
  #elif (defined U102)
    #define CONFIGLABEL "U102"
  #elif (defined U101)
    #define CONFIGLABEL "U101"
  #elif (defined H105)
    #define CONFIGLABEL "H105"
  #elif (defined C101)
    #define CONFIGLABEL "C101"
  #else
    #error
  #endif
#elif (defined B346)
  #define BETALABEL "Beta 3.46"
  #if (defined B450)
    #define SYMMETRIC
    #define CONFIGLABEL "B450"
  #elif(defined D450)
    #define CONFIGLABEL "D450"
  #else
    #error
  #endif
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
#elif (defined B385)
  #define SYMMETRIC
  #define BETALABEL "Beta 3.85"
  #define CONFIGLABEL "J500"
#else
  #error
#endif

#ifdef LIGHT_ONLY
  #undef SYMMETRIC
#endif

static void
integrate( struct input_params *Input ,
	   FILE *txtfile ,
	   const double a )
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
      #ifdef PUT_ZERO
      add( &Int[k] , Int[0] ) ;
      #endif
      fprintf( txtfile , "Integral %e %e %e \n" ,
	       Input -> Data.x[shift+k].avg ,
	       Int[k].avg , Int[k].err ) ;  
    }

    // write out a flat distribution
    char str[256] ;
    sprintf( str , "Integral_%zu.flat" , i ) ;
    write_flat_dist( Int , Input->Data.x+shift ,
		     Input -> Data.Ndata[i] , str ) ;

    #ifdef PUT_ZERO
    Int[0] = Nint_pt( Input->Data.x+shift , Input->Data.y+shift ,
		      Input->Data.Ndata[i] , lerp_pt , true ) ;
    #else
    Int[0] = Nint_pt( Input->Data.x+shift , Input->Data.y+shift ,
		      Input->Data.Ndata[i] , lerp_pt , false ) ;
    #endif
    
    printf( "Lerped int %f %e %e\n" , lerp_pt ,
	    Int[0].avg , Int[0].err ) ;
    sprintf( str , "Lerp_%g.flat" , lerp_pt ) ;
    struct resampled LP = init_dist( NULL , Input->Data.x[shift].NSAMPLES ,
				     Input->Data.x[shift].restype ) ;
    equate_constant( &LP , a , Input->Data.x[shift].NSAMPLES ,
		     Input->Data.x[shift].restype ) ;
    
    write_flat_dist( Int , &LP , 1 , str ) ;

    free( LP.resampled ) ;
    
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
#ifdef STRANGE
  #ifdef CONNECTED
  const double NUM[2] = { 1 , 1 } ;
  #else
  const double NUM[2] = { -1 , -1 } ;
  #endif
#elif (defined CHARM)
  #ifdef CONNECTED
  const double NUM[2] = { 16 , 16 } ;
  #elif (defined LC)
  const double NUM[2] = { -24 , -24 } ;
  #else
  const double NUM[2] = { -16 , -16 } ;
  #endif
#else
  #ifdef SYMMETRIC
    #ifdef CONNECTED
    const double NUM[2] = { 18 , 18 } ;
    #elif defined CPLUSD
    const double NUM[2] = { 18 , -36 } ;
    #else
    const double NUM[2] = { -36 , -36 } ;
    #endif
  #else
    #ifdef CONNECTED
    const double NUM[2] = { 17 , 17 } ;
    #elif defined CPLUSD
    const double NUM[2] = { 17 , -25 } ;
    #else
      #ifdef LS // in the code I add the two contributions
      const double NUM[2] = { -5 , -5 } ;
      #elif (defined LC)
      const double NUM[2] = { -20 , -20 } ;
      #elif (defined CC)
      const double NUM[2] = { -16 , -16 } ;
      #else
      const double NUM[2] = { -25 , -25 } ;
      #endif
    #endif
  #endif
#endif
  const double DEN = 81. ;

  FILE *txtfile = fopen( "HLBL_data.txt" , "w" ) ;
    
  //const double Q = 17/81. ; // -> 25/81. for the disconnected!
  // if we do the symmetric point this factor is 18/81. & 36/81.
  // values from the g-2 paper and the ZV paper
  // https://arxiv.org/pdf/1904.03120.pdf
  // https://arxiv.org/pdf/1811.08209.pdf
#ifdef B334
  const double amu = 0.05316465193824131 ;
  const double a = 0.09929 ;
  #if (defined A654)
  const double Z = 0.69789 ;
  #elif (defined A653)
    #ifdef CHARMZV
  //const double Z = 8.415401e-01 ;
  //const double Z = 8.781319e-01 ;
  //const double Z = 9.153006e-01 ;
  //const double Z = 9.530678e-01 ;
  //const double Z = 9.914530e-01 ;
  //const double Z = 1.069983e+00 ;
  //const double Z = 1.151280e+00 ;
  const double Z = 1.322647e+00 ; // "physical Ds"
    #else
    const double Z = 0.703507 ;
    #endif
  #endif
#elif (defined B340)
  const double amu = 0.04624130625372472 ;
  const double a = 0.08636 ;
  #if (defined H101)
  //const double fpirat = 1 ;
  //const double fpirat = 1.1350853574396977 ;
  const double fpirat = 1.288225 ;
  //const double fpirat = 1.659523650625 ;
    #ifdef CHARMZV
  //const double Z = 8.270813e-01 ;
  //const double Z = 8.858063e-01 ;
  //const double Z = 9.461044e-01 ;
  //const double Z = 1.007931e+00 ;
  //const double Z = 1.071593e+00 ;
    const double Z = 1.20324 ; // g-2 PAP value
    #else
    const double Z = 0.71562 ;
    #endif
  #elif (defined U103)
  const double Z = 0.71562 ;
  #elif (defined H102) || (defined U102)
  const double Z = 0.71226 ;
  #elif (defined U102)
  const double Z = 0.71226 ;
  #elif (defined H105) || (defined U101)
  const double fpirat = 1 ;
  //const double fpirat = 1.0197996401467917 ;
  //const double fpirat = 1.0399913060435257 ;
  const double Z = 0.70908 ;
  #elif (defined C101)
  const double fpirat = 1 ;
  //const double fpirat = 0.9728698738732318 ;
  //const double fpirat = 0.9464757914901178 ;
  const double Z = 0.70717 ;
  #endif
#elif (defined B346)
  const double amu = 0.040876115324332385 ;
  const double a = 0.07634 ;
  #if (defined B450)
    #ifdef CHARMZV
    //const double Z = 8.197728e-01 ;
    //const double Z = 8.685114e-01 ;
    //const double Z = 9.184234e-01 ;
    //const double Z = 9.694756e-01 ;
    //const double Z = 1.021638e+00 ;
    const double Z = 1.12972 ; // g-2 paper
    #else
    const double Z = 0.72647 ;
    #endif
  #elif (defined D450)
  const double Z = 0.719209 ;
  #endif
#elif (defined B355)
  const double amu = 0.034407901110055004 ;
  const double a = 0.06426 ;
  #if (defined H200)
    #ifdef CHARMZV
    const double Z = 1.04843 ;
    #else
    const double Z = 0.74028 ;
    #endif
  #elif (defined N202)
    #ifdef CHARMZV
    //const double Z = 8.132188e-01 ;
    //const double Z = 8.509598e-01 ;
    //const double Z = 8.895498e-01 ;
    //const double Z = 9.288259e-01 ;
    //const double Z = 9.677787e-01 ;
    const double Z = 1.04843 ; // "physical Ds" g-2 paper
    #else
    const double Z = 0.74028 ;
    #endif
  #elif (defined N203)
    #ifdef CHARMZV
    const double Z = 1.04534 ;
    #else
    const double Z = 0.73792 ;
    #endif
  #elif (defined N200)
    #ifdef CHARMZV
    const double Z = 1.03587 ;
    #else
    const double Z = 0.73614 ;
    #endif
  #elif (defined D200)
    #ifdef CHARMZV
    const double Z = 1.03310 ;
    #else
    const double Z = 0.73429 ;
    #endif
  #else
  fprintf( stdout , "Ensemble not recognised for Z_V\n") ;
  return FAILURE ;
  #endif
#elif (defined B370)
  const double amu = 0.02667067466996327 ;
  const double a = 0.04981 ;
  #ifdef CHARMZV
  //const double Z = 8.111842e-01 ;
  //const double Z = 8.379414e-01 ;
  //const double Z = 8.650583e-01 ;
  //const double Z = 8.925378e-01 ;
  //const double Z = 9.204328e-01 ;
  const double Z = 0.97722 ; // physical Ds g-2 paper
  #else
  const double Z = 0.75909 ;
  #endif
#elif (defined B385)
  const double amu = 0.020936025425789363 ;
  const double a = 0.0391 ;
    #ifdef CHARMZV
  //const double Z = 8.126588e-01 ;
  //const double Z = 8.324016e-01 ;
  const double Z = 8.522841e-01 ;
  //const double Z = 8.724131e-01 ;
  //const double Z = 8.927380e-01 ;
    #endif
#else
  fprintf( stderr , "I do not understand your selected beta\n" ) ;
  return FAILURE ;
#endif
  const double ZV = pow(Z,PL) ;
  const double alpha_QED = 1.0/137.035999 ;
#ifdef FPIRESCALE
  const double prefac = ZV*amu*pow(4*M_PI*alpha_QED,3)/3.*2*M_PI*M_PI*fpirat ;
#else
  const double prefac = ZV*amu*pow(4*M_PI*alpha_QED,3)/3.*2*M_PI*M_PI ;
#endif

  size_t i , j , shift = 0 ;
  fprintf( txtfile , "Config is %s, %s\n" , CONFIGLABEL , BETALABEL ) ;
  fprintf( txtfile , "a = %g [fm], amu %g\n" , a , amu ) ;
  fprintf( txtfile , "Renormalisation (%g)^%d\n" , Z , PL ) ;
  fprintf( txtfile , "          |y| [fm]\tf(|y|)[fm^-1]\terror\n" ) ;
  #ifndef CPLUSD
  fprintf( txtfile , "Using %g/%g charge factor\n" , NUM[0] , DEN ) ;
  #endif
  
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    const double Q = NUM[i%2]/DEN ;
    // multiply by |x|^3*prefac
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      // x is |y| when comparing to note
      root( &Input -> Data.x[j] ) ;

      // multiply by Q-factors
      mult_constant( &Input -> Data.y[j] , Q*prefac ) ;

      #ifdef MULT_X3
        #ifdef IN_FERMI
        mult_constant( &Input -> Data.y[j] , pow( Input -> Data.x[j].avg , 3 )/a ) ;
        mult_constant( &Input -> Data.x[j] , a ) ;
        #else
        mult_constant( &Input -> Data.y[j] , pow( Input -> Data.x[j].avg , 3 ) ) ;
        #endif
      #endif
      
      fprintf( txtfile , "Integrand %e %e %e\n" ,
	       Input->Data.x[j].avg ,
	       Input -> Data.y[j].avg ,
	       Input -> Data.y[j].err ) ;
    }    
    shift += Input -> Data.Ndata[i] ;
  }

#ifdef CPLUSD
  // list is sorted so we just go through and add into shortest
  const size_t sidx = Input->Data.Ndata[0] < Input->Data.Ndata[1] ? 0 : 1 ;
  if( sidx == 0 ) {
    for( j = 0 ; j < Input->Data.Ndata[0] ; j++ ) {
      size_t k ;
      for( k = j + Input->Data.Ndata[0] ; k < Input->Data.Ndata[1] + Input->Data.Ndata[0] ; k++ ) {
	if( fabs( Input -> Data.x[k].avg - Input -> Data.x[j].avg ) < 1E-14 ) {
	  printf("CHECK %f == %f\n" , Input -> Data.x[j].avg , 
		 Input -> Data.x[k].avg ) ;
	  add( &Input -> Data.y[j] , Input -> Data.y[k] ) ;
	  equate( &Input -> Data.x[j] , Input -> Data.x[k] ) ;
	}	
      }
    }
  } else {
    for( j = Input->Data.Ndata[0] ; j < Input->Data.Ndata[0]+Input->Data.Ndata[1] ; j++ ) {
      size_t k ;
      for( k = 0 ; k < Input->Data.Ndata[0] ; k++ ) {
	if( fabs( Input -> Data.x[k].avg - Input -> Data.x[j].avg ) < 1E-14 ) {
	  printf("CHECK %f == %f\n" , Input -> Data.x[j].avg , 
		 Input -> Data.x[k].avg ) ;
	  add( &Input -> Data.y[j] , Input -> Data.y[k] ) ;
	  equate( &Input -> Data.x[j] , Input -> Data.x[k] ) ;
	}
      }
    }     
  }
#endif

  // write out integrands
  shift=0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    #ifdef PUT_ZERO
    fprintf( txtfile , "Integrand %e %e %e\n" , 0. , 0. , 0. ) ;
    #endif    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      fprintf( txtfile , "Integrand %e %e %e\n" ,
	       Input -> Data.x[j].avg ,
	       Input -> Data.y[j].avg ,
	       Input -> Data.y[j].err ) ;
    }
    shift += Input -> Data.Ndata[i] ;
  }
  #ifndef NINT_HLBL
  integrate( Input , txtfile , a ) ;
  #endif
  fclose( txtfile ) ;

  // perform a fit
  double Chi ;
#ifdef NINT_HLBL
  struct resampled *Fit = fit_and_plot_and_Nint( *Input , &Chi ) ;
#else
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
#endif
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}
