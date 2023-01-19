/**
   @file cornell.c
   @brief cornell fit

   big global fit to A exp( -( B/r + C + D*r)*t )
   
 */
#include "gens.h"

static const double tmap[ ] = { 11.5 , 12.5 , 13.5 , 14.5 , 15.5, 16.5} ;
//static const double tmap[ ] = { 11, 12 , 13, 14 , 15 , 16 , 17 } ;

#define PAR3

double
fcornellv2( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef PAR3
  return
    fparams[0]*exp( -( fparams[3]/X.X + fparams[2] + fparams[1]*X.X)*(tmap[Npars]) ) ;
#else
  return
    fparams[0]*exp( -( M_PI/(12*X.X) + fparams[2] + fparams[1]*X.X)*(tmap[Npars]) ) ;
#endif
}

void
cornellv2_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fcornellv2( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

void
cornellv2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double r = DATA -> x[i] ;
    const double t = tmap[ DATA -> map[i].bnd ] ; //+0.5;

#ifdef PAR3
    const double ee = exp( -(fparams[3]/r + fparams[2] + fparams[1]*r)*t ) ;
    // derivatives wrt the various fit params
    df[0][i] = ee ;
    df[3][i] = -(t/r)*fparams[0]*ee ;
    df[2][i] = -(t  )*fparams[0]*ee ;
    df[1][i] = -(t*r)*fparams[0]*ee ;
#else
    const double ee = exp( -(M_PI/(12*r) + fparams[2] + fparams[1]*r)*t ) ;
    // derivatives wrt the various fit params
    df[0][i] = ee ;
    df[1][i] = -(t*r)*fparams[0]*ee ;
    df[2][i] = -(t  )*fparams[0]*ee ;
#endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
cornellv2_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
cornellv2_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit )
{
  fparams[0] = 1000 ;
  fparams[1] = 1.0 ;
  fparams[2] = 0.1 ;
}
