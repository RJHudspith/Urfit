/**
   @file poles.c
   @brief fit some poles
 */
#include "gens.h"

#include "Nder.h"

//#define NOPOLE

double
fpoles( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef NOPOLE
  return fparams[0] + X.X * fparams[1] ;
#else

  switch( Npars ) {
  case 0 : return  fparams[0] + X.X * ( fparams[3] + X.X * fparams[1] ) ;
  case 1 : return  fparams[0] + X.X * ( fparams[5] + X.X * fparams[1] ) ;
  case 2 : return -fparams[0] + X.X * ( fparams[4] + X.X * fparams[2] ) ;
  case 3 : return -fparams[0] + X.X * ( fparams[6] + X.X * fparams[2] ) ;
  }
#endif
}

void
poles_f( double *f , const void *data , const double *fparams )
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
    f[i] = fpoles( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
poles_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const double X = DATA -> x[i] ;
    const double XX = X * X ;
    
    switch( DATA -> map[i].bnd ) {
    case 0 :
      df[0][i] = 1.0 ;
      df[1][i] = XX ;
      df[2][i] = 0.0 ;
      df[3][i] = X ;
      df[4][i] = 0.0 ;
      df[5][i] = 0.0 ;
      df[6][i] = 0.0 ;
      break ;
    case 1 :
      df[0][i] = 1.0 ;
      df[1][i] = XX ;
      df[2][i] = 0.0 ;
      df[3][i] = 0.0 ;
      df[4][i] = 0.0 ;
      df[5][i] = X ;
      df[6][i] = 0.0 ;
      break ;
    case 2 :
      df[0][i] = -1.0 ;
      df[1][i] = 0.0 ;
      df[2][i] = XX ;
      df[3][i] = 0.0 ;
      df[4][i] = X ;
      df[5][i] = 0.0 ;
      df[6][i] = 0.0 ;
      break ;
    case 3 :
      df[0][i] = -1.0 ;
      df[1][i] = 0.0 ;
      df[2][i] = XX ;
      df[3][i] = 0.0 ;
      df[4][i] = 0.0 ;
      df[5][i] = 0.0 ;
      df[6][i] = X ;
      break ;
    }

  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
poles_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
poles_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
#ifdef NOPOLE
  fparams[0] = 0.1 ;
  fparams[1] = 1 ;
  fparams[2] = 0.1 ;
#else
  fparams[0] = 0.008 ;
  fparams[2] = -0.6 ;
  fparams[1] = -0.88 ;
  fparams[3] = -0.1 ;
  fparams[4] = -0.88 ;
  fparams[5] = -0.1 ;
#endif
}
