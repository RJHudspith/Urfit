#include "gens.h"
#include "Nder.h"

enum{ udbb , udcb , lsbb , lscb } ;

#define INVERSE

double
fudcb_heavy( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double r = X.X ;
#ifdef INVERSE
  switch( Npars ) {
  case udbb : return fparams[0]/(2*r) + fparams[1] + fparams[2]*2.*r + 22.8*(r) ;
  case udcb :
    if( r < 1 ) {
      return fparams[0]/(1+r) + fparams[1] + fparams[2]*(1+r) + 34.1 - 11.4*r ;
    } else {
      return fparams[0]/(1+r) + fparams[1] + fparams[2]*(1+r) + 34.1*r - 11.4 ;
    }
  case lsbb : return fparams[0]/(2*r) + fparams[4] + fparams[3]*2.*r + 24.4*r ;
  case lscb :
    if( r < 1 ) {
      return fparams[0]/(1+r) + fparams[4] + fparams[3]*(1+r) + 34.1 - 11.9*r ;
    } else {
      return fparams[0]/(1+r) + fparams[4] + fparams[3]*(1+r) + 35.7*r - 11.4 ;
    }
  }
#else
  switch( Npars ) {
  case udbb : return fparams[0]*r/2 + fparams[1] + fparams[2]*2./r + 22.8/(r) ;
  case udcb :
    if( r > 1 ) {
      return fparams[0]*r/(1+r) + fparams[1] + fparams[2]*(1+1./r) + 34.1 - 11.4/r ;
    } else {
      return fparams[0]*r/(1+r) + fparams[1] + fparams[2]*(1+1./r) + 34.1/r - 11.4 ;
    }
  case lsbb : return fparams[0]*r/2 + fparams[4] + fparams[3]*2./r + 24.4/r ;
  case lscb :
    if( r > 1 ) {
      return fparams[0]*r/(1+r) + fparams[4] + fparams[3]*(1+1./r) + 34.1 - 11.9/r ;
    } else {
      return fparams[0]*r/(1+r) + fparams[4] + fparams[3]*(1+1./r) + 35.7/r - 11.4 ;
    }
  }
#endif
  return 0 ;
}

void
udcb_heavy_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = fudcb_heavy( X , fparams , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
udcb_heavy_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    // initialise everything to zero
    for( j = 0 ; j < 5 ; j++ ) df[j][i] = 0.0 ;

    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    size_t j ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( fudcb_heavy , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }
  }
}

void
udcb_heavy_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
udcb_heavy_guesses( double *fparams ,
		    const struct data_info Data ,
		    const struct fit_info Fit )
{
  fparams[0] = -11 ;
  fparams[1] = -200 ;
  fparams[2] = 300 ;
  fparams[3] = -200 ;
  fparams[4] = 300 ;
  return ;
}
