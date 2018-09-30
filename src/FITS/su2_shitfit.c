/**
   @file poly.c
   @brief polynomial fit
 */
#include "gens.h"

#include "Nder.h"

double
fsu2_shitfit( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return sqrt( fparams[0] + fparams[1]*X.X ) ;
  //return fparams[0] + fparams[1]*X.X*( 1 + fparams[2]*log(X.X) ) ;
  //return ( fparams[0] + fparams[1]*X.X + fparams[2]*X.X*X.X)*exp(-fparams[3]*X.X) + ( fparams[4] + fparams[5]*X.X)*(1-exp(-fparams[6]*X.X)) ;
}

void
su2_shitfit_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , p ; 
  for (i = 0; i < DATA -> n ; i++) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    double par[ DATA -> Npars ] ;
    for( p = 0 ; p < DATA -> Npars ; p++ ) {
      par[ p ] = fparams[ DATA -> map[i].p[p] ] ;
    }
    f[i] = fsu2_shitfit( X , par , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
su2_shitfit_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i  ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    size_t j ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = Nder( fsu2_shitfit , X ,
					   DATA -> map[i].bnd ,
					   fparams , j , DATA -> Npars ) ;
    }
    
  }
  return ;
}

// second derivatives are all zero
void
su2_shitfit_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// polynomial guesses
void
su2_shitfit_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit )
{

  return ;
}
