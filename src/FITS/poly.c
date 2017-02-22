/**
   exponential Model f = ( p(x_i) - y_i )
 */
#include "gens.h"

#include "gls_bootfit.h"

double
fpoly( const struct x_desc X , const double *fparams , const size_t Npars )
{
  size_t i ;
  if( X.N < 1 ) return fparams[0] ;
  register double poly = X.X * fparams[ X.N ] ;
  for( i = X.N-1 ; i > 0 ; i-- ) {
    poly = X.X * ( fparams[i] + poly ) ;
  }
  return poly + fparams[0] ;
}

void
poly_f( double *f , const void *data , const double *fparams )
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
    f[i] = fpoly( X , par , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
poly_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double xloc = 1.0 ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      df[ DATA -> map[i].p[j] ][i] = xloc ;
      xloc *= DATA->x[i] ;
    }
  }
  return ;
}

// second derivatives are all zero
void
poly_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// compute the linearised matrix for the GLS
void
poly_linmat( double **U ,
	     const void *data ,
	     const size_t Nparam ,
	     const size_t Nlogic )
{
  struct data *Data = (struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < Data -> n ; i++ ) {
    for( j = 0 ; j < Nlogic ; j++ ) {
      U[i][j] = 0.0 ;
    }
    register double x0 = 1 ;
    for( j = 0 ; j < Nparam ; j++ ) {
      U[i][ Data -> map[i].p[j] ] = x0 ;
      x0 *= Data -> x[i] ;
    }
  }
  return ;
}

// polynomial guesses
void
poly_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit )
{
  double chisq = 0.0 ;
  size_t i ;

  single_gls( fparams , &chisq , Data , Fit , 1 , true ) ;
  
  // tell us about the guesses
  printf( "\n" ) ;
  for( i = 0 ; i < Fit.Nlogic ; i++ ) {
    printf( "[GUESS] Fit param guess %zu -> %f \n" , i , fparams[i] ) ; 
  }
  printf( "\n" ) ;
  
  
  return ;
}
