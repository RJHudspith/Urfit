/**
   @file qslab2.c
   @brief fit the slab formula Y = \chi * X ( 1 - X )
 */
#include "gens.h"
#include "Nder.h"

#define WITHOUT_FSE

double
fqslab2( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef WITHOUT_FSE
  return fparams[0] * X.X*( 1. - X.X ) - fparams[1] ;
#else
  return fparams[0] * X.X*( 1. - X.X/X.LT ) 
    - X.LT*fparams[1]*( 1.
			- exp( -fparams[2]*X.X )
			- exp( -fparams[2]*(X.LT-X.X ) )
			+ exp( -fparams[2]*(X.LT) )
		   ) ;
#endif
}

void
qslab2_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    f[i] = fqslab2( X , p , DATA -> Npars ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
qslab2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double x = DATA->x[i] ;
    const double LT = DATA -> LT[i] ;
    const double xp = x ;

    df[ DATA -> map[i].p[0] ][ i ] = x*(1-xp) ;

#ifdef WITHOUT_FSE
    df[ DATA -> map[i].p[1] ][ i ] = -1 ;
#else
    df[ DATA -> map[i].p[1] ][ i ] = 
      -LT*(1.- exp( -fparams[2] * x )
	   -exp( -fparams[2] *(LT-x))
	   +exp( -fparams[2] *LT)
	) ;
    df[ DATA -> map[i].p[2] ][ i ] =
      fparams[1]*LT*(
		     -x*exp( -fparams[2] * x )
		     -(LT-x)*exp( -fparams[2] *(LT-x) )
		     +LT*exp( -fparams[2] *(LT))
		     ) ;
#endif
  }
  return ;
}

// second derivatives
void
qslab2_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// guesses using the data, we don't need to be too accurate
void
qslab2_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{ 
  fparams[0] = 7E-4 ;
  fparams[1] = 1E-6 ;
  //fparams[3] = 0.3 ;
}
