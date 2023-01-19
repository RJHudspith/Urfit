/**
   multiple fit to SUN continuum limit
 */
#include "gens.h"

#include "Nder.h"

#define V2

const double NMAP[4] = { 3,4,5,6 } ;

double
fLargeNB( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef V1 // not terrible
  return (fparams[0] + ( fparams[1] + (fparams[2] + fparams[3]/X.X)/X.X)/X.X) *( NMAP[Npars] ) ;
#elif (defined V2)
  return (fparams[0] + ( fparams[1] + (fparams[2] + fparams[3]/X.X)/X.X)/X.X) *( 1+fparams[4]*NMAP[Npars] ) ;
#else
  return (fparams[0] + fparams[1]/X.X + fparams[2]*exp( -fparams[3]/X.X ) )*( NMAP[Npars] ) ;
#endif
}

void
LargeNB_f( double *f , const void *data , const double *fparams )
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
    f[i] = fLargeNB( X , p , DATA->map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
LargeNB_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;  
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double x  = DATA -> x[i] ;

#ifdef V1
    const double N = NMAP[ DATA->map[i].bnd ] ;
    df[ DATA->map[i].p[0] ][i] = N ;
    df[ DATA->map[i].p[1] ][i] = N/x ;
    df[ DATA->map[i].p[2] ][i] = N/(x*x) ;
    df[ DATA->map[i].p[3] ][i] = N/(x*x*x) ;
#elif (defined V2)
    const double Nm = NMAP[ DATA->map[i].bnd ] ;
    const double N = 1.+fparams[ DATA->map[i].p[4]] * Nm ;
    const double P0 = fparams[ DATA->map[i].p[0] ] ;
    const double P1 = fparams[ DATA->map[i].p[1] ] ;
    const double P2 = fparams[ DATA->map[i].p[2] ] ;
    const double P3 = fparams[ DATA->map[i].p[3] ] ;
    
    df[ DATA->map[i].p[0] ][i] = N ;
    df[ DATA->map[i].p[1] ][i] = N/x ;
    df[ DATA->map[i].p[2] ][i] = N/(x*x) ;
    df[ DATA->map[i].p[3] ][i] = N/(x*x*x) ;
    df[ DATA->map[i].p[4] ][i] = Nm*( P0 + P1/x + P2/(x*x) + P3/(x*x*x) ) ;
#else
    const double N = NMAP[ DATA->map[i].bnd ] ;
    df[ DATA->map[i].p[0] ][i] = (N) ;
    df[ DATA->map[i].p[1] ][i] = (N)/x ;
    df[ DATA->map[i].p[2] ][i] = (N)*exp( -fparams[ DATA->map[i].p[3] ]/x ) ;
    df[ DATA->map[i].p[3] ][i] = (-N*fparams[ DATA->map[i].p[2]]/x)*exp( -fparams[ DATA->map[i].p[3] ]/x ) ;
#endif
  }
  return ;
}

void
LargeNB_linmat( double **U ,
		const void *data ,
		const size_t N ,
		const size_t M ,
		const size_t Nlogic )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    // set matrix to zero
    for( j = 0 ; j < Nlogic ; j++ ) {
      U[i][j] = 0.0 ;
    }
    const double x  = DATA -> x[i] ;
    const double N = NMAP[ DATA->map[i].bnd ] ;

    U[i][ DATA -> map[i].p[0] ] = (N) ;
    U[i][ DATA -> map[i].p[1] ] = (N)/x ;
    U[i][ DATA -> map[i].p[2] ] = (N)/(x*x) ;
    U[i][ DATA -> map[i].p[3] ] = (N)/(x*x*x) ;
  }
}

// second derivatives? Will we ever use them - J?
void
LargeNB_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
LargeNB_guesses( double *fparams ,
		 const struct data_info Data ,
		 const struct fit_info Fit )
{
}
