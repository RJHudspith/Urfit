/**
   @file qsusc_su2.c
   @brief fit the slab formula Y = ( C + \chi*t ) * ( 1 + A / L_t ) * ( m - m_0 )
 */
#include "gens.h"
#include "Nder.h"

// have the L_t = 32 first and then the L_t = 48
const double masses[] = { -0.845 , -0.855 , -0.865 , -0.875 ,
			  -0.845 , -0.855 , -0.865 , -0.875 } ;

double
fqsusc_su2( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return ( fparams[0] * ( X.X - fparams[1] ) ) * ( 1 + fparams[2] / X.LT ) ;
}

void
qsusc_su2_f( double *f , const void *data , const double *fparams )
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
    f[i] = fqsusc_su2( X , p , DATA -> map[i].bnd ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
qsusc_su2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {


    const double x = DATA -> x[i] ;
    const double LT = DATA -> LT[i] ;
    const size_t bnd = DATA -> map[i].bnd ;

    //return ( fparams[0] * ( X.X - fparams[1] ) ) * ( 1 + fparams[2] / X.LT ) ;
    df[DATA -> map[i].p[0]][i] = ( x - fparams[DATA -> map[i].p[1]] ) * ( 1 + fparams[DATA -> map[i].p[2]] / LT ) ;
    df[DATA -> map[i].p[1]][i] = -fparams[DATA -> map[i].p[0]] * ( 1 + fparams[DATA -> map[i].p[2]] / LT ) ;
    df[DATA -> map[i].p[2]][i] = ( fparams[DATA -> map[i].p[0]] * ( x - fparams[DATA -> map[i].p[1]] ) ) / LT ;
  }
  return ;
}

void
qsusc_su2_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

// guesses using the data, we don't need to be too accurate
void
qsusc_su2_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit )
{ 
  fparams[0] = -10 ;
  fparams[1] = -0.9 ;
  fparams[2] = 10 ;
  //fparams[3] = -0.901 ;
  //fparams[4] = 10 ;
  //fparams[5] = 0.0 ;
}
