/**
   multiple fit to c4 and c7
 */
#include "gens.h"

/*
static const double c4map[ 12 ] = { 1 , 1.23 , 1.4 ,
				    1 , 1.23 , 1.4 ,
				    1 , 1.23 , 1.4 ,
				    1 , 1.23 , 1.4 } ;
*/

static const double c7map[ 12 ] = { 1 , 1 , 1 ,
				    1.23 , 1.23 , 1.23 ,
				    1.4 , 1.4 , 1.4 ,
				    0 , 0 , 0 } ;

double
fc4c7( const struct x_desc X , const double *fparams , const size_t Npars )
{
  return fparams[0] + X.X*X.X*(fparams[1])  + c7map[Npars]*(fparams[2] ) ;
}

void
c4c7_f( double *f , const void *data , const double *fparams )
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
    f[i] = fc4c7( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
c4c7_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;  
  for( i = 0 ; i < DATA -> n ; i++ ) {    
    const double c4 = DATA -> x[i] ;
    const double c7 = c7map[i] ;
    df[ DATA -> map[ i ].p[0] ][i] = 1 ;
    df[ DATA -> map[ i ].p[1] ][i] = c4*c4 ;
    df[ DATA -> map[ i ].p[2] ][i] = c7 ;
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
c4c7_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
c4c7_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit )
{
  fparams[0] = 0.0 ;
  fparams[1] = 0.04 ;
  fparams[2] = 0. ;
  fparams[3] = -0.03 ;
}
