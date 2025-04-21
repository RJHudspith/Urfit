/**
   Fuck you
 */
#include "gens.h"

#include "Nder.h"

double
ftest( const struct x_desc X , const double *fparams , const size_t Npars )
{
  if( Npars > 5 ) {
    return fparams[2]*( 1 + fparams[3]*(X.X - fparams[4]) ) ;
  } else if( Npars > 2 ) {
    return fparams[1]*( 1 + fparams[3]*(X.X - fparams[4]) ) ;
  } else {
    return fparams[0]*( 1 + fparams[3]*(X.X - fparams[4]) ) ;
  }
}

void
test_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;
    f[i] = ftest( X , fparams , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
test_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {

    const struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			      DATA -> N , DATA -> M } ;

    if( i > 5 ) {
      df[0][i] = 0 ;
      df[1][i] = 0 ;
      df[2][i] = ( 1 + fparams[3]*(X.X - fparams[4]) ) ;
      df[3][i] = fparams[2]*(X.X-fparams[4]) ;
      df[4][i] = -fparams[2]*fparams[3] ;
    } else if( i > 2 ) {
      df[0][i] = 0 ;
      df[1][i] = ( 1 + fparams[3]*(X.X - fparams[4]) ) ;
      df[2][i] = 0.0 ;
      df[3][i] = fparams[1]*(X.X-fparams[4]) ;
      df[4][i] = -fparams[1]*fparams[3] ;
    } else {
      df[0][i] = ( 1 + fparams[3]*(X.X - fparams[4]) ) ;
      df[1][i] = 0 ;
      df[2][i] = 0.0 ;
      df[3][i] = fparams[0]*(X.X-fparams[4]) ;
      df[4][i] = -fparams[0]*fparams[3] ;
    }
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
test_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
test_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 0 ;
  fparams[1] = 0 ;
  fparams[2] = 0 ;
  fparams[3] = 0 ;

	 
  return ;
}
