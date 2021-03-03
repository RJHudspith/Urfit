/**
   multiple fit to c4 and c7
 */
#include "gens.h"

#include "Nder.h"

//#define LINEARA

static const size_t fparam_map[12] = { 2 , 2 , 3 , 4 , 4 , 5 ,
				       2 , 2 , 3 , 4 , 4 , 5 } ;

double
fHLBL_cont( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const size_t idx = fparam_map[ Npars ] ;
#ifdef LINEARA
  return fparams[0] + fparams[1]*fparams[idx] ;
#else
  return fparams[0] + fparams[1]*fparams[idx]*fparams[idx] ;
#endif
}

void
HLBL_cont_f( double *f , const void *data , const double *fparams )
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
    f[i] = fHLBL_cont( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
HLBL_cont_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;  
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const size_t idx = fparam_map[ i ] ;
    df[ DATA -> map[i].p[0] ][i] = 1 ;
#ifdef LINEARA
    df[ DATA -> map[i].p[1] ][i] = fparams[idx] ;
    df[ idx ][i] = fparams[DATA -> map[i].p[1]] ;
#else
    df[ DATA -> map[i].p[1] ][i] = fparams[idx]*fparams[idx] ;
    df[ idx ][i] = 2*fparams[DATA -> map[i].p[1]]*fparams[idx] ;
#endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
HLBL_cont_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
HLBL_cont_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit )
{
}
