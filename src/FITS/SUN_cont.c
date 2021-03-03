/**
   multiple fit to SUN continuum limit
 */
#include "gens.h"

#include "Nder.h"

#define FIT7

double
fSUN_cont( const struct x_desc X , const double *fparams , const size_t Npars )
{
#ifdef FIT1
  return fparams[0]*( 1 + fparams[1]*X.X + fparams[2]*X.X/(X.LT*X.LT) ) ;
#elif defined FIT2
  return fparams[0]*( 1 + fparams[1]*X.X + fparams[2]*X.X/(X.LT*X.LT)
		      + fparams[3]*X.X*X.X/(X.LT*X.LT)) ;
#elif defined FIT3
  return fparams[0]*( 1 + fparams[1]*X.X + fparams[2]*X.X/(X.LT*X.LT)
		      + fparams[3]/(X.LT*X.LT)) ;
#elif defined FIT4
  return fparams[0]*( 1 + fparams[1]*X.X + fparams[2]/(X.LT*X.LT)) ;
#elif defined FIT5
  return fparams[0]*( 1 + fparams[1]*X.X + fparams[2]/(X.LT*X.LT)
		      + fparams[3]*X.X*X.X/(X.LT*X.LT*X.LT*X.LT) ) ;
#elif defined FIT6
  return fparams[0]*( 1 + fparams[1]*X.X + fparams[2]/(X.LT*X.LT)
		      + fparams[3]*X.X*X.X/(X.LT*X.LT*X.LT*X.LT)
		      + fparams[4]/(X.LT*X.LT*X.LT)
		      ) ;
#elif defined FIT7
  return fparams[0]*( 1 + fparams[1]*X.X
		      + fparams[2]*X.X*X.X/(X.LT*X.LT*X.LT*X.LT*X.LT*X.LT) ) ;
#elif defined FIT8
  return fparams[0]*( 1 + fparams[1]*X.X
		      + fparams[2]*X.X*X.X/(X.LT*X.LT*X.LT*X.LT*X.LT) ) ;
#elif defined FIT9
  return fparams[0]*( 1 + fparams[1]*X.X
		      + fparams[2]*X.X*X.X*X.X/(X.LT*X.LT*X.LT*X.LT*X.LT*X.LT) ) ;  
#else
  return fparams[0]*( 1 + fparams[3]/(X.LT*X.LT) + fparams[1]*X.X + fparams[2]*X.X/(X.LT*X.LT) ) ;
#endif
}

void
SUN_cont_f( double *f , const void *data , const double *fparams )
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
    f[i] = fSUN_cont( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
SUN_cont_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i ;  
  for( i = 0 ; i < DATA -> n ; i++ ) {
    const double x  = DATA -> x[i] ;
    const size_t NC = DATA -> LT[i] ;
    const double p0 = fparams[ DATA -> map[i].p[0] ] ;
    const double p1 = fparams[ DATA -> map[i].p[1] ] ;
    const double p2 = fparams[ DATA -> map[i].p[2] ] ;
#ifdef FIT1
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2*x/(NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x/(NC*NC);
#elif defined FIT2
    const double p3 = fparams[ DATA -> map[i].p[3] ] ;
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2*x/(NC*NC)+p3*x*x/(NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x/(NC*NC);
    df[ DATA -> map[i].p[3] ][i] = p0*x*x/(NC*NC);
#elif defined FIT3
    const double p3 = fparams[ DATA -> map[i].p[3] ] ;
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2*x/(NC*NC)+p3/(NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x/(NC*NC);
    df[ DATA -> map[i].p[3] ][i] = p0/(NC*NC);
#elif defined FIT4
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2/(NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0/(NC*NC);
#elif defined FIT5
    const double p3 = fparams[ DATA -> map[i].p[3] ] ;
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2/(NC*NC)+p3*x*x/(NC*NC*NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0/(NC*NC);
    df[ DATA -> map[i].p[3] ][i] = p0*x*x/(NC*NC*NC*NC);
#elif defined FIT6
    const double p3 = fparams[ DATA -> map[i].p[3] ] ;
    const double p4 = fparams[ DATA -> map[i].p[3] ] ;
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2/(NC*NC)
				    +p3*x*x/(NC*NC*NC*NC)
				    +p4/(NC*NC*NC*NC)
				    ) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0/(NC*NC);
    df[ DATA -> map[i].p[3] ][i] = p0*x*x/(NC*NC*NC*NC);
    df[ DATA -> map[i].p[4] ][i] = p0/(NC*NC*NC);
#elif defined FIT7
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2*x*x/(NC*NC*NC*NC*NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x*x/(NC*NC*NC*NC*NC*NC);
#elif defined FIT8
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2*x*x/(NC*NC*NC*NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x*x/(NC*NC*NC*NC*NC);
#elif defined FIT9
    df[ DATA -> map[i].p[0] ][i] = (1 + p1*x+p2*x*x*x/(NC*NC*NC*NC*NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x*x*x/(NC*NC*NC*NC*NC*NC);
#else
    const double p3 = fparams[ DATA -> map[i].p[3] ] ;
    df[ DATA -> map[i].p[0] ][i] = (1+p3/(NC*NC) + p1*x+p2*x/(NC*NC)) ;
    df[ DATA -> map[i].p[1] ][i] = p0*x ;
    df[ DATA -> map[i].p[2] ][i] = p0*x/(NC*NC);
    df[ DATA -> map[i].p[3] ][i] = p0/(NC*NC);
#endif
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
SUN_cont_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
SUN_cont_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit )
{
  fparams[0] = 7E-4 ;
  fparams[1] = 0.1 ;
  fparams[2] = -1 ;
#if (defined FIT5)
  fparams[3] = 1 ;
#elif (defined FIT6)
  fparams[3] = 1 ;
  fparams[4] = 1 ;
#endif
}
