#include "gens.h"

#include "fit_and_plot.h"
#include "init.h"
#include "momenta.h"
#include "resampled_ops.h"
#include "stats.h"
#include "write_flat.h"

static double
b2sgl( const double Q4 , const double Q2 )
{
  return (Q4+3*Q2*Q2)/(12*Q2) ;
}

static double
b4sgl( const double Q6 , const double Q4 , const double Q2 )
{
  return (Q6 - 15*Q2*Q4 + 30*Q2*Q2*Q2)/(180*360.*Q2) ;
}

static void
compute_b2( struct resampled *b2 ,
	    const struct resampled Q4 ,
	    const struct resampled Q2 )
{
  size_t i ;
  for( i = 0 ; i < Q2.NSAMPLES ; i++ ) {
    b2->resampled[i] = b2sgl( Q4.resampled[i] , Q2.resampled[i] ) ;
  }
  b2->avg = b2sgl( Q4.avg , Q2.avg ) ;
  compute_err( b2 ) ;
}

static void
compute_b4( struct resampled *b4 ,
	    const struct resampled Q6 ,
	    const struct resampled Q4 ,
	    const struct resampled Q2 )
{
  size_t i ;
  for( i = 0 ; i < Q2.NSAMPLES ; i++ ) {
    b4->resampled[i] = b4sgl( Q6.resampled[i] ,
			      Q4.resampled[i] ,
			      Q2.resampled[i] ) ;
  }
  b4->avg = b4sgl( Q6.avg , Q4.avg , Q2.avg ) ;
  compute_err( b4 ) ;
}

// topological moments of Q^2
int
Qmoments( struct input_params *Input )
{
  size_t j ;

  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    printf( "[QMOM] QMOM_%zu %f %f \n" , j ,
	    Input -> Data.y[j].avg ,
	    Input -> Data.y[j].err ) ;
  }
  
  for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
    if( j!=2 ) {
      divide( &Input -> Data.y[j] , Input -> Data.y[2] ) ;
    }
  }

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  return SUCCESS ;
}

// using the moments of the Q^2 distribution we compute the cumulants
int
CumFromMom( struct input_params *Input )
{
  struct resampled tmp = init_dist( NULL ,
				    Input -> Data.y[0].NSAMPLES ,
				    Input -> Data.y[0].restype ) ;

  compute_b4( &tmp , Input->Data.y[6] , Input->Data.y[4] , Input->Data.y[2] ) ;

  printf( "Q2 Q4 Q6 :: %f %f %f\n" , Input -> Data.y[0].avg , Input -> Data.y[2].avg , Input -> Data.y[4].avg ) ;

  
  printf( "Result b4 :: %e +/- %e\n" , tmp.avg , tmp.err ) ;
  
  compute_b2( &tmp , Input->Data.y[4] , Input->Data.y[2] ) ;
  printf( "Result b2 :: %e +/- %e\n" , tmp.avg , tmp.err ) ;

  printf( "Result chi :: %e +/- %e\n" ,
	  Input->Data.y[0].avg , Input->Data.y[0].err ) ;
  
  return SUCCESS ;
}


// topological moments of Q^2
int
TraditionalQ( struct input_params *Input )
{
  struct resampled tmp = init_dist( NULL ,
				    Input -> Data.y[0].NSAMPLES ,
				    Input -> Data.y[0].restype ) ;

  compute_b4( &tmp , Input->Data.y[6] , Input->Data.y[4] , Input->Data.y[2] ) ;
  printf( "Result b4 :: %e +/- %e\n" , tmp.avg , tmp.err ) ;

  compute_b2( &tmp , Input->Data.y[4] , Input->Data.y[2] ) ;
  printf( "Result b2 :: %e +/- %e\n" , tmp.avg , tmp.err ) ;

  equate( &tmp , Input -> Data.y[2] ) ;
  printf( "Result Chi :: %e +/- %e\n" , Input->Data.y[2].avg ,
	  Input->Data.y[2].err ) ;
  
  return SUCCESS ;
}
