/**
   @file hvp_pade.c
   @brief fit a pade to the hvp
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "fit_chooser.h"
#include "ffunction.h"
#include "make_xmgrace.h"
#include "resampled_ops.h"
#include "stats.h"

void
pade_derivative( struct resampled *adler ,
		 const double q2 , 
		 const struct resampled *fit ,
		 const struct input_params *Input )
{
  // do the individual bootstraps
  size_t j ;
  
  for( j = 0 ; j < fit[0].NSAMPLES ; j++ ) {
    double denom = 0.0 , num = 0.0 ;
    size_t n , m ;
    // left hand side
    for( n = 0 ; n < Input -> Fit.N ; n++ ) {
      num += fit[ n ].resampled[j] * pow( q2 , n ) ;
    }
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      denom += fit[ Input -> Fit.N + m ].resampled[j] * pow( q2 , (m+1) ) ;
    }
    denom += 1.0 ;
    // left hand side
    double left_sum = 0.0 , right_sum = 0.0 ;
    for( n = 1 ; n < Input -> Fit.N ; n++ ) {
      left_sum += n * fit[ n ].resampled[j] * pow( q2 , n ) ;
    }
    left_sum = left_sum / denom ;
    // right hand side
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      right_sum += (m+1) * fit[ Input -> Fit.N + m ].resampled[j] * pow( q2 , m+1 ) ;
    }
    right_sum = right_sum * num / ( denom * denom ) ;

    //printf( "%f -> %f %f \n" , num/denom , left_sum , right_sum ) ;
    adler -> resampled[j] = left_sum - right_sum ;
  }

  // do the average
  {
    double denom = 0.0 , num = 0.0 ;
    size_t n,m ;
    // left hand side
    for( n = 0 ; n < Input -> Fit.N ; n++ ) {
      num += fit[ n ].avg * pow( q2 , n ) ;
    }
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      denom += fit[ Input -> Fit.N + m ].avg * pow( q2 , (m+1) ) ;
    }
    denom += 1.0 ;
    // left hand side
    double left_sum = 0.0 , right_sum = 0.0 ;
    for( n = 1 ; n < Input -> Fit.N ; n++ ) {
      left_sum += n * fit[ n ].avg * pow( q2 , n ) ;
    }
    left_sum = left_sum / denom ;
    // right hand side
    for( m = 0 ; m < Input -> Fit.M ; m++ ) {
      right_sum += (m+1) * fit[ Input -> Fit.N + m ].avg * pow( q2 , m+1 ) ;
    }
    right_sum = right_sum * num / ( denom * denom ) ;

    //printf( "%f -> %f %f \n" , num/denom , left_sum , right_sum ) ;
    adler -> avg = left_sum - right_sum ;
  }
  
  compute_err( adler ) ;
  mult_constant( adler , 12*M_PI*M_PI*5./9. ) ;
  
  return ;
}

// numerical derivative
void
pade_derivative2( struct resampled *adler ,
		  const double q2 , 
		  const struct resampled *fit ,
		  const struct input_params *Input ,
		  const size_t shift )
{
  // set up the fit again
  struct fit_descriptor fdesc = init_fit( Input -> Data , Input -> Fit ) ;
  double fparams[ fdesc.Nparam ] ;
  size_t j , p ;
  const double h = 1E-6 ;

  // x+h and x-h
  struct x_desc xph = { q2 + h , 0 , Input -> Fit.N , Input -> Fit.M } ;
  struct x_desc xmh = { q2 - h , 0 , Input -> Fit.N , Input -> Fit.M } ;

  for( j = 0 ; j < fit[0].NSAMPLES ; j++ ) {
    // evaluate the fitfunc
    for( p = 0 ; p < fdesc.Nparam ; p++ ) {
      fparams[ p ] = fit[p].resampled[j] ;
    }
    adler -> resampled[j] =
      ( fdesc.func( xph , fparams , 0 ) -
	fdesc.func( xmh , fparams , 0 ) ) / ( 2 * h ) ;
  }

  // and the average
  for( p = 0 ; p < fdesc.Nparam ; p++ ) {
    fparams[ p ] = fit[p].avg ;
  }
  adler -> avg = ( fdesc.func( xph , fparams , 0 ) -
		   fdesc.func( xmh , fparams , 0 ) ) / ( 2 * h ) ;

  compute_err( adler ) ;
  mult_constant( adler , q2*12*M_PI*M_PI*5./9. ) ;

  // free the fitfunction
  free_ffunction( &fdesc.f , fdesc.Nlogic ) ;

  return ;
}

int
fit_hvp( struct input_params *Input )
{
  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  // compute the derivative q^2 d/dq^2
  struct resampled adler ;
  adler = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;

  make_xmgrace_graph( "der_coarse.agr" , "Q\\S2\\N [GeV\\S2\\N]" , "D(Q\\S2\\N)" ) ;
  
  size_t i = 0 ;
  const size_t Ndata = (size_t)( ( Input -> Traj[i].Fit_High - Input -> Traj[i].Fit_Low ) / 0.1 + 0.5 ) ;
  double *Q = malloc( Ndata * sizeof( double ) ) ;
  double *Ave = malloc( Ndata * sizeof( double ) ) ;
  double *Hi = malloc( Ndata * sizeof( double ) ) ;
  double *Lo = malloc( Ndata * sizeof( double ) ) ;

  size_t idx = 0 ;
  double q2 ;
  for( q2 = Input -> Traj[i].Fit_Low ; q2 < Input -> Traj[i].Fit_High ; q2 += 0.1 ) {

    pade_derivative2( &adler , q2 , fit , Input , i*Input -> Fit.Nparam ) ;

    Q[idx] = q2 ; Ave[idx] = adler.avg ; Hi[idx] = adler.err_hi ; Lo[idx] = adler.err_lo ;
    
    printf( "ADLER :: %f %f %f \n" , q2 , adler.avg , adler.err ) ;
    idx++ ;
  }

  draw_line( Q , Hi , Ndata ) ;
  draw_line( Q , Ave , Ndata ) ;
  draw_line( Q , Lo , Ndata ) ;

  free( Q ) ; free( Ave ) ; free( Hi ) ; free( Lo ) ;

  free( adler.resampled ) ;
  
  size_t k ;
  for( k = 0 ; k < Input -> Fit.Nlogic ; k++ ) {
    free( fit[k].resampled ) ;
  }
  free( fit ) ;
  
  return SUCCESS ;
}



#if 0

#include "gen_ders.h"

int
fit_hvp( struct input_params *Input )
{
  // compute the generic derivatives of the data
  size_t i , j , k , l , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    // set the derivatives
    struct resampled *der = malloc( Input -> Data.Ndata[i] *
				    sizeof( struct resampled ) ) ;
    for( j = 0 ; j < Input -> Data.Ndata[i] ; j++ ) {
      der[j] = init_dist( NULL ,
			  Input -> Data.x[j+shift].NSAMPLES ,
			  Input -> Data.x[j+shift].restype ) ;
    }

    // loop bootstraps
    for( k = 0 ; k < Input -> Data.x[shift].NSAMPLES ; k++ ) {
      
      double X[ Input -> Data.Ndata[i] ] , Y[ Input -> Data.Ndata[i] ] ;

      size_t idx = 0 ;
      X[0] = Input -> Data.x[shift].resampled[k] ;
      Y[0] = Input -> Data.y[shift].resampled[k] ;
     
      for( j = 1 ; j < Input -> Data.Ndata[i] ; j++ ) {
	if( fabs( Input -> Data.x[j+shift].resampled[k] - Input -> Data.x[j+shift-1].resampled[k] ) < 1E-12 ) continue ;
	X[idx] = Input -> Data.x[j+shift].resampled[k] ;
	Y[idx] = Input -> Data.y[j+shift].resampled[k] ;
	idx++ ;
      }
      double **ders = get_ders( Y , X , idx , 12 ) ;

      // set the derivatives
      for( l = 0 ; l < idx ; l++ ) {
	der[l].resampled[k] = ders[l][1] ;
	//printf( "HERE :: %f \n" , ders[l][0] ) ;
      }
	
      // 
      for( l = 0 ; l < idx ; l++ ) {
	free( ders[l] ) ;
      }
      free( ders ) ;
      //
    }

    // do the average
    double X[ Input -> Data.Ndata[i] ] , Y[ Input -> Data.Ndata[i] ] ;

    size_t idx = 0 ;
    X[0] = Input -> Data.x[shift].avg ;
    Y[0] = Input -> Data.y[shift].avg;
     
    for( j = 1 ; j < Input -> Data.Ndata[i] ; j++ ) {
      if( fabs( Input -> Data.x[j+shift].avg - Input -> Data.x[j+shift-1].avg ) < 1E-12 ) continue ;
      X[idx] = Input -> Data.x[j+shift].avg ;
      Y[idx] = Input -> Data.y[j+shift].avg ;
      idx++ ;
    }
    double **ders = get_ders( Y , X , idx , 4 ) ;

    // set the derivatives
    for( l = 0 ; l < idx ; l++ ) {
      der[l].avg = ders[l][1] ;
      //printf( "HERE :: %f \n" , ders[l][0] ) ;
    }
	
    // 
    for( l = 0 ; l < idx ; l++ ) {
      free( ders[l] ) ;
    }
    free( ders ) ;

    // print out the derivatives
    for( j = 0 ; j < idx ; j++ ) {
      compute_err( &der[j] ) ;
      mult_constant( &der[j] , 12*M_PI*M_PI*X[j]*5./9. ) ;
      printf( "%f %f %f \n" , X[j] , der[j].avg , der[j].err ) ;
    }
    
    shift += Input -> Data.Ndata[i] ;
  }
  
  return SUCCESS ;
}

#endif
