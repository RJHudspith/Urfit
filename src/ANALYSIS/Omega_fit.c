#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams

#include "fvol_delta.h"
#include "write_flat.h"

#include "resampled_ops.h"
#include "correlation.h"

#define NA (6)
#define NPARS (13)

// just a linear fit
int
omega_analysis( struct input_params *Input )
{
  init_phi3( 1000 ) ;
  set_phi3( 0 , true ) ;

  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  // convert to fermi
  struct resampled tmp ;
  for( size_t a = 0 ; a < NA+NPARS ; a++ ) {
    if( a == 0 || a > NPARS ) {
      tmp = init_dist( &fit[a] , fit[a].NSAMPLES , fit[a].restype ) ;
      raise( &tmp , -1 ) ;
      fprintf( stdout , "In GeV -> %e %e || %g (percentage) relerr\n" , tmp.avg , tmp.err , 100*tmp.err/tmp.avg ) ;
    }
  }
  free( tmp.resampled ) ;

  // convert to fm
  const double h = 0.197326980 ;
  for( size_t a = 0 ; a < NA+NPARS ; a++ ) {
    if( a == 0 || a > NPARS ) {
      mult_constant( &fit[a] , h ) ;
      fprintf( stdout , "In fermi -> %e %e || %g \% relerr\n" , fit[a].avg , fit[a].err , 100*fit[a].err/fit[a].avg ) ;
    }
  }

  // compute correlation
  struct resampled *corr = malloc( (NA+1) * sizeof( struct resampled ) ) ;
  double **relation = malloc( (NA+1) * sizeof( double* ) );
  size_t idx = 0 ;
  for( size_t a = 0 ; a < NA+NPARS ; a++ ) {
    if( a == 0 || a > NPARS ) {
      relation[idx] = malloc( 6*sizeof( double ) ) ;
      corr[idx] = init_dist( &fit[a] , fit[a].NSAMPLES , fit[a].restype ) ;
      idx ++ ;
    }    
  }

  // do latst one
  relation[idx] = malloc( 6*sizeof( double ) ) ;
  corr[idx] = init_dist( &fit[13] , fit[13].NSAMPLES , fit[13].restype ) ;
  raise( &corr[idx] , 0.5 ) ;
  mult_constant( &corr[idx] , 0.197326980 ) ;
  idx ++ ; 

  compute_upper_correlation( relation , corr , NA+1 , CORRELATED ) ;
  fill_lower_triangular( relation , NA+1 ) ;
  write_corrmatrix( (const double **)relation , NA+1 , CORRELATED ) ;

  modified_covariance( relation , NA+1 ) ;
  fill_lower_triangular( relation , NA+1 ) ;
  write_corrmatrix( (const double **)relation , NA+1 , CORRELATED ) ;

  /*
  // computes inverse by cholesky decomp
  printf( "Inverted covariance\n" ) ;
  
  gsl_matrix *A = gsl_matrix_calloc( 7 , 7 ) ;
  for( size_t i = 0 ; i < 7 ; i++ ) {
    size_t j ;
    for( j = 0 ; j < 7 ; j++ ) {
      gsl_matrix_set( A , i , j , relation[i][j] ) ;
    }
  }    
  gsl_linalg_cholesky_decomp( A ) ;
  gsl_linalg_cholesky_invert( A ) ;
  for( size_t i = 0 ; i < 7 ; i++ ) {
    size_t j ;
    for( j = 0 ; j < 7 ; j++ ) {
      relation[i][j] = gsl_matrix_get( A , i , j ) ;
    }
  }
  write_corrmatrix( (const double **)relation , NA+1 , CORRELATED ) ;
  */
  
  for( size_t a = 0 ; a < NA+1 ; a++ ) {
    free( corr[a].resampled ) ;
    free( relation[a] ) ;
  }
  free( corr ) ;
  free( relation ) ;


  // conversions for Matthias all in units of powers of GeV
  const double Momega = 1.6695 ; //1.67245 ; // GeV

  struct resampled t08_Momega = init_dist( &fit[13] , fit[13].NSAMPLES , fit[13].restype ) ;
  mult_constant( &t08_Momega , 8*Momega ) ;

  // sqrt( sqrt(8t0)*Momega )
  struct resampled fac1 = init_dist( &fit[13] , fit[13].NSAMPLES , fit[13].restype ) ;
  raise( &fac1 , 0.5 ) ;
  mult_constant( &fac1 , sqrt(8)*Momega ) ;
  raise( &fac1 , 0.5 ) ;

  struct resampled fac2 = init_dist( &fit[13] , fit[13].NSAMPLES , fit[13].restype ) ;
  mult_constant( &fac2 , 8*Momega ) ;
  mult( &fac2 , fit[13] ) ;

  mult( &fit[12] , t08_Momega ) ;
  fprintf( stdout , "d_0 :: %e +/- %e GeV^-1\n" , fit[12].avg , fit[12].err ) ;
  
  // d_\omega
  mult( &fit[1] , t08_Momega ) ;
  fprintf( stdout , "d_\Omega :: %e +/- %e GeV^-1\n" , fit[1].avg , fit[1].err ) ;

  // d_\omega prime
  mult( &fit[2] , t08_Momega ) ;
  fprintf( stdout , "d_\Omega^prime :: %e +/- %e GeV^-1\n" , fit[2].avg , fit[2].err ) ;

  // e_\omega^{\eta}
  mult( &fit[3] , fac2 ) ;
  fprintf( stdout , "e_\Omega :: %e +/- %e GeV^-3\n" , fit[3].avg , fit[3].err ) ;

  // c_\Omega
  mult( &fit[6] , fac1 ) ;
  fprintf( stdout , "c_\Omega :: %e +/- %e \n" , fit[6].avg , fit[6].err ) ;

  // h_\Omega
  mult( &fit[7] , fac1 ) ;
  fprintf( stdout , "h_\Omega :: %e +/- %e \n" , fit[7].avg , fit[7].err ) ;

  // G_\OmegaQ^{(s)}
  mult( &fit[4] , t08_Momega ) ;
  fprintf( stdout , "g_{\Omega K}^{(s)} :: %e +/- %e GeV^-1\n" , fit[4].avg , fit[4].err ) ;

  // G_\OmegaK^{(V)}
  mult( &fit[5] , t08_Momega ) ;
  fprintf( stdout , "g_{\Omega K}^{(v)} :: %e +/- %e GeV^-1\n" , fit[5].avg , fit[5].err ) ;

  // G_\Omega\pi^{(s)}
  mult( &fit[8] , t08_Momega ) ;
  fprintf( stdout , "g_{\Omega pi}^{(s)} :: %e +/- %e GeV^-1\n" , fit[8].avg , fit[8].err ) ;

  // G_\Omega\pi^{(v)}
  mult( &fit[9] , t08_Momega ) ;
  fprintf( stdout , "g_{\Omega pi}^{(v)} :: %e +/- %e GeV^-1\n" , fit[9].avg , fit[9].err ) ;

  // G_\Omega\eta^{(s)}
  mult( &fit[10] , t08_Momega ) ;
  fprintf( stdout , "g_{\Omega eta}^{(s)} :: %e +/- %e GeV^-1\n" , fit[10].avg , fit[10].err ) ;

  // G_\Omega\eta^{(s)}
  mult( &fit[11] , t08_Momega ) ;
  fprintf( stdout , "g_{\Omega eta}^{(v)} :: %e +/- %e GeV^-1\n" , fit[11].avg , fit[11].err ) ;

  // t0 in fm
  raise( &fit[13] , 0.5 ) ;
  mult_constant( &fit[13] , 0.197326980 ) ;
  fprintf( stdout , "rt0 :: %f +/- %f\n" , fit[13].avg , fit[13].err ) ;
  
  return SUCCESS ;
}
