#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h" // free_fitparams

#include "fvol_delta.h"
#include "write_flat.h"

#include "resampled_ops.h"
#include "correlation.h"

#define NA (6)
#define NPARS (12)

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
  struct resampled *corr = malloc( NA * sizeof( struct resampled ) ) ;
  double **relation = malloc( NA * sizeof( double* ) );
  size_t idx = 0 ;
  for( size_t a = 0 ; a < NA+NPARS ; a++ ) {
    if( a == 0 || a > NPARS ) {
      relation[idx] = malloc( 6*sizeof( double ) ) ;
      corr[idx] = init_dist( &fit[a] , fit[a].NSAMPLES , fit[a].restype ) ;
      idx ++ ;
    }
  }

  compute_upper_correlation( relation , corr , NA , CORRELATED ) ;
  fill_lower_triangular( relation , NA ) ;
  write_corrmatrix( relation , NA , CORRELATED ) ;
  
  modified_covariance( relation , NA ) ;
  fill_lower_triangular( relation , NA ) ;
  write_corrmatrix( relation , NA , CORRELATED ) ;

  for( size_t a = 0 ; a < NA ; a++ ) {
    free( corr[a].resampled ) ;
    free( relation[a] ) ;
  }
  free( corr ) ;
  free( relation ) ;


  // conversions for Matthias all in units of powers of GeV
  const double t08 = 8*0.5391144626032651 ;
  const double rtt08 = sqrt( t08 ) ; // GeV^-1
  const double Momega = 1.67245 ; // GeV
  const double f = 0.0924 ; // GeV

  // d_\omega
  mult_constant( &fit[1] , t08*Momega ) ;
  fprintf( stdout , "d_\Omega :: %e +/- %e GeV^-1\n" , fit[1].avg , fit[1].err ) ;

  // d_\omega prime
  mult_constant( &fit[2] , t08*Momega ) ;
  fprintf( stdout , "d_\Omega^prime :: %e +/- %e GeV^-1\n" , fit[2].avg , fit[2].err ) ;

  // e_\omega^{\eta}
  mult_constant( &fit[3] , t08*t08*Momega ) ;
  fprintf( stdout , "e_\Omega :: %e +/- %e GeV^-3\n" , fit[3].avg , fit[3].err ) ;

  // c_\Omega
  mult_constant( &fit[6] , sqrt(rtt08*Momega) ) ;
  fprintf( stdout , "c_\Omega :: %e +/- %e \n" , fit[6].avg , fit[6].err ) ;

  // h_\Omega
  mult_constant( &fit[7] , sqrt(rtt08*Momega) ) ;
  fprintf( stdout , "h_\Omega :: %e +/- %e \n" , fit[7].avg , fit[7].err ) ;

  // G_\OmegaQ^{(s)}
  mult_constant( &fit[4] , t08*Momega ) ;
  fprintf( stdout , "g_{\Omega K}^{(s)} :: %e +/- %e GeV^-1\n" , fit[4].avg , fit[4].err ) ;

  // G_\OmegaK^{(V)}
  mult_constant( &fit[5] , t08*Momega ) ;
  fprintf( stdout , "g_{\Omega K}^{(v)} :: %e +/- %e GeV^-1\n" , fit[5].avg , fit[5].err ) ;

  // G_\Omega\pi^{(s)}
  mult_constant( &fit[8] , t08*Momega ) ;
  fprintf( stdout , "g_{\Omega pi}^{(s)} :: %e +/- %e GeV^-1\n" , fit[8].avg , fit[8].err ) ;

  // G_\Omega\pi^{(v)}
  mult_constant( &fit[9] , t08*Momega ) ;
  fprintf( stdout , "g_{\Omega pi}^{(v)} :: %e +/- %e GeV^-1\n" , fit[9].avg , fit[9].err ) ;

  // G_\Omega\eta^{(s)}
  mult_constant( &fit[10] , t08*Momega ) ;
  fprintf( stdout , "g_{\Omega eta}^{(s)} :: %e +/- %e GeV^-1\n" , fit[10].avg , fit[10].err ) ;

  // G_\Omega\eta^{(s)}
  mult_constant( &fit[11] , t08*Momega ) ;
  fprintf( stdout , "g_{\Omega eta}^{(v)} :: %e +/- %e GeV^-1\n" , fit[11].avg , fit[11].err ) ;
  
  return SUCCESS ;
}
