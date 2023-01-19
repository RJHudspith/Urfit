#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "gevp.h"
#include "make_xmgrace.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"
#include "stats.h"
#include "write_flat.h"

static int
write_fitmass_graph( FILE *file , 
		     const struct resampled mass ,
		     const double lo ,
		     const double hi ,
		     const int t0 )
{
  fprintf( file , "%e %e\n%e %e\n\n" , lo+t0 , mass.err_hi , hi+t0 , mass.err_hi ) ; 
  fprintf( file , "%e %e\n%e %e\n\n" , lo+t0 , mass.avg , hi+t0 , mass.avg ) ;
  fprintf( file , "%e %e\n%e %e\n\n" , lo+t0 , mass.err_lo , hi+t0 , mass.err_lo ) ;
  return SUCCESS ;
}

int
binding_corr_analysis( struct input_params *Input )
{
  size_t i , j ;
  const size_t N = Input -> Fit.N ;
  const size_t LT = Input -> Data.Ndata[0] ;

  for( j = 0 ; j < LT ; j++ ) {
    mult( &Input -> Data.y[ j+(Input -> Data.Nsim-2)*LT ] ,
	  Input -> Data.y[ j+(Input -> Data.Nsim-1)*LT ] ) ;
  }
  
  for( i = 0 ; i < Input -> Data.Nsim-2 ; i++ ) {
    const size_t shift = Input -> Data.Nsim-2 ;
    for( j = 0 ; j < LT ; j++ ) {  
      divide( &Input -> Data.y[ j+i*LT ] ,
	      Input -> Data.y[ j+shift*LT ] ) ;
    }
  }

  if( ( Input -> Fit.N*Input-> Fit.M ) != (Input -> Data.Nsim-2) ) {
    fprintf( stderr , "GEVP N,M not equal to Nsim-2\n" ) ;
    goto end ;
  }
  
  // compute evalues
  const size_t t0 = 2 ;
  struct resampled *evalues = solve_GEVP( Input -> Data.y ,
					  LT ,
					  Input -> Fit.N ,
					  Input -> Fit.M ,
					  t0 , t0+1 ) ;

  for( i = 0 ; i < Input -> Fit.N ; i++ ) {
    for( j = 0 ; j < Input -> Data.Ndata[0] ; j++ ) {
      equate( &Input -> Data.y[j+i*Input->Data.Ndata[0]] ,
	      evalues[ j + i*Input->Data.Ndata[0] ] ) ;
    }
  }

  // free the eigenvalues
  for( i = 0 ; i < Input ->Data.Ndata[0] * Input -> Fit.N  ; i++ ) {
    free( evalues[i].resampled ) ;
  }
  free( evalues ) ;

  // subtract the shift automatically for the fit?
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    subtract_constant( &Input -> Data.x[ i ] , (double)t0 ) ;
  }

  // change the data size
  const size_t Nsim_prev = Input -> Data.Nsim ;
  const size_t Ntot_prev = Input -> Data.Ntot ;
  Input -> Data.Nsim = N ;
  Input -> Data.Ntot = N*Input->Data.Ndata[0] ;
    
  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;

  // set this stuff
  Input -> Data.Nsim = 1 ;
  Input -> Data.Ntot = Input->Data.Ndata[0] ;
  Input -> Fit.N = Input -> Fit.M = 1 ;
  Input -> Fit.Nlogic = 2 ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  // write out the fit result for reasons
  if( Input -> Fit.Fitdef == EXP ||
      Input -> Fit.Fitdef == COSH ) {
    FILE *massfile = fopen( "massfits.dat" , "w" ) ;
    size_t j , shift = 0 ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      for( j = 0 ; j < 2*Input -> Fit.N ; j+= 2 ) {
	write_fitmass_graph( massfile , Fit[j+1] ,
			     Input -> Traj[i].Fit_Low ,
			     Input -> Traj[i].Fit_High ,
			     t0 ) ;
      }
      shift += Input -> Data.Ndata[i] ;
    }
    fclose( massfile ) ;
  }
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;

  Input -> Data.Nsim = Nsim_prev ;
  Input -> Data.Ntot = Ntot_prev ;

 end:
  
  return SUCCESS ;
}

int
binding_corr_analysis2( struct input_params *Input )
{
  size_t i , j ;
  const size_t N = Input -> Fit.N ;
  const size_t LT = Input -> Data.Ndata[0] ;
  
  for( j = 0 ; j < LT ; j++ ) {
    struct resampled temp = init_dist( &Input -> Data.y[ j+LT ] ,
				       Input -> Data.y[ j+LT ].NSAMPLES ,
				       Input -> Data.y[ j+LT ].restype ) ;
    raise( &temp , 2./2 ) ;
    divide( &Input -> Data.y[ j ] , temp ) ;
  }

  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

    // write out a flat file
  if( Input -> Fit.Fitdef == EXP ||
      Input -> Fit.Fitdef == COSH ||
      Input -> Fit.Fitdef == SINH ) {
    
    struct resampled mpi2 = init_dist( NULL ,
				       Fit[1].NSAMPLES ,
				       Fit[1].restype ) ;
    write_flat_dist( &Fit[1] , &mpi2 , 1 , "Mass_0.flat" ) ;
    FILE *massfile = fopen( "massfits.dat" , "w+a" ) ;

    write_fitmass_graph( massfile , Fit[1] ,
			 Input -> Traj[0].Fit_Low ,
			 Input -> Traj[0].Fit_High , 0 ) ;
    
    fclose( massfile ) ;
  }
  
  return SUCCESS ;
}
