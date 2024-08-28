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
  const size_t t0 = 1 ;
  struct resampled *evalues = solve_GEVP( Input -> Data.y ,
					  LT ,
					  Input -> Fit.N ,
					  Input -> Fit.M ,
					  4 , 6 ) ;

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
    /*
    struct resampled temp = init_dist( &Input -> Data.y[ j+2*LT ] ,
				       Input -> Data.y[ j+2*LT ].NSAMPLES ,
				       Input -> Data.y[ j+2*LT ].restype ) ;
    mult( &Input -> Data.y[j] , Input -> Data.y[ j + LT ] ) ;

    raise( &temp , 2. ) ;
    divide( &Input -> Data.y[ j ] , temp ) ;

    divide( &Input -> Data.y[ j ] , Input->Data.y[j+LT] ) ;
    */
    raise( &Input -> Data.y[j+LT] , 4./2 ) ;
    divide( &Input -> Data.y[ j ] , Input->Data.y[j+LT] ) ;
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

static void
deriv( struct resampled *res ,
       struct resampled *tempb , 
       const struct resampled *y ,
       const size_t i ,
       const size_t LT ,
       const size_t offset )
{
  // compute the derivative and put it in tempf
  if( i == 0 ) {
    // f(x+1)-f(x)
    equate( res , y[1+offset] ) ;
    equate( tempb , y[offset] ) ;
    subtract( res , *tempb ) ; 
  } else if( i == LT-1 ) {
    // f(x) - f(x-1)
    equate( res , y[i+offset] ) ;
    equate( tempb , y[i-1+offset] ) ;
    subtract( res , *tempb ) ; 
  } else {
    equate( res , y[i+1+offset] ) ;
    equate( tempb , y[i-1+offset] ) ;
    subtract( res , *tempb ) ;
    mult_constant( res , 0.5 ) ;
  }
}

// form from https://arxiv.org/pdf/1502.04999.pdf
const double getCA( const double beta )
{
  const double g0sq = 6/beta ;
  return -0.006033*g0sq*( 1 + exp( 9.2056 - 13.9847/g0sq ) ) ;
}

int
PCAC_analysis( struct input_params *Input )
{
  size_t i , j = 0 ;
  const size_t N = Input -> Fit.N ;
  const size_t LT = Input -> Data.Ndata[0] ;

  struct resampled tempf = init_dist( &Input -> Data.y[ j ] ,
				      Input -> Data.y[ j ].NSAMPLES ,
				      Input -> Data.y[ j ].restype ) ;
  struct resampled tempb = init_dist( &Input -> Data.y[ j ] ,
				      Input -> Data.y[ j ].NSAMPLES ,
				      Input -> Data.y[ j ].restype ) ;

  // first we O(a) improve the axial current
  // ASSUMES first element is PP and second set is  <At P>
  const double beta = 3.34 ;
  const double CA = getCA( beta ) ;

  printf( "CA(%f) :: %e\n" , beta , CA ) ;

  for( j = 0 ; j < LT ; j++ ) {
    deriv( &tempf , &tempb , Input -> Data.y , j , LT , 0 ) ;
    mult_constant( &tempf , CA ) ;
    add( &Input -> Data.y[j+LT] , tempf ) ; 
  }
  
  // ASSUMES first element is PP and second set is  <At P>
  for( j = 0 ; j < LT ; j++ ) {
    
    deriv( &tempf , &tempb , Input -> Data.y , j , LT , LT ) ;

    // divide by PP correlator
    divide( &tempf , Input -> Data.y[j] ) ;
    mult_constant( &tempf , -0.5 ) ;
    equate( &Input -> Data.y[j] , tempf ) ;
  }
  free( tempf.resampled ) ;
  free( tempb.resampled ) ;
    
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
