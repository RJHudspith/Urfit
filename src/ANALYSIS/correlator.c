/**
   @file correlator.c
   @brief correlator analysis
 */
#include "gens.h"

#include "blackbox.h"
#include "decays.h"
#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "make_xmgrace.h"
#include "pade_laplace.h"
#include "resampled_ops.h"
#include "stats.h"
#include "write_flat.h"

//#define MATRIX_PRONY
//#define FIT_EFFMASS
//#define PADE_LAPLACE

#define CDIV

void 
my_little_prony( const struct input_params *Input )
{
  make_xmgrace_graph( "MProny.agr" , "t/a" , "am\\seff" ) ;
  
  const size_t Nstates = Input -> Fit.N ;
  size_t i , j , k , l , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

    const size_t Ndata = Input -> Data.Ndata[i] ;
    double data[ Ndata ] ;
    double masses[ Nstates ][ Ndata ] ;

    struct resampled *effmass =
      malloc( Nstates * Ndata * sizeof( struct resampled ) ) ;
    for( j = 0 ; j < Nstates * Ndata ; j++ ) {
      effmass[j] = init_dist( NULL , Input -> Data.x[shift].NSAMPLES ,
			      Input -> Data.x[shift].restype ) ;
    }
    
    for( k = 0 ; k < Input -> Data.x[shift].NSAMPLES ; k++ ) {
      for( j = shift ; j < shift + Ndata ; j++ ) {
	data[ j - shift ] = Input -> Data.y[j].resampled[k] ;
      }
      
      blackbox( data , Ndata , Nstates , masses ) ;
      
      for( l = 0 ; l < Nstates ; l++ ) {
	for( j = 0 ; j < Ndata ; j++ ) {
	  effmass[j + l*Ndata].resampled[k] = masses[l][j] ;
	}
      }
      printf( "[PRONY] Samples %1.2f percent done\n" ,
	      100. * k/(double)Input->Data.x[shift].NSAMPLES ) ;
    }

    // same for the average
    for( j = shift ; j < shift + Ndata ; j++ ) {
      data[ j - shift ] = Input -> Data.y[j].avg ;
    }
    blackbox( data , Ndata , Nstates , masses ) ;
    for( l = 0 ; l < Nstates ; l++ ) {
      for( j = 0 ; j < Ndata ; j++ ) {
	effmass[j + l*Ndata].avg = masses[l][j] ;
      }
    } 
    
    for( j = 0 ; j < Nstates * Ndata ; j++ ) {
      compute_err( &effmass[j] ) ;
    }

    for( l = 0 ; l < Nstates ; l++ ) {
      plot_data( Input -> Data.x + shift , effmass + l*Ndata ,
		 Input -> Data.Ndata[i] ) ;
    }

    for( j = 0 ; j < Nstates * Ndata ; j++ ) {
      free( effmass[j].resampled ) ;
    }
    free( effmass ) ;
    
    shift += Ndata ;
  }
  close_xmgrace_graph() ;

  return ;
}

static int
write_fitmass_graph( FILE *file , 
		     const struct resampled mass ,
		     const double lo ,
		     const double hi )
{
  fprintf( file , "%e %e\n%e %e\n\n" , lo , mass.err_hi , hi , mass.err_hi ) ; 
  fprintf( file , "%e %e\n%e %e\n\n" , lo , mass.avg , hi , mass.avg ) ;
  fprintf( file , "%e %e\n%e %e\n\n" , lo , mass.err_lo , hi , mass.err_lo ) ;
  return SUCCESS ;
}

void 
boot_pade_laplace( const struct input_params *Input )
{
  // do the average first
  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    double p0 ;
    const size_t Ndata = Input -> Data.Ndata[i] ;
    double y[ Ndata ] , x[ Ndata ] ;
    double fparams[ 2*Input -> Fit.N ] ;
    struct resampled *poles = malloc( Input -> Fit.N * 2 * sizeof( struct resampled ) ) ;

    // loop taylor expansion point
    for( p0 = 0.0 ; p0 < 5.0 ; p0 += 0.25 ) {

      for( j = 0 ; j < Input -> Fit.N * 2 ; j++ ) {
	poles[ j ] = init_dist( NULL ,
				Input -> Data.x[shift].NSAMPLES ,
				Input -> Data.x[shift].restype ) ;
      }
      
      size_t idx = 0 ;
      // loop N samples
      for( k = 0 ; k < Input -> Data.x[shift].NSAMPLES ; k++ ) {
	idx = 0 ;
	for( j = shift ; j < shift + Ndata ; j++ ) {
	  x[ idx ] = Input -> Data.x[ j ].resampled[k] ;
	  y[ idx ] = Input -> Data.y[ j ].resampled[k] ;
	  idx++ ;
	}
	pade_laplace( fparams , x , y , Ndata , Input -> Fit.N , p0 ) ;
	// equate fparams to our resampled data
	for( j = 0 ; j < 2*Input -> Fit.N ; j++ ) {
	  poles[j].resampled[k] = fparams[j] ;
	}
      }
      // do the average too
      idx = 0 ;
      for( j = shift ; j < shift + Ndata ; j++ ) {
	x[ idx ] = Input -> Data.x[ j ].avg ;
	y[ idx ] = Input -> Data.y[ j ].avg ;
	idx++ ;
      }
      pade_laplace( fparams , x , y , Ndata , Input -> Fit.N , p0 ) ;
      for( j = 0 ; j < 2*Input -> Fit.N ; j++ ) {
	poles[j].avg = fparams[j] ;
	compute_err( &poles[j] ) ;

	printf( "[PLAP] %f PARAM_%zu %e %e \n" , p0 , j , poles[j].avg , poles[j].err ) ;
      }
      printf( "\n" ) ;
    }

    shift += Ndata ;
    for( j = 0 ; j < 2*Input -> Fit.N ; j++ ) {
      free( poles[j].resampled ) ;
    }
    free( poles ) ;
  }
  return ;
}

int
correlator_analysis( struct input_params *Input )
{
  size_t i , shift = 0 ;

#ifdef MATRIX_PRONY
  fprintf( stdout , "[Correlator] Matrix prony\n" ) ;
  my_little_prony( Input ) ;
#endif

#ifdef PADE_LAPLACE
  boot_pade_laplace( Input ) ;
#endif
  
  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ATANH_EFFMASS ) ;

#ifdef FIT_EFFMASS
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    equate( &Input -> Data.y[i] , effmass[i] ) ;
  }
#endif

  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    free( effmass[i].resampled ) ;
  }
  free( effmass ) ;


  #ifdef CDIV
  for( i = 0 ; i < Input -> Data.Ntot ; i++ ) {
    divide_constant( &Input -> Data.y[i] , 1E14 ) ;
  }
  #endif


  shift = 0 ;
  
  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

  printf( "P2 %f \n" ,
	  Input -> Traj[0].mom[0]*Input -> Traj[0].mom[0] +
	  Input -> Traj[0].mom[1]*Input -> Traj[0].mom[1] +
	  Input -> Traj[0].mom[2]*Input -> Traj[0].mom[2] ) ;
  
  // compute a decay constant
  if( Input -> Fit.Fitdef == PP_AA_WW ||
      Input -> Fit.Fitdef == PP_AA ||
      Input -> Fit.Fitdef == PPAA ) {

    FILE *massfile = fopen( "massfits.dat" , "w" ) ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      write_fitmass_graph( massfile , Fit[0] ,
			   Input -> Traj[i].Fit_Low ,
			   Input -> Traj[i].Fit_High ) ;
    }
    fclose( massfile ) ;


    struct resampled Mass = init_dist( &Fit[0] ,
				       Fit[0].NSAMPLES ,
				       Fit[0].restype ) ;


    write_flat_dist( &Fit[0] , &Fit[0] , 1 , "Mass.flat" ) ;
    
    struct resampled dec = decay( Fit , *Input , 0 , 2 ) ;

    write_flat_dist( &dec , &dec , 1 , "Decay.flat" ) ;

    divide( &Fit[0] , dec ) ;

    raise( &Fit[0] , 2 ) ;
    
    printf( "(M/F)^2 %e %e \n" , Fit[0].avg , Fit[0].err ) ;

    write_flat_dist( &Fit[0] , &Fit[0] , 1 , "MovFsq.flat" ) ;

    //////////////// PCAC ? /////////////////////
    // is d_t A_t^P P^W / 2P^L P^W
    // which I make
    //
    // m_\pi/2 | A^L/P^L |
    equate( &dec , Fit[2] ) ;

    printf( "PCAC %e %e \n" , dec.avg , dec.err ) ;
    divide( &dec , Fit[1] ) ;
    printf( "PCAC %e %e \n" , dec.avg , dec.err ) ;
    mult( &dec , Mass ) ;
    printf( "PCAC %e %e \n" , dec.avg , dec.err ) ;
    mult_constant( &dec , 0.5 ) ;

    printf( "PCAC %e %e \n" , dec.avg , dec.err ) ;
    
    free( dec.resampled ) ;
  }

  // compute a decay constant
  if( Input -> Fit.Fitdef == COSH ) {
    struct resampled dec = decay( Fit , *Input , 1 , 0 ) ;
    free( dec.resampled ) ;
  }
    
  // write out a flat file
  if( Input -> Fit.Fitdef == EXP ||
      Input -> Fit.Fitdef == COSH ||
      Input -> Fit.Fitdef == SINH ) {
    //write_flat_single( &Fit[1] , "Mass.flat" ) ;
    struct resampled mpi2 = init_dist( NULL ,
				       Fit[1].NSAMPLES ,
				       Fit[1].restype ) ;

    size_t shift = 0 , j ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      for( j = 0 ; j < 2*Input -> Fit.N ; j+= 2 ) {
	double psq = 0 ;
	size_t mu ;
	for( mu = 0 ; mu < 3 ; mu++ ) {
	  const double ptilde = sin( Input -> Traj[0].mom[mu]*2.*M_PI/
				     Input -> Traj[0].Dimensions[mu] ) ;
	  psq += ptilde*ptilde ;
	}
	equate_constant( &mpi2 , psq ,
			 Fit[1].NSAMPLES , Fit[1].restype ) ;
	char str[256] ;
	sprintf( str , "Mass_%zu.flat" , j+i*2*Input->Fit.N ) ;
	write_flat_dist( &Fit[j+1+i*2*Input -> Fit.N] , &mpi2 , 1 , str ) ;
      }
    }
    free( mpi2.resampled ) ;

    FILE *massfile = fopen( "massfits.dat" , "w+a" ) ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      for( j = 0 ; j < 2*Input -> Fit.N ; j+= 2 ) {
	write_fitmass_graph( massfile , Fit[j+1+i*2*Input -> Fit.N] ,
			     Input -> Traj[i].Fit_Low ,
			     Input -> Traj[i].Fit_High ) ;
      }
      shift += Input -> Data.Ndata[i] ;
    }
    fclose( massfile ) ;

    struct resampled Omega = init_dist( NULL ,
					Fit[0].NSAMPLES ,
					Fit[0].restype ) ;

    equate_constant( &Omega , 1.67245 , Fit[0].NSAMPLES , Fit[0].restype ) ;
    divide( &Omega , Fit[1] ) ;
    fprintf( stdout , "ainverse %f %f\n" , Omega.avg , Omega.err ) ; 
    
    free( Omega.resampled ) ;
  }

  /*
  subtract( &Fit[3] , Fit[1] ) ;
  printf( "%f +/- %f \n" , Fit[3].avg , Fit[3].err ) ;
  */

    // write out a flat file
  if( Input -> Fit.Fitdef == PEXP ) {
    struct resampled mpi2 = init_dist( NULL ,
				       Fit[1].NSAMPLES ,
				       Fit[1].restype ) ;
    FILE *massfile = fopen( "massfits.dat" , "w+a" ) ;
    size_t i , j ;
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      for( j = 0 ; j < 2 ; j++ ) {
	write_fitmass_graph( massfile , Fit[j*2] ,
			     Input -> Traj[i].Fit_Low ,
			     Input -> Traj[i].Fit_High ) ;
	char str[256] ;
	sprintf( str , "Mass_%zu.flat" , j*2 ) ;
	write_flat_dist( &Fit[j*2] , &mpi2 , 1 , str ) ;
      }
      shift += Input -> Data.Ndata[i] ;
    }
    fclose( massfile ) ;
    free( mpi2.resampled ) ;
  }
  
  free_fitparams( Fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
