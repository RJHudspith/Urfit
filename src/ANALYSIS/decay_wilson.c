#include "gens.h"

#include "effmass.h"
#include "fit_and_plot.h"
#include "init.h"
#include "decays.h"
#include "make_xmgrace.h"
#include "plot_fitfunc.h"
#include "resampled_ops.h"
#include "stats.h"
#include "write_flat.h"

#include "fake.h"

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
const double getCA( const int beta )
{
  const double g0sq = 600./beta ;
  return -0.006033*g0sq*( 1 + exp( 9.2056 - 13.9847/g0sq ) ) ;
}

static struct resampled
getZAdist( const int beta , const double kappal , const double kappas , struct resampled *fit )
{
  struct resampled kappa_crit , ZA ;
  // values for kappa_crit from terrible regensburg paper
  // values for ZA are ZAlsub in Tab7 of https://arxiv.org/pdf/1808.09236
  // ba and babar from other terrible paper https://arxiv.org/pdf/2106.05398
  switch( beta ) {
  case 334 :
    kappa_crit = generate_fake_single( 0.1366938 , 0.000045 , fit[0].NSAMPLES ) ;
    break ;
  case 340 :
    kappa_crit = generate_fake_single( 0.1369153 , 0.000009 , fit[0].NSAMPLES ) ;
    ZA = generate_fake_single( 0.75642 , 0.00072 , fit[0].NSAMPLES ) ;
    break ;
  case 346 :
    kappa_crit = generate_fake_single( 0.1370613 , 0.000010 , fit[0].NSAMPLES ) ;
    ZA = generate_fake_single( 0.76169 , 0.00093 , fit[0].NSAMPLES ) ;
    break ;
  case 355 :
    kappa_crit = generate_fake_single( 0.1371715 , 0.000010 , fit[0].NSAMPLES ) ;
    ZA = generate_fake_single( 0.76979 , 0.00043 , fit[0].NSAMPLES ) ;
    break ;
  case 370 :
    kappa_crit = generate_fake_single( 0.1371530 , 0.000009 , fit[0].NSAMPLES ) ;
    ZA = generate_fake_single( 0.78378 , 0.00047 , fit[0].NSAMPLES ) ;
    break ;
  case 380 :
    kappa_crit = generate_fake_single( 0.1369767 , 0.000026 , fit[0].NSAMPLES ) ;
    ZA = generate_fake_single( 0.79667 , 0.00047 , fit[0].NSAMPLES ) ;
    break ;
  default :
    exit(1) ;
    break ;
  }
  printf( "kappa_crit %e %e\n" , kappa_crit.avg , kappa_crit.err ) ;
  printf( "ZA %e %e\n" , ZA.avg , ZA.err ) ;
  raise( &kappa_crit , -1 ) ; mult_constant( &kappa_crit , 0.5 ) ;
  printf( "amc %e %e\n" , kappa_crit.avg , kappa_crit.err ) ;
  struct resampled aml = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;
  equate_constant( &aml , kappal , fit[0].NSAMPLES , fit[0].restype ) ;
  raise( &aml , -1 ) ; mult_constant( &aml , 0.5 ) ;
  subtract( &aml , kappa_crit ) ;
  
  struct resampled ams = init_dist( NULL , fit[0].NSAMPLES , fit[0].restype ) ;
  equate_constant( &ams , kappas , fit[0].NSAMPLES , fit[0].restype ) ; 
  raise( &ams , -1 ) ; mult_constant( &ams , 0.5 ) ;
  subtract( &ams , kappa_crit ) ;
  
  printf( "aml %e %e\n" , aml.avg , aml.err ) ;
  printf( "ams %e %e\n" , ams.avg , ams.err ) ;

  const double g0sq = 600./beta ;
  //const double ba = 1 + g0sq*( 0.0881*(4/3.) + 0.0113*g0sq ) ;

  struct resampled ba = generate_fake_single( 0.0113 , 0.0044 , fit[0].NSAMPLES ) ;
  mult_constant( &ba , g0sq*g0sq ) ;
  add_constant( &ba , 1 + 0.0881*4/3*g0sq ) ;

  add( &aml , ams ) ;
  mult_constant( &aml , 0.5 ) ;
  mult( &ba , aml ) ;
  printf( "ba part %e %e \n" , ba.avg , ba.err ) ;

  add_constant( &ba , 1.0 ) ;
  mult( &ZA , ba ) ;

  printf( "ZA -> %e %e \n" , ZA.avg , ZA.err ) ;
  return ZA ;
}

const double
getZA( const int beta , const double kappal , const double kappas )
{
  double kappa_crit = 0.125 , ZA = 1 , babar = 0 ;
  // values for kappa_crit from terrible regensburg paper
  // values for ZA are ZAlsub in Tab7 of https://arxiv.org/pdf/1808.09236
  // ba and babar from other terrible paper https://arxiv.org/pdf/2106.05398
  switch( beta ) {
  case 334 :
    kappa_crit = 0.1366938 ;
    break ;
  case 340 :
    kappa_crit = 0.1369153 ;
    ZA = 0.75642 ;
    //babar = -0.11 ;
    break ;
  case 346 :
    kappa_crit = 0.1370613 ;
    ZA = 0.76169 ;
    //babar = 0.10 ;
    break ;
  case 355 :
    kappa_crit = 0.1371715 ;
    ZA = 0.76979 ;
    //babar = -0.04 ;
    break ;
  case 370 :
    kappa_crit = 0.1371530 ;
    ZA = 0.78378 ;
    //babar = -0.05 ;
    break ;
  case 380 :
    kappa_crit = 0.1369767 ;
    break ;
  default :
    exit(1) ;
    break ;
  }
  const double amql = 1/(2*kappal) - 1/(2*kappa_crit), amqs = 1/(2*kappas) - 1/(2*kappa_crit) ;
  printf( "What %e %e\n" , amql , amqs ) ;
  const double amqav = 0.0 ; //(2*amql + amqs)/3. ;
  const double amqrr = (amql+amqs)/2. ; // only care about light-quark f_\pi for now

  
  const double g0sq = 600./beta ;
  const double ba = 1 + g0sq*( 0.0881*(4/3.) + 0.0113*g0sq ) ;
  return ZA*(1 + 3*babar*amqav + ba*amqrr ) ;
}

int
fpi_analysis( struct input_params *Input )
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
  const int beta = 355 ;
  const double kappal = 0.1372067 , kappas = 0.136546844 ;
  const double CA = getCA( beta ) ;

  printf( "CA(%f) :: %e\n" , beta/100. , CA ) ;
  printf( "kl :: %e || ks :: %e\n" , kappal , kappas ) ;

  for( j = 0 ; j < LT ; j++ ) {
    deriv( &tempf , &tempb , Input -> Data.y , j , LT , 0 ) ;
    mult_constant( &tempf , CA ) ;
    add( &Input -> Data.y[j+LT] , tempf ) ; 
  }
  free( tempf.resampled ) ;
  free( tempb.resampled ) ;
 
    
  // compute an effective mass 
  struct resampled *effmass = effective_mass( Input , ACOSH_ITERATIVE_EFFMASS ) ;

  // perform a fit
  double Chi ;
  struct resampled *Fit = fit_and_plot( *Input , &Chi ) ;

    // write out a flat file
  if( Input -> Fit.Fitdef == PPAA) {
    struct resampled mpi2 = init_dist( NULL ,
				       Fit[1].NSAMPLES ,
				       Fit[1].restype ) ;
    write_flat_dist( &Fit[0] , &mpi2 , 1 , "Mass_0.flat" ) ;
    FILE *massfile = fopen( "massfits.dat" , "w+a" ) ;
    write_fitmass_graph( massfile , Fit[0] ,
			 Input -> Traj[0].Fit_Low ,
			 Input -> Traj[0].Fit_High , 0 ) ;
    fclose( massfile ) ;
  }

  
  // improved current is in Data.y[j+LT]
  struct resampled fpi = decay( Fit , *Input , 0 , 2 ) ;

  fprintf( stdout , "ZA %f\n" , getZA( beta , kappal , kappas ) ) ; 

  mult_constant( &fpi , getZA( beta , kappal , kappas ) ) ;

  fprintf( stdout , "Renormalised f_pi %e +/- %e\n" , fpi.avg , fpi.err ) ;
  free( fpi.resampled ) ;

  fpi = decay( Fit , *Input , 0 , 2 ) ;
  struct resampled ZA = getZAdist( beta , kappal , kappal , Fit ) ;
  mult( &fpi , ZA ) ;
  fprintf( stdout , "Renormalised f_pi %e +/- %e\n" , fpi.avg , fpi.err ) ;

  fprintf( stdout , "Stat error rat %e\n" , 100*fpi.err/fpi.avg ) ;
  
  return SUCCESS ;
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
  const int beta = 355 ;
  const double CA = getCA( beta ) ;

  printf( "CA(%f) :: %e\n" , beta/100. , CA ) ;

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

