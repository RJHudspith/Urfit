/**
   @file Nint.c
   @brief trapezoid and generalized Simpson's numerical integrators
 */
#include "gens.h"

#include "fit_chooser.h"
#include "ffunction.h"
#include "resampled_ops.h"
#include "stats.h"

// simpson evaluation
static double
simp_eval( const double *fparams ,
	   const struct data_info Data ,
	   const struct fit_info Fit ,
	   struct fit_descriptor fdesc ,
	   const double x ,
	   const double h ,
	   const size_t shift )
{
  struct x_desc xdesc = { x , Data.LT[shift] , Fit.N , Fit.M } ;

  double sum = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

  xdesc.X = x+h/2. ;
  
  sum += 4*fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

  xdesc.X = x+h ; 
  sum += fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

  return h*sum ;
}

// driver for the adaptive simpson's
// usual two half-step approach to estimate the error
static double
adaptive_simpsons( const double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit ,
		   struct fit_descriptor fdesc ,
 		   const double low , 
		   const double upp ,
		   const double eps ,
		   const size_t shift )
{
  const double inveps = 1.0 / eps ; 
  register double fx , fxph_4 , fxph_2 , fxp3h_4 , fxph ; // some evaluations I use a lot
  register double sum = 0.0 , diff ;
  double h = ( upp - low ) / 100. ; // initial guess
  double x = low ;
  const size_t max_steps = 150 ;
  while( x < upp ) {    
    register double step = 0.0 , step2 = 0.0 ;
    struct x_desc xdesc = { x , Data.LT[shift] , Fit.N , Fit.M } ;
    
    int nsteps = 0 ;
    while( nsteps < max_steps ) {

      xdesc.X = x ;
      fx      = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

      xdesc.X = x+h/4. ;
      fxph_4  = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

      xdesc.X = x+h/2. ;
      fxph_2  = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

      xdesc.X = x+3*h/4. ;
      fxp3h_4 = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

      xdesc.X = x+h ;
      fxph    = fdesc.func( xdesc , fparams , Fit.map[shift].bnd ) ;

      step  = h * ( fx + 4. * fxph_2 + fxph ) ;
      step2 = h * ( fx + 4. * ( fxph_4 + fxp3h_4 ) + 2. * fxph_2 + fxph ) / 2.0 ;
 
      diff = fabs( step - step2 ) * inveps ;
      
      // break if the diff is within tolerance
      if( diff < 1.0 ) 
	break ;

      // if we run into too much trouble we return a NAN
      if( nsteps == max_steps-1 ) {
	fprintf( stderr , "[ADS] ill convergence met\n" ) ;
	fprintf( stderr , "[ADS] %e \n" , diff ) ;
	return ( sum + simp_eval( fparams , Data , Fit ,
				  fdesc , x , upp-x , shift ) ) / 6. ;
      }

      h *= 0.9 * pow( diff , -0.25 ) ; 
      nsteps++ ;
    }
    sum = sum + ( step2 + ( step2 - step ) / 15. ) ;
    x += h ;
    if( diff > 1E-15 ) {
      h *= 0.9 * pow( diff , -0.2 ) ;
    }
  }
  // small correction exactly to "x"
  return ( sum + simp_eval( fparams , Data , Fit ,
			    fdesc , x , upp-x , shift ) ) / 6. ;
}

static double
general_simpsons_array_O3( const double *x , 
			   const double *y , 
			   const size_t i )
{
  const register double a = x[ i ] ;
  const register double b = x[ i + 1 ] ;
  const register double c = x[ i + 2 ] ;
  const register double fa = y[ i ] / 6.0 ;
  const register double fb = y[ i + 1 ] / 6.0 ;
  const register double fc = y[ i + 2 ] / 6.0 ;

  // up to some tolerance we just return the trapezoid evaluation
  if( ( b - a ) < 1E-20 ) {
    return ( fc - fb ) * ( c - b ) * 0.5 ;
  } else if( ( c - b ) < 1E-20 ) {
    return ( fb - fa ) * ( b - a ) * 0.5 ;
  } else {
    register double loc_sum = a * a * ( fb - fc ) ;
    loc_sum += c * c * ( fb - fa ) ;
    loc_sum += -3 * b * b * ( fa + fc ) ;
    loc_sum += 2 * b * c * ( 2 * fa + fc ) ;
    loc_sum += 2 * a * b * ( fa + 2 * fc ) ;
    loc_sum += -2 * a * c * ( fa + fb + fc ) ;
    loc_sum *= ( c - a ) / ( ( a - b ) * ( b - c ) ) ;
    return loc_sum ;
  }
}

static double
general_simpsons_array_O4( const double *x , 
			   const double *y , 
			   const size_t i )
{
  const double a = x[ i ] ;
  const double b = x[ i + 1 ] ;
  const double c = x[ i + 2 ] ;
  const double d = x[ i + 3 ] ;

  const double fa = y[ i ] ;
  const double fb = y[ i + 1 ] ;
  const double fc = y[ i + 2 ] ;
  const double fd = y[ i + 3 ] ;

  long double loc_sumfa , loc_sumfb , loc_sumfc , loc_sumfd ;

  // f(a) multiplier
  loc_sumfa  = 3 * a * a + 6 * b * c - 2 * ( b + c ) * d ;
  loc_sumfa += d * d + 2 * a * ( -2 * ( b + c ) + d ) ;
  loc_sumfa *= -fa / ( ( a - b ) * ( a - c ) ) ;

  // f(b) multiplier
  loc_sumfb  = ( a - d ) * ( a - d ) * ( a - 2. * c + d ) ;
  loc_sumfb *= fb / ( ( b - a ) * ( b - c ) * ( b - d ) ) ;
  
  // f(c) multiplier
  loc_sumfc  = ( a - d ) * ( a - d ) * ( a - 2. * b + d ) ;
  loc_sumfc *= fc / ( ( c - a ) * ( c - b ) * ( c - d ) ) ;

  // f(d) multiplier
  loc_sumfd  = a * a + 6 * b * c - 2 * a * ( b + c - d ) ;
  loc_sumfd += -4. * b * d - 4. * c * d + 3 * d * d ;
  loc_sumfd *= fd / ( ( b - d ) * ( d - c ) ) ;

  return ( loc_sumfa + loc_sumfb + loc_sumfc + loc_sumfd ) * ( a - d ) / 12. ;
}

// numerical integrator
static double
simpsons_arr5( const double *y , 
	       const double *x ,
	       const size_t N )
{
  register double sum = 0.0 ;

  const double lo = general_simpsons_array_O4( x , y , 0 )
    - general_simpsons_array_O3( x , y , 1 ) ;

  size_t i ;
  for( i = 0 ; i < N - 2 ; i++ ) {
    sum += general_simpsons_array_O3( x , y , i ) ;
  }

  const double hi = general_simpsons_array_O4( x , y , N-4 )
    - general_simpsons_array_O3( x , y , N - 4 ) ;

  return ( lo + sum + hi )/2. ;
}

// simple trapezoid rule
static double
run_the_trap( const double *y ,
	      const double *x ,
	      const size_t N )
{
  register double sum = 0.0 ;
  size_t i ;
  for( i = 0 ; i < N-1 ; i++ ) {
    sum += ( y[i] + y[i+1] )*( x[i+1] - x[i] )/2. ;
  }
  return sum ;
}

// generic numerical integrator
struct resampled
Nint( const struct resampled *dataX ,
      const struct resampled *dataY ,
      const size_t Ndata ,
      const bool is_trap )
{
  static double (*integrator) ( const double *y ,
				const double *x ,
				const size_t N ) ;
  if( is_trap ) {
    integrator = run_the_trap ;
  } else {
    integrator = simpsons_arr5 ;
  }
  
  struct resampled Int = init_dist( NULL ,
				    dataX[0].NSAMPLES ,
				    dataX[0].restype ) ;
  
  // loop number of samples
  size_t k ;
  for( k = 0 ; k < Int.NSAMPLES ; k++ ) {
    double xloc[ Ndata ] , yloc[ Ndata ] ;
    size_t j ;
    for( j = 0 ; j < Ndata ; j++ ) {
      xloc[j] = dataX[j].resampled[k] ;
      yloc[j] = dataY[j].resampled[k] ;
    }
    Int.resampled[k] = integrator( yloc , xloc , Ndata ) ;
  }
  // and the average
  double xloc[ Ndata ] , yloc[ Ndata ] ;
  size_t j ;
  for( j = 0 ; j < Ndata ; j++ ) {
    xloc[j] = dataX[j].avg ;
    yloc[j] = dataY[j].avg ;
  }
  Int.avg = integrator( yloc , xloc , Ndata ) ;

  compute_err( &Int ) ;
  
  return Int ;
}

// placeholder for numerically integrating a fit
struct resampled
Nint_fit( struct resampled *f ,
	  const struct data_info Data ,
	  const struct fit_info Fit ,
	  const double upp ,
	  const double low ,
	  const double eps ,
	  const size_t shift )
{
  struct resampled Int = init_dist( NULL ,
				    f[0].NSAMPLES ,
				    f[0].restype ) ;
  struct fit_descriptor fdesc = init_fit( Data , Fit ) ;

  // could be done in parallel
  size_t j , p ;
  for( j = 0 ; j < f[0].NSAMPLES ; j++ ) {
    double fparams[ fdesc.Nparam ] ;
    for( p = 0 ; p < fdesc.Nparam ; p++ ) {
      fparams[ p ] = f[ Fit.map[shift].p[p] ].resampled[j] ;
    }
    Int.resampled[j] = adaptive_simpsons( fparams , Data , Fit ,
					  fdesc , low , upp , eps ,
					  shift ) ;
  }
  // do the average
  double fparams[ fdesc.Nparam ] ;
  for( p = 0 ; p < fdesc.Nparam ; p++ ) {
    fparams[ p ] = f[ Fit.map[shift].p[p] ].avg ;
  }
  Int.avg = adaptive_simpsons( fparams , Data , Fit ,
			       fdesc , low , upp , eps ,
			       shift ) ;

  fprintf( stdout , "\n[INT] Integrated fit parameters\n" ) ;
  fprintf( stdout , "[INT] Integration range %e -> %e\n" ,
	   low , upp ) ;
  fprintf( stdout , "[INT] Integral %e %e\n\n" , Int.avg , Int.err ) ;
  
  return Int ;
}
