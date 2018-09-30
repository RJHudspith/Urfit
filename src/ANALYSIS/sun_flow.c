/**
   @file sun_flow.c
   @brief gradient flow analysis
 */
#include "gens.h"

#include "fit_chooser.h"
#include "fit_and_plot.h"
#include "ffunction.h"
#include "GLU_bswap.h"
#include "init.h"
#include "resampled_ops.h"
#include "stats.h"

//#define LT

#define Ntau (5)

static struct resampled
read_single( const char *filename )
{
  struct resampled sample ;
  sample.resampled = NULL ;
  FILE *file = fopen( filename , "rb" ) ;
  if( file == NULL ) {
    return sample ;
  }
  uint32_t num_mom[1] ;
  if( fread( num_mom , sizeof(uint32_t) , 1 , file ) != 1 ) {
    printf( "[IO] num_mom reading failure %u\n" , num_mom[0] ) ;
  }
#ifndef WORDS_BIGENDIAN
  bswap_32( 1 , num_mom ) ;
#endif
  // read a single momentum
  uint32_t mom[4] ;
  if( fread( mom , sizeof(uint32_t) , 4 , file ) != 4 ) { 
    printf( "[IO] mom reading failure \n" ) ;
  }
#ifndef WORDS_BIGENDIAN
  bswap_32( 4 , mom ) ;
#endif
  
  uint32_t vals[2] ;
  if( fread( vals , sizeof(uint32_t) , 2 , file ) != 2 ) {
    fprintf( stderr , "[IO] Vals reading failure \n" ) ;
    return sample ;
  }
#ifndef WORDS_BIGENDIAN
  bswap_32( 2 , vals ) ;
#endif

  switch( vals[0] ) {
  case 0 : sample.restype  = Raw ; break ;
  case 1 : sample.restype  = JackKnife ; break ;
  case 3 : sample.restype  = BootStrap ; break ;
  }
  sample.NSAMPLES = vals[1] ;

  sample.resampled = malloc( sample.NSAMPLES * sizeof( double ) ) ;

  if( fread( sample.resampled , sizeof(double) ,
	     sample.NSAMPLES , file ) != sample.NSAMPLES ) {
    fprintf( stderr , "[SUN_T0] Sample read error\n" ) ;
    return sample ;
  }
#ifndef WORDS_BIGENDIAN
    bswap_64( sample.NSAMPLES , sample.resampled ) ;
#endif
  
   double average[1] ;
  if( fread( average , sizeof(double) , 1 , file ) != 1 ) {
    fprintf( stderr , "[SUN_T0] Average read failure \n" ) ;
    return sample ;
  }
#ifndef WORDS_BIGENDIAN
  bswap_64( 1 , average ) ;
#endif
  const double compar = average[ 0 ] ;
  compute_err( &sample ) ;
  sample.avg = compar ;

  fprintf( stdout , "[SUN_T0] %f %f \n" , sample.avg , sample.err ) ;

  fclose(file) ;

  return sample ;
}

// do cont -> write to flat as per usual
static int
do_cont( const struct input_params *Input ,
	 const struct resampled *fit ,
	 const struct resampled sub ,
	 const size_t NC ,
	 const size_t shift )
{
  const size_t Lt[ Ntau ] = { 6 , 7 , 8 , 9 , 10 } ;

  char str[ 256 ] ;
  sprintf( str , "flat_cont_SU%zu.dat" , NC ) ;
    
  FILE *outfile = fopen( str , "w" ) ;

  fprintf( outfile , "%d\n" , fit[0].restype ) ;
  fprintf( outfile , "%d\n" , Ntau ) ;
  
  size_t nt , p , i ;
  for( nt = 0 ; nt < Ntau ; nt++ ) {

    if( NC < 6 ) {
      sprintf( str , "NS%zu_NT%zu_crossing_SU%zu.bin" , 3*Lt[nt] , Lt[nt] , NC ) ;
    } else {
      sprintf( str , "NS%zu_NT%zu_crossing_SU%zu.bin" , 2*Lt[nt] , Lt[nt] , NC ) ;
    }
  
    // read in some data
    struct resampled this = read_single( str ) ;

    if( this.resampled == NULL ) {
      fprintf( stderr , "FILE %s not found \n" , str ) ;
      continue ;
    }

    // set up the fit again
    struct fit_descriptor fdesc = init_fit( Input -> Data , Input -> Fit ) ;
    double fparams[ fdesc.Nparam ] ;

    printf( "fit input \n" ) ;

    subtract( &this , sub ) ;
    divide_constant( &this , NC*NC ) ;

    printf( "Shifted %f %f\n" , this.avg , this.err ) ;
    struct resampled data = init_dist( NULL ,
				       fit[0].NSAMPLES ,
				       fit[0].restype ) ;

    for( i = 0 ; i < fit[0].NSAMPLES ; i++ ) {
      
      // evaluate the fitfunc
      struct x_desc xdesc = { this.resampled[i] ,
			      Input -> Data.LT[0] ,
			      Input -> Fit.N , Input -> Fit.M } ;

      for( p = 0 ; p < fdesc.Nparam ; p++ ) {
	fparams[ p ] = fit[p].resampled[i] ;
      }
      data.resampled[i] = fdesc.func( xdesc , fparams , 0 ) ;
    }

    struct x_desc xdesc = { this.avg ,
			    Input -> Data.LT[0] ,
			    Input -> Fit.N , Input -> Fit.M } ;
    // and the average
    for( p = 0 ; p < fdesc.Nparam ; p++ ) {
      fparams[ p ] = fit[p].avg ;
    }
    data.avg = fdesc.func( xdesc , fparams , 0 ) ;

    compute_err( &data ) ;

    struct resampled t0tc = init_dist( &data ,
				       data.NSAMPLES ,
				       data.restype ) ;
	  
    mult_constant( &t0tc , 1/(double)Lt[nt] ) ;
    raise( &data , -2 ) ;
 
    printf( "Result :: %f %f %f %f \n" , data.avg , t0tc.avg , data.err , t0tc.err ) ;

    // write out a flat file
    fprintf( outfile , "%zu\n" , fit[0].NSAMPLES ) ;
    for( i = 0 ; i < data.NSAMPLES ; i++ ) {
      fprintf( outfile , "%1.12e %1.12e\n" , data.avg , t0tc.resampled[i] ) ;
    }
    fprintf( outfile , "AVG %1.12e %1.12e\n" , data.avg , t0tc.avg ) ;

    free( t0tc.resampled ) ;
  
    // free the data
    free( this.resampled ) ;
    free( data.resampled ) ;
  }

  // appended to flat file
  fclose( outfile ) ;

  return SUCCESS ;
}

int
sun_wflow_analysis( struct input_params *Input )
{
  // shift the data x = x - x[0]
  size_t i , j , shift = 0 ;

  const double NC[ 6 ] = { 3 , 4 , 5 , 6 , 7 , 8 } ;
  const double NCSQUARED[ 6 ] = { (3*3.) , (4*4.) , (5*5.) ,
				  (6*6.) , (7*7.) , (8*8.) } ;
  
  struct resampled *sub = malloc( Input -> Data.Nsim * sizeof( struct resampled ) ) ;
  
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {

     sub[i] = init_dist( &(Input -> Data.x[shift+Input -> Fit.M]) ,
			 Input -> Data.x[shift+Input -> Fit.M].NSAMPLES , 
			 Input -> Data.x[shift+Input -> Fit.M].restype ) ;

    #ifndef LT
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      subtract( &(Input -> Data.x[j]) , sub[i] ) ;
      divide_constant( &(Input -> Data.x[j]) , NCSQUARED[i] ) ;
    }
    
    #else
    
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {
      raise( &(Input -> Data.y[j]) , -1 ) ;
      mult_constant( &(Input -> Data.y[j]) , Input -> Traj[i].Dimensions[3] ) ;
      printf( "LT :: %f %f \n" , Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
    }
    
    #endif
    
    shift = j ;
  }


  // perform the fit
  double chi = 0.0 ;
  struct resampled *fit = fit_and_plot( *Input , &chi ) ;

  shift = 0 ;
  for( i = 0 ; i < 4 ; i++ ) {
    do_cont( Input , fit , sub[i] , NC[i] , shift ) ;
    shift += Input -> Data.Ndata[i] ;
  }
  
  // free the fit
  free_fitparams( fit , Input -> Fit.Nlogic ) ;

  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    free( sub[i].resampled ) ;
  }
  free( sub ) ;
  
  return SUCCESS ;
}
