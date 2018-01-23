/**
   @file beta_crit.c
   @brief critical beta evaluation
 */
#include "gens.h"

#include "fit_and_plot.h"
#include "GLU_bswap.h"
#include "histogram.h"
#include "init.h"
#include "NR.h"
#include "stats.h"
#include "rng.h"

#define Q (5)

// quicksort comparator
static int
cmp( const void *a ,
     const void *b )
{
  return *(double*)a < *(double*)b ? -1 : \
    ( *(double*)a > *(double*)b ? 1 : 0 ) ;
}

// the separatrix
static double
compute_sep( const double *data ,
	     const int Ndata ,
	     const double cutoff )
{
  double *sdata = malloc( Ndata * sizeof( double ) ) ;
  memcpy( sdata , data , Ndata * sizeof( double ) ) ;
  
  // sort the data
  qsort( sdata , Ndata , sizeof( double ) , cmp ) ;
  
  size_t i , Ncount = 0 ;
  for( i = 0 ; i < Ndata ; i++ ) {
    Ncount++ ;
    if( sdata[i] > cutoff ) break ;
  }

  free( sdata ) ;
  
  return ( (double)( Q + 1 ) * Ncount - (double)Ndata )
    / ( (double)( Q - 1 ) * Ncount + (double)Ndata ) ;
}

// test whether various evaluations are equivalent within bin widths
static bool
are_overlapping( const double mid1 , const double err1 ,
		 const double mid2 , const double err2 )
{
  if( mid1 > mid2 ) {
    if( ( mid1 - err1 ) < ( mid2 + err2 ) ) {
      return true ;
    }
  } else {
    if( ( mid1 + err1 ) > ( mid2 - err2 ) ) {
      return true ;
    }
  }
  return false ;
}

static double
histmin( double *width ,
	 const struct resampled y ,
	 const size_t Nbin )
{
  struct histogram *hist = histogram( y.resampled ,
				      y.NSAMPLES ,
				      Nbin ) ;

  size_t min ;
  const double x = hist_min( &min , hist , Nbin ) ;

  *width = hist[min].width ;

  free( hist ) ;
  
  return x ;
}

int
beta_crit( struct input_params *Input )
{
  // initialise the rng to some value
  init_rng( 123456 ) ;

  double sep_best = 1.0 , beta_best = Input -> Data.x[0].avg ;
  
  size_t i , j , k , shift = 0 ;
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    for( j = shift ; j < shift + Input -> Data.Ndata[i] ; j++ ) {

      double width , midprev = 0.0 , errprev = 0.0 ;
      double mid = histmin( &width , Input -> Data.y[j] , 20 ) ;
      double best = mid , widthbest = width ;

      size_t nconsistent = 0 , bestconsistent = 0 , Nbin ;
      for( Nbin = 15 ; Nbin < 27 ; Nbin++ ) {

        mid = histmin( &width , Input -> Data.y[j] , Nbin ) ;

	if( are_overlapping( mid , width/2. , midprev , errprev/2. ) ) {
	  nconsistent++ ;
	} else {
	  if( nconsistent > bestconsistent ) {
	    best = midprev ;
	    widthbest = errprev ;
	    bestconsistent = nconsistent ;
	  }
	}
	midprev = mid ;
	errprev = width ;
      }

      if( nconsistent > bestconsistent ) {
	best = midprev ;
	widthbest = errprev ;
	bestconsistent = nconsistent ;
      }

      fprintf( stdout , "\n%zu %f %f \n" , Nbin , best , widthbest/2. ) ;

      const double sep = compute_sep( Input -> Data.y[j].resampled  ,
				      Input -> Data.y[j].NSAMPLES ,
				      best ) ;
      const double seperr = 0.5 * ( compute_sep( Input -> Data.y[j].resampled  ,
						 Input -> Data.y[j].NSAMPLES ,
						 best + widthbest/2. ) -
				    compute_sep( Input -> Data.y[j].resampled  ,
						 Input -> Data.y[j].NSAMPLES ,
						 best - widthbest/2. ) ) ;

      fprintf( stdout , "SEP %f %f %f \n" , Input -> Data.x[j].avg , sep , seperr ) ;


      // reallocate x and y and create gaussian distributions
      Input -> Data.x[j].resampled =		\
	realloc( Input -> Data.x[j].resampled ,
		 Input -> Data.Nboots * sizeof( double ) ) ;
      Input -> Data.y[j].resampled =		\
	realloc( Input -> Data.y[j].resampled ,
		 Input -> Data.Nboots * sizeof( double ) ) ;

      // set y to be some gaussian width
      Input -> Data.y[j].avg = sep ;
      
      Input -> Data.y[j].NSAMPLES = Input -> Data.Nboots ;
      Input -> Data.y[j].restype = BootStrap ;
      
      Input -> Data.x[j].NSAMPLES = Input -> Data.Nboots ;
      Input -> Data.x[j].restype = BootStrap ;
      
      //rng_reseed() ;
      for( k = 0 ; k < Input -> Data.Nboots ; k++ ) {
        Input -> Data.y[j].resampled[k] = sep + rng_gaussian( seperr ) ;
      }

      compute_err( &Input -> Data.y[j] ) ;
      compute_err( &Input -> Data.x[j] ) ;

      #ifdef VERBOSE
      fprintf( stdout , "%f %f %f %f \n" ,
	Input -> Data.x[j].avg , Input -> Data.x[j].err ,
	Input -> Data.y[j].avg , Input -> Data.y[j].err ) ;
      #endif
      
      if( fabs( sep ) < sep_best ) {
	sep_best  = Input -> Data.y[j].avg ;
	beta_best = Input -> Data.x[j].avg ;
      }
    }
    shift = j ;
  }

  fprintf( stdout , "[BC] SEP closest to zero %f %f \n", beta_best , sep_best ) ;
  
  free_rng( ) ;

  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;

  if( Input -> Fit.Fitdef != NOFIT ) {
    struct resampled zero = fit_zero( fit , Input , beta_best ) ;
    
    fprintf( stdout , "[BC] ZERO prediction :: %f %f \n" , zero.avg , zero.err ) ;

    // write out a crossing file
    char str[ 256 ] ;
    sprintf( str , "NS%zu_NT%zu_crossing.bin" ,
	     Input -> Traj[0].Dimensions[0] ,
	     Input -> Traj[0].Dimensions[3] ) ;
    FILE *file = fopen( str , "wb" ) ;

    uint32_t num_mom[1] = { 1 } ;
    int32_t mom[ 4 ] = { 3 , 0 , 0 , 0 } ;
#ifndef WORDS_BIGENDIAN
    bswap_32( 1 , num_mom ) ;
#endif
    fwrite( num_mom , sizeof( uint32_t ) , 1 , file ) ;
#ifndef WORDS_BIGENDIAN
    bswap_32( 4 , mom ) ;
#endif
    fwrite( mom , sizeof( int32_t ) , 4 , file ) ;

    
    int32_t vals[2] = { 3 , zero.NSAMPLES } ;
#ifndef WORDS_BIGENDIAN
    bswap_32( 2 , vals ) ;
#endif
    fwrite( vals , sizeof( uint32_t ) , 2 , file ) ;

#ifndef WORDS_BIGENDIAN
    bswap_64( zero.NSAMPLES , zero.resampled ) ;
#endif
    fwrite( zero.resampled , sizeof( double ) , zero.NSAMPLES , file ) ;
    
    double avg[1] = { zero.avg } ;
#ifndef WORDS_BIGENDIAN
    bswap_64( 1 , avg ) ;
#endif
    fwrite( avg , sizeof( double ) , 1 , file ) ;
    
    fclose( file ) ;
    
    free( zero.resampled ) ;
  }
  
  free_fitparams( fit , Input -> Fit.Nlogic ) ;
  
  return SUCCESS ;
}
