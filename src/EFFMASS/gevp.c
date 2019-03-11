/**
   @file gevp.c
   @brief solves a generalised eigenvalue problem

   Solves the system

   A.v = \lambda B.v
   
   where A and B are real matrices

   Uses GSL's eigenvalue sort routines for the eigenvalues
 */
#include "gens.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "stats.h"

// #define ABS_EVALUES

struct GEVP_temps {
  gsl_matrix *a ;
  gsl_matrix *b ;
  gsl_eigen_genv_workspace *work ;
  gsl_vector_complex *alpha ;
  gsl_vector *beta ;
  gsl_matrix_complex *evec ;
  double complex *ev ;
  double *evec_norm ;
  size_t N ;
} ;

#ifdef ABS_EVALUES

static int
comp_desc( const double complex ev1 ,
	   const double complex ev2 )
{
  return cabs(ev1) < cabs(ev2) ;
}

static int
comp_asc( const double complex ev1 ,
	  const double complex ev2 )
{
  return cabs(ev1) > cabs(ev2) ;
}

#else

static int
comp_desc( const double complex ev1 ,
	   const double complex ev2 )
{
  return creal(ev1) < creal(ev2) ;
}

static int
comp_asc( const double complex ev1 ,
	  const double complex ev2 )
{
  return creal(ev1) > creal(ev2) ;
}

#endif

// sort by whatever I specify in comp at the moment just eigenvalues
static void
insertion_sort( struct GEVP_temps *G ,
		int (*comp)( const double complex ev1 ,
			     const double complex ev2 ) )
{
  gsl_vector_complex *tevec = gsl_vector_complex_alloc( G -> N ) ;
  
  size_t i , j ;
  for( i = 0 ; i < G -> N ; i++ ) {

    const double complex ev1 = G -> ev[i] ;
    for( j = 0 ; j < G -> N ; j++ ) {
      gsl_vector_complex_set( tevec , j ,
			      gsl_matrix_complex_get( G->evec , j , i ) ) ;
    }
    
    int hole = (int)i - 1 ;
    while( hole >= 0 && comp( G -> ev[hole] , ev1 ) ) {
      // copy data
      G -> ev[hole+1] = G -> ev[hole] ;
      for( j = 0 ; j < G -> N ; j++ ) {
        gsl_matrix_complex_set( G->evec , j , hole+1 ,
				gsl_matrix_complex_get( G->evec , j , hole ) ) ;
      }
      hole-- ;
    }
    G -> ev[hole+1] = ev1 ;
    for( j = 0 ; j < G -> N ; j++ ) {
      gsl_matrix_complex_set( G -> evec , j , i ,
			      gsl_vector_complex_get( tevec , j ) ) ;
    }
  }

  gsl_vector_complex_free( tevec ) ;
  
  return ;
}

// uses GSL to solve a generalised eigenvalue problem
// returns the real part of the eigenvalues
// solves A.v = \lambda B.v
// where A and B are real matrices
static int
gevp( struct GEVP_temps *G ,
      const size_t n ,
      const bool before ,
      const size_t t ,
      const bool write_evalues ) 
{
  int flag = SUCCESS ;
  
  // perform decomposition
  const int err = gsl_eigen_genv( G->a , G->b , G->alpha ,
				  G->beta , G->evec , G->work ) ;

  if( err != 0 ) {
    fprintf( stderr , "%s\n" , gsl_strerror( err ) ) ;
    fprintf( stderr , "GEVP Aborting\n" ) ;
    flag = FAILURE ;
  }

  // set the ev
  size_t i ;
  for( i = 0 ; i < G -> N ; i++ ) {
    G -> ev[i] = ( gsl_vector_complex_get( G -> alpha , i ).dat[0] 
		   + I*gsl_vector_complex_get( G -> alpha , i ).dat[1] )
      / gsl_vector_get( G->beta , i ) ;
  }

  // sort by the real part
  if( before ) {
    insertion_sort( G , comp_asc ) ;
  } else {
    insertion_sort( G , comp_desc ) ;
  }
  
  size_t j ;
  for( j = 0 ; j < G -> N ; j++ ) {
    if( write_evalues ) {
      // if we want them written
      fprintf( stdout , "STATE_%zu %zu evalue %e %e \n" ,
	       j , t , creal( G -> ev[j]) , cimag( G -> ev[j]) ) ;
      size_t nevec ;
      for( nevec = 0 ; nevec < G->N ; nevec++ ) {
	fprintf( stdout , "STATE_%zu %zu evec %zu %e %e\n" ,
		 j , t , nevec ,
		 ( gsl_matrix_complex_get( G -> evec , nevec , j ).dat[0] )
		 ,
		 ( gsl_matrix_complex_get( G -> evec , j , nevec ).dat[1] )
		 ) ;
      }
      //
    }
  }
  
  return flag ;
}

// take the SVD
static int 
svd_TLS( gsl_matrix *a ,
	 gsl_matrix *b ,
	 const double *C0 ,
	 const double *C1 ,
	 const size_t NROWS ,
	 const size_t NCOLS )
{
  // allocations
  gsl_matrix *U  = gsl_matrix_alloc( 2*NROWS , NCOLS ) ;
  gsl_matrix *Q  = gsl_matrix_alloc( NCOLS , NCOLS ) ;
  gsl_vector *S  = gsl_vector_alloc( NCOLS ) ;
  gsl_vector *WORK  = gsl_vector_alloc( NCOLS ) ;

  size_t i , j ;
  int FLAG = SUCCESS ;

  // initialize our matrix
  for( i = 0 ; i < NROWS ; i++ ){
    for( j = 0 ; j < NCOLS ; j++ ) {
      gsl_matrix_set( U , i , j , C0[j+i*NCOLS] ) ;
      gsl_matrix_set( U , i + NROWS , j , C1[j+i*NCOLS] ) ;
    }
  }

  // perform the decomposition
  if( gsl_linalg_SV_decomp_jacobi( U , Q , S ) ) {
    printf( "[SVD] GSL SVD comp failure \n" ) ;
    FLAG = FAILURE ;
    goto FREE ;
  }

  for( i = 0 ; i < NCOLS ; i++ ) {
    for( j = 0 ; j < NCOLS ; j++ ) {
      gsl_matrix_set( a , i , j , gsl_matrix_get( U , i , j ) ) ;
      gsl_matrix_set( b , i , j , gsl_matrix_get( U , i + NROWS , j ) ) ;
    }
  }

 FREE :
  gsl_matrix_free( U ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( WORK ) ;

  return 0 ;
}

// set the a and b matrices
static void
set_ab( gsl_matrix *a , gsl_matrix *b ,
	const double *C0 , const double *C1 ,
	const size_t N , const size_t M )
{
  size_t i , j ;
  if( N == M ) {
    for( i = 0 ; i < N ; i++ ) {
      for( j = 0 ; j < N ; j++ ) {
	gsl_matrix_set( a , i , j , *C0 ) ; C0++ ;
	gsl_matrix_set( b , i , j , *C1 ) ; C1++ ;
      }
    }
  } else {
    svd_TLS( a , b , C0 , C1 , M , N ) ;
  }
  return ;
}

static void
allocate_G( struct GEVP_temps *G ,
	    const size_t N )
{
  // set the data
  G -> a = gsl_matrix_alloc( N , N ) ;
  G -> b = gsl_matrix_alloc( N , N ) ;
  // allocations for GSL
  G -> work  = gsl_eigen_genv_alloc( N ) ;
  G -> alpha = gsl_vector_complex_alloc( N ) ;
  G -> beta  = gsl_vector_alloc( N ) ;
  G -> evec  = gsl_matrix_complex_alloc( N , N ) ;
  G -> evec_norm = malloc( N*sizeof(double) ) ;
  G -> ev = malloc( N*sizeof( double complex ) ) ;
  G -> N = N ;
}

static void
free_G( struct GEVP_temps *G )
{
  gsl_matrix_complex_free( G -> evec ) ;
  gsl_vector_complex_free( G -> alpha ) ;
  gsl_vector_free( G -> beta ) ;
  gsl_eigen_genv_free( G-> work ) ;
  gsl_matrix_free( G -> a ) ;
  gsl_matrix_free( G -> b ) ;
  free( G -> ev ) ;
  free( G -> evec_norm ) ;
}

// solve the GEVP
struct resampled *
solve_GEVP( const struct resampled *y ,
	    const size_t Ndata ,
	    const size_t N ,
	    const size_t M ,
	    const size_t t0 )
{
  if( N > M ) {
    fprintf( stderr , "[GEVP] cannot solve when N states "
	     "are less than M correlators\n" ) ;
    return NULL ;
  }
  
  // initialise the generalised eigenvalues
  struct resampled *evalues = malloc( Ndata*N*
				      sizeof( struct resampled ) ) ;
  size_t i , j , k ;
  for( j = 0 ; j < Ndata*N ; j++ ) {
    evalues[j].resampled = malloc( y[0].NSAMPLES *
				   sizeof( double ) ) ;
    evalues[j].restype = y[0].restype ;
    evalues[j].NSAMPLES = y[0].NSAMPLES ;
  }

  struct GEVP_temps *G = malloc( Ndata*sizeof( struct GEVP_temps ) ) ;

  double C0[ N*M ] , C1[ N*M ] ;
  
  // ugh loop order is all weird
  for( j = 0 ; j < Ndata ; j++ ) {
    allocate_G( &G[j] , N );
  }

  // loop resamples
  for( k = 0 ; k < y[0].NSAMPLES ; k++ ) {
    
    for( j = 0 ; j < Ndata ; j++ ) {

      // put these in the linearised matrices
      size_t shift = 0 ;
      for( i = 0 ; i < N*M ; i++ ) {
	C0[i] = y[ j  + shift ].resampled[k] ;
	C1[i] = y[ t0 + shift ].resampled[k] ;
	shift += Ndata ;
      }

      // set AB by squaring the matrices, possibly
      set_ab( G[j].a , G[j].b , C0 , C1 , N , M ) ;
      
      // compute eigenvalues
      if( gevp( &G[j] , N , j<t0 , j , false ) == FAILURE ) {
	fprintf( stderr , "[GEVP] GEVP solve failed \n" ) ;
        return NULL ;
      }
    }

    for( j = 0 ; j < Ndata ; j++ ) {
      // poke into solution
      for( i = 0 ; i < N ; i++ ) {
	evalues[j+Ndata*i].resampled[ k ] = creal( G[j].ev[i] ) ;
      }
      //
    }
  }

  for( j = 0 ; j < Ndata ; j++ ) {
    // redo for the average
    size_t shift = 0 ;
    for( i = 0 ; i < N*M ; i++ ) {
      C0[i] = y[ j  + shift ].avg ;
      C1[i] = y[ t0 + shift ].avg ;
      shift += Ndata ;
    }

    // set AB by squaring the matrices possibly
    set_ab( G[j].a , G[j].b , C0 , C1 , N , M ) ;

    // compute eigenvalues
    if( gevp( &G[j] , N , j<t0 , j , true ) == FAILURE ) {
      fprintf( stderr , "[GEVP] GEVP solve failed \n" ) ;
      return NULL ;
    }
  }

  for( j = 0 ; j < Ndata ; j++ ) {
    // poke into solution
    for( i = 0 ; i < N ; i++ ) {
      evalues[j+Ndata*i].avg = creal( G[j].ev[i] ) ;
      compute_err( &evalues[j+Ndata*i] ) ;
    }
  }
  
  for( j = 0 ; j < Ndata ; j++ ) {
    free_G( &G[j] ) ;
  }
  free( G ) ;
  
  return evalues ;
}
