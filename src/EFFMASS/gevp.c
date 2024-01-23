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
#include <gsl/gsl_version.h>

#include "stats.h"

//#define ABS_EVALUES

#define OPTIMISED_CORRELATOR

//#define SYMMETRIC_GEVP

struct GEVP_temps {
  gsl_matrix *a ;
  gsl_matrix *b ;
#ifdef SYMMETRIC_GEVP
  gsl_eigen_gensymmv_workspace *work ;
  double *ev ;
  gsl_matrix *evec ;
  gsl_vector *eval ;
#else
  gsl_eigen_genv_workspace *work ;
  gsl_vector_complex *alpha ;
  gsl_vector *beta ;
  gsl_matrix_complex *evec ;
  double complex *ev ;
#endif
  double *evec_norm ;
  size_t N ;
} ;

#ifdef SYMMETRIC_GEVP
  static int
  comp_desc( const double ev1 ,
	     const double ev2 )
  {
    return ev1 < ev2 ;
  }
  static int
  comp_asc( const double ev1 ,
	    const double ev2 )
  {
    return ev1 > ev2 ;
  }
#else
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
#endif

// sort by whatever I specify in comp at the moment just eigenvalues
static void
insertion_sort( struct GEVP_temps *G ,
		#ifdef SYMMETRIC_GEVP
		int (*comp)( const double ev1 ,
			     const double ev2 )
		#else
		int (*comp)( const double complex ev1 ,
			     const double complex ev2 )
		#endif
		)
{
  #ifdef SYMMETRIC_GEVP
  gsl_vector *tevec = gsl_vector_alloc( G -> N ) ;
  #else
  gsl_vector_complex *tevec = gsl_vector_complex_alloc( G -> N ) ;
  #endif
  
  size_t i , j ;
  for( i = 0 ; i < G -> N ; i++ ) {

    const double complex ev1 = G -> ev[i] ;
    for( j = 0 ; j < G -> N ; j++ ) {
#ifdef SYMMETRIC_GEVP
      gsl_vector_set( tevec , j , gsl_matrix_get( G->evec , j , i ) ) ;

#else
      gsl_vector_complex_set( tevec , j ,
			      gsl_matrix_complex_get( G->evec , j , i ) ) ;
#endif
    }
    
    int hole = (int)i - 1 ;
    while( hole >= 0 && comp( G -> ev[hole] , ev1 ) ) {
      // copy data
      G -> ev[hole+1] = G -> ev[hole] ;
      for( j = 0 ; j < G -> N ; j++ ) {
#ifdef SYMMETRIC_GEVP
        gsl_matrix_set( G->evec , j , hole+1 ,
			gsl_matrix_get( G->evec , j , hole ) ) ;
#else
        gsl_matrix_complex_set( G->evec , j , hole+1 ,
				gsl_matrix_complex_get( G->evec , j , hole ) ) ;
#endif
      }
      hole-- ;
    }
    G -> ev[hole+1] = ev1 ;
    for( j = 0 ; j < G -> N ; j++ ) {
#ifdef SYMMETRIC_GEVP
      gsl_matrix_set( G -> evec , j , hole+1 , gsl_vector_get( tevec , j ) ) ;
#else
      gsl_matrix_complex_set( G -> evec , j , hole+1 ,
			      gsl_vector_complex_get( tevec , j ) ) ;
#endif
    }
  }

#ifdef SYMMETRIC_GEVP
  gsl_vector_free( tevec ) ;
#else
  gsl_vector_complex_free( tevec ) ;
#endif
  return ;
}

#ifdef OPTIMISED_CORRELATOR

// compute the "optimised correlator" from the eigenvectors
static void
optimised_correlator( struct GEVP_temps *G ,
		      const struct resampled *y ,
		      const size_t k ,
		      const bool is_avg ,
		      const size_t N ,
		      const size_t M ,
		      const size_t Ndata ,
		      const size_t td )
{
  // check diagonalisation with V^\dag C(t) V
  double C0[ N*M ] ;
  double complex V2[N][N][N];
  for( size_t a = 0 ; a < N ; a++ ) {
    for( size_t b = 0 ; b < N ; b++ ) {
      const double complex v1 = conj( gsl_matrix_complex_get( G[td].evec , b , a ) ) ;
      for( size_t c = 0 ; c < N ; c++ ) {
	const double complex v2 = gsl_matrix_complex_get( G[td].evec , c , a ) ;
	V2[a][b][c] = v1*v2 ;
      }
    }
  }
    
  for( size_t j = 0 ; j < Ndata ; j++ ) {    
    // put these in the linearised matrices
    size_t shift = 0 ;
    for( size_t i = 0 ; i < N*M ; i++ ) {
      if( is_avg == true ) {
	C0[i] = y[ j  + shift ].avg ;
      } else {
	C0[i] = y[ j  + shift ].resampled[k] ;
      }
      shift += Ndata ;
    }
    // diagonalise
    double complex *p1 = (double complex*)V2 ;
    for( size_t a = 0 ; a < N ; a++ ) {
      double *p2 = (double*)C0 ;
      register double sum = 0.0 ;
      for( size_t i = 0 ; i < N*N ; i++ ) {
	sum += creal((*p1))*(*p2) ;
	p1++ ; p2++ ;
      }
      G[j].ev[a] = creal(sum) ;
    }
  }
  return ;
}

// compute the "optimised correlator" from the eigenvectors
static void
optimised_correlator2( struct GEVP_temps *G1 ,
		       const struct GEVP_temps G2 ,
		       const double *C0 ,
		       const size_t N )
{
  size_t a , b , c ;
  for( a = 0 ; a < N ; a++ ) {
#ifdef SYMMETRIC_GEVP
    register double Sum1 = 0. ;
#else
    register double complex Sum1 = 0. ;
#endif
    for( b = 0 ; b < N ; b++ ) {
      #ifdef SYMMETRIC_GEVP
      register double Sum1 = 0. ;
      const double A = gsl_matrix_get( G2.evec , b , a ) ;
      const double AC = A ;
      #else
      register double complex Sum1 = 0. ;
      const gsl_complex A = gsl_matrix_complex_get( G2.evec , b , a ) ;
         #if GSL_MINOR_VERSION > 6
         const double complex AC = conj( A ) ;
         #else
         const double complex AC = A.dat[0] - I*A.dat[1] ;
         #endif
      #endif
      for( c = 0 ; c < N ; c++ ) {
	#ifdef SYMMETRIC_GEVP
	const double B = gsl_matrix_get( G2.evec , c , a ) ;
	const double BC = B ;
        #else
	const gsl_complex B = gsl_matrix_complex_get( G2.evec , c , a ) ;
           #if GSL_MINOR_VERSION > 6
	   const double complex BC = B ;
           #else
	   const double complex BC = B.dat[0] + I*B.dat[1] ;
           #endif
	#endif
	Sum1 += AC * BC * C0[ c + N*b ] ;
      }
    }
    #ifdef SYMMETRIC_GEVP
    G1 -> ev[a] = Sum1 ;
    #else
    G1 -> ev[a] = creal( Sum1 ) ;
    #endif
  }
}
#endif

// uses GSL to solve a generalised eigenvalue problem
// returns the real part of the eigenvalues
// solves A.v = \lambda B.v
// where A and B are real matrices
static int
gevp2( struct GEVP_temps *G ,
       const size_t n ,
       const bool before ,
       const size_t t ,
       const bool write_evalues ) 
{
  int flag = SUCCESS ;
  
  // perform decomposition
#ifdef SYMMETRIC_GEVP
  const int err = gsl_eigen_gensymmv( G[t].a , G[t].b ,
				      G[t].eval , G[t].evec , G[t].work ) ;
#else
  const int err = gsl_eigen_genv( G[t].a , G[t].b , G[t].alpha ,
				  G[t].beta , G[t].evec , G[t].work ) ;
#endif
  if( err != 0 ) {
    fprintf( stderr , "%s\n" , gsl_strerror( err ) ) ;
    fprintf( stderr , "GEVP Aborting\n" ) ;
    flag = FAILURE ;
  }

#ifndef OPTIMISED_CORRELATOR
  // set the ev
  size_t i ;
  for( i = 0 ; i < G -> N ; i++ ) {
    #ifdef SYMMETRIC_GEVP
    G[t].ev[i] = gsl_vector_get( G[t].eval , i ) ;
    #else
    G[t].ev[i] = gsl_vector_complex_get( G[t].alpha , i )/gsl_vector_get( G[t].beta , i ) ;
    if( creal( G[t].ev[i] ) < 0.0 ) {
      G[t].ev[i] = 1E-16 ;
    }
    #endif
  }

  // sort by the real part
  if( before ) {
    insertion_sort( &G[t] , comp_asc ) ;
  } else {
    insertion_sort( &G[t] , comp_desc ) ;
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
		 #ifdef SYMMETRIC_GEVP
		 ( gsl_matrix_get( G -> evec , nevec , j ) ) , 0.0 
		 #else
		 ( creal( gsl_matrix_complex_get( G -> evec , nevec , j ) ) )
		 ,
		 ( cimag( gsl_matrix_complex_get( G -> evec , j , nevec ) ) )
		 #endif
		 ) ;
      }
      //
    }
  }
#endif

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
#ifdef SYMMETRIC_GEVP
  G -> work  = gsl_eigen_gensymmv_alloc( N ) ;
  G -> evec  = gsl_matrix_alloc( N , N ) ;
  G -> evec_norm = malloc( N*sizeof(double) ) ;
  G -> eval = gsl_vector_alloc( N ) ;
  G -> ev = malloc( N*sizeof(double) ) ;
#else
  // allocations for GSL
  G -> work  = gsl_eigen_genv_alloc( N ) ;
  G -> alpha = gsl_vector_complex_alloc( N ) ;
  G -> beta  = gsl_vector_alloc( N ) ;
  G -> evec  = gsl_matrix_complex_alloc( N , N ) ;
  G -> evec_norm = malloc( N*sizeof(double) ) ;
  G -> ev = malloc( N*sizeof( double complex ) ) ;
#endif
  G -> N = N ;
}

static void
free_G( struct GEVP_temps *G )
{
#ifdef SYMMETRIC_GEVP
  gsl_matrix_free( G -> evec ) ;
  gsl_eigen_gensymmv_free( G-> work ) ;
  gsl_vector_free( G->eval ) ;
  free( G -> evec_norm ) ;
#else
  gsl_matrix_complex_free( G -> evec ) ;
  gsl_vector_complex_free( G -> alpha ) ;
  gsl_vector_free( G -> beta ) ;
  gsl_eigen_genv_free( G-> work ) ;
  gsl_matrix_free( G -> a ) ;
  gsl_matrix_free( G -> b ) ;
  free( G -> ev ) ;
  free( G -> evec_norm ) ;
#endif
}

// solve the GEVP
struct resampled *
solve_GEVP( const struct resampled *y ,
	    const size_t Ndata ,
	    const size_t N ,
	    const size_t M ,
	    const size_t t0 ,
	    const size_t td )
{
  fprintf( stdout , "GEVP solver ----> " ) ;
  
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
      if( gevp2( G , N , j<t0 , j , false ) == FAILURE ) {
	fprintf( stderr , "[GEVP] GEVP solve failed \n" ) ;
        return NULL ;
      }
    }

#ifdef OPTIMISED_CORRELATOR
    optimised_correlator( G , y , k , false , N , M , Ndata , td ) ;
#endif
    
    for( j = 0 ; j < Ndata ; j++ ) {
      // poke into solution
      for( i = 0 ; i < N ; i++ ) {
	#ifdef SYMMETRIC_GEVP
	evalues[j+Ndata*i].resampled[ k ] = gsl_vector_get( G[j].eval , i ) ;
	#else
	evalues[j+Ndata*i].resampled[ k ] = creal( G[j].ev[i] ) ;
	#endif
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

    // print the matrices
    for( i = 0 ; i < N ; i++ ) {
      for( size_t k = 0 ; k < M ; k++ ) {
	printf( "%e ", C0[k+i*M] ) ;
      }
      printf( "\n" ) ;
    }
    printf( "\n" ) ;
    // print the matrices
    for( i = 0 ; i < N ; i++ ) {
      for( size_t k = 0 ; k < M ; k++ ) {
	printf( "%e ", C1[k+i*M] ) ;
      }
      printf( "\n" ) ;
    }
    printf( "\n" ) ;
    
    // set AB by squaring the matrices possibly
    set_ab( G[j].a , G[j].b , C0 , C1 , N , M ) ;

    // compute eigenvalues
    if( gevp2( G , N , j<t0 , j , true ) == FAILURE ) {
      fprintf( stderr , "[GEVP] GEVP solve failed \n" ) ;
      return NULL ;
    }
  }
#ifdef OPTIMISED_CORRELATOR
  optimised_correlator( G , y , 0 , true , N , M , Ndata , td ) ;
#endif
  
  for( j = 0 ; j < Ndata ; j++ ) {
    // poke into solution
    for( i = 0 ; i < N ; i++ ) {
      #ifdef SYMMETRIC_GEVP
      evalues[j+Ndata*i].avg = gsl_vector_get( G[j].eval , i ) ;
      #else
      evalues[j+Ndata*i].avg = creal( G[j].ev[i] ) ;
      #endif
      compute_err( &evalues[j+Ndata*i] ) ;
    }
  }
  
  for( j = 0 ; j < Ndata ; j++ ) {
    free_G( &G[j] ) ;
  }
  free( G ) ;
  
  return evalues ;
}

// solve the GEVP
struct resampled *
solve_GEVP_fixed( const struct resampled *y ,
		  const size_t Ndata ,
		  const size_t N ,
		  const size_t M ,
		  const size_t t0 ,
		  const size_t td )
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
	C0[i] = y[ j + shift ].resampled[k] ;
	C1[i] = y[ (j+t0)%Ndata + shift ].resampled[k] ;
	shift += Ndata ;
      }

      // set AB by squaring the matrices, possibly
      set_ab( G[j].a , G[j].b , C0 , C1 , N , M ) ;
      
      // compute eigenvalues
      if( gevp2( G , N , true , j , false ) == FAILURE ) {
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
      C1[i] = y[ (j+t0)%Ndata + shift ].avg ;
      shift += Ndata ;
    }

    // set AB by squaring the matrices possibly
    set_ab( G[j].a , G[j].b , C0 , C1 , N , M ) ;

    // compute eigenvalues
    if( gevp2( G , N , true , j , true ) == FAILURE ) {
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
