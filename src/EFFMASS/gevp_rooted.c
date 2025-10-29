/**
   @file gevp_rooted.c
   @brief solves a rooted (Colin-like) Correlator problem by first computing the matrix root

   C1 = \sqrt(C[t_0])

   dressing C

   C[t] = C1.C[t].C1

   computin the evecs for this matrix at some time t_D "U"

   and diagonalising to get the Evs

   \lambda[t][a] = (U.C[t].U^\dagger)_{aa}

   This performs very similarly to the solve_GEVP routine but is probably more efficient

   Seems to work fine for non-symmetric matrices too!
 */
#include "gens.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_version.h>

#include "stats.h"

struct GEVP_temps {
  gsl_matrix *a ;
  gsl_matrix *b ;
  gsl_eigen_nonsymmv_workspace *work ;
  gsl_vector_complex *eval ;
  gsl_matrix_complex *evec ;
  double *ev ;
  size_t N ;
} ;

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
  G -> work = gsl_eigen_nonsymmv_alloc( N ) ;
  G -> evec = gsl_matrix_complex_alloc( N , N ) ;
  G -> eval = gsl_vector_complex_alloc( N ) ;
  G -> ev = (double*)malloc( N*sizeof( double ) ) ;
  G -> N = N ;
}

static void
free_G( struct GEVP_temps *G )
{
  gsl_matrix_free( G -> a ) ;
  gsl_matrix_free( G -> b ) ;
  gsl_vector_complex_free( G -> eval ) ;
  gsl_matrix_complex_free( G -> evec ) ;
  gsl_eigen_nonsymmv_free( G -> work ) ;
  free( G -> ev ) ;
}

// simple matrix multiplication ( atomic left multiply ) a = b * a 
static void 
multab_atomic_left( const int NC ,
		    double a[ NC*NC ] , 
		    const double b[ NC*NC ] )
{
  // standard looped version
  size_t i , j , m ;
  double R[ NC ] ;
  register double sum ;
  for( i = 0 ; i < NC ; i++ ) { // loop cols
    for( j = 0 ; j < NC ; j ++ ) { // loop rows
      sum = 0.0 ;
      for( m = 0 ; m < NC ; m ++  ) { // loop elements in row or column
        sum += b[ m + j*NC ] * a[ i + m*NC ] ;
      }	
      R[j] = sum ;
    }
    // copy back over to a ...
    for( m = 0 ; m < NC ; m++ ) {
      a[ i + m*NC ] = R[m] ;
    }
  }
  return ;
}

// column elimination
static void
eliminate_column( const int NC ,
		  double *a , 
		  double *inverse , 
		  const size_t i , 
		  const size_t j )
{
  const double fac1 = a[ i + j*NC ] / a[ i*(NC+1) ] ;
  double *pA = a + i*NC ;
  size_t k ;
  // such a pattern elimintates cancelling zeros
  for( k = i + 1 ; k < NC ; k++ ) {
    a[ k + j*NC ] -= fac1 * pA[k] ;
  }
  pA = inverse + i*NC ;
  for( k = 0 ; k < NC ; k++ ) {
    // whatever we do to a, we do to the identity
    inverse[ k + j*NC ] -= fac1 * pA[k] ;
  }
  return ;
}

// row swapper
static void
swap_rows( const int NC ,
	   double *a , 
	   double *inverse , 
	   const size_t row_idx , 
	   const size_t piv )
{
  size_t l ;
  for( l = 0 ; l < NC ; l++ ) {
    register double tmp = a[l+row_idx*NC] ;
    a[l+row_idx*NC] = a[l+piv*NC] ;
    a[l+piv*NC] = tmp ;
    tmp = inverse[l+row_idx*NC] ;
    inverse[l+row_idx*NC] = inverse[l+piv*NC] ;
    inverse[l+piv*NC] = tmp ;
  }
  return ;
}

static int
gauss_jordan( const int NC ,
	      double M_1[ NC*NC ] , 
	      const double M[ NC*NC ] )
{
  double a[NC*NC] , ainv[NC*NC] ;
  double row_best , col_best ;
  size_t i, j , col_piv ;
  for( i = 0 ; i < NC*NC ; i++ ) {
    a[ i ] = M[ i ] ;
    ainv[ i ] = 0.0 ;
  }
  for( i = 0 ; i < NC ; i++ ) {
    ainv[ i*(NC+1) ] = 1.0 ;
  }
  for( i = 0 ; i < NC-1 ; i++ ) {
    const size_t piv_idx = i*(NC+1) ;
    const double pivot = a[piv_idx] * a[piv_idx] ;
    row_best = col_best = pivot ;
    col_piv = piv_idx ;
    size_t j ;
    for( j = i+1 ; j < NC ; j++ ) {
      const double col_attempt = a[i+j*NC] * a[i+j*NC] ;
      if( col_attempt > col_best ) {
	col_best = col_attempt ;
	col_piv = j ;
      }
    }
    if( col_best > row_best && col_best > pivot ) {
      swap_rows( NC , a , ainv , col_piv , i ) ;
    }
    // perform gaussian elimination to obtain the upper triangular
    for( j = NC-1 ; j > i ; j-- ) { // go up in other columns
      eliminate_column( NC , a , ainv , i , j ) ;
    }
  }
  // a is upper triangular, do the same for the upper half
  // no pivoting to be done here
  for( i = NC-1 ; i > 0 ; i-- ) {
    for( j = 0 ; j < i ; j++ ) {
      eliminate_column( NC , a , ainv , i , j ) ;
    }
  }
  // multiply each row by its M_1 diagonal
  double am1 ;
  for( j = 0 ; j < NC ; j++ ) {
    if( a[ j*( NC+1 ) ] == 0.0 ) { 
      fprintf( stderr , "[INVERSE] Matrix is singular!!\n" ) ; 
      return FAILURE ;
    }
    am1 = 1.0 / a[ j*( NC+1 ) ] ;
    for( i = 0 ; i < NC ; i++ ) {
      ainv[ i + j*NC ] = ainv[ i + j*NC ] * am1 ;
    }
  }
  // push the double complex precision matrix back into whatever precision
  for( i = 0 ; i < NC*NC ; i++ ) {
    M_1[i] = ainv[i] ;
  }
  return SUCCESS ;
}

// calculates the inverse of the matrix M outputs to M_1 //
static inline int 
inverse( const int NC , 
	 double M_1[ NC*NC ] , 
	 double M[ NC*NC ] )
{
  return gauss_jordan( NC , M_1 , M ) ;
}

// inverse square root
static int
denman_rootZ( int NC , double Z[NC*NC] )
{
  double M[ NC*NC ] , INVM[ NC*NC ] ;
  double tol = 1.0 ;
  size_t i ,iters = 0 ;
  for( i = 0 ; i < NC*NC ; i++ ) {
    M[i] = Z[i] ;
    Z[i] = ( i%(NC+1) ) ? 0.0 : 1.0 ;
  }
  // computes the inverse square root of A (Z), and the square root Y.
  while( tol > 1E-13 && iters < 25 ) {
    // invert M into INVM
    inverse( NC , INVM , M ) ;
    // add the identity
    for( i = 0 ; i < NC ; i++ ) {
      INVM[ i*(NC+1) ] += 1.0 ;
      M[ i*(NC+1) ]    += 1.0 ;
    }
    multab_atomic_left( NC , Z , INVM ) ;
    // Z <- 0.5 ( 1 + M^-1 ) Z
    // M <- 0.25 ( M + 2I + M^-1 )
    for( i = 0 ; i < NC*NC ; i++ ) {
      Z[i] *= 0.5 ; 
      M[i] = 0.25 * ( M[i] + INVM[i] ) ;
    }
    // make sure it is positive definite
    register double c ;
    tol = 0. ;
    for( i = 0 ; i < NC ; i++ ) {
      c = M[ i*(NC+1) ] ;
      tol += c*c ;
    }
    tol = fabs( tol - NC ) ;
    // should complain here
    if( iters == 24 ) {
      fprintf( stderr , "[TAYLOR_LOGS] Denman root Z ill convergence %e \n" , 
	       tol ) ;
      return FAILURE ;
    }
    iters++ ;
  }
  // that works!
  return SUCCESS ;
}

// compute the "optimised correlator" from the eigenvectors
static void
optimised_correlator( struct GEVP_temps *G ,
		      const size_t N ,
		      const size_t Ndata ,
		      double evals[Ndata][N] ,
		      const double C[Ndata][N*N] )
{
  // check diagonalisation with V^\dag C(t) V
  double complex V2[N][N][N];
  for( size_t a = 0 ; a < N ; a++ ) {
    for( size_t b = 0 ; b < N ; b++ ) {
      const double complex v1 = conj( gsl_matrix_complex_get( G -> evec , b , a ) ) ;
      for( size_t c = 0 ; c < N ; c++ ) {
	const double complex v2 = gsl_matrix_complex_get( G -> evec , c , a ) ;
	V2[a][b][c] = v1*v2 ;
      }
    }
  }
  for( size_t j = 0 ; j < Ndata ; j++ ) {    
    // diagonalise
    double complex *p1 = (double complex*)V2 ;
    for( size_t a = 0 ; a < N ; a++ ) {
      double *p2 = (double*)C[j] ;
      register double sum = 0.0 ;
      for( size_t i = 0 ; i < N*N ; i++ ) {
	sum += creal((*p1))*(*p2) ;
	p1++ ; p2++ ;
      }
      evals[j][a] = creal(sum) ;
    }
  }
  return ;
}

static void
set_C( const struct resampled *y ,
       const size_t N ,
       const size_t t0 ,
       const bool is_avg ,
       const size_t k ,
       const size_t Ndata ,
       double C[ N*N ] )
{
  // set c1 and root here
  size_t shift = 0 ;
  for( size_t i = 0 ; i < N*N ; i++ ) {
    if( is_avg == true ) {
      C[i] = y[ t0 + shift ].avg ;
    } else {
      C[i] = y[ t0 + shift ].resampled[k] ;
    }
    shift += Ndata ;
  }
}

// does the slow loopy matrix product C0 -> C1.C0.C1
static void
normalize_mat( const size_t N ,
	       const double C1[N*N] ,
	       double C0[N*N] )
{
  // overwrite C0 with C1.C0.C1 -> C1_ab C0_bc C1_cd
  double res[ N*N ] = {} ;
  for( int a = 0 ; a < N ; a++ ) {
    for( int b = 0 ; b < N ; b++ ) {
      for( int c = 0 ; c < N ; c++ ) {
	for( int d = 0 ; d < N ; d++ ) {
	  res[ d + N*a ] += \
	    C1[ b + N*a ] * C0[ c + N*b ] * C1[ d + N*c ] ;
	}
      }
    }
  }
  for( int a = 0 ; a < N ; a++ ) {
    for( int b = 0 ; b < N ; b++ ) {
      C0[b+a*N] = res[b+a*N] ;
    }
  }
}

static void
diagonalise( struct resampled *evalues ,
	     const size_t N ,
	     const size_t Ndata ,
	     const size_t td ,
	     const size_t k ,
	     const bool is_avg ,
	     const double C0[Ndata][N*N] )
{
  // ugh loop order is all weird
  struct GEVP_temps G ;
  allocate_G( &G , N );
  // diagonalise C once more for every j
  double evals[Ndata][ N ] = {} ;
  // set AB -> here we only use matrix a
  set_ab( G.a , G.b , C0[td] , C0[td] , N , N ) ;  
  // only need rotation evecs at some diagonalisation time
  gsl_eigen_nonsymmv( G.a, G.eval, G.evec, G.work) ;
  optimised_correlator( &G , N , Ndata , evals , C0 ) ;
  for( size_t j = 0 ; j < Ndata ; j++ ) {
    for( size_t i = 0 ; i < N ; i++ ) {
      if( is_avg == true ) {
	evalues[j+Ndata*i].avg = evals[j][i] ;
      } else {
	evalues[j+Ndata*i].resampled[k] = evals[j][i] ;
      }
    }
  }
  free_G( &G ) ;
}

// solve the GEVP of modified correlator matrix C(t_0)^-1/2 C(t) C(t_0)^{-1/2}
struct resampled *
solve_GEVP_rooted( const struct resampled *y ,
		   const size_t Ndata ,
		   const size_t N ,
		   const size_t M ,
		   const size_t t0 ,
		   const size_t td )
{
  fprintf( stdout , "GEVP solver ----> " ) ;
  
  if( N > M || N != M ) {
    fprintf( stderr , "[GEVP] cannot solve when N states "
	     "are less than M correlators or matrices aren't square\n" ) ;
    return NULL ;
  }

  // initialise the generalised eigenvalues
  struct resampled *evalues = malloc( Ndata*N*
				      sizeof( struct resampled ) ) ;
  for( size_t j = 0 ; j < Ndata*N ; j++ ) {
    evalues[j].resampled = malloc( y[0].NSAMPLES *
				   sizeof( double ) ) ;
    evalues[j].restype = y[0].restype ;
    evalues[j].NSAMPLES = y[0].NSAMPLES ;
  }
  
  // loop resamples
#pragma omp parallel for
  for( size_t k = 0 ; k < y[0].NSAMPLES ; k++ ) {
    double C0[Ndata][ N*M ] = {} , C1[ N*M ] = {} ;
    set_C( y , N , t0 , false , k , Ndata , C1 ) ;
    denman_rootZ( N , C1 ) ;
    for( size_t j = 0 ; j < Ndata ; j++ ) {
      set_C( y , N , j , false , k , Ndata , C0[j] ) ;
      normalize_mat( N , C1 , C0[j] ) ;
    }
    diagonalise( evalues , N , Ndata , td , k , false , C0 ) ;
  }

  {
    // do average
    double C0[Ndata][ N*M ] = {} , C1[ N*M ] = {} ;
    set_C( y , N , t0 , true , 0 , Ndata , C1 ) ;
    denman_rootZ( N , C1 ) ;
    for( size_t j = 0 ; j < Ndata ; j++ ) {
      set_C( y , N , j , true , 0 , Ndata , C0[j] ) ;
      normalize_mat( N , C1 , C0[j] ) ;
    }
    diagonalise( evalues , N , Ndata , td , 0 , true , C0 ) ;
  }
  
  return evalues ;
}
