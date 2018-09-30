/**
   @file pade_laplace.c
   @brief laplace-pade method
 */
#include "gens.h"

#include <gsl/gsl_poly.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "pade_coefficients.h"

//#define VERBOSE

// LUT for the factorial
static const double fac[169] = { 1.000000e+00 , 1.000000e+00 , 2.000000e+00 , 6.000000e+00 , 2.400000e+01 , 1.200000e+02 , 7.200000e+02 , 5.040000e+03 , 4.032000e+04 , 
3.628800e+05 , 3.628800e+06 , 3.991680e+07 , 4.790016e+08 , 6.227021e+09 , 8.717829e+10 , 1.307674e+12 , 2.092279e+13 , 
3.556874e+14 , 6.402374e+15 , 1.216451e+17 , 2.432902e+18 , 5.109094e+19 , 1.124001e+21 , 2.585202e+22 , 6.204484e+23 , 
1.551121e+25 , 4.032915e+26 , 1.088887e+28 , 3.048883e+29 , 8.841762e+30 , 2.652529e+32 , 8.222839e+33 , 2.631308e+35 , 
8.683318e+36 , 2.952328e+38 , 1.033315e+40 , 3.719933e+41 , 1.376375e+43 , 5.230226e+44 , 2.039788e+46 , 8.159153e+47 , 
3.345253e+49 , 1.405006e+51 , 6.041526e+52 , 2.658272e+54 , 1.196222e+56 , 5.502622e+57 , 2.586232e+59 , 1.241392e+61 , 
6.082819e+62 , 3.041409e+64 , 1.551119e+66 , 8.065818e+67 , 4.274883e+69 , 2.308437e+71 , 1.269640e+73 , 7.109986e+74 , 
4.052692e+76 , 2.350561e+78 , 1.386831e+80 , 8.320987e+81 , 5.075802e+83 , 3.146997e+85 , 1.982608e+87 , 1.268869e+89 , 
8.247651e+90 , 5.443449e+92 , 3.647111e+94 , 2.480036e+96 , 1.711225e+98 , 1.197857e+100 , 8.504786e+101 , 6.123446e+103 , 
4.470115e+105 , 3.307885e+107 , 2.480914e+109 , 1.885495e+111 , 1.451831e+113 , 1.132428e+115 , 8.946182e+116 , 7.156946e+118 , 
5.797126e+120 , 4.753643e+122 , 3.945524e+124 , 3.314240e+126 , 2.817104e+128 , 2.422710e+130 , 2.107757e+132 , 1.854826e+134 , 
1.650796e+136 , 1.485716e+138 , 1.352002e+140 , 1.243841e+142 , 1.156773e+144 , 1.087366e+146 , 1.032998e+148 , 9.916779e+149 , 
9.619276e+151 , 9.426890e+153 , 9.332622e+155 , 9.332622e+157 , 9.425948e+159 , 9.614467e+161 , 9.902901e+163 , 1.029902e+166 , 
1.081397e+168 , 1.146281e+170 , 1.226520e+172 , 1.324642e+174 , 1.443860e+176 , 1.588246e+178 , 1.762953e+180 , 1.974507e+182 , 
2.231193e+184 , 2.543560e+186 , 2.925094e+188 , 3.393109e+190 , 3.969937e+192 , 4.684526e+194 , 5.574586e+196 , 6.689503e+198 , 
8.094299e+200 , 9.875044e+202 , 1.214630e+205 , 1.506142e+207 , 1.882677e+209 , 2.372173e+211 , 3.012660e+213 , 3.856205e+215 , 
4.974504e+217 , 6.466855e+219 , 8.471581e+221 , 1.118249e+224 , 1.487271e+226 , 1.992943e+228 , 2.690473e+230 , 3.659043e+232 , 
5.012889e+234 , 6.917786e+236 , 9.615723e+238 , 1.346201e+241 , 1.898144e+243 , 2.695364e+245 , 3.854371e+247 , 5.550294e+249 , 
8.047926e+251 , 1.174997e+254 , 1.727246e+256 , 2.556324e+258 , 3.808923e+260 , 5.713384e+262 , 8.627210e+264 , 1.311336e+267 , 
2.006344e+269 , 3.089770e+271 , 4.789143e+273 , 7.471063e+275 , 1.172957e+278 , 1.853272e+280 , 2.946702e+282 , 4.714724e+284 , 
7.590705e+286 , 1.229694e+289 , 2.004402e+291 , 3.287219e+293 , 5.423911e+295 , 9.003692e+297 , 1.503617e+300 , 2.526076e+302 };

// compute the coefficients
static void
dLdp( double *d ,
      const size_t Nders ,
      const double *x ,
      const double *y ,
      const size_t N ,
      const double p0 )
{
  double epx[ N ] ;
  size_t i , j ;
  // initialise the exponential factors and the powers of t
  for( i = 0 ; i < N ; i++ ) {
    epx[ i ] = exp( -p0 * x[i] ) * y[i] ;
  }
  for( i = 0 ; i < Nders ; i++ ) {
    d[i] = 0.0 ;
    for( j = 0 ; j < N-1 ; j++ ) {
      d[i] += ( epx[j+1] + epx[j] ) * ( x[j+1] - x[j] ) ;
      epx[j] *= -x[j] ;
    }
    epx[j] *= -x[j] ;
    d[i] /= ( 2. * fac[i] ) ;
    #ifdef VERBOSE
    printf( "D[ %zu ] %e %e \n" , i , d[i] , fac[i] ) ;
    #endif
  }
  return ;
}

// compute the coefficients using simpson's rule
static void
dLdp_simp( double *d ,
	   const size_t Nders ,
	   const double *x ,
	   const double *y ,
	   const size_t N ,
	   const double p0 )
{
  double epx[ N ] ;
  size_t i , j ;
  // initialise the exponential factors and the powers of t
  for( i = 0 ; i < N ; i++ ) {
    epx[ i ] = exp( -p0 * x[i] ) * y[i] ;
  }
  for( i = 0 ; i < Nders ; i++ ) {
    // perform the integration
    d[i] = 0.0 ;
    for( j = 1 ; j < N/2 ; j++ ) {
      d[i] += ( epx[ 2*j - 2 ] + 4*epx[ 2*j-1 ] + epx[ 2*j ] ) ;
    }
    // update the t-factors
    for( j = 0 ; j < N ; j++ ) {
      epx[j] *= -x[j] ;
    }
    d[i] /= ( 3. * fac[i] ) ;
    #ifdef VERBOSE
    printf( "D[ %zu ] %e %e \n" , i , d[i] , fac[i] ) ;
    #endif
  }
  return ;
}

// qsort comparison
static int 
comp( const void *elem1 , 
      const void *elem2 ) 
{
  const double f = *( (double*)elem1 ) ;
  const double s = *( (double*)elem2 ) ;
  if (f < s) return  1 ;
  if (f > s) return -1 ;
  return 0 ;
}

typedef enum { GOOD , BAD , UGLY } GBU ;

static GBU
filter( const double complex z )
{
  if( creal( z ) > 0 ) {
    return BAD ;
  } else if( fabs( cimag( z ) ) < 1E-12 ) {
    return GOOD ;
  } else {
    return UGLY ;
  }
}

static size_t
stable_pade( double *masses ,
	     const double *d ,
	     const size_t Nexp ,
	     const double p0 )
{
  size_t i ;
  
  const size_t N = 4*Nexp ;
  double *pade = malloc( 2*N * sizeof( double ) ) ;
  
  pades_from_poly( pade , d , N , N ) ;
#ifdef VERBOSE
  for( i = 0 ; i < 2*N ; i++ ) {
    printf( "[PLAP] pade_%zu %e \n" , i , pade[i] ) ;
  }
#endif

  // use a root-finding algorithm to determine the poles
  double *denom = malloc( ( N + 1 ) * sizeof( double ) ) ;
  denom[ 0 ] = 1 ;
  for( i = 0 ; i < N  ; i++ ) {
    denom[ i + 1 ] = pade[ N + i ] ;
  }

  gsl_poly_complex_workspace *w =		\
    gsl_poly_complex_workspace_alloc( N+1 ) ;
  double complex *z = malloc( N * sizeof( double complex ) ) ;
  
  gsl_poly_complex_solve( denom , N+1 , w , (double*)z ) ;

  // check usual possible solutions
  size_t Ngood = 0 , Nbad = 0 , Nugly = 0 ;
  for( i = 0 ; i < N ; i++ ) {
    z[i] += p0 ;
    #ifdef VERBOSE
    fprintf( stdout , "[PLAP] Z_%zu %f %f \n" ,
	     i , creal( z[i] ) , cimag( z[i] ) ) ;
    #endif
    switch( filter( z[i] ) ) {
    case GOOD : Ngood++ ; break ;
    case BAD : Nbad++ ; break ;
    case UGLY : Nugly++ ; break ;
    }
  }
  #ifdef VERBOSE
  fprintf( stdout , "[PLAP] Ngood %zu, Nbad %zu , Nugly %zu \n" ,
	   Ngood , Nbad , Nugly ) ;
  #endif

  size_t gidx = 0 , bidx = 0 , uidx = 0 ;
  double Good[ Ngood ] , Bad[ Nbad ] , Ugly[ Nugly ] ;
  for( i = 0 ; i < N ; i++ ) {
    switch( filter( z[i] ) ) {
    case GOOD : Good[ gidx++ ] = creal( z[i] ) ; break ;
    case BAD : Bad[ bidx++ ] = creal( z[i] ) ; break ;
    case UGLY : Ugly[ uidx++ ] = creal( z[i] ) ; break ;
    }
  }

  qsort( Good , Ngood , sizeof(double) , comp ) ;
  qsort( Ugly , Nugly , sizeof(double) , comp ) ;
  qsort( Bad , Nbad , sizeof(double) , comp ) ;

  size_t idx = 0 ;
  for( i = 0 ; i < Ngood ; i++ ) {
    if( idx >= Nexp ) goto end ;
    masses[ idx++ ] = -Good[i] ;
  }

  for( i = 0 ; i < Nugly ; i++ ) {
    if( idx >= Nexp ) goto end ;
    masses[ idx++ ] = -Ugly[i] ;
  }

  for( i = 0 ; i < Nbad ; i++ ) {
    if( idx >= Nexp ) goto end ;
    masses[ idx++ ] = -Bad[i] ;
  }


end :

  #ifdef VERBOSE
  for( i = 0 ; i < Nexp ; i++ ) {
    printf( "Mass_%zu %e \n" , i , masses[i] ) ;
  }
  #endif
  
  gsl_poly_complex_workspace_free( w ) ;
  free( denom ) ;
  free( pade ) ;
  free( z ) ;

  return Ngood ;
}

// get amplitudes by SVD
static int
get_amplitudes( double *Amps ,
		const double *x ,
		const double *y ,
		const size_t Ndata ,
		const double *masses ,
		const double Nexps )
{
  gsl_matrix *alpha = gsl_matrix_alloc( Ndata-1 , Nexps ) ;
  gsl_matrix *Q     = gsl_matrix_alloc( Nexps , Nexps ) ;
  gsl_vector *S     = gsl_vector_alloc( Nexps ) ;
  gsl_vector *Work  = gsl_vector_alloc( Nexps ) ;
  
  gsl_vector *beta  = gsl_vector_alloc( Ndata-1 ) ;
  gsl_vector *delta = gsl_vector_alloc( Nexps ) ;

  size_t i , j ;
  for( i = 1 ; i < Ndata ; i++ ) {
    // set alpha
    for( j = 0 ; j < Nexps ; j++ ) {
      gsl_matrix_set( alpha , i-1 , j , exp( -masses[j] * x[i] ) ) ;
    }
    gsl_vector_set( beta , i-1 , y[i] ) ;
  }

  if( gsl_linalg_SV_decomp( alpha , Q , S , Work ) != GSL_SUCCESS ) {
    return FAILURE ;
  }
  if( gsl_linalg_SV_solve( alpha , Q , S , beta , delta ) != GSL_SUCCESS ) {
    return FAILURE ;
  }

  // tell us the result
  for( i = 0 ; i < Nexps ; i++ ) {
    #ifdef VERBOSE
    fprintf( stdout , "%f e ^ { %f * x } \n" ,
	     gsl_vector_get( delta , i ) , masses[i] ) ;
    #endif
    Amps[i] = gsl_vector_get( delta , i ) ;
  }

  // free the gsl vectors
  gsl_matrix_free( alpha ) ;
  gsl_matrix_free( Q ) ;
  gsl_vector_free( S ) ;
  gsl_vector_free( Work ) ;
  
  gsl_vector_free( beta ) ;
  gsl_vector_free( delta ) ;

  return SUCCESS ;
}

int
pade_laplace( double *fparams ,
	      const double *x ,
	      const double *y ,
	      const size_t Ndata ,
	      const size_t Nexps ,
	      const double p0 )
{
  const size_t Nders = 4 * ( (Nexps + 8) ) + 1 ;
  if( Nders > 169 ) {
    fprintf( stderr , "[PLAP] too many exps for our simple measurement\n" ) ;
    return FAILURE ;
  }
  
  double *d = malloc( Nders * sizeof( double ) ) ;
  double Masses[ Nexps ] , Amps[ Nexps ] ;
  size_t i ;

  // compute the derivatives of the amplitudes
  //dLdp( d , Nders , x , y , Ndata , p0 ) ;
  dLdp_simp( d , Nders , x , y , Ndata , p0 ) ;
  
  // get the poles of the pade representation
  int Flag = SUCCESS ;
  size_t Trials = Nexps , Nbest = 0 ;
  for( Trials = Nexps ; Trials < Nexps + 8 ; Trials++ ) {
    double Mtrials[ Trials ] ;
    const size_t Ngood = stable_pade( Mtrials , d , Trials , p0 ) ;
    if( Ngood >= Nexps ) {
      size_t j ;
      for( j = 0 ; j < Nexps ; j++ ) {
	Masses[ j ] = Mtrials[ j ] ;
      }
      break ;
    } else if( Ngood > Nbest ) {
      Nbest = Ngood ;
      size_t j ;
      for( j = 0 ; j < Nexps ; j++ ) {
	Masses[ j ] = Mtrials[ j ] ;
      }
    }
  }
  if( Flag == FAILURE ) goto free ;
  
  // compute the amplitudes from the masses
  Flag = get_amplitudes( Amps , x , y , Ndata , Masses , Nexps ) ;
  if( Flag == FAILURE ) goto free ;
  
  // set the fit parameters
  for( i = 0 ; i < Nexps ; i++ ) {
    fparams[ 2*i + 0 ] = Amps[ i ] ;
    fparams[ 2*i + 1 ] = Masses[ i ] ;
  }

 free :
  free( d ) ;

  return Flag ;
}
