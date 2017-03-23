/**
   @file gen_ders.c
   @brief general N-th order finite difference code
 */
#include "gens.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>

static const double fac[ 101 ] = {1.000000000000000e+00 ,
				  1.000000000000000e+00 , 2.000000000000000e+00 , 
				  6.000000000000000e+00 , 2.400000000000000e+01 , 
				  1.200000000000000e+02 , 7.200000000000000e+02 , 
				  5.040000000000000e+03 , 4.032000000000000e+04 , 
				  3.628800000000000e+05 , 3.628800000000000e+06 , 
				  3.991680000000000e+07 , 4.790016000000000e+08 , 
				  6.227020800000000e+09 , 8.717829120000000e+10 , 
				  1.307674368000000e+12 , 2.092278988800000e+13 , 
				  3.556874280960000e+14 , 6.402373705728000e+15 , 
				  1.216451004088320e+17 , 2.432902008176640e+18 , 
				  5.109094217170944e+19 , 1.124000727777608e+21 , 
				  2.585201673888498e+22 , 6.204484017332394e+23 , 
				  1.551121004333099e+25 , 4.032914611266057e+26 , 
				  1.088886945041835e+28 , 3.048883446117138e+29 , 
				  8.841761993739701e+30 , 2.652528598121910e+32 , 
				  8.222838654177922e+33 , 2.631308369336935e+35 , 
				  8.683317618811886e+36 , 2.952327990396041e+38 , 
				  1.033314796638614e+40 , 3.719933267899012e+41 , 
				  1.376375309122634e+43 , 5.230226174666010e+44 , 
				  2.039788208119744e+46 , 8.159152832478977e+47 , 
				  3.345252661316380e+49 , 1.405006117752880e+51 , 
				  6.041526306337383e+52 , 2.658271574788449e+54 , 
				  1.196222208654802e+56 , 5.502622159812088e+57 , 
				  2.586232415111682e+59 , 1.241391559253607e+61 , 
				  6.082818640342675e+62 , 3.041409320171338e+64 , 
				  1.551118753287382e+66 , 8.065817517094388e+67 , 
				  4.274883284060025e+69 , 2.308436973392414e+71 , 
				  1.269640335365828e+73 , 7.109985878048635e+74 ,
				  4.052691950487722e+76 , 2.350561331282879e+78 , 
				  1.386831185456899e+80 , 8.320987112741392e+81 , 
				  5.075802138772248e+83 , 3.146997326038794e+85 , 
				  1.982608315404440e+87 , 1.268869321858842e+89 , 
				  8.247650592082472e+90 , 5.443449390774431e+92 , 
				  3.647111091818868e+94 , 2.480035542436831e+96 , 
				  1.711224524281413e+98 , 1.197857166996989e+100 , 
				  8.504785885678622e+101 , 6.123445837688608e+103 , 
				  4.470115461512683e+105 , 3.307885441519386e+107 , 
				  2.480914081139539e+109 , 1.885494701666050e+111 , 
				  1.451830920282858e+113 , 1.132428117820629e+115 , 
				  8.946182130782973e+116 , 7.156945704626378e+118 , 
				  5.797126020747366e+120 , 4.753643337012840e+122 , 
				  3.945523969720657e+124 , 3.314240134565352e+126 , 
				  2.817104114380549e+128 , 2.422709538367272e+130 , 
				  2.107757298379527e+132 , 1.854826422573984e+134 , 
				  1.650795516090845e+136 , 1.485715964481761e+138 , 
				  1.352001527678402e+140 , 1.243841405464130e+142 , 
				  1.156772507081641e+144 , 1.087366156656742e+146 , 
				  1.032997848823905e+148 , 9.916779348709491e+149 , 
				  9.619275968248206e+151 , 9.426890448883242e+153 , 
				  9.332621544394410e+155 } ;


// selection sort for the shifts
static void
set_shifts( double *shifts ,
	    double *ydata , 
	    const double *data ,
	    const double *x ,
	    const size_t order ,
	    const size_t Ndata ,
	    const size_t idx )
{
  // precompute shifts
  size_t i , j ;
  for( i = 0 ; i < Ndata ; i++ ) {
    shifts[ i ] = x[i] - x[idx] ;
    ydata[ i ] = data[i] ;
  }
  // selection sort
  for( i = 0 ; i < order+1 ; i++ ) {
    size_t min = i ;
    for( j = i+1 ; j < Ndata ; j++ ) {
      if( fabs( shifts[j] ) < fabs( shifts[min] ) ) {
	min = j ;
      }
    }
    //
    if( min != i ) {
      double tmp = shifts[i] ;
      shifts[i] = shifts[min] ;
      shifts[min] = tmp ;
      tmp = ydata[i] ;
      ydata[i] = ydata[min] ;
      ydata[min] = tmp ;
    }
    //
  }
  return ;
}

double **
get_ders( const double *y ,
	  const double *x ,
	  const size_t N ,
	  const size_t order )
{
  if( order > N ) {
    return NULL ;
  }
  
  double **ders = malloc( N * sizeof( double * ) ) ;
  
  double *shifts = malloc( N * sizeof( double ) ) ;
  double *yshift = malloc( N * sizeof( double ) ) ;

  gsl_matrix *alpha     = gsl_matrix_alloc( order+1 , order+1 ) ;
  gsl_vector *beta      = gsl_vector_alloc( order+1 ) ;
  gsl_vector *delta     = gsl_vector_alloc( order+1 ) ;
  gsl_permutation *perm = gsl_permutation_alloc( order+1 ) ;

  size_t i , j , k ;
			  
  for( i = 0 ; i < N ; i++ ) {
    
    // allocate derivatives
    ders[i] = malloc( ( order+1 ) * sizeof( double ) ) ;
      
    // set the shifts
    set_shifts( shifts , yshift , y , x , order , N , i ) ;

    // the rest get filled from the index "lo"
    for( k = 0 ; k < order+1 ; k++ ) {
      for( j = 0 ; j < order+1 ; j++ ) {
	gsl_matrix_set( alpha , k , j , pow( shifts[k] , j ) / fac[j] ) ;
      }
      gsl_vector_set( beta , k , yshift[k] ) ;
    }

    // write out the matrix we are solving
    #ifdef VERBOSE
    for( k = 0 ; k < order+1 ; k++ ) {
      for( j = 0 ; j < order+1 ; j++ ) {
	printf( " %f " , gsl_matrix_get( alpha , k , j ) ) ;
      }
      printf( "\n" ) ;
    }
    printf( "\n" ) ;
    #endif

    int signum , Flag = SUCCESS ;
    if( gsl_linalg_LU_decomp( alpha , perm , &signum ) != GSL_SUCCESS ) {
      fprintf( stderr , "[FINITE] LU decomp failed\n" ) ;
      Flag = FAILURE ;
    }
    if( gsl_linalg_LU_solve( alpha , perm , beta , delta ) != GSL_SUCCESS ) {
      fprintf( stderr , "[FINITE] LU solve failure \n" ) ;
      Flag = FAILURE ;
    }

    for( j = 0 ; j < order+1 ; j++ ) {
      ders[i][j] = gsl_vector_get( delta , j ) ;
    }
  }

  // free the temporary shifts and data
  free( shifts ) ; free( yshift ) ;

  // free the gsl vectors
  gsl_vector_free( beta ) ;
  gsl_matrix_free( alpha ) ;
  gsl_permutation_free( perm ) ;
  gsl_vector_free( delta ) ;

  return ders ;
}
