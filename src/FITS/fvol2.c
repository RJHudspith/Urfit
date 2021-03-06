/**
   @file fvol1.c
   @brief finite volume fit y = A + B*X + e^{-\sqrt{X}*L }
 */
#include "gens.h"

#include "Nder.h"

//#define STRANGE

//#define INA

//#define MPISYST
//#define VSYST
//#define ASYST

//#define POLE
//#define XSQ

#ifdef STRANGE
int
fvol2_NMAX( void )
{
  return 11 ;
}

static const double MPIL[11] = { 5.312 , 4.029 ,
				 4.354 , 5.818 , 3.917 , 4.642 ,
				 5.146 ,
				 4.358 , 6.408 ,
				 5.109 ,
				 1E8 } ;

static const double a2[11] = { 0.25319 , 0.25319 ,
			       0.191536 , 0.191536 ,  0.191536 , 0.191536 ,
			       0.14967 ,
			       0.10605 , 0.10605 ,
			       0.06372 ,
			       0.0 } ;

#else

#ifdef ASYST

int
fvol2_NMAX( void )
{
  return 13 ;
}
// U103 , H101 , U102 , H105 , C101
// B450
// H200 , N202 , N200 , D200
// N300
static const double MPIL[13] = { 4.354 , 5.818 , 3.743 , 3.917 , 4.642 ,
				 5.146 , 5.377 ,
				 4.358 , 6.408 , 4.416 , 4.154 ,
				 5.109 ,
				 1E8 } ;

static const double a2[13] = { 0.191536 , 0.191536 , 0.191536 , 0.191536 , 0.191536 ,
			       0.14967 , 0.14967 ,
			       0.10605 , 0.10605 , 0.10605 , 0.10605 ,
			       0.06372 ,
			       0.0 } ;
#elif (defined VSYST)

int
fvol2_NMAX( void )
{
  return 13 ;
}

// A653 , A654
// U103 , H101 , C101
// B450
// H200 , N202 , N200 , D200
// N300
static const double MPIL[13] = { 5.312 , 4.029 ,
				 4.354 , 5.818 , 4.642 ,
				 5.146 , 5.377 ,
				 4.358 , 6.408 , 4.416 , 4.154 ,
				 5.109 ,
				 1E8 } ;

static const double a2[13] = { 0.25319 , 0.25319 ,
			       0.191536 , 0.191536 , 0.191536 ,
			       0.14967 , 0.14967 ,
			       0.10605 , 0.10605 , 0.10605 , 0.10605 ,
			       0.06372 ,
			       0.0 } ;

#elif (defined MPISYST)

int
fvol2_NMAX( void )
{
  return 8 ;
}

// A654
// U102 , H105 , C101
// D450
// N200 , D200
static const double MPIL[8] = { 4.029 ,
				3.743 , 3.917 , 4.642 ,
				5.377 ,
				4.416 , 4.154 ,
				1E8 } ;

static const double a2[8] = { 0.25319 ,
			      0.191536 , 0.191536 , 0.191536 ,
			      0.14967 ,
			      0.10605 , 0.10605 ,
			      0.0 } ;
#else

int
fvol2_NMAX( void )
{
  return 16 ;
}

// A653 , A654
// U103 , H101 , U102 , H105 , C101
// B450 , D450
// H200 , N202 , N200 , D200
// N300
static const double MPIL[16] = { 5.312 , 4.029 ,
				 4.354 , 5.818 , 3.743 , 3.917 , 4.642 ,
				 5.146 , 5.377 ,
				 4.358 , 6.408 , 4.416 , 4.154 ,
				 5.109 ,
				 1E8 ,
				 1E8 } ;

static const double a2[16] = { 0.25319 , 0.25319 ,
			       0.191536 , 0.191536 , 0.191536 , 0.191536 , 0.191536 ,
			       0.14967 , 0.14967 ,
			       0.10605 , 0.10605 , 0.10605 , 0.10605 ,
			       0.06372 ,
			       0.191536 ,
			       0.0 } ;
#endif

#endif 

static double
map_a( const int Npars ,
       const double *fparams )
{
  switch( Npars) {
  case 0 : case 1 : return fparams[4] ;
  case 2 : case 3 : case 4 : case 5 : case 6 : return fparams[5] ;;
  case 7 : case 8 : return fparams[6] ;
  case 9 : case 10 : case 11 : case 12 : return fparams[7] ;
  case 13 : return fparams[8] ;
  case 14 : return fparams[5] ;
  case 15 : return 0 ;
  }
}

double
ffvol2( const struct x_desc X , const double *fparams , const size_t Npars )
{
  const double asq = map_a( Npars , fparams ) ;
  return fparams[0] + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars]/2 ) + fparams[3]*asq ;
  /*
#ifdef INA
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars]/2 ) + fparams[3]*sqrt(a2[Npars]) ) ;
#elif (defined POLE)
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[4]*log(X.X) + fparams[2]*exp( -MPIL[Npars]/2 ) + fparams[3]*sqrt(a2[Npars]) ) ;
#elif (defined XSQ)
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[4]*X.X*X.X + fparams[2]*exp( -MPIL[Npars]/2 ) + fparams[3]*a2[Npars] ) ;
#else // a^2
  return fparams[0] * ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[Npars]/2 ) + fparams[3]*a2[Npars] ) ;
#endif
  */
}

void
fvol2_f( double *f , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ; 
  for( i = 0 ; i < DATA -> n ; i++ ) {
    double p[ DATA -> Npars ] ;
    for( j = 0 ; j < DATA -> Npars ; j++ ) {
      p[ j ] = fparams[ DATA -> map[ i ].p[ j ] ] ;
    }
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;
    f[i] = ffvol2( X , p , i ) - DATA -> y[i] ;
  }
  return ;
}

// derivatives
void
fvol2_df( double **df , const void *data , const double *fparams )
{
  const struct data *DATA = (const struct data*)data ;
  size_t i , j ;
  for( i = 0 ; i < DATA -> n ; i++ ) {
    struct x_desc X = { DATA -> x[i] , DATA -> LT[i] ,
			DATA -> N , DATA -> M } ;

    const double asq = map_a( i , fparams ) ;
    df[0][i] = 1 ;
    df[1][i] = X.X ;
    df[2][i] = exp( -MPIL[i]/2 ) ;
    df[3][i] = asq ;

    // hmmmm
    df[4][i] = df[5][i] = df[6][i] = df[7][i] = df[8][i] = 0.0 ;

    switch( i ) {
    case 0 : case 1 :
      df[4][i] = fparams[3] ;
      break ;
    case 2 : case 3 : case 4 : case 5 : case 6 :
      df[5][i] = fparams[3] ;
      break ;
    case 7 : case 8 :
      df[6][i] = fparams[3] ;
      break ;
    case 9 : case 10 : case 11 : case 12 :
      df[7][i] = fparams[3] ;
      break ;
    case 13 :
      df[8][i] = fparams[3] ;
      break ;
    }
      
    /*
    #ifdef INA
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[i]/2 ) + fparams[3]*sqrt(a2[i]) ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i]/2 ) ;
    df[3][i] = fparams[0] * sqrt(a2[i]) ;
    #elif (defined POLE)
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[4]*log(X.X) + fparams[2]*exp( -MPIL[i]/2 ) + fparams[3]*sqrt(a2[i]) ) ;
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i]/2 ) ;
    df[3][i] = fparams[0] * sqrt(a2[i]) ;
    df[4][i] = fparams[0] * log( X.X ) ;
    #elif (defined XSQ)
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[4]*X.X*X.X + fparams[2]*exp( -MPIL[i]/2 ) + fparams[3]*a2[i] ) ; 
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i]/2 ) ;
    df[3][i] = fparams[0] * a2[i] ;
    df[4][i] = fparams[0] * X.X * X.X ;
    #else
    df[0][i] = ( 1 + fparams[1]*X.X + fparams[2]*exp( -MPIL[i]/2 ) + fparams[3]*a2[i] ) ;
    df[1][i] = fparams[0] * X.X ;
    df[2][i] = fparams[0] * exp( -MPIL[i]/2 ) ;
    df[3][i] = fparams[0] * a2[i] ;
    #endif
    */
  }
  return ;
}

// second derivatives? Will we ever use them - J?
void
fvol2_d2f( double **d2f , const void *data , const double *fparams )
{
  return ;
}

void
fvol2_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit )
{
  fparams[0] = 80 ;
  fparams[1] = -0.1 ;
  fparams[2] = 1 ;
  fparams[3] = -0.1 ;
  return ;
}
