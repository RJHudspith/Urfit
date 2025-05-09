/**
   @file read_inputs.c
   @brief read an input file
 */
#include "gens.h"

#include <string.h>

#include "init.h"      // free_fit && free_data
#include "read_fit.h"
#include "read_graph.h"
#include "read_inputs.h"
#include "read_stats.h"
#include "read_traj.h"

// set these to something reasonable
#define STR1_LENGTH (64)
#define STR2_LENGTH (256)

#define VERBOSE

// get the file length
static size_t
get_file_length( FILE *Infile )
{
  size_t Ntags = 0 ;
  char str1[ STR1_LENGTH ] , str2[ STR2_LENGTH ] ;
  while( true ) {
    const int test = fscanf( Infile , "%s = %s" , str1 , str2 ) ;
    if( test == EOF ) break ; // if the file is finished we leave
    if( test != 2 ) continue ; // if there are comments and stuff
    Ntags ++ ;
  }
  rewind( Infile ) ;
  return Ntags ;
}

// allocate the input file struct
static struct flat_file *
pack_inputs( FILE *Infile ,
	     const size_t Ntags )
{
  struct flat_file *Flat = malloc( Ntags * sizeof( struct flat_file ) ) ;
  char str1[ STR1_LENGTH ] , str2[ STR2_LENGTH ] ;
  size_t i = 0 ;
  while( true ) {
    
    const int test = fscanf( Infile , "%s = %s" , str1 , str2 ) ;
    
    if( test == EOF ) break ; // if the file is finished we leave
    if( test != 2 ) continue ; // if there are comments and stuff
    
    Flat[i].Token_Length = 64 ; //strlen( str1 ) ;
    Flat[i].Token = malloc( Flat[i].Token_Length * sizeof( char ) ) ;
    
    sprintf( Flat[i].Token , "%s" , str1 ) ;
    
    Flat[i].Value_Length = 256 ; //strlen( str2 ) ;
    Flat[i].Value = malloc( Flat[i].Value_Length * sizeof( char ) ) ;
    
    sprintf( Flat[i].Value , "%s" , str2 ) ;
    
    i++ ;
  }
  return Flat ;
}

// frees up the struct
static void
unpack_inputs( struct flat_file *Flat ,
	       const size_t Ntags )
{
  size_t i ;
  for( i = 0 ; i < Ntags ; i++ ) {
    free( Flat[i].Token ) ;
    free( Flat[i].Value ) ;
  }
  free( Flat ) ;
  return ;
}

// test if two strings are equivalent
int
are_equal( const char *str_1 , const char *str_2 ) 
{
  return ( strcmp( str_1 , str_2 ) != 0 ) ? 0 : 1 ;
}

// free the input file data
void
free_inputs( struct input_params *Input )
{ 
  size_t i ;
  if( Input -> Traj != NULL ) {
    for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
      free( Input -> Traj[i].Dimensions ) ;
      if( Input -> Traj[i].FileX != NULL ) {
	free( Input -> Traj[i].FileX ) ;
      }
      if( Input -> Traj[i].FileY ) {
	free( Input -> Traj[i].FileY ) ;
      }
    }
    free( Input -> Traj ) ;
  }
  free( Input -> Graph.Name ) ;
  free( Input -> Graph.Xaxis ) ;
  free( Input -> Graph.Yaxis ) ;
  return ;
}

// read the input file
int
read_inputs( struct input_params *Input ,
	     const char *filename )
{
  // open the file
  FILE *infile = fopen( filename , "r" ) ;
  int Flag = SUCCESS ;

  // initialise any arrays of Input to NULL
  Input -> Data.x = NULL ;
  Input -> Data.y = NULL ;
  Input -> Data.Ndata = NULL ;
  Input -> Data.Cov.W = NULL ;    
  Input -> Data.Ntot = 0 ;
  Input -> Data.Nsim = 0 ;
  Input -> Data.LT = NULL ;

  Input -> Traj = NULL ;
  
  Input -> Fit.Sims = NULL ;
  Input -> Fit.Prior = NULL ;
  Input -> Fit.map = NULL ;
  Input -> Fit.Guess = NULL ;

  Input -> Graph.Name = NULL ;
  Input -> Graph.Xaxis = NULL ;
  Input -> Graph.Yaxis = NULL ;
  
  if( infile == NULL ) {
    printf( "[INPUT] Cannot find input file %s \n" , filename ) ;
    return FAILURE ;
  }

  // get the flat file length
  const size_t Ntags = get_file_length( infile ) ;

  fprintf( stdout , "[INPUTS] Ntags = %zu \n" , Ntags ) ;

  // pack the input file
  struct flat_file *Flat = pack_inputs( infile , Ntags ) ;

  fprintf( stdout , "[INPUTS] inputs packed\n" ) ;

  // get the filetype tag
  size_t io_tag = 0 ;
  if( ( io_tag = tag_search( Flat , "FileType" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUT] FileType tag not found in input file \n" ) ;
    Flag = FAILURE ;
  } else {
    if( are_equal( Flat[ io_tag ].Value , "Corr_File" ) ) {
      Input -> FileType = Corr_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "Distribution_File" ) ) {
      Input -> FileType = Distribution_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "Fake_File" ) ) {
      Input -> FileType = Fake_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "Flat_File" ) ) {
      Input -> FileType = Flat_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "GLU_Tcorr_File" ) ) {
      Input -> FileType = GLU_Tcorr_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "GLU_Qmoment_File" ) ) {
      Input -> FileType = GLU_Qmoment_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "GLU_File" ) ) {
      Input -> FileType = GLU_File ;
    } else if( are_equal( Flat[ io_tag ].Value , "Adler_File" ) ) {
      Input -> FileType = Adler_File ;
    } else {
      fprintf( stderr , "[INPUT] FileType %s not recognised\n" ,
	       Flat[ io_tag ].Value ) ;
      Flag = FAILURE ;
    }
  }
  
  // get the filetype tag
  size_t an_tag = 0 ;
  if( ( an_tag = tag_search( Flat , "Analysis" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUT] Analysis tag not found in input file \n" ) ;
    Flag = FAILURE ;
  } else {
    if( are_equal( Flat[ an_tag ].Value , "Alphas" ) ) {
      Input -> Analysis = Alphas ;
    } else if( are_equal( Flat[ an_tag ].Value , "Adler" ) ) {
      Input -> Analysis = Adler ;
    } else if( are_equal( Flat[ an_tag ].Value , "Beta_crit" ) ) {
      Input -> Analysis = Beta_crit ;
    } else if( are_equal( Flat[ an_tag ].Value , "Binding_Corr" ) ) {
      Input -> Analysis = Binding_Corr ;
    } else if( are_equal( Flat[ an_tag ].Value , "Correlator" ) ) {
      Input -> Analysis = Correlator ;
    } else if( are_equal( Flat[ an_tag ].Value , "Exceptional" ) ) {
      Input -> Analysis = Exceptional ;
    } else if( are_equal( Flat[ an_tag ].Value , "HLBL" ) ) {
      Input -> Analysis = HLBL ;
    } else if( are_equal( Flat[ an_tag ].Value , "HVP" ) ) {
      Input -> Analysis = HVP ;
    } else if( are_equal( Flat[ an_tag ].Value , "Wflow" ) ) {
      Input -> Analysis = Wflow ;
    } else if( are_equal( Flat[ an_tag ].Value , "Fit" ) ) {
      Input -> Analysis = Fit ;
    } else if( are_equal( Flat[ an_tag ].Value , "Fpi_CLS" ) ) {
      Input -> Analysis = Fpi_CLS ;
    } else if( are_equal( Flat[ an_tag ].Value , "Nrqcd" ) ) {
      Input -> Analysis = Nrqcd ;
    } else if( are_equal( Flat[ an_tag ].Value , "PCAC" ) ) {
      Input -> Analysis = PCAC ;
    } else if( are_equal( Flat[ an_tag ].Value , "Pof" ) ) {
      Input -> Analysis = Pof ;
    } else if( are_equal( Flat[ an_tag ].Value , "Qcorr" ) ) {
      Input -> Analysis = Qcorr ;
    } else if( are_equal( Flat[ an_tag ].Value , "Qsusc" ) ) {
      Input -> Analysis = Qsusc ;
    } else if( are_equal( Flat[ an_tag ].Value , "Qslab" ) ) {
      Input -> Analysis = Qslab ;
    } else if( are_equal( Flat[ an_tag ].Value , "QslabFix" ) ) {
      Input -> Analysis = QslabFix ;
    } else if( are_equal( Flat[ an_tag ].Value , "Ren_Rats" ) ) {
      Input -> Analysis = Ren_Rats ;
    } else if( are_equal( Flat[ an_tag ].Value , "KKops" ) ) {
      Input -> Analysis = KKops ;
    } else if( are_equal( Flat[ an_tag ].Value , "KK_BK" ) ) {
      Input -> Analysis = KK_BK ;
    } else if( are_equal( Flat[ an_tag ].Value , "Sol" ) ) {
      Input -> Analysis = Sol ;
    } else if( are_equal( Flat[ an_tag ].Value , "General" ) ) {
      Input -> Analysis = General ;
    } else if( are_equal( Flat[ an_tag ].Value , "TetraGEVP" ) ) {
      Input -> Analysis = TetraGEVP ;
    } else if( are_equal( Flat[ an_tag ].Value , "TetraGEVP_Fixed" ) ) {
      Input -> Analysis = TetraGEVP_Fixed ;
    } else if( are_equal( Flat[ an_tag ].Value , "SpinOrbit" ) ) {
      Input -> Analysis = SpinOrbit ;
    } else if( are_equal( Flat[ an_tag ].Value , "StaticPotential" ) ) {
      Input -> Analysis = StaticPotential ;
    } else if( are_equal( Flat[ an_tag ].Value , "ZV" ) ) {
      Input -> Analysis = ZV ;
    } else {
      fprintf( stderr , "[INPUT] Analysis %s not recognised\n" ,
	       Flat[ an_tag ].Value ) ;
      Flag = FAILURE ;
    }
  }

  printf( "Here get_traj\n" ) ;
  
  if( get_traj( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }

  printf( "Here get fit\n" ) ;
  
#ifdef VERBOSE
  printf( "FLAG %d \n" , Flag ) ;
#endif
  if( get_fit( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }
#ifdef VERBOSE
  fprintf( stdout , "Fit FLAG %d \n" , Flag ) ;
#endif
  if( get_stats( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }
#ifdef VERBOSE
  fprintf( stdout , "FLAG %d \n" , Flag ) ;
#endif
  if( get_graph( Input , Flat , Ntags ) == FAILURE ) {
    Flag = FAILURE ;
  }
#ifdef VERBOSE
  fprintf( stdout , "FLAG %d \n" , Flag ) ;
#endif
  // close and unpack
  unpack_inputs( Flat , Ntags ) ;
  
  fclose( infile ) ;

  return Flag ;
}

// search for a specific tag
size_t
tag_search( const struct flat_file *Flat ,
	    const char *tag ,
	    const size_t Start ,
	    const size_t Ntags ) 
{
  size_t i ;
  for( i = Start ; i < Ntags ; i++ ) {
    if( are_equal( Flat[i].Token , tag ) ) {
      return i ;
    }
  }
  return i ;
}
