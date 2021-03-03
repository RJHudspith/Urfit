/**
   @file read_traj.c
   @brief read in the trajectory information from the flat input file

   Expected trajectory information is (in this order)

   TrajName = %s
   TrajStep = Begin,Increment,End
   TrajStat = Bin
   TrajFitr = Fit_Low,Fit_High
   TrajDims = %d,%d,%d,%d
   TrajGsGK = Gamma1,Gamma_2,Tfold
   TrajMom = 0.,0.,0.
   TrajRW = %s
 */
#include "gens.h"

#include <string.h>

#include "read_inputs.h"

// block size for the traj information
#define Nblock (8)

// little enum for getting the map right
enum { TrajName , TrajStep , TrajStat ,
       TrajFitr , TrajDims , TrajGsGk , TrajMom , TrajRW } TrajBlock ;

// linked list for the dimensions
struct node {
  int dim ;
  struct node *next ;
} ;

struct node_dbl {
  double dim ;
  struct node_dbl *next ;
} ;

// set the record for each traj file
static char **
set_record( void )
{
  char **Record = malloc( Nblock * sizeof( char* ) ) ;
  Record[0] = "TrajXY" ;
  Record[1] = "TrajStep" ;
  Record[2] = "TrajStat" ;
  Record[3] = "TrajFitr" ;
  Record[4] = "TrajDims" ;
  Record[5] = "TrajGsGk" ;
  Record[6] = "TrajMom" ;
  Record[7] = "TrajRW" ;
  return Record ;
}

// free the record for each traj file
static void
free_record( char **Record )
{
  free( Record ) ;
}

// get the number of trajectories, returns 0 if there are incomplete records
static size_t*
get_Ntraj( size_t *Ntraj ,
	   const struct flat_file *Flat ,
	   const char **Record ,
	   const size_t Ntags )
{ 
  size_t block_idx = 0 ;
  
  size_t *Block = NULL ;
  *Ntraj = 0 ;

  // loop the indexes sanity-checking the records
  while( ( block_idx = tag_search( Flat , Record[0] , block_idx , Ntags ) )
	 != Ntags ) {
    // sanity check that we have enough tags
    if( ( block_idx + TrajDims ) >= Ntags ) {
      fprintf( stderr , "[INPUTS] insufficient record length for traj %zu\n" ,
	       *Ntraj ) ;
      fprintf( stderr , "[INPUTS] expected TrajName\nTrajStep\n" 
	       "TrajStat\nTrajFitr\nTrajDims\n" ) ;
      return NULL ;
    }
    size_t j ;
    for( j = 1 ; j < Nblock ; j++ ) {
      if( !are_equal( Flat[ block_idx + j ].Token , Record[j] ) ) {
	fprintf( stderr , "[INPUTS] %s not found for traj %zu\n" ,
		 Record[j] , *Ntraj ) ;
	fprintf( stderr , "[INPUTS] is the ordering correct?\n" ) ;
        return NULL ;
      }
    }
    // increment our number of measurements and the block index
    *Ntraj = *Ntraj + 1 ; block_idx += Nblock ;
  }

  if( *Ntraj == 0 ) return Block ;
  
  Block = malloc( *Ntraj * sizeof( size_t ) ) ;
  *Ntraj = 0 ; block_idx = 0 ;
  while( ( block_idx = tag_search( Flat , Record[ TrajName ] , block_idx , Ntags ) )
	 != Ntags ) {
    Block[ *Ntraj ] = block_idx ;
    *Ntraj = *Ntraj + 1 ; block_idx += Nblock ;
  }
  
  return Block ;
}

// get the number of dimensions
static struct node *
get_Ndims( size_t *Nd , const char *Value )
{
  struct node *head = NULL , *curr ;
  *Nd = 1 ;
  char *tok = strtok( (char*)Value , "," ) , *endptr ;
  curr = (struct node*)malloc( sizeof( struct node ) ) ;
  curr -> dim = strtol( tok , &endptr , 10 ) ;
  curr -> next = head ;
  head = curr ;
  
  while( ( tok = strtok( NULL , "," ) ) != NULL ) {
    curr = (struct node*)malloc( sizeof( struct node ) ) ;
    curr -> dim = strtol( tok , &endptr , 10 ) ;
    curr -> next = head ;
    head = curr ;
    
    *Nd = *Nd + 1 ;
  }
  return head ;
}

// get the number of dimensions
static struct node_dbl*
get_Moms( size_t *Nd , const char *Value )
{
  struct node_dbl *head = NULL , *curr ;
  *Nd = 1 ;
  char *tok = strtok( (char*)Value , "," ) , *endptr ;
  curr = (struct node_dbl*)malloc( sizeof( struct node_dbl ) ) ;
  curr -> dim = strtod( tok , &endptr ) ;
  curr -> next = head ;
  head = curr ;
  while( ( tok = strtok( NULL , "," ) ) != NULL ) {
    curr = (struct node_dbl*)malloc( sizeof( struct node_dbl ) ) ;
    curr -> dim = strtod( tok , &endptr ) ;
    curr -> next = head ;
    head = curr ;
    
    *Nd = *Nd + 1 ;
  }
  return head ;
}

// Gamma index map
static size_t
GsGk_map( const char *tok )
{
  char *endptr ;
  if( are_equal( tok , "Vi" ) ) {
    return Vi ;
  } else if( are_equal( tok , "Tit" ) ) {
    return Tit ;
  } else if( are_equal( tok , "Ai" ) ) {
    return Ai ;
  } else if( are_equal( tok , "Tij" ) ) {
    return Tij ;
  } 
  return strtol( tok , &endptr , 10 ) ;
}

// set the gammas and the time folding
static int
set_GsGk( struct traj *Traj ,
	  char *Tok )
{
  // tokenize gamma source and gamma sink
  char *tok = strtok( Tok , "," ) ;
  Traj -> Gs = GsGk_map( tok ) ; tok = strtok( NULL , "," ) ;
  Traj -> Gk = GsGk_map( tok ) ; tok = strtok( NULL , "," ) ;
  // set the fold type
  if( are_equal( tok , "PLUS_PLUS" ) ) {
    Traj -> Fold = PLUS_PLUS ;
  } else if( are_equal( tok , "PLUS_MINUS" ) ) {
    Traj -> Fold = PLUS_MINUS ;
  } else if( are_equal( tok , "MINUS_PLUS" ) ) {
    Traj -> Fold = MINUS_PLUS ;
  } else if( are_equal( tok , "MINUS_MINUS" ) ) {
    Traj -> Fold = MINUS_MINUS ;
  } else if( are_equal( tok , "NOFOLD" ) ) {
    Traj -> Fold = NOFOLD ;
  } else if( are_equal( tok , "NOFOLD_MINUS" ) ) {
    Traj -> Fold = NOFOLD_MINUS ;
  } else if( are_equal( tok , "TDER" ) ) {
    Traj -> Fold = NOFOLD_MINUS ;
  } else if( are_equal( tok , "NOFOLD_SWAPT" ) ) {
    Traj -> Fold = NOFOLD_SWAPT ;
  } else if( are_equal( tok , "NOFOLD_MINUS_SWAPT" ) ) {
    Traj -> Fold = NOFOLD_MINUS_SWAPT ;
  } else if( are_equal( tok , "ZV_SUB" ) ) {
    Traj -> Fold = ZV_SUB ;
  } else {
    fprintf( stderr , "[INPUT] I don't understand the fold %s\n" ,
	     tok ) ;
    return FAILURE ;
  }
  return SUCCESS ;
}

// set the trajectory structs
static struct traj *
set_trajs( const struct flat_file *Flat ,
	   const size_t *Block ,
	   const size_t Ntraj )
{
  struct traj *Traj = malloc( Ntraj * sizeof( struct traj ) ) ;
  char *tok , *endptr = NULL ;
  size_t i ;
  for( i = 0 ; i < Ntraj ; i++ ) {
    
    // set the filename for the x data
    tok = strtok( Flat[ Block[i] + TrajName ].Value , "," ) ;
    if( are_equal( tok , "NULL" ) ) {
      Traj[i].FileX = NULL ;
    } else {
      Traj[i].FileX = malloc( strlen( tok ) * sizeof( char ) ) ;
      sprintf( Traj[i].FileX , "%s" , tok ) ;
    }
    printf( "%s \n" , tok ) ;
    tok = strtok( NULL , "," ) ;

    // set the filename for the y data
    if( are_equal( tok , "NULL" ) ) {
      Traj[i].FileY = NULL ;
    } else {
      Traj[i].FileY = malloc( strlen( tok ) * sizeof( char ) ) ;
      sprintf( Traj[i].FileY , "%s" , tok ) ;
    }
    printf( "%s \n" , tok ) ;

    // set the file start, increment, and end. Sanity check them
    tok = strtok( Flat[ Block[i] + TrajStep ].Value , "," ) ;
    Traj[i].Begin     = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;
    Traj[i].Increment = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;
    Traj[i].End       = strtol( tok , &endptr , 10 ) ;
    if( Traj[i].Begin >= Traj[i].End ) {
      fprintf( stderr , "[INPUT] non-sensical trajectory beginning and end"
	       "%zu -> %zu\n" , Traj[i].Begin , Traj[i].End ) ;
      return NULL ;
    }

    // set the binning
    Traj[i].Bin = strtol( Flat[ Block[i] + TrajStat ].Value , &endptr , 10 ) ;

    // set the fit range
    tok = strtok( Flat[ Block[i] + TrajFitr ].Value , "," ) ;
    Traj[i].Fit_Low  = strtod( tok , &endptr ) ; tok = strtok( NULL , "," ) ;
    Traj[i].Fit_High = strtod( tok , &endptr ) ;
    
    // set the dimensions, linked list is backwards
    size_t Ndims = 0 , j ;
    struct node *Node = get_Ndims( &Ndims , Flat[ Block[i] + TrajDims ].Value ) ;
    Traj[i].Nd = Ndims ;
    Traj[i].Dimensions = malloc( Traj[i].Nd * sizeof( size_t ) ) ;
    for( j = Traj[i].Nd ; j != 0 ; j-- ) {
      Traj[i].Dimensions[j-1] = (size_t)Node -> dim ;
      //free( Node ) ;        // free this node
      Node = Node -> next ; // point to the next node
    }

    // set the GSGK
    if( set_GsGk( &Traj[i] , Flat[ Block[i] + TrajGsGk ].Value ) == FAILURE ) {
      return NULL ;
    }

    // set the momenta linked list is backwards
    struct node_dbl *ndbl = get_Moms( &Ndims , Flat[ Block[i] + TrajMom ].Value ) ;
    Traj[i].mom = malloc( 4 * sizeof( double ) ) ;
    for( j = Ndims ; j != 0 ; j-- ) {
      Traj[i].mom[j-1] = (double)ndbl -> dim ;
      //free( ndbl ) ;        // free this node
      ndbl = ndbl -> next ; // point to the next node
    }

    // set the RW name
    tok = Flat[ Block[i] + TrajRW ].Value ;
    if( are_equal( tok , "NULL" ) ) {
      Traj[i].RW = NULL ;
    } else {
      Traj[i].RW = malloc( strlen( tok ) * sizeof( char ) ) ;
      sprintf( Traj[i].RW , "%s" , tok ) ;
    }
  }
  return Traj ;
}

// get all of the trajectory information from the input file
int
get_traj( struct input_params *Input ,
	  const struct flat_file *Flat ,
	  const size_t Ntags )
{
  // set the record for what we are supposed to look for
  char **Record = set_record( ) ;
  int flag = SUCCESS ;
  size_t i ;
  Input -> Data.Nsim = 0 ;

  // check the consistency of the trajectories
  const size_t *Block = get_Ntraj( &( Input -> Data.Nsim ) , Flat ,
				   (const char**)Record , Ntags ) ;
  if( Block == NULL ) {
    flag = FAILURE ;
    goto memfree ;
  }

  printf( "In traj set_trajs\n" ) ;

  // set the trajectory information into the Input struct
  if( ( Input -> Traj = set_trajs( Flat , Block , Input -> Data.Nsim ) )
      == NULL ) {
    fprintf( stderr , "[INPUT] set_trajs failure\n" ) ;
    flag = FAILURE ;
    goto memfree ;
  }

  // summarise the trajectories
  for( i = 0 ; i < Input -> Data.Nsim ; i++ ) {
    printf( "\n[INPUT] Summary for Traj %zu \n" , i ) ;
    printf( "(Begin, Inc , End ) -> ( %zu , %zu , %zu ) \n" ,
	    Input -> Traj[i].Begin , Input -> Traj[i].Increment ,
	    Input -> Traj[i].End ) ;
    printf( "(Fit Low, Fit High ) -> ( %f , %f ) \n" ,
	    Input -> Traj[i].Fit_Low , Input -> Traj[i].Fit_High ) ;
    printf( "(Binning) -> %zu \n" , Input -> Traj[i].Bin ) ;
    printf( "(FileX) -> %s \n" , Input -> Traj[i].FileX ) ;
    printf( "(FileY) -> %s \n" , Input -> Traj[i].FileY ) ;
    printf( "(RW) -> %s \n" , Input -> Traj[i].RW ) ;
    size_t j ;
    printf( "(Nd,Dims) ->  ( %zu , " , Input -> Traj[i].Nd ) ;
    for( j = 0 ; j < Input -> Traj[i].Nd ; j++ ) {
      printf( " %zu " , Input -> Traj[i].Dimensions[j] ) ;
    }
    printf( ") \n" ) ;
    printf( "(Mom) ->  ( " ) ;
    for( j = 0 ; j < 3 ; j++ ) {
      printf( " %1.15f " , Input -> Traj[i].mom[j] ) ;
    }
    printf( ") \n" ) ;
  }
  
  // memory free
 memfree :
  
  if( Record != NULL ) {
    free_record( Record ) ;
  }
  if( Block != NULL ) {
    free( (void*)Block ) ;
  }
  
  return flag ;
}

#undef Nblock
