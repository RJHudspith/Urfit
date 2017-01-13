/**
   @file read_traj.c
   @brief read in the trajectory information from the flat input file

   Expected trajectory information is (in this order)

   TrajName = %s
   TrajStep = Begin,Increment,End
   TrajStat = Bin
   TrajFitr = Fit_Low,Fit_High
   TrajDims = %d,%d,%d,%d...
 */
#include "gens.h"

#include <string.h>

#include "read_inputs.h"

// block size for the traj information
#define Nblock (5)

enum { TrajName , TrajStep , TrajStat , TrajFitr , TrajDims } TrajBlock ;

// linked list for the dimensions
struct node {
  size_t dim ;
  struct node *next ;
} ;

// set the record for each traj file
static char **
set_record( void )
{
  char **Record = malloc( Nblock * sizeof( char* ) ) ;
  Record[0] = "TrajName" ;
  Record[1] = "TrajStep" ;
  Record[2] = "TrajStat" ;
  Record[3] = "TrajFitr" ;
  Record[4] = "TrajDims" ;
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

  size_t *Block = malloc( *Ntraj * sizeof( size_t ) ) ;
  *Ntraj = 0 ; block_idx = 0 ;
  while( ( block_idx = tag_search( Flat , Record[ TrajName ] , block_idx , Ntags ) )
	 != Ntags ) {
    Block[ *Ntraj ] = block_idx ;
    *Ntraj = *Ntraj + 1 ; block_idx += Nblock ;
  }
  
  return Block ;
}

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
    // set the filename
    Traj[i].Filename = malloc( Flat[ Block[i] ].Value_Length * sizeof( char ) ) ;
    sprintf( Traj[i].Filename , "%s" , Flat[ Block[i] ].Value ) ;
    
    // set the file start, increment, and end
    tok = strtok( Flat[ Block[i] + TrajStep ].Value , "," ) ;
    Traj[i].Begin     = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;
    Traj[i].Increment = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;
    Traj[i].End       = strtol( tok , &endptr , 10 ) ; tok = strtok( NULL , "," ) ;

    // set the binning
    Traj[i].Bin = strtol( Flat[ Block[i] + TrajStat ].Value , &endptr , 10 ) ;

    // set the fit range
    tok = strtok( Flat[ Block[i] + TrajFitr ].Value , "," ) ;
    Traj[i].Fit_Low  = strtod( tok , &endptr ) ; tok = strtok( NULL , "," ) ;
    Traj[i].Fit_High = strtod( tok , &endptr ) ; tok = strtok( NULL , "," ) ;
    
    // set the dimensions, linked list is backwards
    size_t Ndims = 0 , j ;
    struct node *Node = get_Ndims( &Ndims , Flat[ Block[i] + TrajDims ].Value ) ;
    Traj[i].Nd = Ndims ;
    Traj[i].Dimensions = malloc( Traj[i].Nd * sizeof( size_t ) ) ;
    for( j = Traj[i].Nd ; j != 0 ; j-- ) {
      Traj[i].Dimensions[j-1] = Node -> dim ;
      free( Node ) ;        // free this node
      Node = Node -> next ; // point to the next node
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

  // set the trajectory information into the Input struct
  if( ( Input -> Traj = set_trajs( Flat , Block , Input -> Data.Nsim ) )
      == NULL ) {
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
    printf( "(Filename) -> %s \n" , Input -> Traj[i].Filename ) ;
    size_t j ;
    printf( "(Nd,Dims) ->  ( %zu , " , Input -> Traj[i].Nd ) ;
    for( j = 0 ; j < Input -> Traj[i].Nd ; j++ ) {
      printf( " %zu " , Input -> Traj[i].Dimensions[j] ) ;
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
