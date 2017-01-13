/**
   @file read_graph.c
   @brief read the input file for the graph information

   Graph = %s
   Graph_X = xaxis
   Graph_Y = yaxis
 */
#include "gens.h"

#include <string.h>

#include "read_inputs.h"

int
get_graph( struct input_params *Input ,
	   const struct flat_file *Flat ,
	   const size_t Ntags )
{
  size_t tag = 0 ;
  if( ( tag = tag_search( Flat , "Graph" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] Graph not found in input file!\n" ) ;
    return FAILURE ;
  }
  Input -> Graph.Name = malloc( Flat[tag].Value_Length * sizeof( char ) ) ;
  strcpy( Input -> Graph.Name , Flat[tag].Value ) ;
  
  if( ( tag = tag_search( Flat , "Graph_X" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] Graph_X not found in input file!\n" ) ;
    return FAILURE ;
  }
  Input -> Graph.Xaxis = malloc( Flat[tag].Value_Length * sizeof( char ) ) ;
  strcpy( Input -> Graph.Xaxis , Flat[tag].Value ) ;

  if( ( tag = tag_search( Flat , "Graph_Y" , 0 , Ntags ) ) == Ntags ) {
    fprintf( stderr , "[INPUTS] Graph_Y not found in input file!\n" ) ;
    return FAILURE ;
  }
  Input -> Graph.Yaxis = malloc( Flat[tag].Value_Length * sizeof( char ) ) ;
  strcpy( Input -> Graph.Yaxis , Flat[tag].Value ) ;

  fprintf( stdout , "\n[INPUTS] summary for graph information\n" ) ;
  fprintf( stdout , "[INPUTS] Graph name -> %s \n" , Input -> Graph.Name ) ;
  fprintf( stdout , "[INPUTS] Graph xaxis -> %s \n" , Input -> Graph.Xaxis ) ;
  fprintf( stdout , "[INPUTS] Graph yaxis -> %s \n" , Input -> Graph.Yaxis ) ;
  
  return SUCCESS ;
}
