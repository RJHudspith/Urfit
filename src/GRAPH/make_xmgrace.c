/**
   @file make_xmgrace.c
   @brief draw an xmgrace graph
 */
#include <stdlib.h>
#include <stdio.h>

#include "gens.h"
#include "fitfunc.h"

static FILE *file ;
static size_t dataset = 0 , colorset = 1 ;

static void
PrintHeader( void )
{
  const size_t symbol = ( dataset + 1 ) % 11 ;
  fprintf( file , "@\ts%zu hidden false\n" , dataset ) ;  
  fprintf( file , "@\ts%zu symbol %zu\n" , dataset , symbol ) ;  
  fprintf( file , "@\ts%zu symbol size 0.62\n" , dataset ) ;  
  fprintf( file , "@\ts%zu symbol color %zu\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%zu symbol fill color %zu\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%zu symbol fill pattern 1\n" , dataset ) ;  
  fprintf( file , "@\ts%zu line type 0\n" , dataset ) ;  
  fprintf( file , "@\ts%zu line color %zu\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%zu errorbar on\n" , dataset ) ;  
  fprintf( file , "@\ts%zu errorbar place both\n" , dataset ) ;  
  fprintf( file , "@\ts%zu errorbar color %zu\n" , dataset , colorset ) ;  
  fprintf( file , "@\ts%zu errorbar size 1\n" , dataset ) ; 
  return ;
}

static void
initialise_graph( const char *x_axis , 
		  const char *y_axis )
{
  fprintf( file , "@map font 0 to \"Times-Roman\", \"Times-Roman\"\n" ) ;
  fprintf( file , "@map font 1 to \"Symbol\", \"Symbol\"\n" ) ;
  fprintf( file , "@map font 2 to \"Helvetica\", \"Helvetica\"\n" ) ;
  fprintf( file , "@map color 0 to (255,255,255), \"white\"\n" ) ;
  fprintf( file , "@map color 1 to (0,0,0), \"black\"\n" ) ;
  fprintf( file , "@map color 2 to (255,0,0), \"red\"\n" ) ;
  fprintf( file , "@map color 3 to (0,100,0), \"dark green\"\n" ) ;
  fprintf( file , "@map color 4 to (0,0,255), \"blue\"\n" ) ;
  fprintf( file , "@map color 5 to (188,143,143), \"rosy brown\"\n" ) ;
  fprintf( file , "@map color 6 to (178,34,34), \"firebrick\"\n" ) ;
  fprintf( file , "@map color 7 to (205,92,92), \"indian red\"\n" ) ;
  fprintf( file , "@map color 8 to (46,139,87), \"sea green\"\n" ) ;
  fprintf( file , "@map color 9 to (255,69,0), \"orange red\"\n" ) ;
  fprintf( file , "@map color 10 to (0,128,128), \"teal\"\n" ) ;
  fprintf( file , "@map color 11 to (128,128,0), \"olive\"\n" ) ;
  fprintf( file , "@map color 12 to (184,134,11), \"dark golden rod\"\n" ) ;
  fprintf( file , "@map color 13 to (47,79,79), \"dark slate gray\"\n" ) ;
  fprintf( file , "@map color 14 to (128,0,128), \"purple\"\n" ) ;
  fprintf( file , "@map color 15 to (160,82,45), \"sienna\"\n" ) ;
  fprintf( file , "@default linestyle 1\n" ) ;
  fprintf( file , "@default linewidth 1.0\n" ) ;
  fprintf( file , "@default color 1\n" ) ;
  fprintf( file , "@default pattern 1\n" ) ;
  fprintf( file , "@default font 0\n" ) ;
  fprintf( file , "@background color 0\n" ) ;
  fprintf( file , "@page background fill on\n" ) ;
  fprintf( file , "@g0 on\n" ) ;
  fprintf( file , "@g0 type XY\n" ) ;
  fprintf( file , "@with g0\n" ) ;
  fprintf( file , "@\tview 0.180000, 0.180000, 1.250000, 0.950000\n" ) ;
  fprintf( file , "@\txaxes invert off\n" ) ;  
  fprintf( file , "@\tyaxes invert off\n" ) ;  
  fprintf( file , "@\txaxes scale Normal\n" ) ;  
  fprintf( file , "@\tyaxes scale Normal\n" ) ;  
  fprintf( file , "@\txaxis on\n" ) ;  
  fprintf( file , "@\txaxis\tlabel char size 1.490000\n" ) ;  
  fprintf( file , "@\txaxis\ttick minor ticks 0\n" ) ; 
  fprintf( file , "@\txaxis\ttick out\n" ) ;  
  fprintf( file , "@\tyaxis on\n" ) ;  
  fprintf( file , "@\tyaxis\tlabel char size 1.490000\n" ) ; 
  fprintf( file , "@\tyaxis\ttick minor ticks 0\n" ) ; 
  fprintf( file , "@\tyaxis\ttick out\n" ) ; 
  // X axis label
  fprintf( file , "@\txaxis label \"%s\" \n" , x_axis ) ;  
  // Y axis label
  fprintf( file , "@\tyaxis label \"%s\" \n" , y_axis ) ;  
  fprintf( file , "@\txaxis label place auto\n" ) ;  
  fprintf( file , "@\tyaxis label place auto\n" ) ;  
  fprintf( file , "@\taltxaxis\toff\n" ) ;  
  fprintf( file , "@\taltyaxis\toff\n" ) ;  
  fprintf( file , "@\tlegend on\n" ) ;  
  fprintf( file , "@\tlegend 0.75, 0.94\n" ) ;  
  fprintf( file , "@\tlegend char size 1.50000\n") ;  
  return ;
}

// and close
void
close_xmgrace_graph( void )
{
  dataset = 0 ; colorset = 1 ;
  fclose( file ) ;
  return ;
}

// assumes graph is already open
void
draw_line( const double *x , const double *y , const size_t n )
{
  fprintf( file , "@\ts%zu hidden false\n" , dataset ) ;  
  fprintf( file , "@\ts%zu symbol 0\n" , dataset ) ;  
  fprintf( file , "@\ts%zu line color 1\n" , dataset ) ;  
  fprintf( file , "@\ts%zu errorbar color 1\n" , dataset ) ;  
  fprintf( file , "@\ts%zu legend \"\"\n" , dataset ) ;  
  fprintf( file , "@target G0.S%zu\n" , dataset ) ;
  fprintf( file , "@type xy\n" ) ;

  size_t i ;
  for( i = 0 ; i < n ; i++ ) {
    fprintf( file , "%e %e\n" , x[i] , y[i] ) ;
  }
  fprintf( file , "&\n" ) ;

  dataset ++ ;

  return ;
}

// open
void
make_xmgrace_graph( const char *filename ,
		    const char *x_axis , 
		    const char *y_axis )
{
  printf( "\n[GRAPH] name :: %s \n" , filename ) ;
  printf( "[GRAPH] axes (x) %s (y) %s \n" , x_axis , y_axis ) ;
  file = fopen( filename , "w" ) ;
  initialise_graph( x_axis , y_axis ) ;
  return ;
}

// plot the data
void
plot_data( const struct resampled *x ,
	   const struct resampled *y ,
	   const size_t Ndata )
{
  // plot the data
  PrintHeader( ) ;

  fprintf( file , "@\ts%zu legend \"\" \n" , dataset ) ;  
  fprintf( file , "@target G0.S%zu\n" , dataset ) ;
  fprintf( file , "@type xydxdxdydy\n" ) ;

  size_t i ;
  for( i = 0 ; i < Ndata ; i++ ) {
    fprintf( file , "%e %e %e %e %e %e\n" ,
	     x[i].avg , 
	     y[i].avg ,
	     x[i].err_hi - x[i].avg , 
	     x[i].avg - x[i].err_lo ,
	     y[i].err_hi - y[i].avg ,
	     y[i].avg - y[i].err_lo ) ;
  }
  fprintf( file , "&\n" ) ;
  dataset ++ ;
  colorset = ( colorset + 1 ) % 16 ;
  return ;
}
