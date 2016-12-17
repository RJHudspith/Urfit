#ifndef MAKE_XMGRACE_H
#define MAKE_XMGRACE_H

void
close_xmgrace_graph( void ) ;

void
draw_line( const double *x , 
	   const double *y , 
	   const size_t n ) ;

void
make_xmgrace_graph( const char *filename ,
		    const char *x_axis , 
		    const char *y_axis ) ;

void
plot_data( const struct resampled *x ,
	   const struct resampled *y ,
	   const size_t Ndata ) ;

#endif
