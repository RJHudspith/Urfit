#ifndef FIT_AND_PLOT_H
#define FIT_AND_PLOT_H

struct resampled *
fit_and_plot( struct input_params Input ,
	      double *Chi ) ;

struct resampled *
fit_and_plot_and_Nint( struct input_params Input ,
		       double *Chi ) ;

#endif
