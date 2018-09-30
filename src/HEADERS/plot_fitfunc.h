#ifndef PLOT_FITFUNC_H
#define PLOT_FITFUNC_H

// derivative evaluated at xpos
struct resampled
deriv_fitfunc( const struct resampled *f ,
	       const struct data_info Data ,
	       const struct fit_info Fit ,
	       const double xpos ,
	       const size_t shift ) ;

struct resampled
extrap_fitfunc( const struct resampled *f ,
		const struct data_info Data ,
		const struct fit_info Fit ,
		const double xpos ,
		const size_t shift ) ;

void
plot_fitfunction( const struct resampled *f ,
		  const struct data_info Data ,
		  const struct fit_info Fit ) ;

#endif
