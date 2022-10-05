#ifndef PLOT_FITFUNC_H
#define PLOT_FITFUNC_H

struct resampled
extrap_fitfunc( const struct resampled *f ,
		const struct data_info Data ,
		const struct fit_info Fit ,
		const double xpos ,
		const size_t shift ) ;

double
bisect_fitfunc( const struct resampled *f ,
		const struct data_info Data ,
		const struct fit_info Fit ,
		const double xlo ,
		const double target ,
		const double xhi ,
		const size_t shift ) ;

void
plot_fitfunction_HACK( const struct resampled *f ,
		       const struct data_info Data ,
		       const struct fit_info Fit ) ;

int
plot_feffmass( const struct resampled *f ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

void
plot_fitfunction( const struct resampled *f ,
		  const struct data_info Data ,
		  const struct fit_info Fit ) ;

#endif
