#ifndef URFIT_NINT_H
#define URFIT_NINT_H

// numerical integrator over the distributions x and y
struct resampled
Nint( const struct resampled *dataX ,
      const struct resampled *dataY ,
      const size_t Ndata ,
      const bool is_trap ) ;

// numerical integrator over the distributions x and y
struct resampled
Nint_pt( const struct resampled *dataX ,
	 const struct resampled *dataY ,
	 const size_t Ndata ,
	 const double pt ,
	 const bool put_zero ) ;

// numerically integrate using the fit parameters
struct resampled
Nint_fit( struct resampled *f ,
	  const struct data_info Data ,
	  const struct fit_info Fit ,
	  const double upp ,
	  const double low ,
	  const double eps ,
	  const size_t shift ) ;

#endif
