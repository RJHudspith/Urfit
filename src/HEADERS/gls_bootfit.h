#ifndef GLS_BOOTFIT_H
#define GLS_BOOTFIT_H

int
single_gls( double *coeffs ,
	    double *chisq ,
	    const struct data_info Data ,
	    const struct fit_info Fit , 
	    const size_t sample_idx ,
	    const bool is_average ) ;

struct resampled *
GLS( const struct data_info Data ,
     const struct fit_info Fit ) ;

#endif
