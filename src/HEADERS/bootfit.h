#ifndef BOOTFIT_H
#define BOOTFIT_H

struct resampled *
perform_bootfit( const struct resampled *x ,
		 const struct resampled *y ,
		 const double **W ,
		 const size_t Ndata ,
		 const size_t LT ,
		 const fittype fit ,
		 const corrtype CORRFIT ) ;

#endif
