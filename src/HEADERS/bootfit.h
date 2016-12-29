#ifndef BOOTFIT_H
#define BOOTFIT_H

struct resampled *
perform_bootfit( size_t *Npars ,
		 const struct resampled *x ,
		 const struct resampled *y ,
		 const double **W ,
		 const size_t *Ndata ,
		 const size_t Nsims ,
		 const bool *sims , // simultaneous fit indices
		 const size_t LT ,
		 const fittype fit ,
		 const corrtype CORRFIT ) ;

#endif
