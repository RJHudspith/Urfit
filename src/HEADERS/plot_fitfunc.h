#ifndef PLOT_FITFUNC_H
#define PLOT_FITFUNC_H

void
plot_fitfunction( const struct resampled *f ,
		  const fittype fit ,
		  const struct resampled *x ,
		  const size_t *Ndata ,
		  const size_t LT ,
		  const corrtype CORRFIT ,
		  const size_t Nsims ,
		  const bool *sims ) ;

#endif
