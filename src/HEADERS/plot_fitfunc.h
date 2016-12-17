#ifndef PLOT_FITFUNC_H
#define PLOT_FITFUNC_H

void
plot_fitfunction( const struct fitparams *f ,
		  const fittype fit ,
		  const struct resampled *x ,
		  const size_t Ndata ,
		  const size_t LT ,
		  const corrtype CORRFIT ) ;

#endif
