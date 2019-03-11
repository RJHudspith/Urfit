#ifndef URFIT_NINT_H
#define URFIT_NINT_H

// numerical integrator over the distributions x and y
struct resampled
Nint( const struct resampled *dataX ,
      const struct resampled *dataY ,
      const size_t Ndata ,
      const bool is_trap ) ;

#endif
