/**
   @file nonexceptional.c
   @brief fit a lambda and extrapolate to various values
 */
#include "gens.h"

int
fit_nonexceptional( struct input_params *Input )
{
  double chisq ;
  struct resampled *fit = fit_and_plot( *Input , &chisq ) ;
  return SUCCESS ;
}
