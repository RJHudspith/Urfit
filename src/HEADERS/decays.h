#ifndef DECAYS_H
#define DECAYS_H

struct resampled
decay( const struct resampled *fitparams ,
       const struct input_params Input ,
       const size_t Mass_idx ,
       const size_t Amp_idx ) ;

struct resampled
decay_AP( const struct resampled *fitparams ,
	  const struct input_params Input ,
	  const size_t Mass_idx  ,
	  const size_t ALPW_idx ,
	  const size_t PWPW_idx ) ;

#endif
