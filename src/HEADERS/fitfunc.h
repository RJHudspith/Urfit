#ifndef FITFUNC_H
#define FITFUNC_H

// error types
enum { ERR , HI , LO , AVE } ;

// this is the same order as UKhadron
// where did the 1 go?
typedef enum { RAWDATA , JACKDATA , BOOTDATA } resample_type ;

// struct containing our statistics
struct resampled {
  double *resampled ;
  double avg ;
  double err_hi ;
  double err_lo ;
  double err ;
  size_t NSAMPLES ;
  resample_type restype ;
} ;

#endif
