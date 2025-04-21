#ifndef FAKE_H
#define FAKE_H

struct resampled
generate_fake_single( const double y,
		      const double dy ,
		      const size_t Nsamples ) ;

struct resampled *
generate_fake_boot( const size_t N,
		    const size_t Nsamples,
		    const double yarr[N] ,
		    const double dyarr[N] ) ;

int
generate_fake_data( struct data_info *Data ,
		    struct fit_info Fit ,
		    struct traj *Traj ,
		    const double xsigma ,
		    const double ysigma ) ;

#endif
