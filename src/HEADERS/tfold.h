#ifndef TFOLD_H
#define TFOLD_H

double complex *
map_correlator( struct traj Traj ,
		const char *str ,
		const uint32_t *mompoint ,
		const size_t Nlt ) ;

int
time_fold( struct resampled *sample ,
	   const double complex *C ,
	   const size_t LT ,
	   const fold_type fold ,
	   const size_t meas ) ;

#endif
