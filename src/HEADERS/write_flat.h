#ifndef WRITE_FLAT_H
#define WRITE_FLAT_H

int
write_flat_file( const struct input_params Input ,
		 const char *name ) ;

int
write_flat_dist( const struct resampled *y ,
		 const struct resampled *x ,
		 const size_t Ndata ,
		 const char *name ) ;


int
write_flat_single( const struct resampled *y ,
		   const char *name ) ;

#endif
