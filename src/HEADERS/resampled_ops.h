#ifndef RESAMPLED_OPS_H
#define RESAMPLED_OPS_H

void
add( struct resampled *a , 
     const struct resampled b ) ;

void
add_constant( struct resampled *a , 
	      const double b ) ;

void
divide( struct resampled *a , 
	const struct resampled b ) ;

void
divide_constant( struct resampled *a , 
		 const double b ) ;

void
equate( struct resampled *a , 
	const struct resampled b ) ;
void
equate_constant( struct resampled *a ,
		 const double constant ,
		 const size_t NSAMPLES ,
		 const resample_type restype ) ;

struct resampled
init_dist( const struct resampled *d , 
	   const size_t NSAMPLES , 
	   const resample_type restype ) ;

void
mult_constant( struct resampled *a , 
	       const double b ) ;

void
mult( struct resampled *a , 
      const struct resampled b ) ;

void
raise( struct resampled *a , 
       const double b ) ;

void
rapby( struct resampled *a ,
       const struct resampled b ,
       const double y ) ;

void
res_acosh( struct resampled *a ) ;

void
res_asinh( struct resampled *a ) ;

void
res_atanh( struct resampled *a ) ;

void
res_log( struct resampled *a ) ;

void
res_exp( struct resampled *a ) ;

void
root( struct resampled *a ) ;

void
spin_average( struct resampled *a , 
	      const struct resampled b ) ;

void
subtract( struct resampled *a , 
	  const struct resampled b ) ;

void
subtract_constant( struct resampled *a , 
		   const double b ) ;

void
phi4_comb( struct resampled *a , 
	   const struct resampled b ) ;

#endif
