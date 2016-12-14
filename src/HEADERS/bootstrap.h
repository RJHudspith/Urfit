#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

void
free_boot_order( void ) ;

void
init_boot_order( const size_t NBOOTS , 
		 const size_t NRAW ) ;

void
bootstrap_error( struct resampled *replicas ) ;

void
bootstrap( struct resampled *Bootstrap ,
	   const struct resampled Raw ) ;

#endif
