#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

void
bootstrap_error( struct resampled *replicas ) ;

void
bootstrap_single( struct resampled *data ,
		  const size_t Nboots ) ;

void
bootstrap_full( struct input_params *Input ) ;

#endif
