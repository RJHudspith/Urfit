#ifndef JACKNIFE_H
#define JACKNIFE_H

void
jackknife_error( struct resampled *replicas ) ;

void
jackknife( struct resampled *Jackknife ,
	   const struct resampled Raw ) ;

#endif
