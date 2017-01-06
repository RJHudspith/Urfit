#ifndef INIT_H
#define INIT_H

void
free_Data( struct data_info *Data ) ;

void
free_Fit( struct fit_info *Fit ,
	  const struct data_info Data ) ;

#endif
