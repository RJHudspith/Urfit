#ifndef INIT_H
#define INIT_H

void
free_Data( struct data_info *Data ,
	   const struct fit_info Fit ) ;

void
free_Fit( struct fit_info *Fit ,
	  const struct data_info Data ) ;

int
init_LT( struct data_info *Data ,
	 const struct traj *Traj ) ;

#endif
