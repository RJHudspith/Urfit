#ifndef FIT_CHOOSER_H
#define FIT_CHOOSER_H

size_t
get_Nparam( const struct fit_info Fit ) ;

struct fit_descriptor
init_fit( const struct data_info Data ,
	  const struct fit_info Fit ) ;

#endif
