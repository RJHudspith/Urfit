#ifndef READ_INPUTS_H
#define READ_INPUTS_H

int
are_equal( const char *str_1 ,
	   const char *str_2 ) ;

void
free_inputs( struct input_params *Input ) ;

int
read_inputs( struct input_params *Input ,
	     const char *filename ) ;

size_t
tag_search( const struct flat_file *Flat ,
	    const char *tag ,
	    const size_t Start ,
	    const size_t Ntags ) ;

#endif
