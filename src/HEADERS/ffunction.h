#ifndef FFUNCTION_H
#define FFUNCTION_H

struct ffunction
allocate_ffunction( const size_t NPARAMS ,
		    const size_t NDATA ) ;

void
copy_ffunction( struct ffunction *f1 ,
		const struct ffunction f ) ;

void
free_ffunction( struct ffunction *f , 
		const size_t NPARAMS ) ;

#endif
