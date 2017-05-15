#ifndef NDER_H
#define NDER_H

double
Nder( double (*f) ( const struct x_desc X ,
		    const double *fparams ,
		    const size_t Npars ) ,
      const struct x_desc X ,
      const size_t Npars ,
      const double *fparams ,
      const size_t idx ,
      const size_t Nparams ) ;

#endif
