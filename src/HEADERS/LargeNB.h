#ifndef LARGENB_H
#define LARGENB_H

double
fLargeNB( const struct x_desc X , const double *fparams , const size_t Npars );
void
LargeNB_f( double *f , const void *data , const double *fparams );
void
LargeNB_df( double **df , const void *data , const double *fparams );
void
LargeNB_linmat( double **U ,
		const void *data ,
		const size_t N ,
		const size_t M ,
		const size_t Nlogic );
void
LargeNB_d2f( double **d2f , const void *data , const double *fparams );
void
LargeNB_guesses( double *fparams ,
		 const struct data_info Data ,
		 const struct fit_info Fit );
#endif
