#ifndef SVD_H
#define SVD_H

int 
svd_inverse( double **Ainv , 
	     const double **A ,
	     const size_t NCOLS ,
	     const size_t NROWS ,
	     const double tolerance ,
	     const bool col_balance ) ;

#endif
