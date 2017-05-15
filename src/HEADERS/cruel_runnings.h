#ifndef CRUEL_RUNNINGS_H
#define CRUEL_RUNNINGS_H

double
RUN( double mu ,
     const double alpha_mu ,
     const double muprime , 
     const size_t nf ,
     const size_t loops ) ;

double 
run_nf3_2MZ( const double alpha ,
	     const double mu ,
	     const size_t loops ) ;

double 
run_MZ_2nf3( const double alpha ,
	     const double mu ,
	     const size_t loops ) ;

struct resampled
run_distribution_nf3_2MZ( struct resampled alpha ,
			  const double mu ,
			  const size_t loops ) ;
  
void
test_running( void ) ;

#endif
