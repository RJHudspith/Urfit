#ifndef AUTOCORR_H
#define AUTOCORR_H

int
autocorrelation( const struct resampled RAW ,
		 const size_t NSEP ,
		 const char *output ) ;

int
ACmeasure( const struct input_params Input ) ;

#endif
