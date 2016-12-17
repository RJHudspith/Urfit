#ifndef GLU_BSWAP_H
#define GLU_BSWAP_H

void
bswap_16( const size_t n , void *u ) ;

void
bswap_32( const size_t n , void *u ) ;

void
bswap_64( const size_t n , void *u ) ;

int 
is_big_endian( void ) ;

#endif
