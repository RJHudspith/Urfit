#ifndef CRC32C_H
#define CRC32C_H

void 
DML_checksum_accum_crc32c( uint32_t *checksuma , 
			   uint32_t *checksumb , 
			   const uint32_t rank , 
			   const void *buf , 
			   const size_t size ) ;

#endif
