TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libURFIT.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

FITS_FILES=./FITS/exp.c

MINIMIZE_FILES=./MINIMIZE/CG.c ./MINIMIZE/GA.c ./MINIMIZE/line_search.c \
	./MINIMIZE/LM.c ./MINIMIZE/SD.c

UTILS_FILES=./UTILS/chisq.c ./UTILS/ffunction.c

## all the source files apart from ./Run/Mainfile.c
libURFIT_a_SOURCES = \
	${FITS_FILES} \
	${MINIMIZE_FILES} \
	${UTILS_FILES}

## and their respective headers
include_HEADERS = 
	$(FITS_FILES:./FITS/%.c=./HEADERS/%.h) \
	$(MINIMIZE_FILES:./MINIMIZE/%.c=./HEADERS/%.h) \
	$(UTILS_FILES:./UTILS/%.c=./HEADERS/%.h) \
	./HEADERS/gens.h \
	../config.h

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = URFIT

URFIT_SOURCES = Mainfile.c
URFIT_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
URFIT_LDADD = libURFIT.a ${LDFLAGS}

endif