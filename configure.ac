AC_INIT([URFIT], [0.1], [renwick.james.hudspith@gmail.com])
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AM_PROG_AR
AC_C_BIGENDIAN 

## SET up the m4
AC_CONFIG_MACRO_DIR([m4])

## prefix!!
AC_DEFINE_UNQUOTED([HAVE_PREFIX],["$prefix"],[Look for where we prefix to])
	  
## Tells us where we are prefixed to
AC_MSG_NOTICE([Configuring for location "$prefix"])
## If we prefix away from the default /usr/ then create binary GLU
## will still try and create shared object libraries there
AM_CONDITIONAL([PREF],[ test "$prefix" == "NONE" ])
   if test "$prefix" == "NONE" ; then
       AC_MSG_NOTICE([Not creating binary "URFIT"])
   else
       AC_MSG_NOTICE([Creating binary "URFIT" in $prefix/bin])
       mkdir $prefix  
       mkdir $prefix/bin
   fi
AC_SUBST(PREF)		   

## C standard libraries
AC_CHECK_HEADERS([math.h],[
LDFLAGS="$LDFLAGS -lm" ])

## code needs gsl
gsl=false 
AC_ARG_WITH([gsl],
	AC_HELP_STRING([--with-gsl=<name>],
		       [Specify GSL location for the GF]),
            [GSL=${with_gsl}
	    ## test to see if something is there ...
	    if [ test -a ${GSL}include/gsl/gsl_rng.h ]; then
	    	AC_MSG_NOTICE([User defined Gnu Scientific Library (GSL) @ $GSL])
	  	## check if the gsl static lib is available
	    	if [ test -a ${GSL}lib/libgsl.a ]; then 
		   AC_MSG_NOTICE([Statically linking to the user's GSL lib])
	    	   CFLAGS="$CFLAGS -I${GSL}include"
	    	   LDFLAGS="$LDFLAGS ${GSL}lib/libgsl.a"	
		   AC_DEFINE([HAVE_GSL],,[Do we have the GSL libs?])
		elif [ test -a ${GSL}lib/libgsl.so ]; then
	    	   AC_MSG_NOTICE([Static link to GSL not found!! Shared object (.so) linking to the user's GSL link ])	
	    	   CFLAGS="$CFLAGS -I${GSL}include -L${GSL}lib/ "
	    	   LDFLAGS="$LDFLAGS -lgsl -lgslcblas"
		   AC_DEFINE([HAVE_GSL],,[Do we have the GSL libs?])
	    	else
	           AC_MSG_ERROR([gsl_rng.h found, but cannot find the static or shared libraries needed to link against!])
	    	fi
		echo " "
	    fi
	    ## have a look for the ones in /usr/
	    AC_CHECK_HEADERS([gsl/gsl_eigen.h],[LDFLAGS="$LDFLAGS -lgsl -lgslcblas"])
		],[ AC_MSG_ERROR([GSL needed (need --with-gsl) ]) ])
AM_CONDITIONAL([DEFGSL],[$gsl])
## Include where it is just in case people want to do something in Makefile.am
AC_SUBST(GSL)

## Look for fftw in the useer specified place
fftw=false 
AC_ARG_WITH([fftw],
	AC_HELP_STRING([--with-fftw=<name>],
		       [Specify FFTW location for the GF]),
            [
	    ## see if we already have it	
	    AC_CHECK_HEADERS([fftw3.h],[local_fftw=true],[local_fftw=false])
	    
	    FFTW=${with_fftw}

	    ## are we compiling for single precision or not?
	    if test x$single = xtrue ; then
	       float_version=f
	    else	
	       float_version=""
	    fi

	    AC_CHECK_HEADERS([${FFTW}include/fftw3.h],[
	    AC_MSG_NOTICE([User defined fftw library @ $FFTW])
	    AC_MSG_NOTICE([Attempting linking to the user's FFTW lib ${FFTW}lib/libfftw3${float_version}])
	
	    ## if we use the user one, only link to the static lib
	    if [ test -a ${FFTW}lib/libfftw3${float_version}.a ]; then
		AC_MSG_NOTICE([Statically linking to the user's FFTW lib])
	    	CFLAGS="$CFLAGS -I${FFTW}include"
	    	LDFLAGS="$LDFLAGS ${FFTW}lib/libfftw3${float_version}.a"
		if test x$OMP_FFTW = xtrue ; then
		   if test -a ${FFTW}lib/libfftw3${float_version}_omp.a ; then
		      LDFLAGS="$LDFLAGS ${FFTW}lib/libfftw3${float_version}_omp.a"
		      AC_DEFINE([OMP_FFTW],[],[OpenMP-parallel FFT routines.])
		   fi
		fi
	        AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    ## Have a look for the shared object library ...
            elif [ test -a ${FFTW}lib/libfftw3${float_version}.so ]; then
		AC_MSG_NOTICE([Static link not found !! Shared object (.so) linking to the user's FFTW lib])
	    	CFLAGS="$CFLAGS -I${FFTW}include -L${FFTW}lib/"
	    	LDFLAGS="$LDFLAGS -lfftw3${float_version}"
		if test x$OMP_FFTW = xtrue ; then
		   ## check to see if it is there otherwise don't include it
		   if test -a ${FFTW}lib/libfftw3${float_version}_omp.so ; then
		      LDFLAGS="$LDFLAGS -lfftw3${float_version}_omp"
		      AC_DEFINE([OMP_FFTW],[],[OpenMP-parallel FFT routines.])
		   fi
		fi
	        AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
            else
	    	AC_MSG_ERROR([fftw3.h found, but cannot find the shared or static libraries!])
	    fi
	    ]	    
	    , [
	    ## see if we already have it	
	    AC_CHECK_HEADERS([fftw3.h],[local_fftw=true],[local_fftw=false])

	    ## get desperate here and check for the header in the usual places
	    ## as omp routines are special I doubt they will be installed, maybe should do a special check
	    if test x$local_fftw = "xtrue" ; then
	       LDFLAGS="$LDFLAGS -lfftw3${float_version}"
	       AC_MSG_NOTICE([Using the system's fftw])	       
	       AC_DEFINE([HAVE_FFTW3_H],[],[FFTW library.])
	    else
	       AC_MSG_NOTICE([No local fftw3.h found either!])
	    fi
     ]) ## closes the check headers ...
]
## closes with-FFTW
)
AM_CONDITIONAL([DEFFFTW],[$fftw])
AC_SUBST(FFTW)

## Some compiler checks
## My code is littered with consts
AC_C_CONST([])

## Also use the restrict flag a lot
AC_C_RESTRICT([])

## Many of my functions are "inlined" {compiler's choice so maybe not} ,
## hence the large binary
AC_C_INLINE([])

## Set up config.h, one of the most important headers of them all
AC_CONFIG_HEADERS([config.h])

## where are my files?
AC_CONFIG_FILES([
Makefile
src/Makefile
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
