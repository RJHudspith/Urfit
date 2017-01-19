TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libURFIT.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

EFFMASS_FILES=./EFFMASS/blackbox.c ./EFFMASS/effmass.c ./EFFMASS/gevp.c

INIT_FILES=./INIT/fake.c ./INIT/init.c

FITS_FILES=./FITS/exp_plusc.c ./FITS/exp.c ./FITS/cosh.c ./FITS/sinh.c \
	./FITS/poly.c ./FITS/pade.c ./FITS/pmap.c

GRAPH_FILES=./GRAPH/make_xmgrace.c ./GRAPH/plot_fitfunc.c

IO_FILES=./IO/GLU_bswap.c ./IO/read_corr.c ./IO/tfold.c

INPUT_FILES=./IO/INPUT/read_inputs.c ./IO/INPUT/read_traj.c \
	./IO/INPUT/read_fit.c ./IO/INPUT/read_graph.c \
	./IO/INPUT/read_stats.c

RUN_FILES=./RUN/bootfit.c ./RUN/fit_chooser.c ./RUN/fit_and_plot.c

MINIMIZE_FILES=./MINIMIZE/CG.c ./MINIMIZE/GA.c ./MINIMIZE/line_search.c \
	./MINIMIZE/LM.c ./MINIMIZE/SD.c

STATS_FILES=./STATS/bootstrap.c ./STATS/jacknife.c ./STATS/stats.c \
	./STATS/resampled_ops.c ./STATS/correlation.c ./STATS/autocorr.c

UTILS_FILES=./UTILS/chisq.c ./UTILS/ffunction.c ./UTILS/rng.c ./UTILS/svd.c \
	./UTILS/poly_coefficients.c ./UTILS/pade_coefficients.c \
	./UTILS/summation.c

## all the source files apart from ./Run/Mainfile.c
libURFIT_a_SOURCES = \
	${EFFMASS_FILES} \
	${FITS_FILES} \
	${GRAPH_FILES} \
	${INIT_FILES} \
	${INPUT_FILES} \
	${IO_FILES} \
	${MINIMIZE_FILES} \
	${RUN_FILES} \
	${STATS_FILES} \
	${UTILS_FILES}

## and their respective headers
include_HEADERS = 
	$(EFFMASS_FILES:./EFFMASS/%.c=./HEADERS/%.h) \
	$(FITS_FILES:./FITS/%.c=./HEADERS/%.h) \
	$(GRAPH_FILES:./GRAPH/%.c=./HEADERS/%.h) \
	$(INIT_FILES:./INIT/%.c=./HEADERS/%.h) \
	$(INPUT_FILES:./IO/INPUT/%c=./HEADERS/%.h) \
	$(IO_FILES:./IO/%.c=./HEADERS/%.h) \
	$(MINIMIZE_FILES:./MINIMIZE/%.c=./HEADERS/%.h) \
	$(RUN_FILES:./RUN/%.c=./RUN/%.h) \
	$(STATS_FILES:./STATS/%.c=./HEADERS/%.h) \
	$(UTILS_FILES:./UTILS/%.c=./HEADERS/%.h) \
	./HEADERS/gens.h \
	./HEADERS/fitfunc.h \
	./HEADERS/fits.h \
	../config.h

if !PREF

Bindir = "${prefix}"/bin

bin_PROGRAMS = URFIT

URFIT_SOURCES = Mainfile.c
URFIT_CFLAGS = ${CFLAGS} -I${TOPDIR}/src/HEADERS/
URFIT_LDADD = libURFIT.a ${LDFLAGS}

endif
