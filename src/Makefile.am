TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

lib_LIBRARIES = libURFIT.a

## just include the header directory
AM_CFLAGS=-I${TOPDIR}/src/HEADERS/

ANALYSIS_FILES=./ANALYSIS/alphas.c ./ANALYSIS/an_wrapper.c \
	./ANALYSIS/beta_crit.c ./ANALYSIS/correlator.c ./ANALYSIS/hvp_pade.c \
	./ANALYSIS/exceptional.c ./ANALYSIS/Qcorr.c ./ANALYSIS/sun_flow.c \
	./ANALYSIS/inverse.c ./ANALYSIS/KKops.c ./ANALYSIS/tetra_gevp.c \
	./ANALYSIS/KK_BK.c ./ANALYSIS/general_ops.c ./ANALYSIS/Ren_rats.c \
	./ANALYSIS/CPCV.c ./ANALYSIS/Wall_Local.c ./ANALYSIS/Qmoments.c \
	./ANALYSIS/statpot.c ./ANALYSIS/adler.c ./ANALYSIS/udcb.c \
	./ANALYSIS/nrqcd.c ./ANALYSIS/nrqcd_old.c ./ANALYSIS/nrqcd_slope.c \
	./ANALYSIS/nrqcd_baremass.c ./ANALYSIS/c4c7_analysis.c \
	./ANALYSIS/su2_shit.c ./ANALYSIS/binding_corr.c ./ANALYSIS/HAL.c \
	./ANALYSIS/pof.c ./ANALYSIS/HLBL.c ./ANALYSIS/ZV.c ./ANALYSIS/sol.c \
	./ANALYSIS/Qcorr_fixed.c ./ANALYSIS/SpinOrbit.c

EFFMASS_FILES=./EFFMASS/blackbox.c ./EFFMASS/effmass.c ./EFFMASS/gevp.c

INIT_FILES=./INIT/fake.c ./INIT/init.c

FITS_FILES=./FITS/alpha_D0.c ./FITS/alpha_D0_multi.c ./FITS/fit_chooser.c \
	./FITS/cosh.c ./FITS/exp_plusc.c ./FITS/exp.c \
	./FITS/HLBL_cont.c ./FITS/Nder.c \
	./FITS/poly.c ./FITS/pade.c ./FITS/pmap.c ./FITS/Pexp.c ./FITS/pp_aa.c \
	./FITS/pp_aa_ww.c ./FITS/ppaa.c ./FITS/pp_aa_ww_r2.c \
	./FITS/Qcorr_bessel.c ./FITS/sinh.c ./FITS/poles.c ./FITS/qslab.c \
	./FITS/qslab2.c ./FITS/exp_xinv.c \
	./FITS/cornell.c ./FITS/pp_aa_exp.c ./FITS/qsusc_su2.c \
	./FITS/adler_alpha_D0.c ./FITS/adler_alpha_D0_multi.c ./FITS/fvol1.c \
	./FITS/alpha_D0_multi_new.c ./FITS/udcb_heavy.c ./FITS/nrqcd_exp.c \
	./FITS/nrqcd_exp2.c  ./FITS/c4c7.c ./FITS/su2_shitfit.c \
	./FITS/SUN_cont.c ./FITS/HALexp.c ./FITS/cosh_plusc.c \
	./FITS/tanh.c ./FITS/ZV_exp.c ./FITS/fvol2.c ./FITS/fvol3.c ./FITS/fsol.c	\
	./FITS/fvolcc.c

GRAPH_FILES=./GRAPH/make_xmgrace.c ./GRAPH/plot_fitfunc.c

IO_FILES=./IO/distribution.c ./IO/GLU_bswap.c ./IO/io_wrapper.c \
	./IO/read_flat.c ./IO/read_corr.c ./IO/read_GLU.c \
	./IO/read_GLU_Qmoment.c \
	./IO/read_GLU_tcorr.c ./IO/tfold.c ./IO/write_flat.c

INPUT_FILES=./IO/INPUT/read_inputs.c ./IO/INPUT/read_traj.c \
	./IO/INPUT/read_fit.c ./IO/INPUT/read_graph.c \
	./IO/INPUT/read_stats.c

MINIMIZE_FILES=./MINIMIZE/CG.c ./MINIMIZE/GA.c ./MINIMIZE/GLS.c \
	./MINIMIZE/GLS_pade.c ./MINIMIZE/line_search.c \
	./MINIMIZE/LM.c ./MINIMIZE/SD.c

PHYSICS_FILES=./PHYSICS/cruel_runnings.c ./PHYSICS/decays.c ./PHYSICS/momenta.c\
	./PHYSICS/sort.c

RUN_FILES=./RUN/bootfit.c ./RUN/fit_and_plot.c ./RUN/gls_bootfit.c

STATS_FILES=./STATS/bootstrap.c ./STATS/jacknife.c ./STATS/stats.c \
	./STATS/resampled_ops.c ./STATS/correlation.c ./STATS/autocorr.c \
	./STATS/raw.c ./STATS/bin.c ./STATS/reweight.c

UTILS_FILES=./UTILS/chisq.c ./UTILS/crc32c.c ./UTILS/ffunction.c \
	./UTILS/gen_ders.c ./UTILS/histogram.c ./UTILS/Nint.c \
	./UTILS/NR.c ./UTILS/poly_coefficients.c \
	./UTILS/pade_coefficients.c ./UTILS/pade_laplace.c \
	./UTILS/rng.c ./UTILS/svd.c ./UTILS/summation.c

## all the source files apart from ./Run/Mainfile.c
libURFIT_a_SOURCES = \
	${ANALYSIS_FILES} \
	${EFFMASS_FILES} \
	${FITS_FILES} \
	${GRAPH_FILES} \
	${INIT_FILES} \
	${INPUT_FILES} \
	${IO_FILES} \
	${MINIMIZE_FILES} \
	${PHYSICS_FILES} \
	${RUN_FILES} \
	${STATS_FILES} \
	${UTILS_FILES}

## and their respective headers
include_HEADERS =
	$(ANALYSIS_FILES:./ANALYSIS/%.c=./HEADERS/%.h) \
	$(EFFMASS_FILES:./EFFMASS/%.c=./HEADERS/%.h) \
	$(FITS_FILES:./FITS/%.c=./HEADERS/%.h) \
	$(GRAPH_FILES:./GRAPH/%.c=./HEADERS/%.h) \
	$(INIT_FILES:./INIT/%.c=./HEADERS/%.h) \
	$(INPUT_FILES:./IO/INPUT/%c=./HEADERS/%.h) \
	$(IO_FILES:./IO/%.c=./HEADERS/%.h) \
	$(MINIMIZE_FILES:./MINIMIZE/%.c=./HEADERS/%.h) \
	$(PHYSICS_FILES:./PHYSICS/%.c=./HEADERS/%.h) \
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

