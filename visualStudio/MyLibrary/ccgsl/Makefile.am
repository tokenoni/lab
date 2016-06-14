## $Id: Makefile.am,v 1.25 2010/04/07 13:43:30 cms124 Exp $
## Process this file with automake to produce Makefile.in

ldflags = -version-info 0:1:0 

sources = 

headers = exception.hpp \
 complex.hpp \
 complex_float.hpp \
 complex_long_double.hpp \
 block.hpp \
 vector.hpp \
 matrix.hpp \
 block_float.hpp \
 vector_float.hpp \
 matrix_float.hpp \
 block_long_double.hpp \
 vector_long_double.hpp \
 matrix_long_double.hpp \
 block_int.hpp \
 vector_int.hpp \
 matrix_int.hpp \
 block_uint.hpp \
 vector_uint.hpp \
 matrix_uint.hpp \
 block_long.hpp \
 vector_long.hpp \
 matrix_long.hpp \
 block_ulong.hpp \
 vector_ulong.hpp \
 matrix_ulong.hpp \
 block_short.hpp \
 vector_short.hpp \
 matrix_short.hpp \
 block_ushort.hpp \
 vector_ushort.hpp \
 matrix_ushort.hpp \
 block_char.hpp \
 vector_char.hpp \
 matrix_char.hpp \
 block_uchar.hpp \
 vector_uchar.hpp \
 matrix_uchar.hpp \
 block_complex_float.hpp \
 vector_complex_float.hpp \
 matrix_complex_float.hpp \
 block_complex_long_double.hpp \
 vector_complex_long_double.hpp \
 matrix_complex_long_double.hpp \
 blas.hpp \
 permutation.hpp \
 permute.hpp \
 mode.hpp \
 linalg.hpp \
 sf.hpp \
 sf_result.hpp \
 sf.hpp \
 sf_airy.hpp \
 sf_bessel.hpp \
 sf_clausen.hpp \
 sf_coulomb.hpp \
 sf_coupling.hpp \
 sf_dawson.hpp \
 sf_debye.hpp \
 sf_dilog.hpp \
 sf_ellint.hpp \
 sf_elljac.hpp \
 sf_erf.hpp \
 sf_exp.hpp \
 sf_expint.hpp \
 sf_fermi_dirac.hpp \
 sf_gamma.hpp \
 sf_gegenbauer.hpp \
 sf_hyperg.hpp \
 sf_laguerre.hpp \
 sf_lambert.hpp \
 sf_legendre.hpp

locallibs = 

if WITH_DEBUG
AM_CXXFLAGS = -DDEBUG
endif
pkginclude_HEADERS = $(headers)

clean-local: clean-src
	rm -f *~

clean-src:
	cd $(srcdir); rm -f *~
