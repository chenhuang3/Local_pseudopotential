#
# make auxiliary library for pseudopotential package
# sources extracted from LAPACK 
# see eg WWW: http://www.netlib.org/lapack/index.html
#

FC      = ifort
FC      = ifort
FFLAGS  = -c -O

LIB     = libFREE.a
DIR     = src.libFREE/

lib_obj = $(LIB)($(DIR)dasum.o)  $(LIB)($(DIR)daxpy.o)  $(LIB)($(DIR)dcopy.o) \
          $(LIB)($(DIR)dgemm.o)  $(LIB)($(DIR)dgemv.o)  $(LIB)($(DIR)dger.o) \
          $(LIB)($(DIR)dgesv.o)  $(LIB)($(DIR)dgetf2.o) $(LIB)($(DIR)dgetrf.o) \
          $(LIB)($(DIR)dgetrs.o) $(LIB)($(DIR)dlae2.o)  $(LIB)($(DIR)dlaebz.o) \
          $(LIB)($(DIR)dlaev2.o) $(LIB)($(DIR)dlagtf.o) $(LIB)($(DIR)dlagts.o) \
          $(LIB)($(DIR)dlamch.o) $(LIB)($(DIR)dlansp.o) $(LIB)($(DIR)dlanst.o) \
          $(LIB)($(DIR)dlapy2.o) $(LIB)($(DIR)dlarf.o)  $(LIB)($(DIR)dlarfg.o) \
          $(LIB)($(DIR)dlarnv.o) $(LIB)($(DIR)dlartg.o) $(LIB)($(DIR)dlaruv.o) \
          $(LIB)($(DIR)dlascl.o) $(LIB)($(DIR)dlaset.o) $(LIB)($(DIR)dlasr.o) \
          $(LIB)($(DIR)dlasrt.o) $(LIB)($(DIR)dlassq.o) $(LIB)($(DIR)dlaswp.o) \
          $(LIB)($(DIR)dnrm2.o)  $(LIB)($(DIR)dopgtr.o) $(LIB)($(DIR)dopmtr.o) \
          $(LIB)($(DIR)dorg2l.o) $(LIB)($(DIR)dorg2r.o) $(LIB)($(DIR)dscal.o) \
          $(LIB)($(DIR)dspevx.o) $(LIB)($(DIR)dspmv.o)  $(LIB)($(DIR)dspr2.o) \
          $(LIB)($(DIR)dsptrd.o) $(LIB)($(DIR)dstebz.o) $(LIB)($(DIR)dstein.o) \
          $(LIB)($(DIR)dsteqr.o) $(LIB)($(DIR)dsterf.o) $(LIB)($(DIR)dswap.o) \
          $(LIB)($(DIR)dtrsm.o)  $(LIB)($(DIR)idamax.o) $(LIB)($(DIR)ilaenv.o) \
          $(LIB)($(DIR)ddot.o)   $(LIB)($(DIR)lsame.o)  $(LIB)($(DIR)xerbla.o) \
          $(LIB)($(DIR)dlag2.o)

$(LIB): $(lib_obj)
	ranlib $(LIB)

%.a : %.f
	$(FC) $(FFLAGS) -o $% $*.f
	ar rv $@ $%
	rm $%
