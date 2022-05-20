# $Header:$
#
# makefile (GNU make) for
#
# pseudopotential generating program      make gncpp.x
# pseudopotential testing program         make pslp.x
# removing garbage                        make clean
#
# ---- begin settings ----

FC      = fort
FFLAGS  = -p -g -qxref=full 
FFLAGS  = 

LD      = $(FC)
LDFLAGS = $(FFLAGS) 

# make and use library of ESSL, or (LAPACK) public domain routines
LIB     = essl
LIB     = lapack

# ---- end settings ----

ifeq ($(LIB),essl)
	MODE    = ESSL
	LIBS    = -lessl
else
	MODE    = LAPACK
	LIBS    = ../lib/libFREE.a
#     LIBS    = -L../lib -lFREE
endif

OBJ_GNCPP         = gncpp.o ncpp.o hamann.o tromar$(MODE).o dnlcc7$(MODE).o defrtm.o 

OBJ_PSLP          = pslp.o psatom.o rcovalent.o ppcheck.o derlkb.o vexcos.o

OBJ_SPECTRAL_TEST = klbyii$(MODE).o laguerre.o gaussq.o 

OBJ_KINETIC_TEST  = kinkon.o dfbt_srt.o bessel.o

OBJ_SHARED        = outwfct.o labelmap.o ggarad.o gga91_sr.o bpgcxc.o ldaxc.o corlyp.o pbe.o \
                    logmesh.o atomini.o sratom_n.o vestat.o vexcor.o dftseq.o \
                    atoaux.o darraux.o anderson.o moment.o spinorbit.o relxc.o \
                    pbe-modify.o stat_orb.o fxc_mgga_pk.o\
                    vklix.o invermat.o arhf.o acgc.o gaunt.o ecp.o cepvwn.o overlap.o 

OBJ_ESSL          = dgef.o dges.o dspev.o errsav.o errset.o einfo.o

all : gncpp.x pslp.x

OBJ_GNCPP_X = $(OBJ_GNCPP) $(OBJ_SHARED)

OBJ_PSLP_X  = $(OBJ_PSLP) $(OBJ_SHARED) $(OBJ_SPECTRAL_TEST) $(OBJ_KINETIC_TEST)

%.o : %.f
	$(FC) $(FFLAGS) -c $*.f

gncpp.x : $(OBJ_GNCPP_X)
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

pslp.x : $(OBJ_PSLP_X)
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

$(OBJ_SPECTRAL_TEST) : basis.h gauss.h

gncpp.o ncpp.o pslp.o dftseq.o logmesh.o : default.h

$(OBJ_PSLP) $(OBJ_GNCPP) : $(LIBS) parameter.h

$(LIBS) : 
ifeq ($(MODE),LAPACK)
	cd ../lib; $(MAKE) -f make.libFREE
endif

clean :
	rm -f *.o 

# contents
#
#   gncpp.f           main
#   ncpp.f            pseudopotentials
#   hamann.f          hamann generator
#   tromar.f          troullier-martins generator
#   dnlcc7.f          pseudocore generator
#   defrtm.f          default cutoff radii for troullier-martins 
#   invermat.f        matrix invertor

#   pslp.f            main
#   psatom.f          pseudoatom
#   rcovalent.f       radii for logarithmic derivatives
#   ppcheck.f         log derivatives & ghost state analysis
#   derlkb.f          radial integrator for separable case

#   labelmap.f        labels for elements and XC functionals
#   outwfct.f         output subroutines
#   ggarad.f          gradient-corrected xc
#   gga91_sr.f        Perdew/Wang GGA91 xc
#   bpgcxc.f          Becke/Perdew xc
#   ldaxc.f           various LDA xc 
#   corlyp.f          Lee/Yang/Parr c
#   pbe.f             Burke/Perdew/Ernzerhof GGA
#   atomini.f         potential initialisation
#   logmesh.f         radial mesh initialisation
#   sratom_n.f        all-electron atom
#   vestat.f          Poisson & Hartree
#   vexcor.f          exchange-correlation 
#   dftseq.f          radial deq solver
#   atoaux.f          predictor-corrector auxiliaries
#   darraux.f         array manipulation auxiliary
#   anderson.f        potential mixer
#   moment.f          radial moments etc.

#   parameter.h       array dimensions
#   default.h         mesh & accuracy settings

#   klbyii.f          nonlocal hamiltonian in polynomial representation
#   laguerre.f        auxiliary: associated laguerre polynomials
#   gauss.f           gauss quadrature 
#   gauss.h           gauss quadrature weights and absissas
#   basis.h           scaling parameters for laguerre polynomials

#   kinkon.f          fourier analysis of wavefunctions
#   dfbt_srt.f        fourier bessel transform
#   bessel.f          bessel functions
