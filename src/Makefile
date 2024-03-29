# $Header:$
#mf 28-06-03 cleaned, deleted ESSL
#
# Makefile_linux (GNU make) for theorem@fhi
#
# pseudopotential generating program      make fhipp.x
# pseudopotential testing program         make pslp.x
# removing garbage                        make clean
#
# ---- begin settings ----

FC      = ifort
FFLAGS  = 

LD      = $(FC)

# cannot find static libraries on HPC and have errors "cannot find -la -lm", remove "-static"
#LDFLAGS = $(FFLAGS) -static
LDFLAGS = $(FFLAGS)

# (LAPACK) public domain routines, ESSL no longer supported
LIB     = LAPACK

# ---- end settings ----

ifeq ($(LIB), LAPACK)
	MODE    = LAPACK
	LIBS    = -L../lib -lFREE
	LIBS    = ../lib/libFREE.a 
endif

OBJ_GNCPP         = fhipp.o ncpp.o hamann.o tromar$(MODE).o dnlcc7$(MODE).o defrtm.o

OBJ_PSLP          = pslp.o psatom.o rcovalent.o ppcheck.o derlkb.o \
                    satf.o outlps.o inverse.o psatom_chen.o cost_chen.o \
					legendre_chen.o legendre_p_chen.o legendre_pp_chen.o psatom_etot.o \
					cgfam.o cgsearch.o lbfgsb.o linpack.o timer.o  # my additional codes 

OBJ_SPECTRAL_TEST = klbyii$(MODE).o laguerre.o gaussq.o 

OBJ_KINETIC_TEST  = kinkon.o dfbt_srt.o bessel.o

OBJ_SHARED        = outwfct.o labelmap.o ggarad.o gga91_sr.o bpgcxc.o ldaxc.o corlyp.o pbe.o \
                    logmesh.o atomini.o sratom_n.o vestat.o vexcor.o dftseq.o \
                    atoaux.o darraux.o anderson.o moment.o spinorbit.o relxc.o \
                    pberev.o stat_orb.o fxc_mgga_pk.o\
                    vklix.o invermat.o arhf.o acgc.o gaunt.o ecp.o cepvwn.o \
                    vexcos.o

OBJ_ESSL          = dgef.o dges.o dspev.o errsav.o errset.o einfo.o

all : fhipp.x pslp.x

OBJ_GNCPP_X = $(OBJ_GNCPP) $(OBJ_SHARED)

OBJ_PSLP_X  = $(OBJ_PSLP) $(OBJ_SHARED) $(OBJ_SPECTRAL_TEST) $(OBJ_KINETIC_TEST)

%.o : %.f
	$(FC) $(FFLAGS) -c $*.f

%.o : %.f90
	$(FC) $(FFLAGS) -c $*.f90

fhipp.x : $(OBJ_GNCPP_X)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

pslp.x : $(OBJ_PSLP_X)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

$(OBJ_SPECTRAL_TEST) : basis.h gauss.h

fhipp.o ncpp.o pslp.o dftseq.o logmesh.o : default.h

$(OBJ_PSLP) $(OBJ_GNCPP) : $(LIBS) parameter.h

$(LIBS) : 
ifeq ($(MODE),LAPACK)
	cd ../lib; $(MAKE) -f make.libFREE
endif

clean :
	rm -f *.o 

# contents
#
#   fhipp.f           main
#   ncpp.f            pseudopotentials
#   hamann.f          hamann generator
#   tromar.f          troullier-martins generator
#   dnlcc7.f          pseudocore generator
#   defrtm.f          default cutoff radii for troullier-martins 

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
#   cepvwn.f          Vosko/Wilk/Nussair LDA <inactive>
#   corlyp.f          Lee/Yang/Parr c
#   pbe.f             Burke/Perdew/Ernzerhof (PBE) GGA
#   pberev.f          revised PBE GGAs
#   acgc.f            (kli-x) for gaunt coefficients 
#   arhf.f            (kli-x) fock matrix elements
#   vklix.f           (kli-x) kli equation solver
#   fxc_mgga_pk.f     Perdew/Kurth/Blaha/Zupan meta GGA
#   stat_orb.f        kinetic energy densit for meta GGA
#   atomini.f         potential initialisation
#   logmesh.f         radial mesh initialisation
#   sratom_n.f        all-electron atom
#   vestat.f          Poisson & Hartree
#   vexcor.f          exchange-correlation 
#   vexcos.f          exchange-correlation spin-polarized
#   dftseq.f          radial deq solver
#   atoaux.f          predictor-corrector auxiliaries
#   darraux.f         array manipulation auxiliary
#   anderson.f        potential mixer
#   moment.f          radial moments etc.
#   spinorbit.f       spin-orbit couplings <inactive>

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
