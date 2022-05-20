c Header:$
cmf 27-06-03 added spin polarization for selected xc functionals
cmf 04-09-02 added kli exchange functional
c
      program pslp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c FHI pseudopotential tool's
c
c              beta version as of 27-06-2003
c
c Except for some cleanup and fine tuning this may be about the final
c version. If you encounter bugs or have some comments, I would like to
c to hear about them. Of course you may use this program freely, but
c please do not redistribute it. Until the final version is released
c I would appreciate if interested parties obtained their copies from
c me instead. Please give proper credit as indicated below.
c
c Martin Fuchs     fuchs@fhi-berlin.mpg.de
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------
c Pseudopotential analyzing routine for semilocal pseudopotentials
c - self-consistent electronic structure for pseudo atom
c - local density/ generalized gradient approximation,
c   some useable for spin-polarized calculation
c - kinetic energy-density dependent functionals (experimental)
c - exact kohn-sham exchange within kli approximation (experimental)
c - linearized/ nonlinear core-valence exchange-correlation
c - valuation of separable Kleinman-Bylander (KB) representation 
c   * calculation of logarithmic derivatives
c   * search for ghost states (cf Ref [1])
c   * bound state spectrum w/ basis of orthogonal polynomials
c    (a crude method for any excited state, yet sufficient here)
c
c Familiarity with the references below should be regarded as a 
c prerequisite to a profitable use of this program. 
c Please cite these references when publishing results obtained with
c pseudopotentials constructed with the help of this program.
c
c [1] X. Gonze, R. Stumpf, M. Scheffler, Phys. Rev. B 44, 8503 (1991)
c     "Analysis of separable potentials"
c
c [2] M. Fuchs, M. Scheffler, Comput. Phys. Commun. 119, 67-98 (1999)
c
c Martin Fuchs
c Fritz-Haber-Institut der Max-Planck-Gesellschaft
c Abteilung Theorie
c Faradayweg 4-6
c D-14195 Berlin - Dahlem
c Germany
c E-mail: fuchs@fhi-berlin.mpg.de
c-----------------------------------------------------------------------
c
c input files
c
c fort.21    control parameters 
c fort.22    atomic parameters 
c fort.31    semilocal ionic pseudopotential (fhi94md format)
c fort.37    all-electron (reference) potential 
c
c output
c
c log. derivatives
c fort.6[l = 0,1,2...]	all-electron case
c fort.5[l = 0,1,2...]	semilocal case
c fort.7[l = 0,1,2...]	separable case (sbrt. klbyan)
c
c fort(iu)    monitoring data, ghost state analysis (-> parameter.h)
c fort.27     modelcore density
c fort.38     pseudo wavefunctions
c fort.4[0-4] ionic pseudopotentials
c fort.4[5-9] screened pseudopotentials
c fort.80     info file for kli case
c-----------------------------------------------------------------------
c
c input description fort.21 and fort.22
c
c fort.21 - line 1
c l_loc  angular momentum for local potential
c lbeg   dto. lowest  for evaluating logarithmic derivatives
c lend   dto. highest for evaluating logarithmic derivatives
c lmax   dto. maximum represented by pseudopotentials
c rld    diagnostic radius for evaluating logarithmic derivatives
c
c fort.21 - line 2
c tlgd   .true.  compute logarithmic derivatives
c tkb    .true.  perform ghost state analysis
c tiwf   .true.  use input wavefunctions to make nonlocal potentials
c tnrl   .true.  assume pseudopotentials are nonrelativistic
c        default is .false. (scalar-relativistic pseudopotentials)
c
c fort.21 - line 3
c tcut_global    .true.  before calculating anything, enforce a
c                 Coulomb tail in pseudopotential for radii beyond
c rcut_global    global cutoff radius. Default is .false. . The
c                 modified pseudopotential is written out.
c
c fort.21 - line 4
c tspin  .true.  perform a spin-polarized calculation.
c        Default is .false. . Implemented for selected xc functionals
c        only. Requires modified input file fort.22 .
c
c fort.22 - line 1
c zfull  atomic number
c nc     number of core states 
c nv     number of valence states
c iexc   exchange-correlation functional
c        an S indicates valid spin-polarized functionals
c         1  LDA   Wigner
c         2  LDA   Hedin/Lundqvist
c         3  LDA   Ceperley/Alder Perdew/Zunger (1980)           
c         4 SGGA   Perdew/Wang (1991)
c         5 SGGA   Becke (1988) X, Perdew (1986) C
c         6 SGGA   Perdew/Burke/Ernzerhof (1996)
c         7  LDA   Zhao/Parr
c         8 SLDA   Ceperley/Alder Perdew/Wang (1991)
c         9  GGA   Becke (1988) X, Lee/Yang/Parr (1988) C
c        10  GGA   Perdew/Wang (1991) X, Lee/Yang/Parr (1988) C
c        11  LDA   exchange only
c        12  KLI   excact Kohn-Sham exchange in the Krieger et al.
c                  approximation
c        13  KLI   KLI+effective core potential
c        14  GGA   Hammer/Norskov revised PBE GGA (RPBE)
c        15  GGA   Zhang/Wang revised PBE GGA (revPBE)
c        16  MGGA  Perdew/Kurth/Zupan/Blaha MetaGGA (1999)
c        17  GGA   GGA exchange + LDA correlation, PBE X, LDA C
c        18  KLI   exact exchange + LDA correlation
c        19  KLI   exact exchange + PBE GGA correlation
c        default is 8
c rnlc   pseudocore radius 
c
c fort.22 - line 2 etc.
c n(i)   principal quantum number
c l(i)   angular momentum
c f(i)   occupation number (two entries for spin-polarized mode)
c
c sample input fort.21 for silicon (3s3p3d potentials)
c > 2 0 2 2 0.0     : l_loc lbeg lend lmax rld
c >.t. .t. .t. .f.  : tlgd  tkb  tiwf tnrl
c
c sample input fort.22 for silicon (LDA, linearized CV XC)
c > 14.0 3 2 8 0.0  : z nc nv iexc rnlc
c > 1   0  2.0      : n l f
c > 2   0  2.0
c > 2   1  6.0
c > 3   0  2.0
c > 3   1  2.0
c excess lines are not read
c
c sample input fort.22 for silicon (LSDA, linearized CV XC)
c > 14.0 3 2 8 0.0  : z nc nv iexc rnlc
c > 1   0  1.0 1.0  : n l f
c > 2   0  1.0 1.0
c > 2   1  3.0 3.0
c > 3   0  1.0 1.0
c > 3   1  2.0 0.0
c excess lines are not read
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical tkb,tins,terr,tnlc,tlgd,tnrl,tiwf,tcut_global,tspin
      character*8 symz,symxc*40,srev*10
      character*1 sn,sl,sm

      include 'parameter.h'
      include 'default.h'

      parameter(rmx=80.d0)

c &&&&&&&&&&&&& REVISION &&&&&&&&&&&&&
      parameter(srev='rev270603B')

      integer   ir,iu,mmaxtmp
      integer   np(ms),lp(ms),n(ms),m(ms),l(ms),ntl(ms),ninu(ms)
      real*8    f(ms),fp(ms),e(ms),rdv(100),vps(mx,ms)
      real*8    r(mx),vi(mx),ves(mx),vxc(mx,ms),vee(mx,ms),ups(mx,ms)
      real*8    vae(mx),ep(ms),igrmap(20),dummy(ms),uout(mx,ms)
      real*8    uu(mx,ms),vbare(mx,ms),vorb(mx,ms),vsl(mx)
      real*8    dc(mx),dcp(mx),dcpp(mx),rho(mx),rhop(mx),rhopp(mx)
      real*8    rhos(mx,ms),rhosp(mx,ms),rhospp(mx,ms)
      real*8    rxx(mx),vvo(mx),rc(ms)
      real*8    x_l(mx),vbare_l(mx,ms)

      real*8  dexc
      common/xced/dexc(mx)

c kli
      logical tkli
      real*8  svm(ms),svm_roks(ms),svkli(ms)
c meta gga
      logical tmgga
      real*8  stat,stat_core
      common/use_mgga/stat(mx),stat_core(mx)

c===========================
c chen: vlocal related 
      parameter (npdim=mx)
      integer vtype,npar,icut_end(ms),icut,icut2,
     &        iorb_start,norb2,mag_fit
      logical scf_mode,test_config,five_point
      real*8  param_new(npdim),param(npdim),param1(npdim),wfopt(ms)
      real*8  vbare_v(mx),vks(mx),vloc_new(mx),vloc(mx),vlocs(mx,ms),
     &        zncpp,zval,rcut_coul,r_cut,lam,fac,eref(ms),
     &        step,pen,pen1,fcost,fcost1,fcost2,wf_match,
     &        vtmp(mx,ms),ups_ref(mx,ms),ndiff,pcoeff,
     &        h0,h1,f0,f1,f2,rat,vs,vsp,vspp,
     &        norm_tmp,norm0(ms),norm(ms),
     &        rho_ref(mx),rhop_ref(mx),rhopp_ref(mx),
     &        occ2(mx),vorb_nl(mx,ms),
     &        mag_energy,mag_weight,occ_sp(ms,2),mag_de

      ! lbfgs 
      character*60 task,csave
      integer isave(44),iprint,lbfgs_m
      logical lsave(4)
      real*8  dsave(29),lbfgs_factr,pgtol
      real*8,allocatable :: lbfgs_wa(:),grad(:),lbfgs_l(:),lbfgs_u(:)
      integer,allocatable :: lbfgs_iwa(:),lbfgs_nbd(:)
c===========================

      external  dmelm,gltfmv,fmom,rcovalent

      data igrmap/0,0,0,1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,1,0/
c
c begin 
c output unit
      iu=iopslp
c defaults
      iexc=8
      lmax=2
      lbeg=0
      lend=lmax
      l_loc=lmax
      rld=0.d0
c do log derivatives tlgd=.true. 
      irl=1
      tcut_global=.false.
      rcut_global = 1e3
      tiwf=.true.
      tkb=.true.
      tkli=.false.
      tlgd=.true.
      tnlc=.false.
      tnrl=.false.
      tmgga=.false.
      tspin=.false.
c log file
      write(iu,619) srev
c initialize
      do i=1,ms
        fp(i)=0.d0
        f(i)=0.d0
        f(i)=0.d0
        n(i)=0
        l(i)=0
      enddo
      ispinmx=1

c header file
      read(21,*,err=10,end=10) l_loc,lbeg,lend,lmax,rld
   10 read(21,*,err=11,end=11) tlgd,tkb,tiwf,tnrl
   11 read(21,*,err=12,end=12) tcut_global,rcut_global
   12 read(21,*,err=13,end=13) tspin
   13 if(tspin) ispinmx=2 
c debug
!d     write(ie,*) '%pslp - ispin,tspin=',ispinmx,tspin

      read(22,*,err=14,end=14) zfull,nc,nv,iexc,rnlc
      zfull=abs(zfull)
      nc=abs(nc)
      nv=abs(nv)
      imx=nc+nv
      do i=1,imx
        if(ispinmx.eq.1) then
          read(22,*,err=14,end=14) np(i),lp(i),fp(i)
        else
          read(22,*,err=14,end=14) np(i),lp(i),fp(i),fp(i+imx)
        endif
        fp(i)=abs(fp(i))
        fp(i+imx)=abs(fp(i+imx))
        if(np(i) .le. lp(i)) stop 'pslp - data: l .le. n'
        if(lp(i) .gt. 3)     stop 'pslp - data: l .gt. 3'
        npmax=max(np(i),npmax)
      enddo
   14 imx=i-1
      if(ispinmx.gt.1) then
        do i=1,imx
c for occupation of spin states ...
          if(fp(i+imx).eq.0.d0) then
c ... follow 2nd Hund rule
            fup=min(fp(i),2.d0*lp(i)+1.d0)
            fdn=fp(i)-fup
          else
            if(fp(i).gt.fp(i+imx)) then
c ... follow fixed majority spin
              fup=fp(i+imx)
              fdn=fp(i)-fup
            else
c ... use occupancies of input
              fup=fp(i)
              fdn=fp(i+imx)
            endif
          endif
          fp(i)=fup
          np(i+imx)=np(i)
          lp(i+imx)=lp(i)
          fp(i+imx)=fdn
        enddo
      endif  
c debug
!d     write(ie,*) '%pslp - state count: imx=', imx,' of',(nc+nv)*ispinmx
!d     do i=1,imx
!d       write(ie,*) np(i),lp(i),fp(i),np(i+imx),lp(i+imx),fp(i+imx)
!d     enddo

      if(tnrl) irl=2
      if(abs(rnlc) .gt. 0.0) tnlc=.true.
      iexc=abs(iexc)
      if(iexc.eq.12 .or. iexc.eq.13 .or. iexc.eq.18 .or.
     1   iexc.eq.19 .or. iexc.eq.16) tkli=.true.
      igr=igrmap(iexc)
      igr=1
      if(iexc .eq. 16) tmgga = .true. ! meta gga
c map for spin polarized calculation
      if(ispinmx.gt.1) then
        iexcs=3
        if(iexc.eq.3 .or. iexc.eq.8) then
          iexcs=3                !LDA
        else if(iexc.eq.5) then
          iexcs=5                !BP GGA
        else if(iexc.eq.4) then  
          iexcs=6                !PW91 GGA
        else if(iexc.eq.6) then
          iexcs=9                !PBE GGA
        else
          write(ie,'(a/,a,i3,a)')
     1      '%pslp - ERROR: spin polarized calculation not implemented',
     1      ' for input iexc=',iexc,'. ACTION: set iexc=3,4,5,6,8.'
          stop
        endif
      endif

c open files
      if(tkli) open(unit=80,file='kliinfoPS.dat',status='unknown')

      lmax=max(0,lmax)
      l_loc=max(0,min(l_loc,lmax))
      lbeg=max(0,lbeg)
      lend=max(lbeg,min(lend,5))
      llmx=lmax+1
      llbeg=lbeg+1
      llend=lend+1
      ll_loc=l_loc+1

      do ll=1,llmx
        n(ll)=ll
        l(ll)=ll-1
      enddo

      cval=0.d0
      norb=0
      do ispin=1,ispinmx
        norbstart=norb+1
        istart=nc*ispin+nv*(ispin-1)+1
        iend=istart+nv-1
c debug
!d       write(ie,*) '%pslp - istart=',istart,' iend=',iend
        do i=istart,iend
          norb=norb+1
          n(norb)=lp(i)+1
          do j=norbstart,norb-1
            if(l(j) .eq. lp(i)) n(norb)=n(norb)+1
          enddo
          l(norb)=lp(i)
          f(norb)=fp(i)
          m(norb)=ispin
          cval=cval+fp(i)
c debug
!d         write(ie,*) 'norb,n,l,m,f',
!d    1      norb,n(norb),l(norb),m(norb),f(norb)
        enddo
      enddo

c pseudopotentials, radial meshs assumed to be the same throughout !!!
c format same as for fhi94md 
      mmax=mx
      do ll=1,llmx
        terr=.true.
        read(31,*,end=58,err=58) mmax,amesh
c fort.31: semilocal ionic pseudopotential (fhi94md format)
        do i=1,mmax
          read(31,*,end=58,err=58) iunused,r(i),uu(i,ll),vbare(i,ll)
        enddo
        terr=.false.
   58   if(i .le. mmax) terr=.true.
        if(abs(amesh-r(2)/r(1)) .gt. 1e-10) terr=.true.
        if(terr) then
          write(ie,*) '& pslp - stop: bad pseudopotential file fort.31'
          write(ie,*) '&        current l,mmax,i ',ll-1,mmax,i
          write(ie,*) '& amesh as input          ',amesh
          write(ie,*) '& amesh set to r(2)/r(1)  ',r(2)/r(1) 
          amesh = r(2)/r(1)
        endif
      enddo
      do i=1,mmax
        read(31,*,end=60,err=60)  r(i),dc(i),dcp(i),dcpp(i)
      enddo
   60 if(tnlc .and. i .eq. 1) then 
        write(ie,'(a/,a)')
     1    '%pslp - ERROR: bad pseudopotential file fort.31, cannot',
     1    ' read core density. ACTION: check pseudopotential file.'
        stop
      endif
      if(.not. tnlc .and. i .gt. 1) then
        write(ie,*) '% pslp - pseudocore present (fort.31)'
        tnlc=.true.
      endif
      if(tnlc.and.tkli) then
          write(ie,'(a/,a,i3,a)')
     1      '%pslp - ERROR: core density not implemented',
     1      ' for iexc=',iexc,'. ACTION: set iexc=3,4,5,6,8.'
          stop
      endif

c meta gga read in core kinetic energy density
      do i=1,mmax
        stat(i) = 0.d0
        stat_core(i) = 0.d0
      enddo
      if(tnlc .and. tmgga) then
        do i=1,mmax
          read(17,*,end=81,err=81) rxx(i),stat_core(i)
        enddo
   81   if(i .eq. 1) then
           write(ie,*) '& pslp - warning: check length of file fort.17'
           write(ie,*) '&        core kinetic energy density'
        endif
        write(ie,*) '& pslp - warning: t_core = t_core_vW'
        do i=1,mmax
          stat_core(i) = 0.125*dcp(i)**2/max(1.d-12,dc(i))
        enddo
      endif  

c restrict range of radial mesh
      amesh=r(2)/r(1)
      al=log(amesh)
      do i=mmax,1,-1
        if(r(i) .lt. rmx) goto 69
      enddo
  69  mmax=i
      if(mod(mmax,2) .eq. 0) mmax=mmax-1

c apply global cutoff to potentials (set to coulomb tail)
      if(tcut_global) then
        write(ie,'(a,/,a,f8.4,a)') 
     1    ' %pslp - WARNING: the input pseudopotential is modified',
     1    '  by imposing a global cutoff at radius = ',rcut_global,'.'  
        do ll=1,llmx
          write(ie,'(a,1x,i1)') 
     1    '        .cutoff imposed for l =',ll-1
          do i=mmax,1,-1
            v_coulomb = -cval/r(i)
            if(r(i) .lt. rcut_global) exit
            vbare(i,ll)=v_coulomb+exp(-3.5*(r(i)-rcut_global)**2)
     1              *(vbare(i,ll)-v_coulomb)
          enddo
        enddo
        write(ie,*)
      endif

c self-consistent pseudoatom 
c initial eigenvalues
      iend=norb/ispinmx
      do i=1,iend
        fp(i)=f(i)
        if(tspin) fp(i)=fp(i)+f(i+iend)
c debug
!d       write(ie,*) '%pslp - i,fp',i,fp(i)
      enddo
      call atomini(iend,mmax,cval,npmax,fp,r,vi,e)
      call dnuv(mx,1,mmax,0.d0,vi)

      do i=1,norb
        e(i)=-0.5d0*max(1.d0,sqrt(cval))
        do ll=1,llmx
          if(l(i)+1 .eq. ll ) then
            call dcpv(mx,1,mmax,vbare(1,ll),vorb(1,i))
          else if(l(i)+1 .gt. llmx ) then
            call dcpv(mx,1,mmax,vbare(1,ll_loc),vorb(1,i))
          endif
        enddo
      enddo
      do i=1,mmax  !new
c       vi(i)=-max(1.d0,cval)/r(i)
        vi(i)=-max(1.d0,sqrt(cval))/r(i)
        do ispin=1,ispinmx
          vee(i,ispin)=vi(i)
        enddo
      enddo !new

c tell psatom when to do spin polarized calculation
      iexcnow=iexc
      if(ispinmx.gt.1) iexcnow=-iexcs 






c================================================
c chen: build vlocal psp 
c================================================
      if (nspinmx>1) then 
        write(iu,*)'code does not work for spin polarized case. stop'
        stop 
      endif 
      print *,''
      print *,' ***** use fine mesh in default.h *****'
      print *,''

      ! read KS potential from psgen run
      print *,' *** read total system KS potential from vks.dat *** '
      open(file='vks.dat',unit=111,action='read',form='formatted')
      do i=1,mmax
        read(111,*)dtmp,vks(i)
      enddo
      close(111)

      ! read unscreened valence AE potential (written by fhi98.f)
      open(file='vbare_valence.dat',action='read',
     &     unit=111,form='formatted')
      do i=1,mmax
        read(111,*)dtmp,vbare_v(i)
      enddo
      close(111)

      open(file='sys.dat',unit=111,action='read',form='formatted')
      read(111,*)z     ! valence electron number 
      read(111,*)zval  ! valence electron number 
      close(111)

      !-------------------------
      ! input file: vlocal.ini
      !-------------------------
      open(file='vlocal.ini',unit=111,action='read',form='formatted')
      read(111,*)i
      if (i>0) mk_vlocal=.true.
      if (i<0) mk_vlocal=.false.
      read(111,*)rcut 
      read(111,*)vtype ! 1: cutoff function (Starkloff and Joannopoulos PRB 1977, 16, 5212)
                       ! 2: polynomial 
                       ! 3: polynomial (consider linear term)
                       ! 4: vlocal on grid 
                       ! 5: polynomial (no linear term, x=r/rc)
                       ! 6: polynomial (no linear term, legendre)
      read(111,*)iorb_start ! (iorb_start .. norb) orbitals are matched for eigenvalues
      read(111,*)pcoeff        ! penalty coefficient 
      read(111,*)wfopt(1:norb) ! flag for optimizing each orbital
      read(111,*)mag_fit    ! fit spin polarization? 
      read(111,*)mag_weight ! weight for the magnetic fitting 
      read(111,*)mag_energy ! high spin (energy) - low spin (energy)
      read(111,*)occ_sp(1:nv,1) ! occ number for spin up for magnetic fitting
      read(111,*)occ_sp(1:nv,2) ! occ number for spin down for magnetic fitting 
      read(111,*)npar   

      if (vtype/=4) then 
        do i=1,npar       
          read(111,*)param(i) 
        enddo 
      endif 
      close(111)

      ! check parameters 
      if (iorb_start<1 .or. iorb_start>norb) then 
        print *,'incorrect value for iorb_start, iorb_start:',iorb_start
        stop
      endif 
      
      scf_mode=.true. ! true: compute vlocal 
                      ! false: compute vks 
 
      if (mk_vlocal) then 

      !====================================================
      ! find cutoff position, compute dv/dx and dv^2/dx^2
      !====================================================
      do i=1,mmax-1
        if (r(i)<rcut.and.rcut<r(i+1)) then 
          icut=i
          if (scf_mode) then 
            ! get gradient and curvature of vbare_v at r(icut)
            h0=r(i)-r(i-1)
            h1=r(i+1)-r(i)
            f0=vbare_v(i-1)
            f1=vbare_v(i)
            f2=vbare_v(i+1)
            rat = h1/h0

            vs = vbare_v(icut)
            ! using f(i-1),f(i),f(i+1) [five-point stencil]
            vsp = (f2-rat**2*f0-(1.d0-rat**2)*f1)/(h1*(1.d0+rat)) ! dv/dx
            vspp = 2.d0*(f2+rat*f0-(1.d0+rat)*f1)/(h0*h1*(1.d0+rat)) ! dv^2/dx^2
          else
            ! get gradient and curvature of vks at r(icut)
            h0=r(i)-r(i-1)
            h1=r(i+1)-r(i)
            f0=vks(i-1)
            f1=vks(i)
            f2=vks(i+1)
            rat = h1/h0

            vs = vks(icut)
            ! using f(i-1),f(i),f(i+1)
            vsp = (f2-rat**2*f0-(1.d0-rat**2)*f1)/(h1*(1.d0+rat)) ! dv/dx
            vspp = 2.d0*(f2+rat*f0-(1.d0+rat)*f1)/(h0*h1*(1.d0+rat)) ! dv^2/dx^2
          endif 
          exit
        endif 
      enddo
      print *,'rcut:      ',rcut
      print *,'icut:      ',icut
      print *,'pcoeff:    ',pcoeff
      print *,'iorb_start:',iorb_start
      print *,'mag_fit:    ',mag_fit 
      print *,'mag_weight: ',mag_weight
      print *,'mag_energy: ',mag_energy
      print *,'occ_sp (up):',occ_sp(1:nv,1)
      print *,'occ_sp (dn):',occ_sp(1:nv,2)
      if (vtype==6) 
     &  print *,"=== legendre polynomial, dv/dx(x=0)=0 ==="
      if (vtype==7) 
     &  print *,"=== legendre polynomial, v'(0)=0 & v''(0)=0 ==="

      ! vlocal is on mesh 
      if (vtype==4) then 
        npar=icut 
        if (scf_mode) then 
          param=vbare_v
        else
          param=vks
        endif 
      endif 

      ! load AE wfval.dat (written by fhipp.f)
      ups_ref = 0.d0 
      open(file='wfval.dat',unit=111,action='read',form='unformatted')
      read(111)ups_ref(:,1:nv)
      close(111)
      print *,' *** loaded AE-valence wf from wfval.dat ***'

      ! load AE valence rho (written by fhipp.f)
      open(file='rhoval.dat',unit=111,action='read',form='unformatted')
      read(111)rho_ref
      read(111)rhop_ref
      read(111)rhopp_ref
      close(111)
      print *,' *** loaded AE-valence density from rhoval.dat ***'

      ! find the vanishing r for each orbital by tesing norm of each orbital 
      do i=1,norb
        do icut2=1,mmax-1
          norm_tmp=0.d0 
          do j=1,icut2
            norm_tmp=norm_tmp+al*r(j)*ups_ref(j,i)**2
          enddo
          norm_tmp = norm_tmp
     &            + al*(23.0d0*r(icut2-2)*ups_ref(icut2-2,i)**2
     &            + 28.0d0*r(icut2-1)*ups_ref(icut2-1,i)**2
     &            +  9.0d0*r(icut2)*ups_ref(icut2,i)**2)/24.0d0

          if (abs(norm_tmp-1.d0)<1e-8) then 
            icut_end(i)=icut2
            write(6,'(a,i3,a,f10.4)')
     &      ' orbital ',i,' vanish at ',r(icut2)
            exit 
          endif 
        enddo
      enddo
      
      ! compute norm for r<rcut
      norm0=0.d0
      do i=1,norb
        ! integeration over [0,r_cut]
        ndim=icut
        do j=1,ndim-3
          norm0(i)=norm0(i)+al*r(j)*ups_ref(j,i)**2
        enddo
        norm0(i)=norm0(i) + al*(23.0d0*r(ndim-2)*ups_ref(ndim-2,i)**2
     &              + 28.0d0*r(ndim-1)*ups_ref(ndim-1,i)**2
     &              +  9.0d0*r(ndim)*ups_ref(ndim,i)**2)/24.0d0
      enddo

      ! load eigval.dat 
      open(file='eigval.dat',unit=111,action='read')
      do i=1,nv
        read(111,*)eref(i)
      enddo
      close(111)

      write(6,'(/,a)')'---- reference eigenvalues and wavefunction ----'
      do i=1,norb
        write(6,'(a,2i2,a,f7.2,a,f8.4,a,f9.3,a,a,es10.3,a,f7.4)') 
     &  'n,l:',n(i),l(i),' => occ:',f(i), 
     &  ' eref:',eref(i),' (',eref(i)*ry2,' eV)',
     &  ' u_ref(icut): ',ups_ref(icut,i),' norm(r<rc):',norm0(i)
      enddo
      print *,''
      print *,' --- orbital opt. weights (wfopt) ---'
      do i=1,norb
        write(6,'(a,i3,a,es9.2)')'    orb',i,': ',wfopt(i) 
      enddo
      print *,''

      !========================================
      ! optimize vloc to reproduce eigenvalue 
      !========================================
      task='START' 
      lbfgs_factr=0.d0
      pgtol=1e-16
      iprint=90
      lbfgs_m=10
      allocate(lbfgs_wa(2*lbfgs_m*npar+5*npar+11*lbfgs_m**2+8*lbfgs_m))
      allocate(lbfgs_nbd(npar),lbfgs_l(npar),lbfgs_u(npar))
      allocate(grad(npar),lbfgs_iwa(3*npar))
      lbfgs_nbd = 0
      iter=0
      
      do while (.true.)   
        call setulb(npar,lbfgs_m,param(1:npar), 
     &              lbfgs_l,lbfgs_u,lbfgs_nbd,
     &              fcost,grad(1:npar),lbfgs_factr,
     &              pgtol,lbfgs_wa,lbfgs_iwa,task,iprint,
     &              csave,lsave,isave,dsave)
        
        if (task(1:2) .eq. 'FG') then
          ! compute cost
          call cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,
     &     z,r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &     rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &     icut,eref,ups_ref,npdim,npar,param,iorb_start,scf_mode,zval,
     &     mag_fit,mag_weight,mag_energy,occ_sp,mag_de,
     &     vtype,0,pcoeff,pen,fcost,wf_match,vloc)
          
          if (vtype/=6) then 
            five_point=.true.  
            step=1e-4
          endif
          if (vtype==6) then
            five_point=.true.
            step=1e-4  ! Legendre polynomial 
          endif 

          ! compute gradient (finite diff)
          do i=1,npar
            if (five_point) then 
              param1 = param
              param1(i) = param(i)+2.d0*step
              call cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,
     &         z,r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &         rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &         icut,eref,ups_ref,npdim,npar,param1,iorb_start,
     &         scf_mode,zval,
     &         mag_fit,mag_weight,mag_energy,occ_sp,mag_de,
     &         vtype,0,pcoeff,pen1,fcostA,wf_match,vtmp)
            endif 

            param1 = param
            param1(i) = param(i)+step
            call cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,
     &       z,r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &       rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &       icut,eref,ups_ref,npdim,npar,param1,iorb_start,
     &       scf_mode,zval,
     &       mag_fit,mag_weight,mag_energy,occ_sp,mag_de,
     &       vtype,0,pcoeff,pen1,fcost1,wf_match,vtmp)

            param1 = param
            param1(i) = param(i)-step
            call cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,
     &       z,r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &       rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &       icut,eref,ups_ref,npdim,npar,param1,iorb_start,
     &       scf_mode,zval, 
     &       mag_fit,mag_weight,mag_energy,occ_sp,mag_de,
     &       vtype,0,pcoeff,pen1,fcost2,wf_match,vtmp)
            
            if (five_point) then 
              param1 = param
              param1(i) = param(i)-2.d0*step
              call cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,
     &         z,r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &         rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &         icut,eref,ups_ref,npdim,npar,param1,iorb_start,
     &         scf_mode,zval,
     &         mag_fit,mag_weight,mag_energy,occ_sp,mag_de,
     &         vtype,0,pcoeff,pen1,fcostB,wf_match,vtmp)
            endif 
            
            ! five-point stencil (much more accurate than central finite diff.)
            if (five_point) then 
             grad(i)=(-fcostA+8.d0*fcost1-8.d0*fcost2+fcostB)/12.d0/step
            else 
             grad(i)=(fcost1-focst2)/2.d0/step
            endif 
          enddo

        elseif (task(1:5) .eq. 'NEW_X') then 
          if (mag_fit<=0) then 
           write(6,'(a,i3,a,es10.3,a,es9.2,a,es8.2)')'iter: ',iter,
     &     ' NEW_X, fcost:',fcost,'  eig_match:',fcost-pen,
     &     ' |g|: ',sqrt(sum(grad**2))
          else 
           write(6,'(a,i3,a,es10.3,a,es9.2,a,es8.2,
     &           a,f11.6,a,f11.6,a)')'iter: ',iter,
     &     ' NEW_X, fcost:',fcost,'  eig_match:',fcost-pen,
     &     ' |g|: ',sqrt(sum(grad**2)),' dE_mag:',mag_de,
     &     ' (',mag_energy,')'
          endif 
          vloc_new = vloc   ! back up vloc
          param_new = param ! back up param 
          iter = iter + 1
          if (fcost<1e-14) then 
            write(6,'(/,a,/)')'** L-BFGS reaches minimum, finished **'
            exit 
          endif 
        else
          print *,'exit BFGS. task: ',trim(task),' last call to cost():'
          ! last time call cost() 
          call cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,
     &     z,r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &     rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &     icut,eref,ups_ref,npdim,npar,param_new,iorb_start,
     &     scf_mode,zval,mag_fit,mag_weight,mag_energy,occ_sp,mag_de,
     &     vtype,0,pcoeff,pen,fcost,wf_match,vloc)
          write(6,'(/,a,es12.4,a)')' *** fcost:    ',fcost,' ***'
          write(6,'(/,a,es12.4,a)')' *** penalty:  ',pen,' ***'
          write(6,'(a,es12.4,a,/)')' *** wf_match: ',wf_match,' ***'
          exit
        endif 
      enddo
      vloc = vloc_new 
      param = param_new 

      print *,'------ final param ------'
      if (vtype==6) then 
        do i=1,npar+4
          write(6,'(a,i2,a,es14.6)')'c(',i,'):',param(i)
        enddo
        print *,''
      elseif (vtype==7) then 
        do i=1,npar+5
          write(6,'(a,i2,a,es14.6)')'c(',i,'):',param(i)
        enddo
        print *,''
      else
        do i=1,npar+3
          write(6,'(a,i2,a,es14.6)')'c(',i,'):',param(i)
        enddo
        print *,''
      endif 

      deallocate(lbfgs_wa,lbfgs_nbd,lbfgs_l,lbfgs_u,grad,lbfgs_iwa)
      print *,'>>> optimization is done <<< '

      !==================================
      ! final calculation using vlocal
      !==================================
      do i=1,mx
        do j=1,ms 
          vlocs(i,j)=vloc(i)
        enddo
      enddo
      write(6,'(/,a,/)')'*** final psatom_chen calc. with vlpsp *** '
      vtmp=0.d0
      zncpp=0.d0
      it=1 
      if (scf_mode) it=itmx
      call psatom(it,igr,iexcnow,zncpp,0,norb,n,l,f,m,
     &  e,nin,de,mmax,r,vtmp,vlocs,rhos,rhosp,rhospp,
     &  dc,dcp,dcpp,ups,tkli,svm,svm_roks)

      open(file='vks.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es14.6)')r(i),vlocs(i,1)+vtmp(i,1)
      enddo
      close(111)
      print *,' >>> vks.dat is written <<<'

      ! compute norm for r<rcut
      norm=0.d0
      do i=1,norb
        do j=1,icut-3
          norm(i)=norm(i)+al*r(j)*ups(j,i)**2
        enddo
        norm(i)=norm(i) + al*(23.0d0*r(icut-2)*ups(icut-2,i)**2
     &              + 28.0d0*r(icut-1)*ups(icut-1,i)**2
     &              +  9.0d0*r(icut)*ups(icut,i)**2)/24.0d0
      enddo

      open(file='rho_lps.dat',action='write',unit=111)
      do i=1,mmax
        write(111,*)r(i),rhos(i,1)
      enddo
      close(111)
      
      ! check eigenvalues 
      do i=1,norb
        write(6,'(a,2i2,a,f9.3,a,f10.5,a,f10.5,a,es10.2,a)') 
     &  'n,l: ',n(i),l(i),' => occ:',f(i),' eig:',e(i),' eref:',eref(i),
     &  ' dE:',(e(i)-eref(i))*27.2114,' eV'
      enddo
      print *,''
      ! check wf norms 
      do i=1,norb
        write(6,'(a,i2,a,es10.3,a,es10.3,a,f8.4,a,f7.4,a)')'orb:',i,
     &    '  u(icut): ',ups(icut,i),' u(ref): ',ups_ref(icut,i),
     &    '  norm(r<rc):',norm(i),' (ref:',norm0(i),')'
      enddo     
      print *,''
    
      ! unscreen vloc
      if (.not. scf_mode) then 
        print *,' >>> unscreen vlocal <<< '
        do i=1,mmax
          rho(i)=rhos(i,1)
          rhop(i)=rhosp(i,1)
          rhopp(i)=rhospp(i,1)
        enddo
        do ispin=2,ispinmx
          do i=1,mmax
            rho(i)=rho(i)+rhos(i,ispin)
            rhop(i)=rhop(i)+rhosp(i,ispin)
            rhopp(i)=rhopp(i)+rhospp(i,ispin)
          enddo
        enddo
        call vestat(mmax,cval,eeel,r,rho,ves,.true.)
        call vexcor(iexc,mmax,r,rho,rhop,rhopp,vxc,eexc,ex,ec,.true.)
        do i=1,mmax
          vloc(i)=vloc(i)-vxc(i,1)-ves(i)
        enddo
      endif 

      ! find Coulomb tail
      do ir=1,mmax
        vvo(ir)=-zval/r(ir)
      enddo
      epstail=1.d-3
      do i=mmax,1,-1
        if(abs(vloc(i)-vvo(i)) .gt. epstail) then
          icoul=i
          write(6,'(a,f7.2,a)')'*** rcut: ',r(icut),'  ***'
          write(6,'(a,f7.2,a,f7.2,a)') 
     &     '*** Coulomb tail starts at ',r(i),'   rmax:',r(mmax),' ***'
          exit 
        endif 
      enddo

      open(file='vlps.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es14.6)')r(i),vloc(i)
      enddo
      close(111)
      print *,' >>> vlps.dat is written <<<'
      open(file='vhxc.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es14.6)')r(i),vxc(i,1)+ves(i)
      enddo
      close(111)
      open(file='rho_val.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es14.6)')r(i),rhos(i,1)
      enddo
      close(111)

      open(file='orb1.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(3es14.6)')r(i),ups(i,1),ups_ref(i,1)
      enddo
      close(111)
      open(file='orb2.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(3es14.6)')r(i),ups(i,2),ups_ref(i,2)
      enddo
      close(111)
      open(file='orb3.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(3es14.6)')r(i),ups(i,3),ups_ref(i,3)
      enddo
      close(111)
      open(file='orb4.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(3es14.6)')r(i),ups(i,4),ups_ref(i,4)
      enddo
      close(111)
      open(file='orb5.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(3es14.6)')r(i),ups(i,5),ups_ref(i,5)
      enddo
      close(111)

      ! set vlocal to vorb 
      write(6,'(a)')'  >>> set vorb=vloc, vbare=vloc, uu=ups <<<'
      vorb_nl = vorb 
      do j=1,ms
        vorb(:,j)=vloc
        vbare(:,j)=vloc
        uu(:,j)=ups(:,j)
      enddo

      call outlps(mmax,z,zval,norb,npdim,npar,param,wfopt,
     &            pcoeff,vtype,iorb_start,rcut,icut,icoul,r,vloc)
      endif ! mk_vlocal
      
      !============================================
      !======== test other configurations =========
      !============================================
      test_config=.false.
      if (test_config) then 
        print *,''
        print *,'  >> other configurations <<<'
        print *,''
        occ2(1)=1.9d0
        occ2(2)=0.1d0
        occ2(3)=0.d0 
        norb2=3
        vtmp=0.d0
        zncpp=0.d0 
        it=itmx 
        call psatom(it,igr,iexcnow,zncpp,0,norb2,n,l,occ2,m,
     &  e,nin,de,mmax,r,vtmp,vlocs,rhos,rhosp,rhospp,
     &  dc,dcp,dcpp,ups,tkli,svm,svm_roks)
        print *,'e:',e(1:norb2)

        it=itmx 
        vtmp=0.d0
        call psatom(it,igr,iexcnow,zncpp,0,norb2,n,l,occ2,m,
     &  e,nin,de,mmax,r,vtmp,vorb_nl,rhos,rhosp,rhospp,
     &  dc,dcp,dcpp,ups,tkli,svm,svm_roks)
        print *,'e:',e(1:norb2)
      endif 

      print *,''
      print *,'  >>> continue pswatch <<<'
      print *,''
c=====================================
c chen: end of solving vlocal 
c=====================================


      it=itmx
      call psatom(it,igr,iexcnow,0.d0,0,norb,n,l,f,m,e,nin,de,mmax,r
     &, vee,vorb,rhos,rhosp,rhospp,dc,dcp,dcpp,ups,tkli,svm,svm_roks)
      iexc=abs(iexc)

c compute total density to retain compatibility with existing code
      do i=1,mmax
        vi(i)=vee(i,1)
        rho(i)=rhos(i,1)
        rhop(i)=rhosp(i,1)
        rhopp(i)=rhospp(i,1)
      enddo
      do ispin=2,ispinmx
        do i=1,mmax
          rho(i)=rho(i)+rhos(i,ispin)
          rhop(i)=rhop(i)+rhosp(i,ispin)
          rhopp(i)=rhopp(i)+rhospp(i,ispin)
        enddo
      enddo

c meta gga
      ek_from_stat = fmom(0,mmax,al,1.d0,r,stat)
      if(tmgga .and. tnlc) then
        call dadv(mx,1,mmax,stat,stat_core,stat)
        ek_core_from_stat = fmom(0,mmax,al,1.d0,r,stat_core)
        ek_from_stat = ek_core_from_stat + ek_from_stat
        do i=1,mmax
          rho(i)=rho(i)+dc(i)
          rhop(i)=rhop(i)+dcp(i)
          rhopp(i)=rhopp(i)+dcpp(i)
        enddo
        open(11,file='alpha-ps.ekin_core')
        open(12,file='alpha-ps.s_eff')   ! effective scaled gradient
        open(13,file='alpha-ps.t_ratio') ! effective t_w/t_s
        do i=1,mmax
          write(11,'(e13.6,1x,e15.8)') r(i),stat(i)+1e-12
          vara = fx_xvar(rho(i),rhop(i),stat(i))
          write(12,'(e13.6,1x,e15.8)') r(i),vara+1e-12
          write(13,'(e13.6,1x,e15.8)')
     1      r(i),rhop(i)**2/( 8d0* rho(i)*stat(i) )+1e-12
        enddo
        write(11,*)'&'
        write(12,*)'&'
        write(13,*)'&'
        do i=1,mmax
          write(11,'(e13.6,1x,e15.8)') r(i),stat(i)-stat_core(i)+1e-12
          vara = rho(i) - dc(i)
          varb = rhop(i) - dcp(i)
          varc = stat(i) - stat_core(i)
          vard = fx_xvar(vara,varb,varc)
          write(12,'(e13.6,1x,e15.8)') r(i),vard+1e-12
          write(13,'(e13.6,1x,e15.8)') r(i),(rhop(i)-dcp(i))**2/
     1      ( 8d0* (rho(i)-dc(i)) * (stat(i)-stat_core(i)) )+1e-12
        enddo
        write(11,*)'&'
        write(12,*)'&'
        write(13,*)'&'
        do i=1,mmax
          write(11,'(e13.6,1x,e15.8)') r(i),stat_core(i)+1e-12
          vara = fx_xvar(dc(i),dcp(i),stat_core(i))
          write(12,'(e13.6,1x,e15.8)') r(i),vara+1e-12
          write(13,'(e13.6,1x,e15.8)') r(i),
     1       dcp(i)**2/( 8d0* dc(i)*stat_core(i) )+1e-12
        enddo
        close(11)
        close(12)
        close(13)
        do i=1,mmax
          rho(i)=rho(i)-dc(i)
          rhop(i)=rhop(i)-dcp(i)
          rhopp(i)=rhopp(i)-dcpp(i)
        enddo
      endif

c total charges
      tc_val=fmom(0,mmax,al,1.d0,r,rho)
      if(tnlc) then
        tc_core=fmom(0,mmax,al,1.d0,r,dc)
        if(igr .eq. 1) then
          t1c=-fmom(1,mmax,al,tc_core,r,dcp)/3.d0
          t2c=fmom(2,mmax,al,tc_core,r,dcpp)/12.d0
          if(tmgga) then
            do ir=1,mmax
              vvo(ir) = 0.125*dcp(ir)**2/max(1.d-10,dc(ir))
            enddo
            t3c=fmom(0,mmax,al,1.d0,r,vvo)
          endif
        endif
      endif

c total energy components
c hartree energy
      call vestat(mmax,cval,eeel,r,rho,ves,.true.)

c xc energy
      if(tnlc) then
        if(ispinmx.gt.1) then
          do i=1,mmax
            dc(i)=0.5d0*dc(i)
            dcp(i)=0.5d0*dcp(i)
            dcpp(i)=0.5d0*dcpp(i)
          enddo
        endif
        do ispin=1,ispinmx
          do i=1,mmax
            rhos(i,ispin)=rhos(i,ispin)+dc(i)
            rhosp(i,ispin)=rhosp(i,ispin)+dcp(i)
            rhospp(i,ispin)=rhospp(i,ispin)+dcpp(i)
          enddo
        enddo
      endif
      if(ispinmx.eq.1) then
        call vexcor(iexc,mmax,r,rhos(1,1),rhosp(1,1),rhospp(1,1),
     1              vxc(1,1),eexc,ex,ec,.true.)
      else
        call vexcos(iexcs,mmax,mmax,r,rhos,rhosp,rhospp,vxc,
     1              eexc,ex,ec,.true.)
      endif
      if(tnlc) then
        do ispin=1,ispinmx
          do i=1,mmax
            rhos(i,ispin)=rhos(i,ispin)-dc(i)
            rhosp(i,ispin)=rhosp(i,ispin)-dcp(i)
            rhospp(i,ispin)=rhospp(i,ispin)-dcpp(i)
          enddo
        enddo
        if(ispinmx.gt.1) then
          do i=1,mmax
            dc(i)=2.d0*dc(i)
            dcp(i)=2.d0*dcp(i)
            dcpp(i)=2.d0*dcpp(i)
          enddo
        endif
      endif
c xc potential energy
      evxc=dmelm(mmax,al,r,vxc(1,1),rhos(1,1))
     1    +dmelm(mmax,al,r,vxc(1,2),rhos(1,2))
c valence xc energy
      vexc=dmelm(mmax,al,r,dexc,rho)

      if(tkli) then
        ex=0.d0
        ex_roks=0.d0
        do i=1,norb
         ex=ex+f(i)*svm(i)
         ex_roks=ex_roks+f(i)*svm_roks(i)
        enddo
        ex=0.5d0*ex
        ex_roks=0.5d0*ex_roks
c       vi holds the hartree + xc potential, vxc the correlation potential
        evxc=dmelm(mmax,al,r,vi,rho)+dmelm(mmax,al,r,vxc,rho)-eeel
        vexc=ex+dmelm(mmax,al,r,dexc,rho)
      endif
      cexc=ex+ec-vexc

c ionic pseudopotential energy 
      evps=0.d0
      do i=1,norb
        if(f(i) .gt. 0.0) then
          ep(i)=gltfmv(mmax,al,r,ups(1,i),vorb(1,i),ups(1,i))
          evps=evps+f(i)*ep(i)
        endif
      enddo
c local pseudopotential energy
      evps_loc=dmelm(mmax,al,r,vbare(1,ll_loc),rho)
c screened pseudopotential energy
      epsp=evps
      do ispin=1,ispinmx
        epsp=epsp+dmelm(mmax,al,r,vee(1,ispin),rhos(1,ispin))
      enddo
c kinetic energy from eigenvalues
      ekin=dspv(ms,1,norb,f,e)-epsp
      etot=ekin+.5*eeel+evps+ex+ec

      call labelmap(int(zfull),symz,iexc,symxc)
      write(iu,'(a30,2x,a8)')   'chemical symbol', symz
      write(iu,'(a30,2x,f5.2)') 'nuclear charge', zfull
      write(iu,'(a30,2x,f5.2)') 'number of valence electrons',cval
      write(iu,'(a30,2x,i2)')   'number of valence states', nv
      write(iu,'(a30,2x,i2,2x,a40)')
     &    'exchange-correlation model', iexc,symxc
      if(tspin) write(iu,'(a30)')'spin-polarized calculation'
      if(tnlc) write(iu,'(a30)') 'nonlinear core-valence XC'
      write(iu,'(a30,2x,i4,f12.6,2x,e13.6)')
     &   'parameters radial mesh',mmax,amesh,r(1)
      call labels(1,llmx-1,1,sn,sl,sm)
      write(iu,'(a30,3x,a1)') 'input pseudopotentials up to',sl
      if(tcut_global)
     &    write(iu,'(a30,1x,f8.3)') 'global cutoff radius',rcut_global
      
      write(iu,620)
  619 format(/'fhi pseudopotential tool pslp - version ',a10,/)
  620 format(/10x,'=== pseudo atom (Hartree a.u.) ===',//
     &  '<        n     l   occupation  eigenvalue(eV)  ',
     &  'potential energy')
      do i=1,norb
        if(ispinmx.gt.1) then
          if(m(i).eq.1) then
            write(iu,631) i,n(i),l(i),f(i),e(i)*ry2,ep(i)
          else
            write(iu,632) i,n(i),l(i),f(i),e(i)*ry2,ep(i)
          endif
        else
            write(iu,630) i,n(i),l(i),f(i),e(i)*ry2,ep(i)
        endif
      enddo
  630 format('< ',i2,4x,2(i2,4x),f8.4,2x,f12.4,4x,f12.5)
  631 format('< ',i2,4x,i2,'_u',2x,i2,4x,f8.4,2x,f12.4,4x,f12.5)
  632 format('< ',i2,4x,i2,'_d',2x,i2,4x,f8.4,2x,f12.4,4x,f12.5)

  608 format(a30,2x,f12.5)
  610 format(a30,2x,f12.5,2x,f12.5)
      write(iu,*) 
      if(tkli) then
        write(iu,610) 'total energy [rks,roks]',etot,etot-ex+ex_roks
      else
        write(iu,608) 'total energy',etot
      endif
      write(iu,608) 'kinetic energy',ekin
      write(iu,608) 'ionic pseudopotential energy',evps
      write(iu,608) 'hartree energy',0.5*eeel
      if(tkli) then
        write(iu,610) 'xc energy [rks,roks]',ex+ec,ex_roks+ec
        write(iu,608) ' c energy           ',ec
      else
        write(iu,608) 'xc energy',ex+ec
      endif
      write(iu,608) 'local potential energy',evps_loc
      write(iu,608) 'xc potential energy',evxc
      if(tmgga) then
        write(iu,608) 'g kinetic energy',ek_from_stat
      endif
      if(tnlc) then
        write(iu,608) 'xc energy core',cexc
        write(iu,608) 'xc energy valence',vexc
        if(tmgga) then
          write(iu,608) 'g kinetic energy core',ek_core_from_stat
          write(iu,608) 'g kinetic energy valence',
     &                   ek_from_stat-ek_core_from_stat
        endif
      endif
      write(iu,608) 'integrated valence density', tc_val
      if(tnlc) then
        write(iu,608) 'integrated core density', tc_core
        if(igr.eq.1) then
          write(iu,608) ' ... 1st derivative test',t1c
          write(iu,608) ' ... 2nd derivative test',t2c
          if(tmgga) then
            write(iu,608) ' ... kinetic energy test',t3c
          endif
        endif
        do i=nin-10,1,-1
          if(rho(i) .le. dc(i)) goto 590
        enddo
 590    write(iu,608) 'estimated equidensity radius >',r(i+1)
      endif
      write(iu,'(a30,2x,i12,2x,a12,1pe9.1)')
     &  'number of iterations', it, 'convergence', de

c for plot x y range
      vl=1.d20
      vh=-1.d20
      do ll=1,llmx
        call dextv(mx,1,mmax,vbare(1,ll),bl,bh)
        vl=min(bl,vl)
        vh=max(bh,vh)
      enddo
      i=int(vl-int(vl/5)*5)-1+int(vl/5)*5
      j=int(vh-int(vh/5)*5)+1+int(vh/5)*5
      ll=max(1,(j-i)/4)
      write(iu,637) 'y range plot',i,j,ll
  637 format(a30,6x,3i4)

c spin polarized calculation makes sense up to here only, so stop
      if(tspin) then
        write(iu,'(/a)') 'pslp - spin-polarized pseudoatom done'
        stop 
      endif
     
c manipulate and output pseudopotential
      do ll=1,llmx
        ntl(ll)=0
        do i=1,norb
          if(n(i) .eq. ll) ntl(ll)=i
        enddo
c pseudopotential output with calculated wavefunctions
        if(ntl(ll) .gt. 0) then
          call dcpv(mx,1,mmax,ups(1,ntl(ll)),uout(1,ll))
        else
          call dnuv(mx,1,mmax,0.d0,uout(1,ll))
        endif
c screened pseudopotentials
        call dadv(mx,1,mmax,vbare(1,ll),vi,vorb(1,ll))

        do ir=1,mmax
          rxx(ir)=0.d0
        enddo
      enddo

c find radius outside of which pseudopotential's tail is coulomb like
      epstail=1.d-3
      do ir=1,mmax
        !vvo(ir)=-cval/r(ir)
        vvo(ir)=-zval/r(ir)
      enddo
      do ll=1,llmx
        rc(ll)=0.d0
        do ir=mmax,1,-1
          if(abs(vbare(ir,ll)-vvo(ir)).gt.epstail) then
            rc(ll)=r(ir)
            exit
          endif
        enddo
      enddo

      if(tcut_global) then
        write(ie,'(/,a,/,a)') 
     1    ' %pslp: WARNING - the output pseudopotential corresponds' ,
     1    '        to the globally cutoff input potential!'
        call outpot(40,mx,mmax,llmx,cval,rc,r,uu,vbare,vorb)
      else
        call outpot(40,mx,mmax,llmx,cval,rc,r,uout,vbare,vorb)
      endif
      if(tnlc) call outcore(27,mmax,1,r,dc,dcp,dcpp)

c radial pseudo-wavefunctions [the u(r)]
      do i=1,nv
        ninu(i)=nin
      enddo
      call out38(mx,1,nv,ninu,n,l,r,ups)

c informal stuff that is under construction
      if(tmgga) then
        if(tnlc) then
          do ir=1,mmax
            rho(ir) = rho(ir) + dc(ir)
            rhop(ir) = rhop(ir) + dcp(ir)
            rhopp(ir) = rhopp(ir) + dcpp(ir)
          enddo
        endif

        write(iu,'(/,3x,a)') 
     &                   '----- under construction ... begin -------'

        call vexcor(6,mmax,r,rho,rhop,rhopp,vxc,eexc,exv,ecv,.true.)
        write(iu,608) 'pbe post total energy',etot-ex-ec+exv+ecv
        call vexcor(16,mmax,r,rho,rhop,rhopp,vxc,eexc,exv,ecv,.true.)
        write(iu,608) 'mgga post total energy',etot-ex-ec+exv+ecv
     
        call vexcor(iexc,mmax,r,rho,rhop,rhopp,vxc,eexc,exv,ecv,.true.)
        write(iu,608) 'valence xc energy',exv+ecv
        write(iu,608) 'total linearized xc energy',epcv+exv+ecv
        write(iu,608) 'post linearized total energy',
     &    etot-ex-ec+epcv+exv+ecv
        if(tnlc) then
          do ir=1,mmax
            rho(ir) = rho(ir) - dc(ir)
            rhop(ir) = rhop(ir) - dcp(ir)
            rhopp(ir) = rhopp(ir) - dcpp(ir)
          enddo
        endif
        write(iu,'(3x,a)') '----- ... end ----------------------------'
      endif

      write(iu,*) 
      write(iu,*) 'pslp - pseudoatom done - now testing'

      if(.not. tlgd .and. .not. tkb) stop 'pslp - pseudoatom done'

c logarithmic derivatives and ghost state test
c read all-electron potential
      if(tlgd) then
        call dnuv(mx,1,mmax,0.d0,vae)
        if(irl .eq. 1) then
          write(iu,638) 
        else if( irl .eq. 2) then
          write(iu,639) 
        endif
        terr=.true.
        read(37,*,end=50,err=50) xxx,mmaxtmp,iexcae,irlae
        if(iexc .ne. iexcae) 
     &    write(ie,*)'& pslp - XC-type mismatch - using IEXC', iexc
        if(irlae .ne. irl) 
     &    write(ie,*)'& pslp - Relativity type mismatch - using IRL',irl
        do i=1,mmaxtmp
          read(37,*,end=48,err=48) rxx(i),vae(i)
        enddo
        terr=.false.
   48   if(i .le. mmaxtmp) terr=.true.
        if( abs(amesh-rxx(2)/rxx(1)) .gt. 5e-6 )
     &    write(ie,*)
     &      '& pslp - warning: grid from unit fort.37 incompatible'
   50   if(terr) then
          write(ie,*) 
     &      '& pslp - warning: bad/missing full potential file fort.37'
        endif
      endif
  638 format(/' --- assuming scalar-relativistic all-electron atom ---')
  639 format(/' --- assuming nonrelativistic all-electron atom ---')

c input or calculated wavefunctions for kb projectors
      call labels(1,l_loc,1,sn,sl,sm)
      write(iu,645) sl
      if(tkb) then
        if(tiwf) then 
          write(iu,641)
        else 
          write(iu,643)
        endif 
      endif
  641 format(' --- input wavefunctions used for kb potentials ---')
  643 format(' --- calculated wavefunctions used for kb potentials ---')
  645 format(/' --- ',a1,' component taken as local potential ---') 

      do ll=1,llmx
        if(.not. tiwf) call dcpv(mx,1,mmax,ups(1,ntl(ll)),uu(1,ll))
        ep(ll)=e(ntl(ll))
      enddo

c diagnostic radius
      if(rld .le. 0.0) rld=rcovalent(int(zfull))*1.3

      call ppcheck(iu,tnrl,tlgd,tkb,l_loc,lmax,lbeg,lend,mmax,rld,ep,r
     &,  vae,vi,uu(1,1),vbare(1,1))

c reciprocal space analysis (occupied states)
      write(iu,652)
 652  format(/' --- kinetic energy convergence in momentum space ---',//
     & 5x,'l  n  bracket   cutoff    norm   kinet. energy   cutoff',/
     & 12x,'(eV)      (Ry)               (Hartree)    (eV)')

      do ll=llbeg,llend
        neff=ll
        do i=1,nv
          if(l(i) .eq. ll-1 .and. f(i) .gt. 0.0) then
            call dadv(mx,1,mmax,vi,vbare(1,ll),vsl)
            call kinkon(iu,mmax,400.d0,0.05d0,neff,ll-1,e(i),r,vsl)
            neff=neff+1
          endif
        enddo
      enddo
c
      write(iu,'(1x,a)') '--- coulomb tail of pseudopotentials ---'
      write(iu,'(5x,a,1p,e8.1,a)') 'Tolerance',epstail,' is met for'
      do ll=1,llmx
       write(iu,'(5x,a,i2,a,f8.3)') 'l=',ll,' at radii >',rc(ll)
      enddo

      write(iu,'(/,a)') ' --- done & exiting ---'
c
c close file
      if(tkli) close(80)
c            
      end
c




