c $Header:$
cmf 27-06-03 added spin polarization for selected xc functionals,
cmf          some cleansing
cmf 04-09-02 added kli exchange functional
c
      program fhipp
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
c Self-consistent electronic structure for atoms 
c based on density-functional theory assuming spherical screening
c - nonrelativistic/ scalar-relativistic 
c - local density/ generalized gradient approximation,
c   some useable for spin-polarized calculation (all-electron atom only)
c - kinetic energy-density dependent functionals (experimental)
c - exact kohn-sham exchange within kli approximation (experimental)
c - all-electron/ valence-only (frozen core) calculation 
c - generalized norm-conserving pseudopotentials of Hamann or 
c   Troullier-Martins type
c - linearized/ nonlinear core-valence exchange-correlation
c 
c Familiarity with the references below should be regarded as a 
c prerequisite to a profitable use of this program. 
c Please cite these references when publishing results obtained with
c pseudopotentials constructed with the help of this program.
c
c [1] M. Fuchs, M. Scheffler, Comput. Phys. Commun. 119, 67-98 (1999)
c
c The pseudopotential generation in this program is based on 
c
c [2] D.R. Hamann, Phys. Rev. B 40, 2980 (1989)
c     "Generalized norm-conserving pseudopotentials"
c
c [3] N. Troullier, J.L. Martins, Phys. Rev. B 43, 1993 (1991)
c     "Efficient pseudopotentials for plane-wave calculations"
c
c Before employing pseudopotentials check their transferability. 
c In particular when using them in fully separable form make sure 
c to avoid ghost-states, see [1] and e.g. 
c
c [4] X. Gonze, R. Stumpf, M. Scheffler, Phys. Rev. B 44, 8503 (1991)
c     "Analysis of separable potentials"
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
c fort.20     control parameters
c fort.22     atomic & pseudopotential parameters
c fort.18     frozen core density (optional)
c fort.36     effective potential (optional)
c
c output files
c fort.10-14  scratch for meta ggas
c fort.19     full core density
c fort.23     monitoring data (-> parameter.h)
c fort.24     monitoring pseudocore 
c fort.25     pseudovalence density
c fort.27     modelcore density
c fort.28     monitoring unscreening
c fort.37     effective potential
c fort.38     all-electron wavefunctions
c fort.39     pseudo and all-electron valence wavefunctions
c fort.4[0-4] ionic pseudopotentials
c fort.4[5-9] screened pseudopotentials
c fort.80     info file for kli case
c-----------------------------------------------------------------------
c
c input description fort.20 and fort.22
c
c fort.20 - line 1
c tdopsp generate pseudo potentials if true (default true)
c tnrl   nonrelativistic atom if true (default false)
c
c fort.20 - line 2
c tspin  .true.  perform a spin-polarized calculation.
c         Default is .false. . Implemented for selected xc functionals
c         only. Requires modified input file fort.22 .
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
c        = 0           linearized core-valence XC, no pseudo core
c        > 0           nonlinear core-valence XC, pseudo core inside
c                      cutoff radius, full core outside
c
c fort.22 - line 2 
c n(i)   principal quantum number
c l(i)   angular momentum
c f(i)   occupation number (two entries for spin-polarized mode)
c
c fort.22 - line 2+nc+nv
c ltmx        maximum angular momentum to generate pseudopotentials for
c spptype     choice of pseudopotential scheme
c             h  Hamann type
c             t  Troullier-Martins type
c
c fort.22 - line 3+nc+nv (optional)
c lt          apply optional paramteres to pseudopotential with this
c             angular momentum
c rct         core cutoff radius
c et          reference energy
c spptype     pseudopotential type 
c
c sample input fort.21 for silicon (3s3p3d potentials)
c > 2 0 2 2 0.0     : l_loc lbeg lend lmax rld
c >.t. .t. .t. .f.  : tlgd  tkb  tiwf tnrl
c
c sample input fort.22 for silicon (Hamann potential, LDA, linearized CV XC)
c > 14.0 3 2 8 0.0  : z nc nv iexc rnlc
c > 1   0  2.0      : n l f
c > 2   0  2.0
c > 2   1  6.0
c > 3   0  2.0
c > 3   1  2.0
c > 2 h             : ltmx spptype
c > 1 1.5 0.0 h     : lt rct et spptype (alternative parameters)
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
      implicit none

      include 'parameter.h'
      include 'default.h'

      logical tcore,tdopsp,tnrl,tsic,tinpot,tspin
      logical tinopt(ms),t_pp_mode(ms)
      character*8 symz,symxc*40,srev*10,spptype*25
      character*20 spp
      character*1 sppopt(ms)
      integer i,ir,iu,iexc,it,l1,lmax,lt,irl,igr,nopt,ishift
      integer nin,mmax,nc,nv,itype,ltmx,i_pp_type_def,j,nstart
      integer iexcnow,iexcs,ispin,nspinmx,norb,norbnow
      integer n(ms),l(ms),ninu(ms),igrmap(20),i_pp_type(ms),np(ms)
      integer inode(ms),ipeak(ms),ltopt(ms),m(ms)
      real*8  al,amesh,csc,csc1,fmom,dmelm,ex,ec,emax,et,rcmax,rct,sf
      real*8  z,eexc,eeel,de,tc,t1,t2,evxc,ecou,dexc,rnlc,zv,zc,ekin
      real*8  gltfmv,etot,eeig,epot_core,ekin_core,dspv,defrtm
      real*8  esic,esic_core,etot_sic,epot_sic,fup,fdn,dcfac,epot
      real*8  f(ms),e(ms),rpk(ms),ves(mx),uu(mx,ms),ups(mx,ms)
      real*8  fp(ms),ep(ms),rc(ms),rc_d(ms,2),rc_dmax(2),vsc(2),vsc1(2)
      real*8  e_sic(ms),vorb(mx,ms),vsic(mx,ms),rxx(mx)
      real*8  r(mx),vi(mx),vae(mx),rctopt(ms),ectopt(ms)
      real*8  dc(mx),dcp(mx),dcpp(mx),rho(mx),rhop(mx),rhopp(mx)
      real*8  ds(mx,ms),dsp(mx,ms),dspp(mx,ms),vxc(mx,ms),vee(mx,ms)
      real*8  vtmp(mx,ms)
      real*8  dcs(mx,ms),dcsp(mx,ms),dcspp(mx,ms)
c kli
      logical tkli
      real*8  svm(ms),svkli(ms)
c meta gga
      logical tmgga
      real*8  ek_from_stat,ek_core_from_stat,stat_orb
      real*8  fx_xvar,vara,varb,varc,vard
      real*8  uup(mx,ms)
      real*8  stat,stat_core
      common/use_mgga/stat(mx),stat_core(mx)
c
c &&&&&&&&&&&&& REVISION &&&&&&&&&&&&&
      parameter(srev='rev270603B')
      parameter(tsic=.false.)
      data igrmap/0,0,0,1,1,1,0,0,1,1,0,0,0,1,1,1,1,1,1,0/
      external defrtm,spp,fmom,dmelm,dspv,stat_orb

ccccccccccccccccccccccccccccccccccc chen ccccccccccccc
      real*8 rhos(mx,ms),rhosp(mx,ms),rhospp(mx,ms)
      real*8 svm_roks(ms)
      real*8 zcore, zval, cval, dtmp
      real*8 epstail,vbare_v(mx),ves_v(mx),vxc_v(mx)    ! KS, hartree, xc potential for valence 
      real*8 rho_v(mx),rho_vp(mx),rho_vpp(mx) ! valence electron and derivatives
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c begin
c defaults
      iexc=8
      iu=iofhipp
      tdopsp=.true.
      tnrl  =.false.
      tcore =.true.
      tkli  =.false.
      tmgga =.false.
      tspin =.false.
      lmax=0
      irl=1
      igr=0
      nspinmx=1

c input 
      read(20,*,end=10,err=10) tdopsp,tnrl
   10 read(20,*,end=11,err=11) tspin
      if(tspin) then 
        nspinmx=2
        tdopsp=.false.
      endif

   11 continue
      read(22,*,end=65,err=65) z,nc,nv,iexc,rnlc
      z=abs(z)
      nc=abs(nc)
      nv=abs(nv)
      norb=nc+nv
      do i=1,norb
        if(nspinmx.eq.1) then
          read(22,*,err=65,end=65) n(i),l(i),f(i)
        else
          read(22,*,err=65,end=65) n(i),l(i),f(i),f(i+norb)
        endif
        n(i)=abs(n(i))
        l(i)=abs(l(i))
        m(i)=1
        f(i)=abs(f(i))
        if(n(i) .le. l(i)) stop 'fhipp - data: l .le. n'
        if(l(i) .gt. 4)    stop 'fhipp - data: l out of range (0 ... 4)'
      enddo
      if(nspinmx.gt.1) then
        do i=1,norb
c for occupation of spin states ...
          if(f(i+norb).eq.0.d0) then
c ... follow 2nd Hund rule
            fup=min(f(i),2.d0*l(i)+1.d0)
            fdn=f(i)-fup
          else
            if(f(i).gt.f(i+norb)) then
c ... follow fixed majority spin
              fup=f(i+norb)
              fdn=f(i)-fup
            else
c ... use occupancies of input
              fup=f(i)
              fdn=f(i+norb)
            endif
          endif
          f(i)=fup
          n(i+norb)=n(i)
          l(i+norb)=l(i)
          m(i+norb)=2
          f(i+norb)=fdn
        enddo
      endif
      norb=norb*nspinmx
      if(nc.lt.1) tcore=.false.
      sf=0.d0
      do i=1,norb
        sf=sf+f(i)
      enddo
c debug 
!     write(ie,*)'%fhipp - state count: norb=',norb,
!    1  ' of',(nc+nv)*nspinmx
!     do i=1,norb
!       write(ie,*) n(i),l(i),f(i) 
!     enddo

      if(sf .gt. z) then
        write(iu,'(a)') '%fhipp - WARNING: negative ion.'
      endif
      iexc=abs(iexc)
      if(iexc.eq.12 .or. iexc.eq.13 .or. iexc.eq.18 .or. iexc.eq.19)
     1   tkli=.true.
      if(iexc .eq. 16) tmgga = .true. ! meta gga
c map for spin polarized calculation
      if(nspinmx.gt.1) then
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
     1      '%fhipp - ERROR: spin polarized calculation not implemented'
     1     ,' for input iexc=',iexc,'. ACTION: set iexc=3,4,5,6,8.'
          stop
        endif
      endif
      rnlc=abs(rnlc)
      if(tsic .and. (iexc.ne.8 .or. tspin)) then
        write(ie,'(a,i3,a/,a/,a)')
     1    '%fhipp - ERROR: SIC not implemented for iexc=',-iexc,' or',
     1    ' spin-polarized calculations. ACTION: set iexc=-8 and',
     1    ' spin-unpolarized mode.'
        stop
      endif
c
c open files
c
      if(tkli) open(unit=80,file='kliinfoAE.dat',status='unknown')
c
c all-electron atom:    irl = 1 scalar-relativistic
c                       irl = 2 non-relativistic
c gradient corrections: igr = 0 gradients not evaluated
c                       igr = 1 gradients evaluated
c
      if(tnrl) irl=2
      igr=igrmap(max(iexc,1))
      if(tcore) igr=1

c grid generation and initial potential 
      call logmesh(mmax,z,amesh,al,r)

c initial potential and eigenvalues
      norbnow=norb/nspinmx
      do i=1,norbnow
        fp(i)=f(i)
        if(tspin) fp(i)=fp(i)+f(i+norbnow)
      enddo
      call atomini(norbnow,mmax,z,n,fp,r,vi,e)

c alternatively do non self-consistent calculation at input potential
      nin=0
      read(36,*,end=48,err=48) rxx(1),nin,i,i
      do i=1,nin
        read(36,*,end=48,err=48) rxx(i),vae(i)
      enddo
   48 close(36)
      if(i .eq. nin+1 .and. nin .eq. mmax) then
        if(tspin) then
          write(iu,'(a/,a/,a/,a)') 
     1      '%fhipp - ERROR: file fort.36 is present and initiates', 
     1      ' non self-consistent calculation. This is not implemented',
     1      ' for the attempted spin-polarized calculation.', 
     1      ' ACTION: modify task and/or input.'
          stop
        endif
        write(iu,'(a/,a)') 
     1    '%fhipp - WARNING: non self-consistent calculation,',
     1    ' input potential read from file fort.36 .'
        call dcpv(mx,1,mmax,vae,vi)
        tinpot = .true.
      else
        tinpot = .false.
      endif

c for frozen core calculation read core density
      nin=0
      read(18,*,end=58,err=58) nin,ekin_core
      do i=1,nin
        read(18,*,end=58,err=58) rxx(i),dc(i),dcp(i),dcpp(i)
      enddo
   58 if(nin .gt. 0) then  
        if(i .ne. nin+1) 
     1    write(ie,*) 
     1      '& fhipp - warning: inconsistent input unit fort.18'
        if( abs(r(2)-rxx(2)) .gt. 1.d-12 )
     1    write(ie,*) 
     1      '& fhipp - warning: grid from unit fort.18 incompatible'
        if(mmax .gt. nin) then
          do i=nin+1,mmax
            dc(i)=0.d0
            dcp(i)=0.d0
            dcpp(i)=0.d0
          enddo
        endif
        tcore=.false.
      endif

c all-electron calculation
c Tell sratom when to do spin polarized calculation: for iexcnow<0
c spin-polarized, for iexcnow>0 unpolarized. 
      do ispin=1,nspinmx
        dcfac=1.d0/float(nspinmx)
        do i=1,mmax
          dcs(i,ispin)=dc(i)*dcfac
          dcsp(i,ispin)=dcp(i)*dcfac
          dcspp(i,ispin)=dcpp(i)*dcfac
          vee(i,ispin)=vi(i)
        enddo
      enddo
      iexcnow=iexc
      if(nspinmx.gt.1) iexcnow=-iexcs
      
      it=itmx
      call sratom(
     1    it,irl,igr,iexcnow,tcore,z,nc,nv,n,l,m,f,e,ninu,de
     1,   mmax,r,rpk,vee,ds,dsp,dspp,dcs,dcsp,dcspp,uu
     1,   tinpot,tsic,e_sic,vorb,tkli,svkli,svm
     1,   tmgga,uup)

c=============================================================
c chen: write KS potential for future use
c=============================================================
      print *,'mmax: ',mmax
      print *,'r(mmax): ',r(mmax)
      open(file='vks.dat',unit=111,action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es18.7)')r(i),vee(i,1)
      enddo
      close(111)
      write(6,'(a)')'>>> vks.dat is written <<<'
      
      ! write valence's eigenvalues 
      open(file='eigval.dat',unit=111,action='write',form='formatted')
      do i=nc+1,nv+nc
        write(111,'(es20.12)')e(i)
      enddo
      close(111)
      write(6,'(a)')'>>> eigval.dat is written <<<'

      ! write wavefunction 
      open(file='wfval.dat',unit=111,action='write',form='unformatted')
      write(111)uu(:,nc+1:nc+nv)
      close(111)
      write(6,'(a)')'>>> wfval.dat is written <<<'

      ! write AE valence rho
      open(file='rhoval.dat',unit=111,action='write',form='unformatted')
      write(111)ds(:,1)-dcs(:,1)
      write(111)dsp(:,1)-dcsp(:,1)
      write(111)dspp(:,1)-dcspp(:,1)
      close(111)
      write(6,'(a,/)')'>>> rhoval.dat is written <<<'

      ! compute number of core electron 
      zcore=0.d0 
      do i=1,nc
        zcore=zcore+f(i)
      enddo 
      zval=z-zcore
      print *,' (chen) z:     ',z
      print *,' (chen) zcore: ',zcore
      print *,' (chen) zval:  ',zval
      open(file='sys.dat',unit=111,action='write',form='formatted')
      write(111,'(f13.6,a)')z,' atomic number '
      write(111,'(f13.6,a)')zval,' valence charge'
      close(111)

c     get valence electron       
      rho_v=ds(:,1)-dcs(:,1)
      rho_vp=dsp(:,1)-dcsp(:,1)
      rho_vpp=dspp(:,1)-dcspp(:,1)
      
      cval=0.d0
      do i=nc+1,nv+nc
        cval=cval+f(i)
      enddo 
      write(6,'(/,a,f8.4,a,/)')
     &  '*** actual valence charge cval: ',cval,' ***'
      call vestat(mmax,cval,eeel,r,rho_v,ves_v,.true.)
      call vexcor(iexc,mmax,r,rho_v,rho_vp,rho_vpp, 
     &            vxc_v,eexc,ex,ec,.true.)

      vbare_v=vee(:,1)-ves_v-vxc_v ! only unscreen valence

      open(file='vbare_valence.dat',unit=111,
     &     action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es18.7)')r(i),vbare_v(i)
      enddo
      close(111)

      ! write valence electron 
      open(file='rho_val_AE.dat',unit=111,
     &     action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es18.7)')r(i),ds(i,1)-dcs(i,1)
      enddo
      close(111)
      write(6,'(a)')'*** written valence density to rho_val_AE.dat ***'

      ! write valence electron 
      open(file='rho_core_AE.dat',unit=111,
     &     action='write',form='formatted')
      do i=1,mmax
        write(111,'(2es18.7)')r(i),dcs(i,1)
      enddo
      close(111)
      write(6,'(a)')'*** written core density to rho_core_AE.dat ***'

      ! find Coulomb tail
      epstail=1.d-3
      do i=mmax,1,-1
        if(abs(vbare_v(i)+zval/r(i)) .gt. epstail) then
          write(6,'(a,f5.2,a,/)') 
     & '*** Coulomb tail (valence) starts at ',r(i),' bohr ***'
          exit 
        endif 
      enddo
c=============================================================

c compute total density to retain compatibility with existing code,
c also copy full potential for pseudopotential construction
      do i=1,mmax
        rho(i)=ds(i,1)
        rhop(i)=dsp(i,1)
        rhopp(i)=dspp(i,1)
        dc(i)=dcs(i,1)
        dcp(i)=dcsp(i,1)
        dcpp(i)=dcspp(i,1)
        vi(i)=vee(i,1)
      enddo
      do ispin=2,nspinmx
        do i=1,mmax
          rho(i)=rho(i)+ds(i,ispin)
          rhop(i)=rhop(i)+dsp(i,ispin)
          rhopp(i)=rhopp(i)+dspp(i,ispin)
          dc(i)=dc(i)+dcs(i,ispin)
          dcp(i)=dcp(i)+dcsp(i,ispin)
          dcpp(i)=dcpp(i)+dcspp(i,ispin)
        enddo
      enddo

      if(it.gt.itmx) 
     1   write(ie,*) '& fhipp - warning: potential not converged'

c some tests of the electron density
c test: 1 = Integral(rho)/charge
      tc = fmom(0,mmax,al,1.d0,r,rho)
      if(igr .eq. 1) then      
c test: 1 = -Integral(r*rho')/(3*charge)
        t1 = -fmom(1,mmax,al,sf,r,rhop)/3.d0
c test: 1 = Integral(r**2 * rho'')/(12*charge)
        t2 = fmom(2,mmax,al,sf,r,rhopp)/12.d0
c meta gga: evaluate kinetic energy density
        if(tmgga) then
        do i=1,mmax
          stat(i) = 0.d0
        enddo
        do i=nc+1,nc+nv
          do j=1,ninu(i)
            stat(j)=stat(j)+f(i)*stat_orb(l(i),al,r(j),uu(j,i),uup(j,i))
          enddo
        enddo
        do i=1,mmax
          stat_core(i) = 0.d0
        enddo
        if(tcore) then
          do i=1,nc
            do j=1,ninu(i)
            stat_core(j)=stat_core(j)
     1         +f(i)*stat_orb(l(i),al,r(j),uu(j,i),uup(j,i))
            enddo
          enddo
        else
          do i=1,mmax
            read(16,*,end=81,err=81) rxx(i),stat_core(i)
          enddo
  81      write(ie,*) '& fhipp - read core kinetic energy density',i-1
        endif
        ek_core_from_stat = fmom(0,mmax,al,1.d0,r,stat)
        do i=1,mmax			! add core and valence parts
          stat(i) = stat(i) + stat_core(i)
        enddo
        ek_from_stat      = fmom(0,mmax,al,1.d0,r,stat)
        ek_core_from_stat = ek_from_stat - ek_core_from_stat
        endif
      endif

c total energy components of all-electron atom
c Coulomb energy, xc-potential energy, etc
      ecou=-z*fmom(-1,mmax,al,1.d0,r,rho)
      call vestat(mmax,sf,eeel,r,rho,ves,.true.)
      ex=0.d0
      ec=0.d0
      if(nspinmx.eq.1) then
        call vexcor(iexc,mmax,r,ds(1,1),dsp(1,1),dspp(1,1),
     1              vxc(1,1),evxc,ex,ec,.true.) 
      else
        call vexcos(iexcs,mmax,mmax,r,ds,dsp,dspp,vxc,eexc,ex,ec,.true.)
      endif
      evxc=0.d0
      do ispin=1,nspinmx
        evxc=evxc+dmelm(mmax,al,r,vxc(1,ispin),ds(1,ispin))
      enddo
      if(tkli) then
       do i=1,nc+nv
        ex=ex+f(i)*svm(i)
       enddo
       ex=0.5d0*ex
       evxc=dmelm(mmax,al,r,vi,rho)-eeel-ecou
      endif
c effective potential energy and kinetic energy
      epot=0.d0
      do ispin=1,nspinmx
        epot=epot+dmelm(mmax,al,r,vee(1,ispin),ds(1,ispin))
      enddo
c core 
      epot_core=0.d0
      do ispin=1,nspinmx
        epot_core=epot_core+dmelm(mmax,al,r,vee(1,ispin),dcs(1,ispin))
      enddo
      esic_core=0.d0
      if(tcore) then
        ekin_core=0.d0
        do ispin=1,nspinmx
          nstart=(nc+nv)*(ispin-1)+1
          norbnow=nstart+nc-1
c debug
!d         write(ie,*)'%fhipp - core nstart=',nstart,' norbnow=',norbnow
!d         tc = fmom(0,mmax,al,1.d0,r,dcs(1,ispin))
!d         write(ie,*)'%fhipp - dcs electrons=',tc
!d         tc = fmom(0,mmax,al,1.d0,r,dc)
!d         write(ie,*)'%fhipp - dc electrons=',tc
!d         tc = fmom(0,mmax,al,1.d0,r,ds(1,ispin))
!d         write(ie,*)'%fhipp - ds electrons=',tc
!d         tc = fmom(0,mmax,al,1.d0,r,rho)
!d         write(ie,*)'%fhipp - rho electrons=',tc
          do i=nstart,norbnow
            ekin_core=ekin_core+f(i)*e(i)
            if(tsic) then 
              ekin_core=ekin_core
     1          +f(i)*gltfmv(mmax,al,r,uu(1,i),vorb(1,i),uu(1,i))
              esic_core=esic_core+f(i)*e_sic(i)
            endif
          enddo
        enddo
        ekin_core=ekin_core-epot_core
      endif
c valence
      ekin=0.d0
      esic=0.d0
      eeig=0.d0
      do ispin=1,nspinmx
        nstart=nc+(nv+nc)*(ispin-1)+1
        norbnow=nstart+nv-1
        do i=nstart,norbnow
          ekin=ekin+f(i)*e(i)
          eeig=eeig+f(i)*e(i)
          if(tsic) then 
            ekin=ekin+f(i)*gltfmv(mmax,al,r,uu(1,i),vorb(1,i),uu(1,i))
            esic=esic+f(i)*e_sic(i)
          endif
        enddo
      enddo
      ekin=(ekin_core+ekin)-epot+epot_core
      esic=esic_core+esic
      if(tcore) then
        norbnow=(nc+nv)*nspinmx
        eeig=dspv(ms,1,norbnow,f,e)
      endif
      etot=ekin+0.5d0*eeel+ecou+ex+ec+esic
      if(iexc.eq.0) etot=ekin+ecou

c output for all-electron calculation
      call labelmap(max(int(abs(z)),1),symz,iexc,symxc)
      write(iu,619) srev 
      write(iu,'(a30,2x,a8)')   'chemical symbol', symz
      write(iu,'(a30,2x,f5.2)') 'nuclear charge', z      
      write(iu,'(a30,2x,f5.2)') 'total charge', z-sf
      write(iu,'(a30,2x,i2)')   'number of core states', nc
      write(iu,'(a30,2x,i2)')   'number of valence states', nv
      write(iu,'(a30,2x,i2,2x,a40)') 
     1    'exchange-correlation model', iexc,symxc 
      if(tspin) write(iu,'(a30)')'spin-polarized calculation'
      if(tsic) write(iu,'(a30)')'self-interaction correction'
      if(tnrl) then
         write(iu,'(a30)')    'non-relativistic mode'
      else
         write(iu,'(a30)')    'scalar-relativistic mode'
      endif
      if(nc.ge.1 .and. .not.tcore) write(iu,'(a30)') 'frozen-core mode'
      write(iu,'(a30,2x,i4,f13.6,2x,e13.6)') 
     1   'parameters radial mesh',mmax,amesh,r(1)
      write(iu,'(/,a32)') '=== all-electron atom ==='
      write(iu,620)
      nstart=nc+1
      if(tcore) nstart=1
      ishift=nc+nv
      do i=nstart,nc+nv
        if(nspinmx.eq.1) then
          write(iu,630) i,n(i),l(i),f(i),e(i)*ry2,e(i)
        else
          write(iu,635) i,n(i),l(i),f(i),e(i)*ry2,e(i)        
          j=i+ishift
          write(iu,636) j,n(j),l(j),f(j),e(j)*ry2,e(j)        
        endif
      enddo

cmf special
c calculate spin-orbit splitting
      if(.false.) then
      do i=nc+1,nc+nv
        if(l(i) > 0) then
          do ir=1,mmax
            rxx(ir)=0.d0
          enddo
          do ir=1,ninu(i)
            rxx(ir)=(uu(ir,i)/r(ir))**2
          enddo
          call spinorbit(mmax,n(i),l(i),ry2,e(i),r,rxx,vi)
        endif
      enddo
      do i=100,nc+nv
        do j=1,i
          do ir=1,mmax
            rxx(ir)=0.d0
          enddo
          if(l(i) == l(j)) then
            do ir=1,ninu(i)
              rxx(ir)=uu(ir,i)*uu(ir,j)
            enddo
            e(ms)=fmom(-2,ninu(i),al,1.d0,r,rxx)
c debug
c           write(ie,*)'%fhipp - SO split: n(i) l(i) n(j) l(j) e',
c    1        n(i),l(i),n(j),l(j),e(ms)
          endif
       enddo
       enddo
       e(ms)=0.d0
      endif
cmf

      write(iu,631) 
c
      write(iu,608)'total energy',etot
      write(iu,608)'kinetic energy',ekin
      write(iu,608)'orbital energy',eeig
      write(iu,608)'coulomb energy',ecou
      write(iu,608)'hartree energy',0.5*eeel
      write(iu,608)'exchange-correlation energy',(ex+ec)
      if(tsic)
     1  write(iu,608)'self-interaction correction',etot_sic
      write(iu,608)'xc potential energy',evxc
      if(tsic)
     1  write(iu,608)'self-interaction correction',epot_sic
      write(iu,'(a30,2x,i16,2x,a12,1pe9.1)')
     1  'number of iterations', it, 'convergence', de
      write(iu,608) 'core kinetic energy',ekin_core
      write(iu,608) 'core orbital energy',ekin_core + epot_core
      write(iu,608) 'core potential energy',epot_core
      if(tmgga) then
        write(iu,608) 'g kinetic energy',ek_from_stat
        write(iu,608) 'g core kinetic energy',ek_core_from_stat
      endif
      write(iu,608) 'integrated density',tc
      if(igr .eq. 1) then
        write(iu,608) '... 1st derivative test 1 = ',t1
        write(iu,608) '... 2nd derivative test 1 = ',t2
      endif
  619 format(/'fhi pseudopotential tool fhipp - version ',a10,/)
  620 format(/'<        n     l      occupation  eigenvalue(eV)'/)
  630 format('< ',i2,4x,2(i2,4x),f10.4,f18.4,1x,f18.4)
  631 format(/34x,'(Hartree a.u.)')
  635 format('< ',i2,4x,i2,'_u',2x,i2,4x,f10.4,2x,f18.4,4x,f18.4)
  636 format('< ',3x,i2,1x,i2,'_d',2x,i2,4x,f10.4,2x,f18.4,4x,f18.4)
  608 format(a30,2x,f16.5)

c spin-polarized calculations implemented up to here
      if(tspin) then
        write(iu,'(/a)') 'fhipp - spin-polarized all-electron atom done'
        stop 
      endif

c all-electron wave functions
      if(tcore) then
        call out38(mx,1,nc+nv,ninu,n,l,r,uu)
      else
        call out38(mx,nc+1,nc+nv,ninu,n,l,r,uu)
      endif

      write(iu,*) 
      write(iu,*) 'fhipp - all-electron atom done'
      if(.not. tdopsp) stop 

c ==== all-electron atom done ===

c all-electron effective potential
      write(37,'(f10.5,1x,i4,1x,i4,1x,i4)') z,mmax,iexc,irl
      do i=1,mmax
        write(37,'(e13.6,1x,e21.14)') r(i),vi(i)
      enddo
      close(37)

c output full core charge density
      if(tcore)then
        write(19,'(i6,1x,e21.14)') ninu(nc),ekin_core
        do i=1,ninu(nc)
          write(19,'(e21.14,3(1x,e21.14))') r(i),dc(i),dcp(i),dcpp(i)
          write(17,'(e21.14,(1x,e21.14))')  r(i),stat_core(i)
        enddo

        if(tmgga) then
          open(10,file='alpha.rho')
          open(11,file='alpha.rho_kin')
          open(12,file='alpha.s_eff')	! effective scaled gradient
          open(13,file='alpha.t_ratio')	! effective t_w/t_s
          open(14,file='alpha.rho_vW')	! effective t_w/t_s
          do i=1,mmax
            write(10,'(e13.6,1x,e15.8)') r(i),rho(i)
            write(11,'(e13.6,1x,e15.8)') r(i),stat(i)+1e-12
            vara = fx_xvar(rho(i),rhop(i),stat(i))
            write(12,'(e13.6,1x,e15.8)') r(i),vara+1e-12
            write(13,'(e13.6,1x,e15.8)') 
     1        r(i),rhop(i)**2/( 8d0* rho(i)*stat(i) )+1e-12
            write(14,'(e13.6,1x,e15.8)') r(i),0.125*rhop(i)**2/rho(i)
          enddo
          write(10,*)'&'
          write(11,*)'&'
          write(12,*)'&'
          write(13,*)'&'
          write(14,*)'&'
          do i=1,mmax
            write(10,'(e13.6,1x,e15.8)') r(i),rho(i)-dc(i)
            write(11,'(e13.6,1x,e15.8)') r(i),stat(i)-stat_core(i)+1e-12
            vara = rho(i) - dc(i)
            varb = rhop(i) - dcp(i)
            varc = stat(i) - stat_core(i)
            vard = fx_xvar(vara,varb,varc)
            write(12,'(e13.6,1x,e15.8)') r(i),vard+1e-12
            write(13,'(e13.6,1x,e15.8)') r(i),(rhop(i)-dcp(i))**2/
     1        ( 8d0* (rho(i)-dc(i)) * (stat(i)-stat_core(i)) )+1e-12
            write(14,'(e13.6,1x,e15.8)') r(i),
     1        0.125*(rhop(i)-dcp(i))**2/(rho(i)-dc(i))
          enddo
          write(10,*)'&'
          write(11,*)'&'
          write(12,*)'&'
          write(13,*)'&'
          write(14,*)'&'
          do i=1,mmax
            write(10,'(e13.6,1x,e15.8)') r(i),dc(i)
            write(11,'(e13.6,1x,e15.8)') r(i),stat_core(i)+1e-12
            vara = fx_xvar(dc(i),dcp(i),stat_core(i))
            write(12,'(e13.6,1x,e15.8)') r(i),vara+1e-12
            write(13,'(e13.6,1x,e15.8)') r(i),
     1         dcp(i)**2/( 8d0* dc(i)*stat_core(i) )+1e-12
            write(14,'(e13.6,1x,e15.8)') r(i),
     1        0.125*dcp(i)**2/dc(i)
          enddo
          close(10)
          close(11)
          close(12)
          close(13)
          close(14)
        endif
        close(17)
        close(19)
      endif

c prepare default pseudopotential input
      read(22,*,end=65,err=65) ltmx,spptype
      if(spptype(1:1) .eq. 'h' .or. spptype(1:1) .eq. 'H') then
        write(iu,'(/4x,a20,2x,a6)') '=== HAMANN mode === ',spptype
        i_pp_type(1)=1
      else if(spptype(1:1) .eq. 't' .or. spptype(1:1) .eq. 'T') then
        write(iu,'(/4x,a30,2x,a6)') 
     1    '=== TROULLIER-MARTINS mode ===',spptype
        i_pp_type(1)=2
      else
        write(ie,*) '& fhipp - stop: bad pseudopotential type',spptype
        stop
      endif
      i_pp_type_def=i_pp_type(1)

      do i=1,ms
        np(i)=0
        ep(i)=0.0d0
        fp(i)=0.0d0
        rc(i)=0.0d0
        inode(i)=1
        ipeak(i)=1
        i_pp_type(i)=i_pp_type_def
      enddo
      
c defaults for Hamann case (i_pp_type = 1), cf Ref [1]
      csc=1.9d0
      csc1=4.0d0
      vsc(1)=0.60d0
      vsc1(1)=0.40d0
c defaults for Troullier-Martins case (i_pp_type = 2)
      vsc(2)=1.55*sqrt(7.d0/z)
      vsc1(2)=1.60*sqrt(z/7.d0)

c Hamann's choice if valence state is unoccupied and unbound
      rcmax=0.0d0
      do i=1,nc
        l1=l(i)+1
        np(l1)=max(np(l1),n(i))
        if(n(i) .eq. l(i)+1) then
          rc_d(l1,1)=csc1*rpk(i)
        else
          rc_d(l1,1)=csc*rpk(i)
        endif
        rcmax=max(rcmax,rc(l1))
        fp(l1)=0.d0
      enddo

c for occupied valence states
      rc_dmax(1)=0.0d0
      rc_dmax(2)=0.0d0
      emax=-100.0d0
      zv=0.0d0
      if(nv .eq. 0) then
        emax = 0.0d0
      else
        do i=nc+nv,nc+1,-1

          zv=zv+f(i)
          emax=max(emax,e(i))
          l1=l(i)+1
          np(l1)=n(i)
          fp(l1)=f(i)
          ep(l1)=e(i)
          if(tsic) call dcpv(mx,1,mmax,vorb(1,i),vsic(1,l1))

c compute default radius
          do itype=1,2
            rc_d(l1,itype)=vsc(itype)*rpk(i)
            if(n(i) .eq. l1) rc_d(l1,itype)=vsc1(itype)*rpk(i)
            if(itype.eq.2 .and. i_pp_type(l1).eq.2) 
     1        rc_d(l1,i_pp_type(l1))=defrtm(np(l1),l1-1,int(z))
            rc_dmax(itype)=max(rc_dmax(itype),rc_d(l1,itype))
          enddo
c find outermost node
          if(n(i) .gt. l1) then
            do j=mmax-1,1,-1
              if(uu(j,i)*uu(j+1,i) .lt. 0.0) goto 700
              if(rpk(i) .le. r(j)) ipeak(l1)=j
            enddo
  700       inode(l1)=j
          else
            do j=mmax-1,1,-1
              if(rpk(i) .ge. r(j)) goto 710
            enddo
  710       ipeak(l1)=j
            inode(l1)=1
          endif
        enddo
      endif

c unoccupied valence states & default assignment 
      lmax=max(ltmx,lmax)
      if(lmax .lt. 0 .or. lmax .gt. 4) 
     1  stop 'fhipp - lmax out of range (0 ... 4)'
      do l1=1,lmax+1
        t_pp_mode(l1)=.false.
        i_pp_type(l1)=i_pp_type(1)
        if(fp(l1) .eq. 0.0 .and. ep(l1) .ge. 0.0) then
          t_pp_mode(l1)=.true.
          ep(l1)=emax
          np(l1)=max(np(l1)+1,l1)
          rc_d(l1,1)=max(rc_dmax(1),rc_d(l1,1))
          rc_d(l1,2)=max(rc_dmax(2),rc_d(l1,2))
        endif
        if(fp(l1).eq.0.0 .and. i_pp_type(l1).eq.2) 
     1    rc_d(l1,i_pp_type(l1))=defrtm(np(l1),l1-1,int(z))

        rc(l1)=rc_d(l1,i_pp_type(l1))
c in case the cutoff radius is not tabulated (troullier-martins) it is computed,
c these radii may need adjustment for transferability/representability reasons
        if(rc(l1) .lt. r(1)) then
          write(ie,*)
     1      '& fhipp - warning: no tabulated default core radius l',l1-1
          rc(l1)=rc_d(l1,i_pp_type(1))/vsc(1)
          if(np(l1) .eq. l1) 
     1      rc(l1)=sqrt(rc_d(l1,1)/vsc1(1)*(l1*(l1-1)+0.5d0))
        endif
      enddo

c print info on valence wavefunctions 
      write(iu,633) 
      do l1=1,lmax+1
        write(iu,'(a,i2,1x,i2,10x,3(5x,f6.3))')
     1    'x', l1-1, np(l1), r(inode(l1)), r(ipeak(l1)), rc(l1)
      enddo
  633 format(/'  l  n     radius:     node      peak       default core'
     1      )

c defaults done

c******************************************************************
c alternative input overriding default values
c
c lt        l value for following parameter changes
c rct       if non-zero, new core radius for angular momentum lt
c et        if non-zero, new energy for lt (effective for unoccupied
c           states only) in eV
c spptype   pseudopotential type 
c*******************************************************************
      spptype='-'
      et=0.0
      rct=0.0
      do i=1,lmax+1
        tinopt(i)=.false.
      enddo
   50 read(22,*,end=60,err=55) lt,rct,et,spptype
        if(lt .gt. lmax) write(ie,*)
     &    '& fhipp - warning: optional input not applied: lmax < lt',lt
        if(lt .gt. lmax) goto 50
        l1=lt+1
        if(lt .lt. 0 .or. lt .gt. 4) 
     &    stop 'fhipp - optional data: l out of range (0 ... 4)'
        if(spptype(1:1) .eq. 'h' .or. spptype(1:1) .eq. 'H') then
          i_pp_type(l1)=1
        else 
     &  if(spptype(1:1) .eq. 't' .or. spptype(1:1) .eq. 'T') then
          i_pp_type(l1)=2
        else
          i_pp_type(l1)=i_pp_type_def
        endif
        if(rct .gt. r(1)) then
          rc(l1)=abs(rct)
          rctopt(l1)=rc(l1)
        else
          rc(l1)=rc_d(l1,i_pp_type(l1))
          rctopt(l1)=0.0d0
        endif
        if(abs(et) .eq. 0.0) then
          ep(l1)=emax
          ectopt(l1)=0.0d0
        else
          t_pp_mode(l1)=.true.
          ep(l1)=et/ry2
          ectopt(l1)=ep(l1)
        endif 
        ltopt(l1)=lt
        sppopt(l1)=spp(i_pp_type(l1))
        tinopt(l1)=.true.
c debug
!d       write(ie,*) '%fhipp - option',
!d    &    l1-1,np(l1),fp(l1),rc(l1),ep(l1),i_pp_type(l1)
      goto 50
   55 write(ie,*) '& fhipp - warning: failed reading optional data'
   60 continue

c verification of cutoff radii
      do l1=1,lmax+1
        if( rc(l1) .lt. r(1) .or. rc(l1) .lt. r(inode(l1)) ) then
          write(ie,*) '& fhipp - stop: cutoff radius too small l',l1-1
          stop
        endif
      enddo

c number of core & valence electrons
      zc=sf-zv

c valence density
      do i=1,mmax
        rho(i)=rho(i)-dc(i)
        rhop(i)=rhop(i)-dcp(i)
        rhopp(i)=rhopp(i)-dcpp(i)
      enddo

c pseudopotentials
      call ncpp(iu,irl,iexc,z,zc,zv,lmax,t_pp_mode,i_pp_type 
     1,  np,ep,fp,rc,vi,mmax,r
     1,  rnlc,dc,dcp,dcpp,rho,rhop,rhopp
     1,  tsic,vsic,ups
     1,  nv,n(nc+1),l(nc+1),e(nc+1),f(nc+1),svkli(nc+1))

      write(iu,*) 
      write(iu,*) 'fhipp - done for input @'
      write(iu,*) 

c output input data
      write(iu,'(a1,1x,f5.2,3i3,1x,e9.2,1x,a19)') '@',z,nc,nv,iexc,rnlc,
     &  ': z nc nv iexc rnlc'
      write(iu,'(a1,5x,i1,2x,i1,2x,f5.2,9x,a7)') '@',n(1),l(1),f(1),
     &  ': n l f'
      do i=2,nc+nv
        write(iu,'(a1,5x,i1,2x,i1,2x,f5.2)') '@',n(i),l(i),f(i)
      enddo
      write(iu,'(a1,1x,i1,2x,a1,19x,a14)') '@',ltmx,spp(i_pp_type_def),
     &  ': ltmx spptype'
      i=0
      do l1=1,lmax+1
        if( tinopt(l1) ) then
          if(i .eq. 0) then 
            i=1
            write(iu,'(a1,1x,i1,2x,f5.2,2x,e9.2,2x,a1,2x,a19)')
     &        '@',ltopt(l1),rctopt(l1),ectopt(l1)*ry2,sppopt(l1),
     &        ': lt rct et spptype'
          else
            write(iu,'(a1,1x,i1,2x,f5.2,2x,e9.2,2x,a1)') 
     &        '@',ltopt(l1),rctopt(l1),ectopt(l1)*ry2,sppopt(l1)
          endif
        else
          spptype=spp(i_pp_type_def)
          write(iu,'(a1,1x,i1,2x,f5.2,2x,e9.2,2x,a1)')
     &        '@',l1,rc(l1),ep(l1)*ry2,spptype(1:1)
        endif
      enddo
c
c close files
c
      if(tkli) close(80)
      close(iu)
      stop
  65  stop 'fhipp - verify input data file fort.22'
      end
c
