c Header:$
cmf 27-06-03 added spin polarization for selected xc functionals 
cmf 04-09-02 added kli exchange functional
c-----------------------------------------------------------------------
c self-consistent (pseudo-) atom calculation
c using log mesh 
c sphericalized, 1st & 2nd radial derivatives 
c
c M. Fuchs, FHI Berlin, 08-1993
c-----------------------------------------------------------------------
c
! Inputs
!  vi:  initial guess for xc+hartree potential
!  vps: semilocal ionic pseudopotential 
!
! Output:
!  vi: final xc+hartree potential
!
      subroutine psatom_etot(it,igr,iexc,z,nc,nv,
     1    n,l,f,m,e,mup,de,mmax,r,vi
     1,   vps,d,dp,dpp,dc,dcp,dcpp,uu,tkli,svm,svm_roks,etot)
c
      implicit real*8 (a-h,o-z)
      logical tconv,tkli
c
      include 'parameter.h'
c
      dimension n(ms),l(ms),f(ms),e(ms),rpk(ms),ves(mx),vps(mx,ms)
     &,         vpw(mx,ms),uu(mx,ms),wa(mx),wb(mx)
      integer m(ms)
      real*8  feps
      real*8  r(mx),r1i(mx),r2i(mx),r3i(mx),r4i(mx),vi(mx,ms)
      real*8  u(mx),up(mx),upp(mx)
      real*8  rho(mx),rhop(mx),rhopp(mx),dc(mx),dcp(mx),dcpp(mx)
      real*8  vo(mx,ms),vi1(mx,ms),vo1(mx,ms),vxc(mx,ms),vee(mx,ms)
      real*8  rhonow(mx),rhonowp(mx),rhonowpp(mx),etot,rhos(mx,2)

c klix
      logical  tstart,tkli_first
      integer  iexckli,iorb,ir,norb_kli
      integer  mp(ms)
      real*8   fp(ms),dens(mx),svkli(ms),svm(ms),svm_roks(ms)
      real*8   svm_tmp(ms),vckli(mx),vsl(mx)

c meta gga
      integer ninu(ms)
c     real*8  uup(mx,ms)
      real*8  stat,stat_core
      common/use_mgga/stat(mx),stat_core(mx)

c spin-polarization variables
      logical tspin
      real*8  ex_s,ec_s,abc,sf_up,sf_dn
      real*8  d(mx,2),dp(mx,2),dpp(mx,2)
      real*8  f_up(ms),f_dn(ms),ns(ms),ls(ms)
      real*8  uup(mx,ms),v_s(mx,2)

c initializations
      feps=1.d-9
      itmax=it
      al=0.1d0*log(r(11)/r(1))
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
        stat(i)=0.d0
        vsl(i)=0.d0
        r1i(i)=1.d0/r(i)
        r2i(i)=r1i(i)*r1i(i)
        r3i(i)=r1i(i)*r2i(i)
        r4i(i)=r2i(i)*r2i(i)
        rhop(i)=0.d0
        rhopp(i)=0.d0
      enddo
c spin initialize
      nspinmx=1
      if(iexc.lt.0) nspinmx=2
      iexc=abs(iexc)
c kli initialize
      iexckli=iexc
      if(iexc.eq.12) iexckli=0
      tkli_first=.true.

c initial potential
      do j=1,nc+nv
        ispin=m(j)
        do i=1,mmax
          vpw(i,j)=vps(i,j)+vi(i,ispin)
        enddo
      enddo

c scale core density to core spin density in spin polarized case
      if(nspinmx.gt.1) then
        do ir=1,mmax
          dc(ir)=0.5d0*dc(ir)
          dcp(ir)=0.5d0*dcp(ir)
          dcpp(ir)=0.5d0*dcpp(ir)
        enddo
      endif
          
c
c bl is the mxing factor for the sc-potential 
c
      it=0
      blmin=0.40d0
      blmax=0.85d0
      als=al*al
      sf=0.d0
      do i=1,nc+nv
        sf=sf+f(i)
      enddo
c
c
c return point for self-consistency loop
c
      do it=1,itmax
        tconv=.true.
        mup=0
        bl=min(blmax,blmin*sqrt(float(it)))
c
c initialize for new iteration
c
        rhos=0.d0 
        do ispin=1,nspinmx
          do j=1,mmax
            rho(j)=0.d0
            d(j,ispin)=0.0d0
            dp(j,ispin)=0.0d0
            dpp(j,ispin)=0.0d0
          enddo
        enddo
c
c solve for pseudo states, rho() is the valence density
c
        do i=1,nc+nv

            if(f(i) .gt. feps .or. tconv) then
              et=e(i)
              call dftseq(8,z,mmax,r,n(i),l(i),1.d0,vpw(1,i),wa,wb,
     1            nin,mch,uld,et,u,up,upp)
              if(e(i) .ne. et) tconv=.false.
              e(i)=et
            endif
c debug
!d           write(ie,*) 'it state',it,i,n(i),l(i),m(i),e(i)
            mup=max(mup,nin)
            ninu(i) = nin
c
c
c accumulate charge & radial derivatives, rhop & rhopp are calculated as
c derivatives w.r.t. r, while up & upp are in terms of the logarithmic
c radial variable
c
ccccccccccccccccccccccc chen: update wave function ccccccccccccccccccccc
            !if(tkli .or. tconv) then
            if(tkli .or. tconv .or. itmax==1) then
              do j=1,nin
                uu(j,i)=u(j)
                uup(j,i)=up(j)
              enddo
            endif
c orbital density and update of spin densities
            ispin=min(nspinmx,m(i))
            do j=1,nin
              rhonow(j)=u(j)*u(j)
              d(j,ispin)=d(j,ispin)+f(i)*rhonow(j)
            enddo
            if(igr.eq.1) then
              do j=1,nin
                rhonowp(j)=2.d0*(up(j)/al-u(j))*u(j)
                rhonowpp(j)=2.d0*( (up(j)/al-u(j))*(up(j)/al
     1                 -3.d0*u(j)) +(upp(j)-al*up(j))*u(j)/als )
                dp(j,ispin)=dp(j,ispin)+f(i)*rhonowp(j)
                dpp(j,ispin)=dpp(j,ispin)+f(i)*rhonowpp(j)
              enddo
            endif

        enddo
c
c densities d() core+valence, rho() valence only
        do ispin=1,nspinmx
          do j=1,mmax
            d(j,ispin)=d(j,ispin)*r2i(j)
            rho(j)=rho(j)+d(j,ispin)
            rhos(j,ispin)=rhos(j,ispin)+d(j,ispin)
            d(j,ispin)=d(j,ispin)+dc(j)
            dp(j,ispin)=dp(j,ispin)*r3i(j)+dcp(j)
            dpp(j,ispin)=dpp(j,ispin)*r4i(j)+dcpp(j)
          enddo
        enddo
c
c output change of electronic part of the effective potential
c vestat(sf... w/out z/r, vestat(sf-z... w/ z/r potential
        call vestat(mmax,sf,eeel,r,rho,ves,tconv)

        if(tkli) then
          if(it < 10) then
            call vexcor(8,mup,r,d(1,1),dp(1,1),dpp(1,1),
     1                  vxc(1,1),eexc,ex,ec,.false.)
          else

            if(tkli_first) then
              do iorb=1,nc+nv
                fp(iorb)=0.5d0*f(iorb)
                mp(iorb)=1
!d               write(iu,*) '%psatom: iorb',iorb,f(iorb),fp(iorb)
!d               write(iu,*) '             ',iorb,n(iorb)
!d               write(iu,*) '             ',iorb,l(iorb)
              enddo
              norb_kli=0
              do iorb=1,nc+nv
                if(fp(iorb) .gt. 1.d-12) norb_kli=norb_kli+1
              enddo
!d             write(iu,'(/a,1x,i2)')
!d    1          ' %psatom: kli - number of occupied orbitals = ',
!d    1          norb_kli
            endif
            do i=1,mmax
              dens(i)=0.5d0*rho(i)
            enddo
 
            call vklix(norb_kli,n,l,mp,fp,mup,r,uu,
     1                 dens,vxc(1,1),vsl,svkli,svm,1,
     1                 tkli_first,tconv,tconv)
            if(iexckli.ne.0) then
              call vexcor(iexckli,mup,r,d(1,1),dp(1,1),dpp(1,1),
     1                    vckli,eexc,ex,ec,.false.)
!d             write(iu,*) '%psatom: iexckli',vckli(1),vckli(100)
              do ir=1,mmax
                vxc(ir,1)=vxc(ir,1)+vckli(ir)
              enddo
            endif
 
c do restricted open shell exchange energy -> roks
            if(tconv) then
              norb_kli=0
              do iorb=1,nc+nv
               fp(iorb)=min(dble(2*l(iorb)+1),f(iorb))
               if(fp(iorb) .gt. 1.d-12) norb_kli=norb_kli+1
              enddo
              call vklix(norb_kli,n,l,mp,fp,mup,r,uu,
     1                 dens,vxc(1,1),vsl,svkli,svm_tmp,1,
     1                 tkli_first,tconv,tconv)
              do iorb=1,norb_kli
               svm_roks(iorb)=fp(iorb)/f(iorb)*svm_tmp(iorb)
              enddo
              norb_kli=0
              do iorb=1,nc+nv
               fp(iorb)=max(0.d0,f(iorb)-fp(iorb))
               if(fp(iorb) .gt. 1.d-12) norb_kli=norb_kli+1
              enddo
              call vklix(norb_kli,n,l,mp,fp,mup,r,uu,
     1                 dens,vxc(1,1),vsl,svkli,svm_tmp,1,
     1                 tkli_first,tconv,tconv)
              do iorb=1,norb_kli
               svm_roks(iorb)=fp(iorb)/f(iorb)*svm_tmp(iorb)
     1                        +svm_roks(iorb)
              enddo
            endif

            if(tkli_first) tkli_first=.false.
          endif

c test case for 2 electrons only
c        do i=1,mmax
c         ves(i)=ves(i)*0.5d0
c        enddo

        else
 
c lda or gga xc potential
          if(nspinmx.eq.1) then
c spin unpolarized
            call vexcor(iexc,mmax,r,d(1,1),dp(1,1),dpp(1,1),vxc(1,1),
     1                  eexc,ex,ec,tconv)
          else
c spin polarized
            call vexcos(iexc,mup,mmax,r,d,dp,dpp,vxc,eexc,ex,ec,tconv)
          endif

        endif

        do ispin=1,nspinmx
          call dadv(mx,1,mmax,ves,vxc(1,ispin),vo(1,ispin))
          call anderson(bl,it,mmax,r,vo1(1,ispin),vo(1,ispin),
     1                  vi1(1,ispin),vee(1,ispin))
        enddo
          
c
c updating the potentials
c
        do j=1,nc+nv
          ispin=m(j)
          do i=1,mmax
            vpw(i,j)=vps(i,j)+vee(i,ispin)
          enddo
        enddo
c
c total energy
c
        eeig=0.d0
        do i=1,nc+nv
          eeig=eeig+f(i)*e(i)
        enddo
        de=abs((eeig-eeig_old)/eeig)
        eeig_old=eeig
c      
        if(tconv) goto 24
      enddo
  24  continue


c=================================================
c total energy 
c local pseudopotential energy
      evps_loc=dmelm(mmax,al,r,vps(1,1),rho)
c screened pseudopotential energy
      epsp=evps_loc
      do ispin=1,nspinmx
        epsp=epsp+dmelm(mmax,al,r,vee(1,ispin),rhos(1,ispin))
      enddo
c kinetic energy from eigenvalues
      ekin=eeig-epsp
      !print *,'ke:',ekin
      !print *,'local:',evps_loc
      !print *,'hartree:',.5*eeel
      !print *,'ex+ec: ',ex+ec
      etot=ekin+.5*eeel+evps_loc+ex+ec
c===================================================      


c output (screening) electronic potential
      do ispin=1,nspinmx
        do i=1,mmax
          vi(i,ispin)=vee(i,ispin)
          d(i,ispin)=d(i,ispin)-dc(i)
          dp(i,ispin)=dp(i,ispin)-dcp(i)
          dpp(i,ispin)=dpp(i,ispin)-dcpp(i)
        enddo
      enddo

c restore core total density
      if(nspinmx.gt.1) then
        do ir=1,mmax
          dc(ir)=2.d0*dc(ir)
          dcp(ir)=2.d0*dcp(ir)
          dcpp(ir)=2.d0*dcpp(ir)
        enddo
      endif

c meta gga kinetic energy density
      do i=1,nc+nv
        do j=1,ninu(i)
          stat(j)=stat(j)+f(i)*stat_orb(l(i),al,r(j),uu(j,i),uup(j,i))
        enddo
      enddo

      return
      end
c
