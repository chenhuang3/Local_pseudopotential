c$Header:$
cmf 27-06-03 added spin polarization for selected xc functionals
cmf 04-09-02 added kli exchange functional
c**********************************************************************
c self-consistent scalar-relativistic spherical all-electron atom 
c 1st & 2nd radial derivatives 
c frozen-core
c non-spinpolarized
c self-interaction correction
c
c Martin Fuchs, FHI der MPG, Berlin, 11-1996
c***********************************************************************
      subroutine sratom(
     1    it,irl,igr,iexc,tcore,z,nc,nv,n,l,m,f,e,ninu,de
     1,   mmax,r,rpk,vee,d,dp,dpp,dcs,dcsp,dcspp,uu
     1,   tinpot,tsic,esic,vsic,tkli,svkli,svm
     1,   tmgga,uup)
c
      implicit none
      include  'default.h'
      include  'parameter.h'

      logical  tconv,tcore,tkli,tspin
      integer  i,iexc,it,j,mch,mmax,nc,nv,nin,igr,irl,itmax
      integer  ir,in,ispin,nspinmx,iexcs,mup
      integer  ninu(ms),n(ms),l(ms),m(ms),nstart(ms),ncore(ms),nend(ms)
      integer  ninmx(ms)
      real*8   al,als,amesh,bl,z,sf,eeel,eexc,ex,ec,et,uld,de,dmelm
      real*8   eeig,eeig_old,feps
      real*8   e(ms),f(ms),rpk(ms),r(mx),rho(mx),rhop(mx),rhopp(mx)
      real*8   dc(mx),dcp(mx),dcpp(mx),vi(mx)
      real*8   u(mx),up(mx),upp(mx),uu(mx,ms)
      real*8   ves(mx),wa(mx),wb(mx),rnode(ms)
      real*8   fmom,sfval,eeel_val,ves_val(mx),rho_val(mx)
      real*8   d(mx,ms),dp(mx,ms),dpp(mx,ms)
      real*8   dcs(mx,ms),dcsp(mx,ms),dcspp(mx,ms)
      real*8   rhonow(mx),rhonowp(mx),rhonowpp(mx)
      real*8   vee(mx,ms),vxc(mx,ms)
      real*8   vi1(mx,ms),vn(mx,ms),vo(mx,ms),vo1(mx,ms)
      real*8   r1i(mx),r2i(mx),r3i(mx),r4i(mx)
c kli
      logical  tstart,tinpot
      integer  mp(ms),norb_kli,iukli,iexckli
      real*8   fp(ms),dens(mx),vckli(mx),vsl(mx),vxc_v(mx)
      real*8   vsl_v(mx),svkli(ms),svm(ms),svmdummy(ms)
      parameter (iukli=80)
c sic
      logical  tsic
      real*8   esic(ms),vsic(mx,ms),dorb(mx,2),dnil(mx,2)
      real*8   blcom
      real*8   vxc_spin(mx,2),vglob(mx),vsic_old(mx,ms)

c meta gga 
      logical  tmgga
      real*8   uup(mx,ms)
      external dmelm
c
      feps=1.d-9
      itmax=it
      amesh=r(2)/r(1)
      al=log(amesh)
      uu=0.d0 
      u= 0.d0 
      up=0.0d0
      upp=0.0d0

c spin initialize
      tspin=.false.
      nspinmx=1
      if(iexc.lt.0) nspinmx=2
      iexc=abs(iexc)
c kli initialize
      iexckli=iexc
      if(iexc.eq.12) iexckli=0
      do i=1,nc+nv
        mp(i) = n(i)
      enddo
      tstart=.true.
c
c sic initialize
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
        dorb(i,2)=0.d0
        dnil(i,1)=0.d0
        dnil(i,2)=0.d0
        r1i(i)=1.d0/r(i)
        r2i(i)=r1i(i)*r1i(i)
        r3i(i)=r1i(i)*r2i(i)
        r4i(i)=r2i(i)*r2i(i)
      enddo
c initialize frozen core
      if(.not. tcore) then
        do ispin=1,nspinmx
          do j=1,mmax
            dcs(j,ispin)=dcs(j,ispin)/r2i(j)
            dcsp(j,ispin)=dcsp(j,ispin)/r3i(j)
            dcspp(j,ispin)=dcspp(j,ispin)/r4i(j)
          enddo
        enddo
      endif
c
      sf=0.d0
      do ispin=1,nspinmx
        nstart(ispin)=(nc+nv)*(ispin-1)+1
        nend(ispin)=nstart(ispin)+nv+nc-1
        do i=nstart(ispin),nend(ispin)
          sf=sf+f(i)
          ninu(i)=1
        enddo
      enddo
c
c begin, bl is the mxing factor for the potential in the
c self-consistency loop
      bl=0.7d0
      blcom=1.d0-bl
      als=al*al
c
c return point for self-consistency loop
      do ispin=1,nspinmx
        nstart(ispin)=(nc+nv)*(ispin-1)+(nc+1)
        ncore(ispin)=0
        nend(ispin)=nstart(ispin)+nv-1
        if(tcore) then
           nstart(ispin)=(nc+nv)*(ispin-1)+1
           ncore(ispin)=nstart(ispin)+nc-1
           nend(ispin)=nstart(ispin)+nv+nc-1
        endif
cdebug
!d       write(ie,*)
!d    1   '%sratom - ispin,nstart,ncore,nend=',
!d    1   ispin,nstart(ispin),ncore(ispin),nend(ispin) 

      enddo
      do it=1,itmax
c
c initialize for new iteration
        tconv=.true.
        mup=0
        do ispin=1,nspinmx
          do j=1,mmax
            rho(j)=0.0d0
            d(j,ispin)=0.0d0
            dp(j,ispin)=0.0d0
            dpp(j,ispin)=0.0d0
          enddo
        enddo
c
        if(tsic) call dcpv(mx,1,mmax,vi,vglob)

c start spin loop ------------------------------------------------------
        do ispin=1,nspinmx
        ninmx(ispin)=0

        do i=nstart(ispin),nend(ispin)
          if(f(i) .gt. feps .or. tconv) then
            et=e(i)
            if(tsic) call dsuv(mx,1,mmax,vglob,vsic(1,i),vee(1,1))

c debug
!d           write(ie,*) 'it state',it,i,n(i),l(i),m(i),e(i)
            call dftseq(irl,z,mmax,r,n(i),l(i),1.d0,vee(1,ispin),wa,wb,
     1                  nin,mch,uld,e(i),u,up,upp)
c debug
!d           write(ie,*) 'it state',it,i,n(i),l(i),m(i),e(i)

            call dcpv(mx,1,nin,u,uu(1,i))
            if(e(i) .ne. et) tconv=.false.
            if(tconv) call dcpv(mx,1,nin,up,uup(1,i)) ! meta gga
            if(nin .gt. ninmx(ispin)) ninmx(ispin)=nin
            if(nin .gt. mup) mup=nin
            ninu(i)=nin
          endif

c accumulate charge & radial derivatives, rhop & rhopp are calculated as
c derivatives w.r.t. r, while up & upp are in terms of the logarithmic
c radial variable
c orbital density and update of spin densities
          do j=1,nin
            rhonow(j)=u(j)*u(j)
            d(j,ispin)=d(j,ispin)+f(i)*rhonow(j)
          enddo
          if(igr.eq.1) then
            do j=1,nin
              rhonowp(j)=2.d0*(up(j)/al-u(j))*u(j)
              rhonowpp(j)=2.d0*( (up(j)/al-u(j))*(up(j)/al
     1               -3.d0*u(j)) +(upp(j)-al*up(j))*u(j)/als )
              dp(j,ispin)=dp(j,ispin)+f(i)*rhonowp(j)
              dpp(j,ispin)=dpp(j,ispin)+f(i)*rhonowpp(j)
            enddo
          endif

c store core density
          if(i.eq.ncore(ispin)) then
            do j=1,ninmx(ispin)
              dcs(j,ispin)=d(j,ispin)
              dcsp(j,ispin)=dp(j,ispin)
              dcspp(j,ispin)=dpp(j,ispin)
              d(j,ispin)=0.d0
              dp(j,ispin)=0.d0
              dpp(j,ispin)=0.d0
            enddo
            do j=ninmx(ispin)+1,mmax
              dcs(j,ispin)=0.d0
              dcsp(j,ispin)=0.d0
              dcspp(j,ispin)=0.d0
              d(j,ispin)=0.d0
              dp(j,ispin)=0.d0
              dpp(j,ispin)=0.d0
            enddo
          endif
c
c find outermost peak of wavefunction (needed for peseudopotentials)
          if(tconv .or. it .ge. itmax) then          
            do j=nin-1,1,-1
              if(up(j)*up(j+1) .lt. 0.0d0) goto 341
            enddo
  341       rpk(i)=r(j)
          endif
c
        enddo ! end orbital loop 

c end spin loop --------------------------------------------------------
        enddo
c
c make total density rho()
        do ispin=1,nspinmx
          do j=1,mmax
            d(j,ispin)=(d(j,ispin)+dcs(j,ispin))*r2i(j)
            rho(j)=rho(j)+d(j,ispin)
            dp(j,ispin)=(dp(j,ispin)+dcsp(j,ispin))*r3i(j)
            dpp(j,ispin)=(dpp(j,ispin)+dcspp(j,ispin))*r4i(j)
co          rho(j)=rho(j)/r(j)**2 +dc(j)
co          rhop(j)=rhop(j)/r(j)**3 +dcp(j)
co          rhopp(j)=rhopp(j)/r(j)**4 +dcpp(j)
          enddo
        enddo
c
c sic potentials
        if(tsic) then
          do i=nstart(ispin),nend(ispin)
            do j=1,mmax
              dorb(j,1)=(uu(j,i)/r(j))**2
            enddo
            call vestat(mmax,1.d0,eeel,r,dorb(1,1),ves,tconv)
c only LSDA no GGA here, 3 corresponds to xc option 8 
c           call vexcos(3,ninu(i),mmax,r,dorb,dnil,dnil,
c    1        vxc_spin,eexc,ex,ec,tconv)
            call dadv(mx,1,mmax,ves,vxc_spin(1,1),vsic(1,i))
            if(it .gt. 1) then
              do j=1,mmax
                vsic(j,i)=bl*vsic(j,i)+blcom*vsic_old(j,i)
              enddo
            endif
            call dcpv(mx,1,mmax,vsic(1,i),vsic_old(1,i))
            esic(i)=-(0.5*eeel+ex+ec)
          enddo
        endif
c
c global effective potential
        call vestat(mmax,sf-z,eeel,r,rho,ves,.false.)
        if(.not. tinpot .and. tkli) then
          if(it < 5) then
co          call vexcor(8,mmax,r,rho,rhop,rhopp,vxc,eexc,ex,ec,.false.)
            call vexcor(8,mmax,r,d(1,1),dp(1,1),dpp(1,1),
     1                  vxc(1,1),eexc,ex,ec,.false.)
          else
            
            do in=1,nc+nv
              fp(in)=0.5d0*f(in)
            enddo
            do ir=1,mmax
              dens(ir)=0.5d0*rho(ir)
            enddo

            if(it == 5) then
             norb_kli=0
             do i=1,nc+nv
              if(fp(i) .gt. feps) norb_kli=norb_kli+1     
             enddo
            endif

            call vklix(norb_kli,n,l,mp,fp,mmax,r,uu,
     1                 dens,vxc(1,1),vsl,svkli,svm,1,tstart,tconv,tconv)
            if(iexckli .ne. 0) then
co            call vexcor(iexckli,mmax,r,rho,rhop,rhopp,vckli,eexc,
              call vexcor(iexckli,mmax,r,d(1,1),dp(1,1),dpp(1,1),
     1                    vckli,eexc,ex,ec,tconv)
!d             write(ie,*) '%sratom_n: ', 
!d    1        it, 'adding vckli',iexckli,vckli(1),vckli(100)
              do ir=1,mmax
                vxc(ir,1)=vxc(ir,1)+vckli(ir)
              enddo
            endif
            tstart=.false.

            if(tconv) then
              write(iukli,'(a)') ' --- kli ae shifts (Ha) --- '
              do i=1,norb_kli
               write(iukli,'(a2,3i3,1x,e14.8)') '%',i,n(i),l(i),svkli(i)
              enddo
              write(iukli,*)
               
              open(11,file='vkli_ae.dat')
              write(11,'(a)') '#ae r(i) vxc(i) vslater(i)'
              do i=1,mmax
                write(11,*) r(i),vxc(i,1),vsl(i)
              enddo
              close(11)
c repeat klix but only for valence orbitals
              do ir=1,mmax
                dens(ir)=0.5d0*(rho(ir)-dc(ir))
              enddo
              i=norb_kli-nc
              call vklix(i,n(nc+1),l(nc+1),mp(nc+1),fp(nc+1),mmax,r
     1,          uu(1,nc+1),dens,vxc_v,vsl_v,svkli(nc+1),svmdummy
     1,          2,.false.,.false.,.false.)

              write(iukli,'(a)') ' --- kli valence shifts (Ha) --- '
              do i=nc+1,norb_kli
               write(iukli,'(a2,3i3,1x,e14.8)') '%',i,n(i),l(i),svkli(i)
              enddo
              write(iukli,*) 

              open(11,file='vkli_valence_ae.dat')
              write(11,'(a)') '#aevalence r(i) vxc(i) vslater(i)'
              do i=1,mmax
                write(11,*) r(i),vxc_v(i),vsl_v(i)
              enddo
              close(11)
            endif
          endif
        else
c lda or gga xc potential
          if(nspinmx.eq.1) then
c spin unpolarized
            call vexcor(iexc,mmax,r,d(1,1),dp(1,1),dpp(1,1),vxc(1,1),
     1                  eexc,ex,ec,.false.)
          else
c spin polarized
           call vexcos(iexc,mmax,mmax,r,d,dp,dpp,vxc,eexc,ex,ec,.false.)
          endif

co        call vexcor(iexc,mmax,r,rho,rhop,rhopp,vxc,eexc,ex,ec,.false.)
        endif

c
c generate next iteration`s potential using anderson method
        if(tsic) call dcpv(mx,1,mmax,vglob,vee(1,1))
        do ispin=1,nspinmx
          call dadv(mx,1,mmax,ves,vxc(1,ispin),vo(1,ispin))
          call anderson(bl,it,mmax,r,vo1(1,ispin),vo(1,ispin),
     1                  vi1(1,ispin),vee(1,ispin))
        enddo
c
c if iexc=0, overwrite with coulomb potential
        if(iexc.eq.0) then
          do ispin=1,nspinmx
            do j=1,mmax
              vee(j,ispin)=-z*r1i(j)
            enddo
          enddo
        endif

co      call dadv(mx,1,mmax,ves,vxc,vo)
co      if(tsic) call dcpv(mx,1,mmax,vglob,vi)
co      call anderson(bl,it,mmax,r,vo1,vo,vi1,vi)
c
c total energy
        eeig=0.d0
        do ispin=1,nspinmx
          do i=nstart(ispin),nend(ispin)
            eeig=eeig+f(i)*e(i)
          enddo
        enddo
        de=abs((eeig-eeig_old)/eeig)
        eeig_old=eeig

      if(tconv) then 
        print *, '*** sratom: all eigens converged ***'
        goto 24
      endif 
      enddo
  24  continue

c proper core density for output
      do ispin=1,nspinmx
        do j=1,mmax
          dcs(j,ispin)=dcs(j,ispin)*r2i(j)
          dcsp(j,ispin)=dcsp(j,ispin)*r3i(j)
          dcspp(j,ispin)=dcspp(j,ispin)*r4i(j)
        enddo
      enddo

c     open(11,file='tmp.coul')
c     write(11,'(i5,1x,e20.14)') mmax,amesh
c     do i=1,mmax
c       write(11,'(i5,1x,e20.14,1x,e20.14,1x,e20.14)') 
c    &     i,r(i),uu(i,1),-z/r(i)
c     enddo
c     close(11)

c     rho_val( : ) = exp(-2*r( : )**2)
c     sfval = fmom(0,mmax,al,1.d0,r,rho_val)
c     print *, 'gaussian total charge',sfval
c     rho_val( : ) = rho_val( : )/sfval
c     sfval = fmom(0,mmax,al,1.d0,r,rho_val)

c     call vestat(mmax,sfval,eeel_val,r,rho_val,ves_val,tconv)
c      
c     open(11,file='tmp.cancel')
c     do i=1,mmax
c     write(11,'(e12.6,1(1x,e12.6))')r(i),
c    1 (uu(i,nc+nv-1))**2
c     enddo
c     close(11)

      return
      end
c
