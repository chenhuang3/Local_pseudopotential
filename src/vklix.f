c $Header: /f21/a/fuchs/Psp/Dkef/Ds/RCS/vklix.f,v 1.1 95/08/22 17:42:15 fuchs Exp Locker: fuchs $
c
c compute Krieger-Li-Iafrate (kli) approximation to the optimized effective
c potential of exchange-interaction, assuming spherical symmetry
c
c input
c iorb ........... number of orbitals in present spin channel
c n() ............ principal quantum number
c l() ............ angular momentum quantum number
c m() ............ spin projection 
c f() ............ orbital occupation numbers
c mmx ............ upper limit index on radial grid
c r() ............ radial mesh
c ud0() .......... radial kohn sham orbitals
c                  normalized to 4*pi
c                  ordered wrt to increasing eigenvalues
c dens() ......... spin density in present channel
c mode ........... 1 ........ evaluate orbital shifts
c                  2 ........ use input orbital shifts
c tene ........... .true. ... evaluate diagonal energy contributions
c                             i.e. exact exchange energy
c tstart ......... .true. ... initialize 
c
c output
c voep() ......... kli approximated oep
c sv()   ......... orbital shifts
c svx()  ......... orbital exchange energy
c
c Martin Fuchs 20-08-1995
c 
      subroutine vklix(iorb,n,l,m,f,mmx,r,ud0,dens,voep,vsl,sv,svx
     1,mode,tstart,tene,tprint)

      implicit real*8 (a-h,o-z)
      include 'parameter.h'

      character*1 sni,sli,smi
      character*1 snj,slj,smj
      logical tene,tprint,tstart
      integer iukli,mode
      real *8 al,amesh,ahat

c     common/arr/r(mx),vi(mx),u(mx),up(mxr,),upp(mx)
c     common/par/z,amesh,al,mmax,nin,mch,iexc
c     common/vslater/vsl(mx)

      parameter(iukli=80,iuslater=81)
      parameter(zero=0.d0,one=1.d0,eps=1.d-12)

      dimension ud0(mx,iorb),dens(mx),voep(mx),dorb(mx,ms),vorb(mx,ms)
     +,   vxi(mx,ms),indx(ms),sm(ms,ms),sv(ms),svx(ms),vav(mx)
     +,   n(ms),l(ms),m(ms),f(ms),svsave(ms),smsave(ms,ms)
     +,   uc(mx),wc(mx),rsqi(mx),r(mx),vsl(mx)
      data icount/0/
      save icount,rsqi

c functions
      ahat(la)=2.d0*la+1

c files
      if(tprint) then
        open(unit=iuslater,file='vslater_accumul.dat',status='unknown')
        write(iukli,'(a,/)') 's List of Fock matrix elements ---'
      endif

c initializations and limiting radial range (overflow precaution)
      amesh=r(2)/r(1)
      al=log(amesh)
      if(tstart) icount = 0
      icount=icount+1
      mup=mmx
      if(icount .eq. 1) then
        write(iukli,*) '%vklix: check orbital parameters'
        write(iukli,*) 'iorb=',iorb
        do i=1,iorb
          write(iukli,*) 'n l f', n(i),l(i),f(i)
        enddo
        write(iukli,*) '%vklix: end printout'
        do j=1,mmx
          rsqi(j)=one/r(j)**2
        enddo
      endif
      do j=1,mmx
        voep(j)=zero
        vav(j)=zero
        if(dens(j) .lt. eps) mup=min(mup,j)
      enddo
      do j=1,iorb+1
        do i=1,iorb+1
          sm(j,i)=zero
      enddo
      enddo

c compute (spherically averaged) orbital densities
      do i=1,iorb
        do j=1,mmx
          dorb(j,i)=ud0(j,i)*ud0(j,i)*rsqi(j)
          vorb(j,i)=zero
        enddo
      enddo

c compute (spherical) exchange potentials
      do i=1,iorb
        if(tprint) then
          write(iukli,'(/,a)') 's Fock matrix elements FME:'
        endif
        do j=1,iorb
          if(i .ge. j .and. tprint) then
            call labels(n(i),l(i),m(i),sni,sli,smi)
            call labels(n(j),l(j),m(j),snj,slj,smj)
            write(iukli,'(a,4x,2i2,30x,a1,a1,1x,a1,a1)')  
     +       's ij',i,j,sni,sli,snj,slj 
            write(iukli,'(a,1x,f8.4,2x,f8.4)')
     +       's f_i f_j    ',f(i),f(j)
          endif
          do ll=abs(l(i)-l(j)),l(i)+l(j)

            cj=-f(j)/ahat(ll)*acgc(l(j),l(i),ll)**2
c configuration averaging for open shells
            if(j .eq. i .and. f(j) .lt. ahat(l(j))) then
c monopole term
              if(ll .eq. 0) then
                cj=-one
c multipole terms
              else
                cj=-max(zero,(f(j)-one))
     +             *ahat(l(j))/(2*l(j)*ahat(ll))*acgc(l(j),l(j),ll)**2
              endif
            endif

cmf
c           cj=-f(j)*ack(l(i),l(j),ll)


            if(cj .lt. zero) then

!d              if(icount .lt. 2) then
!d                  print *,
!d    + i,j,'  ',l(i),l(j),ll,l(i)+l(j)+ll,cj,acgc(l(j),l(i),ll)
!d                endif

              call arhf
     +             (mup,r,l(i),l(j),ll,ud0(1,i),ud0(1,j),uc,wc,vxi(1,j))

              do k=1,mup
                vorb(k,i)=vorb(k,i)+cj*(vxi(k,j)*ud0(k,j)*rsqi(k))
              enddo

            endif

            if(i .ge. j .and. tprint) then
             do k=1,mup
              uc(k)=vxi(k,j)*ud0(k,j)*rsqi(k)
             enddo
             dfock=dmelm(mup,al,r,ud0(1,i),uc)
             write(iukli,'(a,28x,i2,i2,i2)')
     +         's l_i l_j L  ',l(i),l(j),ll
             write(iukli,'(a,f8.4,2x,f8.4)')
     +         's c_j c_j/f_j ',cj,cj/f(j)
             write(iukli,'(a,2(2x,e12.6))')    
     +         's FME FME*c/f ',dfock,dfock*cj/f(j)
             write(iukli,'(a)') 's'
            endif

          enddo

        enddo

c compute density weighted (Slater-like) potential
        do j=1,mup
          vav(j)=vav(j)+f(i)*ud0(j,i)*vorb(j,i)
        enddo
        if(tprint) then
         do j=1,mup
          vsl(j)=vav(j)/dens(j)
          write(iuslater,*) r(j),vsl(j) 
         enddo
         write(iuslater,*) '&',j
        endif

      enddo

      do j=1,mup
        vav(j)=vav(j)/dens(j)
        vsl(j)=vav(j)
      enddo

c average potential option
c     do j=1,mup
c       voep(j)=vav(j)
c     enddo
c     return

      if(mode == 1) then
c potential constants
      do i=1,iorb
        do j=1,mup
          uc(j)=vav(j)*ud0(j,i)*rsqi(j)-vorb(j,i)
        enddo
        sv(i)=dmelm(mup,al,r,ud0(1,i),uc)
!d       svsave(i)=sv(i)
      enddo

c linear equation setup 
      sv(iorb)=zero
      do i=1,iorb-1
        do j=1,mup
          uc(j)=dorb(j,i)/dens(j)
        enddo
        do j=1,i
          sm(i,j)=-f(j)*dmelm(mmx,al,r,uc,dorb(1,j))
          sm(j,i)=sm(i,j)*f(i)/f(j)
          if(i .eq. j) sm(i,j)=sm(i,j)+1.d0
!d         smsave(i,j)=sm(i,j)
        enddo
      enddo

!d     print *, '--- vkli corrector matrix ---'
!d     do i=1,iorb-1
!d       print *, i, (sm(i,j),j=1,iorb-1)
!d     enddo

c solve system --- sv(i) are the kli correctors
      idm=iorb-1
      call ludcmp(sm,idm,ms,indx,xxx)
      call lubksb(sm,idm,ms,indx,sv)

!d     print *, '--- vkli average difference potential ---'
!d     print *, 'iorb',iorb
!d     do i=1,iorb
!d       print *, i, n(i),l(i),m(i),f(i),sv(i)
!d     enddo

!d     print *, '--- vkli verify solution ---'
!d     do i=1,iorb-1
!d       test=zero
!d       do j=1,iorb-1
!d         test=test+smsave(i,j)*sv(j)
!d       enddo
!d       print *, i, test,svsave(i)
!d     enddo
      else if (mode .ne. 2) then
        print *, '& - vklix - warning: unimplemented mode',mode
      endif

!d     print *, '--- vkli orbital correctors ---'
!d     do i=1,iorb-1
!d       print *, i, sv(i)
!d     enddo

c compute kli potential
      do i=1,iorb

!d       iunit=50+i
!d       if(m(i).eq.1) then
!d         open(iunit)
!d         rewind(iunit)
!d         do j=1,mup
!d           write(iunit,*) r(j),f(i)*ud0(j,i)*vorb(j,i)/dens(j)
!d         enddo
!d         close(iunit)
!d       endif

        do j=1,mup
          voep(j)=voep(j)+f(i)*sv(i)*dorb(j,i)
        enddo

      enddo
      do j=1,mup
        voep(j)=vav(j)+voep(j)/dens(j)
      enddo
      do j=mup+1,mmx
        voep(j)=-one/r(j)
      enddo
     
c orbital exchange energies
      if(tene) then
        do i=1,iorb
          svx(i)=dmelm(mup,al,r,ud0(1,i),vorb(1,i))
        enddo
      endif

!d     print *, '--- vkli done ---'
c
      if(tprint) then
        close(iuslater)
        write(iukli,'(a,/)') 's ---------------------------'
      endif
      return
      end
c
c
c ancillary stuff
c
c
      real*8 function ack(l1,l2,ll)
      implicit none
     
      integer l1,l2,ll
      integer arg,factorial,i,num1,den1,parg
      real*8  y1,y2
      external factorial

      parg=l1+l2+ll
      if(mod(parg,2) .eq. 0) then
       arg=-l1+l2+ll
       num1=factorial(arg)
       arg=-l2+l1+ll
       num1=num1*factorial(arg)
       arg=-ll+l1+l2
       num1=num1*factorial(arg)
       arg=l1+l2+ll+1
       den1=factorial(arg)
       ack=dble(num1)/dble(den1)

       parg=parg/2
       num1=factorial(parg)**2
       ack=ack*dble(num1)

       arg=parg-l1
       den1=factorial(arg)**2
       ack=ack/dble(den1)

       arg=parg-l2
       den1=factorial(arg)**2
       ack=ack/dble(den1)

       arg=parg-ll
       den1=factorial(arg)**2
       ack=ack/dble(den1)
      else
       ack=0.d0
      endif
c     print *,'ack',parg,mod(parg,2),ack
      
      return
      end

      integer function factorial(i)
      integer   i
      integer   ftmp,j

      ftmp=1
      do j=2,i
        ftmp=ftmp*j
      enddo
      factorial=ftmp
      return
      end

      
       

      
 
      
