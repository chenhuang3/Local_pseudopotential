c $Header:$
c***********************************************************************
c generalized norm-conserving pseudopotentials of Hamann type
c original version by D.R. Hamann, gncpp
c cite:   D.R. Hamann, Phys. Rev. B 40, 2980 (1989).
c version adopted at Fritz-Haber-Institut Berlin
c Martin Fuchs 12-1995
c
c input
c io ........... I/O unit
c z ............ atomic number
c epspp ........ accuracy parameter for matching reference and pseudo
c                energies
c l ............ angular momentum channel
c rc ........... cutoff radius
c ein .......... reference energy
c mode ......... [modus for radial schroedinger eq]
c mch .......... grid index for reference radius
c mmax ......... maximum grid index
c r() .......... radial grid
c u() .......... radial all-electron wavefunction, changed upon output
c up() ......... dto. 1st derivative (on the transformed grid!)
c upp() ........ dto. 2nd derivative
c vin() ........ all-electron potential
c
c output
c u() .......... radial pseudowavefunction
c up() ......... dto. 1st derivative (on the transformed grid!)
c upp() ........ dto. 2nd derivative
c v() .......... pseudopotential
c 
c parameter
c alam ......... healing parameter for wavefunction 
c 
c***********************************************************************
      subroutine hamann(io,z,epspp,l,rc,ein,mode,mch
     1,  mmax,r,u,up,upp,vin,v)
c
      implicit none
      include  'parameter.h'

      integer  io,mmax,mode,mch,nrc,l,iter,nin,i,mchi
      real*8   z,alam,epspp,ein,rc,al,amesh,etest,umch,upmch
      real*8   uld,cl,texp,ro,sv,sf,sx,dcl,gam,del
      real*8   r(mx),vin(mx),v(mx),fc(mx),vc(mx)
      real*8   wm,wa(mx),wb(mx),u(mx),up(mx),upp(mx)

      parameter (alam=3.5d0)

c initialization
      wm=1.d0
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
      enddo
      amesh=r(2)/r(1)
      al=log(amesh)

      etest=ein
      nrc=log(rc/r(1))/al+1
      upmch=up(mch)
      umch=u(mch)
      uld=upmch/umch

!d     write(ie,*) '& hamann -- alam            ',alam
!d     write(ie,*) '& hamann -- l,e,mode        ',l,etest,mode
!d     write(ie,*) '& hamann -- mch,umch,upmch  ',mch,umch,upmch
!d     write(ie,*) '& hamann -- rc,nrc,vin(nrc) ',rc,nrc,vin(nrc)

c construct pseudopotential 
c form cutoff potential
      do i=1,mmax
        texp=(r(i)/dble(rc))**alam
        if(texp .lt. 700.0d0) then
          fc(i)=exp(-texp)
        else
          fc(i)=0.0d0
        end if
        vc(i)=(1.0d0 - fc(i))*vin(i)
      enddo

c iteration to find correct additive term
      iter=0
      cl = vin(nrc)
c
   30 iter=iter+1
        if(iter .gt. 30) then
          write(ie,*) 
     1      '& hamann - stop: convergence error l,rc,rmch,ein',
     1      l,rc,r(mch),ein
          stop 
        endif
c
        do i=1,mmax
          v(i)=vc(i) + cl*fc(i)
        enddo
c
        etest=ein
        mchi=mch
        call dftseq(mode,z,mmax,r,l+1,l,wm,v,wa,wb,
     1    nin,mchi,uld,etest,u,up,upp)

!d       write(ie,*) '& hamann --- iter,etest   ',iter,etest
!d       write(ie,*) '& hamann --- nin,mchi,uld ',nin,mchi,uld

c
c perform integral times cutoff function
c also integrals for norm correction
        ro=r(1)/sqrt(amesh)
        sv=ro**(2*l+3)*(u(1)**2*fc(1)/r(1)**(2*l+2))/dble(2*l+3)
        sf=ro**(2*l+3)*fc(1)**2/dble(2*l+3)
        sx=ro**(2*l+3)*(u(1)*fc(1)/r(1)**(l+1))/dble(2*l+3)
c
        do i=1,nin-3
          sf=sf + al*(r(i)**(2*l+3))*fc(i)**2
          sx=sx + al*(r(i)**(l+2))*fc(i)*u(i)
          sv=sv + al*r(i)*u(i)**2*fc(i)
        enddo
c
        sv=sv + al*(23.d0*r(nin-2)*u(nin-2)**2*fc(nin-2)
     1            + 28.d0*r(nin-1)*u(nin-1)**2*fc(nin-1)
     1          +  9.d0*r(nin  )*u(nin  )**2*fc(nin  ))/24.d0
        sf=sf + al*(23.d0*(r(nin-2)**(2*l+3))*fc(nin-2)
     1            + 28.d0*(r(nin-1)**(2*l+3))*fc(nin-1)
     1            +  9.d0*(r(nin  )**(2*l+3))*fc(nin  ))/24.d0
        sx=sx + al*(23.d0*(r(nin-2)**(l+2))*u(nin-2)*fc(nin-2)
     1            + 28.d0*(r(nin-1)**(l+2))*u(nin-1)*fc(nin-1)
     1            +  9.d0*(r(nin  )**(l+2))*u(nin  )*fc(nin  ))/24.d0
c
c perturbation theory to improve potential
        dcl=(ein - etest)/sv
        cl=cl + dcl
!d       write(ie,*) '& hamann --- dcl,sv       ',dcl,sv
c
      if(abs(dcl) .gt. epspp) goto 30
c
c end of iteration label 30
c
c find coefficient for norm correction
      gam=abs(umch/u(mch))
      sx=sx*gam**2
      sf=sf*gam**2
      del=(-sx + (sx/abs(sx))
     &    *sqrt(sx*sx - sf*(gam**2-1.0d0)))/sf
c
c construct final wavefunction
      do i=1,nin
        u(i)=gam*(u(i) + del*r(i)**(l+1)*fc(i))
      enddo
c
c construct final pseudopotential
      do i=1,nin
        v(i)=v(i) + (gam*del*r(i)**(l+1)*fc(i)/(2.d0*u(i)))
     1           *((alam**2*(r(i)/rc)**(2*alam)
     1             -(2.d0*alam*l + alam*(alam+1.d0))
     1               *(r(i)/rc)**alam)/r(i)**2
     1            + 2.d0*etest - 2.d0*v(i))
      enddo
      do i=nin+1,mmax
        v(i)=vin(i)
      enddo

      call dftseq(mode,z,mmax,r,l+1,l,wm,v,wa,wb,
     1   nin,mchi,uld,etest,u,up,upp)
      
!d     write(ie,*) '& hamann --- vin(nrc),vin(mch) ',vin(nrc),vin(mch)
!d     write(ie,*) '& hamann --- v(nrc),v(mch)     ',v(nrc),v(mch)
!d     write(ie,*) '& hamann --- gam,uld           ',gam,up(mch)/u(mch)
!d     write(ie,*) '& hamann --- ein,etest         ',ein,etest
!d     write(ie,*) '& hamann --- returning from l= ',l

c
      return
      end
c
