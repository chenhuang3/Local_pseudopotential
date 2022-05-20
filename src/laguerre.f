c $Header:$
c-----------------------------------------------------------------------
c associated laguerre polynomials
c see: G. Arfken, Mathematical methos for physicists, Academic Press,
c      San Diego, 1985, pg. 725
c
c plgr0 ..... L_n^k(x)
c plgr1 ..... 1st derivative wrt x
c plgr2 ..... 2nd derivative wrt x
c
c bas0 ...... related orthogonal polynomial c.f. basis.h
c bas1 ...
c bas2 ...
c
c bas0 etc. need to include the parameter file basis.h .
c
c note: the laguerre functions are programmed in a straightforward
c       (likely: dull) manner, taking no advantage whatsoever of 
c       recursion formulae, not worrying at all about error
c       cancellation issues, etc. --- be aware of that.
c
c Martin Fuchs, FHI, 09-1994
c
c-----------------------------------------------------------------------
      real*8 function plgr0(n,k,x)
c     program plgr

      implicit real*8 (a-h,o-z)
      parameter(nkmx=34)
      dimension fac(0:nkmx)
      logical ifc
      save    ifc,fac
      data    ifc/.true./
c
c     read(5,*) n,k,x
      if(n+k .gt. nkmx) stop 'plgr0'
      if(ifc) then
        ifc=.false.
        fac(0)=1.d0
        do i=1,nkmx
          fac(i)=i*fac(i-1)
        enddo
      endif
c
      plgr0=0.d0
      phase=-1.d0
      y=1.d0
      do m=0,n
        y=y*x
        if(m.eq.0) y=1.d0
        phase=-phase
        plgr0=plgr0
     +        + phase
     +        * fac(n+k)/(fac(n-m)*fac(k+m))
     +        / fac(m)
     +        * y
      enddo
      plgr0=plgr0/sqrt(fac(n+k)/fac(n))
c     print *, 'laguerre: ',n,k,x,plgr0
      return
c
      end
      
      real*8 function plgr1(n,k,x)
      implicit real*8 (a-h,o-z)
      parameter(nkmx=50)
      dimension fac(0:nkmx)
      logical ifc
      save    ifc,fac
      data    ifc/.true./
c
      if(n+k .gt. nkmx) stop 'plgr1'
      if(ifc) then
        ifc=.false.
        fac(0)=1.d0
        do i=1,nkmx
          fac(i)=i*fac(i-1)
        enddo
      endif
c
      if(n .eq. 0) then
        plgr1=0.d0
      else
        plgr1=0.d0
        phase=1.d0
        y=1.d0
        do m=1,n
          y=y*x
          if(m.eq.1) y=1.d0
          phase=-phase
          plgr1=plgr1
     +          + phase
     +          * fac(n+k)/(fac(n-m)*fac(k+m))
     +          / fac(m-1)
     +          * y
        enddo
        plgr1=plgr1/sqrt(fac(n+k)/fac(n))
      endif
c     print *, 'laguerre: ',n,k,x,plgr1
      return
c
      end
      
      
      real*8 function plgr2(n,k,x)
      implicit real*8 (a-h,o-z)
      parameter(nkmx=50)
      dimension fac(0:nkmx)
      logical ifc
      save    ifc,fac
      data    ifc/.true./
c
      if(n+k .gt. nkmx) stop 'plgr1'
      if(ifc) then
        ifc=.false.
        fac(0)=1.d0
        do i=1,nkmx
          fac(i)=i*fac(i-1)
        enddo
      endif
c
      if(n .le. 1) then
        plgr2=0.d0
      else
        plgr2=0.d0
        phase=-1.d0
        y=1.d0
        do m=2,n
          y=y*x
          if(m.eq.2) y=1.d0
          phase=-phase
          plgr2=plgr2
     +          + phase
     +          * fac(n+k)/(fac(n-m)*fac(k+m))
     +          / fac(m-2)
     +          * y
        enddo
        plgr2=plgr2/sqrt(fac(n+k)/fac(n))
      endif
c     print *, 'laguerre: ',n,k,x,plgr2
      return
c
      end

      real*8 function bas0(n,x)
      implicit real*8 (a-h,o-z)
      include 'basis.h'
      external plgr0
      xx=x*a
      bas0=b*exp(-.5d0*xx)*x*plgr0(n,k,xx)
      if(l .gt. 0) bas0=x**l*bas0
      return
      end
      
      real*8 function bas1(n,x)
      implicit real*8 (a-h,o-z)
      include 'basis.h'
      external bas0,plgr1
      xx=x*a
      bas1=(-c+ll/x)*bas0(n,x)+a*b*exp(-.5d0*xx)*x*plgr1(n,k,xx)
      return
      end
      
      real*8 function bas2(n,x)
      implicit real*8 (a-h,o-z)
      include 'basis.h'
      external bas0,bas1,plgr2
      xx=x*a
      bas2=-ll/x**2*bas0(n,x)+(-c+ll/x)*bas1(n,x)
     +     +a*b*x*exp(-.5d0*xx)
     +     *((-c+ll/x)*plgr1(n,k,xx)+a*plgr2(n,k,xx))
      return
      end

      real*8 function geta(i)
      integer i
      include 'basis.h'
      geta = a
      return
      end
c
