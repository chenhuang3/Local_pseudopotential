c $Header:$
c
c logarithmic derivative for separable (Kleinman-Bylander) potential
c at radius r(ich), Hartree a.u.
c
c input
c l ............ angular momentum 
c e ............ energy 
c ich .......... grid index for evaluation
c mmax ......... grid used up to mmax
c r() .......... radial grid
c vkb() ........ u_l(r) * ( v_l(r) - v_loc(r) )
c v() .......... screened local potential
c 
c output
c dlkb ......... r * logarithmic derivative
c dedl ......... energy derivative of dlkb        <disabled>
c nodes ........ number of nodes in wfct.         <disabled>
c
c original version by Gonze & Stumpf
c
      subroutine derlkb(l,e,ich,dlkb,dkb,dedl,nodes,mmax,r,vkb,v)
c
      implicit real*8 (a-h,o-z)
      
      include  'parameter.h'

      dimension cf(mx),cfkb(mx)
      dimension wkb(mx),v1kb(mx),vkb(mmax)
      dimension r(mmax),u(mx),up(mx),upp(mx),v(mmax)

      sls=float(l*(l+1))
      al=0.1d0*log(r(11)/r(1))
      als=al*al
c     coefficient array for u in differential eq.
      do 20 i=1,mmax
   20 cf(i)=als*sls + 2.0*als*(v(i)-e)*r(i)**2

c start wavefunction with series
      ups0=(v(2)*r(2)-v(1)*r(1))/(r(2)-r(1))
      zcorr=-(v(1)-ups0)*r(1)
      c2=-zcorr/float(l+1)
      c3=zcorr**2/((2*l+3)*(l+1))+(ups0-e)/(2*l+3)
c
      do 50 i=1,4
      u(i)=r(i)**(l+1) + c2*r(i)**(l+2) + c3*r(i)**(l+3)
      up(i)=al*((l+1)*r(i)**(l+1) + c2*(l+2)*r(i)**(l+2)
     & + c3*(l+3)*r(i)**(l+3))
   50 upp(i)=al*up(i)+cf(i)*u(i)
c
c     outward integration using predictor once, corrector twice
c 
cmf predictor corrector auxiliaries renamed
      do 70 i=4,ich-1
      u(i+1)=u(i)+aeo(up,i)
      up(i+1)=up(i)+aeo(upp,i)
      do 60 it=1,2
      upp(i+1)=al*up(i+1)+cf(i+1)*u(i+1)
      up(i+1)=up(i)+aio(upp,i)
   60 u(i+1)=u(i)+aio(up,i)
   70 continue
cmf
c
c     uout=u(ich)
c     upout=up(ich)
c     write(ie,998)uout,upout
c998  format('derlkb - uout ',d16.6,'     upout=',d16.6)
c     if(ich.ne.1)goto 200
c
c  normalisation of this first function
c
      do 75 i=1,ich
 75   wkb(i)=u(i)*vkb(i)*r(i)   
      anorm=dsimpson(ich,wkb,al)
c     write(ie,999)anorm,dkb
c999  format('derlkb - norm ',d16.6,'     dkb=',d16.6)
      do 80 i=1,ich
 80   wkb(i)=u(i)*dkb/anorm
      wkbp=up(ich)*dkb/anorm
c
c the second function
c
c     coefficients of the inhomogeneous term
c
      a0=vkb(1)/r(1)**(l+1)  
      do 100 i=1,ich
 100  cfkb(i)=2*al**2*r(i)**2*dkb/anorm/a0*vkb(i)   
c start wavefunction with series
      c3=(ups0-e+zcorr**2/(l+1)+1)/(2*l+3)
c
      do 150 i=1,4
      u(i)=(r(i)**(l+1) + c2*r(i)**(l+2) + c3*r(i)**(l+3))*dkb/anorm
      up(i)=al*((l+1)*r(i)**(l+1) + c2*(l+2)*r(i)**(l+2)
     & + c3*(l+3)*r(i)**(l+3))*dkb/anorm
  150 upp(i)=al*up(i)+cf(i)*u(i)+cfkb(i)
c
c     outward integration using predictor once, corrector twice
c 
      do 170 i=4,ich-1
        u(i+1)=u(i)+aeo(up,i)
        up(i+1)=up(i)+aeo(upp,i)
        do 160 it=1,2
          upp(i+1)=al*up(i+1)+cf(i+1)*u(i+1)+cfkb(i+1)
          up(i+1)=up(i)+aio(upp,i)
 160      u(i+1)=u(i)+aio(up,i)
 170  continue
c     
c     subtraction of function w to give Vv1=0
c
      do 175 i=1,ich
 175  v1kb(i)=u(i)*vkb(i)*r(i)   
      anorm1=dsimpson(ich,v1kb,al)
      do 178 i=ich,ich
 178  v1kb(i)=u(i)-wkb(i)/dkb*anorm1
      do 180 i=ich,ich
 180  v1kb(i)=v1kb(i)*a0/dkb*anorm
      v1kbp=(up(ich)-wkbp/dkb*anorm1)*a0/dkb*anorm
c
c     solution
c
      do 190 i=ich,ich
 190  u(i)=v1kb(i)+wkb(i)
      uout=u(ich)
      upout=v1kbp+wkbp

 200  dlkb=upout/al/uout-1.0d0

c     write(ie,*) '& derlkb:'
c     do 205 i=1,ich
c     write(ie,9996)i,r(i),u(i),v1kb(i),wkb(i)
c9996 format(i3,4f15.7)
c205  continue
c     if(ich.ne.1)stop
c
c   energy derivative of the logarithmic derivative
c
c     do 210 i=1,ich
c210  up(i)=u(i)**2*r(i)
c     dedl=dsimpson(ich,up,al)
c     dedl=-2*dedl*r(ich)/u(ich)**2   
c
c   number of nodes
c
c     nodes=0
c     do 220 i=1,ich-1
c220  if(u(i)*u(i+1).lt.0.0)nodes=nodes+1
c     write(ie,997)dlkb,dedl
c997  format('dlkb =',d16.6,'     dedl =',d16.6)
c
      return
      end

      real*8 function dsimpson(n,f,h)
      implicit real*8 (a-h,o-z)
      include 'parameter.h'
c
c     simpson integration for f with step length h for n points
c     n should be odd !
c
      dimension f(n)
c
      nm12=(n-1)/2
      if(nm12*2.ne.n-1) then 
        write(ie,*) '& simpson - stop: wrong n',n
        stop
      endif
c
      dsimpson=4.*dssum(nm12,f(2),2)+2.*dssum(nm12-1,f(3),2)+f(1)+f(n)
      dsimpson=dsimpson*h/3.
c
      return
      end

      REAL*8 FUNCTION dSSUM ( N, DX, INCX )
      implicit REAL*8 (a-h,o-z)
C
C     SUM OF THE ELEMENTS OF A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78. CJB
C
      REAL*8 DX(*)
      INTEGER I,INCX,N,NINCX,LINCX
C
      dSSUM=0.d0
      IF(N.LE.0)RETURN
C
      NINCX = N*INCX
      LINCX = INCX + 1
      SSUM = DX(1)
      DO 10, I = LINCX,NINCX,INCX
        SSUM = SSUM + DX(I)
   10 CONTINUE
      dSSUM=SSUM
      RETURN
      END


