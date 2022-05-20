c $Header:$
c testing spherical harmonic sum rules
c here: integral_solidangle complexconj[Y_lama] * Y_lbmb * Y_lcmc
c
      real*8 function gaunt(la,ma,lb,mb,lc,mc)

      implicit real*8 (a-h,o-z)
      complex*16 ya,yb,yc,cylm

      parameter (mx=12,zero=0.d0)
      dimension x(mx),y(mx),w(mx),p(mx),yinv(mx),f(mx),
     &          ya(mx),yb(mx),yc(mx)    
      common /testfunction/l,m
      data x,w/
     & -0.981560634246719244,-0.904117256370474909,
     & -0.769902674194304693,-0.587317954286617483,
     & -0.367831498998180184, -0.125233408511468913,
     & 0.125233408511468913, 0.367831498998180184,
     & 0.587317954286617483, 0.769902674194304693,
     & 0.904117256370474909, 0.981560634246719244,
c
     & 0.471753363865118694d-01, 0.106939325995318343,
     & 0.160078328543346360, 0.203167426723065925,
     & 0.233492536538354778, 0.249147045813402884,
     & 0.249147045813402884, 0.233492536538354778,
     & 0.203167426723065925, 0.160078328543346360,
     & 0.106939325995318343, 0.471753363865118694d-01/

      pi=4*atan(1.d0)
10    imx=mx

c
c gauss integration
      do i=1,1

        if(notri(la,lb,lc) .ge. 0) then
          if(-ma+mb+mc .eq. 0) then
            phi=h*(i-1)
            do j=1,imx
              ya(j)=conjg(cylm(la,ma,x(j),phi))
              yb(j)=cylm(lb,mb,x(j),phi)
              yc(j)=cylm(lc,mc,x(j),phi)
            enddo
            p(i)=zero
            do j=1,imx
              p(i)=p(i)+w(j)*ya(j)*yb(j)*yc(j)
            enddo
            gaunt=2*pi*p(i)
          else
c           print *, '-ma+mb+mc != 0'
            gaunt=zero
          endif
        else
c         print *, la,lb,lc
          gaunt=zero
        endif

      enddo
c
      return
      end
c
c verify triangle rule
      FUNCTION NOTRI(K,L,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      NOTRI=-1
      IF(MOD((K+L+M),2).EQ.1) GOTO 10
      IF((K+L-M).LT.0) GOTO 10
      IF((K-L+M).LT.0) GOTO 10
      IF((M+L-K).LT.0) GOTO 10
      NOTRI=1
   10 RETURN
      END
c $Header:$
c
c complex spherical harmonic
c phase convenetion: cylm(l,m,x,phi) =  (-1)^m cylm(l,m,x,phi)
c x = cos(theta)
c
c Martin Fuchs, FHI 16-06-1995
c
      complex*16 function cylm(l,m,x,phi)

      implicit real*8 (a-h,o-z)
      complex*16 a,conjg
      parameter (nmx=16,lmx=nmx-1,pi4=12.5663706143591725d0)
      dimension ll(nmx)
      data ll /1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31/

c normalization factor
      if(l .gt. lmx .or. x*x .gt. 1. .or. m*m .gt. l*l) then
        print *, l,lmx
        print *, x**2,1
        print *, m**2,l**2
        stop 'cylm --- bad argument' 
      endif

      mm=abs(m)
      fin=1
      do i=l-mm+1,l+mm
        fin=fin*i
      enddo
      fnorm=sqrt(ll(l+1)/(fin*pi4))

      cylm=fnorm*plgndr0(l,mm,x)*exp(cmplx(0,mm*phi))

      if(m .lt. 0) then
        cylm=conjg(cylm)
        if(mod(mm,2) .ne. 0) cylm=-cylm
      endif

      return
      end
c
c
c associated legendre polynomial
c taken from
c        numerical recipies,
c        w.h. preuss et al.,
c        cambridge university press, new york, 1986
c

      real*8 function plgndr0(L,M,X)

      implicit real*8 (a-h,o-z)

      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.d0)PAUSE 'bad arguments'
      PMM=1.d0
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.d0-X)*(1.d0+X))
        FACT=1.d0
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.d0
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        plgndr0=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          plgndr0=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          plgndr0=PLL
        ENDIF
      ENDIF
      RETURN
      END
c
c associated legendre polynomial 1st derivative
c
      real*8 function plgndr1(l,m,x)
      
      implicit real*8 (a-h,o-z)
      external plgndr0

      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.d0)PAUSE 'bad arguments'

      plgndr1=(l+1)*(x*plgndr0(l,0,x)-plgndr0(l+1,0,x))/(1.d0-x*x)

c     plgndr1=0.d0
c     if(m+1 .le. l) plgndr1=plgndr0(l,m+1,x)
c     if(m-1 .ge. 0) plgndr1=plgndr1+
c    +  (l+m)*(l-m+1)*plgndr0(l,m-1,x)
c     plgndr1=0.5d0*plgndr1/(1.d0-x*x)**2
      
      return
      end
c
       
        
      
