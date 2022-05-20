c $Header:$
c array manipulations
c
c b = a
      subroutine dcpv(mx,ml,mh,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(mx),b(mx)
      do i=ml,mh
        b(i)=a(i)
      enddo
      return
      end

c
c c = a + b
      subroutine dadv(mx,ml,mh,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(mx),b(mx),c(mx)
      do i=ml,mh
        c(i)=a(i)+b(i)
      enddo
      return
      end

c
c c = a*b + d 
      subroutine dlcv(mx,ml,mh,a,b,d,c)
      implicit real*8 (a-h,o-z)
      dimension b(mx),d(mx),c(mx)
      do i=ml,mh
        c(i)=a*b(i)+d(i)
      enddo
      return
      end
c
c c = a - b
      subroutine dsuv(mx,ml,mh,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(mx),b(mx),c(mx)
      do i=ml,mh
        c(i)=a(i)-b(i)
      enddo
      return
      end
c
c a = 0
      subroutine dnuv(mx,il,ih,const,a)
      implicit real*8 (a-h,o-z)
      dimension a(mx)
      do i=il,ih
        a(i)=const
      enddo
      return
      end
c
c a DOT b
      real*8 function dspv(mx,il,ih,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(mx),b(mx)
      dspv=0.d0
      do i=il,ih
        dspv=dspv+a(i)*b(i)
      enddo
      return
      end
c
c c = a*b
      subroutine dscv(mx,il,ih,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension b(mx),c(mx)
	do i=il,ih
        c(i)=a*b(i)
      enddo
      return
      end
c
      subroutine dsqv(mx,il,ih,a,b,c)
      implicit real*8 (a-h,o-z)
      dimension a(mx),b(mx),c(mx)
	do i=il,ih
        c(i)=a(i)*b(i)
      enddo
      return
      end
c
      subroutine dextv(mx,kl,kh,b,bl,bh)
      implicit real*8 (a-h,o-z)
      dimension b(mx)
      parameter (epsl=-1.d20,epsh=1.d20)

      il=kl
      ih=kh
      is = 1
      if(kl .gt. kh) then
        il = kh 
        ih = kl 
        is = -1
      endif
      bl=epsh
      bh=epsl
      do i=il,ih,is
        bl=min(b(i),bl)
        bh=max(b(i),bh)
      enddo

      return
      end
        
        
      
