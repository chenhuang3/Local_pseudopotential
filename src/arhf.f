c $Header: /f21/a/fuchs/Psp/Dkef/Ds/RCS/arhf.f,v 1.1 95/08/22 17:41:35 fuchs Exp Locker: fuchs $
c
c compute a term of the multipole expansion of the exchange potential
c
c yabc(x) = integral(y^{2}dy)[ R_{la}(y)*R_{lb}(y)* min[y,x]^{lc}/max[y,x]^{lc+1} ]
c where R_{l}(y)=u_{l}(y)/y are radial wavefunctions w/ angular momentum quantum no. l
c
c input
c mmx ............ upper limit index on radial grid
c r() ............ logarithmic radial grid
c la ............. dto
c lb ............. dto
c lc ............. dto
c ua() ........... u_{la}
c ub() ........... dto
c uc(),wc() ...... workspace
c
c output
c yabc() ......... dto
c
c Martin Fuchs 05-08-1995
c
      subroutine arhf(mx,r,la,lb,lc,ua,ub,uc,wc,yabc)
c
      implicit real*8 (a-h,o-z)
c
      dimension r(mx),ua(mx),ub(mx),uc(mx),yabc(mx),wc(mx)
      parameter(zero=0.d0)

      fi(i)=ua(i)*ub(i)*uc(i)
      fo(i)=ua(i)*ub(i)*wc(i)

      al=0.1d0*log(r(11)/r(1))

c evaluate radial matrix element w/ trapezoidal rule 
      k1=lc
      k2=-(lc+1)
      k3=-(lc+1)
      do i=1,mx
        uc(i)=r(i)**k1
        wc(i)=r(i)**k2
      enddo

      do i=1,mx
        yabc(i)=zero
        do j=1,i-1
          yabc(i)=yabc(i)+fi(j)*r(j)
        enddo
        yabc(i)=al*(.5d0*fi(i)*r(i)+yabc(i))*r(i)**k3
        aux=zero
        do j=i+1,mx-1
          aux=aux+fo(j)*r(j)
        enddo
        aux=al*(aux+.5d0*(fo(i)*r(i)+fo(mx)*r(mx)))*r(i)**lc
        yabc(i)=yabc(i)+aux
      enddo

      return
      end

c w/ as deq
c     sb=ub(3)/r(3)**(lb+1)
c     sc=uc(3)/r(3)**(lc+1)
c     ld=la+lb+lc+4
c     wc(1)=sb*sc*r(1)**ld
c     wc(2)=sb*sc*r(2)**ld
c     wc(3)=sb*sc*r(3)**ld

c     do i=1,mx
c       wa(i)=ub(i)*uc(i)*al*r(i)**(1+la)
c     enddo
c     do i=3,mx-1
c       wc(i+1)=wc(i)+aio(wa,i)
c     enddo

c     ld=-la+lb+lc-1
c     wb(1)=sb*sc*r(1)**ld
c     wb(2)=sb*sc*r(2)**ld
c     wb(3)=sb*sc*r(3)**ld

c     do i=1,mx
c       wa(i)=ub(i)*uc(i)*al/r(i)**la
c     enddo
c     do i=3,mx-1
c       wb(i+1)=wb(i)+aio(wa,i)
c     enddo
c     wbs=wb(mx)
c     do i=1,mx
c       wc(i)= wc(i)/r(i)**(la+1)+ (wbs-wb(i))*r(i)**la
c     enddo

c     print *, 'arhf - wbs',wbs,wb(1)
c     do i=1,mx,5
c       write(69,*) r(i),yabc(i),wc(i)
c     enddo
c     close(69)
