c $Header:$
c-----------------------------------------------------------------------
c spectral analysis of radial schroedinger equation  with a fully 
c separable nonlocal potential (given by projector function).
c the integrodifferential equation is transformed into a matrix equation
c using a basis of laguerre polynomials, resulting in a standard 
c eigenvalue problem 
c
c Hartree a.u.
c                        !!! using lapack !!!
c
c input
c l ............ angular momentum
c mmax ......... maximum index on logarithmic radial grid
c r() .......... radial grid
c vsc() ........ local potential (here screening part)
c vsl() ........ semilocal potential 
c nnl .......... number of projectors present
c cnl() ........ inverse of average nonlocal potential (strength parameter)
c vnl() ........ projector functions
c 
c output
c enl() ........ bound state spectrum
c unl() ........ real space eigenfunctions
c
c note: radial grids commonly used in pp generation might be too coarse
c       to ensure a sufficiently accurate transformation into the 
c       polynomial representation 
c
c Martin Fuchs, FHI, 07-1995
c-----------------------------------------------------------------------
c 
      subroutine klbyii(l,mmax,r,vsc,vsl,nnl,cnl,vnl,enl,unl)

      implicit real*8 (a-h,o-z)
      character*8 sav2101

      include 'parameter.h'
      include 'gauss.h'

      logical tfirst,tnl
c nmx: basis size, must be consistent with imx: order of gauss quadrature
      parameter (mp=1,nmx=30,nmx1=nmx*(nmx+1)/2,zero=0.d0)
      dimension r(mx),vsc(mx),vsl(mx),cnl(mp),vnl(mx,mp),vloc(mx),
     &          g0(mx,nmx),f0(imx,nmx),f1(imx,nmx),
     &          hu(nmx,nmx),hessl(nmx1),he(nmx),hz(nmx,nmx),
     &          enl(ms),unl(mx,ms),ht(nmx,nmx),
     &          htl(nmx,nmx),tnl(4),hnl(nmx,4),
     &          rwk1(8*nmx),iwk1(5*nmx),iwk2(nmx)
      external  geta,bas0,bas1,gaussq

      save tfirst,g0,g1,ht,tnl
      data tfirst /.true./
      data tnl    /.true.,.true.,.true.,.true./

      al=.1d0*log(r(11)/r(1))

c construct basis
      if(tfirst) then
!d       t1=cputime(0.d0)
        tfirst=.false.
        a = geta(1)
        do n1=1,nmx
          do i=1,mmax
            g0(i,n1)=bas0(n1-1,r(i))
          enddo
        enddo
!d       t1=cputime(t1)
!d       write(ie,*) '& klbyii - time basis 0       ',t1
!d       t1=cputime(t1)
        do n1=1,nmx
          do i=1,imx
            f0(i,n1)=bas0(n1-1,x(i)/a)*exp(0.5*x(i))/x(i)
            f1(i,n1)=bas1(n1-1,x(i)/a)*exp(0.5*x(i))
          enddo
        enddo

!d       t1=cputime(t1)
!d       write(ie,*) '& klbyii - time basis 1       ',t1
!d       t1=cputime(t1)
            
c kinetic energy operator
        do n2=1,nmx
          do n1=1,n2
            ht(n1,n2)=0.d0
            htl(n1,n2)=0.d0
c explcit gauss quadrature
            do i=1,imx
              ht(n1,n2)=ht(n1,n2)+w(i)*f1(i,n1)*f1(i,n2)
              htl(n1,n2)=htl(n1,n2)+w(i)*f0(i,n1)*f0(i,n2)
            enddo
            ht(n1,n2)=0.5*ht(n1,n2)/a
            htl(n1,n2)=a*htl(n1,n2)
c           ht(n1,n2)=0.5*gaussq(gtlfgq,n1-1,n2-1)
c           htl(n1,n2)=gaussq(gtlfgq,n1-1,n2-1)
          enddo
        enddo
!d       t1=cputime(t1)
!d       write(ie,*) '& klbyii - time initialization',t1
      endif

      ll=l+1
      if(nnl .gt. 0) then
!d       t1=cputime(0.d0)
        if(tnl(ll)) then
          tnl(ll)=.false.
          do i=1,nnl
            do n3=1,nmx
              hnl(n3,ll)=gltfmv(mmax,al,r,g0(1,n3),vnl(1,i),r)
            enddo
          enddo
        endif
!d       t1=cputime(t1)
!d       write(ie,*) '& klbyii - time nonlocal init',t1
      endif
          
c screened (semi-) local potential
      fll=.5d0*l*(l+1)
      do i=1,mmax
        vloc(i)=vsc(i)+vsl(i)
      enddo

c hamiltonian matrix
!d     t1=cputime(0.d0)
      do n2=1,nmx
        do n1=1,n2
          hu(n1,n2)=ht(n1,n2)+fll*htl(n1,n2)
     &             +gltfmv(mmax,al,r,g0(1,n1),vloc,g0(1,n2))
          do i=1,nnl
            hu(n1,n2)=hu(n1,n2)+cnl(i)*hnl(n1,ll)*hnl(n2,ll)
          enddo
          hu(n2,n1)=hu(n1,n2)
        enddo
      enddo
!d     t1=cputime(t1)
!d     write(ie,*) '& klbyii - time matrix setup ',t1

c calculate spectrum and eigenvectors
      n3=0
      do n1=1,nmx
        he(n1)=0.d0
        do n2=n1,nmx
          n3=n3+1
          hessl(n3)=hu(n2,n1)
        enddo
      enddo
c lapack 
      info=0
      call dspevx('N','I','L',nmx,hessl,-1.d4,1.d1,1,4,1.d-6,
     &     i,he,hz,nmx,rwk1,iwk1,iwk2,info)
      if(info .ne. 0) then
        write(ie,*) '& klbyii - warning: eigenvalues not converged'
        write(ie,*) '& klbyii - he(i)',(iwk2(i),i=1,4)
      endif

c real space eigenfunctions
      i=1
10    if(he(i) .lt. 0.d0 .and. i .lt. ms) then
        enl(i)=he(i)
        do j=1,mmax
          unl(j,i)=0.0d0
          do n1=1,nmx
            unl(j,i)=unl(j,i)+hz(n1,i)*g0(j,n1)
          enddo
        enddo
      i=i+1
      goto 10
      endif

      return
      end

c     real*8 function gltfgq(m,n,y)
c     implicit real*8 (a-h,o-z)
c     include 'basis.h'
c     external bas1
c     x=y/a
c     gltfgq=bas1(m,x)*bas1(n,x)*exp(y)/a
c     return
c     end

c     real*8 function gltflq(m,n,y)
c     implicit real*8 (a-h,o-z)
c     include 'basis.h'
c     external bas0
c     x=y/a
c     gltflq=a*bas0(m,x)*bas0(n,x)*exp(y)/y**2
c     return
c     end

!d     REAL*8 FUNCTION CPUTIME(TIME)
!d     REAL*8 TIME,ETIME
!d     integer mclock
!d     REAL*8 TARRAY(2)
!d     external mclock
!d     CPUTIME = mclock()*0.01 -time
!d     END
