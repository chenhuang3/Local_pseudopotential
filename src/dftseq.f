c  $Header:$
c***********************************************************************
c
c integrates radial (pauli-type scalar-relativistic w/ spin-orbit av.) 
c equation on a c logarithmic mesh w/ fixed/variable effective mass 
c
c input
c mode ...... 1 is for full potential bound state (scalar relativistic)
c             2 is for pseudopotential bound state (nonrelativistic)
c             3 is for full potential to find log derivative (sc.rel.)
c             4 is for pseudopotential to find energy which produces
c               specified log derivative (fake bound state)
c             5 is for pseudopotential to produce wavefunction beyond
c               radius used for pseudopotential construction
c             6 like 3 but nonrelativistic
c             7 like 4 
c             8 like 2 but not terminating if no bound state is found
c               returning zero eigenvalue instead
c z ......... nuclear charge
c mmax ...... number of radial grid points
c r() ....... radial grid
c n ......... principal quantum number
c l ......... angular momentum qn
c efm ....... effective mass
c uld ....... logarithmic derivative at mch
c v() ....... potential v*Psi
c v1() ...... potential v1*grad(Psi)
c v2() ...... potential 1/v2 * laplace(Psi)
c
c output
c nin ....... outermost gridpoint used
c mch ....... classical turning point
c e ......... eigenvalue
c u() ....... wave function
c up() ......               1st derivative on log grid
c upp() .....               2nd derivative on log grid
c 
c original version by D.R. Hamann, gncpp
c
c MF FHI Berlin 
c***********************************************************************
c
      subroutine dftseq(mode,z,mmax,r,n,l,efm,v,v1,v2,
     &           nin,mch,uld,e,u,up,upp)
c
      implicit real*8 (a-h,o-z)
c
      include 'parameter.h'
      include 'default.h'

      parameter (itmax=200)
      parameter (fss_i=137.0359896d0)
c
      logical trap,trel
      common/err/trap(10)
c
      dimension r(mx),u(mx),up(mx),upp(mx)
      dimension cf(mx),dv(mx),fr(mx),frp(mx),v(mx),v1(mx),v2(mx),
     1          x1(mx),x2(mx),x4(mx),erec(itmax)

c
c convergence factor for solution of schroedinger eq.  if calculated
c correction to eigenvalue is smaller in magnitude than epsae times
c the magnitude of the current guess, the current guess is not changed.
c
c     epsae=1.0d-9
c
c relativistic - non-relativistic switch
c
      amesh=r(2)/r(1)
      al=log(amesh)
      sls=l*(l+1)

      if(mode .eq. 6 .or. mode .eq. 7) then
        ics=max(1,mch)
        ico=max(1,mch)
      else
        ics=mmax
        ico=mch
      endif
        
      if(mode .eq. 1 .or. mode .eq. 3) then
        trel=.true.
        fss=(1.0d0/fss_i)**2
        if(l .eq. 0) gamma=sqrt(1.0d0-fss*z**2)
        if(l .gt. 0) gamma=(l*sqrt(l**2-fss*z**2) +
     &   (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
      else
        trel=.false.
        fss=1.0d-20
        gamma=l+1
      end if
      imod=mode
      if(imod .eq. 6) imod = 3
      if(imod .eq. 7) imod = 4
      if(imod .eq. 8) imod = 2
!d     if(mode .ne. 1 .and. mode .ne. 2) then
!d       write(ie,*) '& dftseq- mode,imod,fss,gamma',mode,imod,fss,gamma
!d     endif
c   
      if(imod .eq. 1 .or. imod .eq. 2) then
        emax=v(mmax)+0.5d0*sls/r(mmax)**2
        emin=0.0d0
        do 6 i=1,mmax
    6     emin=dmin1(emin,v(i)+0.5d0*sls/r(i)**2)
        if(e .gt. emax) e=1.25d0*emax
        if(e .lt. emin) e=0.75d0*emin
        if(e .gt. emax) e=0.5d0*(emax+emin)
      else if(imod .eq. 4) then
        emax=e +  10.0d0
        emin=e - 10.0d0
      end if
c
c initialize arrays and remove leftover garbage
      do i=1,4
        u(i)=0.d0
        up(i)=0.0d0
        upp(i)=0.0d0
      enddo
      do i=1,mmax
        x1(i)=al*r(i)
        x2(i)=x1(i)*x1(i)
      enddo

      if(trel) then
c calculate dv/dr for darwin correction
        dv(1)=(-50.d0*v(1)+96.d0*v(2)-72.d0*v(3)+32.d0*v(4)
     &         -6.d0*v(5))/(24.d0*x1(1))
        dv(2)=(-6.d0*v(1)-20.d0*v(2)+36.d0*v(3)-12.d0*v(4)
     &         +2.d0*v(5))/(24.d0*x1(2))

        do i=3,ics
          dv(i)=(2.d0*v(i-2)-16.d0*v(i-1)+16.d0*v(i+1)
     &           -2.d0*v(i+2))/(24.d0*x1(i))
        enddo
      else
c nonrelativistic coefficient arrays for u (fr) and up (frp)
        do i=1,ics
          fr(i)=-2*x2(i)*v1(i)*v2(i)
          frp(i)=2*x1(i)*v1(i)*v2(i)
        enddo
      endif

c return point for bound state convergence
      aux=al*al*sls
      nint=0
   10 nint=nint+1
      erec(nint)=e
        if(nint .gt. itmax) then
          if(mode .eq. 8) then
            write(ie,*) 
     1        '& dftseq - no bound state found (iter): e=>0 n l',n,l
            e=0.d0
            return
          endif
          write(ie,*) '& dftseq - stop: convergence error'
          write(ie,*) '&   nint n l                ',nint,n,l
          write(ie,*) '&   eigenvalue log to file err.dftseq_e'
          write(ie,*) '&   potential  log to file err.dftseq_v'
          open(10,file='err.dftseq_e',status='unknown')
          do i=1,10
            write(10,*) i,erec(i)
          enddo
          do i=nint,1,-1
            write(10,*) i,erec(i)
          enddo
          close(10)
          open(10,file='err.dftseq_v',status='unknown')
          do i=1,ics
            write(10,*) r(i),v(i)
          enddo
          close(10)
          trap(4)=.true.
          stop
        endif

c coefficient array for u in differential eq.
        do i=1,ics
          cf(i)=aux+2*(v(i)-e)*v2(i)*x2(i)
        enddo

c relativistic coefficient arrays for u (fr) and up (frp)
        if(trel) then
          do i=1,ics
            fr(i)=x2(i)*fss*(-(v(i)-e)**2 + 0.5d0*dv(i)/
     &      (r(i)*(1.0d0+0.5d0*fss*(e-v(i)))))
            frp(i)=-x1(i)*0.5d0*fss*dv(i)/(1.0d0+0.5d0*fss* (e-v(i)))
          enddo
        endif

c
c find classical turning point for matching
        if(imod .eq. 1 .or. imod .eq. 2) then
          do 30 i=ics,2,-1
            if(cf(i-1) .le. 0.d0 .and. cf(i) .gt. 0.d0) then
              mch=i
              ico=mch
              go to 40
            end if
   30     continue
c
c error trap
          if(mode .eq. 8) then
            write(ie,*) 
     1        '& dftseq - no bound state found (ctp): e=>0 n l',n,l
            e=0.d0
            return
          endif
          write(ie,*) '& dftseq - stop: no classical turning point'
          write(ie,*) '&    nint n l e                      ',nint,n,l,e
          write(ie,*) '&    0-order coefficient to file err.dftseq_tp'
          open(10,file='err.dftseq_tp',status='unknown')
          do i=1,ics
            write(10,*) r(i),cf(i)/r(i)**2
          enddo
          close(10)
          trap(5)=.true.
          stop 
   40     continue
        else
          nin=mch
        end if

c start wavefunction with series
        do 50 i=1,4
          u(i)=r(i)**gamma
          up(i)=al*gamma*r(i)**gamma
   50     upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
c
cmf l=0: hydrogen-like start (improves eigenvalues and wavefunctions)
c
        if(l .eq. 0 .and. imod .eq. 2 .or. imod .eq. 1) then
cmf          z=z*efm
          do i=1,4
            u(i)=r(i)**gamma*(1.d0-z*r(i)+0.5d0*(z*r(i))**2)
            up(i)=al*gamma*r(i)**gamma*((1.d0-z*r(i)+0.5d0*(z*r(i))**2))
     &          +al*r(i)**gamma*(-z*r(i)+(z*r(i))**2)
            upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
          enddo
cmf         z=z/efm
        endif
c
c outward integration using predictor once, corrector
c twice
        node=0
c
        do 70 i=4,ico-1
          u(i+1)=u(i)+aeo(up,i)
          up(i+1)=up(i)+aeo(upp,i)
          do 60 it=1,2
            upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
            up(i+1)=up(i)+aio(upp,i)
   60     u(i+1)=u(i)+aio(up,i)
          if(u(i+1)*u(i) .le. 0.0d0) node=node+1
   70   continue
c
        uout=u(mch)
        upout=up(mch)
c
        if(node-n+l+1 .eq. 0 .or. imod .eq. 3 .or. imod .eq. 5) then
c
          if(imod .eq. 1 .or. imod .eq. 2) then
c
c start inward integration at 15*classical turning
c point with simple exponential
            nin=mch+2.7d0/al
            if(nin .gt. mmax-4) nin=mmax-4
            xkap=sqrt(sls/r(nin)**2 + 2.0d0*(v(nin)-e)*efm)
c
            do 110 i=nin,nin+4
              u(i)=exp(-xkap*(r(i)-r(nin)))
              up(i)=-r(i)*al*xkap*u(i)
  110         upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
c
c integrate inward
c
            do 130 i=nin,mch+1,-1
              u(i-1)=u(i)+aei(up,i)
              up(i-1)=up(i)+aei(upp,i)
              do 130 it=1,2
                upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
                up(i-1)=up(i)+aii(upp,i)
  130           u(i-1)=u(i)+aii(up,i)
c
c scale outside wf for continuity
            sc=uout/u(mch)
c
            do 150 i=mch,nin
              upp(i)=sc*upp(i)
              up(i)=sc*up(i)
  150         u(i)=sc*u(i)
c
            upin=up(mch)
c
          else
c
            upin=uld*uout
c
          end if
c
c perform normalization sum
c
          ro=r(1)/sqrt(amesh)
          sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
c
          do 160 i=1,nin-3
  160       sn=sn+al*r(i)*u(i)**2
c
          sn=sn + al*(23.0d0*r(nin-2)*u(nin-2)**2
     &              + 28.0d0*r(nin-1)*u(nin-1)**2
     &              +  9.0d0*r(nin  )*u(nin  )**2)/24.0d0
c
          cn=1.0d0/sqrt(sn)
          uout=cn*uout
          upout=cn*upout
          upin=cn*upin
c
          do 180 i=1,nin
            upp(i)=cn*upp(i)
            up(i)=cn*up(i)
  180       u(i)=cn*u(i)
          do 190 i=nin+1,mmax
            upp(i)=0.d0
            up(i)=0.d0
  190       u(i)=0.0d0
c
c exit for fixed-energy calculation
c
          if(imod .eq. 3 .or. imod .eq. 5) return
c
c perturbation theory for energy shift
          de=0.5d0*uout*(upout-upin)/(al*r(mch))
c
c convergence test and possible exit
c
          if(abs(de) .ge. max(abs(e),0.2d0)*epsae) then
c
            if(de .gt. 0.0d0) then 
              emin=e
            else
              emax=e
            end if
            e=e+de
            if(e .gt. emax .or. e .lt. emin) e=0.5d0*(emax+emin)
c
c loop back to converge e
c
            go to 10

          else

            return

          endif
c
        else if(node-n+l+1 .lt. 0) then
c too few nodes
          emin=e
          e=0.5d0*(emin+emax)
          go to 10
c
        else
c too many nodes
          emax=e
          e=0.5d0*(emin+emax)
          go to 10
        end if
      continue
c
c
      end
c
