c $Header:$
************************************************************************
c core polarization potential of Igel-Mann et al.
c
c input
c z       atomic number
c mmax    radila mesh range
c r		radial mesh
c
c output
c vcp     polarization potential
c vo1()...  potential arrays, vi() is updated output 
c***********************************************************************
      subroutine ecp(z,mmax,r,vcp)
c
      implicit none
      integer  i,mmax,nzmax,iz
      parameter (nzmax=100)
      real*8   z
      real*8   r(mmax),vcp(mmax),alpha(nzmax),delta(nzmax)
      include  'parameter.h'
c
      do i=1,nzmax
        delta(i) = 1.d1
        alpha(i) = 0.d0
      enddo
      alpha(13)  = 0.2649d0 
      delta(13)  = 1.1600d0 
      alpha(14)  = 0.1624d0 
      delta(14)  = 1.4700d0 
      alpha(31)  = 1.2400d0 
      delta(31)  = 0.5362d0 
      alpha(32)  = 0.7628d0 
      delta(32)  = 0.7105d0 
      alpha(33)  = 0.5096d0 
      delta(33)  = 0.8863d0 
c shirley martin parameters
c     delta(14) = 1.25d0
c     delta(31) = 1.00d0
c     delta(32) = 1.05d0
c     delta(33) = 1.14d0
      iz = aint(z)
      write(iofhipp,*) 'ECP -- effective core potential parameters:'
      write(iofhipp,*) 'ecp - iz alpha delta',iz,alpha(iz),delta(iz)
      do i=1,mmax
      vcp(i)=-0.5d0*alpha(iz)*(1.d0-exp(-delta(iz)*r(i)**2))**2/r(i)**4
      enddo
c
      return
      end
