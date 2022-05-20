c $Header:$
c***********************************************************************
c generate logarithmic mesh 
c
c mmax    index of outmost meshpoint
c amesh   ratio of consecutive points
c al      logarithmic increment
c z       atomic number (scale start point)
c***********************************************************************
c
      subroutine logmesh(mmax,z,amesh,al,r)
c
      implicit none
      include  'parameter.h'
      include  'default.h'
      integer  mmax,i
      real*8   z,amesh,al,cmax
      real*8   r(mx)
c
c  minimum rmin and spacing ainc from include default.h
      r(1)=dble(rmin)/z
      amesh=dble(ainc)
      al=log(amesh)      
      cmax=dble(rmax)/rmin
      mmax=log(cmax*z)/al
      if(mod(mmax,2) .eq. 0) mmax=mmax+1
      if(mmax .gt. mx) stop 'logmesh - mesh too large'
      do i=2,mmax
        r(i)=amesh*r(i-1)
      enddo
      return
      end
c
