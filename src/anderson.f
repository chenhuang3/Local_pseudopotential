c $Header:$
************************************************************************
c updating potentials in self-consitency iteration
c
c input
c bl		mixing parameter
c it		iteration number
c mmax	number of relvant meshpoints points
c r()		here: radial mesh
c vo1()...  potential arrays, vi() is updated output 
c
c original version by D.R. Hamann, gncpp
c***********************************************************************
      subroutine anderson(bl,it,mmax,r,vo1,vo,vi1,vi)
c
      implicit none
      integer  i,it,mmax
      real*8   bl,thl,sn,sd,vn,rl,rl1,dr
      real*8   r(mmax),vo(mmax),vo1(mmax),vi(mmax),vi1(mmax)
c
c
      thl=0.0d0
      if(it .gt. 1) then
        sn=0.0d0
        sd=0.0d0
        do 120 i=1,mmax
          rl=vo(i)-vi(i)
          rl1=vo1(i)-vi1(i)
          dr=rl-rl1
          sn=sn + rl*dr*r(i)**2
          sd=sd + dr*dr*r(i)**2
  120   continue
        thl=sn/sd
      endif
      do 140 i=1,mmax
        vn=(1.0d0-bl)*((1.0d0-thl)*vi(i) + thl*vi1(i))
     &   + bl*((1.0d0-thl)*vo(i) + thl*vo1(i))
        vi1(i)=vi(i)
        vo1(i)=vo(i)
        vi(i)=vn
  140 continue
c
      return
      end
