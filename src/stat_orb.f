      real*8 function stat_orb(l,al,r,u,up)

      implicit none
      integer ll,l
      real*8  al,r,u,up
      real*8  g0,g1,tt

      ll = l*(l+1)
      g0 = u/r
      g1 = (up/(al*r)-g0)/r
      tt = g1*g1
      if(ll .ne. 0) tt = tt + ll/r**2 * g0*g0
      stat_orb = 0.5d0*tt
 
      return
      end
