c $Header:$
c range of sic potential
c LDA-SIC
      subroutine spinorbit(mmax,nn,ll,ry2,eig,r,dorb,vorb)

      implicit none
      include  'parameter.h'
      integer  mmax,i,ll,nn
      real*8   ry2,range,al,dmelm,egauge,delta,converge
      real*8   r(mx),dorb(mx),vorb(mx),vso(mx)
      real*8   eig,fss,ls_up,ls_dn,e_up,e_dn,esplit
      parameter(fss=137.0359896d0**(-2))
      external dmelm

      al=0.1d0*log(r(11)/r(1))
      do i=2,mx
        vso(i)=((vorb(i)-vorb(i-1))/(r(i)-r(i-1)))/r(i)
      enddo
      vso(1)=vso(2)*r(2)/r(1)
      esplit=0.25*fss*dmelm(mmax,al,r,vso,dorb)

      ls_up = ((ll+0.5)*(ll+0.5) - ll*ll - 0.25)
      ls_dn = ((ll-0.5)*(ll-0.5) - ll*ll - 0.25)

      e_up  = eig+esplit*ls_up
      e_dn  = eig+esplit*ls_dn

      write(ie,*) 'n=',nn,' l=',ll, ' eig=',eig*ry2
      write(ie,*) 'j = l+1/2: eig = ', e_up*ry2
      write(ie,*) 'j = l-1/2: eig = ', e_dn*ry2
      write(ie,*) 'so splitting   = ', esplit*(ls_up-ls_dn)*ry2

      return
      end
