c $Header:$
c***********************************************************************
c
c Fourier-Bessel-transform (FBT) of a real function given on a log mesh
c and inverse
c
c input
c mx ........... real space mesh size
c r() .......... real space mesh
c gp() ......... real space function on mesh
c lb ........... angular momentum index for FBT
c ekmx ......... cutoff energy E = k^2 (Ry a.u.) 
c delk ......... mesh increment in reciprocal space 
c                inverse transformation done on a linear mesh 
c mxe .......... maximum dimension for reciprocal space mesh
c iou .......... output unit, if 0 no output
c
c output
c gp() ......... inverse transform 
c fk() ......... reciprocal space mesh
c gk() ......... FBT coefficients
c func() ....... workspace
c
c***********************************************************************
c	
      subroutine dfbt(mx,r,gp,lb,ekmx,delk,mxe,kmx,fk,gk,func,iou,cnorm)
c
      implicit none
      
      integer  mx,lb,mxe,iou,i,k,kmx
      real*8   fmom,bessj,ekmx,delk,cnorm,pi4,gtrf,gint,al,fkmx
      real*8   r(mx),gp(mx),fk(mxe),gk(mxe),func(mx)
      external fmom,bessj

c integrands
      gtrf(i,k)=gp(i)*bessj(lb,r(i)*fk(k))/r(i)
      gint(i,k)=gk(k)*bessj(lb,r(i)*fk(k))/r(i)

      pi4=16.0*atan(1.0d0)
      cnorm=8.0/pi4
      
c cutoff energy in rydberg
      fkmx=sqrt(ekmx)
      fk(1)=0.d0
      do k=2,mxe
        fk(k)=delk*(k-1)
        kmx=k
        if(fk(k) .ge. fkmx) goto 10
      enddo
   10 continue
      if(kmx .gt. mxe) stop '& dfbt - stop: k mesh overflow'

c Fourier-Bessel-Transform 
c radial mesh is assumed logarithmic
      al=log(r(2)/r(1))
      do k=1,kmx
        do i=1,mx
          func(i)=gtrf(i,k)
        enddo
        gk(k)=cnorm*fmom(0,mx,al,1.d0,r,func)
      enddo

c inverse transform
      do i=1,mx
        gp(i)=0.5*(gint(i,1)+gint(i,kmx))
        do k=2,kmx-1
          gp(i)=gp(i)+gint(i,k)
        enddo
        gp(i)=gp(i)*delk
      enddo

      if(iou .gt. 0) then
        do k=1,mx
          write(iou,'(e21.14,1x,e21.14)') r(k),gp(k)
        enddo
        close(iou)
      endif
      
      return
      end
c
