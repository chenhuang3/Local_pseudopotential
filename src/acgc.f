c     program testacgc
c     implicit real*8 (a-h,o-z)
c     read(5,*) l1,l2,l3
c     print *, acgc(l1,l2,l3)
c     end

      real*8 function acgc(l1,l2,l3)

      implicit real*8 (a-h,o-z)
      parameter(mx=60)
      dimension fac(0:mx)
      logical ifc
      data    ifc/.true./
      save    ifc,fac

      t(i)=fac(i/2)/sqrt(fac(i))
      hatl(i)=dble(2*i+1)

      if(ifc) then
        ifc=.false.
        fac(0)=1.d0
        do i=1,mx
          fac(i)=dble(i)*fac(i-1)
        enddo
      endif

      if(l1+l2+l3 .gt. mx) stop 'acgc --- l1+l2+l3'
      if(notri(l1,l2,l3) .gt.0) then
        acgc=sqrt(hatl(l3)/(l1+l2+l3+1))
     &         *t(l1+l2+l3)/
     &           (t(l1+l2-l3)*t(l1-l2+l3)*t(l2+l3-l1))
        if(mod((l1+l2-l3)/2,2) .ne. 0.) acgc=-acgc
      else
        acgc=0.d0
      endif

      return
      end
