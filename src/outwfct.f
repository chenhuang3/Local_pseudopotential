c $Header:$
c***********************************************************************
c output subroutines for pseudopotential program
c***********************************************************************
c
c densities
      subroutine outcore(iu,mx,mstep,r,dc,dcp,dcpp)
      implicit real*8 (a-h,o-z)
      dimension r(mx),dc(mx),dcp(mx),dcpp(mx)
      open(iu)
      do i=1,mx,mstep
        write(iu,700) r(i),dc(i),dcp(i),dcpp(i)
      enddo
      close(iu)
  700 format(e20.14,3(1x,e20.14))

      return
      end
c
c pseudopotentials and wavefunctions
      subroutine outpot(istart,mx,mmax,llmx,cval,rc,r,ups,vps,vscr)
      implicit none
      integer  istart,mx,mmax,llmx,iu,iv,i,ll
      real*8   amesh,cval,rc(llmx),r(mx),ups(mx,*),vps(mx,*),vscr(mx,*)

      amesh=r(2)/r(1)
      do ll=1,llmx
        iu=istart+ll-1
        iv=istart+5+ll-1
        write(iu,'(a,1x,i5,1x,e20.14,1x,i1,1x,f8.4,1x,e20.14)')
     1    '#',mmax,amesh,ll-1,rc(ll),cval
        write(iv,'(a,1x,i5,1x,e20.14,1x,i1,1x,f8.4,1x,e20.14)')
     1    '#',mmax,amesh,ll-1,rc(ll),cval
        do i=1,mmax
c ionic
          write(iu,700) i,r(i),ups(i,ll),vps(i,ll)
c screened
          write(iv,700) i,r(i),ups(i,ll),vscr(i,ll)
        enddo
        close(iu)
        close(iv)
      enddo
  700 format(i4,1x,e20.14,2(1x,e20.14))
  710 format(e20.14,1x,e20.14)

      return
      end
c
c wavefunctions
      subroutine out38(mx,nmn,nmx,ninu,n,l,r,u)
      implicit real*8 (a-h,o-z)
      character*1 sn,sl,sm
      dimension n(nmx),l(nmx),ninu(nmx),r(mx),u(mx,nmx)
c
      open(38)
      do i=nmx,nmn,-1
        call labels(n(i),l(i),1,sn,sl,sm)
        write(38,'(2a2,4(1x,i2),1x,i5)')
     +  '# ',sn//sl,i,n(i),l(i),0,ninu(i)
        do j=1,ninu(i),max(ninu(nmx)/100,1)
          write(38,'(e12.6,1x,e12.6)') r(j),u(j,i)
        enddo
      enddo
      close(38)
c
      return
      end
c
c pseudo and ae valence wavefunctions
      subroutine out39(mx,lmax,ninu,np,rc,r,ups,uae)
      implicit real*8 (a-h,o-z)
      character*1 sn,sl,sm
      dimension np(*),ninu(*),rc(*),r(mx),ups(mx,*),uae(mx,*),ip(5)
      
      open(39)
      do l1=1,lmax+1
        call labels(l1,l1-1,1,sn,sl,sm)
        ip(l1)=1
        if(uae(ninu(l1)-4,l1)*ups(ninu(l1)-4,l1) .le. 0.0) ip(l1)=-1
        write(39,'(2a3,3i2,1x,f6.3)') '##',sn//sl,l1,l1,l1-1,rc(l1)
        do i=1,ninu(l1),max(ninu(lmax+1)/100,1)
          write(39,'(e12.6,1x,e12.6)') r(i),ups(i,l1)
          if(ups(i+1,l1) .eq. 0.0) goto 20
        enddo
  20    call labels(np(l1),l1-1,1,sn,sl,sm)
        write(39,'(2a3,3i2,1x,a2)') '# ',sn//sl,l1,np(l1),l1-1,'ae'
        do i=1,ninu(l1),max(ninu(lmax+1)/100,1)
          write(39,'(e12.6,1x,e12.6)') r(i),ip(l1)*uae(i,l1)
          if(uae(i+1,l1) .eq. 0.0) goto 21
        enddo
  21    continue
      enddo

      return
      end
c 
c sic double counting corrections
      subroutine outsic(istart,mx,mmax,llmx,r,ups,vps)
      implicit none
      integer  istart,mx,mmax,llmx,iu,iv,i,ll
      real*8   amesh,r(mx),ups(mx,*),vps(mx,*)

      amesh=r(2)/r(1)
      do ll=1,llmx
        iu=istart+ll-1
        iv=istart+5+ll-1
        write(iu,'(a,1x,i5,1x,e20.14,1x,i1,1x,f8.4,1x,e20.14)')
     1    '#',mmax,amesh,ll-1
        do i=1,mmax
c ionic
          write(iu,700) i,r(i),ups(i,ll),vps(i,ll)
        enddo
        close(iu)
      enddo
  700 format(i4,1x,e20.14,2(1x,e20.14))
  710 format(e20.14,1x,e20.14)

      return
      end
