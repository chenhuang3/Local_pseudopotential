      !
      ! write local pseudopotential for ABINIT 
      !
      subroutine outlps(mmax,zatom,zion,norb,npdim,npar,param,wfopt,
     &                  pcoeff,vtype,iorb_start,rcut,icut,icoul,r,v)
 
        implicit none
        include 'parameter.h'
       
        logical use_poly

        integer iorb_start,vtype,i,j,ngrid,icoul
        integer icut,norb,npdim,npar,mmax
        real*8 dx,h,t,slope1,slope2,h00,h01,h10,h11
        real*8 param(npdim),rcut,pcoeff,wfopt(ms)
        real*8 r(mx),v(mx),zatom,zion,rmax
        real*8 rr(mx+1),vv(mx+1),vp(mx+1)
        real*8,allocatable :: x(:),v2(:)
       
        print *,''
        print *,'enter outlps()'

        !!dx = 1e-2
        !!rmax = min(r(mmax),20.d0)
        dx = 1e-3 
        rmax = min(r(mmax),10.d0)
        ngrid = int(rmax/dx)
       
        if (r(icoul)>rmax) then 
          print *,'r(icoultail)>rmax, why Coulomb tail is so far? stop'
          stop
        endif 
        !
        ! Can we directly write vlocal using polynomial?
        !
        if (r(icut)>r(icoul) .and. vtype==5) then 
          use_poly=.true. 
        else 
          use_poly=.false.
        endif 
        print *,'outlps() ==> use_poly: ',use_poly,''
       
        allocate(x(ngrid),v2(ngrid))
        
        if (use_poly) then 
          !
          ! directly use the polynomial (for vtype=5 only)
          !
          print *,'outlps() ==> directly use polynomial '
          do i=1,ngrid
            x(i)=dble(i-1)*dx
            if (x(i)>r(icut)) then 
              v2(i)=-zion/x(i)
            else 
              v2(i)=param(1)
              do j=2,npar+3
                v2(i)=v2(i)+param(j)*(x(i)/r(icut))**j
              enddo 
            endif 
          enddo
        else 
          !
          ! interpolate potential based on pchip
          !
          print *,'outlps() ==> interpolate use pchip '
          ! add the first point (0,v(1)) to rr
          rr(1)=0.d0
          rr(2:mmax+1)=r
          vv(1)=v(1)
          vv(2:mmax+1)=v
       
          ! make derivatives 
          vp(1) = (vv(2)-vv(1))/(rr(2)-rr(1))
          vp(mmax+1) = (vv(mmax+1)-vv(mmax))/(rr(mmax+1)-rr(mmax))
       
          do i=2,mmax
            slope1=(vv(i)-vv(i-1))/(rr(i)-rr(i-1))
            slope2=(vv(i+1)-vv(i))/(rr(i+1)-rr(i))
            if (slope1*slope2>0.d0) then 
              ! Harmonic mean 
              vp(i)=2.d0*slope1*slope2/(slope1+slope2)
            else
              vp(i)=0.d0
            endif 
          enddo
       
          ! interpolate potential based on pchip
          do i=1,ngrid
            x(i)=dble(i-1)*dx
            ! get the segment 
            do j=1,mmax
              if (x(i)>=rr(j) .and. x(i)<=rr(j+1)) then 
                h = rr(j+1)-rr(j)
                t = (x(i)-rr(j))/h
                h00=2.d0*t**3-3.d0*t**2+1.d0 
                h10=t**3-2.d0*t**2+t
                h01=-2.d0*t**3+3.d0*t**2
                h11=t**3-t**2
                v2(i)=h00*vv(j)+h10*h*vp(j)+h01*vv(j+1)+h11*h*vp(j+1)
                exit 
              endif 
            enddo 
          enddo
        endif 
       
        !
        ! write CPI file for abinit 
        ! 
        open(file='vlps.cpi',unit=111,action='write',form='formatted')
        write(111,'(a)')'local psp'
        write(111,'(2f10.3,a)')zatom,zion,'   20210830  zatom,zion,pspd'
        write(111,'(a,i6,a)')'8  11  0  0  ',ngrid, 
     &    ' 0  pspcode,pspxc,lmax,lloc,mmax,r2well'
        write(111,'(a)')'0   -1   0    rchrg,fchrg,qchrg'
        write(111,'(a)')'0    0   0   0   0   nproj'
        write(111,'(a)')'0  extension_switch'
        write(111,'(a)')'0'
       
        do i=1,ngrid
          write(111,'(i5,es14.5,es14.6)')i,x(i),v2(i)
        enddo 

        ! print parameter information 
        write(111,'(a,L2)')'# use_poly: ',use_poly
        write(111,'(a,I2)')'# vtype: ',vtype
        write(111,'(a,I2)')'# iorb_start: ',iorb_start
        write(111,'(a,f12.6)') '# r(icut):',r(icut)
        write(111,'(a,f12.6)') '# rcut:   ',rcut
        write(111,'(a,es13.6)')'# pcoeff: ',pcoeff
        do i=1,norb
          write(111,'(a,es20.12)')'# wfopt: ',wfopt(i)
        enddo
        write(111,'(a,i3)')'# npar:',npar
        do j=1,npar+3
          write(111,'(a,i2,a,es20.12)')'# param(',j,')',param(j)
        enddo
        write(111,'(a)')'----------- vlocal.ini ------------'
        close(111)
        call system('cat vlocal.ini >> vlps.cpi')
       
        deallocate(x,v2)

        print *,'leave outlps()'

      end



