!
! For given parameters, construct vlocal and 
! do non-self-consistent calculations on pseudoatom. 
! Compare the eigenvalues with the benchmark.
! 
! inputs:
!   n: main quantum number 
!   l: angular moment quantum number 
!   f: occupation number
!   m: spin_up - spin_down in each orbital?
!   r: radial mesh 
!   vks: KS potential
! 
! output: fcost
!
      subroutine cost(tkli,igr,iexcnow,norb,n,l,f,m,e,mmax,z, 
     &  r,vbare_v,vks,vs,vsp,vspp,lp,nv,nc,
     &  rho_ref,rhop_ref,rhopp_ref,norm0,wfopt,icut_end,
     &  icut,eref,ups_ref,npdim,npar,param,iorb_start,scf_mode,zval,
     &  mag_fit,mag_weight,mag_energy,occ_sp,mag_de, 
     &  vtype,iprint,pcoeff,pen,fcost,wf_match,vloc)

!        implicit real*8 (a-h,o-z)
        implicit none 
        include 'parameter.h'
        include 'default.h'
        
        integer n(ms),l(ms),m(ms),lp(ms),nv,nc
        integer i,j,it,igr,iexcnow,norb,mmax,iprint,npar
        integer icut,npdim,vtype,iorb_start,itest
        integer wfopt_type
        real*8 wfopt(ms),param(npdim)  ! parameters for vlocal 
        real*8 f(ms),e(ms),eref(ms),vks(mx),ups_ref(mx,ms)
        real*8 r(mx),vhxc(mx,ms),vloc0(mx,ms),fac
        real*8 fcost,e_match,wf_match
        real*8 norm0(ms),norm0_sm(ms),norm1(ms),pen,wfdiff(ms)

        ! for vlocal 
        real*8 vbare_v(mx) ! valence-unscreened AE potential 
        real*8 vs,vsp,vspp ! pot, dpot, ddpot at icut 
        real*8 rcut,lam,shift,zval,vks_curv,vks_slop,vks_final(mx)
       
        ! output variables for psatom()
        logical tkli,scf_mode,do_reg
        integer icut_end(ms),icut2,icut3,nin,ndim,k
        real*8 vloc(mx),zncpp
        real*8 de,rhos(mx,ms),rhosp(mx,ms),rhospp(mx,ms) 
        real*8 rho_ref(mx),rhop_ref(mx),rhopp_ref(mx)
        real*8 dc(mx),dcp(mx),dcpp(mx),ups1(mx),ups(mx,ms)
        real*8 svm(ms),svm_roks(ms),svkli(ms)
        real*8 c0,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c(100)
        real*8 xx(npar+1),yy(npar+1),yp(npar+1)
        real*8 sol(3),mat(3,3),invmat(3,3),b(3),amesh,al,t
        real*8 sol4(4),sol5(5),b5(5),bvec(4),mat5(5,5),
     &         mat4(4,4),invmat4(4,4),invmat5(5,5),s1,s2,s3,s4,s5
        real*8 h0,h1,h00,h10,h01,h11,pcoeff,f0,f1,f2,rat,lap,z
        real*8 cheby,d_cheby,d2_cheby,dtdr,norm_rho,norm_rho0,der
        real*8 satf,tolnorm,lder1,lder2,dr,dd,dd_ref,weight,p,p_ref
        real*8 norm_tmp,rcut_low,rcut_high,r_end(norb)
        real*8 legendre,legendre_p,legendre_pp

        integer ispin,ispinmx,iexc_temp,norb2,istart,iend,norbstart
        integer n_sp(ms),l_sp(ms),m_sp(ms) ! for spin polarized 
        real*8  e_sp(ms),occ(ms),vlocs(mx,ms),vtmp(mx,ms)
        
        ! fitting magnetic energy 
        integer mag_fit
        real*8  mag_weight,mag_energy,occ_sp(ms,2),mag_de
        real*8  ups_tmp(mx,ms),etot1,etot2

        rcut=r(icut)
        amesh=r(2)/r(1)
        al=log(amesh)
        pen=0.d0 


        !===============================
        ! make vlocal 
        !===============================
        select case(vtype)
        case (1)
          ! cutoff function 
          lam   = param(1) 
          shift = param(2) 
          rcut  = param(3)
          do i=1,mmax
            fac=(1.d0-exp(-lam*r(i)))/(1.d0+exp(-lam*(r(i)-rcut)))
            do j=1,ms
              vloc0(i,j) = fac*(vks(i)+shift)-shift
            enddo
          enddo
          vloc=vloc0(:,1)

        case (2) ! polynomial function 
          c(1:npar)=param(1:npar)  ! first is c0, then c2,c3,c4,...

          ! compute c4, c5, and c6 based on vs,vsp,vspp
          call get_mat_npar(mat,rcut,npar)
          call inverse(mat,invmat,3)

          b(1)=vs-c(1)
          do i=2,npar
            b(1)=b(1)-c(i)*rcut**i
          enddo

          b(2)=vsp
          do i=2,npar
            b(2)=b(2)-i*c(i)*rcut**(i-1)
          enddo

          b(3)=vspp
          do i=2,npar
            b(3)=b(3)-i*(i-1)*c(i)*rcut**(i-2)
          enddo 

          sol = matmul(invmat,b)
          c(npar+1) = sol(1)
          c(npar+2) = sol(2)
          c(npar+3) = sol(3)
          param(npar+1:npar+3)=c(npar+1:npar+3)

          ! make vlocal
          do i=1,mmax
            do j=1,ms
              if (i<=icut) then 
                if (npar==3) then 
                  vloc0(i,j)=c0+c2*r(i)**2+c3*r(i)**3
     &              +c4*r(i)**4+c5*r(i)**5+c6*r(i)**6
                else 
                  vloc0(i,j)=c(1)
                  do k=2,npar+3
                    vloc0(i,j)=vloc0(i,j)+c(k)*r(i)**k
                  enddo
                endif 
              else 
                ! r>rcut 
                if (scf_mode) then 
                  vloc0(i,j)=vbare_v(i)
                else 
                  vloc0(i,j)=vks(i)
                endif 
              endif 
            enddo
          enddo 
          vloc=vloc0(:,1)

        case (3)  ! polynomial (consider linear term)
          if (.not. scf_mode) then
            print *,'only work with scf_mode=true case'
            stop
          endif 
          c(1:npar)=param(1:npar) 

          ! compute high-order c based on vs,vsp,vspp
          call get_mat_npar_type3(mat,rcut,npar)
          call inverse(mat,invmat,3)

          b(1)=vs
          do i=1,npar
            b(1)=b(1)-c(i)*rcut**(i-1)
          enddo

          b(2)=vsp
          do i=2,npar
            b(2)=b(2)-(i-1)*c(i)*rcut**(i-2)
          enddo

          b(3)=vspp
          do i=3,npar
            b(3)=b(3)-(i-1)*(i-2)*c(i)*rcut**(i-3)
          enddo 

          sol = matmul(invmat,b)
          c(npar+1) = sol(1)
          c(npar+2) = sol(2)
          c(npar+3) = sol(3)
          param(npar+1:npar+3)=c(npar+1:npar+3)

          ! make vlocal
          do i=1,mmax
            do j=1,ms
              if (i<=icut) then 
                vloc0(i,j)=0.d0
                do k=1,npar+3
                  vloc0(i,j)=vloc0(i,j)+c(k)*r(i)**(k-1)
                enddo
              else 
                ! r>rcut 
                if (scf_mode) then 
                  vloc0(i,j)=vbare_v(i)
                else 
                  vloc0(i,j)=vks(i)
                endif 
              endif 
            enddo
          enddo 
          vloc=vloc0(:,1)

        case (4) 
          ! vlocal is on mesh 
          do i=1,mmax
            do j=1,ms
              if (i<=icut) then 
                vloc0(i,j)=param(i)
              else 
                ! r>rcut 
                if (scf_mode) then 
                  vloc0(i,j)=vbare_v(i)
                else 
                  vloc0(i,j)=vks(i)
                endif 
              endif 
            enddo
          enddo 
          vloc=vloc0(:,1)

        case (5) ! polynomial function (no linear term, x=r/rc)
          c(1:npar)=param(1:npar)  ! first is c0, then c2,c3,c4,... do not have c1
          dtdr=1.d0/rcut 

          ! compute c4, c5, and c6 based on vs,vsp,vspp
          call get_mat_npar_t(mat,rcut,npar) ! based on t=r/rc, for r=rc, t=1
          call inverse(mat,invmat,3)

          b(1)=vs-c(1)
          do i=2,npar
            b(1)=b(1)-c(i)
          enddo

          b(2)=vsp
          do i=2,npar
            b(2)=b(2)-i*c(i)*dtdr
          enddo

          b(3)=vspp
          do i=2,npar
            b(3)=b(3)-i*(i-1)*c(i)*dtdr**2
          enddo 

          sol = matmul(invmat,b)
          c(npar+1) = sol(1)
          c(npar+2) = sol(2)
          c(npar+3) = sol(3)
          param(npar+1:npar+3)=c(npar+1:npar+3)
          
          if (iprint==1) then 
            write(6,'(a,i2,a,es20.12)')'c(',npar+1,'):',c(npar+1)
            write(6,'(a,i2,a,es20.12)')'c(',npar+2,'):',c(npar+2)
            write(6,'(a,i2,a,es20.12)')'c(',npar+3,'):',c(npar+3)
          endif 

          ! make vlocal
          do i=1,mmax
            t=r(i)/rcut
            do j=1,ms
              if (i<=icut) then 
                vloc0(i,j)=c(1)
                do k=2,npar+3
                  vloc0(i,j)=vloc0(i,j)+c(k)*t**k
                enddo
              else 
                ! r > rcut 
                if (scf_mode) then 
                  vloc0(i,j)=vbare_v(i)
                else 
                  vloc0(i,j)=vks(i)
                endif 
              endif 
            enddo
          enddo 
          vloc=vloc0(:,1)

        case (6)
          !----------------------------------
          ! legendre polynomial function 
          !----------------------------------
          dtdr=2.d0/rcut  ! t = 2*r/rcut - 1
          c(1:npar)=param(1:npar)  

          ! calculate additional 4 coefficients c(npar+1:npar+4)
          ! based on the continuity of v,vp,vpp at rcut and vpp=0 at r=0
          do i=1,4
            mat4(1,i)=legendre(npar+i-1,1.d0)            ! continuity of v at r=rcut
            mat4(2,i)=legendre_p(npar+i-1,1.d0)*dtdr     ! continuity of vp at r=rcut
            mat4(3,i)=legendre_pp(npar+i-1,1.d0)*dtdr**2 ! continuity of vpp at r=rcut
            mat4(4,i)=legendre_p(npar+i-1,-1.d0)*dtdr    ! vp=0 at r=0
          enddo 

          s1=0.d0;s2=0.d0;s3=0.d0;s4=0.d0
          do i=1,npar
            s1=s1+c(i)*legendre(i-1,1.d0)
            s2=s2+c(i)*legendre_p(i-1,1.d0)*dtdr
            s3=s3+c(i)*legendre_pp(i-1,1.d0)*dtdr**2
            s4=s4+c(i)*legendre_p(i-1,-1.d0)*dtdr
          enddo
          bvec(1) = vs-s1
          bvec(2) = vsp-s2
          bvec(3) = vspp-s3
          bvec(4) = -s4

          call inverse(mat4,invmat4,4)
          sol4 = matmul(invmat4,bvec)
          c(npar+1:npar+4)=sol4
          param(npar+1:npar+4)=c(npar+1:npar+4)

          ! make vlocal
          vloc0=0.d0 
          do i=1,mmax
            t=2.d0*r(i)/rcut-1.d0
            do j=1,ms
              if (i<=icut) then 
                do k=1,npar+4
                  vloc0(i,j)=vloc0(i,j)+c(k)*legendre(k-1,t)
                enddo
              else 
                ! r > rcut 
                if (scf_mode) then 
                  vloc0(i,j)=vbare_v(i)
                else 
                  vloc0(i,j)=vks(i)
                endif 
              endif 
            enddo
          enddo 
          vloc=vloc0(:,1)

        case (7)
          !-------------------------------------------
          ! legendre polynomial function 
          ! condition: v'(0) and v''(0) are zero
          !-------------------------------------------
          dtdr=2.d0/rcut  ! t = 2*r/rcut - 1
          c(1:npar)=param(1:npar)  

          ! calculate additional 4 coefficients c(npar+1:npar+4)
          ! based on the continuity of v,vp,vpp at rcut and vpp=0 at r=0
          do i=1,5
            mat5(1,i)=legendre(npar+i-1,1.d0)             ! continuity of v at r=rcut
            mat5(2,i)=legendre_p(npar+i-1,1.d0)*dtdr      ! continuity of vp at r=rcut
            mat5(3,i)=legendre_pp(npar+i-1,1.d0)*dtdr**2  ! continuity of vpp at r=rcut
            mat5(4,i)=legendre_p(npar+i-1,-1.d0)*dtdr     ! v'=0 at r=0
            mat5(5,i)=legendre_pp(npar+i-1,-1.d0)*dtdr**2 ! v''=0 at r=0
          enddo 

          s1=0.d0;s2=0.d0;s3=0.d0;s4=0.d0;s5=0.d0
          do i=1,npar
            s1=s1+c(i)*legendre(i-1,1.d0)
            s2=s2+c(i)*legendre_p(i-1,1.d0)*dtdr
            s3=s3+c(i)*legendre_pp(i-1,1.d0)*dtdr**2
            s4=s4+c(i)*legendre_p(i-1,-1.d0)*dtdr     ! v'=0 at r=0
            s5=s5+c(i)*legendre_pp(i-1,-1.d0)*dtdr**2 ! v''=0 at r=0
          enddo
          b5(1) = vs-s1
          b5(2) = vsp-s2
          b5(3) = vspp-s3
          b5(4) = -s4
          b5(5) = -s5

          call inverse(mat5,invmat5,5)
          sol5 = matmul(invmat5,b5)
          c(npar+1:npar+5)=sol5
          param(npar+1:npar+5)=c(npar+1:npar+5)

          ! make vlocal
          vloc0=0.d0 
          do i=1,mmax
            t=2.d0*r(i)/rcut-1.d0
            do j=1,ms
              if (i<=icut) then 
                do k=1,npar+5
                  vloc0(i,j)=vloc0(i,j)+c(k)*legendre(k-1,t)
                enddo
              else 
                ! r > rcut 
                if (scf_mode) then 
                  vloc0(i,j)=vbare_v(i)
                else 
                  vloc0(i,j)=vks(i)
                endif 
              endif 
            enddo
          enddo 
          vloc=vloc0(:,1)
        end select 


        !=======================================
        ! self-consistent atom calculation
        !=======================================
        vhxc=0.d0 ! set initial xc+hartree to zero
        if (scf_mode) it=itmx 
        if (.not.scf_mode) it=1
        zncpp = 0.d0 
!        call psatom(it,igr,iexcnow,zncpp,0,norb,
!     &   n,l,f,m,e,nin,de,mmax,r,vhxc,vloc0,rhos,rhosp,rhospp, 
!     &   dc,dcp,dcpp,ups,tkli,svm,svm_roks)
        call psatom_etot(it,igr,iexcnow,zncpp,0,norb,
     &   n,l,f,m,e,nin,de,mmax,r,vhxc,vloc0,rhos,rhosp,rhospp, 
     &   dc,dcp,dcpp,ups,tkli,svm,svm_roks,etot1)
        
        if (iprint==1) then 
          do i=1,norb
            write(6,'(a,2i4,a,f6.2,a,f12.5,a,f12.5)') 
     &      'n,l: ',n(i),l(i),' occ:',f(i),' eig:',e(i),' eref:',eref(i)
          enddo
        endif 


        !====================================
        ! match orbitals in [rcut, inf]
        !====================================
        wfdiff=0.d0 
        do i=1,norb
          ndim = icut_end(i)
          ! copy ups to ups1, and flip it if needed 
          ups1=ups(:,i) 
          if (ups(icut,i)*ups_ref(icut,i)<0.d0) ups1=-ups(:,i)
          do j=icut,ndim
            wfdiff(i) = wfdiff(i) + al*r(j)*(ups1(j)-ups_ref(j,i))**2
          enddo
          wfdiff(i) = wfdiff(i)
     &      + al*(23.0d0*r(ndim-2)*(ups1(ndim-2)-ups_ref(ndim-2,i))**2
     &      + 28.0d0*r(ndim-1)*(ups1(ndim-1)-ups_ref(ndim-1,i))**2
     &      +  9.0d0*r(ndim)*(ups1(ndim)-ups_ref(ndim,i))**2)/24.0d0
        enddo


        !=========================================
        ! compute cost function
        !=========================================
        e_match=0.d0 
        do i=iorb_start,norb
          e_match=e_match+(e(i)-eref(i))**2
        enddo 

        ! optimize wave function
        wf_match = 0.d0 
        wfopt_type = 1  ! 1: match norm inside rcut 
                        ! 2: match u(rcut) and u'(rcut)
        do i=1,norb 
          wf_match = wf_match + wfopt(i)*wfdiff(i)
        enddo
        fcost = e_match + wf_match
        

        !========================================
        ! test spin configuration 
        !========================================
        if (mag_fit>0) then 
          ispinmx=2
          ! map for spin polarized calculation
          if(ispinmx.gt.1) then
            if(iexcnow.eq.3 .or. iexcnow.eq.8) then
              iexc_temp=3                !LDA
            else if(iexcnow.eq.5) then
              iexc_temp=5                !BP GGA
            else if(iexcnow.eq.4) then  
              iexc_temp=6                !PW91 GGA
            else if(iexcnow.eq.6) then
              iexc_temp=9                !PBE GGA
            else
            write(ie,'(a/,a,i3,a)')
     1      '%pslp - ERROR: spin polarized calculation not implemented',
     1      ' for input iexc=',iexcnow,'. ACTION: set iexc=3,4,5,6,8.'
            stop
            endif
          endif
          ! extend lp(i) to consider spin polarization
          do i=1,nc+nv
            lp(i+nc+nv)=lp(i)
          enddo 
          !--- set n, l, m ------------------
          ! n(1:norb/2)      -- spin up
          ! n(norb/2+1:norb) -- spin dn
          ! similar for l, m, and occ
          norb2=0 
          do ispin=1,ispinmx
            norbstart=norb2+1
            istart=nc*ispin+nv*(ispin-1)+1 
            iend=istart+nv-1
            do i=istart,iend
              norb2=norb2+1        ! only count valence orbitals
              n_sp(norb2)=lp(i)+1  ! for pseudo-atom n is based on l 
              do j=norbstart,norb2-1
                if(l_sp(j) .eq. lp(i)) n_sp(norb2)=n_sp(norb2)+1
              enddo
              l_sp(norb2)=lp(i)           ! l is from input file 
              occ(norb2)=occ_sp(i-istart+1,ispin)  ! occ goes with l
              m_sp(norb2)=ispin  ! spin index for valence orbitals
            enddo
          enddo
          do j=1,ms 
            vlocs(:,j)=vloc
          enddo
          e_sp=0.d0; vtmp=0.d0; zncpp=0.d0; it=itmx
          iexc_temp = -abs(iexc_temp) ! activate spin polarized calculation
          call psatom_etot(it,igr,iexc_temp,zncpp,0,norb2,
     &     n_sp,l_sp,occ,m_sp,
     &     e_sp,nin,de,mmax,r,vtmp,vlocs,rhos,rhosp,rhospp,
     &     dc,dcp,dcpp,ups_tmp,tkli,svm,svm_roks,etot2)
          !print *,etot2
          !print *,'e(up):',e_sp(1:norb2/2)
          !print *,'e(dn):',e_sp(norb2/2+1:norb2)
          !print *,'etotal diff:',etot2-etot1
          mag_de = etot2 - etot1 
          fcost = fcost + mag_weight*(mag_de-mag_energy)**2
        endif 
      end 


      subroutine get_mat_npar(mat,rcut,npar)
        implicit none 
        integer npar
        real*8 mat(3,3), rcut

        ! compute c4, c5, and c6 based on vs,vsp,vspp
        mat(1,1)=rcut**(npar+1)
        mat(1,2)=rcut**(npar+2)
        mat(1,3)=rcut**(npar+3)

        mat(2,1)=(npar+1)*rcut**npar
        mat(2,2)=(npar+2)*rcut**(npar+1)
        mat(2,3)=(npar+3)*rcut**(npar+2)

        mat(3,1)=(npar+1)*npar*rcut**(npar-1)
        mat(3,2)=(npar+2)*(npar+1)*rcut**(npar)
        mat(3,3)=(npar+3)*(npar+2)*rcut**(npar+1)
      end 

      subroutine get_mat_npar_t(mat,rcut,npar)
        implicit none 
        integer npar
        real*8 mat(3,3),dtdr,rcut
        dtdr = 1.d0/rcut 

        ! compute c4, c5, and c6 based on vs,vsp,vspp
        mat(1,1)=1.d0
        mat(1,2)=1.d0
        mat(1,3)=1.d0

        mat(2,1)=(npar+1)*dtdr
        mat(2,2)=(npar+2)*dtdr
        mat(2,3)=(npar+3)*dtdr

        mat(3,1)=(npar+1)*npar*dtdr**2
        mat(3,2)=(npar+2)*(npar+1)*dtdr**2
        mat(3,3)=(npar+3)*(npar+2)*dtdr**2
      end 
      
      subroutine get_mat_npar_type3(mat,rcut,npar)
        implicit none 
        integer npar
        real*8 mat(3,3), rcut

        ! compute c4, c5, and c6 based on vs,vsp,vspp
        mat(1,1)=rcut**npar
        mat(1,2)=rcut**(npar+1)
        mat(1,3)=rcut**(npar+2)

        mat(2,1)=npar*rcut**(npar-1)
        mat(2,2)=(npar+1)*rcut**npar
        mat(2,3)=(npar+2)*rcut**(npar+1)

        mat(3,1)=npar*(npar-1)*rcut**(npar-2)
        mat(3,2)=(npar+1)*npar*rcut**(npar-1)
        mat(3,3)=(npar+2)*(npar+1)*rcut**npar
      end 

      ! compute gradient (central finite diff)
      subroutine calc_der(mx,y,r,i,g)
        implicit none 
        integer mx,i
        real*8 y(mx),r(mx),g,h0,h1,f0,f1,f2,rat

        ! get gradient at r(icut)
        h0=r(i)-r(i-1)
        h1=r(i+1)-r(i)
        f0=y(i-1)
        f1=y(i)
        f2=y(i+1)
        rat = h1/h0

        ! using f(i-1),f(i),f(i+1)
        g = (f2-rat**2*f0-(1.d0-rat**2)*f1)/(h1*(1.d0+rat)) ! dv/dx
      endsubroutine 
      
