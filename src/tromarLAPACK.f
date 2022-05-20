c $Header:$
c***********************************************************************
c norm-conserving pseudopotentials of Troullier-Martins type
c cite: [1]  N. Troullier, J.L. Martins, Phys. Rev. B 43, 1993 (1991).
c
c input
c l ............ angular momentum channel
c rc ........... cutoff radius
c ein .......... reference energy
c mmax ......... maxiumum used index on radial grid
c r() .......... lograithmic radial grid
c u() .......... radial all-electron wavefunction, changed upon output 
c up() ......... dto. 1st derivative (on the transformed grid!)
c upp() ........ dto. 2nd derivative
c v() .......... all-electron potential
c
c output
c nrc .......... grid index next to cutoff radius
c u() .......... radial pseudowavefunction up to nrc, beyond ae wfct.
c up() ......... dto. 1st derivative (on the transformed grid!)
c upp() ........ dto. 2nd derivative
c vps() ........ pseudopotential up to nrc, beyond ae-potential
c
c Martin Fuchs, FHI der MPG, Berlin, 08-1996
c***********************************************************************
      subroutine tromar(l,rc,nrc,ein,mmax,r,u,up,upp,v,vps)

      implicit  none
      include   'parameter.h'

      logical   tmod
      integer   nrc,mmax,l,ll,i,j,jmx,k,kmx,m,n
      integer   ipvt(5)
c maximum number of iterations for nonlinear eqns 
c jmx .... zero curvature at orgin 
c kmx .... norm conservation
      parameter (jmx=20,kmx=20)

      real*8    rc,ein,al,truecore,pseucore
      real*8    vir0,vir1,vir2,eq0,eq1,eq2,eq3,eq4,eq5,eq6,eq6mod
      real*8    rc1,rc2,rc3,rc4,rc5,rc6,rc7,rc8,rc9,rc10,rc11,rc12
      real*8    cc0,cc1,cc2,cc3,cc4,cc5,cc6,cc0inc,cc2inc
      real*8    polynom,polynom1,polynom2,xxx,eq0old,eq6old
      real*8    r(mx),u(mx),up(mx),upp(mx),v(mx),vps(mx),rsq(mx),rll(mx)
      real*8    ups(mx),ac(5,5),acsav(5,5),bc(5)

      tmod=.false.
c here we can set the l`s for which the flatness condition in the
c troullier-martins construction should be relaxed
c     if(l .eq. 1) tmod=.true.
      eq6mod=0.d0

      ll=l+1
      al=0.1d0*log(r(11)/r(1))
      nrc=log(rc/r(1))/al+1

c we want a positive wavefunction at infinity
      if(u(nrc) .lt. 0.) then
        do i=1,mmax
          u(i)=-u(i)
          up(i)=-up(i)
          upp(i)=-upp(i)
        enddo
      endif

      do i=1,mmax
        rsq(i)=r(i)*r(i)
        rll(i)=r(i)**ll
      enddo

      cc0=0.d0
      cc2=0.d0

c charge inside, nothing fancy just the trapezoidal rule
      truecore=0.5*(u(1)*u(1)*r(1)+u(nrc)*u(nrc)*r(nrc))
      do i=2,nrc-1
        truecore=truecore+u(i)*u(i)*r(i)
      enddo
      truecore=al*truecore+0.5*u(1)*u(1)*r(1)

      rc1 =r(nrc)
      rc2 =rc1*rc1
      rc3 =rc1*rc2
      rc4 =rc1*rc3
      rc5 =rc1*rc4
      rc6 =rc1*rc5
      rc7 =rc1*rc6
      rc8 =rc1*rc7
      rc9 =rc1*rc8
      rc10=rc1*rc9
      rc11=rc1*rc10
      rc12=rc1*rc11

c derivatives of effective potential
      vir0=v(nrc)

      vir2=(16*(v(nrc-1)+v(nrc+1))-(30*vir0+v(nrc-2)+v(nrc+2)))/12.0d0
      vir1=(v(nrc-2)+8*(v(nrc+1)-v(nrc-1))-v(nrc+2))/12.0d0
        
      vir2=(vir2-al*vir1)/(al*r(nrc))**2
      vir1=vir1/(al*r(nrc))

      eq1=log(u(nrc)/rc1**ll)
      eq2=up(nrc)/(al*rc1)/u(nrc)-ll/rc1
      eq3=2*(vir0-ein-ll*eq2/rc1)-eq2*eq2
      eq4=2*(vir1+ll*(eq2/rc1-eq3)/rc1-eq2*eq3)
      eq5=2*(vir2+ll*(2*(eq3-eq2/rc1)/rc1-eq4)/rc1-eq3*eq3-eq2*eq4)

!d     write(ie,*) '& tromar'
!d     write(ie,*) '& l ein         ',l,ein
!d     write(ie,*) '& rc r(nrc) nrc ',rc,rc1,nrc
!d     write(ie,*) '& truecore      ',truecore
!d     write(ie,*) '& vir0 vir1 vir2',vir0,vir1,vir2
!d     write(ie,*) '& eq1 eq2 eq3   ',eq1,eq2,eq3
!d     write(ie,*) '& eq4 eq5       ',eq4,eq5

      ac(1,1)=rc2
      ac(1,2)=rc6
      ac(1,3)=rc8
      ac(1,4)=rc10
      ac(1,5)=rc12

      ac(2,1)=2*rc1
      ac(2,2)=6*rc5
      ac(2,3)=8*rc7
      ac(2,4)=10*rc9
      ac(2,5)=12*rc11

      ac(3,1)=2.d0
      ac(3,2)=30*rc4
      ac(3,3)=56*rc6
      ac(3,4)=90*rc8
      ac(3,5)=132*rc10
 
      ac(4,1)=0.d0
      ac(4,2)=120*rc3
      ac(4,3)=336*rc5
      ac(4,4)=720*rc7
      ac(4,5)=1320*rc9

      ac(5,1)=0.d0
      ac(5,2)=360*rc2
      ac(5,3)=1680*rc4
      ac(5,4)=5040*rc6
      ac(5,5)=11880*rc8

      do n=1,5
        do m=1,5
          acsav(m,n)=ac(m,n)
        enddo
      enddo

c zero curvature at origin condition
 
      cc2inc=0.25d0
      eq6old=0.d0
      do j=1,jmx

        bc(1)=eq1-rc4*cc2
        bc(2)=eq2-4*rc3*cc2
        bc(3)=eq3-12*rc2*cc2
        bc(4)=eq4-24*rc1*cc2
        bc(5)=eq5-24*cc2
  
c numerical recipies
c       call ludcmp(ac,5,5,ipvt,xxx)
c       call lubksb(ac,5,5,ipvt,bc)

c essl
c       call dgef(ac,5,5,ipvt)
c       call dges(ac,5,5,ipvt,bc,0)

c lapack
        call dgesv(5,1,ac,5,ipvt,bc,5,i)
        if(i .ne. 0) then
          if(i .lt. 0) write(ie,*) '& dnlcc7 - stop: singular matrix'
          if(i .gt. 0) write(ie,*) '& dnlcc7 - stop: bad input'
          stop
        endif
  
        cc1=bc(1)
        cc3=bc(2)
        cc4=bc(3)
        cc5=bc(4)
        cc6=bc(5)

c norm conservation: find cc0 by newton method
  
        cc0inc=0.25d0
        eq0old=0.d0
        do k=1,kmx
  
          pseucore=0.d0
          do i=1,nrc
            polynom=rsq(i)*(cc1+rsq(i)*
     1        (cc2+rsq(i)*(cc3+rsq(i)*(cc4+rsq(i)*(cc5+rsq(i)*cc6))))
     1                     )
            ups(i)=rll(i)*exp(polynom)
            pseucore=pseucore+r(i)*ups(i)*ups(i)
          enddo
          pseucore=al*(pseucore
     1      -0.5*(ups(nrc)*ups(nrc)*r(nrc)+ups(1)*ups(1)*r(1))
     1                )+0.5*ups(1)*ups(1)*r(1)
          eq0=log(truecore/pseucore)-2*cc0

          if(abs(eq0) .lt. 1.d-12) goto 18

          cc0inc=-eq0/((eq0-eq0old)/cc0inc)
          eq0old=eq0
          cc0=cc0+cc0inc

          do n=1,5
            do m=1,5
              ac(m,n)=acsav(m,n)
            enddo
          enddo
         
          bc(1)=eq1-rc4*cc2-cc0
          bc(2)=eq2-4*rc3*cc2
          bc(3)=eq3-12*rc2*cc2
          bc(4)=eq4-24*rc1*cc2
          bc(5)=eq5-24*cc2

c numerical recipies
c         call ludcmp(ac,5,5,ipvt,xxx)
c         call lubksb(ac,5,5,ipvt,bc)

c essl
c         call dgef(ac,5,5,ipvt)
c         call dges(ac,5,5,ipvt,bc,0)

c lapack
          call dgesv(5,1,ac,5,ipvt,bc,5,i)
          if(i .ne. 0) then
            if(i .lt. 0) write(ie,*)'& dnlcc7 - stop: singular matrix 1'
            if(i .gt. 0) write(ie,*)'& dnlcc7 - stop: bad input 1'
            stop
          endif

          cc1=bc(1)
          cc3=bc(2)
          cc4=bc(3)
          cc5=bc(4)
          cc6=bc(5)
 
        enddo
   18   if(k .gt. kmx) then
          write(ie,*) '& tromar - stop: convergence error cc0 (norm)'
          write(ie,*) '&  cc0 cc0inc eq0      ',cc0,cc0inc,eq0
          write(ie,*) '&  cc2 cc2inc eq6      ',cc2,cc2inc,eq6
          write(ie,*) '&  iteration j k       ',j,k
          stop
        endif

!d       write(ie,*) '&  iteration j k',j,k
!d       write(ie,*) '&  cc0 cc4    ',cc0,cc4
!d       write(ie,*) '&  cc1 cc5    ',cc1,cc5
!d       write(ie,*) '&  cc2 cc6    ',cc2,cc6
!d       write(ie,*) '&  cc3        ',cc3

        if(tmod) then
         tmod=.false.
         write(6,*) 'Enter eqmod:'
         read(5,*) eq6mod
        endif
        eq6=cc2*(2*l+5)+cc1*cc1+eq6mod

        if(abs(eq6) .lt. 1.d-6) goto 16

        cc2inc=-eq6/((eq6-eq6old)/cc2inc)
        eq6old=eq6
        cc2=cc2+cc2inc

      enddo 
   16 if(j .gt. jmx) then
        write(ie,*) '& tromar - stop: convergence error cc2 (curvature)'
        write(ie,*) '&  cc0 cc0inc eq0      ',cc0,cc0inc,eq0
        write(ie,*) '&  cc2 cc2inc eq6      ',cc2,cc2inc,eq6
        write(ie,*) '&  iteration j k       ',j,k
        stop
      endif
c do we have a solution?
!d     write(ie,*) '&  cc1**2 + cc2*(2l+5)', cc1**2+(2*ll+3)*cc2

c screened pseudopotential 
      bc(1)=u(nrc)
      bc(2)=up(nrc)
      bc(3)=upp(nrc)
      do i=1,nrc
        polynom=cc0+rsq(i)*(cc1+rsq(i)
     1    *(cc2+rsq(i)*(cc3+rsq(i)
     1    *(cc4+rsq(i)*(cc5+rsq(i)*cc6)))))
        polynom1=r(i)*(2*cc1+rsq(i)*(4*cc2+rsq(i)
     1    *(6*cc3+rsq(i)*(8*cc4+rsq(i)*(10*cc5+rsq(i)*12*cc6)))))
        polynom2=2*cc1+rsq(i)*(12*cc2+rsq(i)*(30*cc3+rsq(i)
     1    *(56*cc4+rsq(i)*(90*cc5+rsq(i)*132*cc6))))
        vps(i)=ein+ll*polynom1/r(i)+0.5*(polynom2+polynom1*polynom1)
        u(i)=ups(i)*exp(cc0)
        up(i)=polynom1+ll/r(i)
        upp(i)=((polynom2-ll/r(i)**2)+up(i)*up(i))*u(i)
        up(i)=up(i)*u(i)
c transform to logarithmic grid
        upp(i)=al*r(i)*(al*up(i)+al*r(i)*upp(i))
        up(i)=al*r(i)*up(i)
      enddo
      do i=nrc+1,mmax
        vps(i)=v(i)
      enddo

!d     write(ie,*) '& done with l                 ',l
!d     write(ie,*) '& input  u, up, upp at nrc    ',bc(1),bc(2),bc(3)
!d     write(ie,*) '& pseudo u, up, upp at nrc    ',
!d    1  u(nrc),up(nrc),upp(nrc)
!d     write(ie,*) '& truecore v(nrc)   v(nrc+1)  ',
!d    1  log(truecore/pseucore)-2*cc0,vir0,v(nrc+1)
!d     write(ie,*) '& pseucore vps(nrc) vps(nrc+1)',
!d    1  pseucore,vps(nrc),vps(nrc+1)

      return
      end
c
