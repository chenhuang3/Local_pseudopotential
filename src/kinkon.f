c $Header:$
c fourier-bessel transform for radial wavefunctions
c   and convergence test for orbital kinetic energy in reciprocal
c   space, cf. e.g. [G. Kresse, J. Hafner, JPCM 6, 8245 (1994)]
c integrals for transforms are evaluated by Simpson's rule using
c   a linear mesh in reciprocal space
c wavefunctions are recalculated from effective potential assuming
c   a non-relativistic Schroedinger eq.
c
c input
c io ........... I/O unit
c mmax ......... maximum used index for radial grid
c ekmx ......... energy (Ry) corresponding to maximum wavevector
c g_del ........ mesh increment in reciprocal space
c n ............ principal quantum number
c l ............ angular momentum quantum number
c e ............ eigenvalue
c vps() ........ effective potential
c
c output to unit iu 
c
      subroutine kinkon(iu,mmax,ekmx,g_del,n,l,e,r,vps)

      implicit none
      include  'parameter.h'

      integer  iu,mmax,n,l,i,kmx,j,mch,ifile
      integer  mxk
      parameter(mxk=2000)
      
      real*8   e,et,ry2,ekmx,g_del,x_norm,x_ekin,x_epot,g_ekin,g_norm
      real*8   delk,al,fmom,dmelm,uld,cnorm,c_conv,c_convc,c_convcend
      real*8   r(mx),vps(mx),fr(mx),rho(mx),fk(mxk),gk(mxk),wk1(mxk)
      real*8   wa(mx),wb(mx),u(mx),up(mx),upp(mx)

      parameter(ry2=27.2116)
      external fmom,dmelm
      
      al=log(r(2)/r(1))
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
      enddo

c transform 
      if(e .ne. 0.) then
        et=e
        call dftseq(2,0.d0,mmax,r,n,l,1.d0,vps,wa,wb,
     1        i,mch,uld,et,u,up,upp)
        do i=1,mmax
          fr(i)=u(i)/r(i)
          rho(i)=fr(i)**2
        enddo
        call dfbt(mmax,r,fr,l,ekmx,g_del,mxk,kmx,fk,gk,wk1,0,cnorm)
      endif

cc fourier coefficients
c     do i=1,kmx
c       write(81,*) fk(i),gk(i)
c     enddo
c     close(81)
         
c real space expectation values
      x_norm=fmom(0,mmax,al,1.d0,r,rho)
      x_epot=dmelm(mmax,al,r,vps,rho)
      x_ekin=et-x_epot

      ifile=80+l
c kinetic energy convergence
      delk=g_del/(2*cnorm)
      g_norm=0.d0
      g_ekin=0.d0
c   convergence threshold 3.675e-5 a.u. = 1 meV 
      c_convc=3.675e-5
      c_conv=c_convc*1000
      c_convcend=c_convc
      do i=2,kmx
        g_ekin=g_ekin*g_norm+0.5*delk*(
     1     (fk(i)*gk(i))**2+(fk(i-1)*gk(i-1))**2
     1    )
        g_norm=g_norm+delk*(
     1     gk(i)*gk(i)+gk(i-1)*gk(i-1)
     1    )
        g_ekin=g_ekin/g_norm
        write(ifile,'(f8.2,1x,e12.5,1x,e12.5)') 
     1    fk(i)**2,(g_ekin-x_ekin)*ry2,gk(i)
        if(abs(g_ekin-x_ekin) .le. c_conv) then
          write(iu,601) l,c_conv*ry2,nint(fk(i)**2),g_norm,g_ekin,
     1     nint(fk(i)**2*ry2/2.d0)
          c_conv=c_conv/10.d0
        endif
        if(c_conv .lt. c_convcend) goto 10
      enddo
      write(iu,603) l,c_conv*ry2,nint(fk(i)**2),'<',g_norm,g_ekin
      write(ie,*) 
     1 '& kinkon - bracket not met for cutoff < ekmx =',
     1 nint(ekmx),' Ry'

   10 continue

      write(iu,605) l,n,x_norm,x_ekin

  601 format('ck  ',i2,4x,1p,e8.1,0p,3x,i4,2x,f9.6,2x,e13.6,2x,i5)
  603 format('ck  ',i2,4x,1p,e8.1,0p,3x,i4,1x,1a,f9.6,2x,e13.6)
  605 format('cx  ',i2,1x,i2,18x,f9.6,2x,e13.6/)

      return
      end
c
