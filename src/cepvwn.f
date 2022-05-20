      subroutine cepvwn(dens,xceng,vxc,nsubpt,npts,npolrz)

      implicit real*8(a-h,o-z)

      dimension dens(npts,3)
      dimension xceng(npts),vxc(npts,2)

      dimension ex(128),ec(128),vx(128,2),vc(128,2)
      dimension grmod(128,3),dlap(128,3),dldlp(128,3)
      dimension rs(128),rx(128),zeta(128)
      dimension dkf(128,2),sperd(128,2),faca(128,2),facb(128,2),
     &          facc(128,2),facd(128,2),fperd(128,2),sdfds(128,2),
     &          dsdfd(128,2),tperd(128,2),uperd(128,2),eperd(128,2)
      dimension ecp(128),decp(128),vcp(128),ecfp(128),decfp(128),
     &          acp(128),dacp(128),vcfp(128)
      dimension d2min1(128),cnperd(128),finc(128),esfi(128),
     &          corrpc(128),denst1(128),densa2(128),densb2(128),
     &          constt(128),work1(128),work2(128),work3(128),
     &          work4(128),work5a(128),work5b(128),work6a(128),
     &          work6b(128),work7a(128),work7b(128)
      
      dimension cap(4),cep(4),cef(4),cjp(4)
      
      data dperdm,dperdn,dperdb,dperdc,alalph,albeta
     $     / 0.06666666666667d0,0.06666666666667d0,
     $     14.d0,0.2d0,0.05d0,5.d0/
      data detol,zetol,ftyl,bkbeta/1.d-20,1.d-10,0.11d0,0.0042d0/
      data cap/-0.0337737,1.13107,13.0045,-0.0047584/
     &    ,cep/ 0.0621814,3.72744,12.9352,-0.1049800/
     &    ,cef/ 0.0310907,7.06042,18.0578,-0.3250000/
     &    ,cjp/0.023266,7.389e-06,8.723,0.472/
      data third,fthrd,fthrd2,fdenom,pi,pi43,pi49
     $     / 0.33333333333333d0,1.33333333333333d0,5.1297633d0,
     $     0.51984204d0,3.14159265358979d0,
     $     4.18879020478638d0,1.63696380104572d0 /
      data twthrd,sevsix,fvthrd,ffifth,svthrd,zero,one,two,three,four,
     $     five,six,dnine,expcut
     $     / 0.66666666666667d0,1.16666666666667d0,1.66666666666667d0,
     $     0.8d0,2.33333333333333d0,
     $     0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,9.d0,80.d0 /
      
c----------------------------------------------------------------------
      
      fspin(z)=(one+z)**fthrd+(one-z)**fthrd
      fspol(z)=(one+z)**third
      fspol1(z)=(fspol(+z)-fspol(-z))*fthrd2/two
      fspol2(z)=(fspin(z)-two)/fdenom
      
      gvwn(x,a,b,c,x0)=a*log(x*x/(x-x0)**2)-a/(x0*x0+b*x0+c)*
     &     ((x0*x0+c)*log((x*x+b*x+c)/(x-x0)**2)+
     &     (x0*x0-c)*atan(sqrt(four*c-b*b)/(two*x+b))*two*b/
     &     sqrt(four*c-b*b))
      dgvwn(x,a,b,c,x0)=two*a*(one/x-(one+b/(x-x0))*x/(x*x+b*x+c))
      
      do ipt=1, nsubpt
        if ( dens(ipt,3).gt.detol ) then
          rs(ipt) = one/(pi43*dens(ipt,3))**third
        else
          rs(ipt) = one/(pi43*detol*detol)**third
        end if
      end do
      
      axperd = -0.75*(three/pi)**third
      
      do ipt=1, nsubpt
        if ( dens(ipt,3).gt.detol ) then
          zeta(ipt) = (dens(ipt,1)-dens(ipt,2))/dens(ipt,3)
          if ( zeta(ipt).gt.+one ) zeta(ipt) = +one
          if ( zeta(ipt).lt.-one ) zeta(ipt) = -one
        else
          zeta(ipt) = zero
        end if
      end do
      
c
c  proceed unto the local correlation functional
c
 3052 continue

c***********************************************************************
c                                                                      *
c  local correlation energy and potentials                             *
c                                                                      *
c***********************************************************************
      do ipt=1, nsubpt
        rx(ipt) = sqrt(rs(ipt))
      end do
      
      do ipt=1, nsubpt
        ecp(ipt) = gvwn(rx(ipt),cep(1),cep(2),cep(3),cep(4))/two
      end do
      
      do ipt=1, nsubpt
        decp(ipt) = dgvwn(rx(ipt),cep(1),cep(2),cep(3),cep(4))/two
      end do
      
      do ipt=1, nsubpt
        vcp(ipt) = ecp(ipt)-decp(ipt)*rx(ipt)/six
      end do
      
      do ipt=1, nsubpt
        ec(ipt) = ecp(ipt)
      end do
      
      do ispin=1, 2
        do ipt=1, nsubpt
          vc(ipt,ispin)=vcp(ipt)
        end do
      end do
      
c     
c  skip the following part if there is no spin polarization
c     
      if ( npolrz.eq.0 ) goto 3053
      
      do ipt=1, nsubpt
        ecfp(ipt) = (gvwn(rx(ipt),cef(1),cef(2),cef(3),cef(4))/two)
     $       - ecp(ipt)
      end do
      
      do ipt=1, nsubpt
        decfp(ipt) = (dgvwn(rx(ipt),cef(1),cef(2),cef(3),cef(4))/two)
     $       - decp(ipt)
      end do
      
      do ipt=1, nsubpt
        acp(ipt) = (gvwn(rx(ipt),cap(1),cap(2),cap(3),cap(4))/two)
     $       * three/fthrd2
      end do
      
      do ipt=1, nsubpt
        dacp(ipt)=(dgvwn(rx(ipt),cap(1),cap(2),cap(3),cap(4))/two)
     $       * three/fthrd2
      end do
      
      do ipt=1, nsubpt
        vcfp(ipt) = fspol2(zeta(ipt))*(ecfp(ipt)*(zeta(ipt)**4)
     &       + acp(ipt)*(one-(zeta(ipt)**4))
     &       - (rx(ipt)/six)*(decfp(ipt)*(zeta(ipt)**4)
     &       + dacp(ipt)*(one-(zeta(ipt)**4))))
      end do
      
      do ipt=1, nsubpt
        ec(ipt) = ec(ipt) + fspol2(zeta(ipt))*(ecfp(ipt)*(zeta(ipt)**4)
     $       + acp(ipt)*(one-zeta(ipt)**4))
      end do
      
      do ipt=1, nsubpt
      vc(ipt,1) = vcp(ipt) + vcfp(ipt) + (one-zeta(ipt))
     $       * ( fspol1(zeta(ipt))
     $       * (ecfp(ipt)*(zeta(ipt)**4) +acp(ipt)*(one-(zeta(ipt)**4)))
     $       + four*fspol2(zeta(ipt))*(ecfp(ipt)*(zeta(ipt)**3)
     $       - acp(ipt)*(zeta(ipt)**3)) )
      end do
      
      do ipt=1, nsubpt
      vc(ipt,2) = vcp(ipt) + vcfp(ipt) - (one+zeta(ipt))
     $       * ( fspol1(zeta(ipt))
     $       * (ecfp(ipt)*(zeta(ipt)**4) +acp(ipt)*(one-(zeta(ipt)**4)))
     $       + four*fspol2(zeta(ipt))*(ecfp(ipt)*(zeta(ipt)**3)
     $       - acp(ipt)*(zeta(ipt)**3)) )
      end do
      
 3053 continue
      
c     
c  proceed unto the end of the subroutine if only a local spin
c  density calculation is being performed
c     
 205  continue
      
      do ipt=1, nsubpt
        xceng(ipt) = ex(ipt) + ec(ipt)
      end do
      
      do ispin=1, 2
        do ipt=1, nsubpt
          vxc(ipt,ispin) = vx(ipt,ispin) + vc(ipt,ispin)
        end do
      end do
      
      return
      end

