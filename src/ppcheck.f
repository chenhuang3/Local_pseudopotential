c $Header:$
c-----------------------------------------------------------------------
c pseudopotential
c analyzing routine for spectrum and logarithmic derivatives
c includes Kleinman-Bylander (KB) separable form 
c spectrum of nonlocal (KB) hamiltonian
c
c input
c iu ........... output unit
c tnrl ......... true  = reference calculation nonrelativistic
c                false = dto scalar-relativistic 
c tkb .......... true  = test separable potentials for ghost states
c                        & evaluate their log derivatives
c tlgd ......... true  = do semilocal/full potential log derivatives
c lmax ......... angular momentum channel present among potentials
c l_loc ........ angular momentum channel for local potential
c lbeg ......... lowest angular momentum channel to check
c lend ......... highest  dto
c mmax ......... maxiumum radial grid index used
c rld .......... diagnostic radius for calculating log derivative
c e() .......... reference spectrum
c r() .......... radial grid
c vae() ........ full potential
c vi() ......... screening potential (from pseudo valence wfct)
c uu() ......... pseudo wavefunctions (projectors for separable case)
c                vkb_l = u_l (v_l - v_loc), psi = u(r)/r * sph. harmonic
c vbare() ...... ionic pseudopotential components (semilocal form)
c 
c output
c Gonze's analysis for ghost states to unit io
c logarithic derivatives to
c   fort.5[l = 0,1,2...]	all-electron case
c   fort.6[l = 0,1,2...]	semilocal case
c   fort.7[l = 0,1,2...]	separable case (sbrt. klbyan)
c
c-----------------------------------------------------------------------
      subroutine ppcheck(iu,tnrl,tlgd,tkb,l_loc,lmax,lbeg,lend,mmax,rld
     1,  e,r,vae,vi,uu,vbare)

      implicit  none
      include   'parameter.h' 
      logical   tprint
      integer   mxe
      parameter (mxe=1000)
      real*8    ry2
      parameter (ry2=27.2116d0)
      logical   tnrl,tkb,tlgd
      character*1 sn,sl,sm
      integer   iu,lbeg,lend,l_loc,lmax,mmax,irl,iep,unit,llbeg,llend
      integer   llmx,i,k,ild,ll,ll_loc,nin,mch,kind,nodes,j
      real*8    ckb,vnorm,ekb,eloc0,eloc1,uld,rld,z
      real*8    es1,es2,eccmx,al,gsum,delta,deltaold,unused,dmelm
      real*8    e(ms),r(mx),vae(mx),vi(mx),uu(mx,ms),vbare(mx,ms)
      real*8    ukb(mx,ms),vsl(mx),u(mx),up(mx),upp(mx),wa(mx),wb(mx)
      real*8    vkb(mx,5),vaverage(5),enl(ms),ecc(mxe),wc(mx)
      real*8    cputime,t1
      external  dmelm

c tprint: print the s wfct. for semilocal/nonlocal potential
      tprint=.false.
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
      enddo
      al=0.1d0*log(r(11)/r(1))
      ll_loc=abs(l_loc)+1
      lbeg=abs(lbeg)
      lend=abs(lend)
      llbeg=lbeg+1
      llend=lend+1
      llmx=lmax+1

c separable potentials
      if(tkb) write(iu,635)

      do ll=llbeg,min(llend,llmx)
        if(tkb .and. ll .ne. ll_loc) then

          do i=mmax,1,-1
            vkb(i,ll)=(vbare(i,ll)-vbare(i,ll_loc))
     1        *uu(i,ll)/r(i)**2
          enddo
          vaverage(ll)=dmelm(mmax,al,r,vkb(1,ll),uu(1,ll))

          ckb=1.d0/vaverage(ll)

          call dsqv(mx,1,mmax,vkb(1,ll),r,vkb(1,ll))
          do i=1,3
            enl(i)=0.d0
          enddo

c spectrum semilocal case w/ polynomial representation
          call klbyii
     1      (ll-1,mmax,r,vi,vbare(1,ll),0,ckb,vkb(1,ll),enl,ukb)
          write(iu,636) ll-1,(ry2*min(enl(i),0.d0),i=1,3)
          if(tprint .and. ll-1 == 0) then
            open(11,file='tmp.ups-s',status='unknown')
            do j=1,ms
              if(enl(j) > -0.001) exit
              do i=1,mmax
                if(r(i) > 10.0) exit
                write(11,'(e12.6,1x,e12.6)') r(i),ukb(i,j)
              enddo
              write(11,*) '&'
            enddo
            close(11)
          endif
c spectrum nonlocal case w/ polynomial representation
          call klbyii
     1      (ll-1,mmax,r,vi,vbare(1,ll_loc),1,ckb,vkb(1,ll),enl,ukb)
          if(tprint .and. ll-1 == 0) then
            open(11,file='tmp.ukb-s',status='unknown')
            do j=1,ms
              if(enl(j) > -0.001) exit
              do i=1,mmax
                if(r(i) > 10.0) exit
                write(11,'(e12.6,1x,e12.6)') r(i),ukb(i,j)
              enddo
              write(11,*) '&'
            enddo
            close(11)
          endif
          write(iu,637) ll-1,(ry2*min(enl(i),0.d0),i=1,3)

        endif

      enddo
635   format(/' --- kb potentials: spectrum of bound states (eV) ---'//
     1        '            l          e0            e1            e2')
636   format('semilocal ',i3,3(2x,f12.4))
637   format('nonlocal  ',i3,3(2x,f12.4))

c Kleinman Bylander analysis following Gonze et al

  629 format(1x,'note: for the semilocal ',a1,' component',
     1  ' no bound ',a1,' state'/
     1  7x,'is found, the ghost state analysis will assume'/
     1  7x,'a zero reference energy (variable eref)'/)
  630 format(1x,'note: for the local potential',
     1  ' no bound ',a1,' state is'/
     1  7x,'found, the ghost state analysis will assume'/
     1  7x,'a zero groundstate energy (variable eloc0)'/)
  632 format(1x,'note: for the local potential', 
     1  ' no bound excited ',a1,' state'/,7x,'is found,',
     1  ' the ghost state analysis will assume'/
     1  7x,'a zero 1st excited state energy (variable eloc1)'/)
  633 format(' --- searching reference state ---')

c screened local potential
      if(tkb) call dadv(mx,1,mmax,vi,vbare(1,ll_loc),vsl)

      do ll=llbeg,min(llend,llmx)
 
        if(tkb .and. ll .ne. ll_loc) then

          call labels(1,ll-1,1,sn,sl,sm)
          write(iu,639) sl
  639  format(/' --- analysis of kb potentials: ',a1,' waves  ---'/)

c check for input reference spectrum (screened semilocal potential)
c if eigenvalue comes as zero assume it needs to be evaluated
          if(e(ll) .ge. 0.0) then
            write(iu,633) 
            call dadv(mx,1,mmax,vi,vbare(1,ll),wc)
            call dftseq(8,0.d0,mmax,r,ll,ll-1,1.d0,wc,wa,wb,
     1          nin,mch,uld,e(ll),u,up,upp)
            if( e(ll) .ge. 0.0) write(iu,629) sl,sl
          endif

          vnorm=dmelm(mmax,al,r,vkb(1,ll),vkb(1,ll))
          ekb=vnorm/vaverage(ll)
          ckb=vaverage(ll)/sqrt(vnorm)

c local potential groundstate
          eloc0=e(ll)
          call dftseq(8,0.d0,mmax,r,ll,ll-1,1.d0,vsl,wa,wb,
     1          nin,mch,uld,eloc0,u,up,upp)
          if( eloc0 .ge. 0.0) write(iu,630) sl
c local potential first excited state
          eloc1=eloc0
          call dftseq(8,0.d0,mmax,r,ll+1,ll-1,1.d0,vsl,wa,wb,
     1          nin,mch,uld,eloc1,u,up,upp)
          if( eloc1 .ge. 0.0) write(iu,632) sl

c ghost state criteria, e(ll) is the bound state (reference) energy
          if(ekb .gt. 0.0) then
            if(e(ll) .ge. eloc0) then
              if(e(ll) .lt. eloc1) then
                write(iu,*) '* no ghost (ekb > 0, eloc0 < eref < eloc1)'
              else
                if(eloc1 .lt. 0.0) then
              write(iu,*) '* one or more ghosts (ekb > 0, eref > eloc1)'
                else
              write(iu,*) '* undetermined (ekb > 0, eref = eloc1 = 0)'
              write(iu,*) '  note: to decide, inspect log derivatives'
                endif
              endif
            else
              write(iu,*) '* illdefined (ekb > 0, eref < eloc0)'
            endif
          else
            if(e(ll) .lt. eloc0) then
              write(iu,*) '* no ghost (ekb < 0, eref < eloc0)'
            else 
              if(eloc0 .lt. 0.0) then
              write(iu,*) '* one or more ghosts (ekb < 0, eref > eloc0)'
              else
            write(iu,*) '* undetermined (ekb < 0, eref = eloc0 = 0)'
            write(iu,*) '  note: to decide, inspect log derivatives'
              endif
            endif
          endif

          write(iu,*) 
          write(iu,640) 'kb cosine',ckb
          write(iu,640) 'kb energy',ekb*ry2 ,'eV','ekb'
          write(iu,640) 'local potential groundstate'
     1,     eloc0*ry2,'eV','eloc0'
          write(iu,640) 'dto. 1st excited state',eloc1*ry2, 'eV','eloc1'
          write(iu,640) 'reference energy',e(ll)*ry2,'eV','eref'
  640  format(a30,2x,f12.4,1x,a2,4x,a5)
        endif

      enddo

      if(.not. tlgd) return
          
c logarithmic derivatives
c diagnostic radius rld
      do i=1,mmax,2
        if(rld .le. r(i)) goto 202
      enddo
      stop 'pslp - radius for log derivative > r(mmax)'
  202 ild=i
      if(mod(2,ild) .eq. 0) ild=ild+1

      write(iu,642) r(ild)
  642 format(/' --- logarithmic derivatives: at radius =',f7.4,' ---'/)

c energy mesh (should be extensive and dense for finding ghost states)
c for full and semilocal potential an adaptive mesh is used (coarser) 
      iep=870
      es1=0.030
      es2=0.003
      ecc(1)=-2.0
      eccmx=2.0
      if(iep .gt. mxe) stop 'ppcheck - dimension mxe too small'
      do i=1,iep-1
        if(ecc(i) .gt. -1.0 .and. ecc(i) .lt. 1.0) then
          ecc(i+1)=ecc(i)+es2
        else
          ecc(i+1)=ecc(i)+es1
        endif
        if(ecc(i+1) .gt. eccmx) goto 301
      enddo
  301 iep=i

c nonlocal (kleinman bylander) potentials
      if(tkb) write(iu,*) '--- nonlocal potentials ---'
      do ll=llbeg,min(llend,llmx)

        if(tkb .and. ll .ne. ll_loc) then
          unit=70+ll-1
          call dsqv(mx,1,mmax,vkb(1,ll),r,vkb(1,ll))
          do i=1,iep
            call derlkb(ll-1,ecc(i),ild,uld,vaverage(ll),unused,
     1        nodes,mmax,r,vkb(1,ll),vsl)
            write(unit,'(2e15.7)') ecc(i),uld
          enddo
          close(unit)
        endif

      enddo

c all-electron potential (kind=1), semilocal potentials (kind=2)
      do kind=1,2

        irl=3
        if(tnrl .or. kind .eq. 2) irl=6
        if(kind .eq. 1) write(iu,*) '--- all-electron potential ---'
        if(kind .eq. 2) write(iu,*) '--- semilocal potentials ---'

        do ll=llbeg,llend

          unit=50+ll-1
          if(kind .eq. 2) then
            unit=60+ll-1
            if(ll .le. llmx) call dadv(mx,1,mmax,vi,vbare(1,ll),vae)
            if(ll .gt. llmx) call dadv(mx,1,mmax,vi,vbare(1,ll_loc),vae)
          endif

          do i=1,iep

            call dftseq(irl,z,mmax,r,ll,ll-1,1.d0,vae,wa,wb,
     1          ild,ild,uld,ecc(i),u,up,upp)
            uld=up(ild)/(u(ild)*al)-1.0
c addaptive stepsize using energy derivative of log-derivative (friedel sum rule)
            gsum=0.
            do k=1,ild-1
              gsum=gsum+u(k)*u(k)*r(k)
            enddo
            gsum=-2*r(k)*al*(gsum+0.5*u(k)*u(k)*r(k))/(u(k)*u(k))
            delta=abs(0.5/gsum)
            if(abs(uld) .gt. 20.) delta=deltaold
            deltaold=delta
c comment out next line if fixed mesh is to be used
c           ecc(i+1)=ecc(i)+delta
            ecc(i+1)=ecc(i)+delta
            write(unit,'(2e15.7)') ecc(i),uld
            if(ecc(i) .gt. eccmx) goto 401
          enddo
  401     close(unit)

        enddo

      enddo
c            
      return
      end
c
