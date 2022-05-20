c $Header:$
cmf 27-06-03 added elements up to 103
c***********************************************************************
c atomic number & xc type - number to symbol conversion map
c
c iz....nuclear charge to elemental symbol symz
c ixc...exchange-correlation type to literal descriptor symxc
c
c Martin Fuchs, FHI der MPG, 11-1996
c***********************************************************************
      subroutine labelmap(iz,symz,ixc,symxc)
      implicit none
      include  'parameter.h'
      integer iz,ixc,ixcmx
      character*8 symz,sz
      character*40 symxc,sxc
      parameter(ixcmx=20)
      dimension sz(mzmx),sxc(0:ixcmx)
      data sxc/
     + 'no Hartree+xc (hydrogenic atom)'
     +,'LDA Wigner'
     +,'LDA Hedin/Lundqvist'
     +,'LDA CA Perdew/Zunger 1980'
     +,'GGA Perdew/Wang 1991'
     +,'GGA Becke/Perdew'
     +,'GGA Perdew/Burke/Ernzerhof'
     +,'LDA CA Perdew/Wang 1991 + X relativistic'
     +,'LDA CA Perdew/Wang 1991'
     +,'GGA X Becke C Lee/Yang/Parr'
     +,'GGA X PW91  C Lee/Yang/Parr'
     +,'LDA exchange only'
     +,'KLI exchange'
     +,'KLI exchange + ECP'
     +,'GGA X Hammer/PBE C PBE'
     +,'GGA X Zhang/Wang/PBE C PBE'
     +,'MGGA Perdew/Kurth post PBE GGA'
     +,'GGA PBE exchange + LDA correlation'
     +,'KLI exchange + LDA correlation'
     +,'KLI exchange + PBE GGA correlation'
     +,'unknown'/
      data sz/
     + 'H ','He'
     +,'Li','Be','B ','C ','N ','O ','F ','Ne'
     +,'Na','Mg','Al','Si','P ','S ','Cl','Ar'
     +,'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn'
     +,     'Ga','Ge','As','Se','Br','Kr'
     +,'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd'
     +,     'In','Sn','Sb','Te','I ','Xe'
     +,'Cs','Ba','La'
     +,               'Ce','Pr','Nd','Pm','Sm','Eu','Gd'
     +,               'Tb','Dy','Ho','Er','Tm','Yb','Lu'
     +,          'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg'
     +,     'Tl','Pb','Bi','Po','At','Rn'
     +,'Fr','Ra','Ac'
     +,               'Th','Pa','U ','Np','Pu','Am','Cm'
     +,               'Bk','Cf','Es','Fm','Md','No','Lw'
     +,'--','--','--','--','--','--','--','--','--','--'
     +,'--','--','--','--','--','--','--'/
       if(iz .le. mzmx .and. iz .gt. 1 ) then
         symz=sz(iz)
       else
         symz='none'
       endif
       if(ixc .le. ixcmx) then
         symxc=sxc(ixc)
       else
         symxc='unknown'
       endif
c
       return
       end
c
c***********************************************************************
c pseudopotential specifier - number to symbol conversion map
c***********************************************************************
       character*20 function spp(ipp)
       implicit none
       integer  ipp,ippmx
       parameter(ippmx=2)
       character*20 symbols_pp(ippmx)
       data symbols_pp/'hamann','troullier-martins'/
       if(ipp .le. ippmx) then
         spp=symbols_pp(ipp)
       else
         spp='what?'
       endif
       return
       end
c
c***********************************************************************
c tags for angular momentum to literal conversion
c***********************************************************************
      subroutine labels(n,l,m,sn,sl,sm)
      integer n,l,m
      character*1 sn,sl,sm,cn(7),cl(6),cm(3)
      data cn/'1','2','3','4','5','6','.'/
      data cl/'s','p','d','f','g','.'/
      data cm/'+','-','.'/
      if(n .lt. 1 .or. n .gt. 6) n=7
      if(l .lt. 0 .or. l .gt. 4) l=5
      if(m .lt. 1 .or. m .gt. 2) m=3
      sn=cn(n)
      sl=cl(l+1)
      sm=cm(m)
      return
      end
c
