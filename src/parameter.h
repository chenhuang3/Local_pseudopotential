c $Header:$
c parameter file for pseudopotential routines 
c mx .......... radial mesh index
c ms .......... states index
      integer   mx,ms,mzmx
      parameter (mx=1200)
      parameter (ms=80)
      parameter (mzmx=120)
c ie .......... error output unit
c iogncpp ..... regular gncpp output unit 
c iopslp ...... regular pslp output unit 
      integer   ie,iofhipp,iopslp
      parameter (ie=6)
      parameter (iofhipp=23)
      parameter (iopslp=6)
