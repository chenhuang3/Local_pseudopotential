c $Header:$
************************************************************************
      real *8 function bessj(l,x)
*
*     Calculates bessj(l,x) = x*j(l,x), where j(l,x) is the Spherical
*     Bessel Function of integer order l.
*     l=0,2 : The analytic expressions in sin(x) & cos(x) are used.
*     l>2   : The functional value is obtained by downward recursion,
*             if x<l, by upward recursion if x>l, calling 'besrec()'.
*     x<xmin,l>0: bessj()=x**(l+1)/(2*l-1)!!
*     If the argument is excessively small, 'bestiny' sets the output
*     equal to the lowest order term in the powerseries expansion.
************************************************************************
      implicit  real*8 (a-h,o-z)
      parameter (xmin=1.0d-10)
      external  besrec,bestiny
      intrinsic sin,cos

      if(l.eq.0)then
         bessj=sin(x)
      else if(x.lt.xmin)then
         bessj=x*bestiny(l,x)
      else if(l.eq.1)then 
         bessj=sin(x)/x-cos(x) 
      else if(l.eq.2)then 
         bessj=(3.0d0/x/x-1.0d0)*sin(x)-3.0d0*cos(x)/x
      else 
         bessj=x*besrec(l,x)
      endif

      return
      end


************************************************************************
      real *8 function besrec(l,x)
*
*     Recursion for Bessel functions. 'acclog' is about the negative
*     log of the desired accuracy.
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter(acclog=1.0d1)
      intrinsic sin,cos,dble

      if(x.gt.dble(l))then
*
* upward recursion
*
         downj=sin(x)/x
         tempj=(downj-cos(x))/x
         do 100 m=2,l
            upj=(2.0*m-1.0)/x*tempj-downj
                 downj=tempj
                 tempj=upj
100      continue
         besrec=upj
      else
*
* downward recursion & renormalization
*
         mstart=2*l+int(acclog*sqrt(dble(l)))
200      upj=0.0d0
         tempj=1.0d0
         do 300 m=mstart-1,0,-1
            downj=(2.0*m+3.0)/x*tempj-upj
            if(m.eq.l)besrec=downj
            upj=tempj
            tempj=downj
300      continue
         besrec=besrec/downj*sin(x)/x
      endif

400   return
      end


************************************************************************
      real *8 function bestiny(l,x)
*
*     Lowest order powerseries term for spherical Bessel function.
************************************************************************
      implicit real*8 (a-h,o-z)

      kdiv=1.0d0
      do 100 k=1,2*l-1,2
         kdiv=kdiv*k
100   continue
      bestiny=x**l/kdiv

      return
      end


