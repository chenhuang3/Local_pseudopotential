c $Header:$
c
c gauss-type quadrature
c here: matrix elements
c
c variables
c imx........order of quadrature
c x(i).......i-th abscissa
c w(i).......i-th weight
c func.......user defined external function
c 
c order, abscissas and weights are initialized by include 'gauss.h'
c hint: do any variable substitutions in func itself
c
      real*8 function gaussq(func,n1,n2)
      implicit real*8 (a-h,o-z)
      include 'gauss.h'
c
      gaussq=0.d0
      do i=1,imx
        gaussq=gaussq+w(i)*func(n1,n2,x(i))
      enddo
c
      return
      end
c
