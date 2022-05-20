      !
      ! satuation function F(x)=(x-xcut)^2 for x>xcut
      !                    F(x)=0 for x<xcut and x>-xcut
      !                    F(x)=(x+xcut)^2 for x<-xcut 
      !
      real*8 function satf(x,xcut)
        implicit none 
        real*8 x, xcut 

        if (x>=xcut) satf=(x-xcut)**2
        if (x<xcut .and. x>-xcut) satf=0.d0
        if (x<=-xcut) satf=(x+xcut)**2

      end function
