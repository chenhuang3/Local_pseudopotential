      !=============================================
      ! second derivative of Legendre polynomial
      !=============================================
      function legendre_pp(n,t) result (p)
        ! from Mathmatica 
        !  n = 3;
        !  f[t_] := LegendreP[n, t]
        !  f[t] // Expand
        !  f'[t]
        !  f''[t]
        implicit none 
        integer :: n
        real(8) :: t,p,y1,y2,y3,step=1e-5,legendre

        y1=legendre(n,t+step)
        y2=legendre(n,t)
        y3=legendre(n,t-step)
        p = (y1-2.d0*y2+y3)/step**2
        return 

      end function 
      
