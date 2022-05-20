      !=============================================
      ! 1st derivative of Legendre polynomial
      !=============================================
      function legendre_p(n,t) result (p)
        ! from Mathmatica 
        !  n = 3;
        !  f[t_] := LegendreP[n, t]
        !  f[t] // Expand
        !  f'[t]
        !  f''[t]
        implicit none 
        integer :: n
        real(8) :: t,p,legendre
        real(8) :: y1,y2,step=1e-5

        y1=legendre(n,t+step)
        y2=legendre(n,t-step)
        p = (y1-y2)/2.d0/step
        return 
        
      end function legendre_p
      
