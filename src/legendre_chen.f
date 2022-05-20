      function legendre(n,t) result (p)
        ! from Mathmatica 
        !  n = 3;
        !  f[t_] := LegendreP[n, t]
        !  f[t] // Expand
        !  f'[t]
        !  f''[t]
        implicit none 
        integer :: n,i
        real(8) :: t,p,p0,p1

        if (n==0) then 
          p=1.d0 
        elseif (n==1) then
          p=t
        else
          ! recursive method 
          p0=1; p1=t
          do i=2,n
            p=((2.d0*(i-1)+1.d0)*t*p1-(i-1)*p0)/dble(i)
            p0=p1
            p1=p
          enddo
        endif 
        return 

      end function legendre
      
