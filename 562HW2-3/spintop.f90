subroutine rkf45step(x,h,y)  !  Runge-Kutta-Fehlberg 4-5  step

	use numtype, only : dp, n_eq
	implicit none
	real(dp), intent(inout) :: x, h
	real(dp), dimension(n_eq), intent(inout) :: y
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, k5, k6, y1, y2
	real(dp), parameter :: alpha = 1.e-8_dp, delta1 = 1.e-20_dp
	real(dp) :: R, delta
	
	interface

        function spin ( x, h, y)

            use numtype, only : dp, n_eq
            implicit none
            real(dp), intent(in) :: x
            real(dp), intent(in) :: h
            real(dp), intent(in), dimension(n_eq) :: y
            real(dp), dimension(n_eq) :: spin

        end function spin

    end interface


	k1 = spin ( x, h, y )
	k2 = spin ( x + h/4, h, y + k1/4)
	k3 = spin ( x + 3*h/8, h, y + (3*k1+9*k2)/32)
	k4 = spin ( x + 12*h/13, h, y+(1932*k1-7200*k2+7296*k3)/2197)
	k5 = spin ( x + h, h, &
	    y + 439*k1/216-8*k2+3680*k3/513-845*k4/4104 )
	k6 = spin ( x + h/2, h, y + &
	    (-8*k1/27 +2*k2-3544*k3/2565 +1859*k4/4104 -11*k5/40))

    y1 = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104  - k5/5
    y2 = y + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 &
        - 9*k5/50 + 2*k6/55

    R = sqrt(dot_product(y1-y2,y1-y2))/h + delta1

    if ( R < alpha ) then
        x = x + h
        y = y1
        delta = 0.92_dp * (alpha/R)**(0.2_dp)
        h = delta*h
    else
        delta = 0.92_dp * (alpha/R)**(0.25_dp)
        h = delta*h
    end if

end subroutine rkf45step