
subroutine rkf45step(x,h,y)  !  Runge-Kutta-Fehlberg 4-5  step

	use setup, only : dp, n_eq
	implicit none
	real(dp), intent(inout) :: x, h
	real(dp), dimension(n_eq), intent(inout) :: y
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, k5, k6, y1, y2
	real(dp), parameter :: epsilon = 1.e-8_dp, tiny = 1.e-20_dp
	real(dp) :: rr, delta

	k1 = kv ( x, h, y )
	k2 = kv ( x + h/4, h, y + k1/4)
	k3 = kv ( x + 3*h/8, h, y + (3*k1+9*k2)/32)
	k4 = kv ( x + 12*h/13, h, y+(1932*k1-7200*k2+7296*k3)/2197)
	k5 = kv ( x + h, h, &
	    y + 439*k1/216-8*k2+3680*k3/513-845*k4/4104 )
	k6 = kv ( x + h/2, h, y + &
	    (-8*k1/27 +2*k2-3544*k3/2565 +1859*k4/4104 -11*k5/40))

    y1 = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104  - k5/5
    y2 = y + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 &
        - 9*k5/50 + 2*k6/55

    rr = sqrt(dot_product(y1-y2,y1-y2))/h + tiny

    if ( rr < epsilon ) then
        x = x + h
        y = y1
        delta = 0.92_dp * (epsilon/rr)**(0.2_dp)
        h = delta*h
    else
        delta = 0.92_dp * (epsilon/rr)**(0.25_dp)
        h = delta*h
    end if

    contains

        function kv (t,h,y)  result(k)

	        use setup, only : dp, n_eq
	        implicit none
	        real(dp), intent(in) :: t, h
	        real(dp), dimension(n_eq), intent(in) :: y
	        real(dp), dimension(n_eq) :: f, k

  

	        k(1:n_eq) = h*f(1:n_eq)

        end function kv

end subroutine rkf45step
