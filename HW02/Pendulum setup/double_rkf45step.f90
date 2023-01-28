subroutine double_rkf45step(x,h,y)  !  Runge-Kutta-Fehlberg 4-5  step

	use double_pendulum, only : dp, n_eq
	implicit none
	real(dp), intent(inout) :: x, h
	real(dp), dimension(n_eq), intent(inout) :: y
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, k5, k6, y1, y2
	real(dp), parameter :: epsilon = 1e-8_dp, tiny = 1.e-20_dp ! tiny --20
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

	        use double_pend 
	        implicit none
	        real(dp), intent(in) :: t, h
	        real(dp), dimension(n_eq), intent(in) :: y
	        real(dp), dimension(n_eq) :: f, k

            real(dp) :: a


            a = 1/(1-((m2*((cos(y(1)-y(2)))**2))/(m1+m2)))

			f(1:2) = y(3:4)
			f(3) = a * ((-m2*(y(3)**2)*sin(2*(y(1)-y(2))))/(2*(m1+m2))+&
                        (m2*gravity*sin(y(2))*cos(y(1)-y(2)))/(l1*(m1+m2))+&
                        (-m2*l2*(y(4)**2)*sin(y(1)-y(2))-(m1+m2)*gravity*sin(y(1)))/&
                        (l1*(m1+m2)))
                        

            f(4) = a * ((m2*(y(4)**2)*sin((y(1)-y(2))))/((m1+m2))+&
                        (gravity*sin(y(1))*cos(y(1)-y(2))/l2)+&
                        (l1*y(3)**2*sin(y(1)-y(2))-gravity*sin(y(2)))/l2)



  

	        k(1:n_eq) = h*f(1:n_eq)


			!----------------------------	
        end function kv

end subroutine double_rkf45step