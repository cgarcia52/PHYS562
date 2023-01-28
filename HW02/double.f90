module setup

    use numtype
    implicit none
    !integer, parameter :: n_eq 
    real(dp), parameter :: g = 10._dp, length = 10._dp
    real(dp), parameter :: q = 0.1_dp !friction coeff

end module setup

program double_pend

    use setup
    implicit none
    real(dp), dimension(n_eq) :: y
    real(dp) :: t, tmax, dt

    t = 0._dp
    tmax = 60._dp
    dt = 0.01_dp

    y(1) = 90*pi/180
    y(2) = 0._dp

    do while (t < tmax)

        write(7,*) y(1), y(2)
        write(8,*) t, y(1)
        write(9,*) t, y(2)

        call rkf45step(t,dt,y)

    end do

end program double_pend

subroutine rkf45step(x,h,y)

    use setup, only : dp, n_eq
    implicit none
    real(dp), intent(inout) :: x, h
    real(dp), dimension(n_eq), intent (inout) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, k5, k6, y1, y2
    real(dp), parameter :: epsilon = 1.e-8_dp, tiny = 1.e-20_dp
    real(dp)  :: rr, delta

    k1 = kv (x,h,y)
    k2 = kv (x + h/4, h, y + k1/4)
    k3 = kv (x + 3*h/8, h ,y + (3*k1 + 9*k2)/32)
    k4 = kv (x + 12*h/13, h, y + (1932*k1 - 7200*k2 + 7296*k3)/2197)
    k5 = kv (x+h, h, y + 439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
    k6 = kv (x+h, h, y + (-8*k1/27 + 2*k2 - 3544*k3/2565 + 1859*k4/4104 - 11*k5/40))

    y1 = y + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5
    y2 = y + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 -9*k5/50 + 2*k6/55

    rr = sqrt(dot_product(y1-y2,y1-y2))/h + tiny

    if (rr < epsilon) then
        x = x + h
        y = y1
        delta = 0.92_dp * (epsilon/rr)**(0.2_dp)
        h = delta*h
    else 
        delta = 0.92_dp * (epsilon/rr)**(0.25_dp)
        h = delta*h
    end if 

    contains

        function kv (t,h,y) result(k)

            use setup, only : dp, n_eq, g, length, q
            implicit none
            real(dp), intent(in) :: t, h
            real(dp), dimension(n_eq), intent(in) :: y
            real(dp), dimension(n_eq) :: f, k

            f(1) = y(2)
            f(2) = -g/length * sin(y(1)) - q*y(2)

            k(1:n_eq) = h*f(1:n_eq)

        end function kv
end subroutine rkf45step