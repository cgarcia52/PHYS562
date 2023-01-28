module setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 2
    real(dp), parameter :: hbar2 = 1._dp, mass = 1._dp, omega = 1._dp
    real(dp) :: energy

end module

program hold

    use setup
    implicit none

    real(dp) :: xmax, x, dx, psi(n_eq)

    xmax = 3.4
    x = 0
    dx = 0.001_dp

    psi(1) = 1
    psi(2) = 0
    energy = 0.52_dp

    do while (x < xmax)
        write (1,*) x, 0._dp, psi(1)
        call rk4step(x, dx, psi)
    end do

end program hold 

subroutine rk4step(x,h,y)

    use setup
    implicit none
    real(dp), intent(inout) :: x
    real(dp), intent(in) :: h
    real(dp), intent(inout), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

    k1 = kv (x,h,y)
    k2 = kv (x+h/2, h, y+k1/2)
    k3 = kv (x+h/2, h, y+k2/2)
    k4 = kv (x+h, h, y+k3)

    dy = (k1 + 2*k2 + 2*k3 + k4) / 6

    y = y + dy
    x = x + h

    contains 

        function kv (x,dx,y) result(k)

            use setup
            real(dp), intent(in) :: x, dx
            real(dp), intent(in), dimension(n_eq) :: y
            real(dp), dimension(n_eq) :: f, k

            f(1) = y(2)
            f(2) = -2*mass/hbar2 * (energy - potential(x) ) * y(1)
            k = dx * f
    
        end function

        function potential(x)

            use setup
            real(dp) :: x, potential

            potential = 1._dp/2 * mass * omega**2 * x**2

        end function potential

end subroutine rk4step 