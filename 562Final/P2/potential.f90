program potential

    use numtype
    implicit none

    integer, parameter :: n_eq = 4
    real(dp) :: r, dr
    real(dp), dimension(n_eq) :: y

    r = 10
    dr = -0.0001_dp

    y(1) = 0._dp
    y(2) = 0._dp
    y(3) = 0._dp
    y(4) = 0._dp

    do while (r > 0)
        write(1,*) r, y(1)
        write(2,*) r, y(2)
        write(3,*) r, y(3)
        write(4,*) r, y(4)
        call rk4step(r,dr,y)
    end do

end program potential

function kv (r,dr,y) result(output)
    
    use numtype
    implicit none

    integer, parameter :: n_eq = 4
    real(dp), intent(in) :: r, dr
    real(dp), intent(in), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: x, output
    real(dp), external :: rho

    x(1) = (1 / (8 * pi)) * exp(-r) 
    x(2) = (1 / (24 * pi)) * exp(-r)
    x(3) = (1 / (2 * pi)) * sin(r) * (exp(-r))
    x(4) = (1 / (2 * pi)) * cos(r) * (exp(-r))

    output  = dr * x

end function kv

function rho(r) result(Vpot)
    
    use numtype
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: Vpot

    Vpot  = 1/(2*pi) * cos(r) * exp(-r)

end function rho