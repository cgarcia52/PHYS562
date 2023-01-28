
program shoot ! solves for Phi and E-field by changing function for each part

    use numtype
    implicit none

    integer, parameter :: n_eq = 2
    real(dp) :: r, dr
    real(dp), dimension(n_eq) :: y 

    r = 10 ! setting up the r boundary from where calculations start
    dr = -0.0001_dp ! integrates it from 6 down to our second bound which we assume to be r = 0

    y(1) = 0._dp
    y(2) = 0._dp
 
    do while ( r > 0 )
        write(1,*) r, y(1)
        write(2,*) r, y(2)

        call rk4step( r, dr, y )
    end do


end program shoot

function kv ( r, dr, y) result(output)

    use numtype
    implicit none

    integer, parameter :: n_eq = 2
    real(dp), intent(in) :: r, dr
    real(dp), intent(in), dimension(n_eq) :: y

    real(dp), dimension(n_eq) :: x, output 
    real(dp), external :: rho

    x(1) = -y(2)
    x(2) = -4 * pi * rho(r) - (2/r) * x(1) 

    output = dr * x
end function kv


function rho(r) result(Vpot)
    use numtype
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: Vpot

    Vpot = 1/(2*pi) * cos(r)* exp(-r) ! change this for a thru d the rho values
end function rho
