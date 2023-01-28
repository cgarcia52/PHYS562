
module sch_setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 3

    real(dp), parameter :: hbar = 1._dp, mass =1._dp, &
            hbar2 = hbar**2, omega = 1._dp

    integer :: iw
    real(dp) :: energy, xmax, norm, y3

end module sch_setup

program schrodinger

    use sch_setup
    use chebyshev
    implicit none
    real(dp), external :: miss
    real(dp) :: aa, bb, en, goal
    integer :: nc, i 

    xmax = 10
    norm = 1

    aa = 0
    bb = 2
    iw = 0
    nc = 13
    call chebyex(miss,nc,cheb,aa,bb) 
    call chebyzero(nc,cheb,aa,bb,z0,iz0) 
    print *,iz0,z0(1:iz0)

    do i = 1, iz0 
        en = z0(i)
        iw = 0 
        norm = 1
        goal = miss(en)
        norm = sqrt(y3)
        print *,norm

        iw = 10 + i
        goal = miss(en)
        print *, en, goal, y3
    end do 

end program schrodinger 



function miss(en)

    use sch_setup
    implicit none
    real(dp) :: en, miss, x, dx 
    real(dp) :: y(n_eq), kappa

    energy = en
    x = xmax
    dx = 0.05_dp
    kappa = sqrt(2*mass/hbar2 * abs(en))

    y(1) = exp(-x**2/2) / norm
    y(2) = -x*exp(-x**2/2) / norm 
    y(3) = 0


    do while ( x > 0  )


        if(iw /= 0 ) write(unit=iw, fmt='(3f10.3)' ) x,0._dp,y(1)

        call rk4step ( x, -dx, y )

    end do

    miss = y(1) * y(2)
    y3 = y(3)
    print *, en, miss, y(1:3) 

end function miss


subroutine rk4step ( x, h, y )

    use sch_setup
    implicit none
    real(dp), intent(in) :: h
    real(dp), intent(inout) :: x
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

    interface

        function kv ( x, h, y )

            use sch_setup
            implicit none
            real(dp), intent(in) :: x, h
            real(dp), dimension(n_eq), intent(in) :: y
            real(dp), dimension(n_eq) :: kv

        end function kv

    end interface

    k1 = kv (x,    h, y)
    k2 = kv (x+h/2,h, y+k1/2)
    k3 = kv (x+h/2,h, y+k2/2)
    k4 = kv (x+h,  h, y+k3)

    dy = (k1 + 2*k2 + 2*k3 + k4)/6

    y = y + dy
    x = x + h

end subroutine rk4step


function kv ( x, dx, y ) result(k)

    use sch_setup
    implicit none
    real(dp), intent(in) :: x, dx
    real(dp), dimension(n_eq), intent(in) :: y
    real(dp), dimension(n_eq) :: f, k
    real(dp), external :: potential
!-----------------------------------------

!---------------------------------------

    f(1) = y(2)                                    ! 
    f(2) = - 2*mass/hbar2 * ( energy - potential(x)) * y(1) !
    f(3) = -2*abs(y(1))**2

!----------------------------------------
    k = dx * f

end function kv


function potential (x)

    use sch_setup
    implicit none
    real(dp) :: x, potential

    potential = 0.5_dp * mass * omega**2 * x**2

end function potential 


