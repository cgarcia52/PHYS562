module setup

    use numtype
    implicit none
    integer, parameter :: n_eq = 4
    real(dp), parameter :: g = 10.0_dp, mass = 1.0_dp, &
        drag = 0.1_dp, distance = 4._dp, wind(2) = (-5._dp, 0._dp)
    real(dp), parameter :: x0 = 0._dp, y0 = 0._dp, &
        v0 = 20._dp
    integer :: iw

end module setup

program targ

    use setup
    use chebyshev
    implicit none
    real(dp), external :: proj
    real(dp) :: alpha, yx, aa, bb, da, eps
    integer :: n, iz, maxf


    aa = 0
    bb = 90
    iw = 0
    n = 5
    call chebyex(proj, n, cheb, aa, bb)
    call chebyzero(n,cheb, aa, bb, z0, iz0)

    eps = 0.01_dp
    maxf = 10
    da = 0.1_dp
    do iz = 1, iz0
        alpha = z0(iz)
        iw = 10 + iz
        !call root_polish(proj, alpha, da, eps, maxf)
        yx = proj(alpha)
        print *, iz, alpha, yx
    end do

end program targ

function proj(alpha)

    use setup
    implicit none
    real(dp) :: alpha, proj
    real(dp) :: t, dt
    real(dp), dimension(n_eq) :: y

    t = 0._dp
    dt = 0.05_dp
    y(1) = 0
    y(2) = mass*v0*cos(alpha/180*pi)
    y(3) = 0.1
    y(4) = mass*v0*sin(alpha/180*pi)

    do while ( y(3) >= 0._dp)
        if(iw /= 0) write(unit=iw, fmt = '(2f15.5)') y(1), y(3)
        call rk4step(t,dt,y)
    end do

    proj = y(1) - distance
    print *, alpha, proj

end function proj

subroutine rk4step(x,h,y)

    use setup, only : dp, n_eq
    implicit none
    real(dp), intent(inout) :: x
    real(dp), intent(in) :: h
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

    k1 = kv(x, h, y)
    k2 = kv(x+h/2, h, y+k1/2)
    k3 = kv(x+h/2, h, y+k2/2)
    k4 = kv(x+h, h, y+k3)

    dy = (k1 + 2*k2 + 2*k3 + k4)/6

    x = x + h
    y = y + dy

    contains 

        function kv (t,dt,y) result(k)

            use setup, only : dp, n_eq, g, mass, drag, wind
            implicit none
            real(dp), intent(in) :: t, dt
            real(dp), dimension(n_eq), intent(in) :: y
            real(dp), dimension(n_eq) :: f, k
            real(dp) :: vv(2), vrel(2), v

            vv(1) = y(2)/mass; vv(2) = y(4)/mass
            vrel = vv - wind
            v = sqrt(vrel(1)**2 + vrel(2)**2)

            f(1) = vv(1)
            f(2) = -drag * v * vrel(2) - mass * g
            k = dt * f
        
        end function kv

end subroutine rk4step

