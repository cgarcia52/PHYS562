program chebytest

use chebyshev
implicit none
real(dp) :: chder2(0:maxch)
integer :: n, i ,np
real(dp), external :: func, dfunc
real(dp) :: x, dx, ya, yb, eps, zz ,dz
n = 13
ya = -pi/4
yb = 9*pi /4
call chebyex(func, n, cheb, ya, yb)

np = 50
dx = (yb-ya)/np
do i = 1, np+1
    x = ya+(i-1)*dx
    write(1,*) x, 0._dp, func(x), cheby(x, cheb, n, ya, yb)
end do

call chebyzero(n, cheb, ya, yb, z0, iz0)
eps = 1.e-15_dp
dz = 0.01_dp
do i =1, iz0
    zz = z0(i)
    call root_polish(func,zz,dz,eps,10)
    print *, 'f(x) = 0', i, z0(i),zz
end do

print *, ' derviative '

call chebyderiv(cheb, n, chder, ya, yb)
call chebyderiv(chder, n-1, chder2, ya, yb)

do i =1, np+1
    x = ya + (i-1) * dx
    write(2,*) x, 0, dfunc(x), cheby(x,chder,n-1,ya,yb)
end do

call chebyzero(n-1,chder,ya,yb,z0,iz0)
do i  = 1, iz0
    zz = z0(i)
    call root_polish(dfunc,zz,dz,eps,10)
    print *, ' dx(x)/dx = 0 ', i , z0(i), zz, cheby(z0(i), chder2, n-2, ya, yb)
end do

end program chebytest

function func(x) result(f)

    use numtype, only : dp
    implicit none
    real(dp) :: x, f
    f = x*cos(x)

end function func

function dfunc(x) result(df)

    use numtype, only : dp
    implicit none
    real(dp) :: x, df
    df = cos(x)-x*sin(x)

end function dfunc