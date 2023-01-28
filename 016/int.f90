program integration

    use numtype
    implicit none
    real(dp), external :: runge
    real(dp) :: a, b, res1, res2
    integer :: n

    a = -5
    b = 5
    n = 20

    call trapezoid(a, b, runge, res1, n)
    call simpson(a, b, runge, res2, n)

    print *, n, res1, res2, 2*atan(b)



end program integration

subroutine trapezoid(a, b, fun, res, n)

    use numtype
    implicit none
    real(dp), external :: fun
    real(dp) :: a, b, res, h, x
    integer :: n, nn, i

    h = (b-a)/n
    res = (fun(a)+fun(b))/2
    do i = 1, n-1
        x = a + i*h
        res = res + fun(x)
    end do

    res = res * h

end subroutine trapezoid  


subroutine simpson(a, b, fun, res, n)

    use numtype
    implicit none
    real(dp), external :: fun
    real(dp) :: a, b, res, h, x
    integer :: n, nn, i

    nn = n
    if(mod(nn,2) /= 0) nn = nn + 1


    h = (b-a)/n
    res = fun(a)+fun(b) + 4*fun(a+h)
    do i = 2, n-2, 2
        x = a + i*h
        res = res + 2*fun(x) + 4*fun(x+h)
    end do

    res = res * h/3

end subroutine simpson

function runge(x)

    use numtype
    implicit none
    real(dp) :: x, runge

    runge = 1/(1+x**2)

end function runge