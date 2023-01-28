program random

    use numtype
    implicit none

    real(dp) :: a, b, integral, integral_error
    !real(dp) :: r(2)
    integer :: i, n, seed

    seed  = 15132604
    call random_seed(seed)


    a = 0
    b = pi

    call montecarlo(n, a, b, fun, integral, integral_error)
    print *, integral, integral_error


    contains 

        subroutine montecarlo(n, a, b, func, integral, integral_error)

            implicit none
            integer, intent(in) :: n
            real(dp), intent(in) :: a, b
            real(dp) :: integral, integral_error
            real(dp), external :: func
            real(dp) :: r, x, fx, f, f2
            integer :: i


    !n = 100000
    f = 0
    f2 = 0
    do i = 1, n
        call random_number(r)
        x = a + r*(b-a)
        fx = func(x)
        fx = func(x)
        f = f + fx
        f2 = f2 + fx**2
        !call random_number(r)
        !write(1,*) r
    end do 
    f = f/n
    f2 = f2/n

    integral = (b-a) * f
    integral_error = (b-a) * sqrt((f2-f**2)/n)

        end subroutine montecarlo

        function fun(x)
            implicit none
            real(dp) :: x, fun
            fun = sin(x)
        end function fun

end program random