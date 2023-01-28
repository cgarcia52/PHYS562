program taylor_fun

    use numtype
    implicit none

    real(dp), dimension(0:100) :: coeff

    integer :: i, imax

    real(dp) :: x, y, res

    imax = 50

    ! exp(x)

    coeff(0) = 1
    do i = 1, imax
        coeff(i) = coeff(i-1)/i 
    end do

    !print *, 1/coeff(0:5)

    x = 5._dp 
    y = coeff(imax)

    do i = imax-1, 0, -1
        y = coeff(i) + x*y 
    end do 

    print *, 'Exponential function ', x, exp(x), y, horner(x, imax, coeff)

    !-----------------------------

    ! log(1+x)

     coeff(0) = 0
    do i = 1, imax
        coeff(i) = -(-1)**(i) *1._dp/i 
    end do 

    !print *, 1/coeff(0:5)

    x = 0.9_dp 
    y = coeff(imax)

    do i = imax-1, 0, -1
        y = coeff(i) + x*y 
    end do 

    print *, 'log (1+x) fucntion ', x, log(1+x), y, horner(x,imax,coeff)




    




end program taylor_fun

function horner(x,nmax,an) result(y) 

    use numtype
    implicit none

    real(dp), intent(in) :: x
    integer, intent(in) :: nmax
    real(dp), dimension(0:namx), intent(in) :: an
    real(dp) :: y
    integer :: i 

    y = coeff(imax)


    do i = imax-1, 0, -1
        y = an(i) + x*y
    end do


end function horner 