
program cf_test

    use numtype
    use taylor_cf_approx
    implicit none
    integer :: i, n
    real(dp) :: x 


    print *,'  exp '

    n = 20
    x = 4

    t_coef(0) = 1
    do i = 1, n
        t_coef(i) = t_coef(i-1)/i
    end do 
    !print *,1/t_coef(0:5)

    call taylor_cfrac(t_coef,n,cf_coef)

    print *, x, exp(x), evalcf(cf_coef,n,x)

end program cf_test



