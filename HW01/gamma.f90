program gamma_function 

    use numtype
    implicit none 
    real(dp) :: x, dx, s

    x = .9_dp ! complex numbers
    s = 1._dp ! real numbers
    dx = 0.1_dp ! step functions

    do while (x<15) ! for x is less than 15

        x = x + dx ! ends when it hits 15

        write(1,*) x, lower(s,x)
            !lower function data, x, lower
        write(2,*) x, upper(s,x)
            !upper function data, x, upper
        write(3,*) x, upper(s,x) + lower(s,x)
            ! checks for addition of gammas against x values

    end do 

    print *, "lower gamma value:", lower(s,x) 
    ! gives out a value for the lower and upper gamma function
    print *, "upper gamma value:", upper(s,x)

    print *, "Total Gamma value:", upper(s,x) + lower(s,x)

        contains

            recursive function upper(s,x) result(s1) 
            !for x to infinity

                implicit none
                real(dp) :: s1 
                real(dp), intent(in) :: x, s 
                    !declare these as already being defined earlier
                integer :: n, i

                n = 500 ! max variable start
                s1 = 0._dp ! starts at 0

                do i = n , 1, -1 !range, ends at 1

                    s1 =  x + (  i-s)/1 + i/(x + (i + 1 - s)/1 + (i + 1))/s1
                        ! upper gamma pattern for both positive and negative

                end do 

                s1 = ( exp(-x) * x ** s)/s1
                ! complete equation

            end function upper


            recursive function lower(s,x) result(s2)
                !for 0 to x of the gamma function

                implicit none
                real(dp) :: s2
                real(dp), intent(in) :: x, s
                    ! these values are defined
                integer :: n, i

                n = 500
                s2 = 0._dp

                do i = n , 0, -1 !range, ends at 0

                    s2 =  (( s + i) * x) / ( i + 1 + s + x - s2)
                    ! equation for lower gamma

                end do

                s2 = (exp(-x) * x**s)/ (s - s2)
                ! full equation


            end function lower


end program gamma_function