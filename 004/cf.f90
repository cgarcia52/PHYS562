
program cf

    use numtype
    implicit none
    real(dp) :: x

    x = 13.73_dp
    print *, x, exp(x), eexp(x)

    contains

        function eexp(x)

            implicit none
            real(dp) :: x, eexp

            eexp = 1/(1 - scf(1,x))

        end function eexp

        recursive function scf(i,x) result(s0)

            implicit none
            real(dp) :: x, s0
            integer :: i
            integer, parameter :: imax = 140

            if ( i > imax ) then
                s0 = 0
            else 
                s0 = x/i / ( 1+ x/i - scf(i+1,x)  )
            end if

        end function scf

end program cf

