program harmonic_oscillator

    use numtype 
    implicit none
    real(dp) :: r, l, d, del, val ! we are setting r = x
    real(dp) :: psi, a, b, c ! the psi function is split up into 3 seperate equations
    integer :: v, n
    

    v = 1 ! these will be changed according to the given values
    n = 2 ! n from 0 to 5, we will test graphs for n= 1, 2, 3
    l = 1 ! l from 0 to 1
    d = 3 ! D from 2 to 3
    r = 0 ! plug in a beginning r value

    del = l + (d/2.0_dp) +n !exponent
    val = l - 1 + (d/(2.0_dp)) !exponent
    
    do while (r < 10)

        r = r + .05 ! step function

        a = ( v**(1.0/4.0) )*( (2*fact(n) ) / GAMMA(del))**(1.0/2.0) ! replace delta
        b = exp((-v*r**2)/2)
        c = (v*r**2)**((l/2) + (d-1)/4)

        psi = a*b*c*Laguerre(n, r, val) !final equation

        print *, a, b, c !check for value consistency

        print *, r, psi ! final values of r (x) and psi (y)
        write(3,*) r, psi

    end do

    contains

        recursive function fact(n) result(s0) !factorial equation

            implicit none
            integer, intent(in) :: n 
            real(dp) :: s0

            if (n<0) then
                stop 'something is wrong'
            else if (n == 0) then
                s0 = 1._dp
            else
                s0 = n * fact(n-1) !simple factorial
            end if

        end function fact


        recursive function Laguerre(n,r,val) result(s0)
        
            implicit none
            integer, intent(in) :: n
            real(dp) :: r
            real(dp) :: val !replace val with actual value

            real(dp) :: s0

            if (n < 0) then
                s0 = 0._dp

            else if (n == 0) then ! note, k = n - 1 so in long equation, k + 1 = n, etc...
                s0 = 1._dp
            else 
                s0 = ( (2*(n-1)+1+val-r**2)*Laguerre(n-1,r,val) - (n-1+val)*Laguerre(n-2,r,val) )/ (n) !laguerre equation for first few terms
            end if

        end function Laguerre
        



end program harmonic_oscillator