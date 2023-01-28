program residue 
    
    use numtype 
    implicit none

    real(dp) :: a ,b ,c, d , r1, r2 
    complex(dp) :: Rt1, Rt2, x1, x2




    ! Poles calculations
    write(*,*) 'Enter values of a, b, c:'

    read*, a, b, c 

    d = b**2 - 4 * a * c

    if ( d == 0 ) then
        print *, 'Roots are real, equal'
        
        r1 = -b / (2*a)

        print *, 'root is:', r1

    elseif (d > 0) then 
        print *, 'Roots are real, equal'

        r1 = (-b + sqrt(d)) / (2 * a)
        r2 = (-b - sqrt(d)) / (2 * a)

        print *, 'roots are unequal but real:' , r1, r2 

    elseif (d < 0) then

        print *, 'roots are complex'
        
        print *, 'roots are:' , -b / (2*a), ' +i ', (sqrt(abs(d))) / (2 * a) , '&' , -b / (2 *a) , '-i', (sqrt(abs(d))) / (2 * a)

        ! roots are by i (or -i) * (sqrt(abs(d))) / (2 * a)

    end if


    ! Residue Calculation
    x1 = (sqrt(abs(d))) / (2 * a) * iic ! root 1
    x2 = -(sqrt(abs(d))) / (2 * a) * iic ! root 2

    Rt1 = ( 1 + exp(-x1**2) ) / (x1 + iic)
    Rt2 = ( 1 + exp(-x2**2) ) / (x2 - iic)

    print *, 'Residue 1:', Rt1
    print *, 'Residue 2:', Rt2

end program residue
