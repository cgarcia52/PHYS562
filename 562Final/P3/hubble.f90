program regression

    use numtype
    implicit none

    real(dp) :: b, m, r, x, y
    real(dp) :: xsum = 0.0_dp, x2sum = 0.0_dp, xysum = 0.0_dp, ysum = 0.0_dp, y2sum = 0.0_dp, n = 0.0_dp
    character (len=80) :: str 
    

    write(unit=*, fmt="(a)") "Perform linear regression"
    write(unit=*, fmt="(a/)") "Enter FIN to stop and compute regression "

    do

        write (unit=*, fmt="(a)", advance="no") "Enter x, y values: "
        read(unit=*, fmt="(a)") str
        if(str == "fin" .or. str == "FIN") exit
        read (unit=str, fmt=*) x, y

        n = n + 1.0_dp

        xsum = xsum + x
        x2sum = x2sum + x * x
        xysum = xysum + x * y
        ysum = ysum + y
        y2sum = y2sum + y * y
    
    end do

    m = (n * xysum - xsum * ysum) / (n * x2sum - xsum**2)
    b = (ysum * x2sum - xsum * xysum) / (n * x2sum - xsum**2)
    r = (xysum - xsum * ysum /n) / &
                        sqrt((x2sum - xsum**2/n) * (y2sum - ysum**2/n))
    
    write (unit=*, fmt="(/a,es15.6)") " Slope       m = ", m
    write (unit=*, fmt="(a, es15.6)") " y-intercept b = ", b
    write (unit=*, fmt="(a, es15.6)") " Correlation r = ", r

end program regression
    

