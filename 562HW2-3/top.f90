module spinning

    use numtype
    implicit none

    real(dp), parameter :: G = 9.81_dp, L = 1._dp, m = .5_dp! 3 kg mass
    real(dp), parameter :: I1 = 1._dp, I2 = I1, I3 = .5_dp

end module spinning

program heavy 

    use spinning
    implicit none
    real(dp), dimension(n_eq) :: y
    real(dp) :: t, tmax, dt, w3

    t = 0
    tmax = 5 ! seconds
    dt = 1._dp

    y(1) = 30 * pi/180         ! phi radians
    y(2) = pi/2       ! phi prime (dot) radians 
    y(3) = 60 * pi/180           ! theta
    y(4) = 45 * pi/180         ! theta prime (dot) radians 

    do while (t < tmax )

        write(1,*) t, y(1) ! time, phi radians, plotting this one
        !write(2,*) t, y(2) ! time, phi dot
        write(3,*) t, y(3) ! time, theta, plotting this one
        !write(4,*) t, y(4) ! time, theta dot

        !write(11,*) y(1), y(2) ! (phi, phi dot)
        !write(22,*) y(3), y(4) ! (theta, theta dot)
        
        !write(33,*) t, dt ! time, time step

        call rkf45step( t, dt, y)

    end do

end program heavy

        function spin( t, dt, y) result(k)

            use spinning
            implicit none
            real(dp), intent(in) :: t ! time
            real(dp), intent(in) :: dt ! time step
            real(dp), intent(in), dimension(n_eq) :: y ! angle measurements
            real(dp), dimension(n_eq) :: f, k ! k in kunga
            real(dp) :: w3, psidot

            psidot = 10

            w3 = y(2)*cos(y(3)) + psidot ! 3rd component of the omega vector (1 x 3 matrix)

            f(1) = y(2) !phi dot, updates

            f(2) = ((I3 * w3 *y(4) ) - 2*I1* y(2) * y(4) *cos(y(3) ))/(I1*sin(y(3))) ! eq of motion where phi double dot is f(2)

            f(3) = y(4) ! theta dot, updates

            f(4) = ((I1* y(2)**2 * cos(y(3)) - I3*w3*y(2) + m*G*L)*sin(y(3)))/I1 ! eq of motion where f(4) = theta double dot
            

            k = dt * f ! time derivative 

        end function spin