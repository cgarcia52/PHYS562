
module pendulum_setup

    use numtype
    implicit none

    integer, parameter :: n_eq = 3
    real(dp) :: g = 0.9, q = 2, omega_d = 1
end module pendulum_setup


program driven_pendulum

    use pendulum_setup
    implicit none
    real(dp), dimension(n_eq) :: y
    real(dp) :: t, tmax, dt


    t = 0
    tmax = 50
    dt = 5

    y(1) = 30 * pi/180 ! \theta
    y(2) = 0           ! \omega = \dot(theta)
    y(3) = 1

    do while ( t < tmax )

        write( 77, * ) t, dt

        write( 21, *  ) t, y(1)
        write( 22, *  ) t, y(2)
        write( 23, *  ) y(1), y(2)

        call rkf45step(t,dt,y)


    end do



end program driven_pendulum



