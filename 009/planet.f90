module planet_setup

    use numtype
    implicit none


    integer, parameter :: n_eq = 6
    real(dp) :: gravity = 6.6743e-11, mass_sun = 1.9885e+30_dp, mass_Earth = 5.97219e+24_dp


end module planet_setup


program kepler

    use planet_setup
    implicit none
    real(dp), dimension(n_eq) :: y
    real(dp) :: t, tmax, dt 


    t = 0
    tmax = 7 * 60 * 60 * 24 * 365
    dt = 30*60*60*24


    y(1:3) = (/ 1.521e+11_dp, 0.0_dp, 0._dp /) !x, y ,z
    y(4:6) = (/ 0._dp, 29.78e+3_dp, 0._dp /) !v_x, v_y, v_z

    do while (t < tmax)

        write(77, *) t, dt

        write(21, *) t, y(1)
        write(22, *) t, y(2)
        write(23, *) t, y(3)


        call rkf45step (t, dt, y)



    end do


end program kepler