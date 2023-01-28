module double_pendulum

    use numtype
    implicit none
    !integer, parameter :: n_eq = 4
    real(dp), parameter :: gravity = 9.81_dp, m1 = 0.2, m2 = 0.4, l1 = 1.2, l2 = 1.5


end module double_pendulum

program main

    use double_pendulum
    use numtype
    implicit none
    real(dp), dimension(n_eq) :: y
    real(dp) :: t, tmax, dt, TE, KE, UE, LG

    t = 0
    tmax = 50
    dt = 0.05

    y(1) = pi/2-0.4
    y(2) = -pi/2-0.5
    y(3) = 10
    y(4) = 0

    do while (t < tmax)

        if (abs(y(1)) > pi) then

            y(2) = -y(2)
        end if

        write(11,*) t, dt
        write(22,*) t, y(1)
        write(33,*) t, y(2)
        write(44,*) t, y(3)
        write(55,*) t, y(4)
        write(66,*) t, mod(y(1),2*pi), l1*sin(y(1)), mod(y(2),2*pi), l1*sin(y(1))+l2*sin(y(2))

        KE = (0.5_dp) * (m1+m2)*(l1**2*y(3)**2) + (0.5_dp)*m2*l2**2*y(4)**2 + m2*l1*l2*y(3)*y(4)*cos(y(1)-y(2))
        UE = -m1*l1*cos(y(1))-m2*gravity*(l1*cos(y(1))+l2*cos(y(2)))
        TE = KE + UE
        LG = KE - UE
        
        write(77,*) t, KE, UE, TE, LG

        call double_rkf45step(t,dt,y)

    end do

end program main
