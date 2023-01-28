module mparameters

    use numtype
    implicit none
    real(dp), parameter :: m1 = 4._dp, m2 = 3._dp, l1 = 1._dp, l2 = 1.5_dp, g = 9.81_dp

end module mparameters

module dynamics


end module dynamics 

program double_pendulum

    use numtype
    use mparameters
    implicit none
    real(dp), dimension(n_eq) :: y
    real(dp) :: t, tmax. dt

    t = 0
    tmax = 20
    dt = 0.1_dp

    do while (t < tmax)

        write(12,*) t, 
        write(13,*) t, 
        write(14,*) t,
        write(15,*) t,

        call rkf45step

    end do

end program double_pendulum

