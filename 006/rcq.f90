
module rcq_setup

    use numtype
    implicit none
    integer, parameter :: n_eq = 2

    real(dp), parameter :: emf = 10._dp, capacity = 2._dp, & 
        resistance = 3._dp 


end module rcq_setup


program rcq

   use rcq_setup
   implicit none
   
   real(dp) :: t, dt, tmax
   real(dp), dimension(n_eq) :: y


   t = 0._dp
   tmax = 30._dp 
   dt = 0.1_dp

   !open( 17, file='f17.r')
   !open(18, file='f18.r')

   qq(1) = 0
   qq(2) = emf/resistance

   do while (t < tmax)

    write(*, *) t, qq(1)
    write(17, *) t, qq(1)
    write(18, *) t, qq(2)
    call rk4step(t, dt, qq)


   end do






end program rcq 


subroutine rk4step (x, h, y)

    use rcq_setup
    implicit none

    real (dp), intent(in) :: h
    real(dp), intent(inout) :: x
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy

    interface 

        function kv (t,dt,y)

            use rcq_setup 
            implicit none
            real(dp), dimension(n_eq), intent(in) :: t, dt
            real(dp), dimension(n_eq), intent(inout) :: y
            real(dp) :: x, h
            real(dp) :: kv



        end function kv


    end interface

    k1 = kv(x, h, y)
    k2 = kv(x+h/2, h, y+k1/2)
    k3 = kv(x+h/2, h, y+k2/2)
    k4 = kv(x+h/2, h, y+k3)

    dy = (k1 +2*k2 + 2*k3 +k4)/6

    y = y + dy
    x = x + h

end subroutine rk4step



function kv (t, dt, y) result(k)
    use rcq_setup
    implicit none
    real(dp), intent(in) :: t, dt
    real(dp), dimension(n_eq), intent(in) :: y
    real(dp), dimension(n_eq) :: k, f 

    !--------------------------------------

    f(1) = (emf  - y(1)/capacity)/resistance
    f(2) = y(2)/(capacity/resistance)

    !--------------------------------------

    k = dt * f 

end function kv 





