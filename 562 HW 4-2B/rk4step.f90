subroutine rk4step(x,h,y)

    use numtype

    
    implicit none

    integer, parameter :: n_eq = 2

    real(dp), intent(inout) :: x
    real(dp), intent(in) :: h
    real(dp), intent(inout), dimension(n_eq) :: y
    real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy 
    interface

        function kv ( x, h, y)

            use numtype

            implicit none

            integer, parameter :: n_eq = 2

            real(dp), intent(in) :: x
            real(dp), intent(in) :: h
            real(dp), intent(in), dimension(n_eq) :: y
            real(dp), dimension(n_eq) :: kv

        end function kv

    end interface

    k1 = kv ( x, h, y)
    k2 = kv ( x + h/2, h, y + k1/2)
    k3 = kv ( x + h/2, h, y + k2/2)
    k4 = kv ( x + h ,  h, y + k3  )

    dy = (k1 + 2*k2 + 2*k3 + k4)/6

    y = y + dy
    x = x + h

end subroutine rk4step