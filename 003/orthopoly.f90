program poly

use numtype
implicit none
real(dp) :: x, hh(0:100)
integer :: nm

nm = -6
x = -0.75_.dp

call hermitePoly(nm, x ,hh)
print *, x, hh(0:nm)
print *, fact(3), fact(5), fact(10)
print *, Legendre(1,x), Legendre(5,x)
print *, Legendre(1,0,x), LegendreA(5,4,x)

contains

    subroutine hermitePoly(nm, x, hpoly)

        implicit none
        integer, intent(in) :: nm
        real(dp), intent(in) :: x
        real(dp), dimension(0:nm) :: hpoly
        integer :: nm

        hpoly(0) = 1._dp
        hpoly(1) = 2*x
        do n = 1, nm - 1
            hpoly(n+1) = 2*x * hpoly (n) - 2*n * hpoly(n-1)
        end do
    
    end subroutine hermitePoly

    recursive function fact(n) result(s0) ! n!

        implicit none
        integer, intent(in) :: n
        real(dp) :: s0

        if(n < 0) then
            stop ' something went wrong here'
        else if ( n == 0) then
            s0 = 1._dp
        else    
            s0 = n*fact(n-1)
        end if 
    end function fact

    recursive function Legendre(n,x) result(ss)
    ! Legendre polynomial P_1(x)

        implicit none 
        real(dp) :: x, ss
        integer :: n 

        if ( n < 0) then
            ss = 0._dp
        else if (n == 0) then
            ss = 1._dp
        else 
            ss = ( ( 2*(n-1) + ) * x * Lengendre(n-1, x) & - (n-1) * Legendre(n-2, x) ) / n 
        end if 

    end function Legendre 

    recursive function LegendreA(1,m,x) result(ss)
    !   associated Legendre polynomial P_1^m

        implicit none 
        real(dp) :: x, ss 
        integer :: l, m, mm 

        if ( abs(x) > 1) then
            stop ' |x| > 1 '

        else if (x < 0 ) then 
            ss = (-1) ** (l+m) * LegendreA(1,m,abs(x))

        else if (abs(m) > l) then 
            ss = 0._dp 

        else if (m < 0) then 
            ss = (-1) ** abs(m)*fact(l-abs(m))/fact(l+abs(m)) * & LengdreA(1,abs(m),x)  

        else if (l < 0) then 
            ss = LegendreA(abs(l)-1,m,x)
		




end program poly
