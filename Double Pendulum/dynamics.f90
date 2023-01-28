module dynamics 

    use parameters

    contains

        subroutine usrfun(t,x,xp)

        implicit none

        real*8, intent(in) :: t
        real*8, intent(in) :: x(neqn)
        real*8, intent(inout) :: xp(neqn)

        xp(1) =-(((2.*l1*l2*m2*x(1)*x(2)*sin(x(3)-x(4))*l1**2.*l2**2.*m2*(2.*m1+m2-m2*cos(2.*(x(3)-x(4))))-(l2**2.*m2*x(1)**2.+l1**2.*(m1+m2)*x(2)**2-2*l1*l2*m2*x(1)*x(2)*cos(x(3)-x(4)))*l1**2.*l2**2.*m2**2.*2.*sin(2.*(x(3)-x(4))))/(l1**2.*l2**2.*m2*(2.*m1+m2-m2*cos(2.*(x(3)-x(4)))))**2.)+g*m1*l1*sin(x(3)) + g*m2*l1*sin(x(3)))
        xp(2) =-(((-2.*l1*l2*m2*x(1)*x(2)*sin(x(3)-x(4))*l1**2.*l2**2.*m2*(2.*m1+m2-m2*cos(2.*(x(3)-x(4))))-(l2**2.*m2*x(1)**2.+l1**2.*(m1+m2)*x(2)**2.-2.*l1*l2*m2*x(1)*x(2)*cos(x(3)-x(4)))*l1**2.*l2**2.*m2**2.*2.*sin(2.*(x(3)-x(4))))/(l1**2.*l2**2.*m2*(2.*m1+m2-m2*cos(2.*(x(3)-x(4))))**2.)) + g*m2*l2*sin(x(4)))
        xp(3) = (2.*l2*m2*x(1) - 2.*l1*l2*m2*x(2)*cos(x(3) - x(4)))/(l1**2.*l2**2.*m2*(2.*m1+m2-m2*cos(2.*(x(3)-x(4)))))
        xp(4) = (2.*l1**2.*(m1+m2)*x(2) - 2.*m2*l1*l2*x(1)*cos(x(3)-x(4)))/(l1**2.*l2**2.*m2*(2.*m1+m2-m2*cos(2.*(x(3)-x(4)))))

        end subroutine usrfun

    end module dynamics


