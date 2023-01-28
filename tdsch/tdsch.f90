
module qm_time

    use numtype
    implicit none
    real(dp), parameter :: hbar = 1, hbar2 = hbar**2, mass = 1
    real(dp), parameter :: v0 = 1, al = 2, x0 = 4

    integer, parameter :: nd = 1000 



end module qm_time


program qm1d

    use qm_time
    implicit none
    real(dp) :: t, dt, tmax
    real(dp) :: xmin, xmax, dx, xm, km, sigma, dp0, energy
    real(dp), dimension(0:nd) :: xi, vxi
    complex(dp), dimension(0:nd) :: psi, phi, dd, du, dl 
    real(dp), external :: potential
    complex(dp), external :: psi0 
    integer :: i, info, nn


    xmax = 100
    xmin = -xmax
    dx = (xmax-xmin)/nd

    xm = -60
    km = 2
    sigma = 2
    dp0 = hbar/(sigma*sq2) 
    energy = hbar2* km**2/(2*mass) + dp0**2/(2*mass) 
    print *,' energy: ', energy

    do i = 0, nd
        xi(i) = xmin + i*dx
        vxi(i) = potential(xi(i))
        psi(i) = psi0(xi(i),xm,km,sigma)
        write(17,*) xi(i), vxi(i)
        write(19,*) xi(i), abs(psi(i))**2
    end do

    tmax = 50
    dt = 0.1
    nn = 0
 
    do while ( t < tmax )

        if( mod(nn,50) == 0 ) then
            do i = 0, nd
                write(100+nn,*) xi(i), abs(psi(i))**2
            end do
        end if

        dd(0:nd) = 0.5_dp + ic*hbar* dt/(4*mass*dx**2) +ic*vxi(i)*dt/(4*hbar)
        du(0:nd) = -ic*hbar*dt/(8*mass*dx**2)
        dl(0:nd) = -ic*hbar*dt/(8*mass*dx**2)
        phi(0:nd) = psi(0:nd)

        call zgtsv(nd+1, 1, dl, dd, du, phi, nd+1, info  )

        psi(0:nd) = phi(0:nd) - psi(0:nd)

        nn = nn + 1
        t = t + dt

    end do 


end program qm1d


function psi0(x,xm,km,sigma)

    use numtype
    implicit none
    real(dp) :: x, xm, km, sigma
    complex(dp) :: psi0 

    psi0 = 1/sqrt(sigma*sqpi) * exp(-(x-xm)**2/(2*sigma**2)) * exp(ic*km*x)

end function psi0


function potential(x)

    use qm_time
    implicit none
    real(dp) :: x, potential

    potential = v0/4 *(1-tanh(al*(x-x0))) * (1+tanh(al*(x+x0)))

end function potential




