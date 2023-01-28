module gamma_setup

    use numtype
    implicit none

    integer, parameter :: npmax = 1200, npar = 8 
    integer, parameter :: nspmin = 230, nspmax = 900
    integer :: yy(0:npmax), icall, iprint, nsp

end module gamma_setup


program gammafit

    use gamma_setup
    implicit none 
    real(dp), external :: chi2 
    integer :: i, itmin, itmax, stat
    real(dp) :: xstart(1:npar), epsf, stepi, fstart

    open(unit=2, file='na22.data')
    i = 0
    do
        read(2, *, iostat=stat) yy(i)
        if (stat /= 0) exit
        i = i + 1
    end do
    nsp = i - 1
    close(2)

    xstart(1:npar) = (/-0.01, 100.0, 3000.0, 320.0, 25.0, &
                        200.0, 720.0, 25.0 /)
    icall = 0
    iprint = 7
    fstart = chi2(xstart)
    stepi = 0.05_dp

    epsf = 0.001_dp
    itmin = 100
    itmax = 2000
    iprint = 0

    call downhill(npar, chi2, xstart, fstart, stepi, epsf, itmin, itmax)
    iprint = 17
    fstart = chi2(xstart)

end program gammafit


function chi2(par) result(s2)

    use gamma_setup
    implicit none
    integer :: i
    real(dp) :: fi, s2, a, b, y1, x1, sig1, y2, x2, sig2, par(npar)

    icall = icall + 1

    a= par(1); b = par(2)
    y1 = par(3); x1 = par(4); sig1 = par(5)
    y2 = par(6); x2 = par(7); sig2 = par(8)
    s2 = 0

    do i = nspmin, nspmax

        fi = a*i + +y1*exp(-(i-x1)**2/sig1**2) &
                + y2**exp(-(i-x2)**2/sig2**2)
        s2 = s2 + (yy(i)-fi)**2 * 1/sqrt(yy(i)+2._dp)

    end do

    s2 = s2 / abs(nspmax-nspmin)
    print '(i4,2x,8f12.4,f20.4)', icall, par(1:npar), s2

    if (iprint /=0) then

        do i = nspmin, nspmax
            fi = a*i + b + y1*exp(-(i-x1)**2/sig1**2) &
                + y2*exp(-(i-x2)**2 / sig2**2)
            write(unit = iprint, fmt = '(i4, i10, f15.2 )') i, yy(i), fi

        end do

    end if

end function chi2