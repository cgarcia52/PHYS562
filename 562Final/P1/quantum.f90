module setup

    use numtype
    implicit none
    real(dp), parameter :: hbar = 1, m = 1, xmin = -25, xmax = 25
    integer, parameter :: dim = 100

end module setup

program quantum
    
    use setup
    implicit none
    real(dp), dimension(0:dim) :: erg, work, wi, wr, x
    real(dp), dimension(4*(dim+1)) :: work2
    real(dp), dimension(0:dim,0:dim) :: One, Two, matrix, V, Haml, L, R

    real(dp), external :: Potential, g
    real(dp) :: dx, lambda

    integer, dimension(0:dim) :: ipiv
    integer :: i, j, info, lwork, lwork2, bound, scatter, zero

    print*, 'What value for lambda?'
    read*, lambda


    do i = 0, dim
        x(i) = xmin + i*dx
        write(20,*) x(i), Potential(lambda, x(i))
    end do

    do i = 0, dim
        do j = 0, dim
            if(i+1==j) then
                One(i,j) = 1.0
                Two(i,j) = 1.0
            else if (i==j+1) then
                One(i,j) = 1.0
                Two(i,j) = 1.0
            else if (i==j) then
                One(i,j) = -2.0
                Two(i,j) = 10.0
            else
                One(i,j) = 0.0
                Two(i,j) = 0.0
            end if
        end do
    end do

    print'(10f7.4)', One(0:9,0:9)
    print*, '------------------------'

    Two = (1.0/12.0) * Two
    print '(10f7.4)', Two(0:9,0:9)
    print*, '-------------------------'

    info = 0
    lwork = dim
    call dgetrf(dim+1, dim+1, Two, dim+1, ipiv, info)

    call dgetri(dim+1, Two, dim+1, ipiv, work+1, lwork+1, info)

        if(info /= 0) then

            stop 'info =/= 0'
        
        end if

    print'(10f7.4)', Two(0:9,0:9)
    print*, '-----------------------------'

    do i = 0, 100
        x(i) = xmin + i*dx
    end do

    do i = 0, dim
        do j = 0, dim
            if (i==j) then
                V(i,j) = Potential(lambda, x(i))
            else
                V(i,j) = 0.0
            end if
        end do
    end do

    matrix = matmul(Two, One)
    Haml = -( (hbar**2))/(2._dp*m) * matrix + V

    print'(10f7.4)', Haml(0:9,0:9)
    print*, '----------------------------'

    info = 0
    lwork = 4 * (dim + 1)

    call dgeev(dim+1, Haml, dim+1, wr, wi, L, dim+1, R, dim+1, work2, lwork, info)

        if(info /=0) then
            stop 'info =/= 0'
        end if

    print '(10f7.4)', Haml(0:9,0:9)
    print*, '---------------------------'

    print '(10f7.4)', Haml(dim-9:dim, dim-9:dim)
    bound = 0; scatter = 0; zero = 0

      do i = 0,dim
        erg(i) = Haml(i,i)
        if(erg(i) < 0) then
        bound = bound + 1
        else if (erg(i) > 0) then
            scatter = scatter + 1
        else
            zero = zero + 1
        end if
    end do

    print*, '-------------------'
    print*, 'scattering state count: ', scatter, ' | ', 'bound state', bound, ' | ', ' E=V: ', zero

end program quantum

function Potential(lambda,x)

    use setup
    implicit none
    real(dp) :: Potential, lambda, x

    if (lambda == 1) then
        Potential = (-1/cosh(x)**2)
    else if (lambda == 2) then
        Potential = (-3/cosh(x)**2)
    else if (lambda == 5) then
        Potential = (-15/cosh(x)**2)
    end if

end function Potential

