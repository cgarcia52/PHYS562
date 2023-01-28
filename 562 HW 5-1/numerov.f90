module setup
    use numtype 
    implicit none
    real(dp), parameter :: hbar = 1, m = 1, xmin = -100, xmax = 100
    integer, parameter :: dim = 500
end module setup

program Numerov
    use setup
    implicit none

    real(dp), dimension(0:dim) ::  Ene, work, wi, wr, x
    real(dp), dimension(4*(dim+1)) :: work2
    real(dp), dimension(0:dim,0:dim) :: One, Two, matrix, V , Ham, L, R
    
    real(dp), external :: Potential, g
    real(dp) :: dx , part, v0, v1
    
    integer, dimension(0:dim) :: IPIV
    integer :: i, j, info, lwork, lwork2, bound, scatter, zero
    
    print *, 'What Part, 1 or 2?, v0, v1.' 
    ! part a = 1, part b = 2, v0 and v1 are constants, only v1 in part a
    read *, part,v0,v1
    
    dx = (xmax-xmin)/dim ! change for interations in x input values

    do i = 0, dim ! 500 point iterations, calls X into potential equation
        x(i) = xmin + i*dx
        write(20,*) x(i), Potential(part,x(i),v0,v1) 
    end do

    do i = 0, dim   !construct constant matrix A and B constructed by element number row by columb (ixj)
        do j = 0, dim
            if (i+1==j) then
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

    print'(10f7.4)', One(0:9,0:9) ! matrix formation check
    print*, '----------------------------'

    Two = (1.0/12.0)*Two ! checks for matrix completion inverse
    print'(10f7.4)', Two(0:9,0:9)
    print*, '----------------------------' 
    
    info = 0
    lwork = dim
    call dgetrf(dim+1, dim+1, Two, dim+1, IPIV, info)             
    ! Decompostion for Lower Upper triangulization
    call dgetri(dim+1, Two, dim+1, IPIV, work+1, lwork+1, info)   
    ! Inverse for Matrix called from Lapack
        if(info /= 0) then
            !print*, info
            stop 'info =/= 0'
            
        end if

    print'(10f7.4)', Two(0:9 ,0:9) ! matrix of 2, 10 by 10 of B, normal B matrix
    print*, '----------------------------' 

    do i = 0, 500
        x(i) = xmin + i*dx
    end do

    do i = 0, dim   !construct V (Potential matrix)
        do j = 0, dim
            if (i==j) then ! diagonal of matrix for eigenvalues
                V(i,j) = Potential(part, x(i), v0, v1)
            else
                V(i,j) = 0.0 !elsewhere in matrix is 0
            end if
        end do
    end do

    matrix = matmul(Two, One) ! multiply matrix
    Ham = -( (hbar)**2 )/(2._dp*m) * matrix + V
    
    print'(10f7.4)',Ham(0:9,0:9) ! checks schrodinger equation
    print*, '----------------------------' 

    info = 0
    lwork = 4 * ( dim + 1 )

    call dgeev('v0','v1', dim+1, Ham, dim+1, wr, wi, L, dim+1, R, dim+1, work2, lwork, info)  
        !Eigenvalues
        if(info /= 0) then
            stop 'info =/= 0'
        end if

    print'(10f7.4)',Ham(0:9,0:9) 
    ! matrix  eigenvalues through lapack function

    print*, '---------------------'

    print'(10f7.4)',Ham(dim-9:dim , dim-9:dim)
    bound = 0 ; scatter = 0; zero = 0 

    do i = 0,dim ! iterates bounds and scatter states higher by 1
        Ene(i) = Ham(i,i)
        if (Ene(i) < 0) then ! negative values means bounded state, positive are scattering states
            bound = bound + 1 
        else if (Ene(i) > 0) then
            scatter = scatter + 1
        else 
            zero = zero + 1
        end if
    end do
    
    print*, '---------------------'
    print *, ' scattering state count: ',scatter,' | ', 'bound state count: ', bound,' | ',' E=V: ',zero

end program Numerov 

function Potential(part,x,v0,v1)
    use setup
    implicit none
    real(dp) ::  Potential, x , part, v0, v1
    
    if (part == 1) then ! part a Potential
        Potential = -10 * (1/cosh(x))**2 + v1 * ( tanh(x) )
    else if (part == 2) then  ! part b Potential
        Potential = -v0 * (1/cosh(x))**2 + v1 * (1/cosh(10 * x))**2
    end if

end function Potential