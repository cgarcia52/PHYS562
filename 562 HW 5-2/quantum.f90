
module qmsetup
    use numtype
    implicit none
    real(dp), parameter :: hsquare = 41.47, mass = 0.5
    integer, parameter :: nsize = 501, lwork = 4*nsize
    real(dp), allocatable :: hammy(:,:), X(:), pot(:), Ene(:), work(:)
end module qmsetup

module identities
    use numtype
    implicit none
    integer, parameter :: maxint = 500, max=5
    !real(dp), external :: q1,q2,q3,q4
    real(dp) :: a, b, c, d, resgl,resgl2,resgl3
    integer :: m, i, ifail, itype
    real(dp), dimension(maxint) :: abscis, weight
end module identities

program quantum
    use numtype
    use identities
    use qmsetup
    implicit none
    integer :: n1, l, pain, n2, info, nn ! where n2 is n1'
    real(dp) :: xmax, xmin, dx
    real(dp), external :: Potential

    b = 2000 ! bounds of the rectangular integration
    a = 0
    c = 0
    d = 0
    itype = 0
    m = 500

    call d01bcf(itype,a,b,c,d,maxint,weight,abscis,ifail) ! program with the types of integration methods
    if( ifail /= 0 ) stop ' ifail '

    resgl = 0
    n1 = 0 ! initial n and n prime values
    n2 = 0
    pain = 3 ! alpha, parameter
    l = 0 ! another variable constant

    do n1 = 0,max ! iterates both n and nprime as well as m (number of points we want)
        do n2 = 1,max
            do i = 1, m ! the three bullet points, the integrals we want to solve for
                resgl=resgl+weight(i)*(Phi(n1,abscis(i) , l, pain))*((1.0/abscis(i)))*(Phi(n2, abscis(i), l, pain))
                resgl2=resgl2+weight(i)*(Phi(n1,abscis(i) , l, pain)) * (Phi(n2, abscis(i), l, pain))
                resgl3 = resgl3+weight(i)*(Phi(n1,abscis(i) , l, pain)) * ((2*pain*(n1+l+1)/abscis(i)-pain**2))&
                    *(Phi(n2, abscis(i), l, pain))
            end do
            write(2, fmt=*) 'integral (Phi*1/r*Phi)',n2,n1,resgl 
            write(3, fmt=*) 'integral (Phi*Phi)',n2,n1,resgl2
            write(4, fmt=*) 'integral Phi*(-d^2/dr^2+l(l+1)/r^2)*phi',n2,n1,resgl3
        end do   
    end do 

    allocate( hammy(nsize,nsize), X(nsize), pot(nsize), Ene(nsize), work(lwork) ) ! distribute the calculated conditions out

    xmax = 6 !range of input values for x, from -6 to 6
    xmin = -xmax
    dx = (xmax - xmin)/(nsize - 1) ! step size equation
    nn = nsize - 1

    !$omp parallel private(i)
    !$omp do
    do i = 1, nn
        X(i) = xmin + (i-1)*dx
        pot(i) = Potential(X(i))
    end do 
    !$omp end do
    !$omp end parallel  

    do i = 1, nn
        write(10, * ) X(i), pot(i) 
    end do    

    hammy(1:nn,1:nn) = 0

    !$omp parallel private(i)
    !$omp do
    do i = 1, nn

        hammy(i,i) = - hsquare**2/(2*mass) * (-2)/dx**2 + pot(i) 
        hammy(i,i+1) = - hsquare**2/(2*mass)  * 1/dx**2
        hammy(i+1,i) = - hsquare**2/(2*mass)  * 1/dx**2

    end do 
    !$omp end do
    !$omp end parallel 

    info = 0
    call dsyev('v','u', nn, hammy, nsize, Ene, work, lwork, info ) ! to get the energy eigenvalues for certain levels of dimensions the matrix is on
    if( info /= 0 ) stop

    print *,Ene(1:25) 
    ! gnuplot requires separate fort files to plot correctly
    ! these represent the bound states for each Ene level of the hamiltonian, which transfers to a plotting tool
    do i = 1, nn
        write(20,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,1)
        write(30,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,2)
        write(40,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,3)
        write(50,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,4)
        write(60,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,5)
        write(70,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,6)
        write(80,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,7)
        write(90,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,8)
        write(100,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,9)
        write(110,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,10)
        write(120,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,11)
        write(130,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,12)
        write(140,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,13)
        write(150,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,14)
        write(160,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,15)
        write(170,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,16)
        write(180,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,17)
        write(190,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,18)
        write(200,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,19)
        write(210,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,20)
        write(220,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,21)
        write(230,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,22)
        write(240,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,23)
        write(250,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,24)
        write(260,fmt=*)  X(i), ( 1/sqrt(dx) )* hammy(i,25)
    end do
      
    contains 

        recursive function fact(n1) result(res)
            ! factorial equation
            implicit none
            integer, intent(in) :: n1
            real(dp) :: res 

            if( n1 < 0 ) then
                stop ' wrong '
            else if ( n1 == 0 ) then
                res = 1._dp
            else 
                res = n1 * fact(n1 - 1)
            end if

        end function fact

        ! here Phi is represented in a Laguerre recursive function
        recursive function Phi(n1, r, l, pain) result(res)
            use numtype
            implicit none
            integer, intent(in) :: n1,l,pain
            real(dp) :: Y, brain, res
            real(dp), intent(inout) :: r
        
            Y = n1 + 1._dp

            brain = 2*l+1._dp

            if(n1 < 0) then ! three integral equations written out
                stop 'n1 greater than 0'
            else if(n1 == 0) then
                res = 1
                res = sqrt(gamma(Y)/gamma(Y + 2*l + 1._dp)) * exp(-pain * r)* (2 * pain * r)**(l + 1._dp) * res
            else if(n1 == 1) then 
                res = 1 + brain - 2 * pain * r
                res = sqrt(gamma(Y)/gamma(Y + 2*l + 1._dp)) * exp(-pain * r)* (2 * pain * r)**(l + 1._dp)*res
            else 
                res = ( (2*(n1-1) + 1 + brain - 2*pain*r)*(Phi(n1-1,r,l,pain)*(sqrt(gamma(Y)/gamma(Y+2*l+1._dp)) * exp(-pain*r) &
                * (2 * pain * r)**(l + 1._dp)))-(n1-1+pain)*(Phi(n1-2,r,l,pain)&
                *(sqrt(gamma(Y)/gamma(Y+2*l+1._dp))* exp(-pain*r)&
                 * (2 * pain * r)**(l + 1._dp))))/(n1)

                res = sqrt(gamma(Y)/gamma(Y+2*l+1._dp)) * exp(-pain*r) * (2 * pain * r)**(l + 1._dp) * res
            end if
       
        end function Phi

end program quantum

function Potential(x)
    ! potential function for V(x) for hamiltonian part 
    use numtype 
    implicit none
    real(dp) :: x, Potential

    Potential = (1438.72 * exp(-3.11 * x))/abs(x + iic) - (626.885 * exp(-1.55 * x))/abs(x + iic)
    ! v1 = 1438.72
    ! B1 = 3.11 
    ! x = r
    ! v2 = -626.885
    ! B2 = 1.55
end function Potential


