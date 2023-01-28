program matrix

    use numtype
    implicit none

    complex(dp), parameter :: & 
        one = (1._dp, 0._dp) , two = (0._dp, 0._dp) , three = (0._dp, 1._dp) ! the set is (real, imaginary)
    integer, parameter :: ndim = 10, lwork = 2*ndim - 1 ! 20, 50, 100 N by N matrix truncated

    integer, parameter :: m = 1, t = 1 ! 1, 2, 3 for w2 and 1 2 for t

    real(dp), dimension(ndim) :: w

    real(dp) :: w2 = 1 ! 1, 2, 3 

    real(dp), dimension(lwork) :: rwork = 3*ndim -2 !for the zheev function, used to calculate eigenvalues
    complex(dp), dimension(lwork) :: work  

    complex(dp), dimension( ndim, ndim) :: d, e, g, f, q, r

    complex(dp), dimension( ndim, ndim) :: X, P, H ! truncating for N by N matrix for X and P matrices, and H is the hamiltonian 

    integer :: a, b, c, n, info, ipiv, LDA ! the position numbers for the matrix elements

!X operator matrix
do a = 1, ndim 
    do b = 1, ndim ! 4 separate checks running thru 2 conditions, a + b
        if (a == b) then
            X(a, b) = two
        else if (abs( a - b) == 1) then
            X(a, b) = min(a,b)
        else
            X(a, b) = two
        end if
            X(a, b) = sqrt(0.5) * sqrt( X(a,b) )
    end do
end do

!P operator matrix
do a = 1, ndim ! 5 separate checks running thru 2 conditions, a + b
    do b = 1, ndim
        if (a == b) then
            P(a, b) = two
        else if ( ( abs( a - b) == 1) .and. (a < b) ) then
            P(a, b) = min(a,b)
            P(a, b) = - sqrt( P(a, b))
        else if ( (abs( a - b) == 1) .and. (a > b) ) then
            P(a, b) = min(a,b)
            P(a, b) = sqrt( P(a, b))
        else
            P(a, b) = two
        end if
            P(a, b) = sqrt(0.5) * ( P(a,b) ) * three
    end do
end do

do n = 1, ndim
    !print '(20 f7 .2)', X(n, 1:ndim) ! to see the matrix in terminal, but its messy
    write(2,*) X(n, 1:ndim) ! prints out in a fort file with easier readability
end do

!print *, '-------------------'

do n = 1, ndim 
    !print '(20 f7 .2)', P(n, 1:ndim) ! to see the matrix in terminal, but its messy
    write(3,*) P(n, 1:ndim) ! prints out in a fort file with easier readability
end do


!print *, '-------------------' ! Hamiltonian equation setup

H = (0.5) * (matmul(P,P) / (m) ) + (0.5) * m * w2 * matmul(X,X) ! matmul takes the X and P and squares it in matrix form

do n = 1, ndim
    !print '(20 f7 .2)', H(n, 1:ndim) ! to see the matrix in terminal, but its messy
    write(1,*) H(1:ndim , n) ! prints out in a fort file with easier readability
end do



!print *, '-------------------' ! part A problem 1, eigenvalues


!print *, ' eigenvalue '

e(1:ndim,1:ndim) = H(1:ndim,1:ndim) ! setting function H equal to arbitrary function

info = 0

call zheev('v','u', ndim, e, ndim, w, work, lwork, rwork, info  ) ! lapack

!print *, info

if( info /= 0 ) stop ' info /= 0 '

!write(9,*) w(1:ndim)

do a = 1, ndim
    !print '(f10.4)',w(a)
    write(10,*) w(a)
    print '(20f7.2)', e(1:ndim,a) !Hamiltonian
    !write(11,*) e(1:ndim,a)
end do

print *, '-------------------' ! part b problem 1


print *,' orthogonality  ' ! finds if eigenvectors are orthonormal

g(1:ndim,1:ndim) = matmul( transpose( conjg( e(1:ndim,1:ndim) ) ), e(1:ndim,1:ndim) )


do a =1, ndim

    print '(20f7.2)', g(1:ndim,a)

    do b = 1, ndim
        write(12,*) b,a, g(b,a), dot_product( e(1:ndim,b)  ,  e(1:ndim ,a) )
    end do
end do

print *, '-------------------'

print *,' completeness  ' ! finds if the matrix forms a complete set

g(1:ndim,1:ndim) = matmul(  e(1:ndim,1:ndim) , transpose(conjg(e(1:ndim, 1:ndim))) )

do a = 1, ndim
    print '(20f7.2)', g(1:ndim,a) 
    write(13,*) b, a, g(1:ndim,a) 
end do


print *, '-------------------' ! part c problem 1

! matrix representation of exp(-iHt)

forall ( a = 1:ndim, b = 1:ndim ) f(a,b) = e(a,b) * exp(- three * w(a) * t) 

g(1:ndim,1:ndim) = matmul(  f(1:ndim,1:ndim) , transpose(conjg(e(1:ndim,1:ndim))) )

!call zgetrf( ndim, ndim, g, ndim, ipiv, info  ) ! 4th variable is Leading dimension of array

!if(info /= 0 ) stop ' info zgetrf /= 0 '



    do b = 1, ndim
        print '(20f7.2)', g(b,1:ndim)
        !write(22,*) g(1:ndim, 1:ndim)
    end do


!q(1:ndim,1:ndim) = matmul( g(1:ndim,1:ndim) , transpose(conjg(g(1:ndim,1:ndim))) )

!do a =1, ndim
    !do b = 1, ndim
        !print *,b,a, q(b,a)  
        !write(23,*) q(1:ndim, 1:ndim)
    !end do
!end do



end program matrix