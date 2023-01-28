program Herm1

    use numtype
    implicit none
    complex(dp), parameter :: one=(1._dp), null=(0._dp)
    integer, parameter :: ndim = 5, lwork = ndim**2
    complex(dp), dimension(ndim,ndim) :: a, b, c, d, e, f, g
    integer :: nn, i, info, j, ipiv(ndim)
    complex(dp) :: ss1, ss2, work(lwork)
    real(dp) :: w(ndim), rwork(3*ndim)


    nn = 4
    a(1:nn, 1:nn) = reshape( (/ 2*one, -one, null, -one, &
     2*one, -one, null, -one, &
     -one, null, 2*one, null, &
      null, -2*one, -one, null/), (/nn,nn/) )
    b(1:nn, 1:nn) = reshape( (/ -one, null, 2*one, null, &
     2*one, null, -one, null, &
      null, 2*one, -one, null, &
      null, -2*one, null, -one/), (/nn,nn/) )


    do i = 1, nn
        print '(10f10.3)',a(i,1:nn)
    end do
    print *,'-----------------------------'
    do i = 1, nn
        print '(10f10.3)',b(i,1:nn)
    end do

    print *,' Hermiticity '
    c(1:nn,1:nn) = transpose(conjg(a(1:nn,1:nn)))
    do i = 1, nn
        print '(10f10.3)',c(i,1:nn) - a(i,1:nn)
    end do
    d(1:nn,1:nn) = transpose(conjg(b(1:nn,1:nn)))
    do i = 1, nn
        print '(10f10.3)',d(i,1:nn) - b(i,1:nn)
    end do

    print *,' [A,B] = ? '
    c(1:nn,1:nn) = matmul( a(1:nn,1:nn) , b(1:nn,1:nn) )
    d(1:nn,1:nn) = matmul( b(1:nn,1:nn) , a(1:nn,1:nn) )
    do i = 1, nn
        print '(10f10.3)',c(i,1:nn) - d(i,1:nn) 
    end do
    print *,'-----' 
    ss1 = null;  ss2 = null
    do i = 1, nn
        ss1 = ss1 + c(i,i)
        ss2 = ss2 + d(i,i)
    end do
    print *,' Tr(AB) & Tr(BA) ', ss1, ss2

    print *, sum( (/ (c(i,i),i = 1,size(c,1)) /) ), sum( (/ (d(i,i),i = 1,size(d,1)) /) )


    print *,' Eigenvalues & Eigenvectors '

    e(1:nn,1:nn) = a(1:nn,1:nn)
    info = 0
    call zheev('v','u', nn, e, ndim, w, work, lwork, rwork, info )
    if ( info .ne. 0 ) stop ' info .ne. 0 '
    print *,' eigenvalues ', w(1:nn)
    do i = 1, nn
        print *,i, w(i)
        print '(5(2f10.4,2x))', e(1:nn,i)
    end do

    print *,' orthogonality '
    d(1:nn,1:nn) = matmul( conjg(transpose(e(1:nn,1:nn))) , e(1:nn,1:nn) )

    do i = 1, nn
        do j = 1, nn
            print *, i, j, dot_product( e(1:nn,i) , e(1:nn,j) ), d(i,j)
        end do 
    end do 

    print *,' completeness $\sum_n | n >  < n | $ '
    d(1:nn,1:nn) = matmul( e(1:nn,1:nn) ,  conjg(transpose(e(1:nn,1:nn)))  )

    do i = 1, nn
        do j = 1, nn
            print *, i, j,  d(i,j)
        end do 
    end do 


    print *,' spectral resolution $\sum_n  | n > w_n < n | $ '

    forall (i=1:nn, j=1:nn) f(i,j) = e(i,j)*w(j)

    d(1:nn,1:nn) = matmul( f(1:nn,1:nn) ,  conjg(transpose(e(1:nn,1:nn))) )

    do i = 1, nn
        do j = 1, nn
            print *, i, j,  d(i,j)
        end do 
    end do 

    print *,' Unitarity of eigenvectors '

    print *,'  U  '
    do i = 1, nn
        print '(10f10.3)', e(i,1:nn)
    end do

    c(1:nn,1:nn) =  transpose(conjg(e(1:nn,1:nn)))
    print *,'  U^+  '
    do i = 1, nn
        print '(10f10.3)', c(i,1:nn)
    end do

    f(1:nn,1:nn) = matmul( e(1:nn,1:nn) , c(1:nn,1:nn) )
    print *,' U U^+  '
    do i = 1, nn
        print '(10f10.3)', f(i,1:nn)
    end do

    print *,'  U^+ U  '
    f(1:nn,1:nn) = matmul( c(1:nn,1:nn) , e(1:nn,1:nn) )
    do i = 1, nn
        print '(10f10.3)', f(i,1:nn)
    end do

    print *,'  U^+ A U  '
    f(1:nn,1:nn) = matmul( a(1:nn,1:nn) , e(1:nn,1:nn) )
    d(1:nn,1:nn) = matmul( c(1:nn,1:nn) , f(1:nn,1:nn) )
    do i = 1, nn
        print '(10f10.3)', d(i,1:nn)
    end do

    print *,'---------- U inverse --------------------------'
    f(1:nn,1:nn) =  e(1:nn,1:nn)

    info = 0
    call zgetrf( nn, nn, f, ndim, ipiv, info  )
    if(info /= 0 ) stop ' info \= 0 in zgetrf '
    call zgetri(nn, f, ndim, ipiv, work, lwork ,info )
    if(info /= 0 ) stop ' info \= 0 in zgetri '
     do i = 1, nn
        print '(10f10.3)', f(i,1:nn)
        print '(10f10.3)', c(i,1:nn)
    end do   



    print *,'--------- A inverse ------------------------'
    f(1:nn,1:nn) =  a(1:nn,1:nn)

    info = 0
    call zgetrf( nn, nn, f, ndim, ipiv, info  )
    if(info /= 0 ) stop ' info \= 0 in zgetrf '
    call zgetri(nn, f, ndim, ipiv, work, lwork ,info )
    if(info /= 0 ) stop ' info \= 0 in zgetri '





    print *,' spectral resolution A^-1 =  | n > (w_n)^-1 < n |  '

    forall (i=1:nn, j=1:nn) g(i,j) = e(i,j) * 1/w(j)

    d(1:nn,1:nn) = matmul( g(1:nn,1:nn) ,  conjg(transpose(e(1:nn,1:nn))) )


    do i = 1, nn
        print '(10f10.3)', f(i,1:nn)
        print '(10f10.3)', d(i,1:nn)
    end do   


end program Herm1 