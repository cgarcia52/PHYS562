program integrate
    implicit none
    integer,parameter :: cp = selected_real_kind(14)
    integer,parameter :: N = 1000
    real(cp),dimension(N) :: f,xc
    real(cp),dimension(N+1) :: x
    real(cp) :: s,xmax,xmin,dx
    integer :: i
    xmin = 0.0_cp
    xmax = 1.0_cp
    dx = (xmax - xmin)/real(N,cp)
    x = (/(xmin + dx*(i-1),i=1,N+1)/)
    ! Define x at center
    do i=1,N
    xc(i) = x(i) + 0.5_cp*dx
    enddo
    ! Define f
    do i=1,N
    f(i) = 1.5 * sqrt(xc(i))
    enddo
    ! Integrate (Midpoint method)
    s = 0.0_cp
    do i=1,N
    s = s + f(i)*dx
    enddo
    write(*,*) 'sum = ',s

end program integrate