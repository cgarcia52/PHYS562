
module taylor_cf_approx

    use numtype
    implicit none
    integer, parameter, public :: max_taylor = 100, ld = max_taylor+1
    real(dp), parameter, private :: eps0 = 1.e-30_dp, eps = 1.e-15_dp 
    real(dp), dimension(0:max_taylor) :: t_coef, cf_coef, c_pade_coef, d_pade_coef   
    
    contains
    
        function taylor_cf_sum(f,n) result(ww)
        !   taylor sum by continued fraction
            
            implicit none
            integer :: n
            real(dp) :: ww, x
            real(dp), dimension(0:max_taylor) :: f, e
            
            call taylor_cfrac(f,n,e)
            x = 1._dp
            ww  = evalcf(e,n,x)
        
        end function taylor_cf_sum
    
        subroutine taylor_cfrac(f,n,e) 
        !   f(n): Taylor coefficients
        !   e(e): continued fraction coeffients
        
            implicit none
            real(dp), dimension(0:max_taylor) :: f, e, r, s
            integer :: n, i, nl, j, k
        
            e(0) = f(0)
            e(1) = f(1)
            e(2) = f(2)/(f(1)+eps0)
            do i = 3, n
                e(i) = f(i)/(f(i-1)+eps0)
                r(i) = e(i)-e(i-1)
            end do
            nl = 1
            do j = 4, n
                if(nl == 1) then
                    do k = j, n
                        s(k) = r(k)/r(k-1)*e(k-1)
                    end do
                else
                    do k = j, n
                        s(k) = r(k)-r(k-1)+e(k-1)
                    end do
                end if
                e(j-1:n) = r(j-1:n)
                r(j:n) = s(j:n)
                nl = -nl
            end do
                
        end subroutine taylor_cfrac
        
        function evalcf(e,n,x) result(g1)
        ! evaluation of continued fraction  
        
            implicit none
            real(dp), dimension(0:max_taylor) :: e
            integer :: n, i
            real(dp) :: x, g0, g1, c0, c1, d0, d1, delta
                   
            g0 = e(0) + eps0 
            c1 = 1+e(1)*x/g0 + eps0 
            g1 = g0*c1
            d0 = 1
            c0 = c1
            g0 = g1
            do i = 2, n
                d1 = 1-e(i)*x*d0 + eps0
                c1 = 1-e(i)*x/c0 + eps0
                d1 = 1/d1
                delta = c1*d1
                g1 = g0*delta
                if (abs(delta-1) < eps) return
                d0 = d1
                c0 = c1
                g0 = g1
            end do
    
        end function evalcf


        
        

        subroutine pade_coef(t_coef, nt, c_pade_coef, d_pade_coef, n )  
        ! sum_{i=0}^nt t_coef_i x^i ~ (sum_{k=0}^{nt/2} c_coef_k x^k)/(1+sum_{k=1}^{nt/2} d_coef_k x^k)
        
            implicit none
            integer :: nt
            real(dp), dimension(0:ld-1) :: t_coef
            real(dp), dimension(0:ld-1) :: c_pade_coef, d_pade_coef
            real(dp), dimension(ld) :: x, y, r, c, ferr, berr
            real(dp), dimension(4*ld) :: work
            !real(dp), dimension(ld,ld) :: q, qf
            real(dp), allocatable, dimension(:,:) :: q, qf
            integer, dimension(ld) :: ipiv, iwork
            character :: equed
            real(dp) :: rcond
            integer :: info, k, n

            allocate ( q(ld,ld), qf(ld,ld))
            n = nt/2 
            x(1:n) = t_coef(n+1:2*n)
            y(1:n) = x(1:n)
            do k = 1, n
                q(1:n,k) = t_coef(n+1-k:2*n-k)
            end do
        
            call dgesvx('e', 'n', n, 1, q, ld, qf, ld, ipiv, equed, r, c, y, ld, x, ld, &
                        rcond, ferr, berr, work, iwork, info)
                        
            do k = 1, n
                y(k) = t_coef(k) - dot_product(x(1:k),t_coef(k-1:0:-1))
            end do
            
            c_pade_coef(0)     =  t_coef(0) 
            c_pade_coef(1:n)   = y(1:n)
            d_pade_coef(0)     = 1._dp
            d_pade_coef(1:n)   = - x(1:n)

            deallocate (q,qf)
    
        end subroutine pade_coef
    
        function pade_f(x, c_pade_coef, d_pade_coef, n) result(ratval)
        
            implicit none
            real(dp) :: x, ratval
            integer :: n
            real(dp), dimension(0:ld-1) :: c_pade_coef, d_pade_coef
            real(dp) :: sumn, sumd

            sumn = horner( c_pade_coef, n, x )
            sumd = horner( d_pade_coef, n, x )

            ratval = sumn/sumd
        
        end function pade_f

              
        function horner(f,n,x) result(y)

            implicit none
            real(dp), dimension(0:max_taylor) :: f
            integer :: n, i
            real(dp) :: x, y
            
            y = f(n)
            do i = n-1, 0, -1
                y = f(i) + x*y
            end do

        end function horner
              
end module taylor_cf_approx
