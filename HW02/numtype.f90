module Numtype

    save
    integer, parameter :: dp = selected_real_kind(15,307)
    !integer , paramter :: qp = selected_real_kind(33,4931)
    real(dp), parameter :: pi = 4*atan(1._dp)
    real(dp), parameter :: iic = (0._dp, 1._dp)
    integer, parameter :: n_eq = 4

end module Numtype