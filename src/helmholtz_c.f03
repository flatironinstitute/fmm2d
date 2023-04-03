module c_helmholtz
  use iso_c_binding
contains
  subroutine c_hfmm2d_s_c_p(eps, zk, ns, sources, charge, pot, ier) &
       bind(c, name='hfmm2d_s_c_p')
    implicit none
    real(c_double), value :: eps
    complex(c_double_complex), value, intent(in) :: zk
    integer(c_int64_t), value :: ns
    real(c_double), intent(in) :: sources(2, ns)
    complex(c_double_complex), intent(in) :: charge(ns)
    complex(c_double_complex), intent(out) :: pot(ns)
    integer(c_int64_t), intent(out) :: ier
    call hfmm2d_s_c_p(eps, zk, ns, sources, charge, pot, ier)
  end subroutine c_hfmm2d_s_c_p
end module c_helmholtz
