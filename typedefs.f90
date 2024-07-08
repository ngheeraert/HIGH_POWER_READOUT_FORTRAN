  module typedefs
  implicit none
  integer, parameter :: sp=kind(1.0)
  integer, parameter :: spc=kind((1.0,1.0))
  integer, parameter :: dp=selected_real_kind(15)
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
#ifdef DP
  integer, parameter :: r_type = dp 
  integer, parameter :: c_type = dpc 
#else
  integer, parameter :: r_type = sp 
  integer, parameter :: c_type = spc 
#endif
end module typedefs
