 subroutine tdefield_update(itim)
 use units,       only : pi
 use precision,   only : dp
 use m_efield,    only : acting_efield, user_specified_field
 use m_variables, only : fieldtype, tadiabatic, dtim, eext
 implicit none
 integer, intent(in) :: itim
 ! Local variables
 real(dp) :: time, fact
!
 if(fieldtype.eq.1) then
   user_specified_field = 0.0_dp
   acting_efield = .false.
 else
   time = dtim * dble(itim)
   fact = 0.0_dp
   if (dabs(time) .le. dabs(tadiabatic)) &
    fact = 0.5_dp*(1.0_dp + dcos(2.0_dp*pi*time/(2.0_dp*tadiabatic)))
   user_specified_field = eext * fact
   call outputfct(itim, fact)
 endif
!
 end subroutine tdefield_update
