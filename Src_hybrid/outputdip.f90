 subroutine outputdip(itim)
 use files,       only : slabel
 use parallel,    only : Node
 use m_variables, only : dtim
 use m_variables, only : dipdft, dipnfr, ddpdft, ddpnfr
 implicit none
 integer, intent(in) :: itim
! Local variables
 integer :: iu
 character*10 :: access_type
 character*99 :: fname
!
 external :: io_assign
!
 fname = slabel
 fname = trim(fname) // '.DIP_ascii'
!
 if(itim.eq.0) then
   access_type = 'sequential'
 else
   access_type = 'append'
 endif
!
 if(Node.eq.0) then
   call io_assign(iu)
   open(unit=iu,file=fname,form='formatted',access=access_type)
   write(iu,'(13e15.7)') dtim*dble(itim), dipdft, dipnfr, ddpdft, ddpnfr
   call io_close(iu)
 endif
!
 end subroutine outputdip
