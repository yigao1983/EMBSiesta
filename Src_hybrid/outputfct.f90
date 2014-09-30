 subroutine outputfct(itim, fact)
 use precision,   only : dp
 use files,       only : slabel
 use parallel,    only : Node
 use m_variables, only : dtim
 implicit none
 integer,  intent(in) :: itim
 real(dp), intent(in) :: fact
! Local variables
 integer :: iu
 character*10 :: access_type
 character*99 :: fname
!
 external :: io_assign
!
 fname = slabel
 fname = trim(fname) // '.FCT_ascii'
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
   write(iu,'(2e15.7)') dtim*dble(itim), fact
   call io_close(iu)
 endif
!
 end subroutine outputfct
