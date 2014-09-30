 subroutine outputpop(itim, natot, nspin, pop)
 use precision,   only : dp
 use files,       only : slabel
 use parallel,    only : Node
 use m_variables, only : dtim
 implicit none
 integer,  intent(in) :: itim
 integer,  intent(in) :: natot, nspin
 real(dp), intent(in) :: pop(natot,nspin)
! Local variables
 integer :: iu, ispin, ia
 character*10 :: access_type
 character*99 :: fname
!
 external :: io_assign
!
 fname = slabel
 fname = trim(fname) // '.POP_ascii'
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
   write(iu,'(e15.7)',advance='no') dtim*dble(itim)
   do ispin = 1,nspin
      do ia = 1,natot
         write(iu,'(e15.7)',advance='no') pop(ia,ispin)
      enddo
   enddo
   write(iu,*)
   call io_close(iu)
 endif
!
 end subroutine outputpop
