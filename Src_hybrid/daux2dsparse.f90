 subroutine daux2dsparse(Daux)
 use precision,       only : dp
 use alloc,           only : re_alloc
 use m_spin,          only : nspin
 use sparse_matrices, only : maxnh, listh, listhptr, numh, Dscf
 use m_variables,     only : nuotot, nuoloc
 implicit none
 real(dp), intent(in) :: Daux(nuotot,nuoloc,nspin)
! Local variables
 integer :: ispin, io, jo, j, ind
!
 if(.not.associated(Dscf)) call re_alloc(Dscf,1,maxnh,1,nspin,name='Dscf',routine='daux2dsparse')
!
 Dscf = 0.0_dp
 do ispin = 1,nspin
    do io = 1,nuoloc
       do j = 1,numh(io)
          ind = listhptr(io) + j
           jo = listh(ind)
          Dscf(ind,ispin) = Daux(jo,io,ispin)
       enddo
    enddo
 enddo
!
 end subroutine daux2dsparse
