! Yi: sparse Hamiltonian to regular Hamiltonian matrix
 subroutine hsparse2haux(Haux)
 use precision,       only : dp
 use m_spin,          only : nspin
 use sparse_matrices, only : maxnh, numh, listh, listhptr, H
 use m_variables,     only : nuotot, nuoloc
 implicit none
 real(dp), intent(out) :: Haux(nuotot,nuoloc,nspin)
! Local variables
 integer :: ispin, io, jo, j, ind
!
 do ispin = 1,nspin
    Haux(:,:,ispin) = 0.0_dp
    do io = 1,nuoloc
       do j = 1,numh(io)
          ind = listhptr(io) + j
           jo = listh(ind)
          Haux(jo,io,ispin) = Haux(jo,io,ispin) + H(ind,ispin)
       enddo
    enddo
 enddo
!
 end subroutine hsparse2haux
