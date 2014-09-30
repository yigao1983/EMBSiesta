 subroutine psi2daux(psi, Daux)
 use precision,    only : dp
 use m_spin,       only : nspin
 use parallel,     only : Node, Nodes
 use parallelsubs, only : LocalToGlobalOrb
 use m_variables,  only : absocctol
 use m_variables,  only : nuotot, nuoloc, nuoocc
 use m_variables,  only : occ
 implicit none
 real(dp), intent(in)  ::  psi(nuotot,nuoloc,nspin)
 real(dp), intent(out) :: Daux(nuotot,nuoloc,nspin)
! Local variables
 logical, save :: FirstCall = .true.
 integer  :: io, jo, iio, ispin
 real(dp), dimension(:,:,:), pointer, save :: psi_red=>null()
!
 if(FirstCall) then
   if(.not.associated(psi_red)) allocate(psi_red(nuotot,nuotot,nspin))
 endif
!
 call collectdata_dp(nuotot, nuoloc, nspin, psi, psi_red)
! Calculate new density matrix
 do ispin = 1,nspin
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       do jo = 1,nuotot
          Daux(jo,io,ispin) = sum(occ(1:nuoocc,ispin)*psi_red(jo,1:nuoocc,ispin)*psi_red(iio,1:nuoocc,ispin))
       enddo
    enddo
 enddo
!
 FirstCall = .false.
!
 end subroutine psi2daux
