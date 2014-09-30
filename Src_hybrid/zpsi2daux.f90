 subroutine zpsi2daux(zpsi, Daux)
 use precision,    only : dp
 use m_spin,       only : nspin
 use parallel,     only : Node, Nodes
 use parallelsubs, only : LocalToGlobalOrb
 use m_variables,  only : nuotot, nuoloc, nuoocc, nuolococc
 use m_variables,  only : occ
 implicit none
 complex(dp), intent(in)  :: zpsi(nuotot,nuolococc,nspin)
 real(dp),    intent(out) :: Daux(nuotot,nuoloc,   nspin)
! Local variables
 logical, save :: FirstCall = .true.
 integer  :: io, jo, iio, ispin
 complex(dp), dimension(:,:,:), pointer, save :: zpsi_red=>null()
!
 if(FirstCall) then
   if(.not.associated(zpsi_red)) allocate(zpsi_red(nuotot,nuoocc,nspin))
   FirstCall = .false.
 endif
!
 call collectdata_dc_wfk(nuotot, nuoocc, nuolococc, nspin, zpsi, zpsi_red)
! Calculate new density matrix
 do ispin = 1,nspin
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       do jo = 1,nuotot
          Daux(jo,io,ispin) = sum(occ(1:nuoocc,ispin)*dreal(zpsi_red(jo,1:nuoocc,ispin)*dconjg(zpsi_red(iio,1:nuoocc,ispin))))
       enddo
    enddo
 enddo
!
 end subroutine zpsi2daux
