 subroutine haux2hbar(Haux, Haux_bar)
 use precision,   only : dp
 use alloc,       only : re_alloc
 use parallel,    only : Node, Nodes, BlockSize
 use m_spin,      only : nspin
 use m_variables, only : simplempi
 use m_variables, only : nuotot, nuoloc
 use m_variables, only : Saux_nhf, Saux_nhf_red
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 integer :: info
 logical :: Serial
 logical, save :: FirstCall = .true.
 real(dp), intent(in)  ::     Haux(nuotot,nuoloc,nspin)
 real(dp), intent(out) :: Haux_bar(nuotot,nuoloc,nspin)
! Local variables
 integer  :: ispin
!real(dp) :: HSnhfaux(nuotot,nuoloc,nspin)
 real(dp), dimension(:,:,:), pointer, save :: HSnhfaux=>null()
 real(dp), dimension(:,:,:), pointer, save :: Haux_red=>null()
#ifdef MPI
 integer :: MPIerror
 integer :: nprow, npcol
 integer :: desch(9)
 integer, save :: ictxt
#endif
!
 call timer( 'haux2hbar', 1 )
!
 Serial = (Nodes.eq.1)
!
 if(FirstCall) then
   if(.not.associated(HSnhfaux)) allocate(HSnhfaux(nuotot,nuoloc,nspin))
   if(simplempi) then
   if(.not.associated(Haux_red)) allocate(Haux_red(nuotot,nuotot,nspin))
   endif
 endif
!
#ifdef MPI
 if(.not.Serial .and. .not.simplempi) then
   nprow = 1
   npcol = Nodes
   if(FirstCall) then
     call blacs_get( -1, 0, ictxt )
     call blacs_gridinit( ictxt, 'C', nprow, npcol )
   endif
!
   call descinit( desch, nuotot, nuotot, BlockSize, BlockSize, 0, 0, ictxt, nuotot, info )
!
 endif
#endif
!
 do ispin = 1,nspin
    HSnhfaux(:,:,ispin) = 0.0_dp
    if(Serial) then
      call dgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Haux(:,:,ispin),nuotot,Saux_nhf,nuotot,0.0_dp,     &
                 HSnhfaux(:,:,ispin),nuotot)
      call dgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux_nhf,nuotot,HSnhfaux(:,:,ispin),nuotot,0.0_dp, &
                 Haux_bar(:,:,ispin),nuotot)
#ifdef MPI
    else
      if(simplempi) then
        call collectdata_dp(nuotot, nuoloc, 1, Haux(:,:,ispin), Haux_red(:,:,ispin))
        HSnhfaux(:,:,ispin) = matmul(Haux_red(:,:,ispin), Saux_nhf)
!       call dgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,Haux_red(:,:,ispin),nuotot,Saux_nhf,nuotot,0.0_dp,     &
!                  HSnhfaux(:,:,ispin),nuotot)
!       call collectdata_dp(nuotot, nuoloc, 1, Saux_nhf, Saux_nhf_red)
        Haux_bar(:,:,ispin) = matmul(Saux_nhf_red, HSnhfaux(:,:,ispin))
!       call dgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,Saux_nhf_red,nuotot,HSnhfaux(:,:,ispin),nuotot,0.0_dp, &
!                  Haux_bar(:,:,ispin),nuotot)
        call MPI_Barrier(MPI_Comm_World,MPIerror)
      else
        call pdgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Haux(:,:,ispin),1,1,desch,Saux_nhf,1,1,desch,0.0_dp,     &
                    HSnhfaux(:,:,ispin),1,1,desch)
        call pdgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux_nhf,1,1,desch,HSnhfaux(:,:,ispin),1,1,desch,0.0_dp, &
                    Haux_bar(:,:,ispin),1,1,desch)
      endif
#endif
    endif
 enddo
!
 FirstCall = .false.
!
 call timer( 'haux2hbar', 2 )
!
 end subroutine haux2hbar
