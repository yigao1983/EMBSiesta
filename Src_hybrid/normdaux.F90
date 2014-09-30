 subroutine normdaux(Daux, TrSDaux)
 use precision,    only : dp
 use parallel,     only : Node, Nodes, BlockSize
 use parallelsubs, only : LocalToGlobalOrb
 use m_spin,       only : nspin
 use m_variables,  only : nuotot, nuoloc
 use m_variables,  only : Saux
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 integer :: info
 logical :: Serial
 logical, save :: FirstCall = .true.
 real(dp), intent(in)  :: Daux(nuotot,nuoloc,nspin)
 real(dp), intent(out) :: TrSDaux(nspin)
! Local variables
 integer :: ispin, io, jo, iio
 real(dp), dimension(Nodes,nspin)         :: TrSDaux_tmp, TrSDaux_red
 real(dp), dimension(nuotot,nuoloc,nspin) :: SDaux
#ifdef MPI
 integer :: nprow, npcol
 integer :: desch(9)
 integer :: MPIerror
 integer, save :: ictxt
#endif
!
 call timer( 'normdaux', 1 )
!
 Serial = (Nodes.eq.1)
!
#ifdef MPI
 if(.not.Serial) then
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
    SDaux(:,:,ispin) = 0.0_dp
    if(Serial) then
      call dgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux,nuotot,Daux(:,:,ispin),nuotot,0.0_dp,SDaux(:,:,ispin),nuotot)
#ifdef MPI
    else
      call pdgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux,1,1,desch,Daux(:,:,ispin),1,1,desch,0.0_dp,SDaux(:,:,ispin),1,1,desch)
#endif
    endif
    TrSDaux_tmp(:,ispin) = 0.0_dp
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       TrSDaux_tmp(Node+1,ispin) = TrSDaux_tmp(Node+1,ispin) + SDaux(iio,io,ispin)
    enddo
#ifdef MPI
    call MPI_AllReduce(TrSDaux_tmp(1,ispin), TrSDaux_red(1,ispin), Nodes, &
         MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerror)
#else
    TrSDaux_red(:,ispin) = TrSDaux_tmp(:,ispin)
#endif
    TrSDaux(ispin) = sum(TrSDaux_red(:,ispin))
 enddo
!
 FirstCall = .false.
!
 call timer( 'normdaux', 2 )
!
 end subroutine normdaux
