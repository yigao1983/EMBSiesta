 subroutine collectdata_dp(nuotot, nuoloc, nspin, Aaux, Aaux_red)
 use precision,    only : dp
 use parallel,     only : Node, Nodes
 use parallelsubs, only : LocalToGlobalOrb
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 integer,  intent(in)  :: nuotot, nuoloc, nspin
 real(dp), intent(in)  ::     Aaux(nuotot,nuoloc,nspin)
 real(dp), intent(out) :: Aaux_red(nuotot,nuotot,nspin)
! Local vaiables
 logical, save :: FirstCall = .true.
 integer  :: ispin, io, iio
!real(dp) :: Aaux_ful(nuotot,nuotot,nspin)
 real(dp), dimension(:,:,:), pointer, save :: Aaux_ful=>null()
#ifdef MPI
 integer :: MPIerror
#endif
!
 if(FirstCall) then
   if(.not.associated(Aaux_ful)) allocate(Aaux_ful(nuotot,nuotot,nspin))
   FirstCall = .false.
 endif
!
 if(nuotot.ne.size(Aaux_ful,1) .or. nuotot.ne.size(Aaux_ful,2) .or. nspin.ne.size(Aaux_ful,3)) then
   if(associated(Aaux_ful)) then
     deallocate(Aaux_ful) ; nullify(Aaux_ful)
   endif
   allocate(Aaux_ful(nuotot,nuotot,nspin))
 endif
!
 Aaux_ful = 0.0_dp
 do ispin = 1,nspin
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       Aaux_ful(:,iio,ispin) = Aaux(:,io,ispin)
    enddo
 enddo
#ifdef MPI
 call MPI_AllReduce(Aaux_ful(1,1,1),Aaux_red(1,1,1),nuotot*nuotot*nspin,MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
 Aaux_red = Aaux_ful
#endif
!
 end subroutine collectdata_dp
