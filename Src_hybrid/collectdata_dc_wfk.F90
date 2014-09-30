 subroutine collectdata_dc_wfk(nuotot, nuoocc, nuolococc, nspin, zwfk, zwfk_red)
 use precision,    only : dp
 use parallel,     only : Node, Nodes
 use parallelsubs, only : LocalToGlobalOrb
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 integer,     intent(in)  :: nuotot, nuoocc, nuolococc, nspin
 complex(dp), intent(in)  ::     zwfk(nuotot,nuolococc,nspin)
 complex(dp), intent(out) :: zwfk_red(nuotot,   nuoocc,nspin)
! Local vaiables
 logical, save :: FirstCall = .true.
 integer  :: ispin, io, iio
 complex(dp), dimension(:,:,:), pointer, save :: zwfk_ful=>null()
#ifdef MPI
 integer :: MPIerror
#endif
!
 if(FirstCall) then
   if(.not.associated(zwfk_ful)) allocate(zwfk_ful(nuotot,nuoocc,nspin))
   FirstCall = .false.
 endif
!
 if(nuotot.ne.size(zwfk_ful,1) .or. nuoocc.ne.size(zwfk_ful,2) .or. nspin.ne.size(zwfk_ful,3)) then
   if(associated(zwfk_ful)) then
     deallocate(zwfk_ful) ; nullify(zwfk_ful)
   endif
   allocate(zwfk_ful(nuotot,nuoocc,nspin))
 endif
!
 zwfk_ful =(0.0_dp,0.0_dp)
 do ispin = 1,nspin
    do io = 1,nuolococc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       zwfk_ful(:,iio,ispin) = zwfk(:,io,ispin)
    enddo
 enddo
#ifdef MPI
 call MPI_AllReduce(zwfk_ful(1,1,1),zwfk_red(1,1,1),nuotot*nuoocc*nspin,MPI_double_complex,MPI_sum,MPI_Comm_World,MPIerror)
#else
 zwfk_red(:,:,:) = zwfk_ful(:,:,:)
#endif
!
 end subroutine collectdata_dc_wfk
