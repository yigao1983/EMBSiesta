 subroutine mullikpop(nuotot, nuoloc, nspin, daux, natot, population)
 use precision,    only : dp
 use parallel,     only : Node, Nodes, BlockSize
 use parallelsubs, only : LocalToGlobalOrb
 use atomlist,     only : lasto
 use m_variables,  only : Saux
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 integer,  intent(in)  :: nuotot, nuoloc, nspin, natot
 real(dp), intent(in)  :: daux(nuotot,nuoloc,nspin)
 real(dp), intent(out) :: population(natot,nspin)
!
 integer :: ispin, io, iio, iomin, iomax, ia, ierror, info
 logical :: Serial
 logical, save :: FirstCall = .true.
!real(dp), dimension(nuotot,nuoloc,nspin) :: Daux_bar
!real(dp), dimension(nuotot,nuotot,nspin) :: Daux_bar_tmp, Daux_bar_red
 real(dp), dimension(:,:,:), pointer, save :: Daux_bar=>null()
 real(dp), dimension(:,:,:), pointer, save :: Daux_bar_tmp=>null(), Daux_bar_red=>null()
#ifdef MPI
 integer :: nprow, npcol, MPIerror
 integer :: desch(9)
 integer, save :: ictxt
#endif
!
 call timer( 'mullikpop', 1 )
!
 Serial = (Nodes.eq.1)
!
 if(FirstCall) then
   if(.not.associated(Daux_bar))     allocate(    Daux_bar(nuotot,nuoloc,nspin))
   if(.not.associated(Daux_bar_tmp)) allocate(Daux_bar_tmp(nuotot,nuotot,nspin))
   if(.not.associated(Daux_bar_red)) allocate(Daux_bar_red(nuotot,nuotot,nspin))
 endif
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
    Daux_bar(:,:,ispin) = 0.0_dp
    if(Serial) then
      call dgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux,nuotot,Daux(:,:,ispin),nuotot,0.0_dp, &
                 Daux_bar(:,:,ispin),nuotot)
#ifdef MPI
    else
      call pdgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux,1,1,desch,Daux(:,:,ispin),1,1,desch,0.0_dp, &
                  Daux_bar(:,:,ispin),1,1,desch)
#endif
    endif
    Daux_bar_tmp(:,:,ispin) = 0.0_dp
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       Daux_bar_tmp(:,iio,ispin) = Daux_bar(:,io,ispin)
    enddo
#ifdef MPI
    call MPI_AllReduce(Daux_bar_tmp(1,1,ispin), Daux_bar_red(1,1,ispin), nuotot*nuotot, &
         MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerror)
#else
    Daux_bar_red(:,:,ispin) = Daux_bar_tmp(:,:,ispin)
#endif
    do ia = 1,natot
       iomin = lasto(ia-1) + 1
       iomax = lasto(ia)
       population(ia,ispin) = 0.0_dp
       do io = iomin,iomax
          population(ia,ispin) = population(ia,ispin) + Daux_bar_red(io,io,ispin)
       enddo
    enddo
 enddo
!
 FirstCall = .false.
!
 call timer( 'mullikpop', 2 )
!
 end subroutine mullikpop
