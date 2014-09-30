 subroutine propagat_zwfk(dt, Haux, zwfk_old, zwfk_new)
 use precision,    only : dp
 use parallel,     only : Node, Nodes, BlockSize
 use parallelsubs, only : LocalToGlobalOrb
 use m_spin,       only : nspin
 use m_variables,  only : simplempi
 use m_variables,  only : i_dc
 use m_variables,  only : nuotot, nuoloc, nuoocc, nuolococc
 use m_variables,  only : Saux_nhf,     Saux_phf
 use m_variables,  only : Saux_nhf_red, Saux_phf_red
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 real(dp),    intent(in)  :: dt
 real(dp),    intent(in)  ::     Haux(nuotot,nuoloc,   nspin)
 complex(dp), intent(in)  :: zwfk_old(nuotot,nuolococc,nspin)
 complex(dp), intent(out) :: zwfk_new(nuotot,nuolococc,nspin)
! Local variables
 integer :: ispin, io, iio, jo, ierror, info
 logical :: Serial
 logical, save :: FirstCall = .true.
 real(dp),    dimension(:,:),   save, pointer:: w=>null()
 real(dp),    dimension(:,:),   save, pointer:: Faux=>null(), Iaux=>null()
 real(dp),    dimension(:,:,:), save, pointer:: Haux_bar=>null()
 real(dp),    dimension(:,:,:), save, pointer:: z=>null()
 complex(dp), dimension(:,:,:), save, pointer:: zw=>null()
!
 real(dp),    dimension(:,:,:), pointer, save:: z_red=>null()
 complex(dp), dimension(:,:,:), pointer, save:: zw_red=>null()
#ifdef MPI
 integer :: MPIerror
 integer :: nprow, npcol
 integer :: desch(9)
 integer, save :: ictxt
#endif
!
 call timer( 'propzwfk', 1 )
!
 Serial = (Nodes.eq.1)
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
 if(FirstCall) then
   allocate(          w(nuotot,nspin))
   allocate(       Faux(nuotot,nuoloc))
   allocate(       Iaux(nuotot,nuoloc))
   allocate(   Haux_bar(nuotot,nuoloc,nspin))
   allocate(          z(nuotot,nuoloc,nspin))
   allocate(         zw(nuotot,nuoloc,nspin))
   if(simplempi) then
     allocate(      z_red(nuotot,nuotot,nspin))
     allocate(     zw_red(nuotot,nuotot,nspin))
   endif
 endif
!
 call haux2hbar(Haux, Haux_bar)
!
 Iaux = 0.0_dp
 do io = 1,nuoloc
    call LocalToGlobalOrb(io,Node,Nodes,iio)
    Iaux(iio,io) = 1.0_dp
 enddo
!
 do ispin = 1,nspin
    Faux(:,:) = Haux_bar(:,:,ispin)
    call rdiag(Faux,Iaux,nuotot,nuoloc,nuotot,w(1,ispin),z(1,1,ispin),nuotot,1,ierror)
 enddo
!
 do ispin = 1,nspin
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       zw(:,io,ispin) = dcmplx(z(:,io,ispin)) * exp(-i_dc*dt*w(iio,ispin))
    end do
 end do
#ifdef MPI
 if(simplempi) call collectdata_dp(nuotot, nuoloc, nspin, z,  z_red)
 if(simplempi) call collectdata_dc(nuotot, nuoloc, nspin, zw, zw_red)
#endif
!
 do ispin = 1,nspin
!
    if(Serial) then
      call zgemm('N','N',nuotot,nuoocc,nuotot,1.0_dp,dcmplx(Saux_phf),nuotot,zwfk_old(:,:,ispin),nuotot,0.0_dp,     &
                 zwfk_new(:,:,ispin),nuotot)
      call zgemm('C','N',nuotot,nuoocc,nuotot,1.0_dp,dcmplx(z(:,:,ispin)),nuotot,zwfk_new(:,:,ispin),nuotot,0.0_dp, &
                 zwfk_new(:,:,ispin),nuotot)
      call zgemm('N','N',nuotot,nuoocc,nuotot,1.0_dp,zw(:,:,ispin),nuotot,zwfk_new(:,:,ispin),nuotot,0.0_dp,        &
                 zwfk_new(:,:,ispin),nuotot)
      call zgemm('N','N',nuotot,nuoocc,nuotot,1.0_dp,dcmplx(Saux_nhf),nuotot,zwfk_new(:,:,ispin),nuotot,0.0_dp,     &
                 zwfk_new(:,:,ispin),nuotot)
#ifdef MPI
    else
      if(simplempi) then
        if(nuolococc.gt.0) zwfk_new(:,1:nuolococc,ispin) = matmul(dcmplx(Saux_phf_red),zwfk_old(:,1:nuolococc,ispin))
        if(nuolococc.gt.0) zwfk_new(:,1:nuolococc,ispin) = matmul(dcmplx(transpose(z_red(:,:,ispin))),zwfk_new(:,1:nuolococc,ispin))
        if(nuolococc.gt.0) zwfk_new(:,1:nuolococc,ispin) = matmul(zw_red(:,:,ispin),zwfk_new(:,1:nuolococc,ispin))
        if(nuolococc.gt.0) zwfk_new(:,1:nuolococc,ispin) = matmul(dcmplx(Saux_nhf_red),zwfk_new(:,1:nuolococc,ispin))
        call MPI_Barrier(MPI_Comm_World,MPIerror)
      else
        call pzgemm('N','N',nuotot,nuoocc,nuotot,1.0_dp,dcmplx(Saux_phf),1,1,desch,dcmplx(zwfk_old(:,:,ispin)),1,1,desch,0.0_dp,&
                    zwfk_new(:,:,ispin),1,1,desch)
        call pzgemm('C','N',nuotot,nuoocc,nuotot,1.0_dp,dcmplx(z(:,:,ispin)),1,1,desch,zwfk_new(:,:,ispin),1,1,desch,0.0_dp,    &
                    zwfk_new(:,:,ispin),1,1,desch)
        call pzgemm('N','N',nuotot,nuoocc,nuotot,1.0_dp,zw(:,:,ispin),1,1,desch,zwfk_new(:,:,ispin),1,1,desch,0.0_dp,           &
                    zwfk_new(:,:,ispin),1,1,desch)
        call pzgemm('N','N',nuotot,nuoocc,nuotot,1.0_dp,dcmplx(Saux_nhf),1,1,desch,zwfk_new(:,:,ispin),1,1,desch,0.0_dp,        &
                    zwfk_new(:,:,ispin),1,1,desch)
      endif
#endif
    endif
 enddo
!
 FirstCall = .false.
!
 call timer( 'propzwfk', 2 )
!
 end subroutine propagat_zwfk
