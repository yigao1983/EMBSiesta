 subroutine propagat_zdaux(dt, Haux, zDaux_old, zDaux_new)
 use precision,    only : dp
 use parallel,     only : Node, Nodes, BlockSize
 use parallelsubs, only : LocalToGlobalOrb
 use m_spin,       only : nspin
 use m_variables,  only : simplempi
 use m_variables,  only : i_dc
 use m_variables,  only : nuotot, nuoloc
 use m_variables,  only : Saux_nhf, Saux_phf
 use m_variables,  only : Saux_nhf_red
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 real(dp),    intent(in)  :: dt
 real(dp),    intent(in)  ::      Haux(nuotot,nuoloc,nspin)
 complex(dp), intent(in)  :: zDaux_old(nuotot,nuoloc,nspin)
 complex(dp), intent(out) :: zDaux_new(nuotot,nuoloc,nspin)
! Local variables
 integer :: ispin, io, iio, jo, ierror, info
 logical :: Serial
 logical, save :: FirstCall = .true.
!real(dp),    dimension(nuotot,nspin)        :: w
!real(dp),    dimension(nuotot,nuoloc)       :: Faux, Iaux
!real(dp),    dimension(nuotot,nuoloc,nspin) :: Haux_bar
!real(dp),    dimension(nuotot,nuoloc,nspin) :: z, zt
 real(dp),    dimension(:,:),   save, pointer:: w=>null()
 real(dp),    dimension(:,:),   save, pointer:: Faux=>null(), Iaux=>null()
 real(dp),    dimension(:,:,:), save, pointer:: Haux_bar=>null()
 real(dp),    dimension(:,:,:), save, pointer:: z=>null(), zt=>null()
!real(dp),    dimension(nuotot,nuotot,nspin) :: z_red
!complex(dp), dimension(nuotot,nuoloc,nspin) :: zw
!complex(dp), dimension(nuotot,nuoloc,nspin) :: expOaux, expOSphfaux
!complex(dp), dimension(nuotot,nuoloc,nspin) :: Uaux, UzDaux, Uaux_dag
 complex(dp), dimension(:,:,:), save, pointer:: zw=>null()
 complex(dp), dimension(:,:,:), save, pointer:: expOaux=>null(), expOSphfaux=>null()
 complex(dp), dimension(:,:,:), save, pointer:: Uaux=>null(), UzDaux=>null(), Uaux_dag=>null()
!complex(dp), dimension(nuotot,nuotot,nspin) :: zw_red
!complex(dp), dimension(nuotot,nuotot,nspin) :: expOaux_red
!complex(dp), dimension(nuotot,nuotot,nspin) :: Uaux_red, UzDaux_red
 real(dp),    dimension(:,:,:), pointer, save:: z_red=>null()
 complex(dp), dimension(:,:,:), pointer, save:: zw_red=>null()
 complex(dp), dimension(:,:,:), pointer, save:: expOaux_red=>null()
 complex(dp), dimension(:,:,:), pointer, save:: Uaux_red=>null(), UzDaux_red=>null()
#ifdef MPI
 integer :: MPIerror
 integer :: nprow, npcol
 integer :: desch(9)
 integer, save :: ictxt
#endif
!
 call timer( 'propzdaux', 1 )
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
   allocate(         zt(nuotot,nuoloc,nspin))
   allocate(         zw(nuotot,nuoloc,nspin))
   allocate(    expOaux(nuotot,nuoloc,nspin))
   allocate(expOSphfaux(nuotot,nuoloc,nspin))
   allocate(       Uaux(nuotot,nuoloc,nspin))
   allocate(     UzDaux(nuotot,nuoloc,nspin))
   allocate(   Uaux_dag(nuotot,nuoloc,nspin))
   if(simplempi) then
     allocate(      z_red(nuotot,nuotot,nspin))
     allocate(     zw_red(nuotot,nuotot,nspin))
     allocate(expOaux_red(nuotot,nuotot,nspin))
     allocate(   Uaux_red(nuotot,nuotot,nspin))
     allocate( UzDaux_red(nuotot,nuotot,nspin))
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
!   do io = 1,nuoloc
!      do jo = 1,nuotot
!         Faux(jo,io) = Haux_bar(jo,io,ispin)
!      enddo
!   enddo
    Faux(:,:) = Haux_bar(:,:,ispin)
    call rdiag(Faux,Iaux,nuotot,nuoloc,nuotot,w(1,ispin),z(1,1,ispin),nuotot,1,ierror)
 enddo
!
#ifdef MPI
 if(simplempi) then
   call collectdata_dp(nuotot, nuoloc, nspin, z, z_red)
   do ispin = 1,nspin
      do io = 1,nuoloc
         call LocalToGlobalOrb(io,Node,Nodes,iio)
         zt(:,io,ispin) = z_red(iio,:,ispin)
      enddo
   enddo
 endif
#endif
!
 do ispin = 1,nspin
    do io = 1,nuoloc
       call LocalToGlobalOrb(io,Node,Nodes,iio)
       zw(:,io,ispin) = dcmplx(z(:,io,ispin)) * exp(-i_dc*dt*w(iio,ispin))
    end do
!
!   do io = 1,nuoloc
!      do jo = 1,nuotot
!             expOaux(jo,io,ispin) =(0.0_dp,0.0_dp)
!         expOSphfaux(jo,io,ispin) =(0.0_dp,0.0_dp)
!                Uaux(jo,io,ispin) =(0.0_dp,0.0_dp)
!              UzDaux(jo,io,ispin) =(0.0_dp,0.0_dp)
!           zDaux_new(jo,io,ispin) =(0.0_dp,0.0_dp)
!      enddo
!   enddo
        expOaux(:,:,ispin) =(0.0_dp,0.0_dp)
    expOSphfaux(:,:,ispin) =(0.0_dp,0.0_dp)
           Uaux(:,:,ispin) =(0.0_dp,0.0_dp)
         UzDaux(:,:,ispin) =(0.0_dp,0.0_dp)
      zDaux_new(:,:,ispin) =(0.0_dp,0.0_dp)
!
    if(Serial) then
      call zgemm('N','C',nuotot,nuotot,nuotot,1.0_dp,zw(:,:,ispin),nuotot,dcmplx(z(:,:,ispin)),nuotot,0.0_dp,      &
                 expOaux(:,:,ispin),nuotot)
      call zgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,expOaux(:,:,ispin),nuotot,dcmplx(Saux_phf),nuotot,0.0_dp,     &
                 expOSphfaux(:,:,ispin),nuotot)
      call zgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,dcmplx(Saux_nhf),nuotot,expOSphfaux(:,:,ispin),nuotot,0.0_dp, &
                 Uaux(:,:,ispin),nuotot)
      call zgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Uaux(:,:,ispin),nuotot,zDaux_old(:,:,ispin),nuotot,0.0_dp,    &
                 UzDaux(:,:,ispin),nuotot)
      call zgemm('N','C',nuotot,nuotot,nuotot,1.0_dp,UzDaux(:,:,ispin),nuotot,Uaux(:,:,ispin),nuotot,0.0_dp,       &
                 zDaux_new(:,:,ispin),nuotot)
#ifdef MPI
    else
      if(simplempi) then
        call collectdata_dc(nuotot, nuoloc, 1, zw(:,:,ispin), zw_red(:,:,ispin))
        expOaux(:,:,ispin)     = matmul(zw_red(:,:,ispin), zt(:,:,ispin))
!       call zgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,zw_red(:,:,ispin),nuotot,dcmplx(zt(:,:,ispin)),nuotot,0.0_dp,     &
!                  expOaux(:,:,ispin),nuotot)
        call collectdata_dc(nuotot, nuoloc, 1, expOaux(:,:,ispin), expOaux_red(:,:,ispin))
        expOSphfaux(:,:,ispin) = matmul(expOaux_red(:,:,ispin), dcmplx(Saux_phf))
!       call zgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,expOaux_red(:,:,ispin),nuotot,dcmplx(Saux_phf),nuotot,0.0_dp,     &
!                  expOSphfaux(:,:,ispin),nuotot)
        Uaux(:,:,ispin) = matmul(dcmplx(Saux_nhf_red), expOSphfaux(:,:,ispin))
!       call zgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,dcmplx(Saux_nhf_red),nuotot,expOSphfaux(:,:,ispin),nuotot,0.0_dp, &
!                  Uaux(:,:,ispin),nuotot)
        call collectdata_dc(nuotot, nuoloc, 1, Uaux(:,:,ispin), Uaux_red(:,:,ispin))
        do io = 1,nuoloc
           call LocalToGlobalOrb(io,Node,Nodes,iio)
           Uaux_dag(:,io,ispin) = dconjg(Uaux_red(iio,:,ispin))
        enddo
        UzDaux(:,:,ispin) = matmul(Uaux_red(:,:,ispin), zDaux_old(:,:,ispin))
!       call zgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,Uaux_red(:,:,ispin),nuotot,zDaux_old(:,:,ispin),nuotot,0.0_dp,    &
!                  UzDaux(:,:,ispin),nuotot)
        call collectdata_dc(nuotot, nuoloc, 1, UzDaux(:,:,ispin), UzDaux_red(:,:,ispin))
        zDaux_new(:,:,ispin) = matmul(UzDaux_red(:,:,ispin), Uaux_dag(:,:,ispin))
!       call zgemm('N','N',nuotot,nuoloc,nuotot,1.0_dp,UzDaux_red(:,:,ispin),nuotot,Uaux_dag(:,:,ispin),nuotot,0.0_dp,   &
!                  zDaux_new(:,:,ispin),nuotot)
        call MPI_Barrier(MPI_Comm_World,MPIerror)
      else
        call pzgemm('N','C',nuotot,nuotot,nuotot,1.0_dp,zw(:,:,ispin),1,1,desch,dcmplx(z(:,:,ispin)),1,1,desch,0.0_dp,      &
                    expOaux(:,:,ispin),1,1,desch)
        call pzgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,expOaux(:,:,ispin),1,1,desch,dcmplx(Saux_phf),1,1,desch,0.0_dp,     &
                    expOSphfaux(:,:,ispin),1,1,desch)
        call pzgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,dcmplx(Saux_nhf),1,1,desch,expOSphfaux(:,:,ispin),1,1,desch,0.0_dp, &
                    Uaux(:,:,ispin),1,1,desch)
        call pzgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Uaux(:,:,ispin),1,1,desch,zDaux_old(:,:,ispin),1,1,desch,0.0_dp,    &
                    UzDaux(:,:,ispin),1,1,desch)
        call pzgemm('N','C',nuotot,nuotot,nuotot,1.0_dp,UzDaux(:,:,ispin),1,1,desch,Uaux(:,:,ispin),1,1,desch,0.0_dp,       &
                    zDaux_new(:,:,ispin),1,1,desch)
      endif
#endif
    endif
 enddo
!
 FirstCall = .false.
!
 call timer( 'propzdaux', 2 )
!
 end subroutine propagat_zdaux
