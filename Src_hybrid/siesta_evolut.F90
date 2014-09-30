 subroutine siesta_evolut
 use precision,     only : dp
 use parallel,      only : Node
 use m_spin,        only : nspin
 use m_dipol,       only : dipol
 use m_variables
 use nearfieldsubs, only : leap_frog, pol2rho3d, pol2dip3d, phi2efd3d
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
! Local variables
 integer  :: itim, idim, ispin, io, jo, ia
 logical  :: newS = .false.
 real(dp) :: TrSDaux(nspin)
 real(dp) :: pop(natot,nspin)
#ifdef MPI
 integer :: MPIerror
#endif
!
 call timer( 'evolution', 1 )
!
 if(Node.eq.0) write(6,'(/a/)') 'siesta_evolut: Evolution type: daux'
!
 stagenfr = 'td'
! 
 dpdft0 = 0.0_dp ; dpnfr0 = 0.0_dp
 dipdft = 0.0_dp ; dipnfr = 0.0_dp
 ddpdft = 0.0_dp ; ddpnfr = 0.0_dp
!
 if(associated(zDaux))     deallocate(zDaux)
 allocate(    zDaux(nuotot,nuoloc,nspin))
 if(associated(zDaux_old)) deallocate(zDaux_old)
 allocate(zDaux_old(nuotot,nuoloc,nspin))
!
 if(nonscfdm.or.nonscfwf) call renormdaux(Daux)
!
 zDaux(:,:,:) = dcmplx(Daux(:,:,:))
!
 itim = 0
 call sauxdaux2haux(itim, newS, Saux, Daux, Haux)
! Initial dipoles
 dpdft0 =-dipol
 dipdft = dpdft0
 ddpdft = dipdft - dpdft0
 if(.not.pseudonfr) then
   call pol2dip3d(poltot, dpnfr0)
   dipnfr = dpnfr0
   ddpnfr = dipnfr - dpnfr0
 endif
!
 call outputdip(0)
!
 if (fieldtype.ne.1) call outputfct(0,1.0_dp)
!
 if(writetdpop) then
!  call lowdinpop(nuotot, nuoloc, nspin, Daux, natot, pop)
   call mullikpop(nuotot, nuoloc, nspin, Daux, natot, pop)
   call outputpop(0, natot, nspin, pop)
 endif
!
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
!
!       acting_efield = .true.
!user_specified_field = eext !/dtim ! user_specified_field is positive
 call tdefield_init()
!
 if(.not.pseudonfr) then
   do idim = 1,3
      efd(:,:,:,idim) = efd(:,:,:,idim) - eext(idim)/dtim
   enddo
 endif
!
 do itim = 1,maxntim-1
!
    zDaux_old(:,:,:) = zDaux(:,:,:)
    ! TDDFT
    if(.not.pseudoatm) call magnus(itim, dtim, zDaux_old, zDaux)
    ! Near-field
    if(.not.pseudonfr) then
      call leap_frog(dtim, efd, curtot, poltot)
      ! Update near-field density and dipole
      call pol2rho3d(poltot, rhonfr)
      call pol2dip3d(poltot, dipnfr)
    endif
    ! Update density matrix
    Daux(:,:,:) = dreal(zDaux(:,:,:))
    ! Normalization
    if(chknormdm) call normdaux(Daux, TrSDaux)
    if(Node.eq.0) then
      write(6,'(/a15,i8)',advance='no') 'evolution: ', itim
      if(.not.pseudoatm .and. chknormdm) write(6,'(e15.7)',advance='no') sum(TrSDaux(:))
      write(6,*)
      call flush(6)
    endif
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
    ! Electric potential and field
    call sauxdaux2haux(itim, newS, Saux, Daux, Haux)
    if(.not.pseudonfr) call phi2efd3d(phi, efd)
    ! Update dipoles
    if(.not.pseudoatm) then
      dipdft =-dipol
      ddpdft = dipdft - dpdft0
    endif
    if(.not.pseudonfr) ddpnfr = dipnfr - dpnfr0
    ! Write to file
    call outputdip(itim)
    if(writetdpop) then
!     call lowdinpop(nuotot, nuoloc, nspin, Daux, natot, pop)
      call mullikpop(nuotot, nuoloc, nspin, Daux, natot, pop)
      call outputpop(itim, natot, nspin, pop)
    endif
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
    ! No external field any more
!   user_specified_field = 0.0_dp
!          acting_efield = .false.
    call tdefield_update(itim)
!
 enddo
!
 call timer( 'evolution', 2 )
!
 end subroutine siesta_evolut
