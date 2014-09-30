 subroutine siesta_spectra
 use precision,   only : dp
 use parallel,    only : Node
 use files,       only : slabel
 use units,       only : eV, pi
 use m_recipes,   only : four1
 use m_variables, only : estrength, sigma, dtim, maxntim
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
! Local variables
 real(dp), parameter :: absestrtol = 1.0e-9_dp
 integer  :: iu
 integer  :: isign, idim
 integer  :: itim, ipowr
 integer  :: maxpowr, maxnfft
 character(len=200) :: fname
 real(dp) :: dfreq, denrg
 real(dp),    pointer :: time(:)
 real(dp),    pointer :: envl(:)
 real(dp),    pointer :: dipoledftt(:,:), ddipoldftt(:,:)
 real(dp),    pointer :: dipolenfrt(:,:), ddipolnfrt(:,:)
 real(dp),    pointer :: frfg(:,:)
 complex(dp), pointer :: ddipoldftf(:,:), ddipolnfrf(:,:)
#ifdef MPI
 integer :: MPIerror
#endif
!
 external :: io_assign, io_close
!
 maxpowr = 100
 maxnfft = 2
 do ipowr = 1,maxpowr
    if(maxnfft.lt.4*maxntim) maxnfft = 2*maxnfft
 enddo
!
 dfreq = 1.d0/(dtim*dble(maxnfft))
 denrg = 2.d0*pi/eV*dfreq
!
 allocate(        time(maxnfft)) ;       time = 0.0_dp
 allocate(        envl(maxnfft)) ;       envl = 0.0_dp
 allocate(dipoledftt(3,maxnfft)) ; dipoledftt = 0.0_dp
 allocate(ddipoldftt(3,maxnfft)) ; ddipoldftt = 0.0_dp
 allocate(dipolenfrt(3,maxnfft)) ; dipolenfrt = 0.0_dp
 allocate(ddipolnfrt(3,maxnfft)) ; ddipolnfrt = 0.0_dp
 allocate(ddipoldftf(3,maxnfft)) ; ddipoldftf =(0.0_dp,0.0_dp)
 allocate(ddipolnfrf(3,maxnfft)) ; ddipolnfrf =(0.0_dp,0.0_dp)
!
 allocate(frfg(2,maxnfft)) ; frfg = 0.0_dp
!
 fname = slabel
 fname = trim(fname) // '.DIP_ascii'
 if(Node.eq.0) then
   call io_assign(iu)
   open(unit=iu,file=fname,form='formatted',status='old')
   do itim = 1,maxntim
      read(iu,*) time(itim), dipoledftt(:,itim), dipolenfrt(:,itim), ddipoldftt(:,itim), ddipolnfrt(:,itim)
   enddo
   call io_close(iu)
 endif
#ifdef MPI
 call MPI_Bcast(time(1),           maxntim, MPI_double_precision, 0, MPI_Comm_World, MPIerror)
 call MPI_Bcast(dipoledftt(1,1), 3*maxntim, MPI_double_precision, 0, MPI_Comm_World, MPIerror)
 call MPI_Bcast(dipolenfrt(1,1), 3*maxntim, MPI_double_precision, 0, MPI_Comm_World, MPIerror)
 call MPI_Bcast(ddipoldftt(1,1), 3*maxntim, MPI_double_precision, 0, MPI_Comm_World, MPIerror)
 call MPI_Bcast(ddipolnfrt(1,1), 3*maxntim, MPI_double_precision, 0, MPI_Comm_World, MPIerror)
#endif
! Envelope function
 envl = dexp(-sigma*time)
! Initialization
 isign = 1
 do idim = 1,3
!   DFT part
    frfg(1,:) = ddipoldftt(idim,:)*envl(:)
    frfg(2,:) = 0.0_dp
    call four1(frfg, maxnfft, isign)
    ddipoldftf(idim,:) = dcmplx(dtim)*dcmplx(frfg(1,:),frfg(2,:))
!   Near-field part
    frfg(1,:) = ddipolnfrt(idim,:)*envl(:)
    frfg(2,:) = 0.0_dp
    call four1(frfg, maxnfft, isign)
    ddipolnfrf(idim,:) = dcmplx(dtim)*dcmplx(frfg(1,:),frfg(2,:))
 enddo
!
 if(dabs(estrength).gt.absestrtol) then
   ddipoldftf = ddipoldftf/estrength
   ddipolnfrf = ddipolnfrf/estrength
 endif
!
 fname = slabel
 fname = trim(fname) // '.POL_ascii'
 if(Node.eq.0) then
   call io_assign(iu)
   open(unit=iu,file=fname,form='formatted',status='unknown')
   do itim = 1,maxntim
      write(iu,'(13e15.7)') denrg*dble(itim-1), &
      dreal(ddipoldftf(:,itim)),aimag(ddipoldftf(:,itim)), &
      dreal(ddipolnfrf(:,itim)),aimag(ddipolnfrf(:,itim))
   enddo
   call io_close(iu)
 endif
!
 deallocate(time)
 deallocate(envl)
 deallocate(dipoledftt)
 deallocate(ddipoldftt)
 deallocate(dipolenfrt)
 deallocate(ddipolnfrt)
 deallocate(ddipoldftf)
 deallocate(ddipolnfrf)
!
 end subroutine siesta_spectra
