 subroutine state_stabiliz
 use precision,         only : dp
 use alloc,             only : de_alloc, re_alloc
 use parallel,          only : Node, Nodes, BlockSize
 use siesta_options,    only : maxsav, wmix, nkick, wmixkick
 use sparse_matrices,   only : Dold, Dscf
 use m_ntm,             only : ntm
 use m_spin,            only : nspin
 use m_variables,       only : pseudonfr, pseudoatm
 use m_variables,       only : chknormdm, chkstabl
 use m_variables,       only : maxnstb, abstolstb
 use m_variables,       only : natot, nuotot, nuoloc
 use m_variables,       only : eig, occ
 use m_variables,       only : wfk
 use m_variables,       only : Saux, Haux, Daux, Daux_old
 use m_variables,       only : Saux_inv
 use m_hyb_pulay_mixer, only : finalize_daux_pulay
 use m_hyb_pulay_mixer, only : daux_pulay_mixer
 use nearfieldsubs,     only : outnfr
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
! Local variables
 integer  :: istb, ispin, io, jo, iio, ia, ierror, info
 logical  :: Serial, newS, cnvg
 real(dp) :: ratio, dDmax, TrSDaux(nspin)
 real(dp) :: abscommutdiff(nspin)
 logical, save :: FirstCall = .true.
 real(dp), dimension(Nodes)               :: dDmax_tmp, dDmax_red
!real(dp), dimension(natot, nspin)        :: population
!real(dp), dimension(nuotot,nspin)        :: eig_new
!real(dp), dimension(nuotot,nuoloc)       :: Faux, Gaux
!real(dp), dimension(nuotot,nuoloc,nspin) :: wfk_new
!real(dp), dimension(nuotot,nuoloc,nspin) :: HDaux
!real(dp), dimension(nuotot,nuoloc,nspin) :: SinvHDaux
!real(dp), dimension(nuotot,nuotot,nspin) :: SinvHDaux_red
 real(dp), dimension(:,:),   pointer, save :: eig_new=>null()
 real(dp), dimension(:,:),   pointer, save :: Faux=>null(), Gaux=>null()
 real(dp), dimension(:,:,:), pointer, save :: wfk_new=>null()
 real(dp), dimension(:,:,:), pointer, save :: HDaux=>null()
 real(dp), dimension(:,:,:), pointer, save :: SinvHDaux=>null()
 real(dp), dimension(:,:,:), pointer, save :: SinvHDaux_red=>null()
! For debug
!real(dp), dimension(nuotot,nuotot,nspin) :: wfk_new_red
#ifdef MPI
 integer :: desch(9)
 integer :: MPIerror
 integer, save :: ictxt
 integer, save :: nprow, npcol
#endif
!
 external :: collectdata_dp
!
 call timer( 'stabiliz', 1 )
!
 Serial = (Nodes.eq.1)
!
 if(.not.associated(eig_new))       allocate(      eig_new(nuotot,nspin))
 if(.not.associated(Faux))          allocate(         Faux(nuotot,nuoloc))
 if(.not.associated(Gaux))          allocate(         Gaux(nuotot,nuoloc))
 if(.not.associated(wfk_new))       allocate(      wfk_new(nuotot,nuoloc,nspin))
 if(.not.associated(HDaux))         allocate(        HDaux(nuotot,nuoloc,nspin))
 if(.not.associated(SinvHDaux))     allocate(    SinvHDaux(nuotot,nuoloc,nspin))
 if(.not.associated(SinvHDaux_red)) allocate(SinvHDaux_red(nuotot,nuotot,nspin))
#ifdef MPI
 if(.not.Serial) then
   nprow = 1
   npcol = Nodes
   if(FirstCall) then
     call blacs_get( -1, 0, ictxt )
     call blacs_gridinit( ictxt, 'C', nprow, npcol )
   endif
   FirstCall = .false.
!
   call descinit( desch, nuotot, nuotot, BlockSize, BlockSize, 0, 0, ictxt, nuotot, info )
!
 endif
#endif
!
 if(associated(Daux_old)) call de_alloc(Daux_old, name='Daux_old')
 call re_alloc(Daux_old, 1, nuotot, 1, nuoloc, 1, nspin, name='Daux_old', routine='stabiliz')
!
 newS = .false.
!
!population = 0.0_dp
!
 abscommutdiff = abstolstb
 if(Node.eq.0) then
   write(6,'(/a)') 'state_stabiliz: convergence criteria'
   if(.not.pseudoatm) write(6,'(a20)',advance='no') '   S^(-1)HD-DHS^(-1)'
   write(6,*)
   if(.not.pseudoatm) write(6,'(3x,e15.7,2x)',advance='no') abscommutdiff 
   write(6,*)
   call flush(6)
 endif
!
#ifdef MPI
 call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
!
 cnvg = .false.
 istb = 1
!
 do while((.not.cnvg) .and. (istb.le.maxnstb))
! Store old density matrix
    do ispin = 1,nspin
       do io = 1,nuoloc
          do jo = 1,nuotot
             Daux_old(jo,io,ispin) = Daux(jo,io,ispin)
          enddo
       enddo
       do jo = 1,size(Dscf,1)
          Dold(jo,ispin) = Dscf(jo,ispin)
       enddo
    enddo
! Construct Hamiltonian
    call sauxdaux2haux(istb, newS, Saux, Daux, Haux)
! Diagonalize Hamiltonian
    do ispin = 1,nspin
       do io = 1,nuoloc
          do jo = 1,nuotot
             Faux(jo,io) = Haux(jo,io,ispin)
             if(ispin.eq.1) &
             Gaux(jo,io) = Saux(jo,io)
          enddo
       enddo
       call rdiag(Faux,Gaux,nuotot,nuoloc,nuotot,eig_new(1,ispin),wfk_new(1,1,ispin),nuotot,istb,ierror)
    enddo
! Calculate commutor of S^(-1)HD-DHS^(-1)
    do ispin = 1,nspin
       HDaux(:,:,ispin) = 0.0_dp
       if(Serial) then
         call dgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Haux(:,:,ispin),nuotot, Daux(:,:,ispin),nuotot,0.0_dp, &
                        HDaux(:,:,ispin),nuotot)
         call dgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux_inv,       nuotot,HDaux(:,:,ispin),nuotot,0.0_dp, &
                    SinvHDaux(:,:,ispin),nuotot)
#ifdef MPI
       else
         call pdgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Haux(:,:,ispin),1,1,desch, Daux(:,:,ispin),1,1,desch,0.0_dp, &
                         HDaux(:,:,ispin),1,1,desch)
         call pdgemm('N','N',nuotot,nuotot,nuotot,1.0_dp,Saux_inv,       1,1,desch,HDaux(:,:,ispin),1,1,desch,0.0_dp, &
                     SinvHDaux(:,:,ispin),1,1,desch)
#endif
       endif
    enddo
    call collectdata_dp(nuotot, nuoloc, nspin, SinvHDaux, SinvHDaux_red)
! Calculate new density matrix
    call psi2daux(wfk_new, Daux)
    call daux2dsparse(Daux)
!
    do ispin = 1,nspin
       abscommutdiff(ispin) = sum(dabs(SinvHDaux_red(:,:,ispin)-transpose(SinvHDaux_red(:,:,ispin))))
    enddo
!
    if(pseudoatm) abscommutdiff = 0.0_dp
!
    cnvg = (sum(abscommutdiff) .lt. abstolstb)
!
    dDmax_tmp = 0.0_dp
    dDmax_tmp(Node+1) = maxval(dabs(Dscf-Dold))
#ifdef MPI
    call MPI_AllReduce(dDmax_tmp(1), dDmax_red(1), Nodes, &
         MPI_double_precision, MPI_sum, MPI_Comm_World, MPIerror)
#else
    dDmax_red = dDmax_tmp
#endif
    dDmax = maxval(dDmax_red)
!
    if(chknormdm) call normdaux(Daux, TrSDaux)
!
    if(Node.eq.0) then
      if(istb.eq.1) then
        write(6,'(/,a12)',advance='no') 'siesta: istb'
        if(.not.pseudoatm) then
         write(6, '(a20)',advance='no') '   S^(-1)HD-DHS^(-1)'
         write(6, '(a20)',advance='no') '         dDmax      '
         if(chknormdm) &
         write(6, '(a20)',advance='no') '        Tr[SD]      '
        endif
        write(6,*)
        call flush(6)
      endif
      write(6,'(a8,i4)',advance='no')   'siesta: ', istb
      if(.not.pseudoatm) then
       write(6,'(3x,e15.7,2x)',advance='no')  sum(abscommutdiff)
       write(6,'(3x,e15.7,2x)',advance='no')              dDmax
       if(chknormdm) &
       write(6,'(3x,e15.7,2x)',advance='no')     sum(TrSDaux(:))
      endif
      write(6,*)
      call flush(6)
    endif
#ifdef MPI
    call MPI_Barrier(MPI_Comm_World,MPIerror)
#endif
! Mixing density matrix
    if(.not.cnvg) then
      if(.not.pseudoatm) &
       call daux_pulay_mixer(istb, maxsav, wmix, nkick, wmixkick, Daux_old, Daux)
    endif
!
    istb = istb + 1
!
!   call lowdinpop(nuotot, nuoloc, nspin, Daux, natot, population)
!   call mullikpop(nuotot, nuoloc, nspin, Daux, natot, population)
 enddo
! For debug
 if(.not.pseudonfr) call outnfr('STB')
! For debug
!if(chkstabl) then
!  if(Node.eq.0) then
!    open(unit=90,file='check_eig.ascii')
!    do ispin = 1,nspin
!       do io = 1,nuotot
!          write(90,'(3e15.7)') eig(io,ispin),  eig_new(io,ispin), &
!                               eig(io,ispin) - eig_new(io,ispin)
!       enddo
!    enddo
!    close(90)
!  endif
!
!  call collectdata_dp(nuotot, nuoloc, nspin, wfk_new, wfk_new_red)
!  if(Node.eq.0) then
!    open(unit=91,file='check_psi.ascii')
!    do ispin = 1,nspin
!       do jo = 1,nuotot
!          do io = 1,nuotot
!             write(91,'(3e15.7)') wfk_red(jo,io,ispin),       wfk_new_red(jo,io,ispin), &
!                              abs(wfk_red(jo,io,ispin)) - abs(wfk_new_red(jo,io,ispin))
!          enddo
!       enddo
!    enddo
!    close(91)
!  endif
!
!  call de_alloc(wfk_red, name='wfk_red')
!endif
!
 do ispin = 1,nspin
    do io = 1,nuotot
       eig(io,ispin) = eig_new(io,ispin)
    enddo
    do io = 1,nuoloc
       do jo = 1,nuotot
          wfk(jo,io,ispin) = wfk_new(jo,io,ispin)
       enddo
    enddo
 enddo
!
 call finalize_daux_pulay()
!
 deallocate(eig_new)       ; nullify(eig_new)
 deallocate(Faux)          ; nullify(Faux)
 deallocate(Gaux)          ; nullify(Gaux)
 deallocate(wfk_new)       ; nullify(wfk_new)
 deallocate(HDaux)         ; nullify(HDaux)
 deallocate(SinvHDaux)     ; nullify(SinvHDaux)
 deallocate(SinvHDaux_red) ; nullify(SinvHDaux_red)
!
 call timer( 'stabiliz', 2 )
!
 end subroutine state_stabiliz
