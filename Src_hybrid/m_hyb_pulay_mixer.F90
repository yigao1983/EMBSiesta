 module m_hyb_pulay_mixer
 use precision,      only : dp
 use siesta_options, only : avoid_first_after_kick
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 integer, save :: n_records_dft
 integer, save :: n_records_nfr
 real(dp), pointer, save ::  dauxmem(:,:,:,:)=>null()
 real(dp), pointer, save ::  rhonmem(:,:,:,:)=>null()
 real(dp), pointer, save :: ddauxmem(:,:,:,:)=>null()
 real(dp), pointer, save :: drhonmem(:,:,:,:)=>null()
!
 contains
!
 subroutine initialz_daux_pulay(npulaymax)
 use alloc,       only : re_alloc, de_alloc
 use m_spin,      only : nspin
 use m_variables, only : nuotot, nuoloc
 implicit none
 integer, intent(in) :: npulaymax
!
 n_records_dft = 0
!
 if(associated( dauxmem)) then
   call de_alloc( dauxmem, name= 'dauxmem')
   nullify( dauxmem)
 endif
 if(associated(ddauxmem)) then
   call de_alloc(ddauxmem, name='ddauxmem')
   nullify(ddauxmem)
 endif
 call re_alloc( dauxmem, 1, nuotot, 1, nuoloc, 1, nspin, 1, npulaymax, &
               name= 'dauxmem', routine='init_hyb_pulay')
 call re_alloc(ddauxmem, 1, nuotot, 1, nuoloc, 1, nspin, 1, npulaymax, &
               name='ddauxmem', routine='init_hyb_pulay')
!
 end subroutine initialz_daux_pulay
!
 subroutine initialz_rhonfr_pulay(npulaymax)
 use alloc,       only : re_alloc, de_alloc
 use m_variables, only : pseudonfr
 use m_variables, only : Nxl, Nyl, Nzl
 implicit none
 integer, intent(in) :: npulaymax
!
 n_records_nfr = 0
!
 if(.not.pseudonfr) then
   if(associated( rhonmem)) then
     call de_alloc( rhonmem, name= 'rhonmem')
     nullify( rhonmem)
   endif
   if(associated(drhonmem)) then
     call de_alloc(drhonmem, name='drhonmem')
     nullify(drhonmem)
   endif
   call re_alloc( rhonmem, 1, Nxl, 1, Nyl, 1, Nzl, 1, npulaymax, &
                  name= 'rhonmem', routine='init_hyb_pulay')
   call re_alloc(drhonmem, 1, Nxl, 1, Nyl, 1, Nzl, 1, npulaymax, &
                  name='drhonmem', routine='init_hyb_pulay')
 endif
!
 end subroutine initialz_rhonfr_pulay
!
 subroutine finalize_daux_pulay()
 use alloc, only : de_alloc
 implicit none
!
 if(associated( dauxmem)) call de_alloc( dauxmem, name= 'dauxmem')
 if(associated(ddauxmem)) call de_alloc(ddauxmem, name='ddauxmem')
!
 end subroutine finalize_daux_pulay
!
 subroutine finalize_rhonfr_pulay()
 use alloc, only : de_alloc
 implicit none
!
 if(associated( rhonmem)) call de_alloc( rhonmem, name= 'rhonmem')
 if(associated(drhonmem)) call de_alloc(drhonmem, name='drhonmem')
!
 end subroutine finalize_rhonfr_pulay
!
 subroutine daux_pulay_mixer(iscf, npulaymax, alpha, nkick, alphakick, Daux_old, Daux_new)
 use alloc,       only : re_alloc, de_alloc
 use parallel,    only : Node
 use m_spin,      only : nspin
 use m_variables, only : nuotot, nuoloc
 implicit none
 integer,  intent(in)    :: iscf
 integer,  intent(in)    :: npulaymax, nkick
 real(dp), intent(in)    :: alpha, alphakick
 real(dp), intent(inout) :: Daux_old(nuotot,nuoloc,nspin)
 real(dp), intent(inout) :: Daux_new(nuotot,nuoloc,nspin)
! Local variables
 integer :: idim, jdim, ipulay, ispin, io, jo, ndim
 integer :: nrecords
 logical :: after_kick
 real(dp) :: sumddaux, sumddaux_red
 real(dp), dimension(:),   pointer :: coef=>null()
 real(dp), dimension(:,:), pointer :: beta=>null()
#ifdef MPI
 integer :: MPIerror
#endif
!
 if(iscf.eq.1) call initialz_daux_pulay(npulaymax)
!
 after_kick = .false.
 if(iscf.eq.1) then
   after_kick = .true.
 elseif(nkick.gt.0) then
   if(mod(iscf,nkick).eq.1) after_kick = .true.
 endif
!
 if(.not.after_kick .or. .not.avoid_first_after_kick) then
   n_records_dft = n_records_dft + 1
   if(n_records_dft .le. npulaymax) then
      dauxmem(:,:,:,n_records_dft) = Daux_old
     ddauxmem(:,:,:,n_records_dft) = Daux_new - Daux_old
   else
     do ipulay = 1,npulaymax-1
         dauxmem(:,:,:,ipulay) =  dauxmem(:,:,:,ipulay+1)
        ddauxmem(:,:,:,ipulay) = ddauxmem(:,:,:,ipulay+1)
     enddo
      dauxmem(:,:,:,npulaymax) = Daux_old
     ddauxmem(:,:,:,npulaymax) = Daux_new - Daux_old
   endif
 endif
!
 if(n_records_dft.le.1) then
   Daux_new = alpha * Daux_new + (1.0_dp-alpha) * Daux_old
   return
 endif
!
 if(nkick.gt.0 .and. mod(iscf,nkick).eq.0) then
   n_records_dft = 0
   Daux_new = alphakick * Daux_new + (1.0_dp-alphakick) * Daux_old
   return
 endif
!
 ndim = min(n_records_dft, npulaymax) + 1
!
 if(associated(beta)) then
   call de_alloc(beta, name='beta')
   nullify(beta)
 endif
 if(associated(coef)) then
   call de_alloc(coef, name='coef')
   nullify(coef)
 endif
 call re_alloc(beta, 1, ndim, 1, ndim, name='beta', routine='hybrid_pulay_mixer')
 call re_alloc(coef, 1, ndim, name='coef', routine='hybrid_pulay_mixer')
!
 do idim = 1,ndim-1
    do jdim = 1,idim
       sumddaux = 0.0_dp
       do ispin = 1,nspin
          do io = 1,nuoloc
             do jo = 1,nuotot
                sumddaux = sumddaux + ddauxmem(jo,io,ispin,jdim) * ddauxmem(jo,io,ispin,idim)
             enddo
          enddo
       enddo
#ifdef MPI
       call MPI_AllReduce(sumddaux, sumddaux_red, 1, MPI_double_precision, &
            MPI_sum, MPI_Comm_World, MPIerror)
#else
       sumddaux_red = sumddaux
#endif
       beta(jdim,idim) = sumddaux_red
       beta(idim,jdim) = sumddaux_red
    enddo
 enddo
!
 beta(1:ndim-1,ndim) = 1.0_dp
 beta(ndim,1:ndim-1) = 1.0_dp
 beta(ndim,ndim)     = 0.0_dp
!
 coef(1:ndim-1) = 0.0_dp
 coef(ndim)     = 1.0_dp
!
 call solve_symlineq(ndim, beta, coef)
!
 Daux_old = 0.0_dp
 Daux_new = 0.0_dp
 do idim = 1,ndim-1
    Daux_old = Daux_old + coef(idim) * dauxmem(:,:,:,idim)
    Daux_new = Daux_new + coef(idim) *(dauxmem(:,:,:,idim) + ddauxmem(:,:,:,idim))
 enddo
!
 Daux_new = alpha * Daux_new + (1.0_dp-alpha) * Daux_old
!
 call de_alloc(beta, name='beta')
 call de_alloc(coef, name='coef')
!
 end subroutine daux_pulay_mixer
!
 subroutine rhonfr_pulay_mixer(iscf, npulaymax, alpha, nkick, alphakick, Rhonfr_old, Rhonfr_new)
 use alloc,       only : re_alloc, de_alloc
 use parallel,    only : Node
 use mesh,        only : nsm
 use m_variables, only : Nxl, Nyl, Nzl, Nx, Ny, Nz
 use m_numeric3d, only : integration
 implicit none
 integer,  intent(in)    :: iscf
 integer,  intent(in)    :: npulaymax, nkick
 real(dp), intent(in)    :: alpha, alphakick
 real(dp), intent(inout) :: Rhonfr_old(Nxl,Nyl,Nzl)
 real(dp), intent(inout) :: Rhonfr_new(Nxl,Nyl,Nzl)
! Local variables
 integer  :: idim, jdim, ipulay, ispin, io, jo, ndim
 integer  :: nmesh(3), nmeshl(3)
 logical  :: after_kick
 real(dp) :: sumdrhonfr
 real(dp), dimension(:),   pointer :: coef=>null()
 real(dp), dimension(:,:), pointer :: beta=>null()
#ifdef MPI
 integer :: MPIerror
#endif
!
 if(iscf.eq.1) call initialz_rhonfr_pulay(npulaymax)
!
 nmesh(1) = Nx  ; nmesh(2) = Ny  ; nmesh(3) = Nz
 nmeshl(1)= Nxl ; nmeshl(2)= Nyl ; nmeshl(3)= Nzl
!
 after_kick = .false.
 if(iscf.eq.1) then
   after_kick = .true.
 elseif(nkick.gt.0) then
   if(mod(iscf,nkick).eq.1) after_kick = .true.
 endif
!
 if(.not.after_kick .or. .not.avoid_first_after_kick) then
   n_records_nfr = n_records_nfr + 1
   if(n_records_nfr .le. npulaymax) then
      rhonmem(:,:,:,n_records_nfr) = Rhonfr_old
     drhonmem(:,:,:,n_records_nfr) = Rhonfr_new - Rhonfr_old
   else
     do ipulay = 1,npulaymax-1
         rhonmem(:,:,:,ipulay) =  rhonmem(:,:,:,ipulay+1)
        drhonmem(:,:,:,ipulay) = drhonmem(:,:,:,ipulay+1)
     enddo
      rhonmem(:,:,:,npulaymax) = Rhonfr_old
     drhonmem(:,:,:,npulaymax) = Rhonfr_new - Rhonfr_old
   endif
 endif
!
 if(n_records_nfr.le.1) then
   Rhonfr_new = alpha * Rhonfr_new + (1.0_dp-alpha) * Rhonfr_old
   return
 endif
!
 if(nkick.gt.0 .and. mod(iscf,nkick).eq.0) then
   n_records_nfr = 0
   Rhonfr_new = alphakick * Rhonfr_new + (1.0_dp-alphakick) * Rhonfr_old
   return
 endif
!
 ndim = min(n_records_nfr, npulaymax) + 1
!
 if(associated(beta)) then
   call de_alloc(beta, name='beta')
   nullify(beta)
 endif
 if(associated(coef)) then
   call de_alloc(coef, name='coef')
   nullify(coef)
 endif
 call re_alloc(beta, 1, ndim, 1, ndim, name='beta', routine='hybrid_pulay_mixer')
 call re_alloc(coef, 1, ndim, name='coef', routine='hybrid_pulay_mixer')
!
 do idim = 1,ndim-1
    do jdim = 1,idim
       call integration(drhonmem(:,:,:,jdim)*drhonmem(:,:,:,idim), nmesh, nmeshl, nsm, sumdrhonfr)
       beta(jdim,idim) = sumdrhonfr
       beta(idim,jdim) = sumdrhonfr
    enddo
 enddo
!
 beta(1:ndim-1,ndim) = 1.0_dp
 beta(ndim,1:ndim-1) = 1.0_dp
 beta(ndim,ndim)     = 0.0_dp
!
 coef(1:ndim-1) = 0.0_dp
 coef(ndim)     = 1.0_dp
!
 call solve_symlineq(ndim, beta, coef)
!
 Rhonfr_old = 0.0_dp
 Rhonfr_new = 0.0_dp
 do idim = 1,ndim-1
    Rhonfr_old = Rhonfr_old + coef(idim) * rhonmem(:,:,:,idim)
    Rhonfr_new = Rhonfr_new + coef(idim) *(rhonmem(:,:,:,idim) + drhonmem(:,:,:,idim))
 enddo
!
 Rhonfr_new = alpha * Rhonfr_new + (1.0_dp-alpha) * Rhonfr_old
!
 call de_alloc(beta, name='beta') ; nullify(beta)
 call de_alloc(coef, name='coef') ; nullify(coef)
!
 end subroutine rhonfr_pulay_mixer
!
 subroutine solve_symlineq(n, a, b)
 use alloc, only : re_alloc, de_alloc
 implicit none
 integer,  intent(in)    :: n
 real(dp), intent(in)    :: a(n,n)
 real(dp), intent(inout) :: b(n)
! Local variables
 integer   :: nrhs, lda, ldb, lwork, info
 integer   :: ipiv(n)
 character :: uplo
 real(dp), pointer :: work(:)
!
 uplo = 'U'
 nrhs = 1
 lda  = n
 ldb  = n
 lwork= n
!
 call re_alloc(work, 1, lwork, name='work', routine='solve_symlineq')
!
 call dsysv( uplo, n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info )
!
 call de_alloc(work, name='work')
!
 end subroutine solve_symlineq
!
 end module m_hyb_pulay_mixer
