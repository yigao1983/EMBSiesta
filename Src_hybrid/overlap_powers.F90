 subroutine overlap_powers
 use precision,    only : dp
 use alloc,        only : re_alloc, de_alloc
 use sys,          only : die
 use parallel,     only : Node, Nodes, BlockSize
 use parallelsubs, only : LocalToGlobalOrb
 use m_variables,  only : simplempi
 use m_variables,  only : nuotot, nuoloc
 use m_variables,  only : Saux, Saux_inv, Saux_nhf, Saux_phf
 use m_variables,  only : Saux_nhf_red, Saux_phf_red
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
! Local variables
 integer :: io, jo, iio, iscf, ierror, info
 logical :: Serial
 logical, save :: FirstCall = .true.
!real(dp), dimension(nuotot)        :: w
!real(dp), dimension(nuotot,nuoloc) :: Faux, Iaux
!real(dp), dimension(nuotot,nuoloc) :: z, zw
 real(dp), dimension(:),   pointer, save :: w=>null()
 real(dp), dimension(:,:), pointer, save :: Faux=>null(), Iaux=>null()
 real(dp), dimension(:,:), pointer, save :: z=>null(), zw=>null()
#ifdef MPI
 integer :: MPIerror
 integer :: nprow, npcol
 integer :: descs(9)
 integer, save :: ictxt
#endif
!
 call timer( 'spowers', 1 )
!
 Serial = (Nodes.eq.1)
! Memory of output: Saux_inv, Saux_phf, Saux_nhf
 if(.not.associated(Saux_inv)) &
  call re_alloc(Saux_inv,1,nuotot,1,nuoloc,name='Saux_inv',routine='overlap_powers')
 if(.not.associated(Saux_nhf)) &
  call re_alloc(Saux_nhf,1,nuotot,1,nuoloc,name='Saux_nhf',routine='overlap_powers')
 if(.not.associated(Saux_phf)) &
  call re_alloc(Saux_phf,1,nuotot,1,nuoloc,name='Saux_phf',routine='overlap_powers')
!
 if(FirstCall) then
   if(.not.associated(w))    allocate(   w(nuotot))
   if(.not.associated(Faux)) allocate(Faux(nuotot,nuoloc))
   if(.not.associated(Iaux)) allocate(Iaux(nuotot,nuoloc))
   if(.not.associated(z))    allocate(   z(nuotot,nuoloc))
   if(.not.associated(zw))   allocate(  zw(nuotot,nuoloc))
 endif
! Identity matrix initialization
 Iaux = 0.0_dp
 do io = 1,nuoloc
    call LocalToGlobalOrb(io,Node,Nodes,iio)
    Iaux(iio,io) = 1.0_dp
 enddo 
! Auxiliary variables initialization
 w(:)   = 0.0_dp
 z(:,:) = 0.0_dp
! Diagonalization of Saux
 iscf = 1
!do io = 1,nuoloc
!   do jo = 1,nuotot
!      Faux(jo,io) = Saux(jo,io)
!   enddo
!enddo
 Faux(:,:) = Saux(:,:)
 call rdiag(Faux,Iaux,nuotot,nuoloc,nuotot,w,z,nuotot,iscf,ierror)
!
 if(ierror.ne.0) call die('Terminating due to failed diagonalisation')
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
   call descinit( descs, nuotot, nuotot, BlockSize, BlockSize, 0, 0, ictxt, nuotot, info )
!
 endif
#endif
! Calculate Saux_inv = S^(-1)
 do io = 1,nuoloc
    call LocalToGlobalOrb(io,Node,Nodes,iio)
    zw(:,io) = z(:,io) / w(iio)
 end do
!
 Saux_inv = 0.0_dp
 if(Serial) then
   call dgemm('N','T',nuotot,nuotot,nuotot,1.0_dp,zw,nuotot,z,nuotot,0.0_dp,Saux_inv,nuotot)
#ifdef MPI
 else
   call pdgemm('N','T',nuotot,nuotot,nuotot,1.0_dp,zw,1,1,descs,z,1,1,descs,0.0_dp,Saux_inv,1,1,descs)
#endif
 endif
! Calculate Saux_nhf = S^(-1/2)
 do io = 1,nuoloc
    call LocalToGlobalOrb(io,Node,Nodes,iio)
    zw(:,io) = z(:,io) / dsqrt(w(iio))
 end do
!
 Saux_nhf = 0.0_dp
 if(Serial) then
   call dgemm('N','T',nuotot,nuotot,nuotot,1.0_dp,zw,nuotot,z,nuotot,0.0_dp,Saux_nhf,nuotot)
#ifdef MPI
 else
   call pdgemm('N','T',nuotot,nuotot,nuotot,1.0_dp,zw,1,1,descs,z,1,1,descs,0.0_dp,Saux_nhf,1,1,descs)
#endif
 endif
! Calculate Saux_phf = S^(1/2)
 do io = 1,nuoloc
    call LocalToGlobalOrb(io,Node,Nodes,iio)
    zw(:,io) = z(:,io) * dsqrt(w(iio))
 end do
!
 Saux_phf = 0.0_dp
 if(Serial) then
   call dgemm('N','T',nuotot,nuotot,nuotot,1.0_dp,zw,nuotot,z,nuotot,0.0_dp,Saux_phf,nuotot)
#ifdef MPI
 else
   call pdgemm('N','T',nuotot,nuotot,nuotot,1.0_dp,zw,1,1,descs,z,1,1,descs,0.0_dp,Saux_phf,1,1,descs)
#endif
 endif
!
!if(.not.associated(Saux_red))     &
! call re_alloc(Saux_red,    1,nuotot,1,nuotot,name='Saux_red',    routine='spowers')
!if(.not.associated(Saux_inv_red)) &
! call re_alloc(Saux_inv_red,1,nuotot,1,nuotot,name='Saux_inv_red',routine='spowers')
!if(.not.associated(Saux_phf_red)) &
! call re_alloc(Saux_phf_red,1,nuotot,1,nuotot,name='Saux_phf_red',routine='spowers')
!call collectdata_dp(nuotot, nuoloc, 1, Saux, Saux_red)
!call collectdata_dp(nuotot, nuoloc, 1, Saux_inv, Saux_inv_red)
!call collectdata_dp(nuotot, nuoloc, 1, Saux_phf, Saux_phf_red)
!
 if(simplempi) then
   if(.not.associated(Saux_nhf_red)) &
   call re_alloc(Saux_nhf_red,1,nuotot,1,nuotot,name='Saux_nhf_red',routine='spowers')
   call collectdata_dp(nuotot, nuoloc, 1, Saux_nhf, Saux_nhf_red)
   if(.not.associated(Saux_phf_red)) &
   call re_alloc(Saux_phf_red,1,nuotot,1,nuotot,name='Saux_phf_red',routine='spowers')
   call collectdata_dp(nuotot, nuoloc, 1, Saux_phf, Saux_phf_red)
 endif
!
 FirstCall = .false.
!
 call timer( 'spowers', 2 )
!
 end subroutine overlap_powers
