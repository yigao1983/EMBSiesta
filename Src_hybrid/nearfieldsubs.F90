 module nearfieldsubs
 use precision, only : dp, grid_p
 use sys,       only : die
 use alloc,     only : re_alloc, de_alloc
#ifdef MPI
 use mpi_siesta
#endif
 implicit none
 real(dp), parameter :: abstolnfr = 1.0e-12_dp
!
 contains
!
 subroutine GetGridTables(nmeshl)
 use parallel,    only : ProcessorY, Node, Nodes
 use m_variables, only : GridTableY, GridTableZ
 implicit none
 integer, intent(in)  :: nmeshl(3)
! Local variables
 integer :: ProcessorZ, Py, Pz
 integer :: YCommunicator, ZCommunicator, YNode, YNodes, ZNode, ZNodes
 integer, pointer :: GridTableY_tmp(:), GridTableZ_tmp(:)
#ifdef MPI
 integer :: MPIerror
#endif
!
 ProcessorZ = Nodes/ProcessorY
 if(ProcessorY*ProcessorZ.ne.Nodes) &
  call die('ERROR: ProcessorY must be a factor of the number of processors!')
!
 Py = (Node/ProcessorZ) + 1
 Pz = Node - (Py - 1)*ProcessorZ + 1
!
 if(associated(GridTableY)) call de_alloc(GridTableY, name='GridTableY')
 call re_alloc(GridTableY, 1, ProcessorY, name='GridTableY', routine='GetGridTables')
 if(associated(GridTableZ)) call de_alloc(GridTableZ, name='GridTableZ')
 call re_alloc(GridTableZ, 1, ProcessorZ, name='GridTableZ', routine='GetGridTables')
!
 if(associated(GridTableY_tmp)) call de_alloc(GridTableY_tmp, name='GridTableY_tmp')
 call re_alloc(GridTableY_tmp, 1, ProcessorY, name='GridTableY_tmp', routine='GetGridTables')
 if(associated(GridTableZ_tmp)) call de_alloc(GridTableZ_tmp, name='GridTableZ_tmp')
 call re_alloc(GridTableZ_tmp, 1, ProcessorZ, name='GridTableZ_tmp', routine='GetGridTables')
!
!  Group processors into subsets by Z and Y
!
#ifdef MPI
 call MPI_Comm_Split(MPI_Comm_World, Pz, Py, ZCommunicator, MPIerror)
 call MPI_Comm_Rank(ZCommunicator,ZNode, MPIerror)
 call MPI_Comm_Size(ZCommunicator,ZNodes,MPIerror)
!
 call MPI_Comm_Split(MPI_Comm_World, Py, Pz, YCommunicator, MPIerror)
 call MPI_Comm_Rank(YCommunicator,YNode, MPIerror)
 call MPI_Comm_Size(YCommunicator,YNodes,MPIerror)
#endif
!
 GridTableY_tmp = 0
 GridTableY_tmp(Py) = nmeshl(2)
#ifdef MPI
 call MPI_AllReduce(GridTableY_tmp(1),GridTableY(1),ProcessorY,MPI_integer,MPI_Sum,ZCommunicator,MPIerror)
#else
 GridTableY = GridTableY_tmp
#endif
!
 GridTableZ_tmp = 0
 GridTableZ_tmp(Pz) = nmeshl(3)
#ifdef MPI
 call MPI_AllReduce(GridTableZ_tmp(1),GridTableZ(1),ProcessorZ,MPI_integer,MPI_Sum,YCommunicator,MPIerror)
#else
 GridTableZ = GridTableZ_tmp
#endif
!
 call de_alloc(GridTableY_tmp, name='GridTableY_tmp')
 call de_alloc(GridTableZ_tmp, name='GridTableZ_tmp')
!
 end subroutine GetGridTables
!
 subroutine LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
 use parallel,    only : ProcessorY
 use m_variables, only : Nyl, Nzl, GridTableY, GridTableZ
 implicit none
 integer, intent(in)  :: Node, Nodes
 integer, intent(in)  :: iix, iiy, iiz
 integer, intent(out) :: ix, iy, iz
! Local variable
 integer :: ProcessorZ, Py, Pz
!
 ProcessorZ = Nodes/ProcessorY
 if(ProcessorY*ProcessorZ.ne.Nodes) &
  call die('ERROR: ProcessorY must be a factor of the number of processors!')
!
 Py = (Node/ProcessorZ) + 1
 Pz = Node - (Py - 1)*ProcessorZ + 1
! X direction
 ix = iix
! Y direction
 iy = sum(GridTableY(1:Py-1)) + iiy
! Z direction
 iz = sum(GridTableZ(1:Pz-1)) + iiz
!
 end subroutine LocalToGlobalGrid
!
 subroutine InitNFmem(ngrid, ngloc)
 use m_variables, only : Nosc, eps_vac
 use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
 use m_variables, only : xbox, ybox, zbox
 use m_variables, only : eps, alph, beta, omeg
 use m_variables, only : phi, efd
 use m_variables, only : pol, cur, poltot, curtot
 use m_variables, only : rhohyb, rhodft, rhonfr, rhonfr_old
 use m_variables, only : TRho
 implicit none
 integer, intent(in) :: ngrid(3), ngloc(3)
! Local variable
 integer :: ierr
!
 Nx = ngrid(1) ; Ny = ngrid(2) ; Nz = ngrid(3)
 Nxl= ngloc(1) ; Nyl= ngloc(2) ; Nzl= ngloc(3)
!
 if(associated(xbox)) call de_alloc(xbox,name='xbox')
 call re_alloc(xbox,1,Nx,name='xbox',routine='initnfmem')
 if(associated(ybox)) call de_alloc(ybox,name='ybox')
 call re_alloc(ybox,1,Ny,name='ybox',routine='initnfmem')
 if(associated(zbox)) call de_alloc(zbox,name='zbox')
 call re_alloc(zbox,1,Nz,name='zbox',routine='initnfmem')
!
 xbox = 0.0_dp ; ybox = 0.0_dp ; zbox = 0.0_dp
!
 if(associated(eps)) call de_alloc(eps,name='eps')
 call re_alloc(eps,1,Nxl,1,Nyl,1,Nzl,name='eps',routine='initnfmem')
!
 eps = eps_vac
!
 if(associated(alph)) call de_alloc(alph,name='alph')
 call re_alloc(alph,1,Nxl,1,Nyl,1,Nzl,1,Nosc,name='alph',routine='initnfmem')
 if(associated(beta)) call de_alloc(beta,name='beta')
 call re_alloc(beta,1,Nxl,1,Nyl,1,Nzl,1,Nosc,name='beta',routine='initnfmem')
 if(associated(omeg)) call de_alloc(omeg,name='omeg')
 call re_alloc(omeg,1,Nxl,1,Nyl,1,Nzl,1,Nosc,name='omeg',routine='initnfmem')
!
 alph = 0.0_dp ; beta = 0.0_dp ; omeg = 0.0_dp
!
 if(associated(phi)) call de_alloc(phi,name='phi')
 call re_alloc(phi,1,Nxl,1,Nyl,1,Nzl,name='phi',routine='initnfmem')
!
 phi = 0.0_dp
!
 if(associated(rhohyb)) call de_alloc(rhohyb,name='rhohyb')
 call re_alloc(rhohyb,1,Nxl,1,Nyl,1,Nzl,name='rhohyb',routine='initnfmem')
 if(associated(rhodft)) call de_alloc(rhodft,name='rhodft')
 call re_alloc(rhodft,1,Nxl,1,Nyl,1,Nzl,name='rhodft',routine='initnfmem')
 if(associated(rhonfr)) call de_alloc(rhonfr,name='rhonfr')
 call re_alloc(rhonfr,1,Nxl,1,Nyl,1,Nzl,name='rhonfr',routine='initnfmem')
 if(associated(rhonfr_old)) call de_alloc(rhonfr_old,name='rhonfr_old')
 call re_alloc(rhonfr_old,1,Nxl,1,Nyl,1,Nzl,name='rhonfr_old',routine='initnfmem')
!
 rhohyb = 0.0_dp
 rhodft = 0.0_dp
 rhonfr = 0.0_dp
 rhonfr_old = 0.0_dp
!
 if(associated(efd)) call de_alloc(efd,name='efd')
 call re_alloc(efd,1,Nxl,1,Nyl,1,Nzl,1,3,name='efd',routine='initnfmem')
!
 efd = 0.0_dp
!
 if(associated(poltot)) call de_alloc(poltot,name='poltot')
 call re_alloc(poltot,1,Nxl,1,Nyl,1,Nzl,1,3,name='poltot',routine='initnfmem')
 if(associated(curtot)) call de_alloc(curtot,name='curtot')
 call re_alloc(curtot,1,Nxl,1,Nyl,1,Nzl,1,3,name='curtot',routine='initnfmem')
!
 poltot = 0.0_dp ; curtot = 0.0_dp
!
 if(associated(pol)) deallocate(pol)
 allocate(pol(Nxl,Nyl,Nzl,3,Nosc), stat=ierr)
 if(associated(cur)) deallocate(cur)
 allocate(cur(Nxl,Nyl,Nzl,3,Nosc), stat=ierr)
!
 pol = 0.0_dp ; cur = 0.0_dp
!
 if(associated(TRho)) call de_alloc(TRho,name='TRho')
 call re_alloc(TRho,1,Nxl*Nyl*Nzl,name='TRho',routine='initnfmem')
!
 TRho = 0.0_grid_p
!
 end subroutine InitNFmem
!
 subroutine ConfigMater
 use m_variables, only : eps_Au, alph_Au, beta_Au, omeg_Au
 use m_variables, only : eps_Ag, alph_Ag, beta_Ag, omeg_Ag
 use m_variables, only : eps_drude, alph_drude, beta_drude, omeg_drude
 use m_variables, only : eps_vac, alph_vac, beta_vac, omeg_vac
 use m_variables, only : Nobj, nfobj
 implicit none
! Local variable
 integer :: iobj
!
 do iobj = 1,Nobj
    select case(nfobj(iobj)%material)
    case('Au') 
        nfobj(iobj)%eps  =  eps_Au
        nfobj(iobj)%alph = alph_Au
        nfobj(iobj)%beta = beta_Au
        nfobj(iobj)%omeg = omeg_Au
    case('Ag')
        nfobj(iobj)%eps  =  eps_Ag
        nfobj(iobj)%alph = alph_Ag
        nfobj(iobj)%beta = beta_Ag
        nfobj(iobj)%omeg = omeg_Ag
    case('drude')
        nfobj(iobj)%eps  =  eps_drude
        nfobj(iobj)%alph = alph_drude
        nfobj(iobj)%beta = beta_drude
        nfobj(iobj)%omeg = omeg_drude
    case('vac')
        nfobj(iobj)%eps  =  eps_vac
        nfobj(iobj)%alph = alph_vac
        nfobj(iobj)%beta = beta_vac
        nfobj(iobj)%omeg = omeg_vac
    endselect
 enddo
!
 end subroutine ConfigMater
!
 subroutine ConfigSpace(ngrid, ngloc, cell)
 use parallel,    only : Node, Nodes
 use m_variables, only : Nosc
 use m_variables, only : Nobj, nfobj
 use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
 use m_variables, only : xbox, ybox, zbox, lbox, dbox
 use m_variables, only : eps_vac, alph_vac, beta_vac, omeg_vac
 use m_variables, only : eps, alph, beta, omeg
 implicit none
 integer,  intent(in) :: ngrid(3), ngloc(3)
 real(dp), intent(in) :: cell(3,3)
! Local variables
 integer  :: iix, iiy, iiz, idim, ix, iy, iz, iobj, iosc
 integer  :: orient, smrdir
 real(dp) :: dx, dy, dz
 real(dp) :: xc, yc, zc, radius, iradius, oradius, length, lengthx, lengthy, lengthz, smear
 logical  :: judge, affected
 real(dp) :: epsl
 real(dp), dimension(Nosc) :: alphl, betal, omegl
 real(dp), dimension(Nosc) :: betas
 real(dp), pointer :: alph_debug_tmp(:,:,:,:), alph_debug(:,:,:,:)
 real(dp), pointer :: beta_debug_tmp(:,:,:,:), beta_debug(:,:,:,:)
 real(dp), pointer :: omeg_debug_tmp(:,:,:,:), omeg_debug(:,:,:,:)
#ifdef MPI
 integer :: MPIerror
#endif
! Lattice constants and spacings
 do idim = 1,3
    lbox(idim) = dsqrt(dot_product(cell(:,idim),cell(:,idim)))
    dbox(idim) = lbox(idim)/dble(ngrid(idim))
 enddo
! Grids
 do ix = 1,Nx
    xbox(ix) = dbox(1)*dble(ix-1)
 enddo
 do iy = 1,Ny
    ybox(iy) = dbox(2)*dble(iy-1)
 enddo
 do iz = 1,Nz
    zbox(iz) = dbox(3)*dble(iz-1)
 enddo
! Vacuum as background
 do iiz = 1,Nzl
    do iiy = 1,Nyl
       do iix = 1,Nxl
           eps(iix,iiy,iiz)   =  eps_vac
          alph(iix,iiy,iiz,:) = alph_vac
          beta(iix,iiy,iiz,:) = beta_vac
          omeg(iix,iiy,iiz,:) = omeg_vac
       enddo
    enddo
 enddo
!
 dx = dbox(1) ; dy = dbox(2) ; dz = dbox(3)
! Objects and grids
 do iobj = 1,Nobj
     epsl = nfobj(iobj)%eps
    alphl = nfobj(iobj)%alph
    betal = nfobj(iobj)%beta
    omegl = nfobj(iobj)%omeg
    do iiz = 1, Nzl
       do iiy = 1, Nyl
          do iix = 1, Nxl
             judge = .false.
             select case(nfobj(iobj)%geometry)
             case('sph')
                 xc = nfobj(iobj)%sph%xc
                 yc = nfobj(iobj)%sph%yc
                 zc = nfobj(iobj)%sph%zc
                 radius = nfobj(iobj)%sph%radius
                 call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
                 judge = isinsphere(ix, iy, iz, dx, dy, dz, xc, yc, zc, radius)
             case('shl')
                 xc = nfobj(iobj)%shl%xc
                 yc = nfobj(iobj)%shl%yc
                 zc = nfobj(iobj)%shl%zc
                 iradius = nfobj(iobj)%shl%iradius
                 oradius = nfobj(iobj)%shl%oradius
                 call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
                 judge = isinshell(ix, iy, iz, dx, dy, dz, xc, yc, zc, iradius, oradius)
             case('wir')
                 xc = nfobj(iobj)%wir%xc
                 yc = nfobj(iobj)%wir%yc
                 zc = nfobj(iobj)%wir%zc
                 orient = nfobj(iobj)%wir%orient
                 length = nfobj(iobj)%wir%length
                 radius = nfobj(iobj)%wir%radius
                 call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
                 judge = isinwire(ix, iy, iz, dx, dy, dz, xc, yc, zc, orient, length, radius)
             case('rct')
                 xc = nfobj(iobj)%rct%xc
                 yc = nfobj(iobj)%rct%yc
                 zc = nfobj(iobj)%rct%zc
                 lengthx = nfobj(iobj)%rct%lengthx
                 lengthy = nfobj(iobj)%rct%lengthy
                 lengthz = nfobj(iobj)%rct%lengthz
                 call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
                 judge = isinrectangle(ix, iy, iz, dx, dy, dz, xc, yc, zc, lengthx, lengthy, lengthz)
             endselect
             if(judge) then
                eps(iix,iiy,iiz)   =  epsl
               alph(iix,iiy,iiz,:) = alphl
               beta(iix,iiy,iiz,:) = betal
               omeg(iix,iiy,iiz,:) = omegl
             endif
             if(nfobj(iobj)%geometry.eq.'sph_smr') then
               xc = nfobj(iobj)%sph_smr%xc
               yc = nfobj(iobj)%sph_smr%yc
               zc = nfobj(iobj)%sph_smr%zc
               radius = nfobj(iobj)%sph_smr%radius
               smear  = nfobj(iobj)%sph_smr%smear
               call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
               call config_sph_smr(ix, iy, iz, dx, dy, dz, xc, yc, zc, radius, smear, &
                                   betal, betas, affected)
               beta(iix,iiy,iiz,:) = beta(iix,iiy,iiz,:) + betas
               if(affected) then
                  eps(iix,iiy,iiz)   =  epsl
                 alph(iix,iiy,iiz,:) = alphl
                 omeg(iix,iiy,iiz,:) = omegl
               endif
             endif
             if(nfobj(iobj)%geometry.eq.'rct_smr') then
               xc = nfobj(iobj)%rct_smr%xc
               yc = nfobj(iobj)%rct_smr%yc
               zc = nfobj(iobj)%rct_smr%zc
               lengthx = nfobj(iobj)%rct_smr%lengthx
               lengthy = nfobj(iobj)%rct_smr%lengthy
               lengthz = nfobj(iobj)%rct_smr%lengthz
               smrdir  = nfobj(iobj)%rct_smr%smrdir
               smear   = nfobj(iobj)%rct_smr%smear
               call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
               call config_rct_smr(ix, iy, iz, dx, dy, dz, xc, yc, zc, &
                                   lengthx, lengthy, lengthz, smrdir, smear, &
                                   betal, betas, affected)
               beta(iix,iiy,iiz,:) = beta(iix,iiy,iiz,:) + betas
               if(affected) then
                  eps(iix,iiy,iiz)   =  epsl
                 alph(iix,iiy,iiz,:) = alphl
                 omeg(iix,iiy,iiz,:) = omegl
               endif
             endif
          enddo
       enddo
    enddo
 enddo
! For debug, symmetrize
 call re_alloc(alph_debug_tmp,1,Nx,1,Ny,1,Nz,1,Nosc,name='beta_debug_tmp',routine='ConfigNFSpace')
 call re_alloc(alph_debug,1,Nx,1,Ny,1,Nz,1,Nosc,name='beta_debug',routine='ConfigNFSpace')
 call re_alloc(beta_debug_tmp,1,Nx,1,Ny,1,Nz,1,Nosc,name='beta_debug_tmp',routine='ConfigNFSpace')
 call re_alloc(beta_debug,1,Nx,1,Ny,1,Nz,1,Nosc,name='beta_debug',routine='ConfigNFSpace')
 call re_alloc(omeg_debug_tmp,1,Nx,1,Ny,1,Nz,1,Nosc,name='omeg_debug_tmp',routine='ConfigNFSpace')
 call re_alloc(omeg_debug,1,Nx,1,Ny,1,Nz,1,Nosc,name='omeg_debug',routine='ConfigNFSpace')
 alph_debug_tmp = 0.0_dp
 beta_debug_tmp = 0.0_dp
 omeg_debug_tmp = 0.0_dp
!
 do iiz = 1,Nzl
    do iiy = 1,Nyl
       do iix = 1,Nxl
          call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
          alph_debug_tmp(ix,iy,iz,:) = alph(iix,iiy,iiz,:)
          beta_debug_tmp(ix,iy,iz,:) = beta(iix,iiy,iiz,:)
          omeg_debug_tmp(ix,iy,iz,:) = omeg(iix,iiy,iiz,:)
       enddo
    enddo
 enddo
#ifdef MPI
 call MPI_AllReduce(alph_debug_tmp(1,1,1,1),alph_debug(1,1,1,1),Nx*Ny*Nz, &
      MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
 call MPI_AllReduce(beta_debug_tmp(1,1,1,1),beta_debug(1,1,1,1),Nx*Ny*Nz, &
      MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
 call MPI_AllReduce(omeg_debug_tmp(1,1,1,1),omeg_debug(1,1,1,1),Nx*Ny*Nz, &
      MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
#else
 alph_debug = alph_debug_tmp
 beta_debug = beta_debug_tmp
 omeg_debug = omeg_debug_tmp
#endif
! Symmetrize along x
 do iz = 1,Nz
    do iy = 1,Ny/2+1
       do ix = 2,Nx/2+1
          alph_debug(Nx-ix+2,iy,iz,:) = alph_debug(ix,iy,iz,:)
          beta_debug(Nx-ix+2,iy,iz,:) = beta_debug(ix,iy,iz,:)
          omeg_debug(Nx-ix+2,iy,iz,:) = omeg_debug(ix,iy,iz,:)
       enddo
    enddo
 enddo
! Symmetriz along y
 do iz = 1,Nz
    do ix = 1,Nx
       do iy = 2,Ny/2+1
          alph_debug(ix,Ny-iy+2,iz,:) = alph_debug(ix,iy,iz,:)
          beta_debug(ix,Ny-iy+2,iz,:) = beta_debug(ix,iy,iz,:)
          omeg_debug(ix,Ny-iy+2,iz,:) = omeg_debug(ix,iy,iz,:)
       enddo
    enddo
 enddo
! Broadcast
 do iiz = 1,Nzl
    do iiy = 1,Nyl
       do iix = 1,Nxl
          call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
          alph(iix,iiy,iiz,:) = alph_debug(ix,iy,iz,:)
          beta(iix,iiy,iiz,:) = beta_debug(ix,iy,iz,:)
          omeg(iix,iiy,iiz,:) = omeg_debug(ix,iy,iz,:)
       enddo
    enddo
 enddo
!
 open(unit=70,file='beta_debug_z.ascii')
 do iz = 1,Nz
    write(70,*) dz*dble(iz-1), alph_debug(Nx/2+1,Ny/2+1,iz,:), beta_debug(Nx/2+1,Ny/2+1,iz,:), omeg_debug(Nx/2+1,Ny/2+1,iz,:)
 enddo
 close(70)
 open(unit=71,file='beta_debug_x.ascii')
 do ix = 1,Nx
    write(71,*) dx*dble(ix-1), alph_debug(ix,Ny/2+1,Nz/2+1,:), beta_debug(ix,Ny/2+1,Nz/2+1,:), omeg_debug(ix,Ny/2+1,Nz/2+1,:)
 enddo
 close(71)
 open(unit=72,file='beta_debug_y.ascii')
 do iy = 1,Ny
    write(72,*) dy*dble(iy-1), alph_debug(Nx/2+1,iy,Nz/2+1,:), beta_debug(Nx/2+1,iy,Nz/2+1,:), omeg_debug(Nx/2+1,iy,Nz/2+1,:)
 enddo
 close(72)
!
 call de_alloc(alph_debug_tmp, name='alph_debug_tmp')
 call de_alloc(alph_debug,     name='alph_debug')
 call de_alloc(beta_debug_tmp, name='beta_debug_tmp')
 call de_alloc(beta_debug,     name='beta_debug')
 call de_alloc(omeg_debug_tmp, name='omeg_debug_tmp')
 call de_alloc(omeg_debug,     name='omeg_debug')
!
 end subroutine ConfigSpace
!
 subroutine lin2mesh(ng, n1, n2, n3, f1d, f3d)
 use parallel, only : Node, Nodes
 implicit none
 integer,      intent(in)  :: ng, n1, n2, n3
 real(grid_p), intent(in)  :: f1d(ng)
 real(dp),     intent(out) :: f3d(n1,n2,n3)
! Local variables
 integer :: ix, iy, iz, ir
!
 if(n1*n2*n3.ne.ng) &
  call die('ERROR: ng != n1*n2*n3')
!
 ir = 1
 do iz = 1,n3
    do iy = 1,n2
       do ix = 1,n1
#ifdef GRID_DP
          f3d(ix,iy,iz) = f1d(ir)
#else
          f3d(ix,iy,iz) = dble(f1d(ir))
#endif
          ir = ir + 1
       enddo
    enddo
 enddo
!
 end subroutine lin2mesh
!
 subroutine mesh2lin(ng, n1, n2, n3, f3d, f1d)
 implicit none
 integer,      intent(in)  :: ng, n1, n2, n3
 real(dp),     intent(in)  :: f3d(n1,n2,n3)
 real(grid_p), intent(out) :: f1d(ng)
! Local variables
 integer :: ix, iy, iz, ir
!
 if(n1*n2*n3.ne.ng) &
  call die('ERROR: ng != n1*n2*n3')
!
 ir = 1
 do iz = 1,n3
    do iy = 1,n2
       do ix = 1,n1
#ifdef GRID_DP
          f1d(ir) = f3d(ix,iy,iz)
#else
          f1d(ir) = real(f3d(ix,iy,iz))
#endif
          ir = ir + 1
       enddo
    enddo
 enddo
!
 end subroutine mesh2lin
!
 subroutine efd2pol(n1, n2, n3, efd_i, pol_o)
 use m_variables, only : Nosc, eps0, beta, omeg, pol
 implicit none
 integer,  intent(in)  :: n1, n2, n3
 real(dp), intent(in)  :: efd_i(n1,n2,n3,3)
 real(dp), intent(out) :: pol_o(n1,n2,n3,3)
! Local variable
 integer :: iosc, idim, ix, iy, iz
!
 pol_o = 0.0_dp
 do iosc = 1,Nosc
    do idim = 1,3
       pol(:,:,:,idim,iosc) = eps0*beta(:,:,:,iosc)*efd_i(:,:,:,idim)
       do iz = 1,n3
          do iy = 1,n2
             do ix = 1,n1
                if(abs(omeg(ix,iy,iz,iosc)).gt.abstolnfr) &
                 pol(ix,iy,iz,idim,iosc) = pol(ix,iy,iz,idim,iosc)/omeg(ix,iy,iz,iosc)**2
             enddo
          enddo
       enddo
    enddo
    pol_o = pol_o + pol(:,:,:,:,iosc)
 enddo
!
 end subroutine efd2pol
!
 subroutine pol2rho(n1, n2, n3, pol_i, rho_o)
 use m_variables, only : nacc, dbox
 use m_numrecipe
 implicit none
 integer,  intent(in)  :: n1, n2, n3
 real(dp), intent(in)  :: pol_i(n1,n2,n3,3)
 real(dp), intent(out) :: rho_o(n1,n2,n3)
! Local variables
 real(dp) :: dx, dy, dz
 real(dp) :: auxx(nacc,n2,n3), auxy(n1,nacc,n3), auxz(n1,n2,nacc)
!
!if(Node.eq.0) write(6,'(a)') 'pol2rho: Calculating induced density...'
!
 dx = dbox(1) ; dy = dbox(2) ; dz = dbox(3)
!
 if(nacc.eq.4) then
   call auxv3d_5pt(n1, n2, n3, pol_i, auxx, auxy, auxz)
   call divg3d_5pt(dx, dy, dz, n1, n2, n3, auxx, auxy, auxz, pol_i, rho_o)
 else
   call auxv3d_3pt(n1, n2, n3, pol_i, auxx, auxy, auxz)
   call divg3d_3pt(dx, dy, dz, n1, n2, n3, auxx, auxy, auxz, pol_i, rho_o)
 end if
 rho_o = -rho_o
!
 end subroutine pol2rho
!
 subroutine phi2efd(n1, n2, n3, phi_i, efd_o)
 use m_variables, only : nacc, dbox
 use m_numrecipe
 implicit none
 integer,  intent(in)  :: n1, n2, n3
 real(dp), intent(in)  :: phi_i(n1,n2,n3)
 real(dp), intent(out) :: efd_o(n1,n2,n3,3)
! Local variables
 real(dp) :: dx, dy, dz
 real(dp) :: auxx(nacc,n2,n3), auxy(n1,nacc,n3), auxz(n1,n2,nacc)
!
!if(Node.eq.0) write(6,'(a)') 'phi2efd: Calculating induced electric field...'
!
 dx = dbox(1) ; dy = dbox(2) ; dz = dbox(3)
!
 if(nacc.eq.4) then
   call auxp3d_5pt(n1, n2, n3, phi_i, auxx, auxy, auxz)
   call grad3d_5pt(dx, dy, dz, n1, n2, n3, auxx, auxy, auxz, phi_i, efd_o)
 else
   call auxp3d_3pt(n1, n2, n3, phi_i, auxx, auxy, auxz)
   call grad3d_3pt(dx, dy, dz, n1, n2, n3, auxx, auxy, auxz, phi_i, efd_o)
 end if
 efd_o = -efd_o
!
 end subroutine phi2efd
!
 subroutine efd2pol3d(efd_i, pol_o)
 use m_variables, only : Nosc, Nxl, Nyl, Nzl, eps0, beta, omeg, pol
 implicit none
 real(dp), intent(in)  :: efd_i(Nxl,Nyl,Nzl,3)
 real(dp), intent(out) :: pol_o(Nxl,Nyl,Nzl,3)
! Local variables
 integer :: n1l, n2l, n3l
 integer :: iosc, idim, ix, iy, iz
!
 n1l = Nxl ; n2l = Nyl ; n3l = Nzl
!
 pol_o = 0.0_dp
 do iosc = 1,Nosc
    do idim = 1,3
       pol(:,:,:,idim,iosc) = eps0*beta(:,:,:,iosc)*efd_i(:,:,:,idim)
       do iz = 1,n3l
          do iy = 1,n2l
             do ix = 1,n1l
                if(abs(omeg(ix,iy,iz,iosc)).gt.abstolnfr) &
                 pol(ix,iy,iz,idim,iosc) = pol(ix,iy,iz,idim,iosc)/omeg(ix,iy,iz,iosc)**2
             enddo
          enddo
       enddo
    enddo
    pol_o = pol_o + pol(:,:,:,:,iosc)
 enddo
!
 end subroutine efd2pol3d
!
 subroutine pol2rho3d(pol_i, rho_o)
 use parallel,    only : Node, Nodes
 use mesh,        only : nsm
 use m_numeric3d, only : divergence
 use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
 implicit none
 real(dp), intent(in)  :: pol_i(Nxl,Nyl,Nzl,3)
 real(dp), intent(out) :: rho_o(Nxl,Nyl,Nzl)
! Local variables
 integer :: nmesh(3), nmeshl(3)
! For debug
!integer :: n1, n2, n3, n1l, n2l, n3l
!integer :: ix, iy, iz, iix, iiy, iiz, idim
!real(dp), pointer :: rho_debug_tmp(:,:,:), rho_debug(:,:,:)
!real(dp), pointer :: pol_serial_tmp(:,:,:,:), pol_serial(:,:,:,:), rho_serial(:,:,:)
!#ifdef MPI
!integer :: MPIerror
!#endif
!
 nmesh(1) = Nx  ; nmesh(2) = Ny  ; nmesh(3) = Nz
 nmeshl(1)= Nxl ; nmeshl(2)= Nyl ; nmeshl(3)= Nzl
!
 call divergence(pol_i, nmesh, nmeshl, nsm, rho_o)
 rho_o = -rho_o
! For debug
!n1 = nmesh(1)  ; n2 = nmesh(2)  ; n3 = nmesh(3)
!n1l= nmeshl(1) ; n2l= nmeshl(2) ; n3l= nmeshl(3)
!
!call re_alloc(rho_debug_tmp,1,n1,1,n2,1,n3,name='rho_debug_tmp',routine='pol2rho3d')
!call re_alloc(rho_debug,1,n1,1,n2,1,n3,name='rho_debug',routine='pol2rho3d')
!call re_alloc(pol_serial_tmp,1,n1,1,n2,1,n3,1,3,name='pol_serial_tmp',routine='pol2rho3d')
!call re_alloc(pol_serial,1,n1,1,n2,1,n3,1,3,name='pol_serial',routine='pol2rho3d')
!call re_alloc(rho_serial,1,n1,1,n2,1,n3,name='rho_serial',routine='pol2rho3d')
!
!rho_debug_tmp = 0.0_dp
!pol_serial_tmp = 0.0_dp
!
!do iiz = 1,n3l
!   do iiy = 1,n2l
!      do iix = 1,n1l
!         call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
!         rho_debug_tmp(ix,iy,iz) = rho_o(iix,iiy,iiz)
!         pol_serial_tmp(ix,iy,iz,:) = pol_i(iix,iiy,iiz,:)
!      enddo
!   enddo
!enddo
!
!#ifdef MPI
!call MPI_AllReduce(rho_debug_tmp(1,1,1),rho_debug(1,1,1),n1*n2*n3, &
!     MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
!call MPI_AllReduce(pol_serial_tmp(1,1,1,1),pol_serial(1,1,1,1),n1*n2*n3*3, &
!     MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
!#else
!rho_debug = rho_debug_tmp
!pol_serial = pol_serial_tmp
!#endif
!
!call pol2rho(n1, n2, n3, pol_serial, rho_serial)
!
!if(Node.eq.0) &
! write(6,*) 'debug: divergence', sum(abs(rho_debug)), sum(abs(rho_serial)), sum(abs(rho_debug))-sum(abs(rho_serial))
!
 end subroutine pol2rho3d
!
!subroutine rho2phi3d(cell, eps_i, rho_i, phi_o)
!use parallel,    only : Node
!use mesh,        only : nsm
!use m_numeric3d, only : conjgrad_mix, conjgrad_mix_precond
!use m_variables, only : eps0_1
!use m_variables, only : cgmaxiter, cgtol, cginfo
!use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
!implicit none
!real(dp), intent(in)    :: cell(3,3)
!real(dp), intent(in)    :: eps_i(Nxl,Nyl,Nzl), rho_i(Nxl,Nyl,Nzl)
!real(dp), intent(inout) :: phi_o(Nxl,Nyl,Nzl)
! Local variables
!integer :: nmesh(3), nmeshl(3)
!
!nmesh(1) = Nx  ; nmesh(2) = Ny  ; nmesh(3) = Nz
!nmeshl(1)= Nxl ; nmeshl(2)= Nyl ; nmeshl(3)= Nzl
!
!if(Node.eq.0) write(6,'(a)') ' Calculating induced potential...'
!
!call conjgrad_mix(cgmaxiter, cgtol, eps_i, eps0_1*rho_i, nmesh, nmeshl, nsm, phi_o, cginfo)
!call conjgrad_mix_precond(cell, cgmaxiter, cgtol, eps_i, eps0_1*rho_i, nmesh, nmeshl, nsm, phi_o, cginfo)
!
!return
!
!end subroutine rho2phi3d
!
 subroutine phi2efd3d(phi_i, efd_o)
 use parallel,    only : Node, Nodes
 use mesh,        only : nsm
 use m_numeric3d, only : gradient
 use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
 implicit none
 real(dp), intent(in)  :: phi_i(Nxl,Nyl,Nzl)
 real(dp), intent(out) :: efd_o(Nxl,Nyl,Nzl,3)
! Local variables
 integer :: nmesh(3), nmeshl(3)
! For debug
!integer :: n1, n2, n3, n1l, n2l, n3l
!integer :: ix, iy, iz, iix, iiy, iiz, idim
!real(dp), pointer :: efd_debug_tmp(:,:,:,:), efd_debug(:,:,:,:)
!real(dp), pointer :: phi_serial_tmp(:,:,:), phi_serial(:,:,:), efd_serial(:,:,:,:)
!#ifdef MPI
!integer :: MPIerror
!#endif
!
 nmesh(1) = Nx  ; nmesh(2) = Ny  ; nmesh(3) = Nz
 nmeshl(1)= Nxl ; nmeshl(2)= Nyl ; nmeshl(3)= Nzl
!
 call gradient(phi_i, nmesh, nmeshl, nsm, efd_o)
 efd_o = -efd_o
! For debug
!n1 = nmesh(1)  ; n2 = nmesh(2)  ; n3 = nmesh(3)
!n1l= nmeshl(1) ; n2l= nmeshl(2) ; n3l= nmeshl(3)
!
!call re_alloc(efd_debug_tmp,1,n1,1,n2,1,n3,1,3,name='efd_debug_tmp',routine='phi2efd3d')
!call re_alloc(efd_debug,1,n1,1,n2,1,n3,1,3,name='efd_debug',routine='phi2efd3d')
!call re_alloc(phi_serial_tmp,1,n1,1,n2,1,n3,name='phi_serial_tmp',routine='phi2efd3d')
!call re_alloc(phi_serial,1,n1,1,n2,1,n3,name='phi_serial',routine='phi2efd3d')
!call re_alloc(efd_serial,1,n1,1,n2,1,n3,1,3,name='efd_serial',routine='phi2efd3d')
!
!efd_debug_tmp = 0.0_dp
!phi_serial_tmp = 0.0_dp
!
!do iiz = 1,n3l
!   do iiy = 1,n2l
!      do iix = 1,n1l
!         call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
!         efd_debug_tmp(ix,iy,iz,:) = efd_o(iix,iiy,iiz,:)
!         phi_serial_tmp(ix,iy,iz) = phi_i(iix,iiy,iiz)
!      enddo
!   enddo
!enddo
!#ifdef MPI
!call MPI_AllReduce(efd_debug_tmp(1,1,1,1),efd_debug(1,1,1,1),n1*n2*n3*3, &
!     MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
!call MPI_AllReduce(phi_serial_tmp(1,1,1),phi_serial(1,1,1),n1*n2*n3, &
!     MPI_double_precision,MPI_Sum,MPI_Comm_World,MPIerror)
!#else
!efd_debug = efd_debug_tmp
!phi_serial = phi_serial_tmp
!#endif
!
!call phi2efd(n1, n2, n3, phi_serial, efd_serial)
!
!if(Node.eq.0) &
! write(6,*) 'debug: gradient', sum(abs(efd_debug)), sum(abs(efd_serial)), sum(abs(efd_debug))-sum(abs(efd_serial))
!
 end subroutine phi2efd3d
!
 subroutine pol2dip3d(pol_i, dip_o)
 use mesh,        only : nsm
 use m_numeric3d, only : integration
 use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
 implicit none
 real(dp), intent(in)  :: pol_i(Nxl,Nyl,Nzl,3)
 real(dp), intent(out) :: dip_o(3)
! Local variables
 integer :: idim
 integer :: nmesh(3), nmeshl(3)
!
 nmesh(1) = Nx  ; nmesh(2) = Ny  ; nmesh(3) = Nz
 nmeshl(1)= Nxl ; nmeshl(2)= Nyl ; nmeshl(3)= Nzl
!
 do idim = 1,3
    call integration(pol_i(:,:,:,idim), nmesh, nmeshl, nsm, dip_o(idim))
 enddo
!
 end subroutine pol2dip3d
!
 subroutine leap_frog(dt, efd_i, cur_o, pol_o)
 use m_variables, only : Nosc, Nxl, Nyl, Nzl, eps0
 use m_variables, only : alph, beta, omeg, cur, pol
 implicit none
 real(dp), intent(in)  :: dt
 real(dp), intent(in)  :: efd_i(Nxl,Nyl,Nzl,3)
 real(dp), intent(out) :: cur_o(Nxl,Nyl,Nzl,3), pol_o(Nxl,Nyl,Nzl,3)
! Local variables
 integer  :: ix, iy, iz, idim, iosc
 real(dp) :: alphl, betal, omegl
 real(dp) ::  efdl,  poll,  curl
 real(dp) :: numer, denom
!
 do iosc = 1,Nosc
    do idim = 1,3
       do iz = 1,Nzl
          do iy = 1,Nyl
             do ix = 1,Nxl
                alphl = alph(ix,iy,iz,iosc)
                betal = beta(ix,iy,iz,iosc)
                omegl = omeg(ix,iy,iz,iosc)
!
                 efdl = efd_i(ix,iy,iz,idim)
                 poll =   pol(ix,iy,iz,idim,iosc)
                 curl =   cur(ix,iy,iz,idim,iosc)
!
                 numer = 1.d0 - 0.5d0*alphl*dt
                 denom = 1.d0 + 0.5d0*alphl*dt
!
                 cur(ix,iy,iz,idim,iosc) = &
                 numer/denom*curl - dt/denom*(omegl*omegl*poll - eps0*betal*efdl)
!
             enddo
          enddo
       enddo
    enddo
 enddo
!
 pol = pol + dt * cur
!
 cur_o = 0.0_dp
 pol_o = 0.0_dp
 do iosc = 1,Nosc
    cur_o = cur_o + cur(:,:,:,:,iosc)
    pol_o = pol_o + pol(:,:,:,:,iosc)
 enddo
!
 end subroutine leap_frog
!
 function isinsphere(ix, iy, iz, dx, dy, dz, xc, yc, zc, radius)
 implicit none
 integer,  intent(in) :: ix, iy, iz
 real(dp), intent(in) :: dx, dy, dz
 real(dp), intent(in) :: xc, yc, zc
 real(dp), intent(in) :: radius
! Local variables
 logical  :: isinsphere
 real(dp) :: x, y, z, dist2
!
 isinsphere = .false.
!
 x = dble(ix-1)*dx
 y = dble(iy-1)*dy
 z = dble(iz-1)*dz
!
 dist2 = (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)
!
 if(dist2 .le. radius*radius) isinsphere = .true.
!
 end function isinsphere
!
 function isinshell(ix, iy, iz, dx, dy, dz, xc, yc, zc, iradius, oradius)
 implicit none
 integer,  intent(in) :: ix, iy, iz
 real(dp), intent(in) :: dx, dy, dz
 real(dp), intent(in) :: xc, yc, zc
 real(dp), intent(in) :: iradius, oradius
! Local variables
 logical  :: isinshell
 real(dp) :: x, y, z, dist2
!
 isinshell = .false.
!
 x = dble(ix-1)*dx
 y = dble(iy-1)*dy
 z = dble(iz-1)*dz
!
 dist2 = (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)
!
 if(dist2 .ge. iradius*iradius .and. dist2 .le. oradius*oradius) &
  isinshell = .true.
!
 end function isinshell
!
 function isinwire(ix, iy, iz, dx, dy, dz, xc, yc, zc, orient, length, radius)
 implicit none
 integer,  intent(in) :: ix, iy, iz
 integer,  intent(in) :: orient
 real(dp), intent(in) :: dx, dy, dz
 real(dp), intent(in) :: xc, yc, zc
 real(dp), intent(in) :: length, radius
! Local variables
 logical  :: isinwire
 logical  :: isincirc, isinleng
 real(dp) :: x, y, z, dist2
!
 isinwire = .false.
 isinleng = .false.
!
 dist2 = 0.0_dp
!
 x = dble(ix-1)*dx
 y = dble(iy-1)*dy
 z = dble(iz-1)*dz
!
 select case(orient)
!
 case(1)
   dist2 = (y-yc)*(y-yc)+(z-zc)*(z-zc)
   isinleng = (dabs(x-xc).le.(0.5d0*length))
 case(2)
   dist2 = (x-xc)*(x-xc)+(z-zc)*(z-zc)
   isinleng = (dabs(y-yc).le.(0.5d0*length))
 case(3)
   dist2 = (x-xc)*(x-xc)+(y-yc)*(y-yc)
   isinleng = (dabs(z-zc).le.(0.5d0*length))
!
 end select
!
 isincirc = (dist2.le.radius*radius)
!
 isinwire = (isincirc.and.isinleng)
!
 end function isinwire
!
 function isinrectangle(ix, iy, iz, dx, dy, dz, xc, yc, zc, lengthx, lengthy, lengthz)
 implicit none
 integer,  intent(in) :: ix, iy, iz
 real(dp), intent(in) :: dx, dy, dz
 real(dp), intent(in) :: xc, yc, zc
 real(dp), intent(in) :: lengthx, lengthy, lengthz
! Local variables
 logical  :: isinrectangle
 logical  :: isinlengthx, isinlengthy, isinlengthz
 real(dp) :: x, y, z
!
 isinrectangle = .false.
!
 x = dble(ix-1)*dx
 y = dble(iy-1)*dy
 z = dble(iz-1)*dz
!
 isinlengthx = (dabs(x-xc).le.(0.5d0*lengthx))
 isinlengthy = (dabs(y-yc).le.(0.5d0*lengthy))
 isinlengthz = (dabs(z-zc).le.(0.5d0*lengthz))
!
 isinrectangle = (isinlengthx.and.isinlengthy.and.isinlengthz)
!
 end function isinrectangle
!
 subroutine config_sph_smr(ix, iy, iz, dx, dy, dz, xc, yc, zc, &
                           radius, smear, betal, betas, affected)
 use m_variables, only : Nosc, eps_vac
 implicit none
 integer,  intent(in)  :: ix, iy, iz
 real(dp), intent(in)  :: dx, dy, dz, xc, yc, zc, radius, smear
 real(dp), intent(in)  :: betal(Nosc)
 real(dp), intent(out) :: betas(Nosc)
 logical,  intent(out) :: affected
! Local variables 
 real(dp) :: x, y, z, dist2
!
 affected = .false.
!
 x = dble(ix-1)*dx
 y = dble(iy-1)*dy
 z = dble(iz-1)*dz
!
 dist2 = (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)
!
 betas = betal / (1.0_dp + dexp((sqrt(dist2)-radius)/smear))
!
 if(maxval(betas) .gt. abstolnfr) affected = .true.
!
 end subroutine config_sph_smr
!
 subroutine config_rct_smr(ix, iy, iz, dx, dy, dz, xc, yc, zc, &
                           lengthx, lengthy, lengthz, smrdir, smear, &
                           betal, betas, affected)
 use m_variables, only : Nosc, eps_vac
 implicit none
 integer,  intent(in)  :: ix, iy, iz, smrdir
 real(dp), intent(in)  :: dx, dy, dz, xc, yc, zc
 real(dp), intent(in)  :: lengthx, lengthy, lengthz, smear
 real(dp), intent(in)  :: betal(Nosc)
 real(dp), intent(out) :: betas(Nosc)
 logical,  intent(out) :: affected
! Local variables
 real(dp) :: x, y, z
 logical  :: isinlengthx, isinlengthy, isinlengthz
!
 affected = .false.
!
 x = dble(ix-1)*dx
 y = dble(iy-1)*dy
 z = dble(iz-1)*dz
!
 isinlengthx = (dabs(x-xc).le.(0.5d0*lengthx))
 isinlengthy = (dabs(y-yc).le.(0.5d0*lengthy))
 isinlengthz = (dabs(z-zc).le.(0.5d0*lengthz))
!
 select case(smrdir)
 case(1)
     if(isinlengthy .and. isinlengthz) then
       betas = betal *(1.0_dp / (1.0_dp + dexp(( x-(xc+0.5_dp*lengthx))/smear)) &
                     + 1.0_dp / (1.0_dp + dexp((-x+(xc-0.5_dp*lengthx))/smear)) &
                     - 1.0_dp)
       if(maxval(betas) .gt. abstolnfr) affected = .true.
     else
       betas = 0.0_dp
     endif
 case(2)
     if(isinlengthx .and. isinlengthz) then
       betas = betal *(1.0_dp / (1.0_dp + dexp(( y-(yc+0.5_dp*lengthy))/smear)) &
                     + 1.0_dp / (1.0_dp + dexp((-y+(yc-0.5_dp*lengthy))/smear)) &
                     - 1.0_dp)
       if(maxval(betas) .gt. abstolnfr) affected = .true.
     else
       betas = 0.0_dp
     endif
 case(3)
     if(isinlengthx .and. isinlengthy) then
       betas = betal *(1.0_dp / (1.0_dp + dexp(( z-(zc+0.5_dp*lengthz))/smear)) &
                     + 1.0_dp / (1.0_dp + dexp((-z+(zc-0.5_dp*lengthz))/smear)) &
                     - 1.0_dp)
       if(maxval(betas) .gt. abstolnfr) affected = .true.
     else
       betas = 0.0_dp
     endif
 case default
     call die('ERROR: rct_smr has and invalid smrdir value!')
 end select
!
 end subroutine config_rct_smr
!
 subroutine outnfr(string)
 use atm_types,   only : species
 use files,       only : slabel
 use parallel,    only : Node, Nodes
 use siesta_geom, only : na_u, isa, xa, ucell
 use m_variables, only : Nx, Ny, Nz, Nxl, Nyl, Nzl
 use m_variables, only : xbox, ybox, zbox
 use m_variables, only : phi, rhonfr
 implicit none
 character*3, intent(in) :: string
! Local variables
 integer :: nmesh(3), nmeshl(3)
 integer :: ix, iy, iz, iix, iiy, iiz, idim, ia
 character*60 :: fname
 real(dp), dimension(:,:,:), pointer :: phi_ful, rhonfr_ful
 real(dp), dimension(:,:,:), pointer :: phi_red, rhonfr_red
#ifdef MPI
 integer :: MPIerror
#endif
!
  nmesh(1) = Nx  ;  nmesh(2) = Ny  ;  nmesh(3) = Nz
 nmeshl(1) = Nxl ; nmeshl(2) = Nyl ; nmeshl(3) = Nzl
!
 call re_alloc(   phi_ful,1,nmesh(1),1,nmesh(2),1,nmesh(3),name=   'phi_ful',routine='outnfr')
 call re_alloc(rhonfr_ful,1,nmesh(1),1,nmesh(2),1,nmesh(3),name='rhonfr_ful',routine='outnfr')
!
 call re_alloc(   phi_red,1,nmesh(1),1,nmesh(2),1,nmesh(3),name=   'phi_red',routine='outnfr')
 call re_alloc(rhonfr_red,1,nmesh(1),1,nmesh(2),1,nmesh(3),name='rhonfr_red',routine='outnfr')
!
    phi_ful = 0.0_dp
 rhonfr_ful = 0.0_dp
!
 do iiz = 1,nmeshl(3)
    do iiy = 1,nmeshl(2)
       do iix = 1,nmeshl(1)
          call LocalToGlobalGrid(iix, iiy, iiz, Node, Nodes, ix, iy, iz)
             phi_ful(ix,iy,iz) =    phi(iix,iiy,iiz)
          rhonfr_ful(ix,iy,iz) = rhonfr(iix,iiy,iiz)
       enddo
    enddo
 enddo
#ifdef MPI
 call MPI_AllReduce(   phi_ful(1,1,1),   phi_red(1,1,1),nmesh(1)*nmesh(2)*nmesh(3), &
      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
 call MPI_AllReduce(rhonfr_ful(1,1,1),rhonfr_red(1,1,1),nmesh(1)*nmesh(2)*nmesh(3), &
      MPI_double_precision,MPI_sum,MPI_Comm_World,MPIerror)
#else
    phi_red =    phi_ful
 rhonfr_red = rhonfr_ful
#endif
!
 fname = trim(slabel) // '.' // 'RHO_NFR_' // trim(string) // '.cube'
 if(Node.eq.0) then
   open(unit=80,file=fname)
   write(80,'(a)') 'HYBSIESTA CUBE FILE'
   write(80,'(a)') 'STATIC NEAR-FIELD DENSITY'
   write(80,'(i4,3f12.6)') na_u, 0.0_dp, 0.0_dp, 0.0_dp
   do idim = 1,3
      write(80,'(i4,3f12.6)') nmesh(idim), ucell(:,idim)/dble(nmesh(idim))
   enddo
   do ia = 1,na_u
      write(80,'(i4,4f12.6)') species(isa(ia))%z, 0.0_dp, xa(:,ia)
   enddo
   do ix = 1,nmesh(1)
      do iy = 1,nmesh(2)
         do iz = 1,nmesh(3)
            write(80,'(e15.7)',advance='no') rhonfr_red(ix,iy,iz)
            if(mod(iz,6).eq.5) write(80,*)
         enddo
         write(80,*)
      enddo
   enddo
   close(80)
 endif
!
 fname = trim(slabel) // '.' // 'PHI_HYB_' // trim(string) // '.cube'
 if(Node.eq.0) then
   open(unit=81,file=fname)
   write(81,'(a)') 'HYBSIESTA CUBE FILE'
   write(81,'(a)') 'STATIC HYBRID POTENTIAL'
   write(81,'(i4,3f12.6)') na_u, 0.0_dp, 0.0_dp, 0.0_dp
   do idim = 1,3
      write(81,'(i4,3f12.6)') nmesh(idim), ucell(:,idim)/dble(nmesh(idim))
   enddo
   do ia = 1,na_u
      write(81,'(i4,4f12.6)') species(isa(ia))%z, 0.0_dp, xa(:,ia)
   enddo
   do ix = 1,nmesh(1)
      do iy = 1,nmesh(2)
         do iz = 1,nmesh(3)
            write(81,'(e15.7)',advance='no') rhonfr_red(ix,iy,iz)
            if(mod(iz,6).eq.5) write(81,*)
         enddo
         write(81,*)
      enddo
   enddo
   close(81)
! For debug
!  open(unit=99,file='debug_nfr.ascii')
!  do iz = 1,nmesh(3)
!     write(99,*) zbox(iz), rhonfr_red(nmesh(1)/2+1,nmesh(2)/2+1,iz)
!  enddo
!  close(99)
 endif
!
 call de_alloc(   phi_ful,name=   'phi_ful')
 call de_alloc(rhonfr_ful,name='rhonfr_ful')
 call de_alloc(   phi_red,name=   'phi_red')
 call de_alloc(rhonfr_red,name='rhonfr_red')
!
 end subroutine outnfr
!
 end module nearfieldsubs
