! 
! This file is part of the SIESTA package.
!
! Copyright (c) UCLA:
! Yi Gao and Daniel Neuhauser, 2012- .
! 
! Use of this software constitutes agreement with the full conditions
! given in the SIESTA license, as signed by all legitimate users.
!
      module m_numeric3d
      use m_variables, only : nacc, dbox
      use m_numrecipe
      implicit none

      contains

      subroutine integration(f, Mesh, Meshl, nsm, v)
C
C  Parallel integration subroutine
C
C  On input :
C
C  real*8 f()      : contains data to be calculated
C  integer Mesh(3) : contains global dimensions of grid
C  integer Meshl(3): contains local dimensions of grid
C  integer nsm     : number of sub-mesh points - needed as mesh
C                  : points are divided over nodes based on
C                  : coarse mesh rather than fine mesh
C  On exit :
C
C  real*8 v        : contains integration of f
C
C  Yi Gao, July 2012
C
C
C
      use precision,   only :dp
      use parallel,    only : Node, Nodes, ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc
      use m_variables, only : GridTableY, GridTableZ
#ifdef MPI
      use mpi_siesta
#endif
      implicit none

C
C  Passed arguments
C
      integer
     .  Mesh(3), Meshl(3), nsm
      real(dp)
     .  f(Meshl(1),Meshl(2),Meshl(3)), v
C
C  Local variables
C

      integer
     .  iix, iiy, iiz, iz
      integer
     .  n1, n2, n3, n2l, n3l, n, ng, ProcessorZ
      integer
     .  jmin, jmax, kmin, kmax
      real(dp)
     . dx, dy, dz
      real(dp), dimension(:),   pointer :: func1d
      real(dp), dimension(:),   pointer :: valu1d
      real(dp), dimension(:,:), pointer :: valu2d
#ifdef MPI
      integer :: MPIerror
      integer, save ::
     .  Py, Pz, YCommunicator, ZCommunicator,
     .  YNode, YNodes, ZNode, ZNodes
      logical, save :: firsttime = .true.
      real(dp), dimension(:),   pointer :: valu1d_tmp, valu1d_red
      real(dp), dimension(:,:), pointer :: valu2d_tmp, valu2d_red
#endif

      n1 = Mesh(1)
      n2 = Mesh(2)
      n3 = Mesh(3)
 
      n2l= Meshl(2)
      n3l= Meshl(3)

      n = n1*n2l*n3l
      ng= n1*n2*n3

      dx = dbox(1) ; dy = dbox(2) ; dz = dbox(3)

      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')

      if(associated(valu2d)) call de_alloc(valu2d, name='valu2d')
      call re_alloc(valu2d, 1, n2l, 1, n3l, name='valu2d',
     .              routine='integration')
      if(associated(valu1d)) call de_alloc(valu1d, name='valu1d')
      call re_alloc(valu1d, 1, n3l,         name='valu1d',
     .              routine='integration')

      if(associated(func1d)) call de_alloc(func1d, name='func1d')
      call re_alloc(func1d, 1, n1, name='func1d', routine='integration')
      do iiz = 1, n3l
         do iiy = 1, n2l
            func1d = f(:,iiy,iiz)
            call intg1d(dx, n1, func1d, valu2d(iiy,iiz))
         enddo
      enddo

#ifdef MPI
      if(firsttime) then
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1

        call MPI_Comm_Split(MPI_Comm_World, Pz, Py, ZCommunicator,
     .       MPIerror)
        call MPI_Comm_Rank(ZCommunicator, ZNode,  MPIerror)
        call MPI_Comm_Size(ZCommunicator, ZNodes, MPIerror)

        call MPI_Comm_Split(MPI_Comm_World, Py, Pz, YCommunicator,
     .       MPIerror)
        call MPI_Comm_Rank(YCommunicator, YNode,  MPIerror)
        call MPI_Comm_Size(YCommunicator, YNodes, MPIerror)

        firsttime = .false.
      endif

      jmin = sum(GridTableY(1:Py-1)) + 1
      jmax = jmin + GridTableY(Py) - 1

      if(associated(valu2d_tmp))
     .call de_alloc(valu2d_tmp, name='valu2d_tmp')
      call re_alloc(valu2d_tmp, 1, n2, 1, n3l, name='valu2d_tmp',
     .              routine='integration')
      valu2d_tmp = 0.0_dp
      valu2d_tmp(jmin:jmax,:) = valu2d(1:n2l,:)

      if(associated(valu2d_red))
     .call de_alloc(valu2d_red, name='valu2d_red')
      call re_alloc(valu2d_red, 1, n2, 1, n3l, name='valu2d_red',
     .              routine='integration')

      call MPI_AllReduce(valu2d_tmp(1,1), valu2d_red(1,1), n2*n3l,
     .     MPI_double_precision, MPI_sum, ZCommunicator, MPIerror)
      call de_alloc(valu2d_tmp, name='valu2d_tmp')

      if(associated(func1d)) call de_alloc(func1d, name='func1d')
      call re_alloc(func1d, 1, n2, name='func1d', routine='integration')

      do iiz = 1, n3l
         func1d = valu2d_red(:,iiz)
         call intg1d(dy, n2, func1d, valu1d(iiz))
      end do
      call de_alloc(valu2d_red, name='valu2d_red')
      call de_alloc(func1d)

      kmin = sum(GridTableZ(1:Pz-1)) + 1
      kmax = kmin + GridTableZ(Pz) - 1

      if(associated(valu1d_tmp))
     .call de_alloc(valu1d_tmp, name='valu1d_tmp')
      call re_alloc(valu1d_tmp, 1, n3, name='valu1d_tmp',
     .              routine='valu1d_tmp')

      valu1d_tmp = 0.0_dp
      valu1d_tmp(kmin:kmax) = valu1d(1:n3l)

      if(associated(valu1d_red))
     .call de_alloc(valu1d_red, name='valu1d_red')
      call re_alloc(valu1d_red, 1, n3, name='valu1d_red',
     .              routine='valu1d_red')

      call MPI_AllReduce(valu1d_tmp(1), valu1d_red(1), n3,
     .     MPI_double_precision, MPI_sum, YCommunicator, MPIerror)
      call de_alloc(valu1d_tmp)

      if(associated(func1d)) call de_alloc(func1d, name='func1d')
      call re_alloc(func1d, 1, n3, name='func1d', routine='integration')

      func1d = valu1d_red
      call intg1d(dz, n3, func1d, v)
      call de_alloc(valu1d_red, name='valu1d_red')

      call de_alloc(func1d, name='func1d')
#else
      if(associated(func1d)) call de_alloc(func1d, name='func1d')
      call re_alloc(func1d, 1, n2, name='func1d', routine='func1d')
      do iz = 1, n3
         func1d = valu2d(:,iz)
         call intg1d(dy, n2, func1d, valu1d(iz))
      enddo
      call de_alloc(func1d, name='func1d')

      if(associated(func1d)) call de_alloc(func1d, name='func1d')
      call re_alloc(func1d, 1, n3, name='func1d')
      func1d = valu1d
      call intg1d(dz, n3, func1d, v)
      call de_alloc(func1d, name='func1d')
#endif
      call de_alloc(valu2d, name='valu2d')
      call de_alloc(valu1d, name='valu1d')

      end subroutine integration

      subroutine divergence(f, Mesh, Meshl, nsm, g)
C
C  Parallel divergence subroutine
C
C  On input :
C
C  real*8 f()      : contains data to be calculated
C  integer Mesh(3) : contains global dimensions of grid
C  integer Meshl(3): contains local dimensions of grid
C  integer nsm     : number of sub-mesh points - needed as mesh
C                  : points are divided over nodes based on
C                  : coarse mesh rather than fine mesh
C  On exit :
C
C  real*8 g()      : contains divergence of f
C
C  Yi Gao, July 2012
C
C
C  Modules
C
      use precision,    only :dp
      use parallel,     only : Node, Nodes, ProcessorY
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C
C  Passed arguments
C
      integer
     .  Mesh(3), Meshl(3), nsm
      real(dp)
     .  f(Meshl(1),Meshl(2),Meshl(3),3), g(Meshl(1),Meshl(2),Meshl(3))
C
C  Local variables
C                                                                                                                                                                                                                  
      integer
     .  n1, n2, n3, n2l, n3l, n, ng, ProcessorZ, i, j, k
      real(dp)
     .  dx, dy, dz
      real(dp), dimension(:), pointer :: aux, func1d, derv1d
#ifdef MPI
      integer
     .  n1lf, nrem, Py, Pz
      real(dp), dimension(:,:,:), pointer :: ft, gt
#endif
      real(dp), dimension(:,:,:), pointer :: gb

C
C  Set mesh size variables
C
      n1 = Mesh(1)
      n2 = Mesh(2)
      n3 = Mesh(3)

      n2l = Meshl(2)
      n3l = Meshl(3)

      n = n1*n2l*n3l
      ng= n1*n2*n3

      dx = dbox(1) ; dy = dbox(2) ; dz = dbox(3)

C
C  Work out processor grid size
C
      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes)
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')

      if(nacc.eq.4) then
        call re_alloc(aux, 1, 4, name='aux', routine='divergence')
      else
        call re_alloc(aux, 1, 2, name='aux', routine='divergence')
      endif
C
C  Differentiation in X direction
C
      call re_alloc(func1d, 1, n1, name='func1d', routine='divergence')
      call re_alloc(derv1d, 1, n1, name='derv1d', routine='divergence')
      do k = 1,n3l
         do j = 1,n2l
            func1d(:) = f(:,j,k,1)
            if(nacc.eq.4) then
              call auxp1d_5pt(n1, func1d, aux)
              call derv1d_5pt(dx, n1, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n1, func1d, aux)
              call derv1d_3pt(dx, n1, aux, func1d, derv1d)
            endif
            g(:,j,k) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
C
      nullify( gb )
      call re_alloc( gb, 1, n1, 1, n2l, 1, n3l, name='gb',
     .               routine='divergence' )
#ifdef MPI
C***********************
C  2-D Processor Grid  *
C***********************
      n1lf = n1/ProcessorY
      nrem = n1 - n1lf*ProcessorY
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1
      if (Py.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )
      call re_alloc( ft, 1, n1lf, 1, n2, 1, n3l, name='ft',
     .               routine='divergence' )
      nullify( gt )
      call re_alloc( gt, 1, n1lf, 1, n2, 1, n3l, name='gt',
     .               routine='divergence' )
C
C  Redistribute data to be distributed by X and Z
C
      call redistribXZ(f(:,:,:,2),n1,n2l,n3l,ft,n1lf,n2,1,nsm,
     .                 Node,Nodes)
C
C  Differentiation in Y direction
C
      call re_alloc(func1d, 1, n2, name='func1d', routine='divergence')
      call re_alloc(derv1d, 1, n2, name='derv1d', routine='divergence')
      do k = 1,n3l
         do i = 1,n1lf
            func1d(:) = ft(i,:,k)
            if(nacc.eq.4) then
              call auxp1d_5pt(n2, func1d, aux)
              call derv1d_5pt(dy, n2, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n2, func1d, aux)
              call derv1d_3pt(dy, n2, aux, func1d, derv1d)
            endif
            gt(i,:,k) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
C
C  Redistribute data back to original form
C
      call redistribXZ(gb,n1,n2l,n3l,gt,n1lf,n2,-1,nsm,Node,Nodes)
      g = g + gb
C
C  Free local memory ready for re-use
C
      call de_alloc( ft, name='ft' )
      call de_alloc( gt, name='gt' )
C
C  Find new distributed x dimension
C
      n1lf = n1/ProcessorZ
      nrem = n1 - n1lf*ProcessorZ
      if (Pz.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )   ! AG
      call re_alloc( ft, 1, n1lf, 1, n2l, 1, n3, name='ft',
     .               routine='divergence' )
      nullify( gt )   ! AG
      call re_alloc( gt, 1, n1lf, 1, n2l, 1, n3, name='gt',
     .               routine='divergence' )
C
C  Redistribute data to be distributed by X and Y
C
      call redistribXY(f(:,:,:,3),n1,n2l,n3l,ft,n1lf,n3,1,nsm,
     .                 Node,Nodes)
C
C  Differentiation in Z direction
C
      call re_alloc(func1d, 1, n3, name='func1d', routine='divergence')
      call re_alloc(derv1d, 1, n3, name='derv1d', routine='divergence')
      do j = 1,n2l
         do i = 1,n1lf
            func1d(:) = ft(i,j,:)
            if(nacc.eq.4) then
              call auxp1d_5pt(n3, func1d, aux)
              call derv1d_5pt(dz, n3, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n3, func1d, aux)
              call derv1d_3pt(dz, n3, aux, func1d, derv1d)
            endif
            gt(i,j,:) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
C
C  Redistribute data to be distributed by Z again
C
      call redistribXY(gb,n1,n2l,n3l,gt,n1lf,n3,-1,nsm,Node,Nodes)
      g = g + gb
C
C  Free local memory
C
      call de_alloc( ft, name='ft' )
      call de_alloc( gt, name='gt' )
#else
C
C  Differentiation in Y direction
C
      call re_alloc(func1d, 1, n2, name='func1d', routine='divergence')
      call re_alloc(derv1d, 1, n2, name='derv1d', routine='divergence')
      do k = 1,n3
         do i = 1,n1
            func1d(:) = f(i,:,k,2)
            if(nacc.eq.4) then
              call auxp1d_5pt(n2, func1d, aux)
              call derv1d_5pt(dy, n2, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n2, func1d, aux)
              call derv1d_3pt(dy, n2, aux, func1d, derv1d)
            endif
            gb(i,:,k) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
      g = g + gb
C
C  Differentiation in Z direction
C
      call re_alloc(func1d, 1, n3, name='func1d', routine='divergence')
      call re_alloc(derv1d, 1, n3, name='derv1d', routine='divergence')
      do j = 1,n2
         do i = 1,n1
            func1d(:) = f(i,j,:,3)
            if(nacc.eq.4) then
              call auxp1d_5pt(n3, func1d, aux)
              call derv1d_5pt(dz, n3, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n3, func1d, aux)
              call derv1d_3pt(dz, n3, aux, func1d, derv1d)
            endif
            gb(i,j,:) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
      g = g + gb
#endif

      call de_alloc(aux, name='aux')
      call de_alloc(gb,  name='gb')

      end subroutine divergence

      subroutine gradient(f, Mesh, Meshl, nsm, g)
C
C  Parallel gradient subroutine
C
C  On input :
C
C  real*8 f()      : contains data to be calculated
C  integer Mesh(3) : contains global dimensions of grid
C  integer Meshl(3): contains local dimensions of grid
C  integer nsm     : number of sub-mesh points - needed as mesh
C                  : points are divided over nodes based on
C                  : coarse mesh rather than fine mesh
C  On exit :
C
C  real*8 g()      : contains gradient of f
C
C  Yi Gao, July 2012
C

C
C  Modules
C
      use precision,    only :dp
      use parallel,     only : Node, Nodes, ProcessorY
      use sys,          only : die
      use alloc,        only : re_alloc, de_alloc
#ifdef MPI
      use mpi_siesta
#endif

      implicit none

C
C  Passed arguments
C
      integer
     .  Mesh(3), Meshl(3), nsm
      real(dp)
     .  f(Meshl(1),Meshl(2),Meshl(3)), g(Meshl(1),Meshl(2),Meshl(3),3)
C
C  Local variables
C
      integer
     .  n1, n2, n3, n2l, n3l, n, ng, ProcessorZ, i, j, k
      real(dp)
     .  dx, dy, dz
      real(dp), dimension(:), pointer :: aux, func1d, derv1d
#ifdef MPI
      integer
     .  n1lf, nrem, Py, Pz
      real(dp), dimension(:,:,:), pointer :: ft, gt
#endif

C
C  Set mesh size variables
C
      n1 = Mesh(1)
      n2 = Mesh(2)
      n3 = Mesh(3)

      n2l = Meshl(2)
      n3l = Meshl(3)

      n = n1*n2l*n3l
      ng = n1*n2*n3

      dx = dbox(1) ; dy = dbox(2) ; dz = dbox(3)

C
C  Work out processor grid size
C
      ProcessorZ = Nodes/ProcessorY
      if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')

      if(nacc.eq.4) then
        call re_alloc( aux, 1, 4, name='aux', routine='gradient')
      else
        call re_alloc( aux, 1, 2, name='aux', routine='gradient')
      endif
C
C  Differentiation in X direction
C
      call re_alloc( func1d, 1, n1, name='func1d', routine='gradient' )
      call re_alloc( derv1d, 1, n1, name='derv1d', routine='gradient' )
      do k = 1,n3l
         do j = 1,n2l
            func1d(:) = f(:,j,k)
            if(nacc.eq.4) then
              call auxp1d_5pt(n1, func1d, aux)
              call derv1d_5pt(dx, n1, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n1, func1d, aux)
              call derv1d_3pt(dx, n1, aux, func1d, derv1d)
            endif
            g(:,j,k,1) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
#ifdef MPI
C***********************
C  2-D Processor Grid  *
C***********************
      n1lf = n1/ProcessorY
      nrem = n1 - n1lf*ProcessorY
      Py = (Node/ProcessorZ) + 1
      Pz = Node - (Py - 1)*ProcessorZ + 1
      if (Py.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )
      call re_alloc( ft, 1, n1lf, 1, n2, 1, n3l, name='ft',
     .               routine='gradient' )
      nullify( gt )
      call re_alloc( gt, 1, n1lf, 1, n2, 1, n3l, name='gt',
     .               routine='gradient' )
C
C  Redistribute data to be distributed by X and Z
C
      call redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,1,nsm,Node,Nodes)
C
C  Differentiation in Y direction
C
      call re_alloc( func1d, 1, n2, name='func1d', routine='gradient' )
      call re_alloc( derv1d, 1, n2, name='derv1d', routine='gradient' )
      do k = 1,n3l
         do i = 1,n1lf
            func1d(:) = ft(i,:,k)
            if(nacc.eq.4) then
              call auxp1d_5pt(n2, func1d, aux)
              call derv1d_5pt(dy, n2, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n2, func1d, aux)
              call derv1d_3pt(dy, n2, aux, func1d, derv1d)
            endif
            gt(i,:,k) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
C
C  Redistribute data back to original form
C
      call redistribXZ(g(:,:,:,2),n1,n2l,n3l,gt,n1lf,n2,-1,nsm,
     .                 Node,Nodes)
C
C  Free local memory ready for re-use
C
      call de_alloc( ft, name='ft' )
      call de_alloc( gt, name='gt' )
C
C  Find new distributed x dimension
C
      n1lf = n1/ProcessorZ
      nrem = n1 - n1lf*ProcessorZ
      if (Pz.le.nrem) n1lf = n1lf + 1
C
C  Allocate local memory
C
      nullify( ft )   ! AG
      call re_alloc( ft, 1, n1lf, 1, n2l, 1, n3, name='ft',
     .               routine='gradient' )
      nullify( gt )   ! AG
      call re_alloc( gt, 1, n1lf, 1, n2l, 1, n3, name='gt',
     .               routine='gradient' )
C
C  Redistribute data to be distributed by X and Y
C
      call redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,1,nsm,Node,Nodes)
C
C  Differentiation in Z direction
C
      call re_alloc( func1d, 1, n3, name='func1d', routine='gradient' )
      call re_alloc( derv1d, 1, n3, name='derv1d', routine='gradient' )
      do j = 1,n2l
         do i = 1,n1lf
            func1d(:) = ft(i,j,:)
            if(nacc.eq.4) then
              call auxp1d_5pt(n3, func1d, aux)
              call derv1d_5pt(dz, n3, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n3, func1d, aux)
              call derv1d_3pt(dz, n3, aux, func1d, derv1d)
            endif
            gt(i,j,:) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
C
C  Redistribute data to be distributed by Z again
C
      call redistribXY(g(:,:,:,3),n1,n2l,n3l,gt,n1lf,n3,-1,nsm,
     .                 Node,Nodes)
C
C  Free local memory
C
      call de_alloc( ft, name='ft' )
      call de_alloc( gt, name='gt' )
#else
C
C  Differentiation in Y direction
C
      call re_alloc( func1d, 1, n2, name='func1d', routine='gradient' )
      call re_alloc( derv1d, 1, n2, name='derv1d', routine='gradient' )
      do k = 1,n3
         do i = 1,n1
            func1d(:) = f(i,:,k)
            if(nacc.eq.4) then
              call auxp1d_5pt(n2, func1d, aux)
              call derv1d_5pt(dy, n2, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n2, func1d, aux)
              call derv1d_3pt(dy, n2, aux, func1d, derv1d)
            endif
            g(i,:,k,2) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
C
C  Differentiation in Z direction
C
      call re_alloc( func1d, 1, n3, name='func1d', routine='gradient' )
      call re_alloc( derv1d, 1, n3, name='derv1d', routine='gradient' )
      do j = 1,n2
         do i = 1,n1
            func1d(:) = f(i,j,:)
            if(nacc.eq.4) then
              call auxp1d_5pt(n3, func1d, aux)
              call derv1d_5pt(dz, n3, aux, func1d, derv1d)
            else
              call auxp1d_3pt(n3, func1d, aux)
              call derv1d_3pt(dz, n3, aux, func1d, derv1d)
            endif
            g(i,j,:,3) = derv1d(:)
         enddo
      enddo
      call de_alloc( func1d, name='func1d' )
      call de_alloc( derv1d, name='derv1d' )
#endif

      call de_alloc( aux, name='aux' )

      return
      end subroutine gradient
!-------------------------------------
      subroutine operation_mix(w, f, Mesh, Meshl, nsm, oprtf)
      use precision, only : dp
      implicit none
      integer
     .  Mesh(3), Meshl(3), nsm
      real(dp)
     .  w(Meshl(1),Meshl(2),Meshl(3)), f(Meshl(1),Meshl(2),Meshl(3))
      real(dp)
     .  oprtf(Meshl(1),Meshl(2),Meshl(3))
!     Local variables
      integer  :: idim
      real(dp) :: gradf(Meshl(1),Meshl(2),Meshl(3),3)
!
      oprtf = 0.0_dp
!
      call gradient(f, Mesh, Meshl, nsm, gradf)
!
      do idim = 1,3
         gradf(:,:,:,idim) = w*gradf(:,:,:,idim)
      enddo
!
      call divergence(gradf,  Mesh, Meshl, nsm, oprtf)
      oprtf = -oprtf
!
      return
!
      end subroutine operation_mix
!-------------------------------------
!     subroutine conjgrad_mix(maxiter, tol, h, g,
!    .                        Mesh, Meshl, nsm, f, info)
!     use precision,   only : dp
!     use parallel,    only : Node
!     use m_variables, only : cgtol
!#ifdef MPI
!     use mpi_siesta
!#endif
!
!     implicit none
!
!     integer
!    .  maxiter, Mesh(3), Meshl(3), nsm
!     logical
!    .  info
!     real(dp)
!    .  tol
!     real(dp)
!    .  h(Meshl(1),Meshl(2),Meshl(3)), g(Meshl(1),Meshl(2),Meshl(3))
!     real(dp)
!    .  f(Meshl(1),Meshl(2),Meshl(3))
! Local variables
!     integer
!    .  iter
!     real(dp)
!    .  alpha, beta, rnorm, value, diff, offset
!     real(dp)
!    .  r(Meshl(1),Meshl(2),Meshl(3)), p(Meshl(1),Meshl(2),Meshl(3))
!     real(dp)
!    .  oprtf(Meshl(1),Meshl(2),Meshl(3)),
!    .  oprtp(Meshl(1),Meshl(2),Meshl(3))
!#ifdef MPI
!     integer :: MPIerror
!#endif
!
!     call operation_mix(h, f, Mesh, Meshl, nsm, oprtf)
!
!     r = -oprtf+g
!     p = r
!     iter = 0
!
!     call integration(r*r, Mesh, Meshl, nsm, rnorm)
!
!     if(Node.eq.0 .and. info) then
!       write(6,'(a)')    ' Conjugate-Gradient begins: '
!       write(6,'(a)')    ' ==========================='
!       write(6,'(a)')    ' iter        || Ax-b ||^2   '
!       write(6,'(a)')    ' ---------------------------'
!     end if
!
!     do while(iter .le. maxiter)
!
!        call operation_mix(h, p, Mesh, Meshl, nsm, oprtp)
!
!        call integration(p*oprtp, Mesh, Meshl, nsm, value)
!        if(dabs(value).lt.cgtol) exit
!        alpha = rnorm/value
!
!        f = f+alpha*p
! The following iteration expression is not numerically stable
!        r = r-alpha*oprtp
! The following gives feedback to the
! steepest descent direction
! but with sacrifice in the speed
!!       call operation_mix(h, f, Mesh, Meshl, nsm, oprtf)
!!       r = -oprtf+g
!        call integration(r*r, Mesh, Meshl, nsm, diff)
!        if(Node.eq.0 .and. info)
!    .    write(6,'(i5,5x,e15.7)') iter, dabs(diff)
!        if(dabs(diff) .le. tol) exit
!
!        beta = diff/rnorm
!        p = r+beta*p
!        rnorm = diff
!
!        iter = iter+1
!     end do
!
!     if(Node.eq.0 .and. info)
!    . write(6,'(a)') ' ---------------------------'
!
!     if(Node.eq.0) offset = f(1,1,1)
!#ifdef MPI
!     call MPI_Bcast(offset, 1, MPI_double_precision, 0,
!    .     MPI_Comm_World, MPIerror)
!#endif
!
!     f = f-offset
!
!     end subroutine conjgrad_mix
!-------------------------------------
!     subroutine conjgrad_mix_precond(cell, maxiter, tol, h, g,
!    .                                Mesh, Meshl, nsm, f, info)
!     use precision,   only : dp
!     use parallel,    only : Node
!     use m_variables, only : cgtol
!#ifdef MPI
!     use mpi_siesta
!#endif
!     implicit none
!
!     integer
!    .  maxiter, Mesh(3), Meshl(3), nsm
!     logical
!    .  info
!     real(dp)
!    .  cell(3,3)
!     real(dp)
!    .  tol
!     real(dp)
!    .  h(Meshl(1),Meshl(2),Meshl(3)), g(Meshl(1),Meshl(2),Meshl(3))
!     real(dp)
!    .  f(Meshl(1),Meshl(2),Meshl(3))
! Local variables
!     integer
!    .  iter, n1, n2, n3
!     real(dp)
!    .  alpha, beta, rnorm, value, diff, offset
!     real(dp)
!    .  r(Meshl(1),Meshl(2),Meshl(3)), 
!    .  z(Meshl(1),Meshl(2),Meshl(3)),
!    .  p(Meshl(1),Meshl(2),Meshl(3))
!     real(dp)
!    .  oprtf(Meshl(1),Meshl(2),Meshl(3)),
!    .  oprtp(Meshl(1),Meshl(2),Meshl(3))
!#ifdef MPI
!     integer :: MPIerror
!#endif
!
!     n1 = Meshl(1) ; n2 = Meshl(2) ; n3 = Meshl(3)
!
!     call operation_mix(h, f, Mesh, Meshl, nsm, oprtf)
!
!     r = -oprtf+g
!     call laplac_invers(cell, n1, n2, n3, Mesh, r, z, nsm)
!     p = z
!     iter = 0
!
!     call integration(r*z, Mesh, Meshl, nsm, rnorm)
!
!     if(Node.eq.0 .and. info) then
!       write(6,'(a)')    ' Conjugate-Gradient begins: '
!       write(6,'(a)')    ' ==========================='
!       write(6,'(a)')    ' iter        || Ax-b ||^2   '
!       write(6,'(a)')    ' ---------------------------'
!     end if
!
!     do while(iter .le. maxiter)
!
!        call operation_mix(h, p, Mesh, Meshl, nsm, oprtp)
!
!        call integration(p*oprtp, Mesh, Meshl, nsm, value)
!        if(dabs(value).lt.cgtol) exit
!        alpha = rnorm/value
!
!        f = f+alpha*p
! The following iteration expression is not numerically stable
!        r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!!       call operation_mix(h, f, Mesh, Meshl, nsm, oprtf)
!!       r = -oprtf+g
!        call laplac_invers(cell, n1, n2, n3, Mesh, r, z, nsm)
!        call integration(r*z, Mesh, Meshl, nsm, diff)
!        if(Node.eq.0 .and. info)
!    .    write(6,'(i5,5x,e15.7)') iter, dabs(diff)
!        if(dabs(diff) .le. tol) exit
!
!        beta = diff/rnorm
!        p = z+beta*p
!        rnorm = diff
!
!        iter = iter+1
!     end do
!
!     if(Node.eq.0 .and. info)
!    . write(6,'(a)') ' ---------------------------'
!
!     if(Node.eq.0) offset = f(1,1,1)
!#ifdef MPI
!     call MPI_Bcast(offset, 1, MPI_double_precision, 0,
!    .     MPI_Comm_World, MPIerror)
!#endif
!
!     f = f-offset
!
!     end subroutine conjgrad_mix_precond
!-------------------------------------
#ifdef MPI
      subroutine redistribXZ(f,n1,n2l,n3l,ft,n1lf,n2,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Z direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Z
C
C  Julian Gale, July 1999
C
      use precision,   only : dp
      use mpi_siesta
      use parallel,    only : ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n2, n2l, n1lf, n3l, Node, Nodes, idir, nsm
      real(dp)
     .  f(n1,n2l,n3l), ft(n1lf,n2,n3l)
C
C  Local variables
C
      integer
     .  BlockSizeY, BlockSizeYMax, jmin, jmax, jloc, n1lmax, NRemY,
     .  i, j, k, jl, kl, BNode, INode, SNode, jminS, jmaxS

      integer
     .  MPIerror, RequestR, RequestS, Status(MPI_Status_Size)

      integer, save ::
     .  ProcessorZ, Py, Pz, ZCommunicator, ZNode, ZNodes

      logical, save :: firsttimeZ = .true.

      real(dp), dimension(:,:,:), pointer  :: ftmp,ftmp2

      if (firsttimeZ) then
C
C  Determine processor grid coordinates
C
        ProcessorZ = Nodes/ProcessorY
        if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1
C
C  Group processors into subsets by Z
C
        call MPI_Comm_Split(MPI_Comm_World, Pz, Py, ZCommunicator,
     .    MPIerror)
        call MPI_Comm_Rank(ZCommunicator,ZNode,MPIerror)
        call MPI_Comm_Size(ZCommunicator,ZNodes,MPIerror)
        firsttimeZ = .false.
      endif

      BlockSizeY = ((n2/nsm)/ProcessorY)*nsm
      NRemY = (n2 - ProcessorY*BlockSizeY)/nsm
      if (NRemY.gt.0) then
        BlockSizeYMax = BlockSizeY + nsm
      else
        BlockSizeYMax = BlockSizeY
      endif
      n1lmax = ((n1-1)/ProcessorY) + 1
C
C  Work out local dimensions
C
      jmin = (Py-1)*BlockSizeY + nsm*min(Py-1,NRemY) + 1
      jmax = jmin + BlockSizeY - 1
      if (Py-1.lt.NRemY) jmax = jmax + nsm
      jmax = min(jmax,n2)
      jloc = jmax - jmin + 1
C
C  Allocate local memory and initialise
C
      nullify( ftmp )
      call re_alloc( ftmp, 1, n1lmax, 1, BlockSizeYMax, 1, n3l,
     &               name='ftmp', routine='redistribXZ' )
      nullify( ftmp2 )
      call re_alloc( ftmp2,1, n1lmax, 1, BlockSizeYMax, 1, n3l,
     &               name='ftmp2', routine='redistribXZ' )

      do i = 1,n3l
        do j = 1,BlockSizeYMax
          do k = 1,n1lmax
            ftmp(k,j,i) = 0.0_dp
            ftmp2(k,j,i) = 0.0_dp
          enddo
        enddo
      enddo

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = 1,n3l
          do j = jmin,jmax
            jl = j - jmin + 1
            kl = 0
            do k = 1+ZNode, n1, ZNodes
              kl = kl + 1
              ft(kl,j,i) = f(k,jl,i)
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorY-1
          BNode = (ZNode+INode)
          BNode = mod(BNode,ProcessorY)
          SNode = (ZNode-INode)
          SNode = mod(SNode+ProcessorY,ProcessorY)
C
C  Collect data to send
C
          do i = 1,n3l
            do jl = 1,jloc
              kl = 0
              do k = 1+BNode, n1, ZNodes
                kl = kl + 1
                ftmp(kl,jl,i) = f(k,jl,i)
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1),n1lmax*BlockSizeYMax*n3l,
     .     MPI_double_precision,SNode,1,ZCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1),n1lmax*BlockSizeYMax*n3l,
     .     MPI_double_precision,BNode,1,ZCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          jminS = SNode*BlockSizeY + nsm*min(SNode,NRemY) + 1
          jmaxS = jminS + BlockSizeY - 1
          if (SNode.lt.NRemY) jmaxS = jmaxS + nsm
          jmaxS = min(jmaxS,n2)
          do i = 1,n3l
            do j = jminS,jmaxS
              jl = j - jminS + 1
              kl = 0
              do k = 1+ZNode, n1, ZNodes
                kl = kl + 1
                ft(kl,j,i) = ftmp2(kl,jl,i)
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = 1,n3l
          do j = jmin,jmax
            jl = j - jmin + 1
            kl = 0
            do k = 1+ZNode, n1, ZNodes
              kl = kl + 1
              f(k,jl,i) = ft(kl,j,i) 
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorY-1
          BNode = (ZNode+INode)
          BNode = mod(BNode,ProcessorY)
          SNode = (ZNode-INode)
          SNode = mod(SNode+ProcessorY,ProcessorY)
C
C  Collect data to send
C
          jminS = SNode*BlockSizeY + nsm*min(SNode,NRemY) + 1
          jmaxS = jminS + BlockSizeY - 1
          if (SNode.lt.NRemY) jmaxS = jmaxS + nsm
          jmaxS = min(jmaxS,n2)
          do i = 1,n3l
            do j = jminS,jmaxS
              jl = j - jminS + 1
              kl = 0
              do k = 1+ZNode, n1, ZNodes
                kl = kl + 1
                ftmp(kl,jl,i) = ft(kl,j,i) 
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1),n1lmax*BlockSizeYMax*n3l,
     .     MPI_double_precision,BNode,1,ZCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1),n1lmax*BlockSizeYMax*n3l,
     .     MPI_double_precision,SNode,1,ZCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          do i = 1,n3l
            do jl = 1,jloc
              kl = 0
              do k = 1+BNode, n1, ZNodes
                kl = kl + 1
                f(k,jl,i) = ftmp2(kl,jl,i) 
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( ftmp2, name='ftmp2' )
      call de_alloc( ftmp, name='ftmp' )

      end subroutine redistribXZ

      subroutine redistribXY(f,n1,n2l,n3l,ft,n1lf,n3,idir,nsm,
     .  Node,Nodes)
C
C  This routine redistributes the data over the Nodes as needed
C  for the FFT routines between the arrays f and ft
C
C  Array f is distributed in the Y/Z direction while ft is distributed
C  in the X/Y direction
C
C  idir = direction of redistribution : > 0 => f->ft
C                                       < 0 => ft->f
C
C  Use the processor grid to divide communicator according to Y
C
C  Currently written in not the most efficient but simple way! 
C  Need to improve communication later.
C
C  Julian Gale, July 1999
C
      use precision,   only: dp
      use mpi_siesta
      use parallel,    only : ProcessorY
      use sys,         only : die
      use alloc,       only : re_alloc, de_alloc

      implicit none

C
C  Passed arguments
C
      integer
     .  n1, n3, n2l, n1lf, n3l, nsm, Node, Nodes, idir
      real(dp)
     .  f(n1,n2l,n3l), ft(n1lf,n2l,n3)
C
C  Local variables
C
      integer
     .  i, j, k, il, kl, BNode, BlockSizeZ, imin, imax, iloc,
     .  INode, n1lmax, SNode, iminS, imaxS, NRemZ, BlockSizeZMax

      integer
     .  MPIerror, RequestR, RequestS, Status(MPI_Status_Size)

      integer, save ::
     .  ProcessorZ, Py, Pz, YCommunicator, YNode, YNodes

      logical, save :: firsttimeY = .true.

      real(dp), dimension(:,:,:), pointer :: ftmp,ftmp2

      if (firsttimeY) then
C
C  Determine processor grid coordinates
C
        ProcessorZ = Nodes/ProcessorY
        if (ProcessorY*ProcessorZ.ne.Nodes) 
     $     call die('ERROR: ProcessorY must be a factor of the' //
     $     ' number of processors!')
        Py = (Node/ProcessorZ) + 1
        Pz = Node - (Py - 1)*ProcessorZ + 1
C
C  Group processors into subsets by Y
C
        call MPI_Comm_Split(MPI_Comm_World, Py, Pz, YCommunicator,
     .    MPIerror)
        call MPI_Comm_Rank(YCommunicator,YNode,MPIerror)
        call MPI_Comm_Size(YCommunicator,YNodes,MPIerror)
        firsttimeY = .false.
      endif

      BlockSizeZ = ((n3/nsm)/ProcessorZ)*nsm
      NRemZ = (n3 - ProcessorZ*BlockSizeZ)/nsm
      if (NRemZ.gt.0) then
        BlockSizeZMax = BlockSizeZ + nsm
      else
        BlockSizeZMax = BlockSizeZ
      endif
      n1lmax = ((n1-1)/ProcessorZ) + 1
C
C  Allocate local memory and initialise
C
      nullify( ftmp )
      call re_alloc( ftmp, 1, n1lmax, 1, n2l, 1, BlockSizeZMax,
     &               name='ftmp', routine='redistribXY' )
      nullify( ftmp2 )
      call re_alloc( ftmp2,1, n1lmax, 1, n2l, 1, BlockSizeZMax,
     &               name='ftmp2', routine='redistribXY' )

      do i = 1,BlockSizeZMax
        do j = 1,n2l
          do k = 1,n1lmax
            ftmp(k,j,i) = 0.0_dp
            ftmp2(k,j,i) = 0.0_dp
          enddo
        enddo
      enddo
C
C  Work out local dimensions
C
      imin = (Pz-1)*BlockSizeZ + nsm*min(Pz-1,NRemZ) + 1
      imax = imin + BlockSizeZ - 1
      if (Pz-1.lt.NRemZ) imax = imax + nsm
      imax = min(imax,n3)
      iloc = imax - imin + 1

      if (idir.gt.0) then
C***********
C  F -> FT *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = imin,imax
          il = i - imin + 1
          do j = 1,n2l
            kl = 0
            do k = 1+YNode, n1, YNodes
              kl = kl + 1
              ft(kl,j,i) = f(k,j,il)
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorZ-1
          BNode = (YNode+INode)
          BNode = mod(BNode,ProcessorZ)
          SNode = (YNode-INode)
          SNode = mod(SNode+ProcessorZ,ProcessorZ)
C
C  Collect data to send
C
          do il = 1,iloc
            do j = 1,n2l
              kl = 0
              do k = 1+BNode, n1, YNodes
                kl = kl + 1
                ftmp(kl,j,il) = f(k,j,il)
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1),n1lmax*n2l*BlockSizeZMax,
     .     MPI_double_precision,SNode,1,YCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1),n1lmax*n2l*BlockSizeZMax,
     .     MPI_double_precision,BNode,1,YCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          iminS = SNode*BlockSizeZ + nsm*min(SNode,NRemZ) + 1
          imaxS = iminS + BlockSizeZ - 1
          if (SNode.lt.NRemZ) imaxS = imaxS + nsm
          imaxS = min(imaxS,n3)
          do i = iminS,imaxS
            il = i - iminS + 1
            do j = 1,n2l
              kl = 0
              do k = 1+YNode, n1, YNodes
                kl = kl + 1
                ft(kl,j,i) = ftmp2(kl,j,il)
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      elseif (idir.lt.0) then
C***********
C  FT -> F *
C***********
C
C  Handle transfer of terms which are purely local
C
        do i = imin,imax
          il = i - imin + 1
          do j = 1,n2l
            kl = 0
            do k = 1+YNode, n1, YNodes
              kl = kl + 1
              f(k,j,il) = ft(kl,j,i)
            enddo
          enddo
        enddo
C
C  Loop over all Node-Node vectors exchanging local data
C
        do INode = 1,ProcessorZ-1
          BNode = (YNode+INode)
          BNode = mod(BNode,ProcessorZ)
          SNode = (YNode-INode)
          SNode = mod(SNode+ProcessorZ,ProcessorZ)
C
C  Collect data to send
C
          iminS = SNode*BlockSizeZ + nsm*min(SNode,NRemZ) + 1
          imaxS = iminS + BlockSizeZ - 1
          if (SNode.lt.NRemZ) imaxS = imaxS + nsm
          imaxS = min(imaxS,n3)
          do i = iminS,imaxS
            il = i - iminS + 1
            do j = 1,n2l
              kl = 0
              do k = 1+YNode, n1, YNodes
                kl = kl + 1
                ftmp(kl,j,il) = ft(kl,j,i)
              enddo
            enddo
          enddo
C
C  Exchange data - send to right and receive from left
C
          call MPI_IRecv(ftmp2(1,1,1),n1lmax*n2l*BlockSizeZMax,
     .     MPI_double_precision,BNode,1,YCommunicator,RequestR,MPIerror)
          call MPI_ISend(ftmp(1,1,1),n1lmax*n2l*BlockSizeZMax,
     .     MPI_double_precision,SNode,1,YCommunicator,RequestS,MPIerror)
C
C  Wait for receive to complete
C
          call MPI_Wait(RequestR,Status,MPIerror)
C
C  Place received data into correct array
C
          do il = 1,iloc
            do j = 1,n2l
              kl = 0
              do k = 1+BNode, n1, YNodes
                kl = kl + 1
                f(k,j,il) = ftmp2(kl,j,il)
              enddo
            enddo
          enddo
C
C  Wait for send to complete
C
          call MPI_Wait(RequestS,Status,MPIerror)
        enddo
      endif
C
C  Free local memory
C
      call de_alloc( ftmp2, name='ftmp2' )
      call de_alloc( ftmp, name='ftmp' )

      return
      end subroutine redistribXY
#endif

      end module m_numeric3d
