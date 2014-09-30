 module m_numrecipe
 use precision, only : dp
 implicit none
!
 contains
!
 subroutine auxp1d_3pt(nx, func, aux)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: func(nx)
 real(dp), intent(out) :: aux(2)

 aux(1) = func(nx)
 aux(2) = func(1)

 end subroutine auxp1d_3pt
!
 subroutine auxp1d_5pt(nx, func, aux)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: func(nx)
 real(dp), intent(out) :: aux(4)

 aux(1:2) = func(nx-1:nx)
 aux(3:4) = func(1:2)

 end subroutine auxp1d_5pt
!
 subroutine auxp3d_3pt(nx, ny, nz, func, auxx, auxy, auxz)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)

 integer :: ix, iy, iz

 do iz = 1, nz
    do iy = 1, ny
       call auxp1d_3pt(nx, func(:,iy,iz), auxx(:,iy,iz))
    end do
 end do
 do iz = 1, nz
    do ix = 1, nx
       call auxp1d_3pt(ny, func(ix,:,iz), auxy(ix,:,iz))
    end do
 end do
 do iy = 1, ny
    do ix = 1, nx
       call auxp1d_3pt(nz, func(ix,iy,:), auxz(ix,iy,:))
    end do
 end do

 end subroutine auxp3d_3pt
!
 subroutine auxp3d_5pt(nx, ny, nz, func, auxx, auxy, auxz)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)

 integer :: ix, iy, iz

 do iz = 1, nz
    do iy = 1, ny
       call auxp1d_5pt(nx, func(:,iy,iz), auxx(:,iy,iz))
    end do
 end do
 do iz = 1, nz
    do ix = 1, nx
       call auxp1d_5pt(ny, func(ix,:,iz), auxy(ix,:,iz))
    end do
 end do
 do iy = 1, ny
    do ix = 1, nx
       call auxp1d_5pt(nz, func(ix,iy,:), auxz(ix,iy,:))
    end do
 end do

 end subroutine auxp3d_5pt
!
 subroutine auxv3d_3pt(nx, ny, nz, vect, auxx, auxy, auxz)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: vect(nx,ny,nz,3)
 real(dp), intent(out) :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)

 integer :: ix, iy, iz

 do iz = 1, nz
    do iy = 1, ny
       auxx(1,iy,iz) = vect(nx,iy,iz,1)
       auxx(2,iy,iz) = vect(1,iy,iz,1)
    end do
 end do
 do iz = 1, nz
    do ix = 1, nx
       auxy(ix,1,iz) = vect(ix,ny,iz,2)
       auxy(ix,2,iz) = vect(ix,1,iz,2)
    end do
 end do
 do iy = 1, ny
    do ix = 1, nx
       auxz(ix,iy,1) = vect(ix,iy,nz,3)
       auxz(ix,iy,2) = vect(ix,iy,1,3)
    end do
 end do

 end subroutine auxv3d_3pt
!
 subroutine auxv3d_5pt(nx, ny, nz, vect, auxx, auxy, auxz)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: vect(nx,ny,nz,3)
 real(dp), intent(out) :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)

 integer :: ix, iy, iz

 do iz = 1, nz
    do iy = 1, ny
       auxx(1:2,iy,iz) = vect(nx-1:nx,iy,iz,1)
       auxx(3:4,iy,iz) = vect(1:2,iy,iz,1)
    end do
 end do
 do iz = 1, nz
    do ix = 1, nx
       auxy(ix,1:2,iz) = vect(ix,ny-1:ny,iz,2)
       auxy(ix,3:4,iz) = vect(ix,1:2,iz,2)
    end do
 end do
 do iy = 1, ny
    do ix = 1, nx
       auxz(ix,iy,1:2) = vect(ix,iy,nz-1:nz,3)
       auxz(ix,iy,3:4) = vect(ix,iy,1:2,3)
    end do
 end do

 end subroutine auxv3d_5pt
!
 subroutine derv1d_3pt(dx, k, aux, func, derv)
 implicit none
 integer,  intent(in)  :: k
 real(dp), intent(in)  :: dx, aux(2)
 real(dp), intent(in)  :: func(k)
 real(dp), intent(out) :: derv(k)
! Local variables
 real(dp), parameter :: coef(3) = (/-1.d0, 0.d0, 1.d0/)
 integer :: ix, iy, ip
 real(dp), allocatable :: ftmp(:)

 if(allocated(ftmp)) deallocate(ftmp)
 allocate(ftmp(k+2))

 ftmp(1) = aux(1)
 ftmp(2:k+1) = func(:)
 ftmp(k+2) = aux(2)

 do ix = 1, k
    derv(ix) = 0.d0
    do ip = 1, 3
       iy = ix+ip-1
       derv(ix) = derv(ix)+coef(ip)*ftmp(iy)
    end do
 end do

 derv = 1.d0/(2.d0*dx)*derv

 if(allocated(ftmp)) deallocate(ftmp)

 end subroutine derv1d_3pt
!
 subroutine derv1d_5pt(dx, k, aux, func, derv)
 implicit none
 integer,  intent(in)  :: k
 real(dp), intent(in)  :: dx, aux(4)
 real(dp), intent(in)  :: func(k)
 real(dp), intent(out) :: derv(k)
! Local variables
 real(dp), parameter :: coef(5) = (/1.d0, -8.d0, 0.d0, 8.d0, -1.d0/)
 integer :: ix, iy, ip
 real(dp), allocatable :: ftmp(:)

 if(allocated(ftmp)) deallocate(ftmp)
 allocate(ftmp(k+4))

 ftmp(1:2) = aux(1:2)
 ftmp(3:k+2) = func(:)
 ftmp(k+3:k+4) = aux(3:4)

 do ix = 1, k
    derv(ix) = 0.d0
    do ip = 1, 5
       iy = ix+ip-1
       derv(ix) = derv(ix)+coef(ip)*ftmp(iy)
    end do
 end do

 derv = 1.d0/(12.d0*dx)*derv

 if(allocated(ftmp)) deallocate(ftmp)

 end subroutine derv1d_5pt
!
 subroutine grad3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, func, grad)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: grad(nx,ny,nz,3)

 integer :: ix, iy, iz
 real(dp), allocatable :: func1d(:), derv1d(:), aux1d(:)

!x part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nx))
 allocate(derv1d(nx))
 allocate(aux1d(2))
!
 do iz = 1, nz
    do iy = 1, ny
       func1d(:) = func(:,iy,iz)
        aux1d(:) = auxx(:,iy,iz)
       call derv1d_3pt(dx, nx, aux1d, func1d, derv1d)
       grad(:,iy,iz,1) = derv1d(:)
    end do
 end do
!y part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(ny))
 allocate(derv1d(ny))
 allocate(aux1d(2))
!
 do iz = 1, nz
    do ix = 1, nx
       func1d(:) = func(ix,:,iz)
        aux1d(:) = auxy(ix,:,iz)
       call derv1d_3pt(dy, ny, aux1d, func1d, derv1d)
       grad(ix,:,iz,2) = derv1d(:)
    end do
 end do
!z part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nz))
 allocate(derv1d(nz))
 allocate(aux1d(2))
!
 do iy = 1, ny
    do ix = 1, nx
       func1d(:) = func(ix,iy,:)
        aux1d(:) = auxz(ix,iy,:)
       call derv1d_3pt(dz, nz, aux1d, func1d, derv1d)
       grad(ix,iy,:,3) = derv1d(:)
    end do
 end do

 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)

 end subroutine grad3d_3pt
!
 subroutine grad3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, func, grad)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: grad(nx,ny,nz,3)

 integer :: ix, iy, iz
 real(dp), allocatable :: func1d(:), derv1d(:), aux1d(:)

!x part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nx))
 allocate(derv1d(nx))
 allocate(aux1d(4))
!
 do iz = 1, nz
    do iy = 1, ny
       func1d(:) = func(:,iy,iz)
        aux1d(:) = auxx(:,iy,iz)
       call derv1d_5pt(dx, nx, aux1d, func1d, derv1d)
       grad(:,iy,iz,1) = derv1d(:)
    end do
 end do
!y part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(ny))
 allocate(derv1d(ny))
 allocate(aux1d(4))
!
 do iz = 1, nz
    do ix = 1, nx
       func1d(:) = func(ix,:,iz)
        aux1d(:) = auxy(ix,:,iz)
       call derv1d_5pt(dy, ny, aux1d, func1d, derv1d)
       grad(ix,:,iz,2) = derv1d(:)
    end do
 end do
!z part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nz))
 allocate(derv1d(nz))
 allocate(aux1d(4))
!
 do iy = 1, ny
    do ix = 1, nx
       func1d(:) = func(ix,iy,:)
        aux1d(:) = auxz(ix,iy,:)
       call derv1d_5pt(dz, nz, aux1d, func1d, derv1d)
       grad(ix,iy,:,3) = derv1d(:)
    end do
 end do

 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)

 end subroutine grad3d_5pt
!
 subroutine divg3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, vect, divg)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)
 real(dp), intent(in)  :: vect(nx,ny,nz,3)
 real(dp), intent(out) :: divg(nx,ny,nz)

 integer :: ix, iy, iz
 real(dp), allocatable :: func1d(:), derv1d(:), aux1d(:)

 divg = 0.d0
!x part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nx))
 allocate(derv1d(nx))
 allocate(aux1d(2))
!
 do iz = 1, nz
    do iy = 1, ny
       func1d(:) = vect(:,iy,iz,1)
        aux1d(:) = auxx(:,iy,iz)
       call derv1d_3pt(dx, nx, aux1d, func1d, derv1d)
       divg(:,iy,iz) = divg(:,iy,iz) + derv1d(:)
    end do
 end do
!y part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(ny))
 allocate(derv1d(ny))
 allocate(aux1d(2))
!
 do iz = 1, nz
    do ix = 1, nx
       func1d(:) = vect(ix,:,iz,2)
        aux1d(:) = auxy(ix,:,iz)
       call derv1d_3pt(dy, ny, aux1d, func1d, derv1d)
       divg(ix,:,iz) = divg(ix,:,iz) + derv1d(:)
    end do
 end do
!z part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nz))
 allocate(derv1d(nz))
 allocate(aux1d(2))
!
 do iy = 1, ny
    do ix = 1, nx
       func1d(:) = vect(ix,iy,:,3)
        aux1d(:) = auxz(ix,iy,:)
       call derv1d_3pt(dz, nz, aux1d, func1d, derv1d)
       divg(ix,iy,:) = divg(ix,iy,:) + derv1d(:)
    end do
 end do

 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)

 end subroutine divg3d_3pt
!
 subroutine divg3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, vect, divg)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)
 real(dp), intent(in)  :: vect(nx,ny,nz,3)
 real(dp), intent(out) :: divg(nx,ny,nz)

 integer :: ix, iy, iz
 real(dp), allocatable :: func1d(:), derv1d(:), aux1d(:)

 divg = 0.d0
!x part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nx))
 allocate(derv1d(nx))
 allocate(aux1d(4))
!
 do iz = 1, nz
    do iy = 1, ny
       func1d(:) = vect(:,iy,iz,1)
        aux1d(:) = auxx(:,iy,iz)
       call derv1d_5pt(dx, nx, aux1d, func1d, derv1d)
       divg(:,iy,iz) = divg(:,iy,iz) + derv1d(:)
    end do
 end do
!y part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(ny))
 allocate(derv1d(ny))
 allocate(aux1d(4))
!
 do iz = 1, nz
    do ix = 1, nx
       func1d(:) = vect(ix,:,iz,2)
        aux1d(:) = auxy(ix,:,iz)
       call derv1d_5pt(dy, ny, aux1d, func1d, derv1d)
       divg(ix,:,iz) = divg(ix,:,iz) + derv1d(:)
    end do
 end do
!z part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nz))
 allocate(derv1d(nz))
 allocate(aux1d(4))
!
 do iy = 1, ny
    do ix = 1, nx
       func1d(:) = vect(ix,iy,:,3)
        aux1d(:) = auxz(ix,iy,:)
       call derv1d_5pt(dz, nz, aux1d, func1d, derv1d)
       divg(ix,iy,:) = divg(ix,iy,:) + derv1d(:)
    end do
 end do

 if(allocated(func1d)) deallocate(func1d)
 if(allocated(derv1d)) deallocate(derv1d)
 if(allocated(aux1d))  deallocate(aux1d)

 end subroutine divg3d_5pt
!
 subroutine lapl1d_3pt(dx, k, aux, func, lapl)
 implicit none

 integer,  intent(in)  :: k
 real(dp), intent(in)  :: dx, aux(2)
 real(dp), intent(in)  :: func(k)
 real(dp), intent(out) :: lapl(k)

 real(dp), parameter :: coef(3) = (/1.d0, -2.d0, 1.d0/)
 integer  :: ix, iy, ip
 real(dp), allocatable :: ftmp(:)

 if(allocated(ftmp)) deallocate(ftmp)
 allocate(ftmp(k+2))

 ftmp(1) = aux(1)
 ftmp(2:k+1) = func(:)
 ftmp(k+2) = aux(2)

 do ix = 1, k
    lapl(ix) = 0.d0
    do ip = 1, 3
       iy = ix+ip-1
       lapl(ix) = lapl(ix)+coef(ip)*ftmp(iy)
    end do
 end do

 lapl = 1.d0/(dx*dx)*lapl

 if(allocated(ftmp)) deallocate(ftmp)

 end subroutine lapl1d_3pt
!
 subroutine lapl1d_5pt(dx, k, aux, func, lapl)
 implicit none

 integer,  intent(in)  :: k
 real(dp), intent(in)  :: dx, aux(4)
 real(dp), intent(in)  :: func(k)
 real(dp), intent(out) :: lapl(k)

 real(dp), parameter :: coef(5) = (/-1.d0, 16.d0, -30.d0, 16.d0, -1.d0/)
 integer :: ix, iy, ip
 real(dp), allocatable :: ftmp(:)

 if(allocated(ftmp)) deallocate(ftmp)
 allocate(ftmp(k+4))

 ftmp(1:2) = aux(1:2)
 ftmp(3:k+2) = func(:)
 ftmp(k+3:k+4) = aux(3:4)

 do ix = 1, k
    lapl(ix) = 0.d0
    do ip = 1, 5
       iy = ix+ip-1
       lapl(ix) = lapl(ix)+coef(ip)*ftmp(iy)
    end do
 end do

 lapl = 1.d0/(12.d0*dx*dx)*lapl

 if(allocated(ftmp)) deallocate(ftmp)

 end subroutine lapl1d_5pt
!
 subroutine lapl3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, func, lapl)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: lapl(nx,ny,nz)

 integer :: ix, iy, iz
 real(dp), allocatable :: func1d(:), lapl1d(:), aux1d(:)

 lapl = 0.d0
!x part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nx))
 allocate(lapl1d(nx))
 allocate(aux1d(2))
!
 do iz = 1, nz
    do iy = 1, ny
       func1d(:) = func(:,iy,iz)
        aux1d(:) = auxx(:,iy,iz)
       call lapl1d_3pt(dx, nx, aux1d, func1d, lapl1d)
       lapl(:,iy,iz) = lapl(:,iy,iz) + lapl1d(:)
    end do
 end do
!y part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(ny))
 allocate(lapl1d(ny))
 allocate(aux1d(2))
!
 do iz = 1, nz
    do ix = 1, nx
       func1d(:) = func(ix,:,iz)
        aux1d(:) = auxy(ix,:,iz)
       call lapl1d_3pt(dy, ny, aux1d, func1d, lapl1d)
       lapl(ix,:,iz) = lapl(ix,:,iz) + lapl1d(:)
    end do
 end do
!z part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nz))
 allocate(lapl1d(nz))
 allocate(aux1d(2))
!
 do iy = 1, ny
    do ix = 1, nx
       func1d(:) = func(ix,iy,:)
        aux1d(:) = auxz(ix,iy,:)
       call lapl1d_3pt(dz, nz, aux1d, func1d, lapl1d)
       lapl(ix,iy,:) = lapl(ix,iy,:) + lapl1d(:)
    end do
 end do

 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)

 end subroutine lapl3d_3pt
!
 subroutine lapl3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, func, lapl)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: lapl(nx,ny,nz)

 integer :: ix, iy, iz
 real(dp), allocatable :: func1d(:), lapl1d(:), aux1d(:)

 lapl = 0.d0
!x part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nx))
 allocate(lapl1d(nx))
 allocate(aux1d(4))
!
 do iz = 1, nz
    do iy = 1, ny
       func1d(:) = func(:,iy,iz)
        aux1d(:) = auxx(:,iy,iz)
       call lapl1d_5pt(dx, nx, aux1d, func1d, lapl1d)
       lapl(:,iy,iz) = lapl(:,iy,iz) + lapl1d(:)
    end do
 end do
!y part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(ny))
 allocate(lapl1d(ny))
 allocate(aux1d(4))
!
 do iz = 1, nz
    do ix = 1, nx
       func1d(:) = func(ix,:,iz)
        aux1d(:) = auxy(ix,:,iz)
       call lapl1d_5pt(dy, ny, aux1d, func1d, lapl1d)
       lapl(ix,:,iz) = lapl(ix,:,iz) + lapl1d(:)
    end do
 end do
!z part
 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)
 allocate(func1d(nz))
 allocate(lapl1d(nz))
 allocate(aux1d(4))
!
 do iy = 1, ny
    do ix = 1, nx
       func1d(:) = func(ix,iy,:)
        aux1d(:) = auxz(ix,iy,:)
       call lapl1d_5pt(dz, nz, aux1d, func1d, lapl1d)
       lapl(ix,iy,:) = lapl(ix,iy,:) + lapl1d(:)
    end do
 end do

 if(allocated(func1d)) deallocate(func1d)
 if(allocated(lapl1d)) deallocate(lapl1d)
 if(allocated(aux1d))  deallocate(aux1d)

 end subroutine lapl3d_5pt
!
 subroutine oprt1d_3pt(dx, nx, auxf, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: dx
 real(dp), intent(in)  :: auxf(2)
 real(dp), intent(in)  :: f(nx)
 real(dp), intent(out) :: oprtf(nx)

 call lapl1d_3pt(dx, nx, auxf, f, oprtf)
 oprtf = -oprtf

 end subroutine oprt1d_3pt
!
 subroutine oprt1d_5pt(dx, nx, auxf, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: dx
 real(dp), intent(in)  :: auxf(4)
 real(dp), intent(in)  :: f(nx)
 real(dp), intent(out) :: oprtf(nx)

 call lapl1d_5pt(dx, nx, auxf, f, oprtf)
 oprtf = -oprtf

 end subroutine oprt1d_5pt
!
 subroutine oprt1d_3pt_mix(dx, nx, w, auxf, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: dx
 real(dp), intent(in)  :: auxf(2)
 real(dp), intent(in)  :: w(nx), f(nx)
 real(dp), intent(out) :: oprtf(nx)

 real(dp) :: aux(2)
 real(dp) :: dervf(nx)

 call derv1d_3pt(dx, nx, auxf, f, dervf)
 dervf = w*dervf

 call auxp1d_3pt(nx, dervf, aux)
 call derv1d_3pt(dx, nx, aux, dervf, oprtf) 

 oprtf = -oprtf

 end subroutine oprt1d_3pt_mix
!
 subroutine oprt1d_5pt_mix(dx, nx, w, auxf, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: dx
 real(dp), intent(in)  :: auxf(4)
 real(dp), intent(in)  :: w(nx), f(nx)
 real(dp), intent(out) :: oprtf(nx)

 real(dp) :: aux(4)
 real(dp) :: dervf(nx)

 call derv1d_5pt(dx, nx, auxf, f, dervf)
 dervf = w*dervf

 call auxp1d_5pt(nx, dervf, aux)
 call derv1d_5pt(dx, nx, aux, dervf, oprtf)

 oprtf = -oprtf

 end subroutine oprt1d_5pt_mix
!
 subroutine oprt3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)
 real(dp), intent(in)  :: f(nx,ny,nz)
 real(dp), intent(out) :: oprtf(nx,ny,nz)

 call lapl3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)
 oprtf = -oprtf

 end subroutine oprt3d_3pt
!
 subroutine oprt3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)
 real(dp), intent(in)  :: f(nx,ny,nz)
 real(dp), intent(out) :: oprtf(nx,ny,nz)

 call lapl3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)
 oprtf = -oprtf

 end subroutine oprt3d_5pt
!
 subroutine oprt3d_3pt_mix(dx, dy, dz, nx, ny, nz, w, auxfx, auxfy, auxfz, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxfx(2,ny,nz), auxfy(nx,2,nz), auxfz(nx,ny,2)
 real(dp), intent(in)  :: w(nx,ny,nz), f(nx,ny,nz)
 real(dp), intent(out) :: oprtf(nx,ny,nz)

 integer  :: ix, iy, iz
 real(dp) :: auxf1d(2)
 real(dp), allocatable :: w1d(:), f1d(:), oprtf1d(:)

 oprtf = 0.d0

! x part
 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)
 allocate(w1d(nx),f1d(nx),oprtf1d(nx))

 do iz = 1, nz
    do iy = 1, ny
       w1d = w(:,iy,iz)
       f1d = f(:,iy,iz)
       auxf1d = auxfx(:,iy,iz)
       call oprt1d_3pt_mix(dx, nx, w1d, auxf1d, f1d, oprtf1d)
       oprtf(:,iy,iz) = oprtf(:,iy,iz)+oprtf1d
    end do
 end do
! y part
 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)
 allocate(w1d(ny),f1d(ny),oprtf1d(ny))

 do iz = 1, nz
    do ix = 1, nx
       w1d = w(ix,:,iz)
       f1d = f(ix,:,iz)
       auxf1d = auxfy(ix,:,iz)
       call oprt1d_3pt_mix(dy, ny, w1d, auxf1d, f1d, oprtf1d)
       oprtf(ix,:,iz) = oprtf(ix,:,iz)+oprtf1d
    end do
 end do
! z part
 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)
 allocate(w1d(nz),f1d(nz),oprtf1d(nz))

 do iy = 1, ny
    do ix = 1, nx
       w1d = w(ix,iy,:)
       f1d = f(ix,iy,:)
       auxf1d = auxfz(ix,iy,:)
       call oprt1d_3pt_mix(dz, nz, w1d, auxf1d, f1d, oprtf1d)
       oprtf(ix,iy,:) = oprtf(ix,iy,:)+oprtf1d
    end do
 end do

 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)

 end subroutine oprt3d_3pt_mix
!
 subroutine oprt3d_5pt_mix(dx, dy, dz, nx, ny, nz, w, auxfx, auxfy, auxfz, f, oprtf)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: auxfx(4,ny,nz), auxfy(nx,4,nz), auxfz(nx,ny,4)
 real(dp), intent(in)  :: w(nx,ny,nz), f(nx,ny,nz)
 real(dp), intent(out) :: oprtf(nx,ny,nz)

 integer  :: ix, iy, iz
 real(dp) :: auxf1d(4)
 real(dp), allocatable :: w1d(:), f1d(:), oprtf1d(:)

 oprtf = 0.d0

! x part
 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)
 allocate(w1d(nx),f1d(nx),oprtf1d(nx))

 do iz = 1, nz
    do iy = 1, ny
       w1d = w(:,iy,iz)
       f1d = f(:,iy,iz)
       auxf1d = auxfx(:,iy,iz)
       call oprt1d_5pt_mix(dx, nx, w1d, auxf1d, f1d, oprtf1d)
       oprtf(:,iy,iz) = oprtf(:,iy,iz)+oprtf1d
    end do
 end do
! y part
 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)
 allocate(w1d(ny),f1d(ny),oprtf1d(ny))

 do iz = 1, nz
    do ix = 1, nx
       w1d = w(ix,:,iz)
       f1d = f(ix,:,iz)
       auxf1d = auxfy(ix,:,iz)
       call oprt1d_5pt_mix(dy, ny, w1d, auxf1d, f1d, oprtf1d)
       oprtf(ix,:,iz) = oprtf(ix,:,iz)+oprtf1d
    end do
 end do
! z part
 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)
 allocate(w1d(nz),f1d(nz),oprtf1d(nz))

 do iy = 1, ny
    do ix = 1, nx
       w1d = w(ix,iy,:)
       f1d = f(ix,iy,:)
       auxf1d = auxfz(ix,iy,:)
       call oprt1d_5pt_mix(dz, nz, w1d, auxf1d, f1d, oprtf1d)
       oprtf(ix,iy,:) = oprtf(ix,iy,:)+oprtf1d
    end do
 end do

 if(allocated(w1d))     deallocate(w1d)
 if(allocated(f1d))     deallocate(f1d)
 if(allocated(oprtf1d)) deallocate(oprtf1d)

 end subroutine oprt3d_5pt_mix
!
 subroutine cjgd1d_3pt(dx, nx, maxiter, tol, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, maxiter
 real(dp), intent(in)    :: dx, tol
 real(dp), intent(in)    :: g(nx)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx)

 integer  :: ix, iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx), p(nx)
 real(dp) :: oprtf(nx), oprtp(nx)
 real(dp) :: aux(2)

 ierr = 1

 call auxp1d_3pt(nx, f, aux)
 call oprt1d_3pt(dx, nx, aux, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg1d(dx, nx, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    '============================'
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    '----------------------------'
 end if

 do while(iter .le. maxiter)
    call auxp1d_3pt(nx, p, aux)
    call oprt1d_3pt(dx, nx, aux, p, oprtp)

    call intg1d(dx, nx, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp1d_3pt(nx, f, aux)
!   call oprt1d_3pt(dx, nx, aux, f, oprtf)
!   r = -oprtf+g

    call intg1d(dx, nx, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') '----------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1)
 f = f-offset

 end subroutine cjgd1d_3pt
!
 subroutine cjgd1d_5pt(dx, nx, maxiter, tol, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, maxiter
 real(dp), intent(in)    :: dx, tol
 real(dp), intent(in)    :: g(nx)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx)

 integer  :: ix, iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx), p(nx)
 real(dp) :: oprtf(nx), oprtp(nx)
 real(dp) :: aux(4)

 ierr = 1

 call auxp1d_5pt(nx, f, aux)
 call oprt1d_5pt(dx, nx, aux, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg1d(dx, nx, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    '============================'
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    '----------------------------'
 end if

 do while(iter .le. maxiter)
    call auxp1d_5pt(nx, p, aux)
    call oprt1d_5pt(dx, nx, aux, p, oprtp)

    call intg1d(dx, nx, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp1d(nx, f, aux)
!   call oprt1d(dx, nx, aux, f, oprtf)
!   r = -oprtf+g

    call intg1d(dx, nx, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') '----------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1)
 f = f-offset

 end subroutine cjgd1d_5pt
!
 subroutine cjgd1d_3pt_mix(dx, nx, maxiter, tol, h, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, maxiter
 real(dp), intent(in)    :: dx, tol
 real(dp), intent(in)    :: h(nx), g(nx)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx)

 integer  :: ix, iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx), p(nx)
 real(dp) :: oprtf(nx), oprtp(nx)
 real(dp) :: auxf(2), auxp(2)

 ierr = 1

 call auxp1d_3pt(nx, f, auxf)
 call oprt1d_3pt_mix(dx, nx, h, auxf, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg1d(dx, nx, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    '============================'
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    '----------------------------'
 end if

 do while(iter .le. maxiter)
    call auxp1d_3pt(nx, p, auxp)
    call oprt1d_3pt_mix(dx, nx, h, auxp, p, oprtp)

    call intg1d(dx, nx, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp1d(nx, f, auxf)
!   call oprt1d_mix(dx, nx, auxh, h, auxf, f, oprtf)
!   r = -oprtf+g

    call intg1d(dx, nx, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') '----------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1)
 f = f-offset

 end subroutine cjgd1d_3pt_mix
!
 subroutine cjgd1d_5pt_mix(dx, nx, maxiter, tol, h, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, maxiter
 real(dp), intent(in)    :: dx, tol
 real(dp), intent(in)    :: h(nx), g(nx)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx)

 integer  :: ix, iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx), p(nx)
 real(dp) :: oprtf(nx), oprtp(nx)
 real(dp) :: auxf(4), auxp(4)

 ierr = 1

 call auxp1d_5pt(nx, f, auxf)
 call oprt1d_5pt_mix(dx, nx, h, auxf, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg1d(dx, nx, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    '============================'
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    '----------------------------'
 end if

 do while(iter .le. maxiter)
    call auxp1d_5pt(nx, p, auxp)
    call oprt1d_5pt_mix(dx, nx, h, auxp, p, oprtp)

    call intg1d(dx, nx, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp1d_5pt(nx, f, auxf)
!   call oprt1d_5pt_mix(dx, nx, h, auxf, f, oprtf)
!   r = -oprtf+g

    call intg1d(dx, nx, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') '----------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1)
 f = f-offset

 end subroutine cjgd1d_5pt_mix
!
 subroutine cjgd3d_3pt(dx, dy, dz, nx, ny, nz, maxiter, tol, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, ny, nz, maxiter
 real(dp), intent(in)    :: dx, dy, dz, tol
 real(dp), intent(in)    :: g(nx,ny,nz)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx,ny,nz)

 integer  :: iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx,ny,nz), p(nx,ny,nz)
 real(dp) :: oprtf(nx,ny,nz), oprtp(nx,ny,nz)
 real(dp) :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)

 ierr = 1

 call auxp3d_3pt(nx, ny, nz, f, auxx, auxy, auxz)
 call oprt3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg3d(dx, dy, dz, nx, ny, nz, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    ' ==========================='
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    ' ---------------------------'
 end if

 do while(iter .le. maxiter)
    
    call auxp3d_3pt(nx, ny, nz, p, auxx, auxy, auxz)
    call oprt3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, p, oprtp)

    call intg3d(dx, dy, dz, nx, ny, nz, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp3d_3pt(nx, ny, nz, f, auxx, auxy, auxz)
!   call oprt3d_3pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)
!   r = -oprtf+g

    call intg3d(dx, dy, dz, nx, ny, nz, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') ' ---------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1,1,1)
 f = f-offset

 end subroutine cjgd3d_3pt
!
 subroutine cjgd3d_5pt(dx, dy, dz, nx, ny, nz, maxiter, tol, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, ny, nz, maxiter
 real(dp), intent(in)    :: dx, dy, dz, tol
 real(dp), intent(in)    :: g(nx,ny,nz)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx,ny,nz)

 integer  :: iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx,ny,nz), p(nx,ny,nz)
 real(dp) :: oprtf(nx,ny,nz), oprtp(nx,ny,nz)
 real(dp) :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)

 ierr = 1

 call auxp3d_5pt(nx, ny, nz, f, auxx, auxy, auxz)
 call oprt3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg3d(dx, dy, dz, nx, ny, nz, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    ' ==========================='
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    ' ---------------------------'
 end if

 do while(iter .le. maxiter)
    
    call auxp3d_5pt(nx, ny, nz, p, auxx, auxy, auxz)
    call oprt3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, p, oprtp)

    call intg3d(dx, dy, dz, nx, ny, nz, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp3d_5pt(nx, ny, nz, f, auxx, auxy, auxz)
!   call oprt3d_5pt(dx, dy, dz, nx, ny, nz, auxx, auxy, auxz, f, oprtf)
!   r = -oprtf+g

    call intg3d(dx, dy, dz, nx, ny, nz, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') ' ---------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1,1,1)
 f = f-offset

 end subroutine cjgd3d_5pt
!
 subroutine cjgd3d_3pt_mix(dx, dy, dz, nx, ny, nz, maxiter, tol, h, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, ny, nz, maxiter
 real(dp), intent(in)    :: dx, dy, dz, tol
 real(dp), intent(in)    :: h(nx,ny,nz), g(nx,ny,nz)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx,ny,nz)

 integer  :: iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx,ny,nz), p(nx,ny,nz)
 real(dp) :: oprtf(nx,ny,nz), oprtp(nx,ny,nz)
 real(dp) :: auxx(2,ny,nz), auxy(nx,2,nz), auxz(nx,ny,2)

 ierr = 1

 call auxp3d_3pt(nx, ny, nz, f, auxx, auxy, auxz)
 call oprt3d_3pt_mix(dx, dy, dz, nx, ny, nz, h, auxx, auxy, auxz, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg3d(dx, dy, dz, nx, ny, nz, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    ' ==========================='
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    ' ---------------------------'
 end if

 do while(iter .le. maxiter)
    
    call auxp3d_3pt(nx, ny, nz, p, auxx, auxy, auxz)
    call oprt3d_3pt_mix(dx, dy, dz, nx, ny, nz, h, auxx, auxy, auxz, p, oprtp)

    call intg3d(dx, dy, dz, nx, ny, nz, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp3d(nx, ny, nz, f, auxx, auxy, auxz)
!   call oprt3d_mix(dx, dy, dz, nx, ny, nz,  auxhx, auxhy, auxhz, h, auxx, auxy, auxz, f, oprtf)
!   r = -oprtf+g

    call intg3d(dx, dy, dz, nx, ny, nz, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') ' ---------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1,1,1)
 f = f-offset

 end subroutine cjgd3d_3pt_mix
!
 subroutine cjgd3d_5pt_mix(dx, dy, dz, nx, ny, nz, maxiter, tol, h, g, f, info, ierr)
 implicit none

 integer,  intent(in)    :: nx, ny, nz, maxiter
 real(dp), intent(in)    :: dx, dy, dz, tol
 real(dp), intent(in)    :: h(nx,ny,nz), g(nx,ny,nz)
 logical,  intent(in)    :: info
 integer,  intent(out)   :: ierr
 real(dp), intent(inout) :: f(nx,ny,nz)

 integer  :: iter
 real(dp) :: alpha, beta, rnorm, value, diff, offset
 real(dp) :: r(nx,ny,nz), p(nx,ny,nz)
 real(dp) :: oprtf(nx,ny,nz), oprtp(nx,ny,nz)
 real(dp) :: auxx(4,ny,nz), auxy(nx,4,nz), auxz(nx,ny,4)

 ierr = 1

 call auxp3d_5pt(nx, ny, nz, f, auxx, auxy, auxz)
 call oprt3d_5pt_mix(dx, dy, dz, nx, ny, nz, h, auxx, auxy, auxz, f, oprtf)

 r = -oprtf+g
 p = r
 iter = 0

 call intg3d(dx, dy, dz, nx, ny, nz, r*r, rnorm)

 if(info) then
   write(6,'(a)')    ' Conjugate-Gradient begins: '
   write(6,'(a)')    ' ==========================='
   write(6,'(a)')    ' iter        || Ax-b ||^2   '
   write(6,'(a)')    ' ---------------------------'
 end if

 do while(iter .le. maxiter)
    
    call auxp3d_5pt(nx, ny, nz, p, auxx, auxy, auxz)
    call oprt3d_5pt_mix(dx, dy, dz, nx, ny, nz, h, auxx, auxy, auxz, p, oprtp)

    call intg3d(dx, dy, dz, nx, ny, nz, p*oprtp, value)
    alpha = rnorm/value

    f = f+alpha*p
! The following iteration expression is not numerically stable
    r = r-alpha*oprtp
! The following gives feedback to the steepest descent direction
! but with sacrifice in the speed
!   call auxp3d_5pt(nx, ny, nz, f, auxx, auxy, auxz)
!   call oprt3d_5pt_mix(dx, dy, dz, nx, ny, nz, h, auxx, auxy, auxz, f, oprtf)
!   r = -oprtf+g

    call intg3d(dx, dy, dz, nx, ny, nz, r*r, diff)
    if(info) write(6,'(i5,f18.6)') iter, dabs(diff)
    if(dabs(diff) .le. tol) then
      ierr = 0
      exit
    endif

     beta = diff/rnorm
        p = r+beta*p
    rnorm = diff

    iter = iter+1
 end do

 if(info) write(6,'(a)') ' ---------------------------'
 if(ierr.gt.0) write(6,'(a)') 'WARNING: Required convergence not reached!'

 offset = f(1,1,1)
 f = f-offset

 end subroutine cjgd3d_5pt_mix
!
 subroutine intg1d(dx, nx, func, valu)
 implicit none

 integer,  intent(in)  :: nx
 real(dp), intent(in)  :: dx
 real(dp), intent(in)  :: func(nx)
 real(dp), intent(out) :: valu

 integer :: ix

 valu = 0.d0
 do ix = 1, nx-1
    valu = valu+0.5d0*(func(ix)+func(ix+1))
 end do
 valu = dx*(valu+0.5d0*(func(nx)+func(1)))

 end subroutine intg1d
!
 subroutine intg3d(dx, dy, dz, nx, ny, nz, func, valu)
 implicit none

 integer,  intent(in)  :: nx, ny, nz
 real(dp), intent(in)  :: dx, dy, dz
 real(dp), intent(in)  :: func(nx,ny,nz)
 real(dp), intent(out) :: valu

 integer :: iy, iz
 real(dp)  :: valuyz(ny,nz), valuz(nz)

 do iz = 1, nz
    do iy = 1, ny
       call intg1d(dx, nx, func(:,iy,iz), valuyz(iy,iz))
    end do
 end do

 do iz = 1, nz
    call intg1d(dy, ny, valuyz(:,iz), valuz(iz))
 end do

 call intg1d(dz, nz, valuz, valu)

 end subroutine intg3d

 end module m_numrecipe
