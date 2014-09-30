 subroutine magnus(itim, dt, zDaux_old, zDaux_new)
 use precision,       only : dp
 use alloc,           only : re_alloc
 use m_spin,          only : nspin
 use m_variables,     only : nuotot, nuoloc
 use m_variables,     only : Saux, Haux
 implicit none
 integer,     intent(in)  :: itim
 real(dp),    intent(in)  :: dt
 complex(dp), intent(in)  :: zDaux_old(nuotot,nuoloc,nspin)
 complex(dp), intent(out) :: zDaux_new(nuotot,nuoloc,nspin)
! Local variables
 logical :: newS = .false.
 logical, save :: FirstCall = .true.
 integer :: ispin, io, jo
!real(dp), dimension(nuotot,nuoloc,nspin) :: Haux_3, Haux_5
!real(dp), dimension(nuotot,nuoloc,nspin) :: Daux_5
 real(dp), dimension(:,:,:), pointer, save :: Haux_3=>null(), Haux_5=>null()
 real(dp), dimension(:,:,:), pointer, save :: Daux_5=>null()
! Saved variables
 real(dp), dimension(:,:,:), pointer, save :: Haux_1a=>null(), Haux_1b=>null()
!
 call timer( 'magnus', 1 )
!
 if(FirstCall) then
   if(.not.associated(Haux_3))  &
    call re_alloc(Haux_3,  1, nuotot, 1, nuoloc, 1, nspin, name='Haux_1a', routine='magnus')
   if(.not.associated(Haux_5))  &
    call re_alloc(Haux_5,  1, nuotot, 1, nuoloc, 1, nspin, name='Haux_1b', routine='magnus')
   if(.not.associated(Daux_5))  &
    call re_alloc(Daux_5,  1, nuotot, 1, nuoloc, 1, nspin, name='Haux_1b', routine='magnus')
   if(.not.associated(Haux_1a)) &
    call re_alloc(Haux_1a, 1, nuotot, 1, nuoloc, 1, nspin, name='Haux_1a', routine='magnus')
   if(.not.associated(Haux_1b)) &
    call re_alloc(Haux_1b, 1, nuotot, 1, nuoloc, 1, nspin, name='Haux_1b', routine='magnus')
   FirstCall = .false.
 endif
!
 if(itim.eq.1) then
!  do ispin = 1,nspin
!     do io = 1,nuoloc
!        do jo = 1,nuotot
!           Daux_5(jo,io,ispin) = dreal(zDaux_old(jo,io,ispin))
!        enddo
!     enddo
!  enddo
   Daux_5(:,:,:) = dreal(zDaux_old(:,:,:))
   call sauxdaux2haux(itim, newS, Saux, Daux_5, Haux_1a)
!  do ispin = 1,nspin
!     do io = 1,nuoloc
!        do jo = 1,nuotot
!           Haux_1b(jo,io,ispin) = Haux_1a(jo,io,ispin)
!        enddo
!     enddo
!  enddo
   Haux_1b(:,:,:) = Haux_1a(:,:,:)
 endif
!
!do ispin = 1,nspin
!   do io = 1,nuoloc
!      do jo = 1,nuotot
!         Haux_3(jo,io,ispin) = -0.75_dp * Haux_1a(jo,io,ispin) + 1.75_dp * Haux_1b(jo,io,ispin)
!      enddo
!   enddo
!enddo
 Haux_3(:,:,:) = -0.75_dp * Haux_1a(:,:,:) + 1.75_dp * Haux_1b(:,:,:)
!
 call propagat_zdaux(0.5_dp*dt, Haux_3, zDaux_old, zDaux_new)
!
!do ispin = 1,nspin
!   do io = 1,nuoloc
!      do jo = 1,nuotot
!         Daux_5(jo,io,ispin) = dreal(zDaux_new(jo,io,ispin))
!      enddo
!   enddo
!enddo
 Daux_5(:,:,:) = dreal(zDaux_new(:,:,:))
!
 call sauxdaux2haux(itim, newS, Saux, Daux_5, Haux_5)
!
 call propagat_zdaux(dt, Haux_5, zDaux_old, zDaux_new)
!
!do ispin = 1,nspin
!   do io = 1,nuoloc
!      do jo = 1,nuotot
!         Haux_1a(jo,io,ispin) = Haux_1b(jo,io,ispin)
!         Haux_1b(jo,io,ispin) =  Haux_5(jo,io,ispin)
!      enddo
!   enddo
!enddo
 Haux_1a(:,:,:) = Haux_1b(:,:,:)
 Haux_1b(:,:,:) =  Haux_5(:,:,:)
!
 call timer( 'magnus', 2 )
!
 end subroutine magnus
