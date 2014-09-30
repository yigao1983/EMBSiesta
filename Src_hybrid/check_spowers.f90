 subroutine check_spowers
 use precision,   only : dp
 use alloc,       only : re_alloc, de_alloc
 use parallel,    only : Node
 use m_variables, only : simplempi
 use m_variables, only : nuotot, nuoloc
 use m_variables, only : Saux,     Saux_inv,     Saux_nhf,     Saux_phf
 use m_variables, only : Saux_red, Saux_inv_red, Saux_nhf_red, Saux_phf_red
 implicit none
!
 integer :: io, jo, iu
 character(len=200) :: fname
 real(dp), dimension(nuotot,nuoloc) :: Maux
 real(dp), dimension(nuotot,nuotot) :: Maux_red
!
 external :: io_assign, io_close
 external :: collectdata_dp
!
 call timer( 'chkspows', 1 )
!
 if(.not.associated(Saux_red))     &
  call re_alloc(Saux_red,1,nuotot,1,nuotot,name='Saux_red',routine='check_spowers')
 if(.not.associated(Saux_inv_red)) &
  call re_alloc(Saux_inv_red,1,nuotot,1,nuotot,name='Saux_inv_red',routine='check_spowers')
 if(.not.associated(Saux_phf_red)) &
  call re_alloc(Saux_phf_red,1,nuotot,1,nuotot,name='Saux_phf_red',routine='check_spowers')
 if(.not.simplempi) then
 if(.not.associated(Saux_nhf_red)) &
  call re_alloc(Saux_nhf_red,1,nuotot,1,nuotot,name='Saux_nhf_red',routine='check_spowers')
 endif
!
 call collectdata_dp(nuotot, nuoloc, 1, Saux,     Saux_red)
 call collectdata_dp(nuotot, nuoloc, 1, Saux_inv, Saux_inv_red)
 call collectdata_dp(nuotot, nuoloc, 1, Saux_phf, Saux_phf_red)
 if(.not.simplempi) &
 call collectdata_dp(nuotot, nuoloc, 1, Saux_nhf, Saux_nhf_red)
!
 do io = 1,nuoloc
    do jo = 1,nuotot
       Maux(jo,io) = sum(Saux_inv_red(jo,:)*Saux(:,io))
    enddo
 enddo
!
 call collectdata_dp(nuotot, nuoloc, 1, Maux, Maux_red)
!
 if(Node.eq.0) then
   call io_assign(iu)
   fname = 'check_SinvS.ascii'
   open(unit=iu,file=fname)
   do io = 1,nuotot
      do jo = 1,nuotot
         write(iu,'(2x,3e15.7)') Maux_red(jo,io)
      enddo
   enddo
   call io_close(iu)
 endif
!
 do io = 1,nuoloc
    do jo = 1,nuotot
       Maux(jo,io) = sum(Saux_phf_red(jo,:)*Saux_phf(:,io))
    enddo
 enddo
!
 call collectdata_dp(nuotot, nuoloc, 1, Maux, Maux_red)
!
 if(Node.eq.0) then
   call io_assign(iu)
   fname = 'check_SphfSphf.ascii'
   open(unit=iu,file=fname)
   do io = 1,nuotot
      do jo = 1,nuotot
         write(iu,'(2x,3e15.7)') Maux_red(jo,io), Saux_red(jo,io), Maux_red(jo,io)-Saux_red(jo,io)
      enddo
   enddo
   call io_close(iu)
 endif
!
 do io = 1,nuoloc
    do jo = 1,nuotot
       Maux(jo,io) = sum(Saux_nhf_red(jo,:)*Saux_nhf(:,io))
    enddo
 enddo
!
 call collectdata_dp(nuotot, nuoloc, 1, Maux, Maux_red)
!
 if(Node.eq.0) then
   call io_assign(iu)
   fname = 'check_SnhfSnhf.ascii'
   open(unit=iu,file=fname)
   do io = 1,nuotot
      do jo = 1,nuotot
         write(iu,'(2x,3e15.7)') Maux_red(jo,io), Saux_inv_red(jo,io), Maux_red(jo,io)-Saux_inv_red(jo,io)
      enddo
   enddo
   call io_close(iu)
 endif
!
 call de_alloc(Saux_red,    name='Saux_red')
 call de_alloc(Saux_inv_red,name='Saux_inv_red')
 call de_alloc(Saux_phf_red,name='Saux_phf_red')
 if(.not.simplempi) &
 call de_alloc(Saux_nhf_red,name='Saux_nhf_red')
!
 call timer( 'chkspows', 2 )
!
 end subroutine check_spowers
