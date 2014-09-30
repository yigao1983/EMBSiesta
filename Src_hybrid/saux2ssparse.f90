 subroutine saux2ssparse(Saux)
 use precision,       only : dp
 use alloc,           only : re_alloc
 use sparse_matrices, only : maxnh, listh, listhptr, numh, S
 use m_variables,     only : nuotot, nuoloc
 implicit none
 real(dp), intent(in) :: Saux(nuotot,nuoloc)
! Local variables
 integer  :: io, jo, j, ind
!
 if(.not.associated(S)) call re_alloc(S,1,maxnh,name='maxnh',routine='saux2ssparse')
!
 S = 0.0_dp
 do io = 1,nuoloc
    do j = 1,numh(io)
       ind = listhptr(io) + j
        jo = listh(ind)
       S(ind) = Saux(jo,io)
    enddo
 enddo
!
 end subroutine saux2ssparse
