 subroutine writetdwave(itim)
 use files,        only : slabel
 use parallel,     only : Node
 use atmfuncs,     only : symfio, cnfigfio, labelfis
 use atomlist,     only : iaorb, iphorb
 use siesta_geom,  only : isa
 use m_spin,       only : nspin
 use m_variables,  only : nuotot
 implicit none
 integer, intent(in) :: itim
! Local variables
 integer :: iu, j
 character(len=200) :: wfs_filename
!
 wfs_filename = trim(slabel)//".WFSX"
!
 if (Node.eq.0) then
    call io_assign( iu )
    open(iu, file=wfs_filename,form="unformatted",status='unknown')
    rewind (iu)
    write(iu) nk, gamma
    write(iu) nspin
    write(iu) nuotot
    write(iu) (iaorb(j),labelfis(isa(iaorb(j))),            &
              iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)), &
              symfio(isa(iaorb(j)),iphorb(j)), j=1,nuotot)
    call io_close(iu)
 endif
!
 do ispin = 1,nspin
!
    if (Node .eq. 0) then
       call io_assign( iu )
       open(iu,file=wfs_filename,form="unformatted",position='append',status='old')
!
       write(iu) ik, k(1:3), kpoint_weight
       write(iu) ispin
       write(iu) nwflist(ik)
!
       write(iu) indwf
       write(iu) eig(indwf)/eV
       write(iu) (aux(1:,j), j=1,ntot)
!
       close (iu)
       call io_close(iu)
!
    endif
!
 enddo
!
 end subroutine writetdwave
