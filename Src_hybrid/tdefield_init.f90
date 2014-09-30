 subroutine tdefield_init()
 use parallel,    only : Node
 use m_efield,    only : acting_efield, user_specified_field
 use m_variables, only : fieldtype, dtim, eext
 implicit none
!
 acting_efield = .true.
 if(fieldtype.eq.1) then
   user_specified_field = eext / dtim ! user_specified_field is positive
   if (Node.eq.0) &
    write(6,'(/a,3f12.6,a/)') 'tdefield_init: Induced external field =', eext, ' Ry/Bohr/e'
 else
   eext = user_specified_field
   if (Node.eq.0) &
    write(6,'(/a,3f12.6,a/)') 'tdefield_init: Begin to release field =', eext, ' Ry/Bohr/e'
 endif
!
 end subroutine tdefield_init
