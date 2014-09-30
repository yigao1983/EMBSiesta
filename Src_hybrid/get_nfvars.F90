 subroutine get_nfvars
 use parallel,    only : Node
 use m_variables, only : Nobj, nacc, nfobj
 use m_variables, only : pseudonfr
 use m_variables, only : infonfr, maxnnfr, maxsavnfr, nkicknfr
 use m_variables, only : alphanfr, alphakicknfr, tolnfr
#ifdef MPI
 use mpi_siesta
#endif 
 implicit none
! Local variables
 integer :: iobj, ierr
#ifdef MPI
 integer :: MPIerror
#endif
 external :: rednfr
!
 if(Node.eq.0) call rednfr()
#ifdef MPI
 call MPI_Bcast(Nobj,1,MPI_integer,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(nacc,1,MPI_integer,0,MPI_Comm_World,MPIerror)
!
 if(Node.ne.0) allocate(nfobj(Nobj))
!
 do iobj = 1,Nobj
    call MPI_Bcast(nfobj(iobj)%geometry,20,MPI_character,0,MPI_Comm_World,MPIerror)
    call MPI_Bcast(nfobj(iobj)%material,20,MPI_character,0,MPI_Comm_World,MPIerror)
!
    if(Node.ne.0) then
      nullify(nfobj(iobj)%sph)
      nullify(nfobj(iobj)%shl)
      nullify(nfobj(iobj)%wir)
      nullify(nfobj(iobj)%rct)
      nullify(nfobj(iobj)%sph_smr)
    endif
!
    select case(nfobj(iobj)%geometry)
    case('sph')
        if(Node.ne.0) allocate(nfobj(iobj)%sph)
        call MPI_Bcast(nfobj(iobj)%sph%xc,    1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph%yc,    1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph%zc,    1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph%radius,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    case('shl')  
        if(Node.ne.0) allocate(nfobj(iobj)%shl)
        call MPI_Bcast(nfobj(iobj)%shl%xc,     1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%shl%yc,     1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%shl%zc,     1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%shl%iradius,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%shl%oradius,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    case('wir')
        if(Node.ne.0) allocate(nfobj(iobj)%wir)
        call MPI_Bcast(nfobj(iobj)%wir%orient,1,MPI_integer,         0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%wir%xc,    1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%wir%yc,    1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%wir%zc,    1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%wir%length,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%wir%radius,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    case('rct')
        if(Node.ne.0) allocate(nfobj(iobj)%rct)
        call MPI_Bcast(nfobj(iobj)%rct%xc,     1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct%yc,     1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct%zc,     1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct%lengthx,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct%lengthy,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct%lengthz,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    case('sph_smr')
        if(Node.ne.0) allocate(nfobj(iobj)%sph_smr)
        call MPI_Bcast(nfobj(iobj)%sph_smr%xc, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph_smr%yc, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph_smr%zc, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph_smr%radius,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%sph_smr%smear, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    case('rct_smr')
        if(Node.ne.0) allocate(nfobj(iobj)%rct_smr)
        call MPI_Bcast(nfobj(iobj)%rct_smr%xc, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%yc, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%zc, 1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%lengthx,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%lengthy,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%lengthz,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%smrdir, 1,MPI_integer,         0,MPI_Comm_World,MPIerror)
        call MPI_Bcast(nfobj(iobj)%rct_smr%smear,  1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
    endselect
 enddo
!
 call MPI_Bcast(pseudonfr,1,MPI_logical,0,MPI_Comm_World,MPIerror)
!
 call MPI_Bcast(     maxnnfr,1,MPI_integer,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(   maxsavnfr,1,MPI_integer,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(    nkicknfr,1,MPI_integer,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(    alphanfr,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(alphakicknfr,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(      tolnfr,1,MPI_double_precision,0,MPI_Comm_World,MPIerror)
 call MPI_Bcast(     infonfr,1,MPI_logical,0,MPI_Comm_World,MPIerror)
#endif
!
 end subroutine get_nfvars
