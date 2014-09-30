 subroutine rednfr
 use fdf
 use units,       only : Ang
 use m_variables, only : Nobj_default, nacc_default
 use m_variables, only : infonfr_default, maxnnfr_default
 use m_variables, only : maxsavnfr_default, nkicknfr_default
 use m_variables, only : alphanfr_default, tolnfr_default
 use m_variables, only : Nobj, nacc, nfobj, pseudonfr
 use m_variables, only : infonfr, maxnnfr
 use m_variables, only : maxsavnfr, nkicknfr
 use m_variables, only : alphanfr, alphakicknfr, tolnfr
 implicit none
! Local variables
 integer :: iu, iobj, iiobj, ierr
!
 Nobj = fdf_integer('NumberOfNFObjects', Nobj_default)
 nacc = fdf_integer('NF.DiffAccuracy',   nacc_default)
 if(nacc.ne.2 .and. nacc.ne.4) nacc = nacc_default
!
 if(Nobj.gt.0) allocate(nfobj(Nobj))
!
 if(fdf_block('NFObjects', iu)) then
   do iobj = 1,Nobj
      read(iu,*) iiobj, nfobj(iiobj)%geometry, nfobj(iiobj)%material
      nullify(nfobj(iiobj)%sph)
      nullify(nfobj(iiobj)%shl)
      nullify(nfobj(iiobj)%wir)
      nullify(nfobj(iiobj)%rct)
      nullify(nfobj(iiobj)%sph_smr)
!
      select case(nfobj(iiobj)%geometry)
      case('sph')
          allocate(nfobj(iiobj)%sph, stat=ierr)
          read(iu,*) nfobj(iiobj)%sph%xc, nfobj(iiobj)%sph%yc, nfobj(iiobj)%sph%zc, &
                     nfobj(iiobj)%sph%radius
          nfobj(iiobj)%sph%xc = nfobj(iiobj)%sph%xc*Ang
          nfobj(iiobj)%sph%yc = nfobj(iiobj)%sph%yc*Ang
          nfobj(iiobj)%sph%zc = nfobj(iiobj)%sph%zc*Ang
          nfobj(iiobj)%sph%radius = nfobj(iiobj)%sph%radius*Ang
      case('shl')
          allocate(nfobj(iiobj)%shl, stat=ierr)
          read(iu,*) nfobj(iiobj)%shl%xc, nfobj(iiobj)%shl%yc, nfobj(iiobj)%shl%zc, &
                     nfobj(iiobj)%shl%iradius, nfobj(iiobj)%shl%oradius
          nfobj(iiobj)%shl%xc = nfobj(iiobj)%shl%xc*Ang
          nfobj(iiobj)%shl%yc = nfobj(iiobj)%shl%yc*Ang
          nfobj(iiobj)%shl%zc = nfobj(iiobj)%shl%zc*Ang
          nfobj(iiobj)%shl%iradius = nfobj(iiobj)%shl%iradius*Ang
          nfobj(iiobj)%shl%oradius = nfobj(iiobj)%shl%oradius*Ang
      case('wir')
          allocate(nfobj(iiobj)%wir, stat=ierr)
          read(iu,*) nfobj(iiobj)%wir%xc, nfobj(iiobj)%wir%yc, nfobj(iiobj)%wir%zc, &
                     nfobj(iiobj)%wir%orient, nfobj(iiobj)%wir%length, nfobj(iiobj)%wir%radius
          nfobj(iiobj)%wir%xc = nfobj(iiobj)%wir%xc*Ang
          nfobj(iiobj)%wir%yc = nfobj(iiobj)%wir%yc*Ang
          nfobj(iiobj)%wir%zc = nfobj(iiobj)%wir%zc*Ang
          nfobj(iiobj)%wir%length = nfobj(iiobj)%wir%length*Ang
          nfobj(iiobj)%wir%radius = nfobj(iiobj)%wir%radius*Ang
      case('rct')
          allocate(nfobj(iiobj)%rct, stat=ierr)
          read(iu,*) nfobj(iiobj)%rct%xc, nfobj(iiobj)%rct%yc, nfobj(iiobj)%rct%zc, &
                     nfobj(iiobj)%rct%lengthx, nfobj(iiobj)%rct%lengthy, nfobj(iiobj)%rct%lengthz
          nfobj(iiobj)%rct%xc = nfobj(iiobj)%rct%xc*Ang
          nfobj(iiobj)%rct%yc = nfobj(iiobj)%rct%yc*Ang
          nfobj(iiobj)%rct%zc = nfobj(iiobj)%rct%zc*Ang
          nfobj(iiobj)%rct%lengthx = nfobj(iiobj)%rct%lengthx*Ang
          nfobj(iiobj)%rct%lengthy = nfobj(iiobj)%rct%lengthy*Ang
          nfobj(iiobj)%rct%lengthz = nfobj(iiobj)%rct%lengthz*Ang
      case('sph_smr')
          allocate(nfobj(iiobj)%sph_smr, stat=ierr)
          read(iu,*) nfobj(iiobj)%sph_smr%xc, nfobj(iiobj)%sph_smr%yc, nfobj(iiobj)%sph_smr%zc, &
                     nfobj(iiobj)%sph_smr%radius, nfobj(iiobj)%sph_smr%smear
          nfobj(iiobj)%sph_smr%xc = nfobj(iiobj)%sph_smr%xc*Ang
          nfobj(iiobj)%sph_smr%yc = nfobj(iiobj)%sph_smr%yc*Ang
          nfobj(iiobj)%sph_smr%zc = nfobj(iiobj)%sph_smr%zc*Ang
          nfobj(iiobj)%sph_smr%radius = nfobj(iiobj)%sph_smr%radius*Ang
          nfobj(iiobj)%sph_smr%smear  = nfobj(iiobj)%sph_smr%smear *Ang
      case('rct_smr')
          allocate(nfobj(iiobj)%rct_smr, stat=ierr)
          read(iu,*) nfobj(iiobj)%rct_smr%xc, nfobj(iiobj)%rct_smr%yc, nfobj(iiobj)%rct_smr%zc, &
                     nfobj(iiobj)%rct_smr%lengthx, nfobj(iiobj)%rct_smr%lengthy, nfobj(iiobj)%rct_smr%lengthz, &
                     nfobj(iiobj)%rct_smr%smrdir, nfobj(iiobj)%rct_smr%smear
          nfobj(iiobj)%rct_smr%xc = nfobj(iiobj)%rct_smr%xc*Ang
          nfobj(iiobj)%rct_smr%yc = nfobj(iiobj)%rct_smr%yc*Ang
          nfobj(iiobj)%rct_smr%zc = nfobj(iiobj)%rct_smr%zc*Ang
          nfobj(iiobj)%rct_smr%lengthx = nfobj(iiobj)%rct_smr%lengthx*Ang
          nfobj(iiobj)%rct_smr%lengthy = nfobj(iiobj)%rct_smr%lengthy*Ang
          nfobj(iiobj)%rct_smr%lengthz = nfobj(iiobj)%rct_smr%lengthz*Ang
          nfobj(iiobj)%rct_smr%smear   = nfobj(iiobj)%rct_smr%smear  *Ang
      endselect
   enddo
 endif
!
 if(Nobj.le.0) then
   pseudonfr = .true.
 else
   pseudonfr = .false.
 endif
!
      maxnnfr = fdf_integer('MaxNFIterations',      maxnnfr_default)
    maxsavnfr = fdf_integer('NF.NumberPulay',     maxsavnfr_default)
     nkicknfr = fdf_integer('NF.NumberKick',       nkicknfr_default)
     alphanfr = fdf_double( 'NF.MixingWeight',     alphanfr_default)
 alphakicknfr = fdf_double( 'NF.KickMixingWeight', alphanfr_default)
       tolnfr = fdf_double( 'NF.Tolerance',          tolnfr_default)
      infonfr = fdf_boolean('NF.Information',       infonfr_default)
!
 end subroutine rednfr
