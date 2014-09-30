! Yi: module for TDDFT and near-field variables
 module m_variables
 use precision, only : dp, grid_p
 use files,     only : label_length
 use units,     only : pi
 implicit none
! Near-field basic parameters
 integer,  parameter ::           Nosc = 1 ! 9
 real(dp), parameter ::        eps0_1  = 2.0_dp * 4.0_dp * pi
 real(dp), parameter ::        eps0    = 1.0_dp / eps0_1
 real(dp), parameter ::        eps_vac = 1.0_dp
 real(dp), parameter :: alph_vac(Nosc) = 0.0_dp
 real(dp), parameter :: beta_vac(Nosc) = 0.0_dp
 real(dp), parameter :: omeg_vac(Nosc) = 0.0_dp
! Near-field metallic parameters
 real(dp), parameter ::  eps_Au = 1.0_dp
 real(dp), parameter :: alph_Au(Nosc) = &
                      (/0.011399455043448_dp/)!0.010877984258316_dp,  0.143004283607872_dp, &
!                       0.102613394133304_dp,  0.086968672396777_dp,  0.144370093236213_dp, &
!                       0.143875622429713_dp,  0.100026907030040_dp,  0.000000000000000_dp/)
 real(dp), parameter :: beta_Au(Nosc) = &
                      (/0.516567746244170_dp/)!-0.067800604864168_dp, -0.220886108538482_dp, &
!                       0.093025415787598_dp,  0.085145044297628_dp,  0.197898312405118_dp, &
!                       0.121802185561148_dp,  0.437787346695808_dp,  0.000000000000000_dp/)
 real(dp), parameter :: omeg_Au(Nosc) = &
                      (/0.017271294327114_dp/)!0.032418267762953_dp,  0.055881483393519_dp, &
!                       0.085336575783313_dp,  0.216560571334773_dp,  0.305871275766289_dp, &
!                       0.422419062083095_dp,  0.581503983914952_dp,  0.000000000000000_dp/)
! Near-field metallic parameters
 real(dp), parameter ::  eps_Ag = 1.0_dp
 real(dp), parameter :: alph_Ag(Nosc) = &
                      (/0.013193007507574_dp/)!0.018389361996630_dp,  0.155376144128204_dp, &
!                       0.119582302032446_dp,  0.133767541302428_dp,  0.077100082871564_dp, &
!                       0.073256103525346_dp,  0.190508498382360_dp,  0.203885252512602_dp/)
 real(dp), parameter :: beta_Ag(Nosc) = &
                      (/0.729276940986316_dp/)!-0.217702672012952_dp, -0.270426693820556_dp, &
!                       0.090376320168156_dp,  0.041331095373972_dp, -0.082975509730000_dp, &
!                       0.097615069063872_dp,  0.218350918182716_dp,  0.167571634884412_dp/)
 real(dp), parameter :: omeg_Ag(Nosc) = &
                      (/0.012465370881808_dp/)!0.026863756234086_dp,  0.046392347291260_dp, &
!                       0.086360912654040_dp,  0.152656694112716_dp,  0.295317571952284_dp, &
!                       0.311854767992418_dp,  0.389763336003722_dp,  0.528969777337128_dp/)
! Near-field metallic parameters
 real(dp), parameter ::  eps_drude = 1.0_dp
 real(dp), parameter :: alph_drude(Nosc) = &
                      (/0.036666666667_dp/) !0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
 real(dp), parameter :: beta_drude(Nosc) = &
                      (/0.640628967642_dp/) !0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
 real(dp), parameter :: omeg_drude(Nosc) = &
                      (/0.023333333333_dp/) !0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
! Sub-type, with geometrical parameters
 type :: sphere
      real(dp) :: xc, yc, zc, radius
 endtype sphere
!
 type :: shell
      real(dp) :: xc, yc, zc, iradius, oradius
 endtype shell
!
 type :: wire
      integer  :: orient
      real(dp) :: xc, yc, zc, length, radius
 endtype wire
!
 type :: tube
      integer  :: orient
      real(dp) :: xc, yc, zc, length, iradius, oradius
 endtype tube
!
 type :: rectangle
      real(dp) :: xc, yc, zc, lengthx, lengthy, lengthz
 endtype rectangle
!
 type :: sphere_smear
      real(dp) :: xc, yc, zc, radius, smear
 endtype sphere_smear
!
 type :: rectangle_smear
      integer  :: smrdir
      real(dp) :: xc, yc, zc, lengthx, lengthy, lengthz, smear
 endtype rectangle_smear
!
 type :: object
      character*20 :: geometry
      character*20 :: material
      real(dp)     :: eps, alph(Nosc), beta(Nosc), omeg(Nosc)
      type(sphere),          pointer :: sph=>null()
      type(shell),           pointer :: shl=>null()
      type(wire),            pointer :: wir=>null()
      type(tube),            pointer :: tub=>null()
      type(rectangle),       pointer :: rct=>null()
      type(sphere_smear),    pointer :: sph_smr=>null()
      type(rectangle_smear), pointer :: rct_smr=>null()
 endtype object
! Occupation threshold
 real(dp), parameter :: absocctol = 1.0e-12_dp
! Near-field parameter
 integer,  parameter ::      Nobj_default = 0
 integer,  parameter ::      nacc_default = 2
! TDDFT parameters
 integer,  parameter ::   maxnnfr_default = 5000
 integer,  parameter ::   maxnstb_default = 100
 integer,  parameter ::   maxntim_default = 16
 integer,  parameter :: tdproptyp_default = 1
 integer,  parameter ::    nuoind_default = 0
 integer,  parameter ::    spnind_default = 0
 integer,  parameter :: maxsavnfr_default = 5
 integer,  parameter ::  nkicknfr_default = 10
 integer,  parameter ::   homolvl_default = 0
 integer,  parameter ::   lumolvl_default = 0
 integer,  parameter :: nuowfkadd_default = 0
 integer,  parameter :: nuooffset_default = 0
 integer,  parameter :: fieldtype_default = 1
!integer,  parameter :: cgmaxiter_default = 9000
 real(dp), parameter ::  alphanfr_default = 1.0e-3_dp
 real(dp), parameter :: abstolstb_default = 1.0e-4_dp
 real(dp), parameter ::    tolnfr_default = 1.0e-4_dp
 real(dp), parameter ::tadiabatic_default = 1.0e12_dp
 real(dp), parameter ::      dtim_default = 1.0e-1_dp
 real(dp), parameter ::     sigma_default = 1.0e-2_dp
 real(dp), parameter ::      eext_default(3) = (/0.0_dp, 0.0_dp, 0.0_dp/)
 real(dp), parameter ::   homochg_default = 0.0_dp
 real(dp), parameter ::   lumochg_default = 0.0_dp
 real(dp), parameter ::    addchg_default = 0.0_dp
!real(dp), parameter ::     cgtol_default = 1.e-6_dp
 logical,  parameter :: simplempi_default = .false.
 logical,  parameter ::   infonfr_default = .true.
 logical,  parameter :: pseudonfr_default = .false.
 logical,  parameter :: pseudoatm_default = .false.
 logical,  parameter ::  nonscfdm_default = .false.
 logical,  parameter ::  nonscfwf_default = .false.
 logical,  parameter ::     tdhyb_default = .false.
 logical,  parameter ::  speconly_default = .false. 
 logical,  parameter ::  chkspows_default = .false.
 logical,  parameter :: chknormdm_default = .false.
 logical,  parameter ::  chkstabl_default = .false.
 logical,  parameter ::writetdpop_default = .false.
!logical,  parameter ::    cginfo_default = .true.
 complex(dp),  parameter :: i_dc =(0.0_dp,1.0_dp)
 character*60, parameter :: slabel_default = 'siesta'
! Near-field and TDDFT stabilization parameters
 integer  :: maxnnfr, maxnstb, maxntim
 integer  :: maxsavnfr, nkicknfr
 integer  :: tdproptyp
! Near-field differentiation accuracy
 integer  :: nacc = 2
! Near-field object number
 integer  :: Nobj
! Near-field grid points
 integer  :: nx, ny, nz, nxl, nyl, nzl
! Near-field and TDDFT stabilization parameters
 real(dp) :: alphanfr, alphakicknfr
 real(dp) :: tolnfr
 real(dp) :: abstolstb, tadiabatic, dtim, sigma
 real(dp) :: lbox(3), dbox(3)
! Near-field stabilization info
 logical  :: infonfr
! Poisson equation method
!character*6 :: poismethod
! Near-field conjugate-gradient parameters
!integer  :: cgmaxiter
!logical  :: cginfo
!real(dp) :: cgtol
! TDDFT variables
 integer  :: natot, nuotot, nuoloc, nuoocc, nuolococc
 integer  :: nuoind, spnind
 integer  :: fieldtype
 real(dp) :: estrength
 real(dp) :: eext(3)
 real(dp) :: dpdft0(3), dipdft(3), ddpdft(3)
 real(dp) :: dpnfr0(3), dipnfr(3), ddpnfr(3)
 logical  :: simplempi
 logical  :: chkspows, chknormdm, chkstabl
 logical  :: pseudonfr, pseudoatm, nonscfdm, nonscfwf, tdhyb, speconly, writetdpop
!
 character*3 :: stagenfr
!
 integer,      dimension(:),         pointer :: GridTableY=>null(), GridTableZ=>null()
 integer,      dimension(:),         pointer :: homolvl=>null(), lumolvl=>null()
 integer,      dimension(:),         pointer :: nuowfkadd=>null(), nuooffset=>null()
 real(dp),     dimension(:),         pointer :: xbox=>null(), ybox=>null(), zbox=>null()
 real(dp),     dimension(:),         pointer :: homochg=>null(), lumochg=>null()
 real(dp),     dimension(:),         pointer :: addchg=>null()
! TDDFT basic variales
 real(dp),     dimension(:,:),       pointer :: eig=>null(), occ=>null()
 real(dp),     dimension(:,:),       pointer :: wfkadd=>null()
 real(dp),     dimension(:,:),       pointer :: Saux=>null(),     Saux_inv=>null(),     Saux_nhf=>null(),     Saux_phf=>null()
 real(dp),     dimension(:,:),       pointer :: Saux_red=>null(), Saux_inv_red=>null(), Saux_nhf_red=>null(), Saux_phf_red=>null()
 real(dp),     dimension(:,:,:),     pointer :: Haux=>null(),     Daux=>null(),     wfk=>null()
 real(dp),     dimension(:,:,:),     pointer :: Haux_red=>null(), Daux_red=>null(), wfk_red=>null()
 real(dp),     dimension(:,:,:),     pointer :: Daux_old=>null()
 real(dp),     dimension(:,:,:),     pointer :: Daux_old_red=>null()
! Near-field basic variables
 real(dp),     dimension(:,:,:),     pointer :: eps=>null(), phi=>null()
 real(dp),     dimension(:,:,:),     pointer :: rhohyb=>null(), rhodft=>null(), rhonfr=>null(), rhonfr_old=>null()
 real(dp),     dimension(:,:,:,:),   pointer :: alph=>null(), beta=>null(), omeg=>null()
 real(dp),     dimension(:,:,:,:),   pointer :: efd=>null(), poltot=>null(), curtot=>null()
 real(dp),     dimension(:,:,:,:,:), pointer :: pol=>null(), cur=>null()
 real(grid_p), dimension(:),         pointer :: TRho=>null()
! TDDFT variables
 complex(dp),  dimension(:,:,:),     pointer :: zDaux=>null(), zDaux_old=>null()
 complex(dp),  dimension(:,:,:),     pointer :: zwfk=>null(),  zwfk_old=>null()
! Near-field objects
 type(object), dimension(:),         pointer :: nfobj=>null()
!
 end module m_variables
