!
! pGFNFF library:
! ---------------
!
! Copyright (C) 2020-2022 Julian Gale
!
! pGFNFF is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! pGFNFF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with pGFNFF.  If not, see <https://www.gnu.org/licenses/>.
!
program main
!
!  Main program to call subroutines from pGFNFF library in order to generate parameters
!  for a material of interest. 
!
!  NB: This example is currently setup for a periodic cubic cell of NaCl but can be
!      edited to run other systems of interest.
!
!  Julian Gale, Curtin University, October 2021
!
  use m_pgfnff_types
  use m_nbr
  use m_io
  implicit none
!
!  User specified
!
  integer(i4)                   :: ndim           ! Number of periodic dimensions
  integer(i4)                   :: numat          ! Number of atoms
  integer(i4)                   :: nnobo          ! Number of bond exclusions
  integer(i4)                   :: maxele         ! Maximum element number supported by GFN-FF
  integer(i4)                   :: maxelein       ! Number of element values for dispersion
  integer(i4)                   :: maxrefin       ! Number of ref values for dispersion
  integer(i4),      allocatable :: nat(:)         ! Atomic number
  integer(i4),      allocatable :: nftype(:)      ! Atom type number
  integer(i4),      allocatable :: nobond(:)      ! Convolution of atomic numbers for bond exclusion (=1000*n1 + n2;n1<n2)
  integer(i4),      allocatable :: nobotyp(:)     ! Convolution of atom types for bond exclusion (=1000*nt1 + nt2)
  logical                       :: ldebug         ! If true then debug printing is output
  logical                       :: lverbose       ! If true then verbose printing is output
  real(dp)                      :: kv(3,3)        ! Reciprocal space lattice vectors (1/Angstroms)
  real(dp)                      :: rv(3,3)        ! Real space lattice vectors (Angstroms)
  real(dp),         allocatable :: q(:)           ! Charges
  real(dp),         allocatable :: x(:)           ! Cartesian x coordinates (Angstroms)
  real(dp),         allocatable :: y(:)           ! Cartesian y coordinates (Angstroms)
  real(dp),         allocatable :: z(:)           ! Cartesian z coordinates (Angstroms)
!
!  Returned quantities
!
  integer(i4)                   :: nfrag        ! Number of fragments
  integer(i4),      allocatable :: nfraglist(:) ! Atomic number
  real(dp),         allocatable :: qfrag(:)     ! Fragment charges
!
!  General variables
!
  integer(i4)                   :: i            ! Looping index
  integer(i4)                   :: ierror       ! Error flag
  integer(i4)                   :: numat2       ! Lower half triangular array size for numat
  logical                       :: lbondsok     ! If true then all atoms have neighbours where they should 
  real(dp)                      :: pi           ! Constant pi
  real(dp)                      :: cut2_nbr     ! Maximum neighbour list cutoff squared
!
!  pGFNFF control parameters 
!
  integer(i4)                   :: pi_change = 0             ! Flag for type of pi electron change : 0 => use charge; + or - for direction
  real(dp)                      :: pi_temp1 = 4000.0_dp      ! Temperature for Huckel pi calculation on first pass
  real(dp)                      :: pi_temp2 =  300.0_dp      ! Temperature for Huckel pi calculation on second pass
  real(dp)                      :: ks = 0.04_dp              ! Maximum spacing between points in reciprocal space
  integer(i4)                   :: maxtoposhell  = 4         ! Maximum number of shells to search for topology
  integer(i4)                   :: maxtoposhell1 = 2         ! Maximum number of shells to search for topology in the first loop
  logical                       :: lnewtopo = .true.         ! This flag controls whether the topological charges use 1d+12 or 1+12.
  logical                       :: lxtbtopo = .false.        ! If true use XTB topology approach
  real(dp)                      :: tdist_thr=12.0_dp         ! R threshold for covalent distance estimated used in apprx EEQ
  real(dp)                      :: atm_alpha1 = 1.0_dp       ! Damping of three-body dispersion parameter
  real(dp)                      :: q_trap = -2.000_dp        ! Largest negative charge allowed for pgfnff_radij
  logical                       :: lhighcn_trap = .false.    ! Controls whether to trap high coordination numbers as per XTB
  logical                       :: lfragment_bond = .false.  ! Controls whether to enforce use of bonds to set fragments
  real(dp)                      :: max_accuracy = 31.622_dp  ! Maximum possible overall accuracy value for GFNFF
  real(dp)                      :: max_acc_disp = 31.622_dp  ! Maximum possible accuracy value for GFNFF for dispersion
  real(dp)                      :: max_acc_rep  = 10000.0_dp ! Maximum possible accuracy value for GFNFF for repulsion
  real(dp)                      :: max_acc_cn   = 100.0_dp   ! Maximum possible accuracy value for GFNFF for coordination numbers
  real(dp)                      :: max_acc_hb1  = 10000.0_dp ! Maximum possible accuracy value for GFNFF for hydrogen bonding 1
  real(dp)                      :: max_acc_hb2  = 1.0d8      ! Maximum possible accuracy value for GFNFF for hydrogen bonding 2
  real(dp)                      :: accuracy = 0.1_dp         ! Overall accuracy value for GFNFF
  real(dp)                      :: accuracy_disp = 0.1_dp    ! Accuracy value for GFNFF for dispersion
  real(dp)                      :: accuracy_rep = 0.1_dp     ! Accuracy value for GFNFF for repulsion
  real(dp)                      :: accuracy_cn = 0.1_dp      ! Accuracy value for GFNFF for coordination numbers
  real(dp)                      :: accuracy_hb1 = 0.1_dp     ! Accuracy value for GFNFF for hydrogen bonding 1
  real(dp)                      :: accuracy_hb2 = 0.1_dp     ! Accuracy value for GFNFF for hydrogen bonding 2
  real(dp)                      :: taper = 0.95_dp           ! Fractional taper range for GFNFF
  real(dp)                      :: wolf_eta = 0.2_dp         ! Eta value for Wolf sum in topology
  real(dp)                      :: scale_r_topo = 1.175_dp   ! Scaling factor for topological distances
  logical                       :: ltopowolf = .true.        ! Use Wolf sum during topology
!
!  pGFNFF force field parameters generated by library
!
  integer(i4),      allocatable :: disp_nref(:)              ! Dispersion number of reference states
  integer(i4)                   :: nABatoms                  ! Number of potential AB atoms for hydrogen bond
  integer(i4)                   :: nHatoms                   ! Number of potential H atoms for hydrogen bond
  integer(i4),      allocatable :: nABatomptr(:)             ! Potential hydrogen bond AB atom pointer
  integer(i4),      allocatable :: nHatomptr(:)              ! Potential hydrogen bond H atom pointer
  integer(i4),      allocatable :: nHatomrptr(:)             ! Potential hydrogen bond H atom reverse pointer
  integer(i4)                   :: nangles                   ! Number of angle bends
  integer(i4),      allocatable :: nangleatomptr(:)          ! Angle bend atom pointers
  integer(i4),      allocatable :: ntorsionatomptr(:)        ! Torsion atom pointers
  integer(i4)                   :: ntorsions                 ! Number of torsions
  integer(i4)                   :: nxbABatoms                ! Number of potential AB atoms for halogen bond
  integer(i4),      allocatable :: nxbABatomptr(:,:)         ! Potential halogen bond AB atom pointer
  real(dp)                      :: cn_kn
  real(dp)                      :: hb_kn
  real(dp)                      :: gfnff_angle_damp          ! Angle damping factor
  real(dp)                      :: gfnff_bend_hb             ! Hydrogen bond parameter
  real(dp)                      :: gfnff_tors_hb             ! Hydrogen bond parameter
  real(dp)                      :: gfnff_hbabmix             ! Hydrogen bond parameter
  real(dp)                      :: gfnff_bondscale           ! Bond scale factor
  real(dp)                      :: gfnff_dispthr             ! Dispersion threshold
  real(dp)                      :: gfnff_cnmax               ! Maximum coordination number parameter
  real(dp)                      :: gfnff_hb_alp              ! Hydrogen bonding cutoff scaling
  real(dp)                      :: gfnff_hb_a_cut            ! Hydrogen bonding cutoff
  real(dp)                      :: gfnff_hb_long_cut         ! Hydrogen bonding cutoff
  real(dp)                      :: gfnff_hb_nb_cut           ! Hydrogen bonding cutoff
  real(dp)                      :: gfnff_hb_s_cut            ! Hydrogen bonding cutoff
  real(dp)                      :: gfnff_hb_scale_gen        ! Hydrogen bonding scale factor
  real(dp)                      :: gfnff_hb_scale_coh        ! Hydrogen bonding scale factor - COH
  real(dp)                      :: gfnff_hbthr1              ! Hydrogen bonding threshold 1
  real(dp)                      :: gfnff_hbthr2              ! Hydrogen bonding threshold 2
  real(dp)                      :: gfnff_repscale_b          ! Repulsion scale factor
  real(dp)                      :: gfnff_repscale_n          ! Repulsion scale factor
  real(dp)                      :: gfnff_repscale_13         ! Repulsion scale factor for 1-3 interactions
  real(dp)                      :: gfnff_repscale_14         ! Repulsion scale factor for 1-4 interactions
  real(dp)                      :: gfnff_repthr              ! Repulsion threshold
  real(dp)                      :: gfnff_torsion_damp        ! Torsion damping factor
  real(dp)                      :: gfnff_xb_a_cut            ! Halogen bonding cutoff
  real(dp)                      :: gfnff_xb_s_cut            ! Halogen bonding cutoff
  real(dp)                      :: gfnff_xb_l_cut            ! Halogen bonding cutoff
  real(dp),         allocatable :: gfnff_xb_scale(:)         ! Halogen bond scale factors for elements
  real(dp),         allocatable :: gfnff_ABhbq(:)            ! Charge-scaling parameters for hydrogen bonds
  real(dp),         allocatable :: gfnff_ABxbq(:,:)          ! Charge-scaling parameters for halogen bonds
  real(dp),         allocatable :: gfnff_rad(:)              ! General radii
  real(dp),         allocatable :: gfnff_rad_cn(:,:)         ! Parameters than control dependence of radii on CN
  real(dp),         allocatable :: gfnff_rcov(:)             ! Covalent radii
  real(dp),         allocatable :: gfnff_alp(:)              ! EEQ parameters
  real(dp),         allocatable :: gfnff_chi(:)              ! EEQ parameters
  real(dp),         allocatable :: gfnff_cnf(:)              ! EEQ parameters
  real(dp),         allocatable :: gfnff_gam(:)              ! EEQ parameters
  real(dp),         allocatable :: gfnff_hb_acid(:)          ! Hydrogen bonding parameters
  real(dp),         allocatable :: gfnff_hb_base(:)          ! Hydrogen bonding parameters
  real(dp),         allocatable :: gfnff_repulsion_a(:)      ! Repulsion parameters
  real(dp),         allocatable :: gfnff_repulsion_p(:)      ! Repulsion parameters
  real(dp),         allocatable :: gfnff_repulsion_z(:)      ! Repulsion parameters
  real(dp),         allocatable :: disp_cn(:,:)              ! Dispersion coordination parameters
  real(dp),         allocatable :: disp_c6(:,:,:,:)          ! Dispersion C6 parameters
  real(dp),         allocatable :: disp_c9(:)                ! Dispersion terms for C9
  real(dp),         allocatable :: disp_r0(:)                ! Dispersion radius parameters
  real(dp),         allocatable :: disp_sqrtZr4r2(:)         ! Dispersion parameters
  real(dp),         allocatable :: disp_zeta(:)              ! Dispersion atom-specific zeta parameters
  real(dp),         allocatable :: par_angle(:,:)            ! Angle bend parameters
  real(dp),         allocatable :: par_bond(:,:,:)           ! Bond parameters
  real(dp),         allocatable :: par_torsion(:,:)          ! Torsion parameters
!
!  Pass values to pGFNFF library routines
!
  call pgfnff_init_param_r('tdist_thr',tdist_thr,ierror)
  call pgfnff_init_param_r('atm_alpha1',atm_alpha1,ierror)
  call pgfnff_init_param_r('q_trap',q_trap,ierror)
  call pgfnff_init_param_r('pi_temp1',pi_temp1,ierror)
  call pgfnff_init_param_r('pi_temp2',pi_temp2,ierror)
  call pgfnff_init_param_r('ks',ks,ierror)
  call pgfnff_init_param_r('max_accuracy',max_accuracy,ierror)
  call pgfnff_init_param_r('max_acc_disp',max_acc_disp,ierror)
  call pgfnff_init_param_r('max_acc_rep',max_acc_rep,ierror)
  call pgfnff_init_param_r('max_acc_cn',max_acc_cn,ierror)
  call pgfnff_init_param_r('max_acc_hb1',max_acc_hb1,ierror)
  call pgfnff_init_param_r('max_acc_hb2',max_acc_hb2,ierror)
  call pgfnff_init_param_r('accuracy_overall',accuracy,ierror)
  call pgfnff_init_param_r('accuracy_disp',accuracy_disp,ierror)
  call pgfnff_init_param_r('accuracy_rep',accuracy_rep,ierror)
  call pgfnff_init_param_r('accuracy_cn',accuracy_cn,ierror)
  call pgfnff_init_param_r('accuracy_hb1',accuracy_hb1,ierror)
  call pgfnff_init_param_r('accuracy_hb2',accuracy_hb2,ierror)
  call pgfnff_init_param_r('wolf_eta',wolf_eta,ierror)
  call pgfnff_init_param_r('taper',taper,ierror)
  call pgfnff_init_param_r('rtopo',scale_r_topo,ierror)
!
  call pgfnff_init_param_i('maxtoposhell1',maxtoposhell1,ierror)
  call pgfnff_init_param_i('maxtoposhell',maxtoposhell,ierror)
  call pgfnff_init_param_i('pi_change',pi_change,ierror)
!
  call pgfnff_init_param_l('newtopo',lnewtopo,ierror)
  call pgfnff_init_param_l('xtbtopo',lxtbtopo,ierror)
  call pgfnff_init_param_l('highcn_trap',lhighcn_trap,ierror)
  call pgfnff_init_param_l('fragment_bond',lfragment_bond,ierror)
  call pgfnff_init_param_l('topowolf',ltopowolf,ierror)
!
!  Output header
!
  write(ioout,'(''********************************************************************************'')')
  write(ioout,'(''*  pGFN-FF Parameter Calculation                                               *'')')
  write(ioout,'(''********************************************************************************'',/)')
!
!  Print levels
!
  ldebug = .true.
  lverbose = .true.
!
!  Set basic system information
!
  numat = 8
  ndim  = 3
!
!  Lattice vectors in Angstroms
!
  rv(1,1) = 5.60
  rv(2,1) = 0.00
  rv(3,1) = 0.00
  rv(1,2) = 0.00
  rv(2,2) = 5.60
  rv(3,2) = 0.00
  rv(1,3) = 0.00
  rv(2,3) = 0.00
  rv(3,3) = 5.60
!
!  Set reciprocal space lattice vectors
!  NB: Below is for special case of a cubic cell!
!
  pi = 4.0_dp*atan(1.0_dp)
  kv(1,1) = 2.0_dp*pi/rv(1,1)
  kv(2,1) = 0.00
  kv(3,1) = 0.00
  kv(1,2) = 0.00
  kv(2,2) = 2.0_dp*pi/rv(2,2)
  kv(3,2) = 0.00
  kv(1,3) = 0.00
  kv(2,3) = 0.00
  kv(3,3) = 2.0_dp*pi/rv(3,3)
!
!  Create arrays for atoms
!
  allocate(nat(numat))
  allocate(nftype(numat))
  allocate(q(numat))
  allocate(x(numat))
  allocate(y(numat))
  allocate(z(numat))
!
!  Atoms
!
  nat(1) = 11
  nftype(1) = 0
  x(1) = 0.0_dp
  y(1) = 0.0_dp
  z(1) = 0.0_dp
!
  nat(2) = 11
  nftype(2) = 0
  x(2) = 0.0_dp
  y(2) = 2.8_dp
  z(2) = 2.8_dp
!
  nat(3) = 11
  nftype(3) = 0
  x(3) = 2.8_dp
  y(3) = 0.0_dp
  z(3) = 2.8_dp
!
  nat(4) = 11
  nftype(4) = 0
  x(4) = 2.8_dp
  y(4) = 2.8_dp
  z(4) = 0.0_dp
!
  nat(5) = 17
  nftype(5) = 0
  x(5) = 2.8_dp
  y(5) = 2.8_dp
  z(5) = 2.8_dp
!
  nat(6) = 17
  nftype(6) = 0
  x(6) = 2.8_dp
  y(6) = 0.0_dp
  z(6) = 0.0_dp
!
  nat(7) = 17
  nftype(7) = 0
  x(7) = 0.0_dp
  y(7) = 2.8_dp
  z(7) = 0.0_dp
!
  nat(8) = 17
  nftype(8) = 0
  x(8) = 0.0_dp
  y(8) = 0.0_dp
  z(8) = 2.8_dp
!
!  Exclude bonds between atoms
!
  nnobo = 2
  allocate(nobond(nnobo))
  allocate(nobotyp(nnobo))
  nobond(1) = 1000*11 + 11   ! Exclude Na-Na bonds
  nobond(2) = 1000*17 + 17   ! Exclude Cl-Cl bonds
!
!  Create arrays for fragments
!
  allocate(nfraglist(numat))
  allocate(qfrag(numat))
!********************************************
!  Initialisation of pGFN-FF library        *
!********************************************
!
!  Call the subroutine to initialise GFN-FF
!
  call pgfnff_init
!
!  Call the subroutine to output input pGFN-FF parameters
!
  call pgfnff_outpar(ioout,ndim)
!********************************************
!  Generate neighbour list for this system  *
!********************************************
!
!  Find the maximum cutoff for the neighbour list 
!
  call pgfnff_get_max_cutoff(numat,nat,cut2_nbr)
!
!  Build the neighbour list
!
  call getnbr(ndim,rv,kv,numat,nat,x,y,z,cut2_nbr,ldebug)
!
!  Check that all atoms have neighbours where they should have
!
  call pgfnff_check_atom_with_nobonds(numat,nat,nnbr,lbondsok)
!
!  If bonding is not OK then try a larger cutoff radius
!
  if (.not.lbondsok) then
    cut2_nbr = 2.0_dp*cut2_nbr
    call getnbr(ndim,rv,kv,numat,nat,x,y,z,cut2_nbr,ldebug)
  endif
!********************************************
!  Generation of parameters in pGFN-FF      *
!********************************************
!
!  Call subroutine to generate parameters
!
  call pgfnff_pargen(ndim,kv,numat,nat,nftype,x,y,z,q,nfrag,nfraglist,qfrag, &
                     nnobo,nobond,nobotyp,maxnbr,nnbr,nbrno,ncnbr,rnbr,xnbr,ynbr,znbr, &
                     nnbr_bond,nbrno_bond,ncnbr_bond,rbnbr,xbnbr,ybnbr,zbnbr,lverbose)
!********************************************
!  Query parameters in pGFN-FF library      *
!********************************************
!
!  At this point the parameters have been generated and stored 
!  Add routines to compute the energy here if required......
!  Or just query the values for information.......
!
  call pgfnff_get_cnpar(cn_kn,hb_kn,gfnff_cnmax)
!
!  Find the maximum element number supported by the library
!
  call pgfnff_get_maxele(maxele)
!
!  Find the maximum atomic number in the actual system
!
  maxelein = 0
  do i = 1,numat
    maxelein = max(maxelein,nat(i))
  enddo
!
!  Covalent radii of elements
!
  allocate(gfnff_rcov(maxele))
  call pgfnff_get_covalent_radii(gfnff_rcov)
!
!  General radii of elements
!
  allocate(gfnff_rad(maxele))
  call pgfnff_get_general_radii(gfnff_rad)
!
!  Get parameters that control the variation of radii with coordination number
!
  allocate(gfnff_rad_cn(5,maxele))
  call pgfnff_get_cn_radii(maxele,gfnff_rad_cn)
!
!  Charge equilibration parameters
!
  allocate(gfnff_alp(numat))
  allocate(gfnff_chi(numat))
  allocate(gfnff_cnf(numat))
  allocate(gfnff_gam(numat))
  call pgfnff_get_eeq(numat,gfnff_alp,gfnff_chi,gfnff_gam,gfnff_cnf)
!
!  Repulsion parameters
!
  call pgfnff_get_repulsion_scale(gfnff_repscale_b,gfnff_repscale_n,gfnff_repscale_13,gfnff_repscale_14)
  call pgfnff_get_repulsion_threshold(gfnff_repthr)
  allocate(gfnff_repulsion_a(maxelein))
  allocate(gfnff_repulsion_z(maxelein))
  numat2 = numat*(numat+1)/2
  allocate(gfnff_repulsion_p(numat2))
  call pgfnff_get_repulsion(numat,maxelein,gfnff_repulsion_a,gfnff_repulsion_z,gfnff_repulsion_p)
!
!  Find the maximum ref value for dispersion in the library
!
  call pgfnff_get_maxref(maxrefin)
!
!  Dispersion parameters
!
  call pgfnff_get_dispersion_threshold(gfnff_dispthr)
  allocate(disp_nref(maxelein))
  allocate(disp_cn(maxrefin,maxelein))
  allocate(disp_c6(maxrefin,maxrefin,maxelein,maxelein))
  allocate(disp_c9(numat))
  allocate(disp_r0(maxelein*(maxelein+1)/2))
  allocate(disp_sqrtZr4r2(maxelein))
  allocate(disp_zeta(numat))
  call pgfnff_get_dispersion(numat,nat,maxelein,maxrefin,disp_nref,disp_cn,disp_c6,disp_c9, &
                             disp_zeta,disp_r0,disp_sqrtZr4r2)
!
!  Bond parameters
!
  allocate(par_bond(3,maxnbr,numat))
  call pgfnff_get_bond_parameters(numat,maxnbr,nnbr_bond,par_bond)
  call pgfnff_get_bond_scale(gfnff_bondscale)
!
!  Angle bend parameters
!
  call pgfnff_get_number_of_angles(nangles)
  allocate(nangleatomptr(nangles))
  allocate(par_angle(2,nangles))
  call pgfnff_get_angles(nangles,nangleatomptr,par_angle,gfnff_angle_damp)
!
!  Torsional parameters
!
  call pgfnff_get_number_of_torsions(ntorsions)
  allocate(ntorsionatomptr(ntorsions))
  allocate(par_torsion(2,ntorsions))
  call pgfnff_get_torsions(ntorsions,ntorsionatomptr,par_torsion,gfnff_torsion_damp)
!
!  Hydrogen bond thresholds
!
  call pgfnff_get_hydrogen_bond_thresholds(gfnff_hbthr1,gfnff_hbthr2)
!
!  Hydrogen bond cutoffs
!
  call pgfnff_get_hydrogen_bond_cutoffs(gfnff_hb_a_cut,gfnff_hb_long_cut,gfnff_hb_nb_cut, &
                                        gfnff_hb_s_cut,gfnff_hb_alp)
!
!  Halogen bond cutoffs
!
  call pgfnff_get_halogen_bond_cutoffs(gfnff_xb_a_cut,gfnff_xb_s_cut,gfnff_xb_l_cut)
!
!  Hydrogen bonding potential AB atoms
!
  call pgfnff_get_number_hydrogen_bond_AB(nABatoms)
  allocate(nABatomptr(nABatoms))
  allocate(gfnff_ABhbq(numat))
  allocate(gfnff_hb_acid(numat))
  allocate(gfnff_hb_base(numat))
  call pgfnff_get_hydrogen_bond_AB(numat,nABatoms,nABatomptr,gfnff_ABhbq,gfnff_hb_acid,gfnff_hb_base)
!
!  Hydrogen bonding potential H atoms
!
  call pgfnff_get_number_hydrogen_bond_H(nHatoms)
  allocate(nHatomptr(nHatoms))
  allocate(nHatomrptr(numat))
  call pgfnff_get_hydrogen_bond_H(numat,nHatoms,nHatomptr,nHatomrptr)
!
!  Hydrogen bond scale factors
!
  call pgfnff_get_hydrogen_bond_scale(gfnff_hb_scale_gen,gfnff_hb_scale_coh)
!
!  Hydrogen bond parameters
!
  call pgfnff_get_hydrogen_bond_parameters(gfnff_bend_hb,gfnff_tors_hb,gfnff_hbabmix)
!
!  Halogen bonding potential AB atoms
!
  call pgfnff_get_number_halogen_bond_AB(nxbABatoms)
  allocate(nxbABatomptr(3,nxbABatoms))
  allocate(gfnff_ABxbq(2,numat))
  call pgfnff_get_halogen_bond_AB(numat,nxbABatoms,nxbABatomptr,gfnff_ABxbq)
!
!  Halogen bond scale factors
!
  allocate(gfnff_xb_scale(maxele))
  call pgfnff_get_halogen_bond_scale(maxele,gfnff_xb_scale)

end program main
