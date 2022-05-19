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
  subroutine pgfnff_get_maxele(max_ele)
!
!  Query the maximum number of elements supported by the library
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: max_ele
!
  max_ele = max_gfnff_ele 
!
  end subroutine pgfnff_get_maxele
!
  subroutine pgfnff_get_maxref(max_ref)
!
!  Query the maximum number of reference states used by the library for dispersion
! 
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  use m_pgfnff_disp,  only : maxRef
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: max_ref
!
  max_ref = maxRef
!
  end subroutine pgfnff_get_maxref
!
  subroutine pgfnff_get_cnpar(cn_kn,hb_kn,gfnff_cnmax)
!
!  Query the coordination number parameters
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                        intent(out)      :: cn_kn
  real(dp),                        intent(out)      :: hb_kn
  real(dp),                        intent(out)      :: gfnff_cnmax
!
  cn_kn = gfnff_kn_cn
  hb_kn = gfnff_kn_hb
  gfnff_cnmax = cnmax
!
  end subroutine pgfnff_get_cnpar
!
  subroutine pgfnff_get_coordination_thresholds(gfnff_cnthr,gfnff_cnhbthr,gfnff_cnerfcut_hb,gfnff_cnerfcut_cn)
!
!  Query the coordination number thresholds
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                        intent(out)      :: gfnff_cnthr
  real(dp),                        intent(out)      :: gfnff_cnhbthr
  real(dp),                        intent(out)      :: gfnff_cnerfcut_hb
  real(dp),                        intent(out)      :: gfnff_cnerfcut_cn
!
  gfnff_cnthr = cnthr
  gfnff_cnhbthr = cnhbthr
  gfnff_cnerfcut_hb = cnerfcut_hb
  gfnff_cnerfcut_cn = cnerfcut_cn
!
  end subroutine pgfnff_get_coordination_thresholds
!
  subroutine pgfnff_get_covalent_radii(max_ele,grcov)
!
!  Query the covalent radii for elements up to max_ele
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)       :: max_ele
  real(dp),                        intent(out)      :: grcov(max_ele)
!
!  Local variables
!
  integer(i4)                                        :: n
!
  do n = 1,max_ele
    grcov(n) = rcov(n)
  enddo
!
  end subroutine pgfnff_get_covalent_radii

  subroutine pgfnff_get_general_radii(max_ele,grad)
!
!  Query the general radii for elements up to max_ele
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)       :: max_ele
  real(dp),                        intent(out)      :: grad(max_ele)
!
!  Local variables
!
  integer(i4)                                        :: n
!
  do n = 1,max_ele
    grad(n) = rad(n)*xtb_autoaa   ! NB Convert units for external use
  enddo
!
  end subroutine pgfnff_get_general_radii

  subroutine pgfnff_get_cn_radii(max_ele,gfnff_rad_cn)
!
!  Query the parameters that control how the radii for elements depend on coordination number
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)       :: max_ele
  real(dp),                        intent(out)      :: gfnff_rad_cn(5,max_ele)
!
!  Local variables
!
  integer(i4)                                        :: n
  integer(i4)                                        :: nrow
  real(dp),                                     save :: rowfct(6,2)
!
  data rowfct/29.84522887,-1.70549806, 6.54013762, 6.39169003, 6.00000000, 5.60000000, & ! First factor for each row
              -8.87843763, 2.10878369, 0.08009374,-0.85808076,-1.15000000,-1.30000000/   ! Second factor for each row
!
!
  do n = 1,max_ele
    nrow = PeriodicTableRow(n)
    gfnff_rad_cn(1,n) = r0_radij(n)
    gfnff_rad_cn(2,n) = cnfak_radij(n)
    gfnff_rad_cn(3,n) = en_radij(n)
    gfnff_rad_cn(4,n) = 0.005_dp*rowfct(nrow,1)
    gfnff_rad_cn(5,n) = 0.005_dp*rowfct(nrow,2)
  enddo
!
  end subroutine pgfnff_get_cn_radii
!*************************
!  Repulsion parameters  *
!*************************
  subroutine pgfnff_get_repulsion(numat,max_ele,gfnff_repulsion_a,gfnff_repulsion_z,gfnff_repulsion_p)
!
!  Query the repulsion parameters for elements
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg,      only : alphanb
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)       :: numat
  integer(i4),                     intent(in)       :: max_ele
  real(dp),                        intent(out)      :: gfnff_repulsion_a(max_ele)
  real(dp),                        intent(out)      :: gfnff_repulsion_z(max_ele)
  real(dp),                        intent(out)      :: gfnff_repulsion_p(numat*(numat+1)/2)
!
!  Local variables
!
  integer(i4)                                        :: n
  integer(i4)                                        :: ind
!
!  Check sizes
!
  if (max_ele.gt.max_gfnff_ele) then
    write(ioout,'('' Error: max_ele input to pgfnff_get_repulsion is too large '')')
    stop
  endif
!
  do n = 1,max_ele
    gfnff_repulsion_a(n) = repa(n)
    gfnff_repulsion_z(n) = repz(n)
  enddo
!
  ind = numat*(numat+1)/2
  do n = 1,ind
    gfnff_repulsion_p(n) = alphanb(n)
  enddo
!
  end subroutine pgfnff_get_repulsion

  subroutine pgfnff_get_repulsion_scale(gfnff_repscale_b,gfnff_repscale_n,gfnff_repscale_13,gfnff_repscale_14)
!
!  Get the repulsion scale factors
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                      intent(out)      :: gfnff_repscale_b
  real(dp),                      intent(out)      :: gfnff_repscale_n
  real(dp),                      intent(out)      :: gfnff_repscale_13
  real(dp),                      intent(out)      :: gfnff_repscale_14
!
  gfnff_repscale_b = repscalb
  gfnff_repscale_n = repscaln
  gfnff_repscale_13 = hh13rep
  gfnff_repscale_14 = hh14rep
!
  end subroutine pgfnff_get_repulsion_scale

  subroutine pgfnff_get_repulsion_threshold(gfnff_repthr)
!
!  Get the repulsion threshold value
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                      intent(out)      :: gfnff_repthr
!
  gfnff_repthr = repthr
!
  end subroutine pgfnff_get_repulsion_threshold
!**************************
!  Dispersion parameters  *
!**************************
  subroutine pgfnff_get_dispersion(numat,nat,maxelein,maxrefin,nref,cn,c6,c9,zeta,r0,rsqrtZr4r2)
!
!  Query the dispersion parameters for elements and atoms
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)       :: numat                                   ! Number of atoms
  integer(i4),                     intent(in)       :: nat(numat)                              ! Atomic numbers of atoms
  integer(i4),                     intent(in)       :: maxelein                                ! Maximum atomic number in system
  integer(i4),                     intent(in)       :: maxrefin                                ! Maximum number of ref values for elements 
  integer(i4),                     intent(out)      :: nref(maxelein)                          ! Number of ref values for each element
  real(dp),                        intent(out)      :: cn(maxrefin,maxelein)                   ! Coefficients of ref values for each element used for CN dependence
  real(dp),                        intent(out)      :: c6(maxrefin,maxrefin,maxelein,maxelein) ! Pairwise coefficients for computing C6 terms
  real(dp),                        intent(out)      :: c9(numat)                               ! Atom-specific coefficients for C9
  real(dp),                        intent(out)      :: zeta(numat)                             ! Atom specific zeta values for C6
  real(dp),                        intent(out)      :: r0(maxelein*(maxelein+1)/2)             ! Element pairwise r0 values 
  real(dp),                        intent(out)      :: rsqrtZr4r2(maxelein)                    ! Element specific sqrtZr4r2 values
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: ind
  integer(i4)                                       :: j
  integer(i4)                                       :: m
  integer(i4)                                       :: n
  integer(i4)                                       :: natm
  real(dp)                                          :: fm
!
!  Check sizes
!
  if (maxRef.gt.maxrefin) then
    write(ioout,'('' Error: Arrays for dispersion parameters are too small : cn, c6 '')')
    stop
  endif
  if (maxEle.gt.maxelein) then
    write(ioout,'('' Error: Arrays for dispersion parameters are too small : nref, cn, c6 '')')
    stop
  endif
!
  nref(1:maxelein) = 0
  cn = 0.0_dp
  c6 = 0.0_dp
!
  do m = 1,numat
    natm = nat(m)
    if (nref(natm).eq.0) then
      nref(natm) = d4_nref(natm)
    endif
  enddo
  ind = 0
  do m = 1,maxelein
    if (nref(m).gt.0) then
      do i = 1,nref(m)
        cn(i,m) = d4_cn(i,m)
      enddo
      do n = 1,maxelein
        if (nref(n).gt.0) then
          do i = 1,nref(m)
            do j = 1,nref(n)
              c6(j,i,n,m) = d4_c6(j,i,n,m)
            enddo
          enddo
        endif
      enddo
    endif
    rsqrtZr4r2(m) = sqrtZr4r2(m)
    do n = 1,m
      ind = ind + 1
      r0(ind) = d3r0(ind)
    enddo
  enddo
!
  do m = 1,numat
    fm = (1.0_dp - 3.0_dp*qf0(m))
    fm = min(max(fm,-4.0_dp),4.0_dp)
    c9(m) = fm*zb3atm(nat(m))*(xtb_autoev)**(1.0_dp/3.0_dp)  ! Charge-corrected C9 term with units of eV for triple product
    zeta(m) = d4_zeta(m)
  enddo
!
  end subroutine pgfnff_get_dispersion

  subroutine pgfnff_get_dispersion_threshold(gfnff_dispthr)
!
!  Get the dispersion threshold value
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                      intent(out)      :: gfnff_dispthr
!
  gfnff_dispthr = dispthr
!
  end subroutine pgfnff_get_dispersion_threshold

!****************
!  Angle bends  *
!****************
  subroutine pgfnff_get_number_of_angles(nangles)
!
!  Query the number of angle bend terms in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: nangles
!
  nangles = nangle
!
  end subroutine pgfnff_get_number_of_angles

  subroutine pgfnff_get_angles(maxangles,nangleatomptr,par_angle,gfnff_angle_damp)
!
!  Query the angle parameters. Arrays should have already been dimensioned
!  correctly otherwise an error will occur.
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types   
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: maxangles                  ! Right-hand dimension of arrays
  integer(i4),                      intent(out)      :: nangleatomptr(3,maxangles) ! Pointer to atoms involved in the angle
  real(dp),                         intent(out)      :: par_angle(2,maxangles)     ! Parameters for angles
  real(dp),                         intent(out)      :: gfnff_angle_damp           ! Damping parameter for angles
!
!  Local variables
!
  integer(i4)                                        :: n
!
!  Check array dimensions
!
  if (nangle.gt.maxangles) then
    write(ioout,'('' Error: Arrays for angles are too small : nangleatomptr, par_angle '')')
    stop
  endif
  do n = 1,nangle
    nangleatomptr(1:3,n) = nangleptr(1:3,n)
    par_angle(1:2,n) = vangle(1:2,n)
  enddo
  gfnff_angle_damp = atcuta
!
  end subroutine pgfnff_get_angles
!*************
!  Torsions  *
!*************
  subroutine pgfnff_get_number_of_torsions(ntorsions)
!
!  Query the number of torsional terms in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types   
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: ntorsions
!
  ntorsions = ntors
!
  end subroutine pgfnff_get_number_of_torsions

  subroutine pgfnff_get_torsions(maxtorsions,ntorsionatomptr,par_torsion,gfnff_torsion_damp)
!
!  Query the torsion parameters. Arrays should have already been dimensioned
!  correctly otherwise an error will occur.
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: maxtorsions                    ! Right-hand dimension of arrays
  integer(i4),                      intent(out)      :: ntorsionatomptr(5,maxtorsions) ! Pointer to atoms involved in the torsion
  real(dp),                         intent(out)      :: par_torsion(2,maxtorsions)     ! Parameters for torsions
  real(dp),                         intent(out)      :: gfnff_torsion_damp             ! Damping parameter for torsions
!
!  Local variables
!
  integer(i4)                                        :: n
!
!  Check array dimensions
!
  if (ntors.gt.maxtorsions) then
    write(ioout,'('' Error: Arrays for torsions are too small : ntorsionatomptr, par_torsion '')')
    stop
  endif
  do n = 1,ntors
    ntorsionatomptr(1:5,n) = ntorsptr(1:5,n)
    par_torsion(1:2,n) = vtors(1:2,n)
  enddo
  gfnff_torsion_damp = atcutt
!
  end subroutine pgfnff_get_torsions

  subroutine pgfnff_get_eeq(numat,gfnff_alp,gfnff_chi,gfnff_gam,gfnff_cnf)
!
!  Query the charge equilibration parameters. Arrays should have already been dimensioned
!  correctly otherwise an error will occur.
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat                ! Number of atoms
  real(dp),                         intent(out)      :: gfnff_alp(numat)     ! Alpha parameters for EEQ
  real(dp),                         intent(out)      :: gfnff_chi(numat)     ! Chi parameters for EEQ
  real(dp),                         intent(out)      :: gfnff_gam(numat)     ! Gamma parameters for EEQ
  real(dp),                         intent(out)      :: gfnff_cnf(numat)     ! Coordination number parameters for EEQ
!
!  Local variables
!
  integer(i4)                                        :: n
!
  do n = 1,numat
    gfnff_alp(n) = alpeeq(n)
    gfnff_chi(n) = chieeq(n)
    gfnff_gam(n) = gameeq(n)
    gfnff_cnf(n) = cnfeeq(n)
  enddo
!
  end subroutine pgfnff_get_eeq

  subroutine pgfnff_get_bond_parameters(numat,maxnbr,nbnbr,par_bond)
!
!  Query the bond parameters. Arrays should have already been dimensioned
!  correctly otherwise an error will occur.
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff_cfg,  only : maxat => maxat_pgfnff
  use m_pgfnff_nbr_lib
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: maxnbr                         ! Array size parameter
  integer(i4),                      intent(in)       :: numat                          ! Number of atoms
  integer(i4),                      intent(in)       :: nbnbr(numat)                   ! Number of neighbours
  real(dp),                         intent(out)      :: par_bond(3,maxnbr,numat)       ! Parameters for bonds
!
!  Local variables
!
  integer(i4)                                        :: i
  integer(i4)                                        :: j
!
!  Check array dimensions
!
  if (maxnbr.gt.maxnbr_lib) then
    write(ioout,'('' Error: Arrays for bond parameters are too small : par_bond '')')
    stop
  endif
  if (numat.gt.maxat) then
    write(ioout,'('' Error: Arrays for bond parameters are too small : nbnbr, par_bond '')')
    stop
  endif
  do i = 1,numat
    do j = 1,nbnbr(i)
      par_bond(1:3,j,i) = vbnbr(1:3,j,i)
    enddo
  enddo
!
  end subroutine pgfnff_get_bond_parameters

  subroutine pgfnff_get_bond_scale(gfnff_bondscale)
!
!  Get the bond scale factor
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                      intent(out)      :: gfnff_bondscale
!
  gfnff_bondscale = vbond_scale
!
  end subroutine pgfnff_get_bond_scale
!***************************
!  Halogen bond AB atoms  *
!***************************
  subroutine pgfnff_get_number_halogen_bond_AB(n_gfnff_xb_AB)
!
!  Query the number of potential AB pairs for halogen bonds in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: n_gfnff_xb_AB
!
  n_gfnff_xb_AB = natxbAB
!
  end subroutine pgfnff_get_number_halogen_bond_AB

  subroutine pgfnff_get_halogen_bond_AB(numat,n_gfnff_xb_AB,n_gfnff_xb_ABptr,gfnff_xb_ABq)
!
!  Query the potential AB pairs for halogen bonds in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat
  integer(i4),                      intent(in)       :: n_gfnff_xb_AB
  integer(i4),                      intent(out)      :: n_gfnff_xb_ABptr(3,n_gfnff_xb_AB)
  real(dp),                         intent(out)      :: gfnff_xb_ABq(2,numat)
!
!  Local variables
!
  integer(i4)                                        :: n
!
!  Trap case if number of terms is zero 
!
  if (n_gfnff_xb_AB.eq.0) then
    gfnff_xb_ABq(1:2,1:numat) = 0.0_dp
  else
    do n = 1,n_gfnff_xb_AB
      n_gfnff_xb_ABptr(1:3,n) = xbatABl(1:3,n)
    enddo
    do n = 1,numat
      gfnff_xb_ABq(1:2,n) = ABxbq(1:2,n)
    enddo
  endif
!
  end subroutine pgfnff_get_halogen_bond_AB
!***************************
!  Hydrogen bond AB atoms  *
!***************************
  subroutine pgfnff_get_number_hydrogen_bond_AB(n_gfnff_hb_AB)
!
!  Query the number of potential AB pairs for hydrogen bonds in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: n_gfnff_hb_AB
!
  n_gfnff_hb_AB = nABat
!
  end subroutine pgfnff_get_number_hydrogen_bond_AB

  subroutine pgfnff_get_hydrogen_bond_AB(numat,n_gfnff_hb_AB,n_gfnff_hb_ABptr,gfnff_hb_ABq,gfnff_hb_acid,gfnff_hb_base)
!
!  Query the potential AB pairs in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat
  integer(i4),                      intent(in)       :: n_gfnff_hb_AB
  integer(i4),                      intent(out)      :: n_gfnff_hb_ABptr(n_gfnff_hb_AB)
  real(dp),                         intent(out)      :: gfnff_hb_ABq(numat)
  real(dp),                         intent(out)      :: gfnff_hb_acid(numat)
  real(dp),                         intent(out)      :: gfnff_hb_base(numat)
!
!  Local variables
!
  integer(i4)                                        :: n
!
!  Trap case if number of terms is zero 
!
  if (n_gfnff_hb_AB.eq.0) then
    gfnff_hb_ABq(1:numat) = 0.0_dp
  else
    do n = 1,n_gfnff_hb_AB
      n_gfnff_hb_ABptr(n) = nABatptr(n)
    enddo
    do n = 1,numat
      gfnff_hb_ABq(n) = ABhbq(n)
      gfnff_hb_acid(n) = hbacid(n)
      gfnff_hb_base(n) = hbbase(n)
    enddo
  endif
!
  end subroutine pgfnff_get_hydrogen_bond_AB

!**************************
!  Hydrogen bond H atoms  *
!**************************
  subroutine pgfnff_get_number_hydrogen_bond_H(n_gfnff_hb_H)
!
!  Query the number of potential H atoms for hydrogen bonds in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(out)      :: n_gfnff_hb_H
!
  n_gfnff_hb_H = nathbH
!
  end subroutine pgfnff_get_number_hydrogen_bond_H

  subroutine pgfnff_get_hydrogen_bond_H(numat,n_gfnff_hb_H,n_gfnff_hb_Hptr,n_gfnff_hb_Hrptr)
!
!  Query the potential H atoms in pGFNFF
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat
  integer(i4),                      intent(in)       :: n_gfnff_hb_H
  integer(i4),                      intent(out)      :: n_gfnff_hb_Hptr(n_gfnff_hb_H)
  integer(i4),                      intent(out)      :: n_gfnff_hb_Hrptr(numat)
!
!  Local variables
!
  integer(i4)                                        :: n
!
  do n = 1,n_gfnff_hb_H
    n_gfnff_hb_Hptr(n) = hbatHl(n)
  enddo
  do n = 1,numat
    n_gfnff_hb_Hrptr(n) = rhbatHl(n)
  enddo
!
  end subroutine pgfnff_get_hydrogen_bond_H

  subroutine pgfnff_get_hydrogen_bond_thresholds(gfnff_hbthr1,gfnff_hbthr2)
!
!  Query the hydrogen bond threshold values
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                         intent(out)      :: gfnff_hbthr1
  real(dp),                         intent(out)      :: gfnff_hbthr2
!
  gfnff_hbthr1 = hbthr1
  gfnff_hbthr2 = hbthr2
!
  end subroutine pgfnff_get_hydrogen_bond_thresholds

  subroutine pgfnff_get_hydrogen_bond_cutoffs(gfnff_hb_a_cut,gfnff_hb_long_cut,gfnff_hb_nb_cut,gfnff_hb_s_cut, &
                                              gfnff_hb_alp)
!
!  Query the hydrogen bond threshold values
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                         intent(out)      :: gfnff_hb_a_cut
  real(dp),                         intent(out)      :: gfnff_hb_long_cut
  real(dp),                         intent(out)      :: gfnff_hb_nb_cut
  real(dp),                         intent(out)      :: gfnff_hb_s_cut
  real(dp),                         intent(out)      :: gfnff_hb_alp
!
  gfnff_hb_a_cut = hbacut
  gfnff_hb_long_cut = hblongcut
  gfnff_hb_nb_cut = hbnbcut
  gfnff_hb_s_cut = hbscut
  gfnff_hb_alp = hbalp
!
  end subroutine pgfnff_get_hydrogen_bond_cutoffs

  subroutine pgfnff_get_hydrogen_bond_scale(gfnff_hb_scale_gen,gfnff_hb_scale_coh)
!
!  Query the hydrogen bond scaling parameters
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                        intent(out)      :: gfnff_hb_scale_gen
  real(dp),                        intent(out)      :: gfnff_hb_scale_coh
!
  gfnff_hb_scale_gen = xhaci_globabh
  gfnff_hb_scale_coh = xhaci_coh
!
  end subroutine pgfnff_get_hydrogen_bond_scale

  subroutine pgfnff_get_hydrogen_bond_parameters(gfnff_bend_hb,gfnff_tors_hb,gfnff_hbabmix)
!
!  Query the hydrogen bond parameters
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                        intent(out)      :: gfnff_bend_hb
  real(dp),                        intent(out)      :: gfnff_tors_hb
  real(dp),                        intent(out)      :: gfnff_hbabmix
!
  gfnff_bend_hb = bend_hb
  gfnff_tors_hb = tors_hb
  gfnff_hbabmix = hbabmix
!
  end subroutine pgfnff_get_hydrogen_bond_parameters

  subroutine pgfnff_get_halogen_bond_cutoffs(gfnff_xb_a_cut,gfnff_xb_s_cut,gfnff_xb_l_cut)
!
!  Query the halogen bond threshold values
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  real(dp),                         intent(out)      :: gfnff_xb_a_cut
  real(dp),                         intent(out)      :: gfnff_xb_s_cut
  real(dp),                         intent(out)      :: gfnff_xb_l_cut
!
  gfnff_xb_a_cut = xbacut
  gfnff_xb_s_cut = xbscut
  gfnff_xb_l_cut = hblongcut_xb
!
  end subroutine pgfnff_get_halogen_bond_cutoffs

  subroutine pgfnff_get_halogen_bond_scale(max_ele,gfnff_xb_scale)
!
!  Query the halogen bond scaling parameters for elements
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)       :: max_ele
  real(dp),                        intent(out)      :: gfnff_xb_scale(max_ele)
!
!  Local variables
!
  integer(i4)                                        :: n
!
!  Check sizes
!
  if (max_ele.gt.max_gfnff_ele) then
    write(ioout,'('' Error: max_ele input to pgfnff_get_halogen_bond_scale is too large '')')
    stop
  endif
!
  do n = 1,max_ele
    gfnff_xb_scale(n) = xbaci(n)
  enddo
!
  end subroutine pgfnff_get_halogen_bond_scale
