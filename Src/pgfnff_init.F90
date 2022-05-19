!
! Original GFN-FF code:
! ---------------------
!
! Copyright (C) 2017-2019 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.
!
! Modifications for pGFNFF library:
! ---------------------------------
!
! Julian Gale, Curtin University, 2020-2022
!
  subroutine pgfnff_init
!
!  Initialises quantities for pGFNFF at start up
!
  use m_pgfnff
  implicit none
!
!  Local variables
!
  integer(i4)              :: i
  integer(i4)              :: j
  integer(i4)              :: k
  real(dp)                 :: tmp
!
!  Basicities
!
  xhbas(1:max_gfnff_ele) = 0.0_dp
!
  xhbas( 6) = 0.80_dp
  xhbas( 7) = 1.68_dp
  xhbas( 8) = 0.67_dp
  xhbas( 9) = 0.52_dp
  xhbas(14) = 4.0_dp
  xhbas(15) = 3.5_dp
  xhbas(16) = 2.0_dp
  xhbas(17) = 1.5_dp
  xhbas(35) = 1.5_dp
  xhbas(53) = 1.9_dp
  xhbas(33) = xhbas(15)
  xhbas(34) = xhbas(16)
  xhbas(51) = xhbas(15)
  xhbas(52) = xhbas(16)
!
!  HB Acidities
!
  xhaci(1:max_gfnff_ele) = 0.0_dp
!
  xhaci( 6) = 0.75_dp
  xhaci( 7) = xhaci_glob + 0.1_dp
  xhaci( 8) = xhaci_glob
  xhaci( 9) = xhaci_glob
  xhaci(15) = xhaci_glob
  xhaci(16) = xhaci_glob
  xhaci(17) = xhaci_glob + 1.0_dp
  xhaci(35) = xhaci_glob + 1.0_dp
  xhaci(53) = xhaci_glob + 1.0_dp
!
!  XB Acidities
!
  xbaci(1:max_gfnff_ele) = 0.0_dp
!
  xbaci(15) = 1.0_dp
  xbaci(16) = 1.0_dp
  xbaci(17) = 0.5_dp
  xbaci(33) = 1.2_dp
  xbaci(34) = 1.2_dp
  xbaci(35) = 0.9_dp
  xbaci(51) = 1.2_dp
  xbaci(52) = 1.2_dp
  xbaci(53) = 1.2_dp
!
!  Generate sqrtZr4r2
!
  do i = 1,max_gfnff_ele
    sqrtZr4r2(i) = sqrt(0.5_dp*r4Overr2(i)*sqrt(dble(i)))
  enddo
!
!  3B bond prefactors and D3 stuff
!
  k = 0
  do i = 1,max_gfnff_ele
    tmp = dble(i)
    zb3atm(i) = - tmp*batmscal**(1.0_dp/3.0_dp)  ! include pre-factor
    do j = 1,i
      k = k + 1
      tmp = sqrtZr4r2(i)*sqrtZr4r2(j)*3.0_dp
      d3r0(k) = (d3a1*dsqrt(tmp) + d3a2)**2   ! save R0^2 for efficiency reasons
    enddo
  enddo
  zb3atm(1) = - 0.25_dp*batmscal**(1.0_dp/3.0_dp) ! slightly better than 1.0
!
!  Bond strength
!
  bstren(1) = 1.00_dp      ! single bond
  bstren(2) = 1.24_dp      ! double bond
  bstren(3) = 1.98_dp      ! triple bond
  bstren(4) = 1.22_dp      ! hyperval bond
  bstren(5) = 1.00_dp      ! M-X
  bstren(6) = 0.78_dp      ! M eta
  bstren(7) = 3.40_dp      ! M-M
  bstren(8) = 3.40_dp      ! M-M
  bstren(9) = 0.5_dp*(bstren(7)+bstren(8))
!
  bsmat(0:3,0:3) = -999.0_dp
  bsmat(0,0) = bstren(1)
  bsmat(3,0) = bstren(1)
  bsmat(3,3) = bstren(1)
  bsmat(2,2) = bstren(2)
  bsmat(1,1) = bstren(3)
  bsmat(1,0) = split0*bstren(1) + split1*bstren(3)
  bsmat(3,1) = split0*bstren(1) + split1*bstren(3)
  bsmat(2,1) = split0*bstren(2) + split1*bstren(3)
  bsmat(2,0) = split0*bstren(1) + split1*bstren(2)
  bsmat(3,2) = split0*bstren(1) + split1*bstren(2)
!
!  Torsion parameters
!
  torsf(1) = 1.00         ! single bond
  torsf(2) = 1.18         ! pi bond
  torsf(3) = 1.05         ! improper
  torsf(5) = 0.50         ! pi part improper
  torsf(6) = -0.90        ! extra sp3 C
  torsf(7) =  0.70        ! extra sp3 N
  torsf(8) = -2.00        ! extra sp3 O
!
!  Bond charge parameters
!
  qfacbm(0)   =1.0d0    ! bond charge dep.gff_srcs += 'gff/gfnff_input.f90'
  qfacbm(1:2) =-0.2d0   !
  qfacbm(  3) =0.70d0   !
  qfacbm(  4) =0.50d0   !
!
!  Huckel 
!
  hdiag(5) =-0.5d0      ! diagonal element relative to c
  hdiag(6) =0.00d0      !
  hdiag(7) =0.14d0      !
  hdiag(8) =-0.38d0     !
  hdiag(9) =-0.29d0     !
  hdiag(16)=-0.30d0     !
  hdiag(17)=-0.30d0     !
  hoffdiag(5)=0.5d0     ! Huckel off-diag constants
  hoffdiag(6)=1.00d0    !
  hoffdiag(7)=0.66d0    !
  hoffdiag(8)=1.10d0    !
  hoffdiag(9)=0.23d0    !
  hoffdiag(16)=0.60d0   !
  hoffdiag(17)=1.00d0   !
!
!  Scaling factor for radii
!
  fat(1:max_gfnff_ele) = 1.0_dp
!
  fat( 1) = 1.02_dp
  fat( 4) = 1.03_dp
  fat( 5) = 1.02_dp
  fat( 8) = 1.02_dp
  fat( 9) = 1.05_dp
  fat(10) = 1.10_dp
  fat(11) = 1.01_dp
  fat(12) = 1.02_dp
  fat(15) = 0.97_dp
  fat(18) = 1.10_dp
  fat(19) = 1.02_dp
  fat(20) = 1.02_dp
  fat(38) = 1.02_dp
  fat(34) = 0.99_dp
  fat(50) = 1.01_dp
  fat(51) = 0.99_dp
  fat(52) = 0.95_dp
  fat(53) = 0.98_dp
  fat(56) = 1.02_dp
  fat(76) = 1.02_dp
  fat(82) = 1.06_dp
  fat(83) = 0.95_dp
!
!  Set periodic table row pointers
!
  PeriodicTableRow(1:2) = 1
  PeriodicTableRow(3:10) = 2
  PeriodicTableRow(11:18) = 3
  PeriodicTableRow(19:36) = 4
  PeriodicTableRow(37:54) = 5
  PeriodicTableRow(55:86) = 6
!
!  Correct the D3 covalent radii by 4/3
!
  do i = 1,max_gfnff_ele
    rcov(i) = 4.0_dp*rcov(i)/3.0_dp
  enddo
!
!  Useful constants
!
  tworootpi = 1.0_dp/sqrt(atan(1.0_dp))
!
!  Set thresholds based on the accuracy (except for cn-hb case)
!
  cnhbthr =  900.0_dp
  dispthr = 1500.0_dp - log10(gfnff_accuracy_disp)*1000.0_dp
  cnthr   =  100.0_dp - log10(gfnff_accuracy_cn)*50.0_dp
  repthr  =  400.0_dp - log10(gfnff_accuracy_rep)*100.0_dp
  hbthr1  =  200.0_dp - log10(gfnff_accuracy_hb1)*50.0_dp
  hbthr2  =  400.0_dp - log10(gfnff_accuracy_hb2)*50.0_dp
!
!  Set start of taper for each threshold
!
  t_cnhbthr = cnhbthr*gfnff_taper
  t_dispthr = dispthr*gfnff_taper
  t_cnthr   = cnthr*gfnff_taper
  t_repthr  = repthr*gfnff_taper
  t_hbthr1  = hbthr1*gfnff_taper
  t_hbthr2  = hbthr2*gfnff_taper
!
!  Set scale factor r0**2 used in coordination numbers so that erf term is non-zero
!
  cnerfcut_hb = (1.0_dp + gfnff_erftol/gfnff_kn_hb)**2
  cnerfcut_cn = (1.0_dp + gfnff_erftol/(-gfnff_kn_cn))
!#######################################################################################
!  Units:                                                                              #
!  Change units of GFNFF parameters from original ones in au to ev and ang for pGFNFF  #
!#######################################################################################

  cnthr                        = cnthr*(xtb_autoaa)**2                             ! Distance**2
  t_cnthr                      = t_cnthr*(xtb_autoaa)**2                           ! Distance**2
  cnhbthr                      = cnhbthr*(xtb_autoaa)**2                           ! Distance**2
  t_cnhbthr                    = t_cnhbthr*(xtb_autoaa)**2                         ! Distance**2
  dispthr                      = dispthr*(xtb_autoaa)**2                           ! Distance**2
  t_dispthr                    = t_dispthr*(xtb_autoaa)**2                         ! Distance**2
  repthr                       = repthr*(xtb_autoaa)**2                            ! Distance**2
  t_repthr                     = t_repthr*(xtb_autoaa)**2                          ! Distance**2
  alp(1:max_gfnff_ele)         = alp(1:max_gfnff_ele)*xtb_autoaa                   ! Distance
  chi(1:max_gfnff_ele)         = chi(1:max_gfnff_ele)*xtb_autoev                   ! Energy
  cnf_gfnff(1:max_gfnff_ele)   = cnf_gfnff(1:max_gfnff_ele)*xtb_autoev             ! Energy
  gam(1:max_gfnff_ele)         = gam(1:max_gfnff_ele)/xtb_autoaa                   ! Inverse distance
  cnfak_radij(1:max_gfnff_ele) = cnfak_radij(1:max_gfnff_ele)*xtb_autoaa           ! Distance
  r0_radij(1:max_gfnff_ele)    = r0_radij(1:max_gfnff_ele)*xtb_autoaa              ! Distance
  repa(1:max_gfnff_ele)        = repa(1:max_gfnff_ele)/(xtb_autoaa**1.5)           ! Inverse distance**(-1.5)
  repscalb                     = repscalb*xtb_autoev*xtb_autoaa                    ! Energy*distance
  repscaln                     = repscaln*xtb_autoev*xtb_autoaa                    ! Energy*distance
  rqshrink                     = rqshrink*xtb_autoaa                               ! Distance
  sqrtZr4r2(1:max_gfnff_ele)   = sqrtZr4r2(1:max_gfnff_ele)*xtb_autoaa             ! Distance
  xbaci(1:max_gfnff_ele)       = xbaci(1:max_gfnff_ele)*xtb_autoev*xtb_autoaa**3   ! Energy*distance**3
  zb3atm(1:max_gfnff_ele)      = zb3atm(1:max_gfnff_ele)*xtb_autoaa**3             ! Distance**3
!
  d3r0 = d3r0*xtb_autoaa**2                                                        ! Distance**2
  hbacut = hbacut*xtb_autoaa                                                       ! Distance
  hblongcut = hblongcut*xtb_autoaa**2                                              ! Distance**2
  hblongcut_xb = hblongcut_xb*xtb_autoaa**2                                        ! Distance**2
  hbnbcut = hbnbcut*xtb_autoaa                                                     ! Distance
  hbscut = hbscut*xtb_autoaa                                                       ! Distance
  hbthr1 = hbthr1*xtb_autoaa**2                                                    ! Distance**2
  hbthr2 = hbthr2*xtb_autoaa**2                                                    ! Distance**2
  t_hbthr1 = t_hbthr1*xtb_autoaa**2                                                ! Distance**2
  t_hbthr2 = t_hbthr2*xtb_autoaa**2                                                ! Distance**2
  mchishift = mchishift*xtb_autoev                                                 ! Energy
  xbscut = xbscut*xtb_autoaa                                                       ! Distance
!
  end subroutine pgfnff_init
