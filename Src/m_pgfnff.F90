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
module m_pgfnff
! 
!  This module contains the variables and subroutine associated
!  with the initialisation of parameters for (p)GFNFF
!
  use m_pgfnff_types
!
  implicit none
!
  integer(i4),     parameter     :: max_gfnff_ele = 86       ! Number of elements supported by GFNFF
!
  real(dp),                 save :: r4Overr2(max_gfnff_ele)  ! Quantities needed to generate dispersion corrections
!
!  Conversion factors for compatability with XTB
!
  real(dp),        parameter     :: xtb_autoaa = 0.52917726_dp
  real(dp),        parameter     :: xtb_autoev = 27.21138505_dp
!
!  Version of GFNFF string
!
  character(len=6),         save :: gfnff_version = 'v6.3.3'        ! Version of GFNN in case of changes
!
!  Pi system flag
!
  integer(i4),              save :: gfnff_pi_change = 0             ! Flag for type of pi electron change : 0 => use charge; + or - for direction
  real(dp),                 save :: gfnff_pi_temp1 = 4000.0_dp      ! Temperature for Huckel pi calculation on first pass
  real(dp),                 save :: gfnff_pi_temp2 =  300.0_dp      ! Temperature for Huckel pi calculation on second pass
  real(dp),                 save :: gfnff_ks = 0.04_dp              ! Maximum spacing between points in reciprocal space
!
!  Trap for negative charges that may cause problems in the topology
!
  real(dp),                 save :: gfnff_q_trap = -2.000_dp        ! Largest negative charge allowed for gfnff_radij
!
!  Trap high coordination numbers
!
  logical,                  save :: lgfnff_highcn_trap = .false.    ! Controls whether to trap high coordination numbers as per XTB
!
!  Flag to control setting of fragments
!
  logical,                  save :: lgfnff_fragment_bond = .false.  ! Controls whether to enforce use of bonds to set fragments
!
!  Accuracy and thresholds
!
  real(dp),                 save :: max_gfnff_accuracy = 31.622_dp  ! Maximum possible overall accuracy value for GFNFF
  real(dp),                 save :: max_gfnff_acc_disp = 31.622_dp  ! Maximum possible accuracy value for GFNFF for dispersion
  real(dp),                 save :: max_gfnff_acc_rep  = 10000.0_dp ! Maximum possible accuracy value for GFNFF for repulsion
  real(dp),                 save :: max_gfnff_acc_cn   = 100.0_dp   ! Maximum possible accuracy value for GFNFF for coordination numbers
  real(dp),                 save :: max_gfnff_acc_hb1  = 10000.0_dp ! Maximum possible accuracy value for GFNFF for hydrogen bonding 1
  real(dp),                 save :: max_gfnff_acc_hb2  = 1.0d8      ! Maximum possible accuracy value for GFNFF for hydrogen bonding 2
  real(dp),                 save :: gfnff_accuracy = 0.1_dp         ! Overall accuracy value for GFNFF
  real(dp),                 save :: gfnff_accuracy_disp = 0.1_dp    ! Accuracy value for GFNFF for dispersion
  real(dp),                 save :: gfnff_accuracy_rep = 0.1_dp     ! Accuracy value for GFNFF for repulsion
  real(dp),                 save :: gfnff_accuracy_cn = 0.1_dp      ! Accuracy value for GFNFF for coordination numbers
  real(dp),                 save :: gfnff_accuracy_hb1 = 0.1_dp     ! Accuracy value for GFNFF for hydrogen bonding 1
  real(dp),                 save :: gfnff_accuracy_hb2 = 0.1_dp     ! Accuracy value for GFNFF for hydrogen bonding 2
  real(dp),                 save :: gfnff_cnc6tol = 1.0d-10         ! Tolerance for reducing the coordination number 2nd derivatives during dispersion
  real(dp),                 save :: gfnff_erftol = 6.0_dp           ! Tolerance that limits the calculation of error functions 
  real(dp),                 save :: gfnff_kn_hb = 27.5_dp           ! Constant used in coordination numbers for hydrogen bonds
  real(dp),                 save :: gfnff_kn_cn = -7.5_dp           ! Constant used in coordination numbers
  real(dp),                 save :: gfnff_taper = 0.95_dp           ! Fractional taper range for GFNFF
  real(dp),                 save :: dispthr                         ! Dispersion threshold
  real(dp),                 save :: cnthr                           ! Coordination number threshold
  real(dp),                 save :: cnhbthr                         ! Coordination number (HB) threshold
  real(dp),                 save :: cnerfcut_hb                     ! Coordination number threshold scale factor for r0 (HB)
  real(dp),                 save :: cnerfcut_cn                     ! Coordination number threshold scale factor for r0 (CN)
  real(dp),                 save :: repthr                          ! Repulsion threshold
  real(dp),                 save :: hbthr1                          ! Hydrogen bond threshold 1
  real(dp),                 save :: hbthr2                          ! Hydrogen bond threshold 2
  real(dp),                 save :: t_dispthr                       ! Dispersion taper
  real(dp),                 save :: t_cnthr                         ! Coordination number taper
  real(dp),                 save :: t_cnhbthr                       ! Coordination number (HB) taper
  real(dp),                 save :: t_repthr                        ! Repulsion taper
  real(dp),                 save :: t_hbthr1                        ! Hydrogen bond taper 1
  real(dp),                 save :: t_hbthr2                        ! Hydrogen bond taper 2
!
  integer(i4),              save :: group(max_gfnff_ele)            ! Group of periodic table
  integer(i4),              save :: normcn(max_gfnff_ele)           ! Normal coordination number
  integer(i4),              save :: metal(max_gfnff_ele)            ! Metal (1) or transition metal (2) flag
  integer(i4),              save :: PeriodicTableRow(max_gfnff_ele) ! Pointer from element to periodic table row
  real(dp),                 save :: alp(max_gfnff_ele)              ! Alpha
  real(dp),                 save :: angl(max_gfnff_ele)             ! Angl
  real(dp),                 save :: angl2(max_gfnff_ele)            ! Angl2
  real(dp),                 save :: atomicRad(max_gfnff_ele)        ! Atomic radii used for molecular fragment setup
  real(dp),                 save :: bond(max_gfnff_ele)             ! Bond
  real(dp),                 save :: chi(max_gfnff_ele)              ! Chi for EEM
  real(dp),                 save :: cnf_gfnff(max_gfnff_ele)        ! CNF 
  real(dp),                 save :: cnfak_radij(max_gfnff_ele)      ! CN factor radij 
  real(dp),                 save :: en(max_gfnff_ele)               ! Pauling electronegativity
  real(dp),                 save :: en_radij(max_gfnff_ele)         ! For radij
  real(dp),                 save :: fat(max_gfnff_ele)              ! Scaling factors for radii
  real(dp),                 save :: gam(max_gfnff_ele)              ! Gamma for EEM
  real(dp),                 save :: rad(max_gfnff_ele)              ! Covalent radii
  real(dp),                 save :: rcov(max_gfnff_ele)             ! Covalent radii from D3
  real(dp),                 save :: repa(max_gfnff_ele)             ! Repa
  real(dp),                 save :: repan(max_gfnff_ele)            ! Repan
  real(dp),                 save :: repz(max_gfnff_ele)             ! Repz
  real(dp),                 save :: r0_radij(max_gfnff_ele)         ! R0 for radij
  real(dp),                 save :: sqrtZr4r2(max_gfnff_ele)        ! sqrtZr4r2
  real(dp),                 save :: tors(max_gfnff_ele)             ! Tors
  real(dp),                 save :: tors2(max_gfnff_ele)            ! Tors2
  real(dp),                 save :: zb3atm(max_gfnff_ele)           ! zb3atm
!
  real(dp),                 save :: d3r0(max_gfnff_ele*(max_gfnff_ele+1)/2)     ! BJ radii 
!
  real(dp),                 save :: xbaci(max_gfnff_ele)            ! xbaci
  real(dp),                 save :: xhaci(max_gfnff_ele)            ! xhaci
  real(dp),                 save :: xhbas(max_gfnff_ele)            ! xhbas
!
!  The following are from param type
!
  real(dp),                 save :: cnmax   = 4.4_dp       ! max. CN considered ie all larger values smoothly set to this val
  real(dp),                 save :: atcuta  = 0.595_dp     ! angle damping
  real(dp),                 save :: atcutt  = 0.505_dp     ! torsion angle damping
  real(dp),                 save :: atcuta_nci  = 0.395_dp ! nci angle damping in HB term
  real(dp),                 save :: atcutt_nci  = 0.305_dp ! nci torsion angle damping in HB term
  real(dp),                 save :: repscalb= 1.7583_dp    ! bonded rep. scaling
  real(dp),                 save :: repscaln= 0.4270_dp    ! non-bonded rep. scaling
  real(dp),                 save :: hbacut   =49.0_dp      ! HB angle cut-off
  real(dp),                 save :: hbscut   =22.0_dp      ! HB SR     "   "
  real(dp),                 save :: xbacut   =70.0_dp      ! same for XB
  real(dp),                 save :: xbscut   = 5.0_dp      !
  real(dp),                 save :: hbsf     = 1.0_dp      ! charge dep.
  real(dp),                 save :: hbst     =15.0_dp      ! 10 is better for S22, 20 better for HCN2 and S30L
  real(dp),                 save :: xbsf     =0.03_dp      !
  real(dp),                 save :: xbst     =15.0_dp      !
  real(dp),                 save :: hbalp    = 6.0_dp      ! damp
  real(dp),                 save :: hblongcut=85.0_dp      ! values larger than 85 yield large RMSDs for P26
  real(dp),                 save :: hblongcut_xb=70.0_dp   ! values larger than 70 yield large MAD for HAL28
  real(dp),                 save :: hbabmix  =0.80_dp      !
  real(dp),                 save :: hbnbcut  =11.20_dp     !
  real(dp),                 save :: tors_hb   =0.94_dp     ! torsion potential shift in HB term
  real(dp),                 save :: bend_hb   =0.20_dp     ! bending potential shift in HB term
  real(dp),                 save :: vbond_scale=0.9_dp     ! vbond(2) scaling for CN(H) = 1
  real(dp),                 save :: xhaci_globabh=0.268_dp ! A-H...B gen. scaling
  real(dp),                 save :: xhaci_coh=0.350_dp     ! A-H...O=C gen. scaling
  real(dp),                 save :: xhaci_glob=1.50_dp     ! acidity
!
!  The following are from gen type
!
  integer(i4),              save :: maxhiter=10           ! the Hückel iterations can diverge so take only a few steps
  real(dp),                 save :: linthr  = 160.0_dp    ! when is an angle close to linear ? (GEODEP) for metals values closer to 170 (than to 160) are better
                                                          ! but this occurs e.g. for Sc in unclear situations. So make it save (160)
  real(dp),                 save :: fcthr   = 1.d-3       ! skip torsion and bending if potential is small
  real(dp),                 save :: rthr     =1.25_dp     ! important bond determination threshold
                                                          ! large values yield more 1.23
  real(dp),                 save :: rthr2    =1.00_dp     ! decrease if a metal is present, larger values yield smaller CN
  real(dp),                 save :: rqshrink =0.23_dp     ! change of R0 for topo with charge qa, larger values yield smaller CN for metals in particular
  real(dp),                 save :: hqabthr  =0.01_dp     ! H charge (qa) threshold for H in HB list 18
  real(dp),                 save :: qabthr  =0.10_dp      ! AB charge (qa) threshold for AB in HB list, avoids HBs with positive atoms,
                                                          ! larger val. better for S30L but worse in PubChem RMSD checks
  real(dp),                 save :: srb1    = 0.3731_dp   ! bond params
  real(dp),                 save :: srb2    = 0.3171_dp   !
  real(dp),                 save :: srb3    = 0.2538_dp   !
  real(dp),                 save :: qrepscal= 0.3480_dp   ! change of non-bonded rep. with q(topo)
  real(dp),                 save :: nrepscal=-0.1270_dp   !   "    "      "       "   CN
  real(dp),                 save :: hhfac   = 0.6290_dp   ! HH repulsion
  real(dp),                 save :: hh13rep = 1.4580_dp   !
  real(dp),                 save :: hh14rep = 0.7080_dp   !
  real(dp),                 save :: bstren(9)             ! Bond strength 
  real(dp),                 save :: qfacBEN =-0.54_dp     ! bend FC change with polarity
  real(dp),                 save :: qfacTOR =12.0_dp      ! torsion FC change with polarity
  real(dp),                 save :: fr3     =0.3_dp       ! tors FC 3-ring
  real(dp),                 save :: fr4     =1.0_dp       ! tors FC 4-ring
  real(dp),                 save :: fr5     =1.5_dp       ! tors FC 5-ring
  real(dp),                 save :: fr6     =5.7_dp       ! tors FC 6-ring
  real(dp),                 save :: torsf(8)              ! torsion 
  real(dp),                 save :: fbs1    =0.50_dp      ! small bend corr.
  real(dp),                 save :: batmscal=0.30_dp      ! bonded ATM scal
  real(dp),                 save :: mchishift=-0.09_dp
  real(dp),                 save :: rabshift    =-0.110_dp! gen shift
  real(dp),                 save :: rabshifth   =-0.050_dp! XH
  real(dp),                 save :: hyper_shift = 0.03_dp ! hypervalent
  real(dp),                 save :: hshift3     = -0.11_dp! heavy
  real(dp),                 save :: hshift4     = -0.11_dp!
  real(dp),                 save :: hshift5     = -0.06_dp!
  real(dp),                 save :: metal1_shift= 0.2_dp  ! group 1+2 metals
  real(dp),                 save :: metal2_shift= 0.15_dp ! TM
  real(dp),                 save :: metal3_shift= 0.05_dp ! main group metals
  real(dp),                 save :: eta_shift   = 0.040_dp! eta bonded
  real(dp),                 save :: qfacbm(0:4)           ! bond charge dep.gff_srcs += 'gff/gfnff_input.f90'
  real(dp),                 save :: qfacbm0  = 0.047_dp   !
  real(dp),                 save :: rfgoed1  = 1.175_dp   ! topo dist scaling
  real(dp),                 save :: htriple  = 1.45_dp    ! decrease Hueckel off-diag for triple bonds because they are less well conjugated 1.4
  real(dp),                 save :: hueckelp2= 1.00_dp    ! increase pot depth depending on P
  real(dp),                 save :: hueckelp3=-0.24_dp    ! diagonal element change with qa
  real(dp),                 save :: hdiag(17)             ! diagonal element relative to C
  real(dp),                 save :: hoffdiag(17)          ! Huckel off-diag constants
  real(dp),                 save :: hiter   =0.700_dp     ! iteration mixing
  real(dp),                 save :: hueckelp=0.340_dp     ! diagonal qa dep.
  real(dp),                 save :: bzref   =0.370_dp     ! ref P value R shift
  real(dp),                 save :: bzref2  =0.315_dp     !  "  "  "    k stretch
  real(dp),                 save :: pilpf   =0.530_dp     ! 2el diag shift
  real(dp),                 save :: d3a1    = 0.58_dp     ! D3, s8 fixed = 2
  real(dp),                 save :: d3a2    = 4.80_dp
  real(dp),                 save :: split0  =0.670_dp     ! mixing of sp^n with sp^n-1
  real(dp),                 save :: split1  =0.330_dp     ! 1 - split0
  real(dp),                 save :: fringbo =0.020_dp     ! str ring size dep.
  real(dp),                 save :: aheavy3 =89.0_dp      ! three coord. heavy eq. angle
  real(dp),                 save :: aheavy4 =100.0_dp     ! four   "       "    "    "
  real(dp),                 save :: bsmat(0:3,0:3)        !
!
  real(dp),                 save :: tworootpi             ! Two divided by square root of pi
  real(dp),                 save :: gfnff_wolf_self       ! Wolf self term for topology
  real(dp),                 save :: gfnff_wolf_eta = 0.2_dp ! Eta value for Wolf sum in topology
  logical,                  save :: lgfnff_topowolf = .false. ! Use Wolf sum during topology
!
  character(len=2),         save :: atsym(max_gfnff_ele)
!
!  Atomic symbols
!
  data atsym/'H ','He','Li','Be','B ','C ','N ','O ','F ', &
             'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar', &
             'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co', &
             'Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
             'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh', &
             'Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &
             'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu', &
             'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
             'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl', &
             'Pb','Bi','Po','At','Rn'/
!
!  The following are parameters set as defaults by XTB
!
  data chi/1.227054_dp, 1.451412_dp, 0.813363_dp, 1.062841_dp, 1.186499_dp, &
           1.311555_dp, 1.528485_dp, 1.691201_dp, 1.456784_dp, 1.231037_dp, &
           0.772989_dp, 1.199092_dp, 1.221576_dp, 1.245964_dp, 1.248942_dp, &
           1.301708_dp, 1.312474_dp, 1.247701_dp, 0.781237_dp, 0.940834_dp, &
           0.950000_dp, 0.974455_dp, 0.998911_dp, 1.023366_dp, 1.047822_dp, &
           1.072277_dp, 1.096733_dp, 1.121188_dp, 1.145644_dp, 1.170099_dp, &
           1.205357_dp, 1.145447_dp, 1.169499_dp, 1.253293_dp, 1.329909_dp, &
           1.116527_dp, 0.950975_dp, 0.964592_dp, 0.897786_dp, 0.932824_dp, &
           0.967863_dp, 1.002901_dp, 1.037940_dp, 1.072978_dp, 1.108017_dp, &
           1.143055_dp, 1.178094_dp, 1.213132_dp, 1.205076_dp, 1.075529_dp, &
           1.206919_dp, 1.303658_dp, 1.332656_dp, 1.179317_dp, 0.789115_dp, &
           0.798704_dp, 1.127797_dp, 1.127863_dp, 1.127928_dp, 1.127994_dp, &
           1.128059_dp, 1.128125_dp, 1.128190_dp, 1.128256_dp, 1.128322_dp, &
           1.128387_dp, 1.128453_dp, 1.128518_dp, 1.128584_dp, 1.128649_dp, &
           1.128715_dp, 1.128780_dp, 1.129764_dp, 1.130747_dp, 1.131731_dp, &
           1.132714_dp, 1.133698_dp, 1.134681_dp, 1.135665_dp, 1.136648_dp, &
           1.061832_dp, 1.053084_dp, 1.207830_dp, 1.236314_dp, 1.310129_dp, &
           1.157380_dp/
!
  data gam/-0.448428_dp, 0.131022_dp, 0.571431_dp, 0.334622_dp,-0.089208_dp, &
           -0.025895_dp,-0.027280_dp,-0.031236_dp,-0.159892_dp, 0.074198_dp, &
            0.316829_dp, 0.326072_dp, 0.069748_dp,-0.120184_dp,-0.193159_dp, &
           -0.182428_dp,-0.064093_dp, 0.061914_dp, 0.318112_dp, 0.189248_dp, &
           -0.104172_dp,-0.082038_dp,-0.059903_dp,-0.037769_dp,-0.015635_dp, &
            0.006500_dp, 0.028634_dp, 0.050768_dp, 0.072903_dp, 0.095037_dp, &
            0.131140_dp, 0.097006_dp,-0.065744_dp,-0.058394_dp, 0.063307_dp, &
            0.091652_dp, 0.386337_dp, 0.530677_dp,-0.030705_dp,-0.020787_dp, &
           -0.010869_dp,-0.000951_dp, 0.008967_dp, 0.018884_dp, 0.028802_dp, &
            0.038720_dp, 0.048638_dp, 0.058556_dp, 0.036488_dp, 0.077711_dp, &
            0.077025_dp, 0.004547_dp, 0.039909_dp, 0.082630_dp, 0.485375_dp, &
            0.416264_dp,-0.011212_dp,-0.011046_dp,-0.010879_dp,-0.010713_dp, &
           -0.010546_dp,-0.010380_dp,-0.010214_dp,-0.010047_dp,-0.009881_dp, &
           -0.009714_dp,-0.009548_dp,-0.009382_dp,-0.009215_dp,-0.009049_dp, &
           -0.008883_dp,-0.008716_dp,-0.006220_dp,-0.003724_dp,-0.001228_dp, &
            0.001267_dp, 0.003763_dp, 0.006259_dp, 0.008755_dp, 0.011251_dp, &
            0.020477_dp,-0.056566_dp, 0.051943_dp, 0.076708_dp, 0.000273_dp, &
           -0.068929_dp/
!
  data cnf_gfnff/0.008904_dp, 0.004641_dp, 0.048324_dp, 0.080316_dp,-0.051990_dp, &
           0.031779_dp, 0.132184_dp, 0.157353_dp, 0.064120_dp, 0.036540_dp, &
          -0.000627_dp, 0.005412_dp, 0.018809_dp, 0.016329_dp, 0.012149_dp, &
           0.021484_dp, 0.014212_dp, 0.014939_dp, 0.003597_dp, 0.032921_dp, &
          -0.021804_dp,-0.022797_dp,-0.023789_dp,-0.024782_dp,-0.025775_dp, &
          -0.026767_dp,-0.027760_dp,-0.028753_dp,-0.029745_dp,-0.030738_dp, &
          -0.004189_dp,-0.011113_dp,-0.021305_dp,-0.012311_dp, 0.049781_dp, &
          -0.040533_dp, 0.012872_dp, 0.021056_dp,-0.003395_dp, 0.000799_dp, &
           0.004992_dp, 0.009186_dp, 0.013379_dp, 0.017573_dp, 0.021766_dp, &
           0.025960_dp, 0.030153_dp, 0.034347_dp,-0.000052_dp,-0.039776_dp, &
           0.006661_dp, 0.050424_dp, 0.068985_dp, 0.023470_dp,-0.024950_dp, &
          -0.033006_dp, 0.058973_dp, 0.058595_dp, 0.058217_dp, 0.057838_dp, &
           0.057460_dp, 0.057082_dp, 0.056704_dp, 0.056326_dp, 0.055948_dp, &
           0.055569_dp, 0.055191_dp, 0.054813_dp, 0.054435_dp, 0.054057_dp, &
           0.053679_dp, 0.053300_dp, 0.047628_dp, 0.041955_dp, 0.036282_dp, &
           0.030610_dp, 0.024937_dp, 0.019264_dp, 0.013592_dp, 0.007919_dp, &
           0.006383_dp,-0.089155_dp,-0.001293_dp, 0.019269_dp, 0.074803_dp, &
           0.016657_dp/
!
  data alp/0.585069_dp, 0.432382_dp, 0.628636_dp, 0.743646_dp, 1.167323_dp, &
           0.903430_dp, 1.278388_dp, 0.905347_dp, 1.067014_dp, 2.941513_dp, &
           0.687680_dp, 0.792170_dp, 1.337040_dp, 1.251409_dp, 1.068295_dp, &
           1.186476_dp, 1.593532_dp, 2.056749_dp, 0.674196_dp, 0.868052_dp, &
           0.575052_dp, 0.613424_dp, 0.651796_dp, 0.690169_dp, 0.728541_dp, &
           0.766913_dp, 0.805285_dp, 0.843658_dp, 0.882030_dp, 0.920402_dp, &
           0.877178_dp, 1.422350_dp, 1.405901_dp, 1.646860_dp, 2.001970_dp, &
           2.301695_dp, 1.020617_dp, 0.634141_dp, 0.652752_dp, 0.668845_dp, &
           0.684938_dp, 0.701032_dp, 0.717125_dp, 0.733218_dp, 0.749311_dp, &
           0.765405_dp, 0.781498_dp, 0.797591_dp, 1.296844_dp, 1.534068_dp, &
           1.727781_dp, 1.926871_dp, 2.175548_dp, 2.177702_dp, 0.977079_dp, &
           0.770260_dp, 0.757372_dp, 0.757352_dp, 0.757332_dp, 0.757313_dp, &
           0.757293_dp, 0.757273_dp, 0.757253_dp, 0.757233_dp, 0.757213_dp, &
           0.757194_dp, 0.757174_dp, 0.757154_dp, 0.757134_dp, 0.757114_dp, &
           0.757095_dp, 0.757075_dp, 0.756778_dp, 0.756480_dp, 0.756183_dp, &
           0.755886_dp, 0.755589_dp, 0.755291_dp, 0.754994_dp, 0.754697_dp, &
           0.868029_dp, 1.684375_dp, 2.001040_dp, 2.067331_dp, 2.228923_dp, &
           1.874218_dp/
!
  data bond/0.417997_dp, 0.258490_dp, 0.113608_dp, 0.195935_dp, 0.231217_dp, &
            0.385248_dp, 0.379257_dp, 0.339249_dp, 0.330706_dp, 0.120319_dp, &
            0.127255_dp, 0.173647_dp, 0.183796_dp, 0.273055_dp, 0.249044_dp, &
            0.290653_dp, 0.218744_dp, 0.034706_dp, 0.136353_dp, 0.192467_dp, &
            0.335860_dp, 0.314452_dp, 0.293044_dp, 0.271636_dp, 0.250228_dp, &
            0.228819_dp, 0.207411_dp, 0.186003_dp, 0.164595_dp, 0.143187_dp, &
            0.212434_dp, 0.210451_dp, 0.219870_dp, 0.224618_dp, 0.272206_dp, &
            0.147864_dp, 0.150000_dp, 0.150000_dp, 0.329501_dp, 0.309632_dp, &
            0.289763_dp, 0.269894_dp, 0.250025_dp, 0.230155_dp, 0.210286_dp, &
            0.190417_dp, 0.170548_dp, 0.150679_dp, 0.192977_dp, 0.173411_dp, &
            0.186907_dp, 0.192891_dp, 0.223202_dp, 0.172577_dp, 0.150000_dp, &
            0.150000_dp, 0.370682_dp, 0.368511_dp, 0.366339_dp, 0.364168_dp, &
            0.361996_dp, 0.359825_dp, 0.357654_dp, 0.355482_dp, 0.353311_dp, &
            0.351139_dp, 0.348968_dp, 0.346797_dp, 0.344625_dp, 0.342454_dp, &
            0.340282_dp, 0.338111_dp, 0.305540_dp, 0.272969_dp, 0.240398_dp, &
            0.207828_dp, 0.175257_dp, 0.142686_dp, 0.110115_dp, 0.077544_dp, &
            0.108597_dp, 0.148422_dp, 0.183731_dp, 0.192274_dp, 0.127706_dp, &
            0.086756_dp/
!
  data repa/2.639785_dp, 3.575012_dp, 0.732142_dp, 1.159621_dp, 1.561585_dp, &
            1.762895_dp, 2.173015_dp, 2.262269_dp, 2.511112_dp, 3.577220_dp, &
            0.338845_dp, 0.693023_dp, 0.678792_dp, 0.804784_dp, 1.012178_dp, &
            1.103469_dp, 1.209798_dp, 1.167791_dp, 0.326946_dp, 0.595242_dp, &
            1.447860_dp, 1.414501_dp, 1.381142_dp, 1.347783_dp, 1.314424_dp, &
            1.281065_dp, 1.247706_dp, 1.214347_dp, 1.180988_dp, 1.147629_dp, &
            0.700620_dp, 0.721266_dp, 0.741789_dp, 0.857434_dp, 0.875583_dp, &
            0.835876_dp, 0.290625_dp, 0.554446_dp, 0.623980_dp, 0.696005_dp, &
            0.768030_dp, 0.840055_dp, 0.912081_dp, 0.984106_dp, 1.056131_dp, &
            1.128156_dp, 1.200181_dp, 1.272206_dp, 0.478807_dp, 0.479759_dp, &
            0.579840_dp, 0.595241_dp, 0.644458_dp, 0.655289_dp, 0.574626_dp, &
            0.560506_dp, 0.682723_dp, 0.684824_dp, 0.686925_dp, 0.689026_dp, &
            0.691127_dp, 0.693228_dp, 0.695329_dp, 0.697430_dp, 0.699531_dp, &
            0.701631_dp, 0.703732_dp, 0.705833_dp, 0.707934_dp, 0.710035_dp, &
            0.712136_dp, 0.714237_dp, 0.745751_dp, 0.777265_dp, 0.808779_dp, &
            0.840294_dp, 0.871808_dp, 0.903322_dp, 0.934836_dp, 0.966350_dp, &
            0.467729_dp, 0.486102_dp, 0.559176_dp, 0.557520_dp, 0.563373_dp, &
            0.484713_dp/
!
  data repan/1.071395_dp, 1.072699_dp, 1.416847_dp, 1.156187_dp, 0.682382_dp, &
             0.556380_dp, 0.746785_dp, 0.847242_dp, 0.997252_dp, 0.873051_dp, &
             0.322503_dp, 0.415554_dp, 0.423946_dp, 0.415776_dp, 0.486773_dp, &
             0.494532_dp, 0.705274_dp, 0.706778_dp, 0.311178_dp, 0.399439_dp, &
             0.440983_dp, 0.475582_dp, 0.510180_dp, 0.544779_dp, 0.579377_dp, &
             0.613976_dp, 0.648574_dp, 0.683173_dp, 0.717772_dp, 0.752370_dp, &
             0.429944_dp, 0.420053_dp, 0.384743_dp, 0.443762_dp, 0.538680_dp, &
             0.472196_dp, 0.423850_dp, 0.385815_dp, 0.249213_dp, 0.285604_dp, &
             0.321995_dp, 0.358387_dp, 0.394778_dp, 0.431169_dp, 0.467560_dp, &
             0.503952_dp, 0.540343_dp, 0.576734_dp, 0.333476_dp, 0.348734_dp, &
             0.358194_dp, 0.351053_dp, 0.404536_dp, 0.389847_dp, 0.302575_dp, &
             0.163290_dp, 0.187645_dp, 0.190821_dp, 0.193998_dp, 0.197174_dp, &
             0.200351_dp, 0.203527_dp, 0.206703_dp, 0.209880_dp, 0.213056_dp, &
             0.216233_dp, 0.219409_dp, 0.222585_dp, 0.225762_dp, 0.228938_dp, &
             0.232115_dp, 0.235291_dp, 0.282937_dp, 0.330583_dp, 0.378229_dp, &
             0.425876_dp, 0.473522_dp, 0.521168_dp, 0.568814_dp, 0.616460_dp, &
             0.242521_dp, 0.293680_dp, 0.320931_dp, 0.322666_dp, 0.333641_dp, &
             0.434163_dp/
!
  data angl/1.661808_dp, 0.300000_dp, 0.018158_dp, 0.029224_dp, 0.572683_dp, &
            0.771055_dp, 1.053577_dp, 2.159889_dp, 1.525582_dp, 0.400000_dp, &
            0.041070_dp, 0.028889_dp, 0.086910_dp, 0.494456_dp, 0.409204_dp, &
            0.864972_dp, 1.986025_dp, 0.491537_dp, 0.050168_dp, 0.072745_dp, &
            0.378334_dp, 0.346400_dp, 0.314466_dp, 0.282532_dp, 0.250598_dp, &
            0.218663_dp, 0.186729_dp, 0.154795_dp, 0.122861_dp, 0.090927_dp, &
            0.140458_dp, 0.653971_dp, 0.528465_dp, 0.420379_dp, 2.243492_dp, &
            0.400000_dp, 0.035341_dp, 0.022704_dp, 0.195060_dp, 0.188476_dp, &
            0.181892_dp, 0.175308_dp, 0.168724_dp, 0.162139_dp, 0.155555_dp, &
            0.148971_dp, 0.142387_dp, 0.135803_dp, 0.169779_dp, 0.265730_dp, &
            0.505495_dp, 0.398254_dp, 2.640752_dp, 0.568026_dp, 0.032198_dp, &
            0.036663_dp, 0.281449_dp, 0.280526_dp, 0.279603_dp, 0.278680_dp, &
            0.277757_dp, 0.276834_dp, 0.275911_dp, 0.274988_dp, 0.274065_dp, &
            0.273142_dp, 0.272219_dp, 0.271296_dp, 0.270373_dp, 0.269450_dp, &
            0.268528_dp, 0.267605_dp, 0.253760_dp, 0.239916_dp, 0.226071_dp, &
            0.212227_dp, 0.198382_dp, 0.184538_dp, 0.170693_dp, 0.156849_dp, &
            0.104547_dp, 0.313474_dp, 0.220185_dp, 0.415042_dp, 1.259822_dp, &
            0.400000_dp/
!
  data angl2/0.624197_dp, 0.600000_dp, 0.050000_dp, 0.101579_dp, 0.180347_dp, &
             0.755851_dp, 0.761551_dp, 0.813653_dp, 0.791274_dp, 0.400000_dp, &
             0.000000_dp, 0.022706_dp, 0.100000_dp, 0.338514_dp, 0.453023_dp, &
             0.603722_dp, 1.051121_dp, 0.547904_dp, 0.000000_dp, 0.059059_dp, &
             0.117040_dp, 0.118438_dp, 0.119836_dp, 0.121234_dp, 0.122632_dp, &
             0.124031_dp, 0.125429_dp, 0.126827_dp, 0.128225_dp, 0.129623_dp, &
             0.206779_dp, 0.466678_dp, 0.496442_dp, 0.617321_dp, 0.409933_dp, &
             0.400000_dp, 0.000000_dp, 0.000000_dp, 0.119120_dp, 0.118163_dp, &
             0.117206_dp, 0.116249_dp, 0.115292_dp, 0.114336_dp, 0.113379_dp, &
             0.112422_dp, 0.111465_dp, 0.110508_dp, 0.149917_dp, 0.308383_dp, &
             0.527398_dp, 0.577885_dp, 0.320371_dp, 0.568026_dp, 0.000000_dp, &
             0.000000_dp, 0.078710_dp, 0.079266_dp, 0.079822_dp, 0.080379_dp, &
             0.080935_dp, 0.081491_dp, 0.082047_dp, 0.082603_dp, 0.083159_dp, &
             0.083716_dp, 0.084272_dp, 0.084828_dp, 0.085384_dp, 0.085940_dp, &
             0.086496_dp, 0.087053_dp, 0.095395_dp, 0.103738_dp, 0.112081_dp, &
             0.120423_dp, 0.128766_dp, 0.137109_dp, 0.145451_dp, 0.153794_dp, &
             0.323570_dp, 0.233450_dp, 0.268137_dp, 0.307481_dp, 0.316447_dp, &
             0.400000_dp/
!
  data tors/0.100000_dp, 0.100000_dp, 0.100000_dp, 0.000000_dp, 0.121170_dp, &
            0.260028_dp, 0.222546_dp, 0.250620_dp, 0.256328_dp, 0.400000_dp, &
            0.115000_dp, 0.000000_dp, 0.103731_dp, 0.069103_dp, 0.104280_dp, &
            0.226131_dp, 0.300000_dp, 0.400000_dp, 0.124098_dp, 0.000000_dp, &
            0.105007_dp, 0.107267_dp, 0.109526_dp, 0.111786_dp, 0.114046_dp, &
            0.116305_dp, 0.118565_dp, 0.120825_dp, 0.123084_dp, 0.125344_dp, &
            0.395722_dp, 0.349100_dp, 0.147808_dp, 0.259811_dp, 0.400000_dp, &
            0.400000_dp, 0.112206_dp,-0.004549_dp, 0.198713_dp, 0.179472_dp, &
            0.160232_dp, 0.140991_dp, 0.121751_dp, 0.102510_dp, 0.083270_dp, &
            0.064029_dp, 0.044789_dp, 0.025548_dp, 0.202245_dp, 0.278223_dp, &
            0.280596_dp, 0.229057_dp, 0.300000_dp, 0.423199_dp, 0.090741_dp, &
            0.076783_dp, 0.310896_dp, 0.309131_dp, 0.307367_dp, 0.305602_dp, &
            0.303838_dp, 0.302073_dp, 0.300309_dp, 0.298544_dp, 0.296779_dp, &
            0.295015_dp, 0.293250_dp, 0.291486_dp, 0.289721_dp, 0.287957_dp, &
            0.286192_dp, 0.284427_dp, 0.257959_dp, 0.231490_dp, 0.205022_dp, &
            0.178553_dp, 0.152085_dp, 0.125616_dp, 0.099147_dp, 0.072679_dp, &
            0.203077_dp, 0.169346_dp, 0.090568_dp, 0.144762_dp, 0.231884_dp, &
            0.400000_dp/
!
  data tors2/1.618678_dp, 1.000000_dp, 0.064677_dp, 0.000000_dp, 0.965814_dp, &
             1.324709_dp, 1.079334_dp, 1.478599_dp, 0.304844_dp, 0.500000_dp, &
             0.029210_dp, 0.000000_dp, 0.417423_dp, 0.334275_dp, 0.817008_dp, &
             0.922181_dp, 0.356367_dp, 0.684881_dp, 0.029210_dp, 0.000000_dp, &
             0.035902_dp, 0.090952_dp, 0.146002_dp, 0.201052_dp, 0.256103_dp, &
             0.311153_dp, 0.366203_dp, 0.421253_dp, 0.476303_dp, 0.531353_dp, &
             0.482963_dp, 1.415893_dp, 1.146581_dp, 1.338448_dp, 0.376801_dp, &
             0.500000_dp, 0.027213_dp,-0.004549_dp, 0.003820_dp, 0.093011_dp, &
             0.182202_dp, 0.271393_dp, 0.360584_dp, 0.449775_dp, 0.538965_dp, &
             0.628156_dp, 0.717347_dp, 0.806538_dp, 0.077000_dp, 0.185110_dp, &
             0.432427_dp, 0.887811_dp, 0.267721_dp, 0.571662_dp, 0.000000_dp, &
             0.000000_dp, 0.122336_dp, 0.131176_dp, 0.140015_dp, 0.148855_dp, &
             0.157695_dp, 0.166534_dp, 0.175374_dp, 0.184214_dp, 0.193053_dp, &
             0.201893_dp, 0.210733_dp, 0.219572_dp, 0.228412_dp, 0.237252_dp, &
             0.246091_dp, 0.254931_dp, 0.387526_dp, 0.520121_dp, 0.652716_dp, &
             0.785311_dp, 0.917906_dp, 1.050500_dp, 1.183095_dp, 1.315690_dp, &
             0.219729_dp, 0.344830_dp, 0.331862_dp, 0.767979_dp, 0.536799_dp, &
             0.500000_dp/
!
  data en/2.200,3.000,0.980,1.570,2.040,2.550,3.040,3.440,3.980, &
          4.500,0.930,1.310,1.610,1.900,2.190,2.580,3.160,3.500, &
          0.820,1.000,1.360,1.540,1.630,1.660,1.550,1.830,1.880, &
          1.910,1.900,1.650,1.810,2.010,2.180,2.550,2.960,3.000, &
          0.820,0.950,1.220,1.330,1.600,2.160,1.900,2.200,2.280, &
          2.200,1.930,1.690,1.780,1.960,2.050,2.100,2.660,2.600, &
          0.79,0.89,1.10,1.12,1.13,1.14,1.15,1.17,1.18,1.20,1.21,1.22, &
          1.23,1.24,1.25,1.26,1.27,1.3,1.5,1.7,1.9,2.1,2.2,2.2,2.2, &   ! value of W-Au modified
          2.00,1.62,2.33,2.02,2.0,2.2,2.2/
!
  data en_radij/2.30085633, 2.78445145, 1.52956084, 1.51714704, 2.20568300, &
                2.49640820, 2.81007174, 4.51078438, 4.67476223, 3.29383610, &
                2.84505365, 2.20047950, 2.31739628, 2.03636974, 1.97558064, &
                2.13446570, 2.91638164, 1.54098156, 2.91656301, 2.26312147, &
                2.25621439, 1.32628677, 2.27050569, 1.86790977, 2.44759456, &
                2.49480042, 2.91545568, 3.25897750, 2.68723778, 1.86132251, &
                2.01200832, 1.97030722, 1.95495427, 2.68920990, 2.84503857, &
                2.61591858, 2.64188286, 2.28442252, 1.33011187, 1.19809388, &
                1.89181390, 2.40186898, 1.89282464, 3.09963488, 2.50677823, &
                2.61196704, 2.09943450, 2.66930105, 1.78349472, 2.09634533, &
                2.00028974, 1.99869908, 2.59072029, 2.54497829, 2.52387890, &
                2.30204667, 1.60119300, 2.00000000, 2.00000000, 2.00000000, &
                2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000, &
                2.00000000, 2.00000000, 2.00000000, 2.00000000, 2.00000000, &
                2.00000000, 2.30089349, 1.75039077, 1.51785130, 2.62972945, &
                2.75372921, 2.62540906, 2.55860939, 3.32492356, 2.65140898, &
                1.52014458, 2.54984804, 1.72021963, 2.69303422, 1.81031095, &
                2.34224386/
!
  data r0_radij/0.55682207, 0.80966997, 2.49092101, 1.91705642, 1.35974851, &
                0.98310699, 0.98423007, 0.76716063, 1.06139799, 1.17736822, &
                2.85570926, 2.56149012, 2.31673425, 2.03181740, 1.82568535, &
                1.73685958, 1.97498207, 2.00136196, 3.58772537, 2.68096221, &
                2.23355957, 2.33135502, 2.15870365, 2.10522128, 2.16376162, &
                2.10804037, 1.96460045, 2.00476257, 2.22628712, 2.43846700, &
                2.39408483, 2.24245792, 2.05751204, 2.15427677, 2.27191920, &
                2.19722638, 3.80910350, 3.26020971, 2.99716916, 2.71707818, &
                2.34950167, 2.11644818, 2.47180659, 2.32198800, 2.32809515, &
                2.15244869, 2.55958313, 2.59141300, 2.62030465, 2.39935278, &
                2.56912355, 2.54374096, 2.56914830, 2.53680807, 4.24537037, &
                3.66542289, 3.19903011, 2.80000000, 2.80000000, 2.80000000, &
                2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000, &
                2.80000000, 2.80000000, 2.80000000, 2.80000000, 2.80000000, &
                2.80000000, 2.34880037, 2.37597108, 2.49067697, 2.14100577, &
                2.33473532, 2.19498900, 2.12678348, 2.34895048, 2.33422774, &
                2.86560827, 2.62488837, 2.88376127, 2.75174124, 2.83054552, &
                2.63264944/
!
  data cnfak_radij/0.17957827, 0.25584045,-0.02485871, 0.00374217, 0.05646607, &
                   0.10514203, 0.09753494, 0.30470380, 0.23261783, 0.36752208, &
                   0.00131819,-0.00368122,-0.01364510, 0.04265789, 0.07583916, &
                   0.08973207,-0.00589677, 0.13689929,-0.01861307, 0.11061699, &
                   0.10201137, 0.05426229, 0.06014681, 0.05667719, 0.02992924, &
                   0.03764312, 0.06140790, 0.08563465, 0.03707679, 0.03053526, &
                  -0.00843454, 0.01887497, 0.06876354, 0.01370795,-0.01129196, &
                   0.07226529, 0.01005367, 0.01541506, 0.05301365, 0.07066571, &
                   0.07637611, 0.07873977, 0.02997732, 0.04745400, 0.04582912, &
                   0.10557321, 0.02167468, 0.05463616, 0.05370913, 0.05985441, &
                   0.02793994, 0.02922983, 0.02220438, 0.03340460,-0.04110969, &
                  -0.01987240, 0.07260201, 0.07700000, 0.07700000, 0.07700000, &
                   0.07700000, 0.07700000, 0.07700000, 0.07700000, 0.07700000, &
                   0.07700000, 0.07700000, 0.07700000, 0.07700000, 0.07700000, &
                   0.07700000, 0.08379100, 0.07314553, 0.05318438, 0.06799334, &
                   0.04671159, 0.06758819, 0.09488437, 0.07556405, 0.13384502, &
                   0.03203572, 0.04235009, 0.03153769,-0.00152488, 0.02714675, &
                   0.04800662/
!
!  COVALENT RADII, used only in neighbor list determination
!  based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!  in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!  edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!  corrected Nov. 17, 2010 for the 92nd edition.
!
  data rad/0.32_dp,0.37_dp,1.30_dp,0.99_dp,0.84_dp,0.75_dp,0.71_dp,0.64_dp,0.60_dp,&
           0.62_dp,1.60_dp,1.40_dp,1.24_dp,1.14_dp,1.09_dp,1.04_dp,1.00_dp,1.01_dp,&
           2.00_dp,1.74_dp,1.59_dp,1.48_dp,1.44_dp,1.30_dp,1.29_dp,1.24_dp,1.18_dp,&
           1.17_dp,1.22_dp,1.20_dp,1.23_dp,1.20_dp,1.20_dp,1.18_dp,1.17_dp,1.16_dp,&
           2.15_dp,1.90_dp,1.76_dp,1.64_dp,1.56_dp,1.46_dp,1.38_dp,1.36_dp,1.34_dp,&
           1.30_dp,1.36_dp,1.40_dp,1.42_dp,1.40_dp,1.40_dp,1.37_dp,1.36_dp,1.36_dp,&
           2.38_dp,2.06_dp,1.94_dp,1.84_dp,1.90_dp,1.88_dp,1.86_dp,1.85_dp,1.83_dp,&
           1.82_dp,1.81_dp,1.80_dp,1.79_dp,1.77_dp,1.77_dp,1.78_dp,1.74_dp,1.64_dp,&
           1.58_dp,1.50_dp,1.41_dp,1.36_dp,1.32_dp,1.30_dp,1.30_dp,1.32_dp,1.44_dp,&
           1.45_dp,1.50_dp,1.42_dp,1.48_dp,1.46_dp/
!
!  Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197) 
!  Values for metals decreased by 10%.
!  NB: Values are corrected by 4/3 in pgfnff_init
!
   data rcov/ &
        0.32_dp,0.46_dp, & ! H,He
        1.20_dp,0.94_dp,0.77_dp,0.75_dp,0.71_dp,0.63_dp,0.64_dp,0.67_dp, & ! Li-Ne
        1.40_dp,1.25_dp,1.13_dp,1.04_dp,1.10_dp,1.02_dp,0.99_dp,0.96_dp, & ! Na-Ar
        1.76_dp,1.54_dp, & ! K,Ca
                        1.33_dp,1.22_dp,1.21_dp,1.10_dp,1.07_dp, & ! Sc-
                        1.04_dp,1.00_dp,0.99_dp,1.01_dp,1.09_dp, & ! -Zn
                        1.12_dp,1.09_dp,1.15_dp,1.10_dp,1.14_dp,1.17_dp, & ! Ga-Kr
        1.89_dp,1.67_dp, & ! Rb,Sr
                        1.47_dp,1.39_dp,1.32_dp,1.24_dp,1.15_dp, & ! Y-
                        1.13_dp,1.13_dp,1.08_dp,1.15_dp,1.23_dp, & ! -Cd
                        1.28_dp,1.26_dp,1.26_dp,1.23_dp,1.32_dp,1.31_dp, & ! In-Xe
        2.09_dp,1.76_dp, & ! Cs,Ba
                1.62_dp,1.47_dp,1.58_dp,1.57_dp,1.56_dp,1.55_dp,1.51_dp, & ! La-Eu
                1.52_dp,1.51_dp,1.50_dp,1.49_dp,1.49_dp,1.48_dp,1.53_dp, & ! Gd-Yb
                        1.46_dp,1.37_dp,1.31_dp,1.23_dp,1.18_dp, & ! Lu-
                        1.16_dp,1.11_dp,1.12_dp,1.13_dp,1.32_dp, & ! -Hg
                        1.30_dp,1.30_dp,1.36_dp,1.31_dp,1.38_dp,1.42_dp/ ! Tl-Rn
!
!  Atomic radii used to define molecular fragments
!
    data atomicRad / &
        0.32_dp, 0.37_dp, 1.30_dp, 0.99_dp, 0.84_dp, 0.75_dp, 0.71_dp, 0.64_dp, &
        0.60_dp, 0.62_dp, 1.60_dp, 1.40_dp, 1.24_dp, 1.14_dp, 1.09_dp, 1.04_dp, &
        1.00_dp, 1.01_dp, 2.00_dp, 1.74_dp, 1.59_dp, 1.48_dp, 1.44_dp, 1.30_dp, &
        1.29_dp, 1.24_dp, 1.18_dp, 1.17_dp, 1.22_dp, 1.20_dp, 1.23_dp, 1.20_dp, &
        1.20_dp, 1.18_dp, 1.17_dp, 1.16_dp, 2.15_dp, 1.90_dp, 1.76_dp, 1.64_dp, &
        1.56_dp, 1.46_dp, 1.38_dp, 1.36_dp, 1.34_dp, 1.30_dp, 1.36_dp, 1.40_dp, &
        1.42_dp, 1.40_dp, 1.40_dp, 1.37_dp, 1.36_dp, 1.36_dp, 2.38_dp, 2.06_dp, &
        1.94_dp, 1.84_dp, 1.90_dp, 1.88_dp, 1.86_dp, 1.85_dp, 1.83_dp, 1.82_dp, &
        1.81_dp, 1.80_dp, 1.79_dp, 1.77_dp, 1.77_dp, 1.78_dp, 1.74_dp, 1.64_dp, &
        1.58_dp, 1.50_dp, 1.41_dp, 1.36_dp, 1.32_dp, 1.30_dp, 1.30_dp, 1.32_dp, &
        1.44_dp, 1.45_dp, 1.50_dp, 1.42_dp, 1.48_dp, 1.46_dp/
!
!  PBE0/def2-QZVP atomic values calculated by S. Grimme in Gaussian (2010)
!  rare gases recalculated by J. Mewes with PBE0/aug-cc-pVQZ in Dirac (2018)
!  He: 3.4698 -> 3.5544, Ne: 3.1036 -> 3.7943, Ar: 5.6004 -> 5.6638,
!  Kr: 6.1971 -> 6.2312, Xe: 7.5152 -> 8.8367
!  not replaced but recalculated (PBE0/cc-pVQZ) were
!   H: 8.0589 ->10.9359, Li:29.0974 ->39.7226, Be:14.8517 ->17.7460
!  also new super heavies Cn,Nh,Fl,Lv,Og
!
  data r4Overr2/&
      8.0589_dp, 3.4698_dp, & ! H,He
     29.0974_dp,14.8517_dp,11.8799_dp, 7.8715_dp, 5.5588_dp, 4.7566_dp, 3.8025_dp, 3.1036_dp, & ! Li-Ne
     26.1552_dp,17.2304_dp,17.7210_dp,12.7442_dp, 9.5361_dp, 8.1652_dp, 6.7463_dp, 5.6004_dp, & ! Na-Ar
     29.2012_dp,22.3934_dp, & ! K,Ca
             19.0598_dp,16.8590_dp,15.4023_dp,12.5589_dp,13.4788_dp, & ! Sc-
             12.2309_dp,11.2809_dp,10.5569_dp,10.1428_dp, 9.4907_dp, & ! -Zn
                     13.4606_dp,10.8544_dp, 8.9386_dp, 8.1350_dp, 7.1251_dp, 6.1971_dp, & ! Ga-Kr
     30.0162_dp,24.4103_dp, & ! Rb,Sr
             20.3537_dp,17.4780_dp,13.5528_dp,11.8451_dp,11.0355_dp, & ! Y-
             10.1997_dp, 9.5414_dp, 9.0061_dp, 8.6417_dp, 8.9975_dp, & ! -Cd
                     14.0834_dp,11.8333_dp,10.0179_dp, 9.3844_dp, 8.4110_dp, 7.5152_dp, & ! In-Xe
     32.7622_dp,27.5708_dp, & ! Cs,Ba
             23.1671_dp,21.6003_dp,20.9615_dp,20.4562_dp,20.1010_dp,19.7475_dp,19.4828_dp, & ! La-Eu
             15.6013_dp,19.2362_dp,17.4717_dp,17.8321_dp,17.4237_dp,17.1954_dp,17.1631_dp, & ! Gd-Yb
             14.5716_dp,15.8758_dp,13.8989_dp,12.4834_dp,11.4421_dp, & ! Lu-
             10.2671_dp, 8.3549_dp, 7.8496_dp, 7.3278_dp, 7.4820_dp, & ! -Hg
                     13.5124_dp,11.6554_dp,10.0959_dp, 9.7340_dp, 8.8584_dp, 8.0125_dp/ ! Tl-Rn

  data metal/ &
   0,                                                                0,&!He
   1,1,                                               0, 0, 0, 0, 0, 0,&!Ne
   1,1,                                               1, 0, 0, 0, 0, 0,&!Ar
   1,1,2,                2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 0, 0, 0, 0, 0,&!Kr
   1,2,2,                2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 0, 0, 0, 0,&!Xe
   1,2,2, 2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
                         2, 2, 2, 2, 2, 2, 2, 2, 2,   1, 1, 1, 1, 0, 0/!Rn
  ! At is NOT a metal, Po is borderline but slightly better as metal
  data group/ &
   1,                                                                   8,&!He
   1,2,                                                  3, 4, 5, 6, 7, 8,&!Ne
   1,2,                                                  3, 4, 5, 6, 7, 8,&!Ar
   1,2,-3,                 -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Kr
   1,2,-3,                 -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8,&!Xe
   1,2,-3,  -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3, &
                           -4,-5,-6,-7,-8,-9,-10,-11,-12,3, 4, 5, 6, 7, 8/!Rn
  data normcn/ &  ! only for non metals well defined
   1,                                                                0,&!He
   4,4,                                               4, 4, 4, 2, 1, 0,&!Ne
   4,4,                                               4, 4, 4, 2, 1, 0,&!Ar
   4,4,4,                4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Kr
   4,4,4,                4, 6, 6, 6, 6, 6, 6, 4, 4,   4, 4, 4, 4, 1, 0,&!Xe
   4,4,4, 4,4,4,4,4,4,4,4,4,4,4,4,4,4, &
                         4, 6, 6, 6, 6, 6, 6, 6, 4,   4, 4, 4, 4, 1, 0/ !Rn
  data repz/ &
   1.,                                                                    2.,&!He
   1.,2.,                                                  3.,4.,5.,6.,7.,8.,&!Ne
   1.,2.,                                                  3.,4.,5.,6.,7.,8.,&!Ar
   1.,2.,3.,                 4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Kr
   1.,2.,3.,                 4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8.,&!Xe
   1.,2.,3.,  3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3., &
                             4.,5.,6.,7.,8.,9.,10.,11.,12.,3.,4.,5.,6.,7.,8./ !Rn

end module m_pgfnff
