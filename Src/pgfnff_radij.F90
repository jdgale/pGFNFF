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
  subroutine pgfnff_radij(nati,natj,cni,cnj,qi,qj,radij,lextra,dradijdcni,dradijdcnj,lgrad1)
!
!  Computes the coordination number dependent radii for a pair of atoms i and j for GFNFF
!
!  On entry : 
!
!  nati            = atomic number for i
!  natj            = atomic number for j
!  cni             = coordination number for i
!  cnj             = coordination number for j
!  qi              = charge for i
!  qj              = charge for j
!  lextra          = if true then compute extra contribution to radij
!  radij           = shift in radii for i and j (in Ang)
!  lgrad1          = if true then compute first derivatives
!
!  On exit :
!
!  radij           = sum of radii for i and j (in Ang)
!  dradijdcni      = derivative of radij w.r.t. cni
!  dradijdcnj      = derivative of radij w.r.t. cnj
!
!   8/20 Created
!  10/20 First derivatives added
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: nati
  integer(i4), intent(in)                          :: natj
  logical,     intent(in)                          :: lextra
  logical,     intent(in)                          :: lgrad1
  real(dp),    intent(in)                          :: cni
  real(dp),    intent(in)                          :: cnj
  real(dp),    intent(in)                          :: qi
  real(dp),    intent(in)                          :: qj
  real(dp),    intent(inout)                       :: radij
  real(dp),    intent(out)                         :: dradijdcni
  real(dp),    intent(out)                         :: dradijdcnj
!
!  Local variables
!
  integer(i4)                                      :: irow
  integer(i4)                                      :: jrow
  real(dp)                                         :: den
  real(dp)                                         :: f12
  real(dp)                                         :: f1
  real(dp)                                         :: f2
  real(dp)                                         :: ri
  real(dp)                                         :: rj
  real(dp),                                   save :: rowfct(6,2)
!
  data rowfct/29.84522887,-1.70549806, 6.54013762, 6.39169003, 6.00000000, 5.60000000, & ! First factor for each row
              -8.87843763, 2.10878369, 0.08009374,-0.85808076,-1.15000000,-1.30000000/   ! Second factor for each row
!
!  Set row of periodic table
!
  irow = PeriodicTableRow(nati)
  jrow = PeriodicTableRow(natj)
!
!  Compute standard radii with coordination corrections
!
  ri = r0_radij(nati) + cnfak_radij(nati)*cni
  rj = r0_radij(natj) + cnfak_radij(natj)*cnj
!
  den = abs(en_radij(nati)-en_radij(natj))
  f1 = 0.005_dp*(rowfct(irow,1) + rowfct(jrow,1))
  f2 = 0.005_dp*(rowfct(irow,2) + rowfct(jrow,2))
  f12 = 1.0_dp - f1*den - f2*den**2
  radij = (radij + ri + rj)*f12
!
  if (lgrad1) then
    dradijdcni = cnfak_radij(nati)*f12
    dradijdcnj = cnfak_radij(natj)*f12
  endif
!
  if (lextra) then
!
!  Now optionally apply corrections for charge/metallic character
!
    f1 = rqshrink
    f2 = rqshrink
    if (metal(nati).gt.0) f1 = 2.0_dp*f1
    if (metal(natj).gt.0) f2 = 2.0_dp*f2
    radij = radij - f1*qi - f2*qj
!
!  Finally scale factors for atoms
!
    radij = fat(nati)*fat(natj)*radij
!
    if (lgrad1) then
      dradijdcni = dradijdcni*fat(nati)*fat(natj)
      dradijdcnj = dradijdcnj*fat(nati)*fat(natj)
    endif
  endif
!
  return
  end
