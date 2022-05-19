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
  subroutine pgfnff_dcn(i,nat,cut,logcn,sumdlogcn,maxnbr,nnbr,nbrno,rnbr,xnbr,ynbr,znbr)
!
!  Computes the coordination number for pGFNFF and the sum of derivatives
!
!  Julian Gale, CIC, Curtin University, August 2020
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                          :: i                   ! Atom whose coordination number is to be computed
  integer(i4), intent(in)                          :: maxnbr              ! Maximum number of neighbours
  integer(i4), intent(in)                          :: nat(*)              ! Atomic numbers
  integer(i4), intent(in)                          :: nnbr(*)             ! Number of neighbours
  integer(i4), intent(in)                          :: nbrno(maxnbr,*)     ! Pointers to neighbours
  real(dp),    intent(in)                          :: cut                 ! Cut-off 
  real(dp),    intent(out)                         :: logcn               ! Log of coordination number for i
  real(dp),    intent(out)                         :: sumdlogcn           ! Sum of log coordination number derivatives for i
  real(dp),    intent(in)                          :: rnbr(maxnbr,*)      ! Distances to neighbours
  real(dp),    intent(in)                          :: xnbr(maxnbr,*)      ! x component of distances to neighbours
  real(dp),    intent(in)                          :: ynbr(maxnbr,*)      ! y component of distances to neighbours
  real(dp),    intent(in)                          :: znbr(maxnbr,*)      ! z component of distances to neighbours
!
!  Local variables
!
  integer(i4)                                      :: j
  integer(i4)                                      :: ni
  real(dp)                                         :: cni
  real(dp)                                         :: dcnidrij
  real(dp)                                         :: dlogcnidcni
  real(dp)                                         :: dlogcni(3)
  real(dp)                                         :: dtrm
  real(dp)                                         :: dr
  real(dp)                                         :: erfCN
  real(dp)                                         :: exptrm
  real(dp)                                         :: logcnmax
  real(dp)                                         :: pi
  real(dp)                                         :: r0
  real(dp)                                         :: rij
  real(dp)                                         :: sumdlogcni
!
  pi = 4.0_dp*atan(1.0_dp)
!************************************************************
!  Loop over pairs of atoms to compute coordination number  *
!************************************************************
  logcnmax = log(1.0_dp + exp(cnmax))
  cni = 0.0_dp
!
!  Compute total coordination number for i, cni
!
  do ni = 1,nnbr(i)
    rij = rnbr(ni,i)
    if (rij.le.cut) then
      j = nbrno(ni,i)
      r0 = rcov(nat(i)) + rcov(nat(j))
      dr = (rij - r0)/r0
      erfCN = 0.5_dp*(1.0_dp + erf(gfnff_kn_cn*dr))
      cni = cni + erfCN
    endif
  enddo
!
!  Now create a function of this coordination number
!
  exptrm = exp(cnmax - cni)
  logcn = logcnmax - log(1.0_dp + exptrm)
!
  sumdlogcn = 0.0_dp
  sumdlogcni = 0.0_dp
  dlogcnidcni = exptrm/(1.0_dp + exptrm)
!
!  Compute derivatives of total coordination number for i
!
  dlogcni(1:3) = 0.0_dp
  do ni = 1,nnbr(i)
    rij = rnbr(ni,i)
    if (rij.le.cut) then
      j = nbrno(ni,i)
      r0 = rcov(nat(i)) + rcov(nat(j))
      dr = (rij - r0)/r0
      dcnidrij = gfnff_kn_cn*exp(-(gfnff_kn_cn*dr)**2)/(sqrt(pi)*r0)
      dtrm = dlogcnidcni*dcnidrij
!
!  i-j contribution
!
      sumdlogcn  = sumdlogcn  + abs(dtrm)
!
!  Self term
!
      dcnidrij = dcnidrij/rij
      dlogcni(1) = dlogcni(1) - dlogcnidcni*dcnidrij*xnbr(ni,i)
      dlogcni(2) = dlogcni(2) - dlogcnidcni*dcnidrij*ynbr(ni,i)
      dlogcni(3) = dlogcni(3) - dlogcnidcni*dcnidrij*znbr(ni,i)
    endif
  enddo
  sumdlogcn  = sumdlogcn  + sqrt(dlogcni(1)**2+dlogcni(2)**2+dlogcni(3)**2)
!
  return
  end
