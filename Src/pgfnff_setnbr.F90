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
  subroutine pgfnff_setnbr(numat,nat,cn,mchar,qf,maxnbr,nnbr,nbrno,rnbr)
!
!  Selects the neighbours needed for each of three neighbour lists
!  from the total neighbour list
!
!  Julian Gale, CIC, Curtin University, March 2021
!
  use m_pgfnff_types
  use m_pgfnff
  use m_pgfnff_nbr_lib
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat                      ! Number of atoms
  integer(i4),                      intent(in)       :: nat(numat)                 ! Atomic number of atoms
  integer(i4),                      intent(in)       :: maxnbr                     ! Maximum number of neighbours in list
  integer(i4),                      intent(in)       :: nnbr(numat)                ! Number of neighbours in list for each atom for each atom
  integer(i4),                      intent(in)       :: nbrno(maxnbr,numat)        ! Pointers to neighbours in list for each atom for each atom
  real(dp),                         intent(in)       :: cn(numat)                  ! Coordination number of atoms
  real(dp),                         intent(in)       :: mchar(numat)               ! Metallic character of atoms
  real(dp),                         intent(in)       :: qf(numat)                  ! Charges of atoms
  real(dp),                         intent(in)       :: rnbr(maxnbr,numat)         ! Distances to neighbours in full list
!
!  Local variables
!
  integer(i4)                                        :: hc_crit
  integer(i4)                                        :: i
  integer(i4)                                        :: j
  integer(i4)                                        :: nati
  integer(i4)                                        :: natj
  integer(i4)                                        :: ni
  logical                                            :: linclude_nohc
  logical                                            :: linclude_nome
  real(dp)                                           :: cni
  real(dp)                                           :: cnj
  real(dp)                                           :: dr0dcni
  real(dp)                                           :: dr0dcnj
  real(dp)                                           :: qi
  real(dp)                                           :: qj
  real(dp)                                           :: radij
  real(dp),          dimension(:), allocatable, save :: rads
  real(dp)                                           :: rij
  real(dp)                                           :: scale
!
!  Initialise specific neighbour lists
!
  nnbr_full(1:numat) = 0
  nnbr_nohc(1:numat) = 0
  nnbr_nome(1:numat) = 0
!
  allocate(rads(maxnbr))
!
  do i = 1,numat
    nati = nat(i)
    cni  = cn(i)
    qi   = qf(i)
!
!  Trap overly negative qi since this will lead to excessive cutoffs
!
    if (qi.lt.gfnff_q_trap) then
      call pgfnff_error('topological charge is below gfnff_qtrap limit',0_i4)
      stop
    endif
!
    do ni = 1,nnbr(i)
      j = nbrno(ni,i)
      natj = nat(j)
!
      cnj  = cn(j)
      qj   = qf(j)
!
!  Compute pairwise radius
!
      radij = 0.0_dp
      call pgfnff_radij(nati,natj,cni,cnj,qi,qj,radij,.true.,dr0dcni,dr0dcnj,.false.)
!
!  Scale by rthr
!
      radij = radij*rthr
      rads(ni) = radij
!
      scale = 1.0_dp
!
!  Full case
!
      if (metal(nati).eq.1) then
        scale = scale*(rthr2 + 0.025_dp)
      elseif (metal(nati).eq.2) then
        scale = scale*rthr2
      endif
      if (metal(natj).eq.1) then
        scale = scale*(rthr2 + 0.025_dp)
      elseif (metal(natj).eq.2) then
        scale = scale*rthr2
      endif
!
      rij = rnbr(ni,i)
      if (rij.lt.radij*scale) then
        nnbr_full(i) = nnbr_full(i) + 1
        nbrno_full(nnbr_full(i),i) = ni
      endif
    enddo
    do ni = 1,nnbr(i)
      j = nbrno(ni,i)
      natj = nat(j)
!
!  No highly coordinated atoms case
!
      linclude_nohc = .true.
      hc_crit = 6
      if (group(nati).le.2) hc_crit = 4
      if (nnbr_full(i).gt.hc_crit) linclude_nohc = .false.
!
      hc_crit = 6
      if (group(natj).le.2) hc_crit = 4
      if (nnbr_full(j).gt.hc_crit) linclude_nohc = .false.
!
!  No metals and unusually coordinated atoms case
!
      linclude_nome = .true.
      if (mchar(i).gt.0.25_dp.or.metal(nati).gt.0) linclude_nome = .false.
      if (mchar(j).gt.0.25_dp.or.metal(natj).gt.0) linclude_nome = .false.
      if (nnbr_full(i).gt.normcn(nati).and.nati.gt.10) linclude_nome = .false.
      if (nnbr_full(j).gt.normcn(natj).and.natj.gt.10) linclude_nome = .false.
!
      rij = rnbr(ni,i)
      if (rij.lt.rads(ni)) then
        if (linclude_nohc) then
          nnbr_nohc(i) = nnbr_nohc(i) + 1
          nbrno_nohc(nnbr_nohc(i),i) = ni
        endif
        if (linclude_nome) then
          nnbr_nome(i) = nnbr_nome(i) + 1
          nbrno_nome(nnbr_nome(i),i) = ni
        endif
      endif
    enddo
  enddo
!
  deallocate(rads)
!
  return
  end
