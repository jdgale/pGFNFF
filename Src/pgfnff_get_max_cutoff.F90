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
  subroutine pgfnff_get_max_cutoff(numat,nat,cut2_nbr)
!
!  Finds the maximum cutoff that is used to generate the neighbour list
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat            ! Number of atoms
  integer(i4),                      intent(in)       :: nat(numat)       ! Atomic number of atoms
  real(dp),                         intent(out)      :: cut2_nbr         ! Maximum cutoff squared for the neighbour list
!
!  Local variables
!
  integer(i4)                                        :: i
  integer(i4)                                        :: j
  integer(i4)                                        :: nati
  integer(i4)                                        :: natj
  integer(i4)                                        :: ierror
  real(dp),    dimension(:),       allocatable, save :: cn
  real(dp)                                           :: cni
  real(dp)                                           :: cnj
  real(dp)                                           :: dr0dcni
  real(dp)                                           :: dr0dcnj
  real(dp)                                           :: r
  real(dp)                                           :: r0
!
!  Allocate local memory
!
  allocate(cn(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : cn '')')
    stop
  endif
!
!  Set coordination number to standard value for estimation of radii
!
  do i = 1,numat
    cn(i) = dble(normcn(nat(i)))
  enddo
!
!  Estimate largest possible cutoff that will be needed
!
  cut2_nbr = 0.0_dp
  cut2_nbr = max(cut2_nbr,cnthr)    ! Coordination number cutoff
!
  do i = 1,numat
    nati = nat(i)
    cni  = cn(i)
    do j = 1,i
      natj = nat(j)
      cnj  = cn(j)
!
!  Compute pairwise radius
!
      r0 = 0.0_dp
      call pgfnff_radij(nati,natj,cni,cnj,0.0_dp,0.0_dp,r0,.false.,dr0dcni,dr0dcnj,.false.)
!
      r = r0*rthr
      if (metal(nat(i)).eq.2) r = r*rthr2
      if (metal(nat(j)).eq.2) r = r*rthr2
      if (metal(nat(i)).eq.1) r = r*(rthr2 + 0.025_dp)
      if (metal(nat(j)).eq.1) r = r*(rthr2 + 0.025_dp)
      cut2_nbr = max(cut2_nbr,r*r) 
    enddo
  enddo
!
  deallocate(cn,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : cn '')')
    stop
  endif
!
  end subroutine pgfnff_get_max_cutoff
