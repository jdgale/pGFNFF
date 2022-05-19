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
module m_pgfnff_topo
! 
!  This module contains the variables that hold the configuration specific
!  topological neighbour list information for pGFNFF in GULP
!
!  Julian Gale, Curtin University, March 2021
!
  use m_pgfnff_types
!
  implicit none
!
  integer(i4),                              save :: maxtopo = 24              ! Maximum number of neighbours
  integer(i4),                              save :: maxtoposhell  = 4         ! Maximum number of shells to search for topology
  integer(i4),                              save :: maxtoposhell1 = 2         ! Maximum number of shells to search for topology in the first loop
  logical,                                  save :: lgfnff_newtopo = .true.   ! This flag controls whether the topological charges use 1d+12 or 1+12.
  logical,                                  save :: lgfnff_xtbtopo = .false.  ! If true use XTB topology approach
  real(dp),                                 save :: cuttopo = 0.0_dp          ! Cutoff used to generate list
  real(dp),                                 save :: tdist_thr=12.0_dp         ! R threshold for covalent distance estimated used in apprx EEQ
                                                                              ! the following two parameters are critical for topo setup
!
!  Topological bonding lists
!
  integer(i4), dimension(:),       pointer, save :: nnbr_topo => null()      ! Number of bonded neighbours
  integer(i4), dimension(:,:),     pointer, save :: nbrno_topo => null()     ! Pointer to atom that is a bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: rtnbr => null()          ! Distance to topologically bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: xtnbr => null()          ! x component of distance to topologically bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: ytnbr => null()          ! y component of distance to topologically bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: ztnbr => null()          ! z component of distance to topologically bonded neighbour

CONTAINS

  subroutine changemaxtopo
!
!  Changes the size of arrays that hold the topological neighbour list for pGFNFF
!
!  Julian Gale, CIC, Curtin University, March 2021
!
  use m_pgfnff_cfg,      only : maxat => maxat_pgfnff
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Local variables
!
  integer(i4)               :: ierror
!
  call pgfnff_realloc(nnbr_topo,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_topo '')')
    stop
  endif
  call pgfnff_realloc(nbrno_topo,maxtopo,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_topo '')')
    stop
  endif
  call pgfnff_realloc(rtnbr,maxtopo,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : rtnbr '')')
    stop
  endif
  call pgfnff_realloc(xtnbr,maxtopo,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xtnbr '')')
    stop
  endif
  call pgfnff_realloc(ytnbr,maxtopo,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ytnbr '')')
    stop
  endif
  call pgfnff_realloc(ztnbr,maxtopo,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ztnbr '')')
    stop
  endif
!
  end subroutine changemaxtopo
!
end module m_pgfnff_topo
