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
module m_pgfnff_nbr_lib
! 
!  This module contains the variables that hold the configuration specific
!  neighbour list information for pGFNFF that is only needed during parameterisation
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
!
  implicit none
!
  integer(i4),                              save :: maxnbr_lib = 6           ! Maximum number of neighbours
!
!  Pointers for specific neighbour lists
!
  integer(i4), dimension(:),       pointer, save :: nnbr_full => null()      ! Number of neighbours - full
  integer(i4), dimension(:),       pointer, save :: nnbr_nohc => null()      ! Number of neighbours - no highly coordinated atoms
  integer(i4), dimension(:),       pointer, save :: nnbr_nome => null()      ! Number of neighbours - no metals or unusual coordination
  integer(i4), dimension(:,:),     pointer, save :: nbrno_full => null()     ! Pointer to nborno for atom that is neighbour - full
  integer(i4), dimension(:,:),     pointer, save :: nbrno_nohc => null()     ! Pointer to nborno for atom that is neighbour - no high coordination
  integer(i4), dimension(:,:),     pointer, save :: nbrno_nome => null()     ! Pointer to nborno for atom that is neighbour - no metals
!
!  Bonding lists
!
  integer(i4), dimension(:),       pointer, save :: nnbr_pi   => null()      ! Number of pi-bonded neighbours
  integer(i4), dimension(:,:),     pointer, save :: nbrno_pi   => null()     ! Pointer to atom that is a pi-bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: xpnbr => null()          ! x component of distance to pi-bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: ypnbr => null()          ! y component of distance to pi-bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: zpnbr => null()          ! z component of distance to pi-bonded neighbour
!
!  Parameters for bonds
!
  real(dp),    dimension(:,:,:),   pointer, save :: vbnbr => null()          ! Parameters for the interaction with a bonded neighbour

CONTAINS

  subroutine changemaxnbrlib
!
!  Changes the size of arrays that hold the neighbour list for pGFNFF
!
!  Julian Gale, CIC, Curtin University, August 2020
!
  use m_pgfnff_cfg,  only : maxat => maxat_pgfnff
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(nnbr_full,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_full '')')
    stop
  endif
  call pgfnff_realloc(nnbr_nohc,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_nohc '')')
    stop
  endif
  call pgfnff_realloc(nnbr_nome,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_nome '')')
    stop
  endif
  call pgfnff_realloc(nbrno_full,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_full '')')
    stop
  endif
  call pgfnff_realloc(nbrno_nohc,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_nohc '')')
    stop
  endif
  call pgfnff_realloc(nbrno_nome,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_nome '')')
    stop
  endif
!
  call pgfnff_realloc(nnbr_pi,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_pi '')')
    stop
  endif
  call pgfnff_realloc(nbrno_pi,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_pi '')')
    stop
  endif
  call pgfnff_realloc(xpnbr,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xpnbr '')')
    stop
  endif
  call pgfnff_realloc(ypnbr,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ypnbr '')')
    stop
  endif
  call pgfnff_realloc(zpnbr,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : zpnbr '')')
    stop
  endif
!
  call pgfnff_realloc(vbnbr,3_i4,maxnbr_lib,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : vbnbr '')')
    stop
  endif
!
  end subroutine changemaxnbrlib

end module m_pgfnff_nbr_lib
