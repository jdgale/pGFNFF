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
module m_pgfnff_cfg
! 
!  This module contains the variables that hold the configuration specific information for pGFNFF
!
!  Julian Gale, Curtin University, May 2021
!
  use m_pgfnff_types
!
  implicit none
!
  integer(i4),                              save :: maxat_pgfnff = 0     ! Maximum number of atoms for GFNFF arrays
!
!  EEM related arrays
!
  real(dp),                                 save :: radeeq               ! Radius for EEQ difference of erf(gam.r) from 1
  real(dp),    dimension(:),       pointer, save :: alpeeq => null()     ! Alpha
  real(dp),    dimension(:),       pointer, save :: chieeq => null()     ! Chi
  real(dp),    dimension(:),       pointer, save :: cnfeeq => null()     ! Coordination number factor
  real(dp),    dimension(:),       pointer, save :: gameeq => null()     ! Gamma
  real(dp),    dimension(:),       pointer, save :: qf0 => null()        ! Reference charges based on topology only
!
!  Dispersion related arrays that are atom specific
!
  real(dp),    dimension(:),       pointer, save :: d4_zeta => null()    ! Charge-dependent zeta value for C6 terms
!
!  Hydrogen bonding lists
!
  real(dp),                                 save :: rvhbref(3,3)         ! Reference cell for HB
  real(dp),    dimension(:),       pointer, save :: hbacid => null()     ! HB acid parameters
  real(dp),    dimension(:),       pointer, save :: hbbase => null()     ! HB base parameters
!
!  Hydrogen bond - pre final search data
!
  integer(i4),                              save :: maxathbH = 0
  integer(i4),                              save :: maxatxbAB = 0
  integer(i4),                              save :: nABat = 0            ! Number of atoms that could be A or B for hydrogen bond
  integer(i4),                              save :: nathbH = 0
  integer(i4),                              save :: natxbAB = 0
  integer(i4), dimension(:),       pointer, save :: nABatptr             ! Pointer to atoms that could be A or B for hydrogen bond
  integer(i4), dimension(:),       pointer, save :: hbatHl               ! Pointer to H atoms potentially in HBs
  integer(i4), dimension(:),       pointer, save :: rhbatHl              ! Reverse pointer to H atoms potentially in HBs
  integer(i4), dimension(:,:),     pointer, save :: xbatABl
  real(dp),    dimension(:),       pointer, save :: ABhbq                ! Charge-related factor for hydrogen bonds
  real(dp),    dimension(:,:),     pointer, save :: ABxbq                ! Charge-related factor for halogen bonds
!
!  Interaction parameters
!
  real(dp),    dimension(:),       pointer, save :: alphanb => null()  ! Non-bond parameter
  real(dp),    dimension(:),       pointer, save :: zetac6 => null()   ! Non-bond parameter
!
!  Bend parameters
!
  integer(i4),                              save :: maxangle = 0       ! Maximum number of angles
  integer(i4),                              save :: nangle = 0         ! Number of angles for bends
  integer(i4), dimension(:,:),     pointer, save :: nangleptr => null()! Pointer to atoms in angle
  real(dp),    dimension(:,:),     pointer, save :: vangle => null()   ! Bending parameters
!
!  Torsion parameters
!
  integer(i4),                              save :: maxtors = 0        ! Maximum number of torsion angles
  integer(i4),                              save :: ntors = 0          ! Number of torsion angles
  integer(i4), dimension(:,:),     pointer, save :: ntorsptr => null() ! Pointer to atoms in torsion
  real(dp),    dimension(:,:),     pointer, save :: vtors => null()    ! Torsion parameters
!
!  Huckel k point sampling
!
  integer(i4),                              save :: maxnkh = 1             ! Maximum number of k points for Huckel
  integer(i4),                              save :: nkh                    ! Number of k points for Huckel
  real(dp),    dimension(:),       pointer, save :: wkh          => null() ! Weight of k point for Huckel
  real(dp),    dimension(:),       pointer, save :: xkh          => null() ! X component of k point for Huckel
  real(dp),    dimension(:),       pointer, save :: ykh          => null() ! Y component of k point for Huckel
  real(dp),    dimension(:),       pointer, save :: zkh          => null() ! Z component of k point for Huckel

CONTAINS

  subroutine changemaxathbH
!
!  Changes the size of arrays that hold information for hydrogen bonds - nathbH
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(hbatHl,maxathbH,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : hbatHl '')')
    stop
  endif
!
  end subroutine changemaxathbH
!
  subroutine changemaxatxbAB
!
!  Changes the size of arrays that hold information for hydrogen bonds - natxbAB
!
!  Julian Gale, CIC, Curtin University, October 2020
!
  use m_io
  use m_pgfnff_reallocate
  implicit none
! 
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(xbatABl,3_i4,maxatxbAB,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xbatABl '')')
    stop
  endif
!
  end subroutine changemaxatxbAB
 
  subroutine changemaxangle
!
!  Changes the size of arrays that hold information for the angle bends
!
!  Julian Gale, CIC, Curtin University, September 2020
!
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(nangleptr,3_i4,maxangle,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nangleptr '')')
    stop
  endif
  call pgfnff_realloc(vangle,2_i4,maxangle,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : vangle '')')
    stop
  endif
!
  end subroutine changemaxangle

  subroutine changemaxtors
!
!  Changes the size of arrays that hold information for the torsions
!
!  Julian Gale, CIC, Curtin University, September 2020
!
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(ntorsptr,5_i4,maxtors,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ntorsptr '')')
    stop
  endif
  call pgfnff_realloc(vtors,2_i4,maxtors,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : vtors '')')
    stop
  endif
!
  end subroutine changemaxtors

  subroutine changemaxnkh
!
!  Changes the size of arrays that hold k points for Huckel solve
!
!  Julian Gale, CIC, Curtin University, April 2021
!
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(wkh,maxnkh,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : wkh '')')
    stop
  endif
  call pgfnff_realloc(xkh,maxnkh,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xkh '')')
    stop
  endif
  call pgfnff_realloc(ykh,maxnkh,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ykh '')')
    stop
  endif
  call pgfnff_realloc(zkh,maxnkh,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : zkh '')')
    stop
  endif
!
  end subroutine changemaxnkh

end module m_pgfnff_cfg
