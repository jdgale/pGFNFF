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
  subroutine pgfnff_check_atom_with_nobonds(numat,nat,nnbr,lbondsok)
!
!  Checks whether there are any atoms without bonds excluding group 8
!
  use m_pgfnff_types
  use m_pgfnff
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: numat            ! Number of atoms
  integer(i4),                      intent(in)       :: nat(numat)       ! Atomic number of atoms
  integer(i4),                      intent(in)       :: nnbr(numat)      ! Numbers of neighbours of atoms
  logical,                          intent(out)      :: lbondsok         ! If true then neighbour list is OK
!
!  Local variables
!
  integer(i4)                                        :: i
!
!  Check for atoms with no bonds that should have them so that we can try an increased cutoff
!
  lbondsok = .true.
  i = 0
  do while (lbondsok.and.i.lt.numat)
    i = i + 1
    if (nnbr(i).eq.0.and.group(nat(i)).ne.8) then
      lbondsok = .false.
    endif
  enddo

  end subroutine pgfnff_check_atom_with_nobonds
