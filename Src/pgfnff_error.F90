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
  subroutine pgfnff_error(errorstring,iline)
!
!  Outputs error message in standard form
!
!  Julian Gale, CIC, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  character(len=*), intent(in) :: errorstring
  integer(i4),      intent(in) :: iline
!******************
!  Output header  *
!******************
  write(ioout,'(/)')
  write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
  write(ioout,'(''!! ERROR : '',a)') trim(errorstring)
  if (iline.gt.0.and.iline.lt.1000000) then
    write(ioout,'(''!!       : Error is apparently on line '',i6)') iline
  endif
  write(ioout,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'',/)')
!
  return
  end
