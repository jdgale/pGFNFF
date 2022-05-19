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
!**************
!  I/O Module *
!**************
!
!  Defines the I/O channels
!
!  Julian Gale, CIC, Curtin University, August 2020
!
  module m_io
    use m_pgfnff_types
    integer(i4),                    save :: ioin  = 5_i4  ! Input channel
    integer(i4),                    save :: ioout = 6_i4  ! Output channel
  end module m_io
