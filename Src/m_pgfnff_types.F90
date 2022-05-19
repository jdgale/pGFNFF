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
  module m_pgfnff_types
!
!  i4 and i2 define the size of integer*4 and integer*2 in the code
!
    use, intrinsic :: iso_c_binding
!
!  Integers
!
    integer, parameter :: i4  = selected_int_kind(5)
    integer, parameter :: i2  = selected_int_kind(3)
!
!  Floats
!
    integer, parameter :: dp  = kind(1.0d0)
    integer, parameter :: dpc = kind((1.0d0,1.0d0))

  end module m_pgfnff_types
