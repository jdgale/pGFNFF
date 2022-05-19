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
  logical function pgfnff_excludebond(nat1,ntyp1,nat2,ntyp2,nnobo,nobond,nobotyp)
!
!  Check to see if a pair of atoms is in the nobond list
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)  :: nat1
  integer(i4),                     intent(in)  :: nat2
  integer(i4),                     intent(in)  :: ntyp1
  integer(i4),                     intent(in)  :: ntyp2
  integer(i4),                     intent(in)  :: nnobo
  integer(i4),                     intent(in)  :: nobond(nnobo)
  integer(i4),                     intent(in)  :: nobotyp(nnobo)
!
!  Local variables
!
  integer(i4)                                  :: ii
  integer(i4)                                  :: indb
  integer(i4)                                  :: nb1
  integer(i4)                                  :: nb2
  integer(i4)                                  :: nti
  integer(i4)                                  :: ntj
!
!  Initialise default value
!
  pgfnff_excludebond = .false.
!
!  If nnobo is zero then return
!
  if (nnobo.eq.0) return
!
!  Check whether bond type is excluded
!
  if (nat1.eq.nat2) then
    indb = nat2 + 1000*nat1
    if (ntyp1.lt.ntyp2) then
      nti = ntyp1
      ntj = ntyp2
    else
      nti = ntyp2
      ntj = ntyp1
    endif
  elseif (nat1.lt.nat2) then
    indb = nat2 + 1000*nat1
    nti = ntyp1
    ntj = ntyp2
  else
    indb = nat1 + 1000*nat2
    nti = ntyp2
    ntj = ntyp1
  endif
!
  ii = 1
  do while (.not.pgfnff_excludebond.and.(ii.le.nnobo))
    if (indb.eq.nobond(ii)) then
      nb1 = nobotyp(ii)/1000
      nb2 = nobotyp(ii) - 1000*nb1
      if ((nb1.eq.nti.or.nb1.eq.0).and.(nb2.eq.ntj.or.nb2.eq.0)) pgfnff_excludebond = .true.
      if (nat1.eq.nat2.and..not.pgfnff_excludebond) then
        if ((nb1.eq.ntj.or.nb1.eq.0).and.(nb2.eq.nti.or.nb2.eq.0)) pgfnff_excludebond = .true.
      endif
    endif
    ii = ii + 1
  enddo
!
  return
  end
