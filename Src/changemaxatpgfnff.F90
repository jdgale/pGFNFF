!
! pGFNFF library:
! ---------------
!
! Copyright (C) 2021-2022 Julian Gale
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
  subroutine changemaxat_pgfnff
!
!  Changes the size of arrays that depend on the maximum number of atoms in pGFNFF arrays
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use m_pgfnff_reallocate
  use m_io
  use m_pgfnff_cfg
  implicit none
!
!  Local variables
!
  integer(i4)              :: ierror
  integer(i4)              :: maxpar
!
  maxpar = maxat_pgfnff*(maxat_pgfnff+1)/2
!
  call pgfnff_realloc(nABatptr,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nABatptr '')')
    stop
  endif
  call pgfnff_realloc(ABhbq,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ABhbq '')')
    stop
  endif
  call pgfnff_realloc(ABxbq,2_i4,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ABxbq '')')
    stop
  endif
  call pgfnff_realloc(d4_zeta,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : d4_zeta '')')
    stop
  endif
  call pgfnff_realloc(qf0,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : qf0 '')')
    stop
  endif
  call pgfnff_realloc(hbacid,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : hbacid '')')
    stop
  endif
  call pgfnff_realloc(hbbase,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : hbbase '')')
    stop
  endif
  call pgfnff_realloc(rhbatHl,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : rhbatHl '')')
    stop
  endif
!
  call pgfnff_realloc(alpeeq,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : alpeeq '')')
    stop
  endif
  call pgfnff_realloc(chieeq,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : chieeq '')')
    stop
  endif
  call pgfnff_realloc(cnfeeq,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : cnfeeq '')')
    stop
  endif
  call pgfnff_realloc(gameeq,maxat_pgfnff,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : gameeq '')')
    stop
  endif
!
!  Arrays that depend on maxpar
!
  call pgfnff_realloc(alphanb,maxpar,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : alphanb '')')
    stop
  endif
  call pgfnff_realloc(zetac6,maxpar,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : zetac6 '')')
    stop
  endif
!
  end subroutine changemaxat_pgfnff
