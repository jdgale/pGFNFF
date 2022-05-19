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
module m_pgfnff_mrec
  private
  public :: pgfnff_mrecgff, pgfnff_mrecgff_pi

contains

  subroutine pgfnff_mrecgff_pi(npi,maxnbr,nnbr,nbrno,molcount,molvec,lcheckpbc,lpbc,lxref)
!
!  Main routine for searching for pi system fragments
!
  use m_pgfnff_types
  use m_pgfnff_nbr_lib
  use m_io
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)     :: npi                 ! Number of pi atoms
  integer(i4),                  intent(out)    :: molvec(npi)         ! Assignment vector of atom to fragment
  integer(i4),                  intent(in)     :: maxnbr              ! Maximum number of neighbours for atoms
  integer(i4),                  intent(in)     :: nnbr(*)             ! Number of neighbours of each atom
  integer(i4),                  intent(in)     :: nbrno(maxnbr,*)     ! Neighbour list of each atom
  integer(i4),                  intent(out)    :: molcount            ! Number of total fragments (increased during search)
  logical,                      intent(in)     :: lcheckpbc           ! If true then check for periodicity
  logical,                      intent(out)    :: lpbc(npi)           ! If true then pi system is periodic
  logical,                      intent(in)     :: lxref               ! If true then cross reference to nbrno
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ni
  integer(i4)                                  :: j
  integer(i4)                                  :: nj
  integer(i4)                                  :: k
  integer(i4)                                  :: ierror
  integer(i4),     dimension(:,:), allocatable :: bond
  integer(i4),     dimension(:,:), allocatable :: bondrptr
  logical,         dimension(:),   allocatable :: taken
  real(dp)                                     :: xyz(3)
  real(dp),        dimension(:,:), allocatable :: xyzpi
!
  allocate(taken(npi),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : taken '')')
    stop
  endif
  allocate(bond(maxnbr,npi),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : bond '')')
    stop
  endif
  allocate(bondrptr(maxnbr,npi),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : bondrptr '')')
    stop
  endif
  allocate(xyzpi(3_i4,npi),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xyzpi '')')
    stop
  endif
!
!  Create pointer for bond i->j in list of i for the same bond in list of j
!
  do i = 1,npi
    do ni = 1,nnbr_pi(i)
      if (lxref) then
        j = nbrno(nbrno_pi(ni,i),i)
      else
        j = nbrno_pi(ni,i)
      endif
!
!  Find i in bonding list of j
!
      bondrptr(ni,i) = 0
      do nj = 1,nnbr_pi(j)
        if (lxref) then
          k = nbrno(nbrno_pi(nj,j),j)
        else
          k = nbrno_pi(nj,j)
        endif
        if (k.eq.i) then
          bondrptr(ni,i) = nj
          exit
        endif
      enddo
!
!  Check that index was found
!
      if (bondrptr(ni,i).eq.0) then
        call pgfnff_error('bond is only present for i to j and not j to i',0_i4)
        stop
      endif
    enddo
  enddo
!
  bond = 0
  do i = 1,npi
    do ni = 1,nnbr_pi(i)
      bond(ni,i) = 1
    enddo
  enddo
  molvec(1:npi) = 0
  lpbc(1:npi) = .false.
  molcount = 1
  taken(1:npi) = .false.
  do i = 1,npi
    if (.not.taken(i)) then
      molvec(i) = molcount
      taken(i) = .true.
      if (lcheckpbc) then
        xyz(1:3)     = 0.0_dp     ! Initialise position of the first atom in the pi fragment
        xyzpi(1:3,i) = 0.0_dp     ! Initialise position of this atom in the pi fragment
      endif
      call mrecgff2_pi(npi,maxnbr,nnbr,nbrno,i,taken,bond,bondrptr,molvec,molcount, &
                       lcheckpbc,lpbc,lxref,xyz,xyzpi)
      molcount = molcount + 1
    endif
  enddo
!
  molcount = molcount - 1
!
  deallocate(xyzpi,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : xyzpi '')')
    stop
  endif
  deallocate(bondrptr,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : bondrptr '')')
    stop
  endif
  deallocate(bond,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : bond '')')
    stop
  endif
  deallocate(taken,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : taken '')')
    stop
  endif

  end subroutine pgfnff_mrecgff_pi

  recursive subroutine mrecgff2_pi(npi,maxnbr,nnbr,nbrno,i,taken,bond,bondrptr,molvec,molcnt, &
                                   lcheckpbc,lpbc,lxref,xyz,xyzpi)
!
  use m_pgfnff_types
  use m_pgfnff_nbr_lib
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)    :: i
  integer(i4),                  intent(in)    :: npi
  integer(i4),                  intent(in)    :: maxnbr
  integer(i4),                  intent(in)    :: nnbr(npi)
  integer(i4),                  intent(in)    :: nbrno(maxnbr,npi)
  integer(i4),                  intent(inout) :: molvec(npi)
  integer(i4),                  intent(in)    :: molcnt
  integer(i4),                  intent(inout) :: bond(maxnbr,npi)
  integer(i4),                  intent(in)    :: bondrptr(maxnbr,npi)
  logical,                      intent(in)    :: lcheckpbc
  logical,                      intent(inout) :: lpbc(npi)
  logical,                      intent(inout) :: taken(npi)
  logical,                      intent(in)    :: lxref
  real(dp),                     intent(inout) :: xyz(3_i4)
  real(dp),                     intent(inout) :: xyzpi(3_i4,npi)
!
!  Local variables
!
  integer(i4)                                 :: icn
  integer(i4)                                 :: j
  integer(i4)                                 :: jj
  integer(i4)                                 :: k
  real(dp)                                    :: diff
  real(dp)                                    :: xyzloc(3)
!
  icn = nnbr(i)
  do k = 1,icn
    jj = maxloc(bond(1:icn,i),1)
    if (lxref) then
      j = nbrno(nbrno_pi(jj,i),i)
    else
      j = nbrno_pi(jj,i)
    endif
    if (lcheckpbc) then
!
!  Set vector from centre of fragment to j
!
      xyzloc(1) = xyz(1) + xpnbr(jj,i)
      xyzloc(2) = xyz(2) + ypnbr(jj,i)
      xyzloc(3) = xyz(3) + zpnbr(jj,i)
    endif
!
    bond(jj,i) = 0
    bond(bondrptr(jj,i),j) = 0
!
    if (i.eq.j) cycle
!
    if (.not.taken(j)) then
      molvec(j) = molcnt
      taken(j) = .true.
      if (lcheckpbc) then
        xyzpi(1:3,j) = xyzloc(1:3)
      endif
      call mrecgff2_pi(npi,maxnbr,nnbr,nbrno,j,taken,bond,bondrptr,molvec,molcnt, &
                       lcheckpbc,lpbc,lxref,xyzloc,xyzpi)
    elseif (lcheckpbc) then
!
!  Check for periodicity
!
      diff = abs(xyzloc(1) - xyzpi(1,j)) + abs(xyzloc(2) - xyzpi(2,j)) + abs(xyzloc(3) - xyzpi(3,j))
      if (diff.gt.1.0d-2) lpbc(molcnt) = .true.
    endif
  enddo
  end subroutine mrecgff2_pi

  subroutine pgfnff_mrecgff(numat,maxnbr,nnbr,nbrno,nbrno2,molcount,molvec,lxref)
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)     :: numat               ! Number of atoms
  integer(i4),                  intent(out)    :: molvec(numat)       ! Assignment vector of atom to fragment
  integer(i4),                  intent(in)     :: maxnbr              ! Maximum number of neighbours for atoms
  integer(i4),                  intent(in)     :: nnbr(numat)         ! Number of neighbours of each atom
  integer(i4),                  intent(in)     :: nbrno(maxnbr,*)     ! Neighbour list of each atom
  integer(i4),                  intent(in)     :: nbrno2(maxnbr,*)    ! Second neighbour list of each atom
  integer(i4),                  intent(out)    :: molcount            ! Number of total fragments (increased during search)
  logical,                      intent(in)     :: lxref               ! If true then cross reference to nbrno
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ni
  integer(i4)                                  :: j
  integer(i4)                                  :: nj
  integer(i4)                                  :: k
  integer(i4)                                  :: ierror
  integer(i4),     dimension(:,:), allocatable :: bond
  integer(i4),     dimension(:,:), allocatable :: bondrptr
  logical,         dimension(:),   allocatable :: taken
!
  allocate(taken(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : taken '')')
    stop
  endif
  allocate(bond(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : bond '')')
    stop
  endif
  allocate(bondrptr(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : bondrptr '')')
    stop
  endif
!
!  Create pointer for bond i->j in list of i for the same bond in list of j
!
  do i = 1,numat
    do ni = 1,nnbr(i)
      if (lxref) then
        j = nbrno2(nbrno(ni,i),i)
      else
        j = nbrno(ni,i)
      endif
!
!  Find i in bonding list of j
!
      bondrptr(ni,i) = 0
      do nj = 1,nnbr(j)
        if (lxref) then
          k = nbrno2(nbrno(nj,j),j)
        else
          k = nbrno(nj,j)
        endif
        if (k.eq.i) then
          bondrptr(ni,i) = nj
          exit
        endif
      enddo
!
!  Check that index was found
!
      if (bondrptr(ni,i).eq.0) then
        call pgfnff_error('bond is only present for i to j and not j to i',0_i4)
        stop
      endif
    enddo
  enddo
!
  bond = 0
  do i = 1,numat
    do ni = 1,nnbr(i)
      bond(ni,i) = 1
    enddo
  enddo
  molvec(1:numat) = 0
  molcount = 1
  taken(1:numat) = .false.
  do i = 1,numat
    if (.not.taken(i)) then
      molvec(i) = molcount
      taken(i) = .true.
      call mrecgff2(numat,maxnbr,nnbr,nbrno,nbrno2,i,taken,bond,bondrptr,molvec,molcount,lxref)
      molcount = molcount + 1
    endif
  enddo
!
  molcount = molcount - 1
!
  deallocate(bondrptr,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : bondrptr '')')
    stop
  endif
  deallocate(bond,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : bond '')')
    stop
  endif
  deallocate(taken,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : taken '')')
    stop
  endif

  end subroutine pgfnff_mrecgff

  recursive subroutine mrecgff2(numat,maxnbr,nnbr,nbrno,nbrno2,i,taken,bond,bondrptr,molvec,molcnt,lxref)
!
  use m_pgfnff_types
  use m_pgfnff_nbr_lib
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)    :: numat
  integer(i4),                  intent(in)    :: nnbr(numat)
  integer(i4),                  intent(in)    :: maxnbr
  integer(i4),                  intent(in)    :: nbrno(maxnbr,numat)
  integer(i4),                  intent(in)    :: nbrno2(maxnbr,numat)
  integer(i4),                  intent(inout) :: molvec(numat)
  integer(i4),                  intent(in)    :: molcnt
  integer(i4),                  intent(inout) :: bond(maxnbr,numat)
  integer(i4),                  intent(in)    :: bondrptr(maxnbr,numat)
  logical,                      intent(inout) :: taken(numat)
  logical,                      intent(in)    :: lxref
!
!  Local variables
!
  integer(i4)                                 :: i
  integer(i4)                                 :: icn
  integer(i4)                                 :: j
  integer(i4)                                 :: jj
  integer(i4)                                 :: k
!
  icn = nnbr(i)
  do k = 1,icn
    jj = maxloc(bond(1:icn,i),1)
    if (lxref) then
      j = nbrno2(nbrno(jj,i),i)
    else
      j = nbrno(jj,i)
    endif
!
    bond(jj,i) = 0
    bond(bondrptr(jj,i),j) = 0
!
    if (i.eq.j) cycle
!
    if (.not.taken(j)) then
      molvec(j) = molcnt
      taken(j) = .true.
      call mrecgff2(numat,maxnbr,nnbr,nbrno,nbrno2,j,taken,bond,bondrptr,molvec,molcnt,lxref)
    endif
  enddo
  end subroutine mrecgff2

end module m_pgfnff_mrec
