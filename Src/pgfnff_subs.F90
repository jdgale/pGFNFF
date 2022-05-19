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
! Modifications for pGFN-FF:
! --------------------------
!
! Julian Gale, Curtin University, 2021
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Condense charges to heavy atoms based on topology  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine qheavy(numat,nat,nnbr,maxnbr,nbrno,q)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in)    ::  numat
  integer(i4),  intent(in)    ::  maxnbr
  integer(i4),  intent(in)    ::  nnbr(numat)
  integer(i4),  intent(in)    ::  nbrno(maxnbr,numat)
  integer(i4),  intent(in)    ::  nat(numat)
  real(dp),     intent(inout) ::  q(numat)
!
!  Local variables
!
  integer(i4)                 :: i
  integer(i4)                 :: j
  integer(i4)                 :: ni
  real(dp)                    :: qtmp(numat)

  qtmp = q
  do i = 1,numat
    if (nat(i).ne.1) cycle
    qtmp(i) = 0.0_dp
    do ni = 1,nnbr(i)
      j = nbrno(ni,i)
      qtmp(j) = qtmp(j) + q(i)/dble(nnbr(i))  ! could be a bridging H
    enddo
  enddo

  q = qtmp

  end subroutine qheavy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! neighbour only version of EEQ model
! included up to 1,4 interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine goedecker_topo(numat,nfrag,nfraglist,qfrag,q,es)
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_topo
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)    :: numat                ! Number of atoms
  integer(i4),       intent(in)    :: nfrag                ! Number of fragments
  integer(i4),       intent(in)    :: nfraglist(numat)     ! Pointer from atom to fragment number
  real(dp),          intent(inout) :: q(*)                 ! Charges
  real(dp),          intent(in)    :: qfrag(numat)         ! Fragment charges
  real(dp),          intent(out)   :: es                   ! ES energy
!
!  Local variables
!
  integer(i4)                      :: i
  integer(i4)                      :: j
  integer(i4)                      :: m
  integer(i4)                      :: ni
  integer(i4)                      :: ierror
  real(dp)                         :: conversion
  real(dp)                         :: es1
  real(dp)                         :: es2
  real(dp)                         :: gammij
  real(dp)                         :: rij
  real(dp)                         :: tsqrt2pi
  real(dp)                         :: tmp
  real(dp)                         :: wtrm
  real(dp),      allocatable, save :: A(:,:)
  real(dp),      allocatable, save :: x(:)
!
!  Parameters
!
  parameter (tsqrt2pi = 0.797884560802866_dp)
!
  conversion = xtb_autoaa*xtb_autoev
!
  m = numat + nfrag ! # atoms frag constrain
!
!  Allocate local memory
!
  allocate(A(m,m),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : A '')')
    stop
  endif
  allocate(x(m),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : x '')')
    stop
  endif
!
  A(1:m,1:m) = 0
!
!  Setup RHS
!
  do i = 1,numat
    x(i) = chieeq(i) ! EN of atom
  enddo
  do i = 1,numat
    A(i,i) = (gameeq(i) + tsqrt2pi/sqrt(alpeeq(i)))*conversion
  enddo
!
!  Setup A matrix by adding Coulomb interactions
!
  do i = 1,numat
    do ni = 1,nnbr_topo(i)
      j = nbrno_topo(ni,i)
      rij = rtnbr(ni,i)
      gammij = 1.0_dp/sqrt(alpeeq(i)+alpeeq(j)) ! squared above
      if (lgfnff_topowolf) then
        wtrm = derfc(gfnff_wolf_eta*rij)/rij
        tmp = erf(gammij*rij)*conversion*(wtrm - gfnff_wolf_self)
      else
        tmp = erf(gammij*rij)*conversion/rij
      endif
      A(j,i) = A(j,i) + tmp
    enddo
  enddo
!
!  Fragment charge constraint
!
  do i = 1,nfrag
    x(numat+i) = qfrag(i)
  enddo
  do i = 1,nfrag
    do j = 1,numat
      if (nfraglist(j).eq.i) then
        A(j,numat+i) = 1
      endif
    enddo
  enddo
  do i = 1,numat
    do j = 1,nfrag
      if (nfraglist(i).eq.j) then
        A(numat+j,i) = 1
      endif
    enddo
  enddo
!
!  Solve for charges
!
  call pgfnff_matsolve(numat,nfrag,m,m,A,x,q,qfrag)
!
!  Energy
!
  es1 = 0.0_dp
  es2 = 0.0_dp
  do i = 1,numat
    do ni = 1,nnbr_topo(i)
      j = nbrno_topo(ni,i)
      rij = rtnbr(ni,i)
      gammij = 1.0_dp/sqrt(alpeeq(i)+alpeeq(j)) ! squared above
      if (lgfnff_topowolf) then
        wtrm = derfc(gfnff_wolf_eta*rij)/rij
        tmp = erf(gammij*rij)*(wtrm - gfnff_wolf_self)
      else
        tmp = erf(gammij*rij)/rij
      endif
      es1 = es1 + 0.5_dp*q(i)*q(j)*tmp
    enddo
    es1 = es1 + q(i)*q(i)*0.5_dp*(gameeq(i)+tsqrt2pi/sqrt(alpeeq(i)))
    es2 = es2 - q(i)*chieeq(i)
  enddo
  es = es1*conversion + es2
!
!  Free local memory
!
  deallocate(x,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : x '')')
    stop
  endif
  deallocate(A,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : A '')')
    stop
  endif

  end subroutine goedecker_topo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! neighbour only version of EEQ model
! included up to 1,4 interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine goedeckera(numat,rtopo,nfrag,nfraglist,qfrag,q,es)
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),       intent(in)    :: numat                ! Number of atoms
  integer(i4),       intent(in)    :: nfrag                ! Number of fragments
  integer(i4),       intent(in)    :: nfraglist(numat)     ! Pointer from atom to fragment number
  real(dp),          intent(in)    :: rtopo(numat,numat)   ! Distances between atoms
  real(dp),          intent(inout) :: q(*)                 ! Charges
  real(dp),          intent(in)    :: qfrag(numat)         ! Fragment charges
  real(dp),          intent(out)   :: es                   ! ES energy
!
!  Local variables
!
  integer(i4)                      :: i
  integer(i4)                      :: j
  integer(i4)                      :: m
  integer(i4)                      :: ierror
  real(dp)                         :: conversion
  real(dp)                         :: es1
  real(dp)                         :: es2
  real(dp)                         :: gammij
  real(dp)                         :: rij
  real(dp)                         :: tsqrt2pi
  real(dp)                         :: tmp
  real(dp),      allocatable, save :: A(:,:)
  real(dp),      allocatable, save :: x(:)
!
!  Parameter
!
  parameter (tsqrt2pi = 0.797884560802866_dp)
!
  conversion = xtb_autoaa*xtb_autoev
!
  m = numat + nfrag ! # atoms frag constrain
!
!  Allocate local memory
!
  allocate(A(m,m),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : A '')')
    stop
  endif
  allocate(x(m),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : x '')')
    stop
  endif
!
  A(1:m,1:m) = 0
!
!  Setup RHS
!
  do i = 1,numat
    x(i) = chieeq(i) ! EN of atom
  enddo
  do i = 1,numat
    A(i,i) = (gameeq(i) + tsqrt2pi/sqrt(alpeeq(i)))*conversion
  enddo
!
!  Setup A matrix
!
  do i = 2,numat
    do j = 1,i-1
      rij = rtopo(j,i)
      gammij = 1.0_dp/sqrt(alpeeq(i)+alpeeq(j)) ! squared above
      tmp = erf(gammij*rij)*conversion/rij
      A(j,i) = A(j,i) + tmp
      A(i,j) = A(i,j) + tmp
    enddo
  enddo
!
!  Fragment charge constraint
!
  do i = 1,nfrag
    x(numat+i) = qfrag(i)
  enddo
  do i = 1,nfrag
    do j = 1,numat
      if (nfraglist(j).eq.i) then
        A(j,numat+i) = 1
      endif
    enddo
  enddo
  do i = 1,numat
    do j = 1,nfrag
      if (nfraglist(i).eq.j) then
        A(numat+j,i) = 1
      endif
    enddo
  enddo
!
!  Solve for charges
!
  call pgfnff_matsolve(numat,nfrag,m,m,A,x,q,qfrag)
!
!  Energy
!
  es1 = 0.0_dp
  es2 = 0.0_dp
  do i = 1,numat
    do j = 1,i-1
      rij = rtopo(j,i)
      gammij = 1.0_dp/sqrt(alpeeq(i)+alpeeq(j)) ! squared above
      tmp = erf(gammij*rij)/rij
      es1 = es1 + q(i)*q(j)*tmp
    enddo
    es1 = es1 + q(i)*q(i)*0.5_dp*(gameeq(i)+tsqrt2pi/sqrt(alpeeq(i)))
    es2 = es2 - q(i)*chieeq(i)
  enddo
  es = es1*conversion + es2
!
!  Free local memory
!
  deallocate(x,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : x '')')
    stop
  endif
  deallocate(A,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : A '')')
    stop
  endif

  end subroutine goedeckera

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function amide(numat,nat,hyb,nnbr,maxnbr,nbrno,gnbrno,pi,ia,lxref)
  use m_pgfnff_types
!
!  Passed variables
!
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: nat(numat)
  integer(i4), intent(in)  :: hyb(numat)
  integer(i4), intent(in)  :: nnbr(numat)
  integer(i4), intent(in)  :: maxnbr
  integer(i4), intent(in)  :: nbrno(maxnbr,numat)
  integer(i4), intent(in)  :: gnbrno(maxnbr,numat)
  integer(i4), intent(in)  :: pi(numat)
  integer(i4), intent(in)  :: ia
  logical,     intent(in)  :: lxref   ! If true then use gnbrno
!
!  Local variables
!
  integer(i4)              :: j
  integer(i4)              :: ni
  integer(i4)              :: no
  integer(i4)              :: nc
  integer(i4)              :: ic
!
  amide = .false. ! don't know
  if (pi(ia).eq.0.or.hyb(ia).ne.3.or.nat(ia).ne.7) return
!
  nc = 0
  no = 0
  if (lxref) then
    do ni = 1,nnbr(ia)
      j = gnbrno(nbrno(ni,ia),ia)
      if (nat(j).eq.6.and.pi(j).ne.0) then  ! a pi C on N?
        nc = nc + 1
        ic = j
      endif
    enddo
  else
    do ni = 1,nnbr(ia)
      j = nbrno(ni,ia)
      if (nat(j).eq.6.and.pi(j).ne.0) then  ! a pi C on N?
        nc = nc + 1
        ic = j
      endif
    enddo
  endif
!
  if (nc.eq.1) then
    if (lxref) then
      do ni = 1,nnbr(ic)
        j = gnbrno(nbrno(ni,ic),ic)
        if (nat(j).eq.8.and.pi(j).ne.0.and.nnbr(j).eq.1) no = no + 1 ! a pi =O on the C?
      enddo
    else
      do ni = 1,nnbr(ic)
        j = nbrno(ni,ic)
        if (nat(j).eq.8.and.pi(j).ne.0.and.nnbr(j).eq.1) no = no + 1 ! a pi =O on the C?
      enddo
    endif
  endif
!
  if (no.eq.1) amide = .true.

  end function amide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function amideH(numat,nat,hyb,nnbr,maxnbr,nbrno,gnbrno,pi,ia,lxref)
  use m_pgfnff_types
!
!  Passed variables
!
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: nat(numat)
  integer(i4), intent(in)  :: hyb(numat)
  integer(i4), intent(in)  :: nnbr(numat)
  integer(i4), intent(in)  :: maxnbr
  integer(i4), intent(in)  :: nbrno(maxnbr,numat)
  integer(i4), intent(in)  :: gnbrno(maxnbr,numat)
  integer(i4), intent(in)  :: pi(numat)
  integer(i4), intent(in)  :: ia
  logical,     intent(in)  :: lxref   ! If true then use gnbrno
!
!  Local variables
!
  integer(i4)              :: j
  integer(i4)              :: ni
  integer(i4)              :: nc
  integer(i4)              :: nn
  logical                  :: amide
!
  amideH = .false.
  if (nnbr(ia).ne.1) return
  if (lxref) then
    nn = gnbrno(nbrno(1,ia),ia)     ! the N
  else
    nn = nbrno(1,ia)
  endif
  if (.not.amide(numat,nat,hyb,nnbr,maxnbr,nbrno,gnbrno,pi,nn,lxref)) return
!
  nc = 0
  if (lxref) then
    do ni = 1,nnbr(nn)
      j = gnbrno(nbrno(ni,nn),nn)
      if (nat(j).eq.6.and.hyb(j).eq.3) then  ! a sp3 C on N?
        nc = nc + 1
      endif
    enddo
  else
    do ni = 1,nnbr(nn)
      j = nbrno(ni,nn)
      if (nat(j).eq.6.and.hyb(j).eq.3) then  ! a sp3 C on N?
        nc = nc + 1
      endif
    enddo
  endif
!
  if (nc.eq.1) amideH = .true.

  end function amideH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ctype                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(i4) function ctype(numat,nat,nnbr,maxnbr,nbrno,pi,ia)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: numat
  integer(i4), intent(in) :: nat(numat)
  integer(i4), intent(in) :: nnbr(numat)
  integer(i4), intent(in) :: maxnbr
  integer(i4), intent(in) :: nbrno(maxnbr,numat)
  integer(i4), intent(in) :: pi(numat)
  integer(i4), intent(in) :: ia
!
!  Local variables
!
  integer(i4)             :: ni
  integer(i4)             :: j
  integer(i4)             :: no
!
  ctype = 0 ! don't know
  no = 0
  do ni = 1,nnbr(ia)
    j = nbrno(ni,ia)
    if (nat(j).eq.8.and.pi(j).ne.0) no = no + 1
  enddo
!
  if (no.eq.1.and.pi(ia).ne.0) ctype = 1 ! a C=O carbon

  end function ctype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Ring analysis routine                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getring36(ndim,numat,nnbr,maxnbr,nbrno,xnbr,ynbr,znbr,lcluster,a0,nra,nrs,rs,nrings)
  use m_pgfnff_types
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Passed variables
!
  integer(i4),                       intent(in)  :: ndim
  integer(i4),                       intent(in)  :: numat
  integer(i4),                       intent(in)  :: nnbr(numat)
  integer(i4),                       intent(in)  :: maxnbr
  integer(i4),                       intent(in)  :: nbrno(maxnbr,numat)
  integer(i4),                       intent(in)  :: a0
  integer(i4),                       intent(out) :: nra(6,20)              ! output: atomlist
  integer(i4),                       intent(out) :: nrs(20)                ! output: ringsize
  logical,                           intent(in)  :: lcluster(numat)
  integer(i4),                       intent(out) :: nrings
  real(dp),                          intent(out) :: rs(3,6,20)             ! output: vectors for ring sides
  real(dp),                          intent(in)  :: xnbr(maxnbr,numat)
  real(dp),                          intent(in)  :: ynbr(maxnbr,numat)
  real(dp),                          intent(in)  :: znbr(maxnbr,numat)
!
!  Local variables
!
  integer(i4)                                    :: i,ir1,ir2,ir3,ir4,ir5,ir6
  integer(i4)                                    :: n0,n1,n2,n3,n4,n5,n6
  integer(i4)                                    :: a1,a2,a3,a4,a5,a6
  integer(i4),                              save :: maxr  = 100
  integer(i4), dimension(:),       pointer, save :: nrsl  => null() ! Ring size - local copy
  integer(i4), dimension(:,:),     pointer, save :: nral  => null() ! Atoms in ring - local copy
  integer(i4), dimension(:),       pointer, save :: same  => null() ! Flag as to whether rings are duplicated
  integer(i4)                                    :: c(6)
  integer(i4)                                    :: ierror
  integer(i4)                                    :: iring
  integer(i4)                                    :: j
  integer(i4)                                    :: kk
  integer(i4)                                    :: m
  integer(i4)                                    :: nn
  logical,                                  save :: first = .true.
  real(dp),    dimension(:),       pointer, save :: cnorm => null()
  real(dp),    dimension(:,:),     pointer, save :: cvec  => null()
  real(dp),    dimension(:,:,:),   pointer, save :: rsl   => null()
  real(dp)                                       :: dum1
  real(dp)                                       :: rv20
  real(dp)                                       :: rv30
  real(dp)                                       :: rv31
  real(dp)                                       :: rv40
  real(dp)                                       :: rv41
  real(dp)                                       :: rv42
  real(dp)                                       :: rv50
  real(dp)                                       :: rv51
  real(dp)                                       :: rv52
  real(dp)                                       :: rv53
  real(dp)                                       :: rv60
  real(dp)                                       :: rv61
  real(dp)                                       :: rv62
  real(dp)                                       :: rv63
  real(dp)                                       :: rv64
  real(dp)                                       :: v(3,6)
  real(dp)                                       :: v20(3)
  real(dp)                                       :: v30(3)
  real(dp)                                       :: v31(3)
  real(dp)                                       :: v40(3)
  real(dp)                                       :: v41(3)
  real(dp)                                       :: v42(3)
  real(dp)                                       :: v50(3)
  real(dp)                                       :: v51(3)
  real(dp)                                       :: v52(3)
  real(dp)                                       :: v53(3)
  real(dp)                                       :: v60(3)
  real(dp)                                       :: v61(3)
  real(dp)                                       :: v62(3)
  real(dp)                                       :: v63(3)
  real(dp)                                       :: v64(3)
!
!  First call allocate arrays
!
  if (first) then
    call pgfnff_realloc(nral,6_i4,maxr,ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : nral '')')
      stop
    endif
    call pgfnff_realloc(nrsl,maxr,ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : nrsl '')')
      stop
    endif
    call pgfnff_realloc(same,maxr,ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : same '')')
      stop
    endif
    call pgfnff_realloc(cnorm,maxr,ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : cnorm '')')
      stop
    endif
    call pgfnff_realloc(cvec,3_i4,maxr,ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : cvec '')')
      stop
    endif
    call pgfnff_realloc(rsl,3_i4,6_i4,maxr,ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : rsl '')')
      stop
    endif
    first = .false.
  endif
!
  nra = 0
  nrs = 0
  rs  = 0.0_dp
  nrings = 0
!
  if ((ndim.eq.0.and.numat.le.2).or.lcluster(a0)) return
!
  nral = 0
  kk = 0
!
!  Loop over neighbours of initial atom, a0
!
  n0 = nnbr(a0)
  do ir1 = 1,n0
    a1 = nbrno(ir1,a0)
    n1 = nnbr(a1)
!
!  If neighbour only has 1 neighbour then this can't lead to a ring and so skip
!
    if (n1.eq.1) cycle
!
    v(1,1) = xnbr(ir1,a0)
    v(2,1) = ynbr(ir1,a0)
    v(3,1) = znbr(ir1,a0)
!
!  Loop over 2nd neighbours
!
    do ir2 = 1,n1
      a2 = nbrno(ir2,a1)
      n2 = nnbr(a2)
!
!  If neighbour only has 1 neighbour then this can't lead to a ring and so skip
!
      if (n2.eq.1) cycle
!
      v(1,2) = xnbr(ir2,a1)
      v(2,2) = ynbr(ir2,a1)
      v(3,2) = znbr(ir2,a1)
!
      v20(1) = v(1,2) + v(1,1)
      v20(2) = v(2,2) + v(2,1)
      v20(3) = v(3,2) + v(3,1)
      rv20 = v20(1)**2 + v20(2)**2 + v20(3)**2
!
!  If a2 is the same as a0 then cycle
!
      if (a2.eq.a0.and.rv20.lt.1.0d-2) cycle
!
!  Loop over 3rd neighbours
!
      do ir3 = 1,n2
        a3 = nbrno(ir3,a2)
        n3 = nnbr(a3)
!
!  If neighbour only has 1 neighbour then this can't lead to a ring and so skip
!
        if (n3.eq.1) cycle
!
        v(1,3) = xnbr(ir3,a2)
        v(2,3) = ynbr(ir3,a2)
        v(3,3) = znbr(ir3,a2)
!
        v31(1) = v(1,3) + v(1,2)
        v31(2) = v(2,3) + v(2,2)
        v31(3) = v(3,3) + v(3,2)
!
        rv31 = v31(1)**2 + v31(2)**2 + v31(3)**2
!
!  If a3 is the same as a1 then cycle
!
        if (a3.eq.a1.and.rv31.lt.1.0d-2) cycle
!
        v30(1) = v(1,3) + v20(1)
        v30(2) = v(2,3) + v20(2)
        v30(3) = v(3,3) + v20(3)
!
        rv30 = v30(1)**2 + v30(2)**2 + v30(3)**2
!
        c(1) = a0
        c(2) = a1
        c(3) = a2
!
        if (a3.eq.a0.and.rv30.lt.1.0d-2) then
!-----------------
!  3-ring found  |
!-----------------
          iring = 3
!
          if (kk.eq.maxr) then
            maxr = maxr + 100_i4
            call pgfnff_realloc(nral,6_i4,maxr,ierror)
            if (ierror.gt.0) then
              write(ioout,'('' Error in Memory Allocation : nral '')')
              stop
            endif
            call pgfnff_realloc(nrsl,maxr,ierror)
            if (ierror.gt.0) then
              write(ioout,'('' Error in Memory Allocation : nrsl '')')
              stop
            endif
            call pgfnff_realloc(same,maxr,ierror)
            if (ierror.gt.0) then
              write(ioout,'('' Error in Memory Allocation : same '')')
              stop
            endif
            call pgfnff_realloc(cnorm,maxr,ierror)
            if (ierror.gt.0) then
              write(ioout,'('' Error in Memory Allocation : cnorm '')')
              stop
            endif
            call pgfnff_realloc(cvec,3_i4,maxr,ierror)
            if (ierror.gt.0) then
              write(ioout,'('' Error in Memory Allocation : cvec '')')
              stop
            endif
            call pgfnff_realloc(rsl,3_i4,6_i4,maxr,ierror)
            if (ierror.gt.0) then
              write(ioout,'('' Error in Memory Allocation : rsl '')')
              stop
            endif
          endif
!
          kk = kk + 1
          nral(1:iring,kk) = c(1:iring)
          cvec(1,kk) = 0.5_dp*(v20(1) + v(1,1))
          cvec(2,kk) = 0.5_dp*(v20(2) + v(2,1))
          cvec(3,kk) = 0.5_dp*(v20(3) + v(3,1))
          cnorm(kk) = cvec(1,kk)**2 + cvec(2,kk)**2 + cvec(3,kk)**2
          if (cnorm(kk).gt.1.0d-6) cnorm(kk) = 1.0_dp/sqrt(cnorm(kk))
          nrsl(kk) = iring
          rsl(1:3,1:iring,kk) = v(1:3,1:iring)
!
!  Having found 3 ring then cycle as subsequent rings would be nested
!
          cycle
        endif
!
!  Loop over 4th neighbours
!
        do ir4 = 1,n3
          a4 = nbrno(ir4,a3)
          n4 = nnbr(a4)
!
!  If neighbour only has 1 neighbour then this can't lead to a ring and so skip
!
          if (n4.eq.1) cycle
!
          v(1,4) = xnbr(ir4,a3)
          v(2,4) = ynbr(ir4,a3)
          v(3,4) = znbr(ir4,a3)
!
          v42(1) = v(1,4) + v(1,3)
          v42(2) = v(2,4) + v(2,3)
          v42(3) = v(3,4) + v(3,3)
!
          rv42 = v42(1)**2 + v42(2)**2 + v42(3)**2
!
!  If a4 is the same as a2 then cycle
!
          if (a4.eq.a2.and.rv42.lt.1.0d-2) cycle
!
          v41(1) = v(1,2) + v42(1)
          v41(2) = v(2,2) + v42(2)
          v41(3) = v(3,2) + v42(3)
!
          rv41 = v41(1)**2 + v41(2)**2 + v41(3)**2
!
!  If a4 is the same as a1 then cycle
!
          if (a4.eq.a1.and.rv41.lt.1.0d-2) cycle
!
          v40(1) = v(1,4) + v30(1)
          v40(2) = v(2,4) + v30(2)
          v40(3) = v(3,4) + v30(3)
!
          rv40 = v40(1)**2 + v40(2)**2 + v40(3)**2
!
          c(4) = a3
!
          if (a4.eq.a0.and.rv40.lt.1.0d-2) then
!-----------------
!  4-ring found  |
!-----------------
            iring = 4
!
            if (kk.eq.maxr) then
              maxr = maxr + 100_i4
              call pgfnff_realloc(nral,6_i4,maxr,ierror)
              if (ierror.gt.0) then
                write(ioout,'('' Error in Memory Allocation : nral '')')
                stop
              endif
              call pgfnff_realloc(nrsl,maxr,ierror)
              if (ierror.gt.0) then
                write(ioout,'('' Error in Memory Allocation : nrsl '')')
                stop
              endif
              call pgfnff_realloc(same,maxr,ierror)
              if (ierror.gt.0) then
                write(ioout,'('' Error in Memory Allocation : same '')')
                stop
              endif
              call pgfnff_realloc(cnorm,maxr,ierror)
              if (ierror.gt.0) then
                write(ioout,'('' Error in Memory Allocation : cnorm '')')
                stop
              endif
              call pgfnff_realloc(cvec,3_i4,maxr,ierror)
              if (ierror.gt.0) then
                write(ioout,'('' Error in Memory Allocation : cvec '')')
                stop
              endif
              call pgfnff_realloc(rsl,3_i4,6_i4,maxr,ierror)
              if (ierror.gt.0) then
                write(ioout,'('' Error in Memory Allocation : rsl '')')
                stop
              endif
            endif
!
            kk = kk + 1
            nral(1:iring,kk) = c(1:iring)
            cvec(1,kk) = (v30(1) + v20(1) + v(1,1))/3.0_dp
            cvec(2,kk) = (v30(2) + v20(2) + v(2,1))/3.0_dp
            cvec(3,kk) = (v30(3) + v20(3) + v(3,1))/3.0_dp
            cnorm(kk) = cvec(1,kk)**2 + cvec(2,kk)**2 + cvec(3,kk)**2
            if (cnorm(kk).gt.1.0d-6) cnorm(kk) = 1.0_dp/sqrt(cnorm(kk))
            nrsl(kk) = iring
            rsl(1:3,1:iring,kk) = v(1:3,1:iring)
!
!  Having found 4 ring then cycle as subsequent rings would be nested
!
            cycle
          endif
!
!  Loop over 5th neighbours
!
          do ir5 = 1,n4
            a5 = nbrno(ir5,a4)
            n5 = nnbr(a5)
!
!  If neighbour only has 1 neighbour then this can't lead to a ring and so skip
!
            if (n5.eq.1) cycle
!
            v(1,5) = xnbr(ir5,a4)
            v(2,5) = ynbr(ir5,a4)
            v(3,5) = znbr(ir5,a4)
!
            v53(1) = v(1,5) + v(1,4)
            v53(2) = v(2,5) + v(2,4)
            v53(3) = v(3,5) + v(3,4)
!
            rv53 = v53(1)**2 + v53(2)**2 + v53(3)**2
!
!  If a5 is the same as a3 then cycle
!
            if (a5.eq.a3.and.rv53.lt.1.0d-2) cycle
!
            v52(1) = v(1,5) + v42(1)
            v52(2) = v(2,5) + v42(2)
            v52(3) = v(3,5) + v42(3)
!
            rv52 = v52(1)**2 + v52(2)**2 + v52(3)**2
!
!  If a5 is the same as a2 then cycle
!
            if (a5.eq.a2.and.rv52.lt.1.0d-2) cycle
!
            v51(1) = v(1,5) + v41(1)
            v51(2) = v(2,5) + v41(2)
            v51(3) = v(3,5) + v41(3)
!
            rv51 = v51(1)**2 + v51(2)**2 + v51(3)**2
!
!  If a5 is the same as a1 then cycle
!
            if (a5.eq.a1.and.rv51.lt.1.0d-2) cycle
!
            v50(1) = v(1,5) + v40(1)
            v50(2) = v(2,5) + v40(2)
            v50(3) = v(3,5) + v40(3)
!
            rv50 = v50(1)**2 + v50(2)**2 + v50(3)**2
!
            c(5) = a4
!
            if (a5.eq.a0.and.rv50.lt.1.0d-2) then
!-----------------
!  5-ring found  |
!-----------------
              iring = 5
!
              if (kk.eq.maxr) then
                maxr = maxr + 100_i4
                call pgfnff_realloc(nral,6_i4,maxr,ierror)
                if (ierror.gt.0) then
                  write(ioout,'('' Error in Memory Allocation : nral '')')
                  stop
                endif
                call pgfnff_realloc(nrsl,maxr,ierror)
                if (ierror.gt.0) then
                  write(ioout,'('' Error in Memory Allocation : nrsl '')')
                  stop
                endif
                call pgfnff_realloc(same,maxr,ierror)
                if (ierror.gt.0) then
                  write(ioout,'('' Error in Memory Allocation : same '')')
                  stop
                endif
                call pgfnff_realloc(cnorm,maxr,ierror)
                if (ierror.gt.0) then
                  write(ioout,'('' Error in Memory Allocation : cnorm '')')
                  stop
                endif
                call pgfnff_realloc(cvec,3_i4,maxr,ierror)
                if (ierror.gt.0) then
                  write(ioout,'('' Error in Memory Allocation : cvec '')')
                  stop
                endif
                call pgfnff_realloc(rsl,3_i4,6_i4,maxr,ierror)
                if (ierror.gt.0) then
                  write(ioout,'('' Error in Memory Allocation : rsl '')')
                  stop
                endif
              endif
!
              kk = kk + 1
              nral(1:iring,kk) = c(1:iring)
              cvec(1,kk) = 0.25_dp*(v40(1) + v30(1) + v20(1) + v(1,1))
              cvec(2,kk) = 0.25_dp*(v40(2) + v30(2) + v20(2) + v(2,1))
              cvec(3,kk) = 0.25_dp*(v40(3) + v30(3) + v20(3) + v(3,1))
              cnorm(kk) = cvec(1,kk)**2 + cvec(2,kk)**2 + cvec(3,kk)**2
              if (cnorm(kk).gt.1.0d-6) cnorm(kk) = 1.0_dp/sqrt(cnorm(kk))
              nrsl(kk) = iring
              rsl(1:3,1:iring,kk) = v(1:3,1:iring)
!
!  Having found 5 ring then cycle as subsequent rings would be nested
!
              cycle
            endif
!
!  Loop over 6th neighbours
!
            do ir6 = 1,n5
              a6 = nbrno(ir6,a5)
              n6 = nnbr(a6)
!
!  If neighbour only has 1 neighbour then this can't lead to a ring and so skip
!
              if (n6.eq.1) cycle
!
              v(1,6) = xnbr(ir6,a5)
              v(2,6) = ynbr(ir6,a5)
              v(3,6) = znbr(ir6,a5)
!
              v64(1) = v(1,6) + v(1,5)
              v64(2) = v(2,6) + v(2,5)
              v64(3) = v(3,6) + v(3,5)
!
              rv64 = v64(1)**2 + v64(2)**2 + v64(3)**2
!
!  If a6 is the same as a4 then cycle
!
              if (a6.eq.a4.and.rv64.lt.1.0d-2) cycle
!
              v63(1) = v(1,6) + v53(1)
              v63(2) = v(2,6) + v53(2)
              v63(3) = v(3,6) + v53(3)
!
              rv63 = v63(1)**2 + v63(2)**2 + v63(3)**2
!
!  If a6 is the same as a3 then cycle
!
              if (a6.eq.a3.and.rv63.lt.1.0d-2) cycle
!
              v62(1) = v(1,6) + v52(1)
              v62(2) = v(2,6) + v52(2)
              v62(3) = v(3,6) + v52(3)
!
              rv62 = v62(1)**2 + v62(2)**2 + v62(3)**2
!
!  If a6 is the same as a2 then cycle
!
              if (a6.eq.a2.and.rv62.lt.1.0d-2) cycle
!
              v61(1) = v(1,6) + v51(1)
              v61(2) = v(2,6) + v51(2)
              v61(3) = v(3,6) + v51(3)
!
              rv61 = v61(1)**2 + v61(2)**2 + v61(3)**2
!
!  If a6 is the same as a1 then cycle
!
              if (a6.eq.a1.and.rv61.lt.1.0d-2) cycle
!
              v60(1) = v(1,6) + v50(1)
              v60(2) = v(2,6) + v50(2)
              v60(3) = v(3,6) + v50(3)
!
              rv60 = v60(1)**2 + v60(2)**2 + v60(3)**2
!
              c(6) = a5
!
              if ((a6.eq.a0.and.rv60.lt.1.0d-2)) then
!-----------------
!  6-ring found  |
!-----------------
                iring = 6
!
                if (kk.eq.maxr) then
                  maxr = maxr + 100_i4
                  call pgfnff_realloc(nral,6_i4,maxr,ierror)
                  if (ierror.gt.0) then
                    write(ioout,'('' Error in Memory Allocation : nral '')')
                    stop
                  endif
                  call pgfnff_realloc(nrsl,maxr,ierror)
                  if (ierror.gt.0) then
                    write(ioout,'('' Error in Memory Allocation : nrsl '')')
                    stop
                  endif
                  call pgfnff_realloc(same,maxr,ierror)
                  if (ierror.gt.0) then
                    write(ioout,'('' Error in Memory Allocation : same '')')
                    stop
                  endif
                  call pgfnff_realloc(cnorm,maxr,ierror)
                  if (ierror.gt.0) then
                    write(ioout,'('' Error in Memory Allocation : cnorm '')')
                    stop
                  endif
                  call pgfnff_realloc(cvec,3_i4,maxr,ierror)
                  if (ierror.gt.0) then
                    write(ioout,'('' Error in Memory Allocation : cvec '')')
                    stop
                  endif
                  call pgfnff_realloc(rsl,3_i4,6_i4,maxr,ierror)
                  if (ierror.gt.0) then
                    write(ioout,'('' Error in Memory Allocation : rsl '')')
                    stop
                  endif
                endif
!
                kk = kk + 1
                nral(1:iring,kk) = c(1:iring)
                cvec(1,kk) = 0.2_dp*(v50(1) + v40(1) + v30(1) + v20(1) + v(1,1))
                cvec(2,kk) = 0.2_dp*(v50(2) + v40(2) + v30(2) + v20(2) + v(2,1))
                cvec(3,kk) = 0.2_dp*(v50(3) + v40(3) + v30(3) + v20(3) + v(3,1))
                cnorm(kk) = cvec(1,kk)**2 + cvec(2,kk)**2 + cvec(3,kk)**2
                if (cnorm(kk).gt.1.0d-6) cnorm(kk) = 1.0_dp/sqrt(cnorm(kk))
                nrsl(kk) = iring
                rsl(1:3,1:iring,kk) = v(1:3,1:iring)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
!------------------
!  Compare rings  !
!------------------
  same = 0
  do i = 1,kk
    do j = i+1,kk
      if (nrsl(i).ne.nrsl(j)) cycle ! different ring size
      if (same(j).eq.1      ) cycle ! already done
!
!  Dot product of normalised centre of ring vectors
!
      dum1 = (cvec(1,i)*cvec(1,j) + cvec(2,i)*cvec(2,j) + cvec(3,i)*cvec(3,j))*cnorm(i)*cnorm(j)
      if (abs(dum1-1.0_dp).lt.1.0d-2) then
        same(j) = 1
      else
        same(j) = 0
      endif
    enddo
  enddo
!-----------------------------------------
!  Reduce to the subset of unique rings  |
!-----------------------------------------
  m = 0
  do i = 1,kk
    if (same(i).eq.0) then
      m = m + 1
      nrs(m) = nrsl(i)     ! number of atoms in ring m
      nn = nrsl(i)
      nra(1:nn,m) = nral(1:nn,i)
      rs(1:3,1:nn,m) = rsl(1:3,1:nn,i)
      if (m.gt.19) then
        m = 19
        goto 999
      endif
    endif
  enddo
999 nrings = m  ! number of rings for this atom

  return
  end subroutine getring36

!!!!!!!!!!
!  Sort  !
!!!!!!!!!!
  subroutine ssort(n,edum,ind)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: n
  real(dp),    intent(inout) :: edum(n)
  integer(i4), intent(inout) :: ind(n)
!
!  Local variables
!
  integer(i4) :: i
  integer(i4) :: ii
  integer(i4) :: j
  integer(i4) :: k
  integer(i4) :: sc1
  real(dp)    :: pp
!
  iiloop: do ii = 2, n
    i = ii - 1
    k = i
    pp = edum(i)
    jloop: do j = ii, n
      if (edum(j).gt.pp) cycle jloop
      k = j
      pp = edum(j)
    enddo jloop
    if (k.eq.i) cycle iiloop
    edum(k) = edum(i)
    edum(i) = pp
    sc1 = ind(i)
    ind(i) = ind(k)
    ind(k) = sc1
  enddo iiloop

  end subroutine ssort

!!!!!!!!!!!!!!!!!!!
!  Zeta function  !
!!!!!!!!!!!!!!!!!!!
  function zeta(nati,qi)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in) :: nati
  real(dp),     intent(in) :: qi
!
!  Local variables
!
  real(dp)                 :: zeta,qmod
  real(dp),      parameter :: zeff(86) = (/ &
   &   1,                                                 2,  & ! H-He
   &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
   &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
   &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
   &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
   &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
   &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26/)  ! Hf-Rn
!
!  Semiempirical Evaluation of the GlobalHardness of the Atoms of 103
!  Elements of the Periodic Table Using the Most Probable Radii as
!  their Size Descriptors DULAL C. GHOSH, NAZMUL ISLAM 2009 in
!  Wiley InterScience (www.inte"rscience.wiley.com).
!  DOI 10.1002/qua.22202
!  values in the paper multiplied by two because
!  (ii:ii)=(IP-EA)=d^2 E/dN^2 but the hardness
!  definition they use is 1/2d^2 E/dN^2 (in Eh)
!
   real(dp),parameter :: c(1:86) = (/ &
  &0.47259288_dp,0.92203391_dp,0.17452888_dp,0.25700733_dp,0.33949086_dp,0.42195412_dp, & ! H-C
  &0.50438193_dp,0.58691863_dp,0.66931351_dp,0.75191607_dp,0.17964105_dp,0.22157276_dp, & ! N-Mg
  &0.26348578_dp,0.30539645_dp,0.34734014_dp,0.38924725_dp,0.43115670_dp,0.47308269_dp, & ! Al-Ar
  &0.17105469_dp,0.20276244_dp,0.21007322_dp,0.21739647_dp,0.22471039_dp,0.23201501_dp, & ! Ca-Cr
  &0.23933969_dp,0.24665638_dp,0.25398255_dp,0.26128863_dp,0.26859476_dp,0.27592565_dp, & ! Mn-Zn
  &0.30762999_dp,0.33931580_dp,0.37235985_dp,0.40273549_dp,0.43445776_dp,0.46611708_dp, & ! Ga-Kr
  &0.15585079_dp,0.18649324_dp,0.19356210_dp,0.20063311_dp,0.20770522_dp,0.21477254_dp, & ! Rb-Mo
  &0.22184614_dp,0.22891872_dp,0.23598621_dp,0.24305612_dp,0.25013018_dp,0.25719937_dp, & ! Tc-Cd
  &0.28784780_dp,0.31848673_dp,0.34912431_dp,0.37976593_dp,0.41040808_dp,0.44105777_dp, & ! In-Xe
  &0.05019332_dp,0.06762570_dp,0.08504445_dp,0.10247736_dp,0.11991105_dp,0.13732772_dp, & ! Cs-Nd
  &0.15476297_dp,0.17218265_dp,0.18961288_dp,0.20704760_dp,0.22446752_dp,0.24189645_dp, & ! Pm-Dy
  &0.25932503_dp,0.27676094_dp,0.29418231_dp,0.31159587_dp,0.32902274_dp,0.34592298_dp, & ! Ho-Hf
  &0.36388048_dp,0.38130586_dp,0.39877476_dp,0.41614298_dp,0.43364510_dp,0.45104014_dp, & ! Ta-Pt
  &0.46848986_dp,0.48584550_dp,0.12526730_dp,0.14268677_dp,0.16011615_dp,0.17755889_dp, & ! Au-Po
  &0.19497557_dp,0.21240778_dp/)

  intrinsic :: exp

  qmod = zeff(nati) + qi
  if (qmod.lt.0.0_dp) then
    zeta = exp( 3.0_dp )
  else
    zeta = exp(3.0_dp*(1.0_dp - exp(c(nati)*(1.0_dp - zeff(nati)/qmod))))
  endif
  end function zeta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Xatom logical function  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function lxatom(nati)
  use m_pgfnff_types
  implicit none
!
!  Passed arguments
!
  integer(i4), intent(in) :: nati
!
  lxatom = .false.
  if (nati.eq.17.or.nati.eq.35.or.nati.eq.53.or. &  ! X in A-X...B
      nati.eq.16.or.nati.eq.34.or.nati.eq.52.or. &
      nati.eq.15.or.nati.eq.33.or.nati.eq.51) lxatom = .true.
  end function lxatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculates the inversion angle  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine pgfnff_improper(ndim,numat,nnbr,maxnbr,rnbr,xnbr,ynbr,znbr,i,n1,n2,n3,phi)
!
!  Compute improper torsion angle
!
!  Uses vectors from the neighbour list to ensure correct images of atoms
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: n1
  integer(i4), intent(in)  :: n2
  integer(i4), intent(in)  :: n3
  integer(i4), intent(in)  :: nnbr(numat)
  integer(i4), intent(in)  :: maxnbr
  real(dp),    intent(in)  :: rnbr(maxnbr,numat)
  real(dp),    intent(in)  :: xnbr(maxnbr,numat)
  real(dp),    intent(in)  :: ynbr(maxnbr,numat)
  real(dp),    intent(in)  :: znbr(maxnbr,numat)
  real(dp),    intent(out) :: phi
!
!  Local variables
!
  logical                  :: ldistOK
  real(dp)                 :: vil(3)
  real(dp)                 :: vjk(3)
  real(dp)                 :: vji(3)
  real(dp)                 :: vn(3)
  real(dp)                 :: ril
  real(dp)                 :: ril2
  real(dp)                 :: rn
  real(dp)                 :: rn2
  real(dp)                 :: rrnil
  real(dp)                 :: rnv
!
!  Check validity of n1, n2 and n3
!
  if (n1.gt.nnbr(i).or.n2.gt.nnbr(i).or.n3.gt.nnbr(i)) then
    call pgfnff_error('invalid arguments for bonds to i in pgfnff_improper',0_i4)
    stop
  endif
!
!  j -> i vector
!
  vji(1) = - xnbr(n1,i)
  vji(2) = - ynbr(n1,i)
  vji(3) = - znbr(n1,i)
!
!  j -> k vector
!
  vjk(1) = xnbr(n2,i) - xnbr(n1,i)
  vjk(2) = ynbr(n2,i) - ynbr(n1,i)
  vjk(3) = znbr(n2,i) - znbr(n1,i)
!
!  i -> l vector
!
  vil(1) = xnbr(n3,i)
  vil(2) = ynbr(n3,i)
  vil(3) = znbr(n3,i)
!
!  Compute cross product of rji and rjk
!
  vn(1) = vji(2)*vjk(3) - vji(3)*vjk(2)
  vn(2) = vji(3)*vjk(1) - vji(1)*vjk(3)
  vn(3) = vji(1)*vjk(2) - vji(2)*vjk(1)
!
  rn2 = vn(1)**2 + vn(2)**2 + vn(3)**2
  rn = sqrt(rn2)
!
  ril = rnbr(n3,i)
  ril2 = ril**2
!
  ldistOK = (rn*ril.gt.1.0d-12)
  if (ldistOK) then
    rrnil = 1.0_dp/(rn*ril)
    rnv = (vn(1)*vil(1) + vn(2)*vil(2) + vn(3)*vil(3))*rrnil
  else
    rnv = 0.0_dp
  endif
!
  phi = asin(rnv)

  end subroutine pgfnff_improper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculates torsion force field value  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp) function valijklff(numat,nnbr,maxnbr,nbrno,xnbr,ynbr,znbr,i,j,n1,n2,n3)
!
!  Computes phi for k-i-j-l
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: numat
  integer(i4), intent(in) :: i
  integer(i4), intent(in) :: j
  integer(i4), intent(in) :: n1                   ! Number of j in list of i
  integer(i4), intent(in) :: n2                   ! Number of k in list of i
  integer(i4), intent(in) :: n3                   ! Number of l in list of j
  integer(i4), intent(in) :: nnbr(numat)
  integer(i4), intent(in) :: maxnbr
  integer(i4), intent(in) :: nbrno(maxnbr,numat)
  real(dp),    intent(in) :: xnbr(maxnbr,numat)
  real(dp),    intent(in) :: ynbr(maxnbr,numat)
  real(dp),    intent(in) :: znbr(maxnbr,numat)
!
!  Local variables
!
  integer(i4)             :: ic
  real(dp)                :: eps
  real(dp),     external  :: getangle
  real(dp)                :: nan
  real(dp)                :: nbn
  real(dp)                :: na(3)
  real(dp)                :: nb(3)
  real(dp)                :: ra(3)
  real(dp)                :: rb(3)
  real(dp)                :: rc(3)
  real(dp)                :: snanb
  real(dp),     external  :: vecnorm
!
  parameter (eps=1.0d-14)
!
!  Check validity of n1, n2 and n3
!
  if (n1.gt.nnbr(i).or.n2.gt.nnbr(i).or.n1.lt.1.or.n2.lt.1) then
    call pgfnff_error('invalid arguments for bonds to i in valijklff',0_i4)
    stop
  endif
  if (n3.gt.nnbr(j).or.n3.lt.1) then
    call pgfnff_error('invalid arguments for bond to j in valijklff',0_i4)
    stop
  endif
!
!  Get torsion coordinate
!
!  l -> i vector
!
  ra(1) = - (xnbr(n3,j) + xnbr(n1,i))
  ra(2) = - (ynbr(n3,j) + ynbr(n1,i))
  ra(3) = - (znbr(n3,j) + znbr(n1,i))
!
!  i -> j vector
!
  rb(1) = xnbr(n1,i)
  rb(2) = ynbr(n1,i)
  rb(3) = znbr(n1,i)
!
!  j -> k vector
!
  rc(1) = xnbr(n2,i) - xnbr(n1,i)
  rc(2) = ynbr(n2,i) - ynbr(n1,i)
  rc(3) = znbr(n2,i) - znbr(n1,i)
!
!  Cross products
!
  na(1) = ra(2)*rb(3) - ra(3)*rb(2)
  na(2) = ra(3)*rb(1) - ra(1)*rb(3)
  na(3) = ra(1)*rb(2) - ra(2)*rb(1)
!
  nb(1) = rb(2)*rc(3) - rb(3)*rc(2)
  nb(2) = rb(3)*rc(1) - rb(1)*rc(3)
  nb(3) = rb(1)*rc(2) - rb(2)*rc(1)
!
  nan = vecnorm(na,3_i4,.true.)
  nbn = vecnorm(nb,3_i4,.true.)
!
  snanb = 0.0_dp
  do ic = 1,3
    snanb = snanb + na(ic)*nb(ic)
  enddo
  if (abs(abs(snanb)-1.0_dp).lt.eps) then
    snanb = sign(1.0_dp,snanb)
  endif

  valijklff = acos(snanb)

  end function valijklff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculates angle value  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp) function getangle(ra,rb)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout) :: ra(3)  ! Vector from i to j
  real(dp),    intent(inout) :: rb(3)  ! Vector from i to k
!
!  Local variables
!
  integer(i4)                :: ic
  real(dp)                   :: eps
  real(dp)                   :: rab
  real(dp)                   :: ran
  real(dp)                   :: rbn
  real(dp),       external   :: vecnorm
!
  parameter (eps=1.d-14)
!
!  Normalise vectors
!
  ran = vecnorm(ra,3_i4,.true.)
  rbn = vecnorm(rb,3_i4,.true.)
!
  rab = 0.0_dp
  do ic = 1,3
    rab = rab + ra(ic)*rb(ic)
  enddo
!
  if (abs(abs(rab)-1.0_dp).lt.eps) then
    rab = sign(1.0_dp,rab)
  endif
  getangle = acos(rab)

  end function getangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculates the vector norm      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(dp) function vecnorm(r,n,lnormalise)
!
!  Computes the vector norm and optionalally normalises the vector
!
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  real(dp),    intent(inout) :: r(n)
  integer(i4), intent(in)    :: n
  logical,     intent(in)    :: lnormalise
!
!  Local variables
!
  integer(i4)                :: i
  real(dp)                   :: dotprod
  real(dp)                   :: or
  real(dp)                   :: rn
!
  dotprod = 0.0_dp
  do i = 1,n
    dotprod = dotprod + r(i)*r(i)
  enddo
  rn = sqrt(dotprod)
  if (lnormalise) then
    if (abs(rn).gt.1.d-14) then
      or = 1.0_dp/rn
      do i = 1,n
        r(i) = or*r(i)
      enddo
    endif
  endif
  vecnorm = rn
  end function vecnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find smallest ring in which bond i-j is located  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ringsbond(ndim,numat,i,j,xij,yij,zij,nring,nringatom,nringsize,ringside,rings)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: nring(numat)
  integer(i4), intent(in)  :: nringatom(6,20,numat)
  integer(i4), intent(in)  :: nringsize(20,numat)
  integer(i4), intent(out) :: rings
  real(dp),    intent(in)  :: ringside(3,6,20,numat)
  real(dp),    intent(in)  :: xij
  real(dp),    intent(in)  :: yij
  real(dp),    intent(in)  :: zij
!
!  Local variables
!
  integer(i4)              :: k
  integer(i4)              :: l
  integer(i4)              :: nr
  integer(i4)              :: rings1
  integer(i4)              :: rings2
  real(dp)                 :: diff
!
  rings1 = 99
  rings2 = 99
!
  if (ndim.eq.0) then
!-----------------------------------------------
!  Non-periodic case : Check atom number only  |
!-----------------------------------------------
    do k = 1,nring(i)          ! all rings of atom i
      do l = 1,nringsize(k,i)  ! all atoms of ring k
        if (nringatom(l,k,i).eq.j.and.nringsize(k,i).lt.rings1) then
          rings1 = nringsize(k,i)
        endif
      enddo
    enddo
!
    do k = 1,nring(j)          ! all rings of atom i
      do l = 1,nringsize(k,j)  ! all atoms of ring k
        if (nringatom(l,k,j).eq.i.and.nringsize(k,j).lt.rings2) then
          rings2 = nringsize(k,j)
        endif
      enddo
    enddo
  else
!-------------------------------------------------
!  Periodic case : Check atom number and vector  |
!-------------------------------------------------
!
!  NB: i and j are bonded and so must be first or last in the vector list
!
    do k = 1,nring(i)          ! all rings of atom i
      nr = nringsize(k,i)
      if (nr.lt.rings1) then
        if (nringatom(2,k,i).eq.j) then
          diff = (ringside(1,1,k,i) - xij)**2 + &
                 (ringside(2,1,k,i) - yij)**2 + &
                 (ringside(3,1,k,i) - zij)**2
          if (diff.lt.1.0d-2) then
            rings1 = nr
          elseif (nringatom(nr,k,i).eq.j) then
            diff = (ringside(1,nr,k,i) + xij)**2 + &
                   (ringside(2,nr,k,i) + yij)**2 + &
                   (ringside(3,nr,k,i) + zij)**2
            if (diff.lt.1.0d-2) then
              rings1 = nr
            endif
          endif
        elseif (nringatom(nr,k,i).eq.j) then
          diff = (ringside(1,nr,k,i) + xij)**2 + &
                 (ringside(2,nr,k,i) + yij)**2 + &
                 (ringside(3,nr,k,i) + zij)**2
          if (diff.lt.1.0d-2) then
            rings1 = nr
          endif
        endif
      endif
    enddo
!
    do k = 1,nring(j)          ! all rings of atom j
      nr = nringsize(k,j)
      if (nr.lt.rings2) then
        if (nringatom(2,k,j).eq.i) then
          diff = (ringside(1,1,k,j) + xij)**2 + &
                 (ringside(2,1,k,j) + yij)**2 + &
                 (ringside(3,1,k,j) + zij)**2 
          if (diff.lt.1.0d-2) then
            rings2 = nr
          elseif (nringatom(nr,k,j).eq.i) then
            diff = (ringside(1,nr,k,j) - xij)**2 + &
                   (ringside(2,nr,k,j) - yij)**2 + &
                   (ringside(3,nr,k,j) - zij)**2
            if (diff.lt.1.0d-2) then
              rings2 = nr
            endif
          endif
        elseif (nringatom(nr,k,j).eq.i) then
          diff = (ringside(1,nr,k,j) - xij)**2 + &
                 (ringside(2,nr,k,j) - yij)**2 + &
                 (ringside(3,nr,k,j) - zij)**2
          if (diff.lt.1.0d-2) then
            rings2 = nr
          endif
        endif
      endif
    enddo
  endif
!
  rings = min(rings1,rings2)
  if (rings.eq.99) rings = 0
!
  end subroutine ringsbond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find smallest ring in which atom i is located  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ringsatom(numat,i,nring,nringsize,rings)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: nring(numat)
  integer(i4), intent(in)  :: nringsize(20,numat)
  integer(i4), intent(out) :: rings
!
!  Local variables
!
  integer(i4)              :: k
!
  rings = 99
!
  do k = 1,nring(i)          ! all rings of atom i
    if (nringsize(k,i).lt.rings) rings = nringsize(k,i)
  enddo
  end subroutine ringsatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find smallest ring in which angle i-j-k is located  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ringsbend(ndim,numat,i,j,k,xij,yij,zij,xik,yik,zik,nring,nringsize,nringatom,ringside,rings)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: k
  integer(i4), intent(in)  :: nring(numat)
  integer(i4), intent(in)  :: nringatom(6,20,numat)
  integer(i4), intent(in)  :: nringsize(20,numat)
  integer(i4), intent(out) :: rings
  real(dp),    intent(in)  :: ringside(3,6,20,numat)
  real(dp),    intent(in)  :: xij
  real(dp),    intent(in)  :: yij
  real(dp),    intent(in)  :: zij
  real(dp),    intent(in)  :: xik
  real(dp),    intent(in)  :: yik
  real(dp),    intent(in)  :: zik
!
!  Local variables
!
  integer(i4)              :: itest
  integer(i4)              :: m
  integer(i4)              :: nr
  integer(i4)              :: l
  integer(i4)              :: rings1
  integer(i4)              :: rings2
  integer(i4)              :: rings3
  real(dp)                 :: diff
!
  rings1 = 99
  rings2 = 99
  rings3 = 99
!
  if (ndim.eq.0) then
!-----------------------------------------------
!  Non-periodic case : Check atom number only  |
!-----------------------------------------------
    do m = 1,nring(i)          ! all rings of atom i
      itest = 0
      do l = 1,nringsize(m,i)  ! all atoms of ring m
        if (nringatom(l,m,i).eq.j.or.nringatom(l,m,i).eq.k) itest = itest + 1
      enddo
      if (itest.eq.2.and.nringsize(m,i).lt.rings1) rings1 = nringsize(m,i)
    enddo
!
    do m = 1,nring(j)          ! all rings of atom j
      itest = 0
      do l = 1,nringsize(m,j)  ! all atoms of ring m
        if (nringatom(l,m,j).eq.i.or.nringatom(l,m,j).eq.k) itest = itest + 1
      enddo
      if (itest.eq.2.and.nringsize(m,j).lt.rings2) rings2 = nringsize(m,j)
    enddo
!
    do m = 1,nring(k)          ! all rings of atom k
      itest = 0
      do l = 1,nringsize(m,k)  ! all atoms of ring m
        if (nringatom(l,m,k).eq.i.or.nringatom(l,m,k).eq.j) itest = itest + 1
      enddo
!
! NB - in the original code it had nringsize(m,j) here, but the code below seems more likely to be correct
!
      if (itest.eq.2.and.nringsize(m,k).lt.rings3) rings3 = nringsize(m,k)
    enddo
    rings = min(rings1,rings2,rings3)
  else
!-------------------------------------------------
!  Periodic case : Check atom number and vector  |
!-------------------------------------------------
!
!  NB: i-j and i-k are bonded and so must be first or last in the vector list
!  NB: Because vectors are being checked, just check for rings of i
!
    do m = 1,nring(i)          ! all rings of atom i
      nr = nringsize(m,i)
      if (nr.lt.rings1) then
        itest = 0
        if (nringatom(2,m,i).eq.j) then
          diff = (ringside(1,1,m,i) - xij)**2 + &
                 (ringside(2,1,m,i) - yij)**2 + &
                 (ringside(3,1,m,i) - zij)**2
          if (diff.lt.1.0d-2) then
            itest = itest + 1
          elseif (nringatom(nr,m,i).eq.j) then
            diff = (ringside(1,nr,m,i) + xij)**2 + &
                   (ringside(2,nr,m,i) + yij)**2 + &
                   (ringside(3,nr,m,i) + zij)**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
            endif
          endif
        elseif (nringatom(nr,m,i).eq.j) then
          diff = (ringside(1,nr,m,i) + xij)**2 + &
                 (ringside(2,nr,m,i) + yij)**2 + &
                 (ringside(3,nr,m,i) + zij)**2
          if (diff.lt.1.0d-2) then
            itest = itest + 1
          endif
        endif
!
!  No point checking k if j wasn't found
!
        if (itest.eq.1) then
          if (nringatom(2,m,i).eq.k) then
            diff = (ringside(1,1,m,i) - xik)**2 + &
                   (ringside(2,1,m,i) - yik)**2 + &
                   (ringside(3,1,m,i) - zik)**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
            elseif (nringatom(nr,m,i).eq.k) then
              diff = (ringside(1,nr,m,i) + xik)**2 + &
                     (ringside(2,nr,m,i) + yik)**2 + &
                     (ringside(3,nr,m,i) + zik)**2
              if (diff.lt.1.0d-2) then
                itest = itest + 1
              endif
            endif
          elseif (nringatom(nr,m,i).eq.k) then
            diff = (ringside(1,nr,m,i) + xik)**2 + &
                   (ringside(2,nr,m,i) + yik)**2 + &
                   (ringside(3,nr,m,i) + zik)**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1 
            endif
          endif
        endif
!
        if (itest.eq.2) rings1 = nr
      endif
    enddo
    rings = rings1
  endif
  if (rings.eq.99) rings = 0

  end subroutine ringsbend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find smallest ring in which torsion i-j-k-l is located  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ringstors(ndim,numat,i,j,k,l,vij,vik,vjl,nring,nringsize,nringatom,ringside,rings)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: k
  integer(i4), intent(in)  :: l
  integer(i4), intent(in)  :: nring(numat)
  integer(i4), intent(in)  :: nringatom(6,20,numat)
  integer(i4), intent(in)  :: nringsize(20,numat)
  integer(i4), intent(out) :: rings
  real(dp),    intent(in)  :: ringside(3,6,20,numat)
  real(dp),    intent(in)  :: vij(3)
  real(dp),    intent(in)  :: vik(3)
  real(dp),    intent(in)  :: vjl(3)
!
!  Local variables
!
  integer(i4)              :: ia
  integer(i4)              :: itest
  integer(i4)              :: m
  integer(i4)              :: nr
  integer(i4)              :: rings1
  integer(i4)              :: rings2
  integer(i4)              :: rings3
  integer(i4)              :: rings4
  real(dp)                 :: diff
!
  if (nring(i).eq.0.or.nring(j).eq.0.or.nring(k).eq.0.or.nring(l).eq.0) then
    rings = 0
    return
  endif
!
  rings1 = 99
  rings2 = 99
  rings3 = 99
  rings4 = 99
!
  if (ndim.eq.0) then
!-----------------------------------------------
!  Non-periodic case : Check atom number only  |
!-----------------------------------------------
    do m = 1,nring(i)          ! all rings of atom i
      itest = 0
      do ia = 1,nringsize(m,i)  ! all atoms of ring m
        if (nringatom(ia,m,i).eq.j.or.nringatom(ia,m,i).eq.k.or.nringatom(ia,m,i).eq.l) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,i).lt.rings1) rings1 = nringsize(m,i)
    enddo
!
    do m = 1,nring(j)          ! all rings of atom j
      itest = 0
      do ia = 1,nringsize(m,j)  ! all atoms of ring m
        if (nringatom(ia,m,j).eq.i.or.nringatom(ia,m,j).eq.k.or.nringatom(ia,m,j).eq.l) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,j).lt.rings2) rings2 = nringsize(m,j)
    enddo
!
    do m = 1,nring(k)          ! all rings of atom k
      itest = 0
      do ia = 1,nringsize(m,k)  ! all atoms of ring m
        if (nringatom(ia,m,k).eq.i.or.nringatom(ia,m,k).eq.j.or.nringatom(ia,m,k).eq.l) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,k).lt.rings3) rings3 = nringsize(m,k)
    enddo
!
    do m = 1,nring(l)          ! all rings of atom l
      itest = 0
      do ia = 1,nringsize(m,l)  ! all atoms of ring m
        if (nringatom(ia,m,l).eq.i.or.nringatom(ia,m,l).eq.j.or.nringatom(ia,m,l).eq.k) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,l).lt.rings4) rings4 = nringsize(m,l)
    enddo
    rings = min(rings1,rings2,rings3,rings4)
  else
!-------------------------------------------------
!  Periodic case : Check atom number and vector  |
!-------------------------------------------------
!
!  NB: bonding order is k-i-j-l
!  NB: Because vectors are being checked, just check for rings of i
!
    do m = 1,nring(i)
      nr = nringsize(m,i)
      if (nr.lt.rings1) then
        itest = 0
        if (nringatom(2,m,i).eq.j) then
          diff = (ringside(1,1,m,i) - vij(1))**2 + &
                 (ringside(2,1,m,i) - vij(2))**2 + &
                 (ringside(3,1,m,i) - vij(3))**2
          if (diff.lt.1.0d-2) then
            itest = itest + 1
!
!  i->j found - test j->l 
!         
            if (nringatom(3,m,i).eq.l) then
              diff = (ringside(1,2,m,i) - vjl(1))**2 + &
                     (ringside(2,2,m,i) - vjl(2))**2 + &
                     (ringside(3,2,m,i) - vjl(3))**2
              if (diff.lt.1.0d-2) then
                itest = itest + 1
              endif
            endif  
          elseif (nringatom(nr,m,i).eq.j) then
            diff = (ringside(1,nr,m,i) + vij(1))**2 + &
                   (ringside(2,nr,m,i) + vij(2))**2 + &
                   (ringside(3,nr,m,i) + vij(3))**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
!
!  i->j found - test j->l 
!         
              if (nr.gt.1) then
                if (nringatom(nr-1,m,i).eq.l) then
                  diff = (ringside(1,nr-1,m,i) + vjl(1))**2 + &
                         (ringside(2,nr-1,m,i) + vjl(2))**2 + &
                         (ringside(3,nr-1,m,i) + vjl(3))**2
                  if (diff.lt.1.0d-2) then
                    itest = itest + 1
                  endif
                endif
              endif
            endif
          endif
        elseif (nringatom(nr,m,i).eq.j) then
          diff = (ringside(1,nr,m,i) + vij(1))**2 + &
                 (ringside(2,nr,m,i) + vij(2))**2 + &
                 (ringside(3,nr,m,i) + vij(3))**2
          if (diff.lt.1.0d-2) then
            itest = itest + 1
!
!  i->j found - test j->l
!
            if (nr.gt.1) then
              if (nringatom(nr-1,m,i).eq.l) then
                diff = (ringside(1,nr-1,m,i) + vjl(1))**2 + &
                       (ringside(2,nr-1,m,i) + vjl(2))**2 + &
                       (ringside(3,nr-1,m,i) + vjl(3))**2
                if (diff.lt.1.0d-2) then
                  itest = itest + 1
                endif
              endif
            endif
          endif
        endif
!
!  No point checking k if j and l weren't found
!
        if (itest.eq.2) then
          if (nringatom(2,m,i).eq.k) then
            diff = (ringside(1,1,m,i) - vik(1))**2 + &
                   (ringside(2,1,m,i) - vik(2))**2 + &
                   (ringside(3,1,m,i) - vik(3))**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
            elseif (nringatom(nr,m,i).eq.k) then
              diff = (ringside(1,nr,m,i) + vik(1))**2 + &
                     (ringside(2,nr,m,i) + vik(2))**2 + &
                     (ringside(3,nr,m,i) + vik(3))**2
              if (diff.lt.1.0d-2) then
                itest = itest + 1
              endif
            endif
          elseif (nringatom(nr,m,i).eq.k) then
            diff = (ringside(1,nr,m,i) + vik(1))**2 + &
                   (ringside(2,nr,m,i) + vik(2))**2 + &
                   (ringside(3,nr,m,i) + vik(3))**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
            endif
          endif
        endif
!
        if (itest.eq.3) rings1 = nr
      endif
    enddo
    rings = rings1
  endif
!
  if (rings.eq.99) then
    rings = 0
  endif

  end subroutine ringstors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Find smallest torsion in which angle i-j-k-l is located  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ringstorl(ndim,numat,i,j,k,l,vij,vik,vjl,nring,nringsize,nringatom,ringside,ringl)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)  :: ndim
  integer(i4), intent(in)  :: numat
  integer(i4), intent(in)  :: i
  integer(i4), intent(in)  :: j
  integer(i4), intent(in)  :: k
  integer(i4), intent(in)  :: l
  integer(i4), intent(in)  :: nring(numat)
  integer(i4), intent(in)  :: nringatom(6,20,numat)
  integer(i4), intent(in)  :: nringsize(20,numat)
  integer(i4), intent(out) :: ringl
  real(dp),    intent(in)  :: ringside(3,6,20,numat)
  real(dp),    intent(in)  :: vij(3)
  real(dp),    intent(in)  :: vik(3)
  real(dp),    intent(in)  :: vjl(3)
!
!  Local variables
!
  integer(i4)              :: ia
  integer(i4)              :: itest
  integer(i4)              :: m
  integer(i4)              :: nr
  integer(i4)              :: rings1
  integer(i4)              :: rings2
  integer(i4)              :: rings3
  integer(i4)              :: rings4
  real(dp)                 :: diff
!
  if (nring(i).eq.0.or.nring(j).eq.0.or.nring(k).eq.0.or.nring(l).eq.0) then
    ringl = 0
    return
  endif
!
  rings1 = -99
  rings2 = -99
  rings3 = -99
  rings4 = -99
!
  if (ndim.eq.0) then
!-----------------------------------------------
!  Non-periodic case : Check atom number only  |
!-----------------------------------------------
    do m = 1,nring(i)          ! all rings of atom i
      itest = 0
      do ia = 1,nringsize(m,i)  ! all atoms of ring m
        if (nringatom(ia,m,i).eq.j.or.nringatom(ia,m,i).eq.k.or.nringatom(ia,m,i).eq.l) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,i).gt.rings1) rings1 = nringsize(m,i)
    enddo
!
    do m = 1,nring(j)          ! all rings of atom j
      itest = 0
      do ia = 1,nringsize(m,j)  ! all atoms of ring m
        if (nringatom(ia,m,j).eq.i.or.nringatom(ia,m,j).eq.k.or.nringatom(ia,m,j).eq.l) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,j).gt.rings2) rings2 = nringsize(m,j)
    enddo
!
    do m = 1,nring(k)          ! all rings of atom k
      itest = 0
      do ia = 1,nringsize(m,k)  ! all atoms of ring m
        if (nringatom(ia,m,k).eq.i.or.nringatom(ia,m,k).eq.j.or.nringatom(ia,m,k).eq.l) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,k).gt.rings3) rings3 = nringsize(m,k)
    enddo
!
    do m = 1,nring(l)          ! all rings of atom l
      itest = 0
      do ia = 1,nringsize(m,l)  ! all atoms of ring m
        if (nringatom(ia,m,l).eq.i.or.nringatom(ia,m,l).eq.j.or.nringatom(ia,m,l).eq.k) itest = itest + 1
      enddo
      if (itest.eq.3.and.nringsize(m,l).gt.rings4) rings4 = nringsize(m,l)
    enddo
    ringl = max(rings1,rings2,rings3,rings4)
  else
!-------------------------------------------------
!  Periodic case : Check atom number and vector  |
!-------------------------------------------------
!
!  NB: bonding order is k-i-j-l
!  NB: Because vectors are being checked, just check for rings of i
!
    do m = 1,nring(i)
      nr = nringsize(m,i)
      if (nr.gt.rings1) then
        itest = 0
        if (nringatom(2,m,i).eq.j) then
          diff = (ringside(1,1,m,i) - vij(1))**2 + &
                 (ringside(2,1,m,i) - vij(2))**2 + &
                 (ringside(3,1,m,i) - vij(3))**2
          if (diff.lt.1.0d-2) then
            itest = itest + 1
!
!  i->j found - test j->l 
!         
            if (nringatom(3,m,i).eq.l) then
              diff = (ringside(1,2,m,i) - vjl(1))**2 + &
                     (ringside(2,2,m,i) - vjl(2))**2 + &
                     (ringside(3,2,m,i) - vjl(3))**2
              if (diff.lt.1.0d-2) then
                itest = itest + 1
              endif
            endif  
          elseif (nringatom(nr,m,i).eq.j) then
            diff = (ringside(1,nr,m,i) + vij(1))**2 + &
                   (ringside(2,nr,m,i) + vij(2))**2 + &
                   (ringside(3,nr,m,i) + vij(3))**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
!
!  i->j found - test j->l 
!         
              if (nr.gt.1) then
                if (nringatom(nr-1,m,i).eq.l) then
                  diff = (ringside(1,nr-1,m,i) + vjl(1))**2 + &
                         (ringside(2,nr-1,m,i) + vjl(2))**2 + &
                         (ringside(3,nr-1,m,i) + vjl(3))**2
                  if (diff.lt.1.0d-2) then
                    itest = itest + 1
                  endif
                endif
              endif
            endif
          endif
        elseif (nringatom(nr,m,i).eq.j) then
          diff = (ringside(1,nr,m,i) + vij(1))**2 + &
                 (ringside(2,nr,m,i) + vij(2))**2 + &
                 (ringside(3,nr,m,i) + vij(3))**2
          if (diff.lt.1.0d-2) then
            itest = itest + 1
!
!  i->j found - test j->l
!
            if (nr.gt.1) then
              if (nringatom(nr-1,m,i).eq.l) then
                diff = (ringside(1,nr-1,m,i) + vjl(1))**2 + &
                       (ringside(2,nr-1,m,i) + vjl(2))**2 + &
                       (ringside(3,nr-1,m,i) + vjl(3))**2
                if (diff.lt.1.0d-2) then
                  itest = itest + 1
                endif
              endif
            endif
          endif
        endif
!
!  No point checking k if j and l weren't found
!
        if (itest.eq.2) then
          if (nringatom(2,m,i).eq.k) then
            diff = (ringside(1,1,m,i) - vik(1))**2 + &
                   (ringside(2,1,m,i) - vik(2))**2 + &
                   (ringside(3,1,m,i) - vik(3))**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
            elseif (nringatom(nr,m,i).eq.k) then
              diff = (ringside(1,nr,m,i) + vik(1))**2 + &
                     (ringside(2,nr,m,i) + vik(2))**2 + &
                     (ringside(3,nr,m,i) + vik(3))**2
              if (diff.lt.1.0d-2) then
                itest = itest + 1
              endif
            endif
          elseif (nringatom(nr,m,i).eq.k) then
            diff = (ringside(1,nr,m,i) + vik(1))**2 + &
                   (ringside(2,nr,m,i) + vik(2))**2 + &
                   (ringside(3,nr,m,i) + vik(3))**2
            if (diff.lt.1.0d-2) then
              itest = itest + 1
            endif
          endif
        endif
!
        if (itest.eq.3) rings1 = nr
      endif
    enddo
    ringl = rings1
  endif
!
  if (ringl.eq.-99) then
    ringl = 0
  endif

  end subroutine ringstorl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Function to return row of periodic table up to six  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(i4) function itabrow6(i)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in) :: i
!
  itabrow6 = 0
!
  if (i.gt.0.and.i.le.2) then
    itabrow6 = 1
  elseif (i.gt.2.and.i.le.10) then
    itabrow6 = 2
  elseif (i.gt.10.and.i.le.18) then
    itabrow6 = 3
  elseif (i.gt.18.and.i.le.36) then
    itabrow6 = 4
  elseif (i.gt.36.and.i.le.54) then
    itabrow6 = 5
  elseif (i.gt.54) then
    itabrow6 = 6
  endif

  return
  end function itabrow6

!!!!!!!!!!!!!
!  AlphaCO  !
!!!!!!!!!!!!!
  logical function alphaCO(numat,nat,hyb,nnbr,maxnbr,nbrno,pi,ia,ib)
  use m_pgfnff_types
!
!  Passed variables
!
  integer(i4), intent(in) :: numat
  integer(i4), intent(in) :: nat(numat)
  integer(i4), intent(in) :: hyb(numat)
  integer(i4), intent(in) :: nnbr(numat)
  integer(i4), intent(in) :: maxnbr
  integer(i4), intent(in) :: nbrno(maxnbr,numat)
  integer(i4), intent(in) :: pi(numat)
  integer(i4), intent(in) :: ia
  integer(i4), intent(in) :: ib
!
!  Local variables
!
  integer(i4)             :: j
  integer(i4)             :: ni
  integer(i4)             :: no
!
  alphaCO = .false.
  if (pi(ia).ne.0.and.hyb(ib).eq.3.and.nat(ia).eq.6.and.nat(ib).eq.6) then
    no = 0
    do ni = 1,nnbr(ia)
      j = nbrno(ni,ia)
      if (nat(j).eq.8.and.pi(j).ne.0.and.nnbr(j).eq.1) no = no + 1 ! a pi =O on the C?
    enddo
    if (no.eq.1) then
      alphaCO = .true.
      return
    endif
  endif
  if (pi(ib).ne.0.and.hyb(ia).eq.3.and.nat(ib).eq.6.and.nat(ia).eq.6) then
    no = 0
    do ni = 1,nnbr(ib)
      j = nbrno(ni,ib)
      if (nat(j).eq.8.and.pi(j).ne.0.and.nnbr(j).eq.1) no = no + 1 ! a pi =O on the C?
    enddo
    if (no.eq.1) then
      alphaCO = .true.
      return
    endif
  endif

  end function alphaCO

