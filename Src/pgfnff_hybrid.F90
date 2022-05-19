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
! --------------------------------
!
! Julian Gale, Curtin University, 2020-2022
!
  subroutine pgfnff_hybrid(numat,nat,qf,hyb,itag,maxnbr,nbrno,rnbr,xnbr,ynbr,znbr)
!
!  Sets hybridisation related parameters
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_nbr_lib
  implicit none
!
!  Passed variables
!
  integer(i4),    intent(in)                       :: numat               ! Number of atoms
  integer(i4),    intent(in)                       :: nat(numat)          ! Atomic number of atoms
  integer(i4),    intent(out)                      :: hyb(numat)          ! Hybridisation state
  integer(i4),    intent(out)                      :: itag(numat)         ! Setup tag
  integer(i4),    intent(in)                       :: maxnbr              ! Maximum number of neighbours
  integer(i4),    intent(in)                       :: nbrno(maxnbr,numat) ! Pointer to neighbours
  real(dp),       intent(in)                       :: qf(numat)           ! Charges
  real(dp),       intent(in)                       :: rnbr(maxnbr,numat)  ! Distances to neighbours
  real(dp),       intent(in)                       :: xnbr(maxnbr,numat)  ! x component of vectors to neighbours
  real(dp),       intent(in)                       :: ynbr(maxnbr,numat)  ! y component of vectors to neighbours
  real(dp),       intent(in)                       :: znbr(maxnbr,numat)  ! z component of vectors to neighbours
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: im
  integer(i4)                                      :: j
  integer(i4)                                      :: jmin
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: ll
  integer(i4)                                      :: nati
  integer(i4)                                      :: nb20i
  integer(i4)                                      :: nbdiff
  integer(i4)                                      :: nbmdiff
  integer(i4)                                      :: ncm
  integer(i4)                                      :: nh
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: nm
  integer(i4)                                      :: nn
  integer(i4)                                      :: nni
  integer(i4)                                      :: no
  integer(i4), dimension(:),   allocatable,   save :: nnbr_dum 
  integer(i4), dimension(:,:), allocatable,   save :: nbrno_dum
  integer(i4)                                      :: ierror
  logical                                          :: etacoord
  real(dp)                                         :: cosphi
  real(dp)                                         :: phi
  real(dp)                                         :: radtodeg
  real(dp)                                         :: r212
  real(dp)                                         :: r2i1
  real(dp)                                         :: r2i2
  real(dp)                                         :: ri1(3)
  real(dp)                                         :: ri2(3)
  real(dp)                                         :: rij
  real(dp)                                         :: rmin
  real(dp)                                         :: rtmp
!
  radtodeg = 45.0_dp/atan(1.0_dp)
!
!  Allocate local memory
!
  allocate(nnbr_dum(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_dum '')')
    stop
  endif
  allocate(nbrno_dum(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_dum '')')
    stop
  endif
!
  do i = 1,numat
    nati = nat(i)
    itag(i) = 0
!
!  Important: determine cases where atom is pi bonded to a metal and thus
!  the hyb must be obtained from the reduced (wo metals) neighbor list
!
    etacoord = .false.
    if (nati.le.10) then
      if (nati.eq.6.and.nnbr_full(i).ge.4.and.nnbr_nome(i).eq.3) etacoord = .true.  ! CP case
      if (nati.eq.6.and.nnbr_full(i).eq.3.and.nnbr_nome(i).eq.2) etacoord = .true.  ! alkyne case
      nm = 0
      do ni = 1,nnbr_full(i)  ! how many metals ? and which
        j = nbrno(nbrno_full(ni,i),i)
        if (metal(nat(j)).ne.0) then
          nm = nm + 1
          im = j
        endif
      enddo
      if (nm.eq.0) then
        etacoord = .false.  ! etacoord makes no sense without metals!
      elseif (nm.eq.1) then  ! distinguish M-CR2-R i.e. not an eta coord.
        ncm = 0
        do ni = 1,nnbr_full(i)
          j = nbrno(nbrno_full(ni,i),i)
          if (j.ne.im) then ! all neighbours that are not the metal im
            do nj = 1,nnbr_full(j)
              if (nbrno(nbrno_full(nj,j),j).eq.im) ncm = ncm + 1 ! ncm=1 is alkyne, =2 is cp
            enddo
          endif
        enddo
        if (ncm.eq.0) etacoord = .false.
      endif
    endif
    if (etacoord) then
      itag(i) = -1
      nnbr_dum(i) = nnbr_nome(i)
      nbrno_dum(1:nnbr_nome(i),i) = nbrno_nome(1:nnbr_nome(i),i)
    else
      nnbr_dum(i) = nnbr_full(i)
      nbrno_dum(1:nnbr_full(i),i) = nbrno_full(1:nnbr_full(i),i)
    endif
  enddo
!
  do i = 1,numat
    nati = nat(i)
    hyb(i) = 0
    nbdiff  = nnbr_full(i) - nnbr_nohc(i)
    nbmdiff = nnbr_full(i) - nnbr_nome(i)
    nb20i   = nnbr_dum(i)
    nh = 0
    no = 0
    do ni = 1,nb20i
      if (nat(nbrno(nbrno_dum(ni,i),i)).eq.1) nh = nh + 1
      if (nat(nbrno(nbrno_dum(ni,i),i)).eq.8) no = no + 1
    enddo
!
! H
!
    if (group(nati).eq.1) then
      if (nb20i.eq.2)               hyb(i)=1 ! bridging H
      if (nb20i.gt.2)               hyb(i)=3 ! M+ tetra coord
      if (nb20i.gt.4)               hyb(i)=0 ! M+ HC
    endif
!
! Be
!
    if (group(nati).eq.2) then
      if(nb20i.eq.2)               hyb(i)=1 ! bridging M
      if(nb20i.gt.2)               hyb(i)=3 ! M+ tetra coord
      if(nb20i.gt.4)               hyb(i)=0 !
    endif
!
! B
!
    if (group(nati).eq.3) then
      if (nb20i.gt.4)                                hyb(i)=3
      if (nb20i.gt.4.and.nati.gt.10.and.nbdiff.eq.0) hyb(i)=5
      if (nb20i.eq.4)                                hyb(i)=3
      if (nb20i.eq.3)                                hyb(i)=2
      if (nb20i.eq.2)                                hyb(i)=1
    endif
!
! C
!
    if (group(nati).eq.4) then
      if (nb20i.ge.4)                                hyb(i)=3
      if (nb20i.gt.4.and.nati.gt.10.and.nbdiff.eq.0) hyb(i)=5
      if (nb20i.eq.3)                                hyb(i)=2
      if (nb20i.eq.2) then
!
!  Compute angle between neighbours
!
        ri1(1) = xnbr(nbrno_dum(1,i),i)
        ri1(2) = ynbr(nbrno_dum(1,i),i)
        ri1(3) = znbr(nbrno_dum(1,i),i)
!
        ri2(1) = xnbr(nbrno_dum(2,i),i)
        ri2(2) = ynbr(nbrno_dum(2,i),i)
        ri2(3) = znbr(nbrno_dum(2,i),i)
!
        r2i1 = ri1(1)**2 + ri1(2)**2 + ri1(3)**2
        r2i2 = ri2(1)**2 + ri2(2)**2 + ri2(3)**2
        r212 = (ri2(1)-ri1(1))**2 + (ri2(2)-ri1(2))**2 + (ri2(3)-ri1(3))**2
!
        rtmp = sqrt(r2i1+r2i2+1.0d-14)
        cosphi = 0.5_dp*(r2i1+r2i2-r212)/rtmp
        cosphi = min(cosphi,1.0_dp)
        cosphi = max(cosphi,-1.0_dp)
        phi = acos(cosphi)
!
        if (phi*radtodeg.lt.150.0_dp) then                         ! geometry dep. setup! GEODEP
          hyb(i) = 2  ! otherwise, carbenes will not be recognized
          itag(i) = 1 ! tag for Hueckel and HB routines
        else
          hyb(i) = 1  ! linear triple bond etc
        endif
        if (qf(i).lt.-0.4) then
          hyb(i) = 2
          itag(i) = 0  ! tag for Hueckel and HB routines
        endif
      endif
      if (nb20i.eq.1) hyb(i)=1  ! CO
    endif
!
! N
!
    if (group(nati).eq.5) then
      if (nb20i.ge.4)                                hyb(i) = 3
      if (nb20i.gt.4.and.nati.gt.10.and.nbdiff.eq.0) hyb(i) = 5
      if (nb20i.eq.3)                                hyb(i) = 3
      if (nb20i.eq.3.and.nati.eq.7)  then
        kk = 0
        ll = 0
        nn = 0
        do ni = 1,3
          j = nbrno(nbrno_dum(ni,i),i)
          if (nat(j).eq. 8.and.nnbr_nohc(j).eq.1) kk = kk + 1 ! check for NO2 or R2-N=O
          if (nat(j).eq. 5.and.nnbr_nohc(j).eq.4) ll = ll + 1 ! check for B-N, if the CN(B)=4 the N is loosely bound and sp2
          if (nat(j).eq.16.and.nnbr_nohc(j).eq.4) nn = nn + 1 ! check for N-SO2-
        enddo
        if (nn.eq.1.and.ll.eq.0.and.kk.eq.0)       hyb(i)=3
        if (ll.eq.1.and.nn.eq.0)                   hyb(i)=2
        if (kk.ge.1) then
          hyb(i)=2
          itag(i)=1  ! tag for Hueckel with no el. for the N in NO2
        endif
        if (nbmdiff.gt.0.and.nn.eq.0)              hyb(i)=2  ! pyridin N coord. to heavy atom
      endif
      if (nb20i.eq.2) then
        hyb(i) = 2
!
!  Compute angle between neighbours
!
        ri1(1) = xnbr(nbrno_dum(1,i),i)
        ri1(2) = ynbr(nbrno_dum(1,i),i)
        ri1(3) = znbr(nbrno_dum(1,i),i)
!
        ri2(1) = xnbr(nbrno_dum(2,i),i)
        ri2(2) = ynbr(nbrno_dum(2,i),i)
        ri2(3) = znbr(nbrno_dum(2,i),i)
!
        r2i1 = ri1(1)**2 + ri1(2)**2 + ri1(3)**2
        r2i2 = ri2(1)**2 + ri2(2)**2 + ri2(3)**2
        r212 = (ri2(1)-ri1(1))**2 + (ri2(2)-ri1(2))**2 + (ri2(3)-ri1(3))**2
!
        rtmp = sqrt(r2i1+r2i2+1.0d-14)
        cosphi = 0.5_dp*(r2i1+r2i2-r212)/rtmp
        cosphi = min(cosphi,1.0_dp)
        cosphi = max(cosphi,-1.0_dp)
        phi = acos(cosphi)
!
        j = nbrno(nbrno_dum(1,i),i)
        k = nbrno(nbrno_dum(2,i),i)
        if (nnbr_dum(j).eq.1.and.nat(j).eq.6)     hyb(i) = 1  ! R-N=C
        if (nnbr_dum(k).eq.1.and.nat(k).eq.6)     hyb(i) = 1  ! R-N=C
        if (nnbr_dum(j).eq.1.and.nat(j).eq.7)     hyb(i) = 1  ! R-N=N in e.g. diazomethane
        if (nnbr_dum(k).eq.1.and.nat(k).eq.7)     hyb(i) = 1  ! R-N=N in e.g. diazomethane
        if (nbrno(nbrno_dum(1,i),i).gt.0.and.metal(nat(nbrno(nbrno_dum(1,i),i))).gt.0) hyb(i) = 1 ! M-NC-R in e.g. nitriles
        if (nbrno(nbrno_dum(2,i),i).gt.0.and.metal(nat(nbrno(nbrno_dum(2,i),i))).gt.0) hyb(i) = 1 ! M-NC-R in e.g. nitriles
        if (nat(j).eq.7.and.nat(k).eq.7.and.nnbr_dum(j).le.2.and.nnbr_dum(k).le.2) hyb(i) = 1  ! N=N=N
        if (phi*radtodeg.gt.linthr)               hyb(i)=1  ! geometry dep. setup! GEODEP
      endif
      if (nb20i.eq.1)                             hyb(i)=1
    endif
!
! O
!
    if (group(nati).eq.6) then
      if (nb20i.ge.3)                                hyb(i) = 3
      if (nb20i.gt.3.and.nati.gt.10.and.nbdiff.eq.0) hyb(i) = 5
      if (nb20i.eq.2)                                hyb(i) = 3
      if (nb20i.eq.2.and.nbmdiff.gt.0) then
!
!  Find the CN of the nearest non-metal atom
!
        nn = 0
        rmin = 1.d+42
        jmin = 0
        do ni = 1,nnbr_nohc(i)
          j = nbrno(nbrno_nohc(ni,i),i)
          if (metal(nat(j)).ne.0) cycle
          rij = rnbr(nbrno_nohc(ni,i),i)
          if (rij.lt.rmin) then
            rmin = rij
            jmin = j
          endif
        enddo
!
        if (jmin.gt.0) nn = nnbr_nohc(jmin)
!
        if (nn.eq.3) hyb(i) = 2 ! M-O-X konj
        if (nn.eq.4) hyb(i) = 3 ! M-O-X non
      endif
      if (nb20i.eq.1) hyb(i) = 2
      if (nb20i.eq.1.and.nbdiff.eq.0) then
        if (nnbr_nohc(nbrno(nbrno_nohc(1,i),i)).eq.1) hyb(i) = 1 ! CO
      endif
    endif
!
! F
!
    if (group(nati).eq.7) then
      if (nb20i.eq.2) hyb(i) = 1
      if (nb20i.gt.2.and.nati.gt.10) hyb(i) = 5
    endif
!
! Ne
!
    if (group(nati).eq.8) then
      hyb(i) = 0
      if (nb20i.gt.0.and.nati.gt.2) hyb(i) = 5
    endif
!
!  Done with main groups
!
    if (group(nati).le.0) then ! TMs
      nni = nb20i
      if (nh.ne.0.and.nh.ne.nni) nni=nni-nh ! don't count Hs
      if (nni.le.2)                        hyb(i) = 1
      if (nni.le.2.and.group(nati).le.-6)  hyb(i) = 2
      if (nni.eq.3)                        hyb(i) = 2
      if (nni.eq.4.and.group(nati).gt.-7)  hyb(i) = 3  ! early TM, tetrahedral
      if (nni.eq.4.and.group(nati).le.-7)  hyb(i) = 3  ! late TM, square planar
      if (nni.eq.5.and.group(nati).eq.-3)  hyb(i) = 3  ! Sc-La are tetrahedral CN=5
    endif
  enddo
!
!  Transfer local dummy neighbour list to nohc data structures
!
  nnbr_nohc(1:numat) = nnbr_dum(1:numat)
  do i = 1,numat
    nbrno_nohc(1:nnbr_dum(i),i) = nbrno_dum(1:nnbr_dum(i),i)
  enddo
!
  j = 0
  do i = 1,numat
    if (nnbr_nohc(i).gt.12) j = j + 1
    do ni = 1,nnbr_nohc(i)
      k = nbrno(nbrno_nohc(ni,i),i)
      if (nat(k).eq.6.and.nat(i).eq.6.and.itag(i).eq.1.and.itag(k).eq.1) then ! check the very special situation of
        itag(i) = 0                                                           ! two carbene C bonded which is an arine
        itag(k) = 0
      endif
    enddo
  enddo
!
  if (lgfnff_highcn_trap) then
    if (dble(j)/dble(numat).gt.0.3_dp) then
      call pgfnff_error('too many atoms in GFNFF with extremely high CN',0_i4)
      stop
    endif
  endif
!
!  Free local memory
!
  deallocate(nbrno_dum,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nbrno_dum '')')
    stop
  endif
  deallocate(nnbr_dum,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nnbr_dum '')')
    stop
  endif
!
  return
  end
