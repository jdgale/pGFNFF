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
  subroutine pgfnff_pargen(ndim,kv,numat,nat,nftype,xclat,yclat,zclat,qf,nfrag,nfraglist,qfrag, &
                           nnobo,nobond,nobotyp,maxnbr,nnbr,nbrno,ncnbr,rnbr,xnbr,ynbr, &
                           znbr,nnbr_bond,nbrno_bond,ncnbr_bond,rbnbr,xbnbr,ybnbr,zbnbr,lverbose)
!
!  Initialises system-specific parameters for pGFNFF at start up
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp,    only : pgfnff_initc6
  use m_pgfnff_mrec
  use m_pgfnff_nbr_lib
  use m_pgfnff_topo
  implicit none
!
!  Passed variables
!
  integer(i4),                      intent(in)       :: ndim                       ! Dimensionality of system (0 -> 3)
  integer(i4),                      intent(in)       :: numat                      ! Number of atoms
  integer(i4),                      intent(in)       :: nat(numat)                 ! Atomic number of atoms
  integer(i4),                      intent(in)       :: nftype(numat)              ! Type number of atoms
  integer(i4),                      intent(out)      :: nfrag                      ! Number of fragments
  integer(i4),                      intent(out)      :: nfraglist(numat)           ! Pointer from atom to fragment number
  integer(i4),                      intent(in)       :: nnobo                      ! Number of bonds to be excluded
  integer(i4),                      intent(in)       :: nobond(nnobo)              ! Convolution of atomic numbers for bond exclusion (=1000*n1 + n2;n1<n2)
  integer(i4),                      intent(in)       :: nobotyp(nnobo)             ! Convolution of atom types for bond exclusion (=1000*nt1 + nt2)
  integer(i4),                      intent(in)       :: maxnbr                     ! Maximum number of neighbours in list
  integer(i4),                      intent(in)       :: nnbr(numat)                ! Number of neighbours in list for each atom for each atom
  integer(i4),                      intent(in)       :: nbrno(maxnbr,numat)        ! Pointers to neighbours in list for each atom for each atom
  integer(i4),                      intent(in)       :: ncnbr(maxnbr,numat)        ! Pointers to cell for neighbours in list for each atom for each atom
  integer(i4),                      intent(out)      :: nnbr_bond(numat)           ! Number of bonded neighbours in list for each atom for each atom
  integer(i4),                      intent(out)      :: nbrno_bond(maxnbr,numat)   ! Pointers to bonded neighbours in list for each atom for each atom
  integer(i4),                      intent(out)      :: ncnbr_bond(maxnbr,numat)   ! Pointers to cell for bonded neighbours in list for each atom for each atom
  real(dp),                         intent(in)       :: kv(3,3)                    ! Cartesian reciprocal lattice vectors 
  real(dp),                         intent(in)       :: xclat(numat)               ! Cartesian coordinates in x
  real(dp),                         intent(in)       :: yclat(numat)               ! Cartesian coordinates in y
  real(dp),                         intent(in)       :: zclat(numat)               ! Cartesian coordinates in z
  real(dp),                         intent(inout)    :: qf(numat)                  ! Atomic charges
  real(dp),                         intent(out)      :: qfrag(numat)               ! Fragment charges
  real(dp),                         intent(in)       :: rnbr(maxnbr,numat)         ! Distances to neighbours in full list
  real(dp),                         intent(out)      :: rbnbr(maxnbr,numat)        ! Distances to neighbours in bond list
  real(dp),                         intent(in)       :: xnbr(maxnbr,numat)         ! x component of vectors to neighbours in full list
  real(dp),                         intent(in)       :: ynbr(maxnbr,numat)         ! y component of vectors to neighbours in full list
  real(dp),                         intent(in)       :: znbr(maxnbr,numat)         ! z component of vectors to neighbours in full list
  real(dp),                         intent(out)      :: xbnbr(maxnbr,numat)        ! x component of vectors to neighbours in bond list
  real(dp),                         intent(out)      :: ybnbr(maxnbr,numat)        ! y component of vectors to neighbours in bond list
  real(dp),                         intent(out)      :: zbnbr(maxnbr,numat)        ! z component of vectors to neighbours in bond list
  logical,                          intent(in)       :: lverbose                   ! If true then provide verbose output
!
!  Local variables
!
  integer(i4)                                        :: bbtyp
  integer(i4)                                        :: ctype
  integer(i4)                                        :: hybi
  integer(i4)                                        :: hybj
  integer(i4), dimension(:),       allocatable, save :: hybrid
  integer(i4), dimension(:),       allocatable, save :: itag
  integer(i4)                                        :: i
  integer(i4)                                        :: ia
  integer(i4)                                        :: idum
  integer(i4)                                        :: ifrag
  integer(i4)                                        :: ii
  integer(i4)                                        :: ind
  integer(i4)                                        :: ind3(3)
  integer(i4), dimension(:),       allocatable, save :: imetal
  integer(i4)                                        :: ip
  integer(i4), dimension(:),       allocatable, save :: ipiadr
  integer(i4), dimension(:),       allocatable, save :: ipiadr2
  integer(i4), dimension(:),       allocatable, save :: ipiadr3
  integer(i4), dimension(:),       allocatable, save :: ipiadr4
  integer(i4)                                        :: ipicount
  integer(i4), dimension(:),       allocatable, save :: ipimvec
  integer(i4), dimension(:),       allocatable, save :: ipis
  integer(i4), dimension(:),       allocatable, save :: itmp
  integer(i4)                                        :: iringl
  integer(i4)                                        :: irings
  integer(i4)                                        :: iringsi
  integer(i4)                                        :: iringsj
  integer(i4)                                        :: iringsk
  integer(i4)                                        :: irings4
  integer(i4)                                        :: itabrow6
  integer(i4)                                        :: ix
  integer(i4)                                        :: j
  integer(i4)                                        :: ji
  integer(i4)                                        :: k
  integer(i4)                                        :: kk
  integer(i4)                                        :: l
  integer(i4)                                        :: m
  integer(i4)                                        :: max_ele
  integer(i4)                                        :: mtyp1
  integer(i4)                                        :: mtyp2
  integer(i4)                                        :: nati
  integer(i4)                                        :: natj
  integer(i4)                                        :: natk
  integer(i4)                                        :: nbonds
  integer(i4), dimension(:,:),     allocatable, save :: nbondtype
  integer(i4), dimension(:,:),     allocatable, save :: nbrno_tmp
  integer(i4)                                        :: nc
  integer(i4)                                        :: ncarbo
  integer(i4)                                        :: nchange
  integer(i4)                                        :: nf
  integer(i4)                                        :: nh
  integer(i4)                                        :: nhi
  integer(i4)                                        :: nhj
  integer(i4)                                        :: nheavy
  integer(i4)                                        :: ni
  integer(i4)                                        :: nj
  integer(i4)                                        :: nk
  integer(i4)                                        :: nl
  integer(i4)                                        :: nloop
  integer(i4)                                        :: nm
  integer(i4)                                        :: nmet
  integer(i4)                                        :: nn
  integer(i4)                                        :: nni
  integer(i4)                                        :: nnn
  integer(i4)                                        :: no
  integer(i4)                                        :: npi
  integer(i4)                                        :: npiall
  integer(i4)                                        :: nr
  integer(i4)                                        :: nra(6,20)
  integer(i4)                                        :: nrot
  integer(i4)                                        :: nrs(20)
  integer(i4), dimension(:),       allocatable, save :: nring
  integer(i4), dimension(:,:,:),   allocatable, save :: nringatom
  integer(i4), dimension(:,:),     allocatable, save :: nringsize
  integer(i4), dimension(:),       allocatable, save :: ntopobond
  integer(i4), dimension(:),       allocatable, save :: ntopolast
  integer(i4)                                        :: nsi
  integer(i4)                                        :: ntoposhell
  integer(i4)                                        :: nts
  integer(i4)                                        :: ntypi
  integer(i4)                                        :: ntypj
  integer(i4)                                        :: oldmaxnbr
  integer(i4)                                        :: pis
  integer(i4)                                        :: ierror
  logical                                            :: alphaCO
  logical                                            :: amide
  logical                                            :: amideH
  logical                                            :: bridge
  logical                                            :: pgfnff_excludebond
  logical                                            :: lanypbc
  logical                                            :: lbondok
  logical                                            :: lccij
  logical                                            :: lfound
  logical                                            :: lfrag_charges_known
  logical                                            :: lnofs
  logical                                            :: lnotpicon
  logical                                            :: lpiatom
  logical                                            :: lpicon
  logical                                            :: lpilist
  logical                                            :: lring
  logical                                            :: lsp3ij
  logical                                            :: lsp3kl
  logical                                            :: ltriple
  logical                                            :: lwarn
  logical                                            :: lxatom
  logical,     dimension(:),       allocatable, save :: lcluster
  logical,     dimension(:),       allocatable, save :: lpimpbc
  real(dp),    dimension(:),       allocatable, save :: cn
  real(dp),    dimension(:),       allocatable, save :: ctmp
  real(dp),    dimension(:),       allocatable, save :: dgam
  real(dp),    dimension(:),       allocatable, save :: dqf
  real(dp),    dimension(:),       allocatable, save :: dxi
  real(dp),    dimension(:),       allocatable, save :: mchar
  real(dp),    dimension(:,:),     allocatable, save :: pibond
  real(dp),    dimension(:),       allocatable, save :: qffrag
  real(dp),    dimension(:),       allocatable, save :: qfh
  real(dp),    dimension(:),       allocatable, save :: qtmp
  real(dp),    dimension(:,:),     allocatable, save :: rtopo            ! Topological distances (rtmp)
  real(dp),    dimension(:,:,:,:), allocatable, save :: ringside
  real(dp)                                           :: alpmax
  real(dp)                                           :: bstrength
  real(dp)                                           :: cni
  real(dp)                                           :: cutcn
  real(dp)                                           :: dr0dcni
  real(dp)                                           :: dr0dcnj
  real(dp)                                           :: dum
  real(dp)                                           :: dum1
  real(dp)                                           :: dum2
  real(dp)                                           :: eem_energy
  real(dp)                                           :: ex1
  real(dp)                                           :: ex2
  real(dp)                                           :: f1
  real(dp)                                           :: f2
  real(dp)                                           :: fbsmall
  real(dp)                                           :: fcn
  real(dp)                                           :: fctot
  real(dp)                                           :: feta
  real(dp)                                           :: ff
  real(dp)                                           :: fheavy
  real(dp)                                           :: fij
  real(dp)                                           :: fijk
  real(dp)                                           :: fkl
  real(dp)                                           :: fn
  real(dp)                                           :: fpi
  real(dp)                                           :: fqq
  real(dp)                                           :: fring
  real(dp)                                           :: fsrb2
  real(dp)                                           :: fxh
  real(dp)                                           :: getangle
  real(dp)                                           :: phi
  real(dp)                                           :: pi
  real(dp)                                           :: qfac
  real(dp)                                           :: qmax
  real(dp)                                           :: qmin
  real(dp)                                           :: r0
  real(dp)                                           :: rdiff
  real(dp)                                           :: rij
  real(dp)                                           :: ril2
  real(dp)                                           :: rkl2
  real(dp)                                           :: rlarge
  real(dp)                                           :: rmax
  real(dp)                                           :: rs(3,6,20)
  real(dp)                                           :: sdum3(3)
  real(dp)                                           :: shift
  real(dp)                                           :: sum
  real(dp)                                           :: sumppi
  real(dp)                                           :: theta
  real(dp)                                           :: valijklff
  real(dp)                                           :: vij(3)
  real(dp)                                           :: vik(3)
  real(dp)                                           :: vjl(3)
  real(dp),    dimension(:,:),     allocatable, save :: xtmp
  real(dp),    dimension(:,:),     allocatable, save :: ytmp
  real(dp),    dimension(:,:),     allocatable, save :: ztmp
  real(dp)                                           :: zeta
!
  pi = 4.0_dp*atan(1.0_dp)
!
!  Check that no element exceeds the maximum supported by GFNFF
!
  max_ele = 0
  do i = 1,numat
    max_ele = max(max_ele,nat(i))
  enddo
!
  if (max_ele.gt.max_gfnff_ele) then
    call pgfnff_error('Element exceeds maximum atomic number allowed for GFNFF',0_i4)
    stop
  endif
!**************************
!  Memory initialisation  *
!**************************
!
!  Check sizes of GFNFF arrays relating to atoms in modules
!
  maxat_pgfnff = numat
  call changemaxat_pgfnff
!
!  Set neighbour list arrays to the correct size
!
  maxnbr_lib = maxnbr
  call changemaxnbrlib
!
  if (lverbose) then
    write(ioout,'(/,''  GFNFF Setup Information: '',/)')
  endif
!
!  Set large distance for topology
!
  if (lgfnff_newtopo) then
    rlarge = 1.0d12
  else
    rlarge = 13.0_dp
  endif
!
!  Set Wolf sum parameters for topology if needed
!
  if (lgfnff_topowolf) then
    gfnff_wolf_self = derfc(gfnff_wolf_eta*tdist_thr)/tdist_thr
  endif
!
!  Allocate local memory
!
  allocate(lcluster(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : lcluster '')')
    stop
  endif
  allocate(cn(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : cn '')')
    stop
  endif
  allocate(ctmp(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ctmp '')')
    stop
  endif
  allocate(dgam(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : dgam '')')
    stop
  endif
  allocate(dqf(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : dqf '')')
    stop
  endif
  allocate(dxi(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : dxi '')')
    stop
  endif
  allocate(hybrid(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : hybrid '')')
    stop
  endif
  allocate(ipiadr(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipiadr '')')
    stop
  endif
  allocate(ipiadr2(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipiadr2 '')')
    stop
  endif
  allocate(ipiadr3(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipiadr3 '')')
    stop
  endif
  allocate(ipiadr4(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipiadr4 '')')
    stop
  endif
  allocate(ipimvec(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipimvec '')')
    stop
  endif
  allocate(ipis(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ipis '')')
    stop
  endif
  allocate(itag(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : itag '')')
    stop
  endif
  allocate(itmp(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : itmp '')')
    stop
  endif
  allocate(lpimpbc(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : lpimpbc '')')
    stop
  endif
  allocate(mchar(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : mchar '')')
    stop
  endif
  allocate(nring(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nring '')')
    stop
  endif
  allocate(nringatom(6,20,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nringatom '')')
    stop
  endif
  allocate(nringsize(20,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nringsize '')')
    stop
  endif
  allocate(ringside(3,6,20,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ringside '')')
    stop
  endif
  allocate(imetal(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : imetal '')')
    stop
  endif
  allocate(qffrag(2*numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : qffrag '')')
    stop
  endif
  allocate(qfh(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : qfh '')')
    stop
  endif
  allocate(qtmp(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : qtmp '')')
    stop
  endif
!
  if (lgfnff_xtbtopo) then
    allocate(rtopo(numat,numat),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : rtopo '')')
      stop
    endif
  else
    allocate(ntopobond(numat),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : ntopobond '')')
      stop
    endif
    allocate(ntopolast(numat),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : ntopolast '')')
      stop
    endif
  endif
!################################################################################
!  Pre-loop setup                                                               #
!################################################################################
!
!  Set number of fragments to 0
!
  nfrag = 0
!
!  Initialise saved maxnbr
!
  oldmaxnbr = 0
!
!  Set coordination number to standard value for estimation of radii
!
  do i = 1,numat
    cn(i) = dble(normcn(nat(i)))
  enddo
!
!  Copy atomic charges from qf to qffrag
!
  do i = 1,numat
    qffrag(i) = qf(i)
  enddo
!
!  Initialise the C6 coefficient parameter
!
  call pgfnff_initc6(numat,nat)
!
  cutcn = sqrt(cnthr)
  do i = 1,numat
    call pgfnff_dcn(i,nat,cutcn,cni,sum,maxnbr,nnbr,nbrno,rnbr,xnbr,ynbr,znbr)
!
!  Convert units of sum for XTB-GULP compatibility
!
    sum = sum*xtb_autoaa
!
    cn(i) = cni
    mchar(i) = exp(-0.005_dp*en(nat(i))**8)*sum/(cn(i)+1.0_dp) ! estimated metallic character as ratio of av. dCN and CN
                                                               ! and an EN cut-off function, used in neigbor routine and for BS estimate
  enddo
!################################################################################
!  Loop until setup is converged                                                #
!################################################################################
  nloop = 0
  do while (nloop.lt.2.and.rqshrink.gt.1.0d-3)
!
!  Set coordination number to standard value for estimation of topology
!
    do i = 1,numat
      ctmp(i) = dble(normcn(nat(i)))
    enddo
!
!  Compute numbers of neighbours
!
    call pgfnff_setnbr(numat,nat,ctmp,mchar,qf,maxnbr,nnbr,nbrno,rnbr)
!
!  Check for atoms with no bonds that should have them
!
    lwarn = .false.
    i = 0
    do while (.not.lwarn.and.i.lt.numat)
      i = i + 1
      if (nnbr_full(i).eq.0.and.group(nat(i)).ne.8) then
        lwarn = .true.
      endif
    enddo
    if (lwarn) then
      call pgfnff_warn('some atoms have no bonds where bonds would be expected',0_i4)
    endif
!
!  Set cluster flag
!
    lcluster(1:numat) = .false.
    do i = 1,numat
      if (nat(i).lt.11.and.nnbr_full(i).gt.2) then
        do ni = 1,nnbr_full(i)
          j = nbrno(nbrno_full(ni,i),i)
          if (metal(nat(j)).ne.0.or.nnbr_nohc(j).gt.4) then
            lcluster(i) = .true.
          endif
        enddo
      endif
    enddo
!
!  Set hybridisation state
!
    call pgfnff_hybrid(numat,nat,qf,hybrid,itag,maxnbr,nbrno,rnbr,xnbr,ynbr,znbr)
!
!  Initialise metallic character of atoms
!
    do i = 1,numat
      imetal(i) = metal(nat(i))
      if (nnbr_nohc(i).le.4.and.group(nat(i)).gt.3) imetal(i) = 0 ! Sn,Pb,Bi, with small CN are better described as non-metals
    enddo
    if (oldmaxnbr.lt.maxnbr) then
!
!  If already allocated then free arrays
!
      if (allocated(nbondtype)) then
        deallocate(pibond,stat=ierror)
        if (ierror.gt.0) then
          write(ioout,'('' Error in Memory Deallocation : pibond '')')
          stop
        endif
        deallocate(nbondtype,stat=ierror)
        if (ierror.gt.0) then
          write(ioout,'('' Error in Memory Deallocation : nbondtype '')')
          stop
        endif
      endif
!
!  Allocate pibond and nbondtype
!
      allocate(nbondtype(maxnbr,numat),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : nbondtype '')')
        stop
      endif
      allocate(pibond(maxnbr,numat),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : pibond '')')
        stop
      endif
!
      oldmaxnbr = maxnbr
    endif
!
!  Initialise bonding lists using the list without high coordination numbers and eliminate pointer to main neighbour list
!
    do i = 1,numat
      nnbr_bond(i) = 0
      nati = nat(i)
      ntypi = nftype(i)
      do ni = 1,nnbr_nohc(i)
        nni = nbrno_nohc(ni,i)
        if (nnobo.gt.0) then
          j = nbrno(nni,i)
          natj = nat(j)
          ntypj = nftype(j)
          if (pgfnff_excludebond(nati,ntypi,natj,ntypj,nnobo,nobond,nobotyp)) then
            lbondok = .false.
          else
            lbondok = .true.
          endif
        else
          lbondok = .true.
        endif
        if (lbondok) then
          nnbr_bond(i) = nnbr_bond(i) + 1
          nbrno_bond(nnbr_bond(i),i) = nbrno(nni,i)
          ncnbr_bond(nnbr_bond(i),i) = ncnbr(nni,i)
          rbnbr(nnbr_bond(i),i) = rnbr(nni,i)
          xbnbr(nnbr_bond(i),i) = xnbr(nni,i)
          ybnbr(nnbr_bond(i),i) = ynbr(nni,i)
          zbnbr(nnbr_bond(i),i) = znbr(nni,i)
        endif
      enddo
!
!  Initialise pibond and nbondtype
!
      pibond(1:nnbr_bond(i),i) = 0.0_dp
      nbondtype(1:nnbr_bond(i),i) = 0
    enddo
!****************************************************
! Hueckel setup for all first-row sp2 and sp atoms  *
!****************************************************
    ipiadr(1:numat) = 0
    ipiadr2(1:numat) = 0
!
    k = 0
    piloop: do i = 1,numat
      lpilist = .false.
      if (nat(i).ge.5.and.nat(i).le.9) lpilist = .true.
      if (nat(i).eq.16.or.nat(i).eq.17) lpilist = .true.
      lpiatom = ((hybrid(i).eq.1.or.hybrid(i).eq.2).and.lpilist) ! sp or sp2 and CNOFS
      kk = 0
      do ni = 1,nnbr_bond(i)
        j = nbrno_bond(ni,i)
        if (nat(i).eq.8.and.nat(j).eq.16.and.hybrid(j).eq.5) then
          lpiatom = .false.
          cycle piloop ! SO3 is not a pi
        endif
        if (hybrid(j).eq.1.or.hybrid(j).eq.2) kk = kk + 1   ! attached to sp2 or sp
      enddo
      lnofs = .false.
      if (nat(i).ge.7.and.nat(i).le.9) lnofs = .true.
      if (nat(i).eq.16.or.nat(i).eq.17) lnofs = .true.
      lpicon = (kk.gt.0.and.lnofs)                          ! an N,O,F (sp3) on sp2
      if (nat(i).eq.7.and.nnbr_bond(i).gt.3) cycle piloop   ! NR3-X is not a pi
      if (nat(i).eq.16.and.hybrid(i).eq.5) cycle piloop     ! SO3   is not a pi
      if (lpicon.or.lpiatom) then
        k = k + 1
        ipiadr(k) = i
        ipiadr2(i) = k
      endif
    enddo piloop

    npiall = k
!
!  Build pi neighbour list
!
    nnbr_pi(1:numat) = 0
    piloop2: do i = 1,numat
      if (ipiadr2(i).eq.0) cycle piloop2
      ii = ipiadr2(i)
      nnbr_pi(ii) = 0
      do ni = 1,nnbr_bond(i)
        j = nbrno_bond(ni,i)
        if (ipiadr2(j).gt.0) then
          nnbr_pi(ii) = nnbr_pi(ii) + 1
          nbrno_pi(nnbr_pi(ii),ii) = ipiadr2(j)
          xpnbr(nnbr_pi(ii),ii) = xbnbr(ni,i)
          ypnbr(nnbr_pi(ii),ii) = ybnbr(ni,i)
          zpnbr(nnbr_pi(ii),ii) = zbnbr(ni,i)
        endif
      enddo
    enddo piloop2

    ipimvec(1:npiall) = 0
    lpimpbc(1:npiall) = .false.
!
!  Assign pi atoms to fragments
!
    call pgfnff_mrecgff_pi(npiall,maxnbr,nnbr_pi,nbrno,ipicount,ipimvec,ndim.gt.0,lpimpbc,.false.)
!********************************
!  Setup xi correction for EEQ  *
!********************************
    dxi = 0 ! default none
    xiloop: do i = 1,numat
      nati = nat(i)
      nn = nnbr_bond(i)
      if (nn.eq.0) cycle xiloop
      ip = ipiadr2(i)
      ji = nbrno_bond(1,i) ! first neighbour
      nh = 0
      nm = 0
      do ni = 1,nn
        if (nat(nbrno_bond(ni,i)).eq.1) nh = nh + 1
        if (imetal(nbrno_bond(ni,i)).ne.0) nm = nm + 1
      enddo
!
!  Hydrogen
!
!      if (nati.eq.1.and.nn.gt.1) dxi(i) = dxi(i) - nn*0.01_dp
!
!  Boron
!
      if (nati.eq.5) dxi(i) = dxi(i) + nh*0.015_dp
!
!  Carbon
!
      if (nati.eq.6.and.nn.eq.2.and.itag(i).eq.1) dxi(i) = -0.15_dp ! make carbene more negative
      if (nati.eq.6.and.nn.eq.1.and.nat(ji).eq.8.and.nnbr_bond(ji).eq.1) dxi(ji) = 0.15_dp ! free CO
!
!  Nitrogen
!
!      if (nati.eq.7.and.nn.eq.1.and.nat(ji).eq.6) dxi(i) = 0.0_dp  !CN
!
!  Oxygen / group 6
!
      if (nati.eq.8.and.nn.eq.1.and.ip.ne.0.and.nat(ji).eq.7.and.ipiadr2(ji).ne.0) dxi(i) =  0.05_dp ! nitro oxygen, otherwise NO2 HBs are too strong
      if (nati.eq.8.and.nn.eq.2.and.nh.eq.2)                                       dxi(i) = -0.02_dp ! H2O
      if (group(nati).eq.6.and.nn.gt.2)                                            dxi(i) = dxi(i) + nn*0.005_dp ! good effect
      if (nati.eq.8.or.nati.eq.16)                                                 dxi(i) = dxi(i) - nh*0.005_dp
!
!   Fluorine / group 7
!
      if (group(nati).eq.7.and.nati.gt.9.and.nn.gt.1) then ! polyvalent Cl,Br ...
        if (nm.eq.0) then
          dxi(i) = dxi(i) - nn*0.021_dp ! good effect
        else
          dxi(i) = dxi(i) + nn*0.05_dp ! good effect for TMs
        endif
      endif
!
!  Convert to GULP units
!
      dxi(i) = dxi(i)*xtb_autoev
    enddo xiloop
!
!  Prepare EEQ xi ATOMIC parameter
!  At this point for the non-geom. dep. charges qa with CN = nb
!
    alpmax = 1.0d-12
    do i = 1,numat
      nati = nat(i)
      dum = min(dble(nnbr_bond(i)),cnmax)  ! limits it
!                   base val  spec. corr.    CN dep.
      chieeq(i) = - chi(nati) + dxi(i) + cnf_gfnff(nati)*sqrt(dum)   ! In eV
      cnfeeq(i) = cnf_gfnff(nati)
      gameeq(i) = gam(nati)
      if (imetal(i).eq.2) then            ! the "true" charges for the TM metals are small (for various reasons)
        chieeq(i) = chieeq(i) - mchishift ! so take for the non-geom. dep. ones less electronegative metals yield more q+
      endif                               ! which reflect better the true polarity used for guessing various
                                          ! potential terms. The positive effect of this is big.
      alpeeq(i) = alp(nati)**2
      alpmax = max(alpmax,alp(nati))
    enddo
!
!  Set radius for which erf should be evaluated as it differs from 1
!
    radeeq = 4.0_dp*alpmax
!***************************
!  Topology based charges  *
!***************************
!
!  Determine topology distances by Floyd-Warshall algorithm
!  as they are used in the EEQ to determine qa (approximate topology charges)
!
    if (.not.lgfnff_xtbtopo) then
      maxtopo = maxnbr
      call changemaxtopo
!
!  Build the topological neighbour list
!
      do i = 1,numat
        nnbr_topo(i) = nnbr_bond(i)
        do ni = 1,nnbr_topo(i)
          nbrno_topo(ni,i) = nbrno_bond(ni,i)
          rtnbr(ni,i) = rbnbr(ni,i)
          xtnbr(ni,i) = xbnbr(ni,i)
          ytnbr(ni,i) = ybnbr(ni,i)
          ztnbr(ni,i) = zbnbr(ni,i)
        enddo
      enddo
!
!  Set first shell of neighbours
!
      do i = 1,numat
        do ni = 1,nnbr_topo(i)
          j = nbrno_topo(ni,i)
          rtnbr(ni,i) = rad(nat(i)) + rad(nat(j))
        enddo
      enddo
!
!  Set number of looping shells
!
      if (nloop.eq.0) then
        ntoposhell = maxtoposhell1
      else
        ntoposhell = maxtoposhell
      endif
!
!  Loop over subsequent neighbouring shells
!
      do nts = 1,ntoposhell
!
!  Copy current nnbr_topo to ensure that shells are treated consistently while distances are added
!
        if (nts.eq.1) then
          ntopolast(1:numat) = 0
        else
          ntopolast(1:numat) = ntopobond(1:numat)
        endif
        ntopobond(1:numat) = nnbr_topo(1:numat)
        nchange = 0
!
        do k = 1,numat
          do ni = ntopolast(k)+1,ntopobond(k)
            i = nbrno_topo(ni,k)
            if (rtnbr(ni,k).gt.tdist_thr) cycle
            rmax = tdist_thr - rtnbr(ni,k)
            do nj = 1,ni-1
              j = nbrno_topo(nj,k)
              if (rtnbr(nj,k).gt.rmax) cycle
              rij = rtnbr(ni,k) + rtnbr(nj,k)
              if (rij.gt.tdist_thr) cycle
!
!  Set up vector from i to j
!
              vij(1) = xtnbr(nj,k) - xtnbr(ni,k)
              vij(2) = ytnbr(nj,k) - ytnbr(ni,k)
              vij(3) = ztnbr(nj,k) - ztnbr(ni,k)
!
!  Search for j in i bonding lists
!
              lfound = .false.
              do nk = 1,nnbr_topo(i)
                if (nbrno_topo(nk,i).eq.j) then
                  rdiff = abs(xtnbr(nk,i) - vij(1)) + abs(ytnbr(nk,i) - vij(2)) + abs(ztnbr(nk,i) - vij(3))
                  if (rdiff.lt.1.0d-12) then
                    lfound = .true.
                    if (rtnbr(nk,i).gt.rij) then
                      rtnbr(nk,i) = rij
                      nchange = nchange + 1
                    endif
                  endif
                endif
              enddo
              if (.not.lfound) then
!
!  No i-j distance found and so add to the topological bonds
!
                nnbr_topo(i) = nnbr_topo(i) + 1
                if (nnbr_topo(i).gt.maxtopo) then
                  maxtopo = nnbr_topo(i) + 24
                  call changemaxtopo
                endif
                nchange = nchange + 1
                nbrno_topo(nnbr_topo(i),i) = j
                rtnbr(nnbr_topo(i),i) = rij
                xtnbr(nnbr_topo(i),i) = vij(1)
                ytnbr(nnbr_topo(i),i) = vij(2)
                ztnbr(nnbr_topo(i),i) = vij(3)
              endif
!
!  Search for i in j bonding lists
!
              lfound = .false.
              do nk = 1,nnbr_topo(j)
                if (nbrno_topo(nk,j).eq.i) then
                  rdiff = abs(xtnbr(nk,j) + vij(1)) + abs(ytnbr(nk,j) + vij(2)) + abs(ztnbr(nk,j) + vij(3))
                  if (rdiff.lt.1.0d-12) then
                    lfound = .true.
                    if (rtnbr(nk,j).gt.rij) then
                      rtnbr(nk,j) = rij
                      nchange = nchange + 1
                    endif
                  endif
                endif
              enddo
              if (.not.lfound) then
!
!  No i-j distance found and so add to the topological bonds
!
                nnbr_topo(j) = nnbr_topo(j) + 1
                if (nnbr_topo(j).gt.maxtopo) then
                  maxtopo = nnbr_topo(j) + 24
                  call changemaxtopo
                endif
                nchange = nchange + 1
                nbrno_topo(nnbr_topo(j),j) = i
                rtnbr(nnbr_topo(j),j) = rij
                xtnbr(nnbr_topo(j),j) = - vij(1)
                ytnbr(nnbr_topo(j),j) = - vij(2)
                ztnbr(nnbr_topo(j),j) = - vij(3)
              endif
            enddo
          enddo
        enddo
!
!  If nothing was done then exit the iteration over shells
!
        if (nchange.eq.0) exit
!
!  End of loop over iterations
!
      enddo
!
!  Final scaling of distances
!
      do i = 1,numat
        do ni = 1,nnbr_topo(i)
          if (rtnbr(ni,i).gt.tdist_thr) rtnbr(ni,i) = rlarge
          rtnbr(ni,i) = rfgoed1*rtnbr(ni,i)
        enddo
      enddo
    else
      rtopo(1:numat,1:numat) = rlarge
      do i = 1,numat
        rtopo(i,i) = 0.0_dp
        do ni = 1,nnbr_bond(i)
          j = nbrno_bond(ni,i)
          rtopo(j,i) = rad(nat(i)) + rad(nat(j))
          rtopo(i,j) = rtopo(j,i)
        enddo
      enddo
!
      do k = 1,numat
        do i = 1,numat
          if (rtopo(i,k).gt.tdist_thr) cycle
          do j = 1,numat
            if (rtopo(k,j).gt.tdist_thr) cycle
            if (rtopo(i,j).gt.(rtopo(i,k)+rtopo(k,j))) then
              rtopo(i,j) = rtopo(i,k) + rtopo(k,j)
            endif
          enddo
        enddo
      enddo
!
      do i = 1,numat
        do j = 1,i-1
          if (rtopo(j,i).gt.tdist_thr) rtopo(j,i) = rlarge ! values not properly considered
          rtopo(j,i) = rfgoed1*rtopo(j,i)
          rtopo(i,j) = rtopo(j,i)
        enddo
      enddo
    endif
!
    lfrag_charges_known = .false.
    if (nfrag.le.1) then
!
!  First check for fragments
!
      if (lgfnff_fragment_bond) then
        call pgfnff_mrecgff(numat,maxnbr,nnbr_bond,nbrno_bond,nbrno,nfrag,nfraglist,.false.)
      else
        call pgfnff_mrecgff(numat,maxnbr,nnbr_full,nbrno_full,nbrno,nfrag,nfraglist,.true.)
      endif
!
!  Use input charges to set fragment charges
!
      qfrag(1:nfrag) = 0.0_dp
      do i = 1,numat
        qfrag(nfraglist(i)) = qfrag(nfraglist(i)) + qffrag(i)
      enddo
!
!  Set qffrag elements for nfrag to zero
!
      do i = 1,nfrag
        qffrag(numat+i) = 0.0_dp
      enddo
    endif
!
!  Make an estimate, topology only EEQ charges from rabd values, including "right" fragment charge
!
    if (lgfnff_xtbtopo) then
      call goedeckera(numat,rtopo,nfrag,nfraglist,qfrag,qffrag,eem_energy)
    else
      call goedecker_topo(numat,nfrag,nfraglist,qfrag,qffrag,eem_energy)
    endif
!
!  Copy atomic charges back to qf
!
    do i = 1,numat
      qf(i) = qffrag(i)
    enddo
!
!  Check on largest charges after initial guess
!
    qmax = -1.0d3
    qmin =  1.0d3
    do i = 1,numat
      qmax = max(qmax,qffrag(i))
      qmin = min(qmin,qffrag(i))
    enddo
!
    if (lverbose) then
      write(ioout,'(2x,''Range of initial charges in GFNFF = '',f8.4,'' to '',f8.4)') qmin,qmax
    endif
!
!  Trap excessive charges
!
    if (qmin.lt.-1.0_dp/qrepscal) then
      call pgfnff_error('GFNFF initial charges are too negative ',0_i4)
      stop
    endif
!
!  Estimate how much of the frag charge is on the pi sub systems
!
    if (ipicount.gt.0.and.nloop.gt.0) then
      ipis(1:ipicount) = 0
      qtmp(1:numat) = qffrag(1:numat) ! save the "right" ones
      qfh(1:numat)  = qffrag(1:numat)
      call qheavy(numat,nat,nnbr_bond,maxnbr,nbrno_bond,qfh) ! heavy atoms only ie H condensed to neighbour
      do pis = 1,ipicount ! loop over pi systems
        do k = 1,npiall
          if (ipimvec(k).eq.pis) then
            kk = ipiadr(k)
            ifrag = nfraglist(kk) !the pi atom of this pi fragment is in EEQ fragment ifrag
            exit
          endif
        enddo
        dum2 = qfrag(ifrag) ! save
        qfrag(ifrag) = 0 ! make only this EEQ fragment neutral
        if (lgfnff_xtbtopo) then
          call goedeckera(numat,rtopo,nfrag,nfraglist,qfrag,qffrag,eem_energy) ! for neutral
        else
          call goedecker_topo(numat,nfrag,nfraglist,qfrag,qffrag,eem_energy)
        endif
!
        qfrag(ifrag) = dum2 ! back
        call qheavy(numat,nat,nnbr_bond,maxnbr,nbrno_bond,qffrag) ! heavy atoms only ie H condensed to neighbour
        dqf(1:numat) = qfh(1:numat) - qffrag(1:numat) ! difference charges upon ionization
        dum1 = 0.0_dp
        dum = 0.0_dp
        do k = 1,npiall
          if (ipimvec(k).eq.pis) dum = dum + dqf(ipiadr(k)) ! only pi atoms
        enddo
        dum = dum*1.1_dp !charges tend to be slightly too small 1.1-1.2
        ipis(pis) = idnint(dum)
        dum1 = dum1 + dum
      enddo
      qffrag(1:numat) = qtmp(1:numat) ! put "right" charges used in FF construction and for HB/XB in place
    else
      ipis(1:ipicount) = 0
    endif

    nloop = nloop + 1

!################################################################################
!  End of loop for setup                                                        #
!################################################################################
  enddo

!**********************************************************************
!  Change EEQ J with estimated q which is a kind of third-order term  *
!**********************************************************************
  do i = 1,numat
    ff = 0.0_dp
    if (nat(i).eq.1) ff = -0.08_dp ! H
    if (nat(i).eq.5) ff = -0.05_dp ! B
    if (nat(i).eq.6) then
      ff = -0.27_dp ! C
      if (hybrid(i).lt.3) ff = -0.45_dp ! unsat
      if (hybrid(i).lt.2) ff = -0.34_dp ! unsat
    endif
    if (nat(i).eq.7) then
      ff = -0.13_dp ! N
      if (ipiadr(i).ne.0) ff = -0.14_dp
      if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,i,.false.)) ff = -0.16
    endif
    if (nat(i).eq. 8) then
      ff = -0.15_dp ! O
      if (hybrid(i).lt.3) ff = -0.08_dp ! unsat
    endif
    if (nat(i).eq. 9)                 ff= 0.10_dp ! F
    if (nat(i).gt.10)                 ff=-0.02_dp ! heavy
    if (nat(i).eq.17)                 ff=-0.02_dp ! Cl
    if (nat(i).eq.35)                 ff=-0.11_dp ! Br
    if (nat(i).eq.53)                 ff=-0.07_dp ! I
    if (imetal(i).eq.1)               ff=-0.08_dp ! M maing
    if (imetal(i).eq.2)               ff=-0.9_dp  ! M TM    ??? too large
    if (group(nat(i)).eq.8)           ff= 0.0_dp  ! RG
!
!  Set dgam including changing units to GULP ones
!
    dgam(i) = qffrag(i)*ff/xtb_autoaa
  enddo
!
!  Prepare true EEQ parameter, they are ATOMIC not element specific!
!
  do i = 1,numat
!                base val   spec. corr.
    chieeq(i) = -chi(nat(i)) + dxi(i)
    gameeq(i) =  gam(nat(i)) + dgam(i)
    if (amideH(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr2,i,.false.)) chieeq(i) = chieeq(i) - 0.02_dp
    ff = 0
    if (nat(i).eq.6)        ff = 0.09_dp
    if (nat(i).eq.7)        ff =-0.21_dp
    if (group(nat(i)).eq.6) ff =-0.03_dp
    if (group(nat(i)).eq.7) ff = 0.50_dp
    if (imetal(i).eq.1)     ff = 0.3_dp
    if (imetal(i).eq.2)     ff =-0.1_dp
    alpeeq(i) = (alp(nat(i)) + xtb_autoaa*ff*qffrag(i))**2
  enddo
!***************************************
!  Get ring info (smallest ring size)  *
!***************************************
!
!  Create a version of nbrno_nome so that referencing to nbrno is not passed to getring36
!
  allocate(nbrno_tmp(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_tmp '')')
    stop
  endif
  allocate(xtmp(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xtmp '')')
    stop
  endif
  allocate(ytmp(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ytmp '')')
    stop
  endif
  allocate(ztmp(maxnbr,numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ztmp '')')
    stop
  endif
!
  do i = 1,numat
    do ni = 1,nnbr_nome(i)
      nbrno_tmp(ni,i) = nbrno(nbrno_nome(ni,i),i)
      xtmp(ni,i) = xnbr(nbrno_nome(ni,i),i)
      ytmp(ni,i) = ynbr(nbrno_nome(ni,i),i)
      ztmp(ni,i) = znbr(nbrno_nome(ni,i),i)
    enddo
  enddo
!
  do i = 1,numat
    call getring36(ndim,numat,nnbr_nome,maxnbr,nbrno_tmp,xtmp,ytmp,ztmp,lcluster,i,nra,nrs,rs,nr)
    nring(i) = nr
    if (nr.gt.0) then
      ringside(1:3,1:6,1:nr,i) = rs(1:3,1:6,1:nr)
      nringatom(1:6,1:nr,i) = nra(1:6,1:nr)
      nringsize(1:nr,i) = nrs(1:nr)
    endif
  enddo
!
  deallocate(ztmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ztmp '')')
    stop
  endif
  deallocate(ytmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ytmp '')')
    stop
  endif
  deallocate(xtmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : xtmp '')')
    stop
  endif
  deallocate(nbrno_tmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nbrno_tmp '')')
    stop
  endif
!******************************
!  Non-bonded pair exponents  *
!******************************
  ind = 0
  do i = 1,numat
    nati = nat(i)
    fn = 1.0_dp + nrepscal/(1.0_dp+dble(nnbr_bond(i))**2)
    dum1 = repan(nati)*(1.0_dp + qffrag(i)*qrepscal)*fn ! a small but physically correct decrease of repulsion with q
!
!  Trap error due to negative charges
!
    if (dum1.lt.0.0_dp) then
      call pgfnff_error('GFNFF has gone unphysical due to initial charges',0_i4)
      stop
    endif
    d4_zeta(i) = zeta(nati,qffrag(i))
    do j = 1,i
      ind = ind + 1
      natj = nat(j)
      fn = 1.0_dp + nrepscal/(1.0_dp+dble(nnbr_bond(j))**2)
      dum2 = repan(natj)*(1.0_dp + qffrag(j)*qrepscal)*fn
      ff = 1.0_dp
      if (nati.eq.1.and.natj.eq.1) then
        ff = ff*hhfac                        ! special H ... H case (for other pairs there is no good effect of this)
      endif
      if ((nati.eq.1.and.metal(natj).gt.0).or.(natj.eq.1.and.metal(nati).gt.0)) ff = 0.85_dp ! M...H
      if ((nati.eq.1.and.natj.eq.6).or.(natj.eq.1.and.nati.eq.6))               ff = 0.91_dp ! C...H, good effect
      if ((nati.eq.1.and.natj.eq.8).or.(natj.eq.1.and.nati.eq.8))               ff = 1.04_dp ! O...H, good effect
      alphanb(ind) = sqrt(dum1*dum2)*ff
!
!  Convert units for GULP
!
      alphanb(ind) = alphanb(ind)/(xtb_autoaa**1.5)    ! Distance**(-1.5)
    enddo
  enddo
!***********************************
!  Make list of HB donor basicity  *
!***********************************
!
!  Atom specific (not element) basicity parameters
!
  do i = 1,numat
    nn = nnbr_bond(i)
    nati = nat(i)
    hbbase(i) = xhbas(nati)
!
!  Carbene:
!
    if (nati.eq.6.and.nn.eq.2.and.itag(i).eq.1) hbbase(i) = 1.46_dp
!
!  Carbonyl R-C=O
!
    if (nati.eq.8.and.nn.eq.1) then
      if (nat(nbrno_bond(nn,i)).eq.6) hbbase(i) = 0.68_dp
    endif
!
!  Nitro R-N=O
!
    if (nati.eq.8.and.nn.eq.1) then
      if (nat(nbrno_bond(nn,i)).eq.7) hbbase(i) = 0.47_dp
    endif
  enddo

!**********************************
!  Make list of HB donor acidity  *
!**********************************
!
!  Atom specific (not element) acidity parameters
!
  do i = 1,numat
    nn = nbrno_bond(1,i)
    hbacid(i) = xhaci(nat(i))
!
!  AmideH:
!
    if (amideH(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr2,i,.false.)) hbacid(nn) = hbacid(nn)*0.80_dp
  enddo
!*****************************
!  Make list of ABs for HAB  *
!*****************************
!
!  Initialise memory if this is the first call
!
  if (maxathbH.eq.0) then
    maxathbH = 1
    call changemaxathbH
  endif
  nathbH = 0
  rhbatHl(1:numat) = 0
!
  do i = 1,numat
    if (nat(i).ne.1)  cycle
    if (hybrid(i).eq.1) cycle      ! exclude bridging hydrogens from HB correction
    ff = hqabthr
    j = nbrno_bond(1,i)
    if (j.le.0) cycle
    if (nat(j).gt.10) ff = ff - 0.20_dp                   ! H on heavy atoms may be negatively charged
    if (nat(j).eq.6.and.hybrid(j).eq.3) ff = ff + 0.05_dp ! H on sp3 C must be really positive 0.05
    if (qffrag(i).gt.ff) then                                 ! make list of HB H atoms but only if they have a positive charge
      nathbH = nathbH + 1
      if (nathbH.gt.maxathbH) then
        maxathbH = nathbH + 20
        call changemaxathbH
      endif
      hbatHl(nathbH) = i
      rhbatHl(i) = nathbH
!
!  H charge-scaling term with energy unit conversion
!
      ex1 = exp(hbst*qffrag(i))
      ex2 = ex1 + hbsf
      ABhbq(i) = (ex1/ex2)*xtb_autoev
    endif
  enddo
!
!  Find possible AB atoms
!
  nABat = 0
  do i = 1,numat
!
!  Exclude hydrogens from being AB
!
    if (nat(i).eq.1)  cycle
!
    if (nat(i).eq.6.and.ipiadr2(i).eq.0) cycle ! C sp or sp2 pi
    ff = qabthr
    if (nat(i).gt.10) ff = ff + 0.2_dp   ! heavy atoms may be positively charged
    if (qffrag(i).gt.ff) cycle
    nABat = nABat + 1
    nABatptr(nABat) = i
!
!  AB charge scaling term
!
    ex1 = exp(-hbst*qffrag(i))
    ex2 = ex1 + hbsf
    ABhbq(i) = ex1/ex2
  enddo
!
!  Initialise memory if this is the first call
!
  if (maxatxbAB.eq.0) then
    maxatxbAB = 1
    call changemaxatxbAB
  endif
!
!  Make ABX list
!
  natxbAB = 0
  do i = 1,numat
    do ia = 1,nnbr_bond(i)
      ix = nbrno_bond(ia,i)
      if (lxatom(nat(ix))) then
        if (nat(ix).eq.16.and.nnbr_bond(ix).gt.2) cycle ! no sulphoxide etc S
        do j = 1,numat
!
!  Exclude i = j and ix = j for non-periodic systems : For PBC need to screen images individually
!
          if (ndim.eq.0.and.(i.eq.j.or.j.eq.ix)) cycle
!
!  Now we check for each image to handle PBC
!
          if (xhbas(nat(j)).lt.1.d-6) cycle  ! B must be O,N,...
          if (group(nat(j)).eq.4) then
            if (ipiadr2(j).eq.0.or.qffrag(j).gt.0.05_dp) cycle   ! must be a (pi)base
          endif
          natxbAB = natxbAB + 1
          if (natxbAB.gt.maxatxbAB) then
            maxatxbAB = natxbAB + 20
            call changemaxatxbAB
          endif
          xbatABl(1,natxbAB) = i
          xbatABl(2,natxbAB) = j
          xbatABl(3,natxbAB) = ia
!
!  Halogen charge scaled term
!
          ex1 = exp(xbst*qffrag(ix))
          ex2 = ex1 + xbsf
          ABxbq(1,ix) = ex1/ex2
!
!  B charge scaled term
!
          ex1 = exp(-xbst*qffrag(j))
          ex2 = ex1 + xbsf
          ABxbq(2,j) = ex1/ex2
        enddo
      endif
    enddo
  enddo
!***************
!  Do Hueckel  *
!***************
  if (ipicount.gt.0) then
!
!  Is any pi system periodic?
!
    lanypbc = .false.
    i = 0
    do while (.not.lanypbc.and.i.lt.ipicount)
      i = i + 1
      if (lpimpbc(i)) lanypbc = .true.
    enddo
!
    if (lanypbc) then
      call pgfnff_sethkpt(ndim,kv,lverbose)
      call pgfnff_huckel_pbc(ndim,kv,numat,nat,xclat,yclat,zclat,qf,ipicount,npiall,ipimvec,lpimpbc, &
                             ipiadr,ipis,itag,hybrid,pibond,maxnbr,nnbr_bond,nbrno_bond,rbnbr, &
                             xbnbr,ybnbr,zbnbr,lverbose)
    else
      call pgfnff_huckel(numat,nat,qf,ipicount,npiall,ipimvec,ipiadr,ipis,itag,hybrid,pibond,maxnbr, &
                         nnbr_bond,nbrno_bond,rbnbr,lverbose)
    endif
  endif
!**********************************************
!  Modify hybridisation due to pi assignment  *
!**********************************************
  do i = 1,numat
    if (hybrid(i).eq.2.and.ipiadr(i).eq.0.and.nnbr_bond(i).eq.3.and.group(nat(i)).eq.4) then ! C,Si,Ge... CN=3, no pi
      j = nbrno_bond(1,i)
      k = nbrno_bond(2,i)
      l = nbrno_bond(3,i)
      call pgfnff_improper(ndim,numat,nnbr_bond,maxnbr,rbnbr,xbnbr,ybnbr,zbnbr,i,1_i4,2_i4,3_i4,phi)
      if (abs(phi)*180.0_dp/pi.gt.40.0_dp) hybrid(i) = 3  ! change to sp^3
    endif
  enddo
!
  if (lverbose) then
    write(ioout,'(2x,"atom   neighbours erfCN metchar sp-hybrid imet pi  qest     coordinates")')
    do i = 1,numat
      j = hybrid(i)
      if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,i,.false.)) j = - hybrid(i)
      if (nat(i).eq.6.and.itag(i).eq.1) j = - hybrid(i)
      write(ioout,'(i5,2x,a2,3x,i4,3x,f5.2,2x,f5.2,8x,i2,3x,i2,3x,i2,2x,f6.3,3f12.6)') &
            i,atsym(nat(i)),nnbr_bond(i),cn(i),mchar(i)*xtb_autoaa,j,imetal(i),ipiadr(i), &
            qffrag(i),xclat(i),yclat(i),zclat(i)
    enddo
!
!  Compute fragments and charges for output (check for CT)
!
    write(ioout,'(/,''molecular fragment  # atoms  topo charge'')')
    do i = 1,nfrag
      dum = 0
      m = 0
      do k = 1,numat
        if (nfraglist(k).eq.i) then
          m = m + 1
          dum = dum + qffrag(k)
        endif
      enddo
      write(ioout,'(5x,i3,10x,i4,10x,f8.3)') i,m,dum
    enddo
  endif

!******************************************************************************************
!  Bonds
!******************************************************************************************
  if (lverbose) then
!
!  Compute total number of bonds for output
!
    nbonds = 0
    do i = 1,numat
      nbonds = nbonds + nnbr_bond(i)
    enddo
    nbonds = nbonds/2
!
    write(ioout,*)
    write(ioout,'(10x,"#atoms :",3x,i0)') numat
    write(ioout,'(10x,"#bonds :",3x,i0)') nbonds
    write(ioout,*)
    write(ioout,'(''  NB: R0 values differ from XTB due to way the shift is handled!'')')
    write(ioout,*)
    write(ioout,*) 'bond atoms        type  in ring    R      R0    piBO    fqq  kbond(tot)  alp'
  endif
!
  do i = 1,numat
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
      nati = nat(i)
      natj = nat(j)
      call ringsbond(ndim,numat,i,j,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),nring,nringatom,nringsize,ringside,irings)
      fxh    = 1.0_dp
      fring  = 1.0_dp
      fqq    = 1.0_dp
      fpi    = 1.0_dp
      fheavy = 1.0_dp
      fcn    = 1.0_dp
      fsrb2  = srb2
      bridge = .false.
      shift  = 0.0_dp
!
!  Assign bond type
!
                                                                    nbondtype(ni,i) = 1 ! single
      if (hybrid(i).eq.2.and.hybrid(j).eq.2)                        nbondtype(ni,i) = 2 ! sp2-sp2 = pi
      if (hybrid(i).eq.3.and.hybrid(j).eq.2.and.nati.eq.7)          nbondtype(ni,i) = 2 ! N-sp2
      if (hybrid(j).eq.3.and.hybrid(i).eq.2.and.natj.eq.7)          nbondtype(ni,i) = 2 ! N-sp2
      if (hybrid(i).eq.1.or. hybrid(j).eq.1)                        nbondtype(ni,i) = 3 ! sp-X i.e. no torsion
      if ((group(nati).eq.7.or.nati.eq.1).and.hybrid(i).eq.1) then
        nbondtype(ni,i) = 3 ! linear halogen i.e. no torsion
        bridge = .true.
      endif
      if ((group(natj).eq.7.or.natj.eq.1).and.hybrid(j).eq.1) then
        nbondtype(ni,i) = 3 ! linear halogen i.e. no torsion
        bridge = .true.
      endif
      if (hybrid(i).eq.5.or.hybrid(j).eq.5)                    nbondtype(ni,i) = 4 ! hypervalent
      if (imetal(i).gt.0.or.imetal(j).gt.0)                    nbondtype(ni,i) = 5 ! metal
      if (imetal(i).eq.2.and.imetal(j).eq.2)                   nbondtype(ni,i) = 7 ! TM metal-metal
      if (imetal(j).eq.2.and.itag(i).eq.-1.and.ipiadr(i).gt.0) nbondtype(ni,i) = 6 ! eta
      if (imetal(i).eq.2.and.itag(j).eq.-1.and.ipiadr(j).gt.0) nbondtype(ni,i) = 6 ! eta
      bbtyp = nbondtype(ni,i)
!
!  Normal bond
!
      if (bbtyp.lt.5) then
        hybi = max(hybrid(i),hybrid(j))
        hybj = min(hybrid(i),hybrid(j))
        if (hybi.eq.5.or.hybj.eq.5) then
          bstrength = bstren(4)                                       ! base value hypervalent
        else
          bstrength = bsmat(hybi,hybj)                                ! base value normal hyb
        endif
        if (hybi.eq.3.and.hybj.eq.2.and.(nati.eq.7.or.natj.eq.7)) bstrength = bstren(2)*1.04   ! N-sp2
        if (bridge) then
          if (group(nati).eq.7)             bstrength = bstren(1)*0.50_dp ! bridging X
          if (group(natj).eq.7)             bstrength = bstren(1)*0.50_dp ! bridging X
          if (nati.eq.1.or.nati.eq.9)       bstrength = bstren(1)*0.30_dp ! bridging H/F
          if (natj.eq.1.or.natj.eq.9)       bstrength = bstren(1)*0.30_dp ! bridging H/F
        endif
        if (bbtyp.eq.4)                        shift = hyper_shift        ! hypervalent
        if (nati.eq.1.or.natj.eq.1)            shift = rabshifth          ! XH
        if (nati.eq.9.and.natj.eq.9)           shift = 0.22_dp            ! f2
        if (hybrid(i).eq.3.and.hybrid(j).eq.0) shift = shift - 0.022_dp   ! X-sp3
        if (hybrid(i).eq.0.and.hybrid(j).eq.3) shift = shift - 0.022_dp   ! X-sp3
        if (hybrid(i).eq.1.and.hybrid(j).eq.0) shift = shift + 0.14_dp    ! X-sp
        if (hybrid(i).eq.0.and.hybrid(j).eq.1) shift = shift + 0.14_dp    ! X-sp
        if ((nati.eq.1.and.natj.eq.6)) then
          call ringsatom(numat,j,nring,nringsize,iringsj)
          if (iringsj.eq.3) fxh = 1.05_dp ! 3-ring CH
          if (ctype(numat,nat,nnbr_bond,maxnbr,nbrno_bond,ipiadr,j).eq.1) fxh = 0.95_dp ! aldehyde CH
        endif
        if ((nati.eq.6.and.natj.eq.1)) then
          call ringsatom(numat,i,nring,nringsize,iringsi)
          if (iringsi.eq.3) fxh = 1.05_dp ! 3-ring CH
          if (ctype(numat,nat,nnbr_bond,maxnbr,nbrno_bond,ipiadr,i).eq.1) fxh = 0.95_dp ! aldehyde CH
        endif
        if ((nati.eq.1.and.natj.eq.5)) fxh = 1.10_dp ! BH
        if ((natj.eq.1.and.nati.eq.5)) fxh = 1.10_dp !
        if ((nati.eq.1.and.natj.eq.7)) fxh = 1.06_dp ! NH
        if ((natj.eq.1.and.nati.eq.7)) fxh = 1.06_dp !
        if ((nati.eq.1.and.natj.eq.8)) fxh = 0.93_dp ! OH
        if ((natj.eq.1.and.nati.eq.8)) fxh = 0.93_dp !
        if (bbtyp.eq.3.and.nati.eq.6.and.natj.eq.8) bstrength = bstren(3)*0.90_dp ! makes CO right and M-CO reasonable
        if (bbtyp.eq.3.and.nati.eq.8.and.natj.eq.6) bstrength = bstren(3)*0.90_dp !
!
!  Modify locally for triple bonds
!
        if (bbtyp.eq.3.and.(hybrid(i).eq.0.or.hybrid(j).eq.0)) bbtyp = 1 ! sp-sp3
        if (bbtyp.eq.3.and.(hybrid(i).eq.3.or.hybrid(j).eq.3)) bbtyp = 1 ! sp-sp3
        if (bbtyp.eq.3.and.(hybrid(i).eq.2.or.hybrid(j).eq.2)) bbtyp = 2 ! sp-sp2
!
!  Pi stuff
!
!  NB: Tolerance used here instead of pibond > 0 as in XTB since otherwise rounding errors can trigger noise
!
        if (pibond(ni,i).gt.1.0d-8) then
          shift = hueckelp*(bzref - pibond(ni,i)) ! ref value = no correction is benzene, P=2/3
          if (bbtyp.ne.3.and.pibond(ni,i).gt.0.1) then
            nbondtype(ni,i) = 2
            bbtyp = 2
          endif
          fpi = 1.0_dp - hueckelp2*(bzref2 - pibond(ni,i)) ! deepness
        endif
        if (nati.gt.10.and.natj.gt.10) then
          fcn = fcn/(1.0_dp + 0.007_dp*dble(nnbr_bond(i))**2)
          fcn = fcn/(1.0_dp + 0.007_dp*dble(nnbr_bond(j))**2)
        endif
!
!  Trap potential overflow due to exponential of qfac
!
        qfac = qffrag(i)*qffrag(j)*70.0_dp*15.0_dp
        if (qfac.gt.80_dp) then
          fqq = 1.0_dp
        elseif (qfac.lt.-80_dp) then
          fqq = 1.0_dp + qfacbm0
        else
          fqq = 1.0_dp + qfacbm0*exp(-qfac)/(1.0_dp + exp(-qfac))
        endif
      else
!
!  Metal involved
!
        shift = 0
        bstrength = bstren(bbtyp)
        if (bbtyp.eq.7) then ! TM-TM
          if (itabrow6(nati).gt.4.and.itabrow6(natj).gt.4) bstrength = bstren(8) ! 4/5d-4/5d
          if (itabrow6(nati).eq.4.and.itabrow6(natj).gt.4) bstrength = bstren(9) ! 3d-4/5d
          if (itabrow6(natj).eq.4.and.itabrow6(nati).gt.4) bstrength = bstren(9) ! 3d-4/5d
          dum = 2.0_dp*(mchar(i) + mchar(j))
          dum = min(dum,0.5_dp)  ! limit the "metallic" correction
          bstrength = bstrength*(1.0_dp - dum)
        endif
        mtyp1 = 0  ! no metal
        mtyp2 = 0
        if (group(nati).eq.1)                    mtyp1 = 1  ! Li...
        if (group(nati).eq.2)                    mtyp1 = 2  ! Be...
        if (group(nati).gt.2.and.imetal(i).eq.1) mtyp1 = 3  ! main group
        if (imetal(i).eq.2)                      mtyp1 = 4  ! TM
        if (group(natj).eq.1)                    mtyp2 = 1  ! Li...
        if (group(natj).eq.2)                    mtyp2 = 2  ! Be...
        if (group(natj).gt.2.and.imetal(j).eq.1) mtyp2 = 3  ! main group
        if (imetal(j).eq.2)                      mtyp2 = 4  ! TM
!
!  Trap potential overflow due to exponential of qfac
!
        qfac = qffrag(i)*qffrag(j)*25.0_dp*15.0_dp
        if (qfac.gt.80_dp) then
          dum = 0.0_dp
        elseif (qfac.lt.-80_dp) then
          dum = 1.0_dp
        else
          dum = exp(-qfac)/(1.0_dp + exp(-qfac))
        endif
!
        fqq = 1.0_dp + dum*(qfacbm(mtyp1) + qfacbm(mtyp2))*0.5_dp   ! metal charge corr.
        if (imetal(i).eq.2.and.natj.gt.10)       fheavy = 0.65_dp ! heavy gen. ligand
        if (imetal(j).eq.2.and.nati.gt.10)       fheavy = 0.65_dp
        if (imetal(i).eq.2.and.natj.eq.15)       fheavy = 1.60_dp ! P ligand
        if (imetal(j).eq.2.and.nati.eq.15)       fheavy = 1.60_dp
        if (imetal(i).eq.2.and.group(natj).eq.6) fheavy = 0.85_dp ! chalcogen ligand
        if (imetal(j).eq.2.and.group(nati).eq.6) fheavy = 0.85_dp
        if (imetal(i).eq.2.and.group(natj).eq.7) fheavy = 1.30_dp ! halogen ligand
        if (imetal(j).eq.2.and.group(nati).eq.7) fheavy = 1.30_dp
        if (imetal(i).eq.2.and.natj.eq.1.and.itabrow6(nati).le.5) fxh = 0.80_dp ! hydrogen 3d/4d
        if (imetal(j).eq.2.and.nati.eq.1.and.itabrow6(natj).le.5) fxh = 0.80_dp ! hydrogen 3d/4d
        if (imetal(i).eq.2.and.natj.eq.1.and.itabrow6(nati).gt.5) fxh = 1.00_dp ! hydrogen 5d
        if (imetal(j).eq.2.and.nati.eq.1.and.itabrow6(natj).gt.5) fxh = 1.00_dp ! hydrogen 5d
        if (imetal(i).eq.1.and.natj.eq.1) fxh = 1.20_dp
        if (imetal(j).eq.1.and.nati.eq.1) fxh = 1.20_dp
        if (imetal(j).eq.2.and.hybrid(i).eq.1) then !CO/CN/NC...
          if (nati.eq.6) then
            fpi   =  1.5_dp
            shift = -0.45_dp
          endif
          if (nati.eq.7.and.nnbr_bond(i).ne.1) then
            fpi   = 0.4_dp
            shift = 0.47_dp
          endif
        endif
        if (imetal(i).eq.2.and.hybrid(j).eq.1) then !CO/CN/NC...
          if (natj.eq.6) then
            fpi   = 1.5_dp
            shift = -0.45_dp
          endif
          if (natj.eq.7.and.nnbr_bond(j).ne.1) then
            fpi   = 0.4_dp
            shift = 0.47_dp
          endif
        endif
        if (imetal(i).eq.2) shift = shift + metal2_shift   ! metal shift TM
        if (imetal(j).eq.2) shift = shift + metal2_shift   !
        if (imetal(i).eq.1.and.group(nati).le.2) shift = shift + metal1_shift   ! metal shift group 1+2
        if (imetal(j).eq.1.and.group(natj).le.2) shift = shift + metal1_shift   !
        if (mtyp1.eq.3)                          shift = shift + metal3_shift   ! metal shift MG
        if (mtyp2.eq.3)                          shift = shift + metal3_shift   !
        if (bbtyp.eq.6.and.metal(nati).eq.2)     shift = shift + eta_shift*nnbr_bond(i)! eta coordinated
        if (bbtyp.eq.6.and.metal(natj).eq.2)     shift = shift + eta_shift*nnbr_bond(j)! eta coordinated
        if (mtyp1.gt.0.and.mtyp1.lt.3) fcn = fcn/(1.0_dp + 0.100_dp*dble(nnbr_bond(i))**2)
        if (mtyp2.gt.0.and.mtyp2.lt.3) fcn = fcn/(1.0_dp + 0.100_dp*dble(nnbr_bond(j))**2)
        if (mtyp1.eq.3)                fcn = fcn/(1.0_dp + 0.030_dp*dble(nnbr_bond(i))**2)
        if (mtyp2.eq.3)                fcn = fcn/(1.0_dp + 0.030_dp*dble(nnbr_bond(j))**2)
        if (mtyp1.eq.4)                fcn = fcn/(1.0_dp + 0.036_dp*dble(nnbr_bond(i))**2)
        if (mtyp2.eq.4)                fcn = fcn/(1.0_dp + 0.036_dp*dble(nnbr_bond(j))**2)
        if (mtyp1.eq.4.or.mtyp2.eq.4)then
          fsrb2 = - srb2*0.22_dp ! weaker, inverse EN dep. for TM metals
        else
          fsrb2 = srb2*0.28_dp   ! "normal" for other metals
        endif
      endif
      if (nati.gt.10.and.natj.gt.10) then  ! both atoms are heavy
        shift = shift + hshift3
        if (nati.gt.18) shift = shift + hshift4
        if (natj.gt.18) shift = shift + hshift4
        if (nati.gt.36) shift = shift + hshift5
        if (natj.gt.36) shift = shift + hshift5
      endif
!
!  Set radius shift
!
      vbnbr(1,ni,i) = rabshift + shift   ! value for all bonds + special part
!
!  Convert units of shift to Angstroms
!
      vbnbr(1,ni,i) = vbnbr(1,ni,i)*xtb_autoaa
!
!  Rings prefactor
!
      if (irings.gt.0) fring = 1.0_dp + fringbo*(6.0_dp - dble(irings))**2  ! max ring size is 6
!
!  Steepness
!
      vbnbr(2,ni,i) = srb1*(1.0_dp + fsrb2*(en(nati) - en(natj))**2 + srb3*bstrength)
!
!  Total prefactor        atoms              spec     typ       qterm    heavy-M  pi   XH(3ring,OH...) CN for M
!
      vbnbr(3,ni,i) = - bond(nati)*bond(natj)*fring*bstrength*fqq*fheavy*fpi*fxh*fcn
!
!  Output
!
      r0 = vbnbr(1,ni,i)
      call pgfnff_radij(nat(i),nat(j),cn(i),cn(j),qffrag(i),qffrag(j),r0,.false.,dr0dcni,dr0dcnj,.false.)
      if (lverbose.and.i.ge.j) then
        write(ioout,'(2a2,2i5,2x,2i5,2x,6f8.3)') &
          atsym(nat(i)),atsym(nat(j)),i,j,bbtyp,irings,rbnbr(ni,i),r0,pibond(ni,i),fqq,vbnbr(3,ni,i),vbnbr(2,ni,i)
      endif
!
!  Convert units of parameters to eV and Angstroms
!
      vbnbr(2,ni,i) = vbnbr(2,ni,i)/(xtb_autoaa**2)
      vbnbr(3,ni,i) = vbnbr(3,ni,i)*xtb_autoev
    enddo
  enddo
!******************************************************************************************
!  Bends
!******************************************************************************************

  if (lverbose) then
    write(ioout,*) 'angle atoms        phi0    phi      FC  pi rings'
  endif

  nangle = 0
  do i = 1,numat
    nn = nnbr_bond(i)                 ! take full set to include M-X-Y
    if (nn.le.1) cycle
    if (nn.gt.6) cycle     ! no highly coordinated atom
    nati = nat(i)
    do ni = 2,nn
      j = nbrno_bond(ni,i)
      natj = nat(j)
      do nj = 1,ni-1
        k = nbrno_bond(nj,i)
        natk = nat(k)
        fijk = angl(nati)*angl2(natj)*angl2(natk)
        if (fijk.lt.fcthr) cycle     ! too small
!
        vij(1) = xbnbr(ni,i)
        vij(2) = ybnbr(ni,i)
        vij(3) = zbnbr(ni,i)
!
        vik(1) = xbnbr(nj,i)
        vik(2) = ybnbr(nj,i)
        vik(3) = zbnbr(nj,i)
!
!  Check whether the k-i-j angle is too close to linear to be a valid torsion
!
        phi = getangle(vij,vik)
!
! DEBUG - below this was originally 60 but changed to 65 to avoid issues with FCC metals!!!
        if (metal(nati).gt.0.and.phi*180.0_dp/pi.lt.65.0_dp) cycle ! skip eta cases even if CN < 6 (e.g. CaCp+)
        feta = 1.0_dp
        if (imetal(i).eq.2.and.itag(j).eq.-1.and.ipiadr(j).gt.0) feta = 0.3_dp       ! eta coord.
        if (imetal(i).eq.2.and.itag(k).eq.-1.and.ipiadr(k).gt.0) feta = feta*0.3_dp  !
        nh = 0
        if (natj.eq.1) nh = nh + 1
        if (natk.eq.1) nh = nh + 1
        nnn = 0
        if (natj.eq.7) nnn = nnn + 1
        if (natk.eq.7) nnn = nnn + 1
        no = 0
        if (natj.eq.8) no = no + 1
        if (natk.eq.8) no = no + 1
        nheavy = 0
        if (natj.gt.14) nheavy = nheavy + 1
        if (natk.gt.14) nheavy = nheavy + 1
        nsi = 0
        if (natj.eq.14) nsi = nsi + 1
        if (natk.eq.14) nsi = nsi + 1
        nc = 0
        if (natj.eq.6) nc = nc + 1
        if (natk.eq.6) nc = nc + 1
        nmet = 0
        if (metal(natj).ne.0) nmet = nmet + 1
        if (metal(natk).ne.0) nmet = nmet + 1
        npi = 0
        if (ipiadr(j).ne.0) npi = npi + 1
        if (ipiadr(k).ne.0) npi = npi + 1
!
        nangle = nangle + 1
!
!  Check size of arrays
!
        if (nangle.gt.maxangle) then
          maxangle = nangle + 10
          call changemaxangle
        endif
!
!  Store pointers to atom i and bonds to j and k in neighbour list of i
!
        nangleptr(1,nangle) = i
        nangleptr(2,nangle) = ni
        nangleptr(3,nangle) = nj
!
        call ringsbend(ndim,numat,i,j,k,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),xbnbr(nj,i),ybnbr(nj,i),zbnbr(nj,i), &
                       nring,nringsize,nringatom,ringside,irings)
        ltriple = (hybrid(i).eq.1.or.hybrid(j).eq.1) .or. &
                  (hybrid(i).eq.1.or.hybrid(k).eq.1)
        if (imetal(i).eq.0.and.imetal(j).eq.0.and.imetal(k).eq.0) then
          fqq = 1.0_dp - (qffrag(i)*qffrag(j) + qffrag(i)*qffrag(k))*qfacBEN      ! weaken it
        else
          fqq = 1.0_dp - (qffrag(i)*qffrag(j) + qffrag(i)*qffrag(k))*qfacBEN*2.5_dp
        endif
        f2 = 1.0_dp
        fn = 1.0_dp

!-------------------------
!  Definitions come here
!-------------------------

!!!!!!!!!!!!!
!  Defaults
!!!!!!!!!!!!!
        r0 = 100.0_dp
        if (hybrid(i).eq.1) r0 = 180.0_dp
        if (hybrid(i).eq.2) r0 = 120.0_dp
        if (hybrid(i).eq.3) r0 = 109.5_dp
        if (hybrid(i).eq.3.and.nati.gt.10) then
          if (nn.le.3) r0 = aheavy3    ! heavy maingroup three coordinated
          if (nn.ge.4) r0 = aheavy4    ! heavy maingroup four  coordinated
          if (nn.eq.4.and.group(nati).eq.5)  r0 = 109.5_dp   ! four coordinated group 5
          if (nn.eq.4.and.group(nati).eq.4.and.nati.gt.49) r0 = 109.5_dp   ! four coordinated Sn, Pb
          if (group(nati).eq.4) r0 = r0 - nh*5.0_dp  ! smaller angles for XHn Si...
          if (group(nati).eq.5) r0 = r0 - nh*5.0_dp  ! smaller angles for XHn P..
          if (group(nati).eq.6) r0 = r0 - nh*5.0_dp  ! smaller angles for XHn S..
        endif
        if (hybrid(i).eq.5) then
          r0 = 90.0_dp
          f2 = 0.11_dp    ! not very important
          if (phi*180.0_dp/pi.gt.linthr) r0 = 180.0_dp    ! hypervalent coordination can be linear GEODEP
        endif
!!!!!!!!!!
!  B
!!!!!!!!!!
        if (nati.eq.5) then
          if (hybrid(i).eq.3) r0 = 115.0_dp
          if (hybrid(i).eq.2) r0 = 115.0_dp
        endif
!!!!!!!!!!
!  C cases
!!!!!!!!!!
        if (nati.eq.6)then
          if (hybrid(i).eq.3.and.nh.eq.2) r0 = 108.6_dp  ! CHH
          if (hybrid(i).eq.3.and.no.eq.1) r0 = 108.5_dp  ! COR
          if (hybrid(i).eq.2.and.no.eq.2) r0 = 122.0_dp  ! COO
          if (hybrid(i).eq.2.and.no.eq.1) f2 = 0.7_dp    ! C=O
          if (hybrid(i).eq.1.and.no.eq.2) then
            ltriple = .false.   ! CO2
            f2 = 2.0_dp
          endif
          if (hybrid(i).eq.3.and.nn.gt.4) then
            if (phi*180.0/pi.gt.linthr) r0 = 180.0_dp    ! hypervalent coordination can be linear GEODEP
          endif
        endif
!!!!!!!!!!
!  O cases
!!!!!!!!!!
        if (nati.eq.8.and.nn.eq.2) then
          r0 = 104.5_dp
!
!  H2O
!
          if (nh.eq.2) then
            r0 = 100.0_dp ! compensate ES of the Hs
            f2 = 1.20_dp ! H2O is better with 1.2-1.3 but H2O in fit behaves differently
          endif
          r0 = r0 +  7.0_dp*nsi   ! O angles widen with Si attached
          r0 = r0 + 14.0_dp*nmet  ! O angles widen with M attached
          if (npi.eq.2) then
            r0 = 109.0_dp ! e.g. Ph-O-Ph
          endif
          if (nmet.gt.0.and.phi*180.0_dp/pi.gt.linthr) then
            r0 = 180.0_dp ! metal coordination can be linear GEODEP
            f2 = 0.3_dp
          endif
        endif
!!!!!!!!!!
! N cases
!!!!!!!!!!
        if (nati.eq.7.and.nn.eq.2) then
          f2 = 1.4_dp
          r0 = 115.0_dp
          if (irings.ne.0) r0 = 105.0_dp
          if (natk.eq.8.or.natj.eq.8) r0 = 103.0_dp
          if (natk.eq.9.or.natj.eq.9) r0 = 102.0_dp
          if (hybrid(i).eq.1) r0 = 180.0_dp  ! NC or NNN
          if (imetal(j).eq.2.and.hybrid(i).eq.1.and.natk.eq.7) r0 = 135.0_dp  ! NN on M
          if (imetal(k).eq.2.and.hybrid(i).eq.1.and.natj.eq.7) r0 = 135.0_dp  ! NN on M
        endif
!
!  NR3
!
        if (nati.eq.7.and.hybrid(i).eq.3) then
!
!  in pi system
!
          if (npi.gt.0) then
            if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,i,.false.)) then
              r0 = 115.0_dp
              f2 = 1.2_dp
            else
              sumppi = pibond(ni,i) + pibond(nj,i)
              r0 = 113.
              f2 = 1.0_dp - sumppi*0.7_dp ! must be -!
            endif
          else
            r0 = 104.0_dp ! sat. pyr. N, steep around 106
            f2 = 0.40_dp  ! 1.0 is better for NH3
            f2 = f2 + nh*0.19_dp
            f2 = f2 + no*0.25_dp
            f2 = f2 + nc*0.01_dp
          endif
        endif
!!!!!!!!!!!!!!
!  Ring < 5
!!!!!!!!!!!!!!
        if (irings.eq.3) r0 = 82.0_dp ! 60 gives too little strain
        if (irings.eq.4) r0 = 96.0_dp
        if (irings.eq.5.and.nati.eq.6) r0 = 109.0_dp
!!!!!!!!!!!!!!
!  Specials
!!!!!!!!!!!!!!
!
!  R-X in 3-rings e.g. cyclopropene
!
        if (irings.eq.0) then
          call ringsatom(numat,i,nring,nringsize,iringsi)
          if (iringsi.eq.3) then
            call ringsatom(numat,j,nring,nringsize,iringsj)
            call ringsatom(numat,k,nring,nringsize,iringsk)
            if (iringsj+iringsk.eq.102) r0 = r0 + 4.0_dp
          endif
        endif
!
!  Triple bonds
!
        if (ltriple) then
          f2 = 0.60_dp  ! complex 7 in S30L makes artificial torsions if this is 0.4 which is
                        ! slightly better for the phenylmethylethyne bending pot.
          if (natj.eq.7.or.natk.eq.7) f2 = 1.0_dp
          if ((imetal(j).eq.2.or.imetal(k).eq.2).and.phi*180.0_dp/pi.gt.linthr) then
            if (nati.eq.6.and.natj.eq.6)        f2 = 3.0_dp  ! M-CC
            if (nati.eq.6.and.natk.eq.6)        f2 = 3.0_dp  ! M-CC
            if (nati.eq.6.and.natj.eq.7)        f2 = 3.0_dp  ! M-CN
            if (nati.eq.6.and.natk.eq.7)        f2 = 3.0_dp  ! M-CN
            if (nati.eq.6.and.group(natj).eq.6) f2 = 14.0_dp ! M-CO or CS
            if (nati.eq.6.and.group(natk).eq.6) f2 = 14.0_dp ! M-CO or CS
            if (nati.eq.7.and.natj.eq.7)        f2 = 10.0_dp ! M-NN
            if (nati.eq.7.and.natj.eq.6)        f2 = 10.0_dp ! M-NC
            if (nati.eq.7.and.natk.eq.6)        f2 = 10.0_dp ! M-NC
            if (nati.eq.7.and.natj.eq.8) then
              r0 = 180.0_dp
              f2 = 12.0_dp 
            endif  ! M-NO
            if (nati.eq.7.and.natk.eq.8) then
              r0 = 180.0_dp
              f2 = 12.0_dp
            endif  ! M-NO
          endif
        endif
!
!  Carbene analogues
!
        if (group(nati).eq.4.and.nn.eq.2.and.itag(i).eq.1)  then
          if (nati.eq.6) r0 = 145.0_dp
          if (nati.gt.6) r0 = 90.0_dp
        endif
!
!  SO3X
!
        if (group(nati).eq.6.and.nn.eq.4.and.no.ge.1) r0 = 115.0_dp
!
!  Halogens CN=2
!
        if (group(nati).eq.7.and.hybrid(i).eq.1) then
          if (nati.eq.9)  r0 = 90.0_dp
          if (nati.eq.17) r0 = 90.0_dp
          if (nati.eq.35) r0 = 90.0_dp
          if (nati.eq.53) r0 = 90.0_dp
          if (nati.gt.9.and.phi*180.0_dp/pi.gt.linthr) r0 = 180.0_dp ! change to linear if linear coordinated, GEODEP
          f2 = 0.6_dp/dble(nati)**0.15_dp
        endif
!
!  Pb or Sn can be pyramidal
!
        if (hybrid(i).eq.3.and.group(nati).eq.4.and.nati.gt.32.and.qffrag(i).gt.0.4_dp) then
          if (phi*180.0_dp/pi.gt.140.0_dp) then
            r0 = 180.0_dp ! change to linear
          endif
          if (phi*180.0_dp/pi.lt.100.0_dp) then
            r0 = 90.0_dp
          endif
          f2 = 1.0_dp
        endif
!
!  Metal
!
        if (imetal(i).gt.0) then
          if (hybrid(i).eq.0) then
            r0 = 90.0_dp
            f2 = 1.35_dp  ! important difference to other bends, big effect 1.15,1.25,1.35
          endif
          if (hybrid(i).eq.1) r0 = 180.0_dp
          if (hybrid(i).eq.2) r0 = 120.0_dp
          if (hybrid(i).eq.3) r0 = 109.5_dp
          if (phi*180.0_dp/pi.gt.linthr) r0 = 180.0_dp ! change to linear
        endif

        fn = 1.0_dp - 2.36_dp/dble(nn)**2

!----------------------
! end of definitions
!----------------------
        vangle(1,nangle) = r0*pi/180.0_dp
        fbsmall = (1.0_dp - fbs1*exp(-0.64_dp*(vangle(1,nangle)-pi)**2))
!
!  Central*neigbour charge spec. met.  small angle corr.
!
        vangle(2,nangle) = fijk*fqq*f2*fn*fbsmall*feta
        if (lverbose) then
          write(ioout,'(3i5,2x,3f8.3,l2,i4)') i,j,k,r0,phi*180./pi,vangle(2,nangle),lpicon,irings
        endif
!
!  Convert units of parameters to eV and Angstroms
!
        vangle(2,nangle) = vangle(2,nangle)*xtb_autoev
      enddo
    enddo
  enddo
!
  if (lverbose) then
    write(ioout,'(10x,"#angl  :",3x,i0,/)') nangle
  endif

!******************************************************************************************
!  Torsions
!******************************************************************************************
  if (lverbose) then
    write(ioout,*) 'torsion atoms        nrot   rings    phi0    phi      FC'
  endif

  ntors = 0
  do i = 1,numat
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
!
!  Because bond list has duplicate bonds only do j.le.i
!
      if (j.gt.i) cycle
!
      if (nbondtype(ni,i).eq.3.or.nbondtype(ni,i).eq.6) cycle    ! metal eta or triple
!
      fij = tors(nat(i))*tors(nat(j))             ! atom contribution, central bond
!
      if (fij.lt.fcthr) cycle
!
      if (tors(nat(i)).lt.0.or.tors(nat(j)).lt.0) cycle ! no negative values
      if (metal(nat(i)).gt.1.and.nnbr_bond(i).gt.4)  cycle ! no HC metals
      if (metal(nat(j)).gt.1.and.nnbr_bond(j).gt.4)  cycle !
!
      fqq = 1.0_dp + abs(qffrag(i)*qffrag(j))*qfacTOR      ! weaken it for e.g. CF-CF and similar
!
      call ringsbond(ndim,numat,i,j,xbnbr(ni,i),ybnbr(ni,i),zbnbr(ni,i),nring,nringatom,nringsize,ringside,irings)
      lring = .false.
      lccij  = .false.
      if (irings.gt.0) lring = .true.
      lsp3ij = (hybrid(i).eq.3.and.hybrid(j).eq.3)
      if (nat(i).eq.6.and.nat(j).eq.6) lccij = .true.
      nhi = 1
      nhj = 1
      do nk = 1,nnbr_bond(i)
        if (nat(nbrno_bond(nk,i)).eq.1) nhi = nhi + 1
      enddo
      do nl = 1,nnbr_bond(j)
        if (nat(nbrno_bond(nl,j)).eq.1) nhj = nhj + 1
      enddo
      fij = fij*(dble(nhi)*dble(nhj))**0.07_dp ! n H term
!
      if (gfnff_version.eq.'v6.3.2') then
        if (alphaCO(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,ipiadr,i,j)) fij = fij*1.15_dp
        if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,i,.false.).and. &
            hybrid(j).eq.3.and.nat(j).eq.6) fij = fij*1.15_dp
        if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,j,.false.).and. &
            hybrid(i).eq.3.and.nat(i).eq.6) fij = fij*1.15_dp
      else
        if (alphaCO(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,ipiadr,i,j)) fij = fij*1.3_dp
        if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,i,.false.).and. &
            hybrid(j).eq.3.and.nat(j).eq.6) fij = fij*1.3_dp
        if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,j,.false.).and. &
            hybrid(i).eq.3.and.nat(i).eq.6) fij = fij*1.3_dp
      endif
!
!  Hypervalent
!
      if (gfnff_version.eq.'v6.3.2') then
        if (nbondtype(ni,i).eq.4) fij = fij*0.4_dp ! good effect
      else
        if (nbondtype(ni,i).eq.4) fij = fij*0.2_dp ! good effect
      endif
!
!   Loop over neighbours of ij
!
      do nk = 1,nnbr_bond(i)
!
!  Exclude case where j and k are the same atom
!
        if (nk.eq.ni) cycle
!
        k = nbrno_bond(nk,i)
!
        do nl = 1,nnbr_bond(j)
          l = nbrno_bond(nl,j)
!
!  Exclude case where i and l / k and l are the same atom - for 0D can be done using atom numbers
!
          if (ndim.eq.0) then
            if (l.eq.i) cycle
            if (l.eq.k) cycle
          endif
!
!  Check for angles within the torsion that are too close to linear
!
          vij(1) = xbnbr(ni,i)
          vij(2) = ybnbr(ni,i)
          vij(3) = zbnbr(ni,i)
!
          vik(1) = xbnbr(nk,i)
          vik(2) = ybnbr(nk,i)
          vik(3) = zbnbr(nk,i)
!
          vjl(1) = xbnbr(nl,j)
          vjl(2) = ybnbr(nl,j)
          vjl(3) = zbnbr(nl,j)
!
!  Exclude case where i and l / k and l are the same atom - for PBC case use distances
!
          if (ndim.gt.0) then
            ril2 = (vij(1) + vjl(1))**2 + (vij(2) + vjl(2))**2 + (vij(3) + vjl(3))**2
            if (ril2.lt.1.0d-2) cycle
            rkl2 = (vij(1) + vjl(1) - vik(1))**2 + (vij(2) + vjl(2) - vik(2))**2 + (vij(3) + vjl(3) - vik(3))**2
            if (rkl2.lt.1.0d-2) cycle
          endif
!
!  Check whether the k-i-j angle is too close to linear to be a valid torsion
!
          theta = getangle(vij,vik)
          if (theta*180.0_dp/pi.gt.170.0_dp) cycle 
!
!  Check whether the i-j-l angle is too close to linear to be a valid torsion
!
          vij = - vij
          theta = getangle(vij,vjl)
          vij = - vij
          if (theta*180.0_dp/pi.gt.170.0_dp) cycle 
!
          fkl = tors2(nat(k))*tors2(nat(l))      ! outer kl term
          if (nat(k).eq.7.and.ipiadr(k).eq.0) fkl = tors2(nat(k))*tors2(nat(l))*0.5_dp
          if (nat(l).eq.7.and.ipiadr(l).eq.0) fkl = tors2(nat(k))*tors2(nat(l))*0.5_dp
          if (fkl.lt.fcthr) cycle
          if (tors(nat(k)).lt.0.or.tors(nat(l)).lt.0) cycle ! no negative values
          f1 = torsf(1)
          f2 = 0.0_dp
          fkl = fkl*(dble(nnbr_bond(k))*dble(nnbr_bond(l)))**(-0.14_dp)  ! CN term

!-------------------------
!  Definitions come here |
!-------------------------
          if (lring) then
!
!  Ring
!
            if (irings.gt.3) then
              call ringstors(ndim,numat,i,j,k,l,vij,vik,vjl,nring,nringsize,nringatom,ringside,irings4) ! find smallest ring containing i,j,k,l
            else
              irings4 = 3 ! the 3-ring is special
            endif
            nrot = 1
            if (nbondtype(ni,i).eq.2) nrot = 2 ! max at 90 for pi and symmetric at 0,-180,180
            phi = 0  ! cis
            if (nbondtype(ni,i).eq.1.and.irings4.gt.0) then
              call ringstorl(ndim,numat,i,j,k,l,vij,vik,vjl,nring,nringsize,nringatom,ringside,iringl) ! find largest ring containing i,j,k,l
              lnotpicon = (ipiadr(k).eq.0.and.ipiadr(l).eq.0)  ! do it only for sat. rings
              if (irings4.eq.3.and.lnotpicon) then
                nrot = 1 
                phi = 0.0_dp 
                f1 = fr3
              endif
              if (irings4.eq.4.and.iringl.eq.irings4.and.lnotpicon) then
                nrot = 6
                phi = 30.0_dp
                f1 = fr4
              endif
              if (irings4.eq.5.and.iringl.eq.irings4.and.lnotpicon) then
                nrot = 6
                phi = 30.0_dp
                f1 = fr5
              endif
              if (irings4.eq.6.and.iringl.eq.irings4.and.lnotpicon) then
                nrot = 3
                phi = 60.0_dp
                f1 = fr6
              endif
            endif
            if (irings4.eq.0.and.nbondtype(ni,i).eq.1.and.nnbr_bond(k).eq.1.and.nnbr_bond(l).eq.1) then
              nrot = 6
              phi = 30.0_dp
              f1 = 0.30_dp
            endif
            if (nbondtype(ni,i).eq.2.and.irings.eq.5.and.nat(i)*nat(j).eq.42) then
              if (amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,i,.false.).or. &
                  amide(numat,nat,hybrid,nnbr_bond,maxnbr,nbrno_bond,nbrno,ipiadr,j,.false.)) f1 = 5.0_dp  ! improving CB7
            endif
          else
!
!  Acyclic
!
            phi = 180.0_dp ! trans
            nrot = 1
            if (hybrid(i).eq.3.and.hybrid(j).eq.3) nrot = 3 ! Me case
            if (nbondtype(ni,i).eq.2) nrot = 2 ! max at 90 for pi and symmetric at 0,-180,180
            if (ipiadr(i).gt.0.and.(ipiadr(j).eq.0.and.hybrid(j).eq.3)) then  ! pi-sp3
              f1 = 0.5_dp
              if (nat(i).eq.7) f1 = 0.2_dp ! important for CB7 conf.
              phi = 180.0_dp
              nrot = 3
            endif
            if (ipiadr(j).gt.0.and.(ipiadr(i).eq.0.and.hybrid(i).eq.3)) then
              f1 = 0.5_dp
              if (nat(j).eq.7) f1 = 0.2_dp ! important for CB7 conf.
              phi = 180.0_dp
              nrot = 3
            endif
          endif
!
!  SP3 specials
!
          if (hybrid(i).eq.3.and.hybrid(j).eq.3) then
!
!  N-N, P-P ...
!
            if (group(nat(i)).eq.5.and.group(nat(j)).eq.5) then
              nrot = 3
              phi = 60.0_dp
              f1 = 3.0_dp
            endif
!
!  5-6
!
            if ((group(nat(i)).eq.5.and.group(nat(j)).eq.6).or. &
     &          (group(nat(i)).eq.6.and.group(nat(j)).eq.5)) then
              nrot = 2
              phi = 90.0_dp
              f1 = 1.0_dp
              if (nat(i).ge.15.and.nat(j).ge.15) f1 = 20.0_dp
            endif
!
!  O-O, S-S ...
!
            if (group(nat(i)).eq.6.and.group(nat(j)).eq.6) then
              nrot = 2
              phi = 90.0_dp
              f1 = 5.0_dp
              if (nat(i).ge.16.and.nat(j).ge.16) f1 = 25.0_dp ! better for h2s2
            endif
          endif
!
!  pi system
!
          if (pibond(ni,i).gt.0) then
            f2 = pibond(ni,i)*exp(-2.5_dp*(1.24_dp-pibond(ni,i))**14)  ! decrease to very small values for P < 0.3
                                                                       ! values of 2.5 instead of 2.4 give larger tangles
                                                                       ! the parameter 1.24 is very sensitive ie 1.25 yield 5 deg more in 1,3cB
            if (ipiadr(k).eq.0.and.nat(k).gt.10) f2 = f2*1.3_dp ! the pi BO becomes more significant if heavies are attached
            if (ipiadr(l).eq.0.and.nat(l).gt.10) f2 = f2*1.3_dp
            f1 = f1*0.55_dp
          endif
!
          if (hybrid(k).eq.5.or.hybrid(l).eq.5) fkl = fkl*1.5_dp ! hypervalent corr.
!-----------------------
!  End of definitions  |
!----------------------
!
!  total FC            sigma       pi             charge central outer kl
!
          fctot = (f1 + 10.0_dp*torsf(2)*f2)*fqq*fij*fkl
!
          if (fctot.gt.fcthr) then ! avoid tiny potentials
            ntors = ntors + 1
            if (ntors.gt.maxtors) then
              maxtors = ntors + 10
              call changemaxtors
            endif
            ntorsptr(1,ntors) = i
            ntorsptr(2,ntors) = ni
            ntorsptr(3,ntors) = nk
            ntorsptr(4,ntors) = nl
            ntorsptr(5,ntors) = nrot
            vtors(1,ntors) = phi*pi/180.0_dp
            vtors(2,ntors) = fctot
!
!  Printout
!
            phi = valijklff(numat,nnbr_bond,maxnbr,nbrno_bond,xbnbr,ybnbr,zbnbr,i,j,ni,nk,nl)

            if (lverbose) then
              write(ioout,'(4i5,2x,i2,5x,i2,4x,3f8.3)') &
     &          i,j,k,l,ntorsptr(5,ntors),irings,vtors(1,ntors)*180.0_dp/pi,phi*180.0_dp/pi,vtors(2,ntors)
            endif
          endif
!
!  Extra rot=1 torsion potential for sp3-sp3 to get gauche conf energies well
!
          lsp3kl = (hybrid(k).eq.3.and.hybrid(l).eq.3)
          if (lsp3kl.and.lsp3ij.and.(.not.lring).and.nbondtype(ni,i).lt.5) then
            ntors = ntors + 1
            if (ntors.gt.maxtors) then
              maxtors = ntors + 10
              call changemaxtors
            endif
            ff = torsf(6)
            if (nat(i).eq.7.or.nat(j).eq.7) ff = torsf(7)
            if (nat(i).eq.8.or.nat(j).eq.8) ff = torsf(8)
            ntorsptr(1,ntors) = i
            ntorsptr(2,ntors) = ni
            ntorsptr(3,ntors) = nk
            ntorsptr(4,ntors) = nl
            ntorsptr(5,ntors) = 1
            vtors(1,ntors) = pi
            vtors(2,ntors)= ff * fij * fkl *fqq
            if (lverbose) then
              write(ioout,'(4i5,2x,i2,5x,i2,4x,3f8.3)') &
     &          i,j,k,l,ntorsptr(5,ntors),irings,vtors(1,ntors)*180.0_dp/pi,phi*180.0_dp/pi,vtors(2,ntors)
            endif
          endif
        enddo ! Loop (jneig) over neighours of j -> l
      enddo ! Loop (ineig) over neighours of i -> k
    enddo ! Loop over neighbours of i -> j
  enddo ! Loop over atoms
!
!  Convert units of parameters to eV and Angstroms
!
  do i = 1,ntors
    vtors(2,i) = vtors(2,i)*xtb_autoev
  enddo

!******************************************************************************************
!  Out of planes
!******************************************************************************************
  if (lverbose) then
    write(ioout,*) 'out-of-plane atoms          phi0    phi      FC'
  endif

  do i = 1,numat
    if (nnbr_bond(i).ne.3) cycle
    if (ipiadr(i).eq.0) then
      if (nat(i).ne.7) cycle
    endif
    ntors = ntors + 1
    if (ntors.gt.maxtors) then
      maxtors = ntors + 10
      call changemaxtors
    endif
!
    nj = 1
    nk = 2
    nl = 3
!
!  Sort atoms according to distance to central atom such that the same inversion angle def. always results
!
    sdum3(1) = rbnbr(1,i)
    sdum3(2) = rbnbr(2,i)
    sdum3(3) = rbnbr(3,i)
!
!  Use sort of references to bond numbers to i
!
    ind3(1) = nj
    ind3(2) = nk
    ind3(3) = nl
!
    call ssort(3_i4,sdum3,ind3)
!
    nj = ind3(1)
    nk = ind3(2)
    nl = ind3(3)
!
    j = nbrno_bond(nj,i)
    k = nbrno_bond(nk,i)
    l = nbrno_bond(nl,i)
!
    ntorsptr(1,ntors) = i
    ntorsptr(2,ntors) = nj
    ntorsptr(3,ntors) = nk
    ntorsptr(4,ntors) = nl
!
    if (ipiadr(i).eq.0.and.nat(i).eq.7) then  ! sat N case
      r0 = 80.0_dp
      ff = 0.60_dp
      ntorsptr(5,ntors) = -1
      vtors(1,ntors) = r0*pi/180.0_dp ! double min at +/- phi0
      vtors(2,ntors) = 0.0_dp
      do m = 1,nnbr_bond(i)
        idum = nbrno_bond(m,i)
        vtors(2,ntors) = vtors(2,ntors) + ff*sqrt(repz(nat(idum)))  ! NX3 has higher inv barr. than NH3
      enddo
    else
      ncarbo = 0
      nf = 0
      do m = 1,nnbr_bond(i)
        idum = nbrno_bond(m,i)
        if (nat(idum).eq.8.or.nat(idum).eq.16) ncarbo = ncarbo + 1
        if (group(nat(idum)).eq.7) nf = nf + 1
      enddo
      fqq = 1.0_dp + qffrag(i)*5.0_dp
      ntorsptr(5,ntors) = 0 ! phi0=0 case (pi)
      vtors(1,ntors) = 0.0_dp     !  "      "
      sumppi = pibond(nj,i) + pibond(nk,i) + pibond(nl,i)
      f2 = 1.0_dp - sumppi*torsf(5)
      vtors(2,ntors) = torsf(3)*f2*fqq  ! Base value x piBO x charge term
!
!  Carbonyl correction
!
      if (nat(i).eq.5.and.ncarbo.gt.0)             vtors(2,ntors) = vtors(2,ntors)*38.0_dp
      if (nat(i).eq.6.and.ncarbo.gt.0)             vtors(2,ntors) = vtors(2,ntors)*38.0_dp
      if (nat(i).eq.6.and.nf.gt.0.and.ncarbo.eq.0) vtors(2,ntors) = vtors(2,ntors)*10.0_dp
      if (nat(i).eq.7.and.ncarbo.gt.0)             vtors(2,ntors) = vtors(2,ntors)*10.0_dp/f2 ! no pi dep
    endif
!
!  Printout
!
    call pgfnff_improper(ndim,numat,nnbr_bond,maxnbr,rbnbr,xbnbr,ybnbr,zbnbr,i,nj,nk,nl,phi)
    if (lverbose) then
      write(ioout,'(4i5,7x,3f8.3)') i,j,k,l,vtors(1,ntors)*180.0_dp/pi,phi*180.0_dp/pi,vtors(2,ntors)
    endif
!
!  Convert units of parameters to eV and Angstroms
!
    vtors(2,ntors) = vtors(2,ntors)*xtb_autoev
  enddo

  if (lverbose) then
    write(ioout,'(10x,"#tors  :",3x,i0)') ntors
    write(ioout,'(10x,"#nmol  :",3x,i0)') nfrag
  endif
!
!  Save reference charges for use later
!
  qf0(1:numat) = qffrag(1:numat)
  qf(1:numat)  = qffrag(1:numat)
!
!  All done
!
  if (lverbose) then
    write(ioout,'(10x,"#optfrag :",3x,i0,/)') nfrag
    write(ioout,'(''  Generation of (p)GFN-FF parameters completed'')')
  endif
!
!  Deallocate local memory
!
  deallocate(pibond,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : pibond '')')
    stop
  endif
  deallocate(nbondtype,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nbondtype '')')
    stop
  endif
!
  deallocate(ipis,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipis '')')
    stop
  endif
  if (lgfnff_xtbtopo) then
    deallocate(rtopo,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : rtopo '')')
      stop
    endif
  else
    deallocate(ntopolast,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : ntopolast '')')
      stop
    endif
    deallocate(ntopobond,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : ntopobond '')')
      stop
    endif
  endif
  deallocate(qtmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : qtmp '')')
    stop
  endif
  deallocate(qfh,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : qfh '')')
    stop
  endif
  deallocate(qffrag,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : qffrag '')')
    stop
  endif
  deallocate(imetal,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : imetal '')')
    stop
  endif
  deallocate(ringside,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nringside '')')
    stop
  endif
  deallocate(nringsize,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nringsize '')')
    stop
  endif
  deallocate(nringatom,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nringatom '')')
    stop
  endif
  deallocate(nring,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : nring '')')
    stop
  endif
  deallocate(mchar,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : mchar '')')
    stop
  endif
  deallocate(lpimpbc,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : lpimpbc  '')')
    stop
  endif
  deallocate(itmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : itmp '')')
    stop
  endif
  deallocate(itag,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : itag '')')
    stop
  endif
  deallocate(ipimvec,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipimvec '')')
    stop
  endif
  deallocate(ipiadr4,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipiadr4 '')')
    stop
  endif
  deallocate(ipiadr3,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipiadr3 '')')
    stop
  endif
  deallocate(ipiadr2,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipiadr2 '')')
    stop
  endif
  deallocate(ipiadr,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ipiadr '')')
    stop
  endif
  deallocate(hybrid,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : hybrid '')')
    stop
  endif
  deallocate(dxi,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : dxi '')')
    stop
  endif
  deallocate(dqf,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : dqf '')')
    stop
  endif
  deallocate(dgam,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : dgam '')')
    stop
  endif
  deallocate(ctmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ctmp '')')
    stop
  endif
  deallocate(cn,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : cn '')')
    stop
  endif
  deallocate(lcluster,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : lcluster '')')
    stop
  endif
!
  end subroutine pgfnff_pargen
