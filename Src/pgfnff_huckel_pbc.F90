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
  subroutine pgfnff_huckel_pbc(ndim,kv,numat,nat,xclat,yclat,zclat,qf,ipicount,npiall,ipimvec,lpimpbc,ipiadr,ipis,itag,hybrid, &
                               pibond,maxnbr,nnbr_bond,nbrno_bond,rbnbr,xbnbr,ybnbr,zbnbr,lverbose)
!
!  Periodic boundary conditions Hueckel theory
!
!  Based on original code of Spicher and Grimme for non-periodic Hueckel solve
!
!  Modifications for PBC:
!
!  Julian Gale, Curtin University, July 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                          intent(in)    :: ndim                      ! Number of dimensions
  integer(i4),                          intent(in)    :: numat                     ! Number of atoms
  integer(i4),                          intent(in)    :: nat(numat)                ! Atomic numbers
  integer(i4),                          intent(in)    :: hybrid(numat)             ! Hybridisation state
  integer(i4),                          intent(in)    :: ipicount
  integer(i4),                          intent(inout) :: ipiadr(numat)
  integer(i4),                          intent(in)    :: ipimvec(numat)
  integer(i4),                          intent(in)    :: ipis(numat)
  integer(i4),                          intent(in)    :: itag(numat)
  integer(i4),                          intent(in)    :: maxnbr                    ! Maximum number of neighbours
  integer(i4),                          intent(in)    :: nnbr_bond(numat)          ! Number of neighbours
  integer(i4),                          intent(in)    :: nbrno_bond(maxnbr,numat)  ! Pointers to neighbours
  integer(i4),                          intent(in)    :: npiall
  logical,                              intent(in)    :: lpimpbc(numat)
  logical,                              intent(in)    :: lverbose                  ! If true then provide verbose output
  real(dp),                             intent(in)    :: kv(3,3)                   ! Reciprocal lattice vectors
  real(dp),                             intent(out)   :: pibond(maxnbr,numat)
  real(dp),                             intent(in)    :: qf(numat)                 ! Charges
  real(dp),                             intent(in)    :: rbnbr(maxnbr,numat)       ! Distances for neighbours
  real(dp),                             intent(in)    :: xbnbr(maxnbr,numat)       ! x component of distances for neighbours
  real(dp),                             intent(in)    :: ybnbr(maxnbr,numat)       ! y component of distances for neighbours
  real(dp),                             intent(in)    :: zbnbr(maxnbr,numat)       ! z component of distances for neighbours
  real(dp),                             intent(in)    :: xclat(numat)              ! Cartesian x coordinates
  real(dp),                             intent(in)    :: yclat(numat)              ! Cartesian y coordinates
  real(dp),                             intent(in)    :: zclat(numat)              ! Cartesian z coordinates
!
!  Local variables
!
  integer(i4)                                         :: hybi
  integer(i4)                                         :: i
  integer(i4)                                         :: ia
  integer(i4)                                         :: ifail
  integer(i4)                                         :: ii
  integer(i4)                                         :: ilaenv
  integer(i4),  dimension(:),       allocatable, save :: ipiadr3
  integer(i4),  dimension(:),       allocatable, save :: ipiadr4
  integer(i4),  dimension(:),       allocatable, save :: ipiel
  integer(i4),  dimension(:),       allocatable, save :: itmp
  integer(i4)                                         :: j
  integer(i4)                                         :: ja
  integer(i4)                                         :: k
  integer(i4)                                         :: lcwrk
  integer(i4)                                         :: n
  integer(i4)                                         :: nalpha
  integer(i4)                                         :: nati
  integer(i4)                                         :: nb
  integer(i4)                                         :: nbeta
  integer(i4)                                         :: nelpi
  integer(i4)                                         :: ni
  integer(i4)                                         :: nk
  integer(i4)                                         :: nkhloc
  integer(i4)                                         :: nn
  integer(i4)                                         :: npi
  integer(i4)                                         :: npicfg
  integer(i4)                                         :: pis
  integer(i4)                                         :: ierror
  logical                                             :: lgammaonly
  logical                                             :: lbiradical
  logical,                                       save :: lwarned = .false.
  complex(dpc)                                        :: Atrm
  complex(dpc)                                        :: phase
  complex(dpc), dimension(:,:),     allocatable, save :: Api
  real(dp),     dimension(:,:),     allocatable, save :: eps
  real(dp),     dimension(:,:),     allocatable, save :: occ
  real(dp),     dimension(:,:),     allocatable, save :: occa
  real(dp),     dimension(:,:),     allocatable, save :: occb
  real(dp),     dimension(:),       allocatable, save :: pisea
  real(dp),     dimension(:),       allocatable, save :: pisip
  real(dp),     dimension(:),       allocatable, save :: pispop
  real(dp),     dimension(:,:),     allocatable, save :: Pnew
  real(dp),     dimension(:,:),     allocatable, save :: Pold
  complex(dpc), dimension(:,:),     allocatable, save :: Ptmp
  complex(dpc), dimension(:,:),     allocatable, save :: Spi
  complex(dpc), dimension(:),       allocatable, save :: cwrk
  real(dp),     dimension(:),       allocatable, save :: w1
  real(dp),     dimension(:),       allocatable, save :: wkhloc
  real(dp)                                            :: bandocc
  real(dp)                                            :: dum
  real(dp)                                            :: dum2
  real(dp)                                            :: eel
  real(dp)                                            :: eafermi
  real(dp)                                            :: ebfermi
  real(dp)                                            :: eahomo
  real(dp)                                            :: ebhomo
  real(dp)                                            :: ealumo
  real(dp)                                            :: eblumo
  real(dp)                                            :: eold
  real(dp)                                            :: epsmax
  real(dp)                                            :: epsmin
  real(dp),                                      save :: kbinev = 8.617333262d-5
  real(dp)                                            :: kbt
  real(dp)                                            :: kr
  real(dp)                                            :: qpi
  real(dp)                                            :: qtot
  real(dp)                                            :: xd
  real(dp)                                            :: yd
  real(dp)                                            :: zd
  real(dp)                                            :: xk
  real(dp)                                            :: yk
  real(dp)                                            :: zk
  real(dp),     dimension(:),       allocatable, save :: xkv
  real(dp),     dimension(:),       allocatable, save :: ykv
  real(dp),     dimension(:),       allocatable, save :: zkv
!
!  Allocate local memory
!
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
  allocate(itmp(numat),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : itmp '')')
    stop
  endif
  allocate(wkhloc(nkh),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : wkhloc '')')
    stop
  endif
  allocate(xkv(nkh),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xkv '')')
    stop
  endif
  allocate(ykv(nkh),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ykv '')')
    stop
  endif
  allocate(zkv(nkh),stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : zkv '')')
    stop
  endif
!***************
!  Do Hueckel  *
!***************
  if (ipicount.gt.0) then
!
!  Generate Cartesian k vectors
!
    do nk = 1,nkh
      xk = xkh(nk)
      yk = ykh(nk)
      zk = zkh(nk)
      if (ndim.eq.3) then
        xkv(nk) = xk*kv(1,1) + yk*kv(1,2) + zk*kv(1,3)
        ykv(nk) = xk*kv(2,1) + yk*kv(2,2) + zk*kv(2,3)
        zkv(nk) = xk*kv(3,1) + yk*kv(3,2) + zk*kv(3,3)
      elseif (ndim.eq.2) then
        xkv(nk) = xk*kv(1,1) + yk*kv(1,2)
        ykv(nk) = xk*kv(2,1) + yk*kv(2,2)
        zkv(nk) = 0.0_dp
      elseif (ndim.eq.1) then
        xkv(nk) = xk*kv(1,1)
        ykv(nk) = 0.0_dp
        zkv(nk) = 0.0_dp
      endif
    enddo
!
    allocate(ipiel(numat),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : ipiel '')')
      stop
    endif
    allocate(pisea(ipicount),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : ipiseap '')')
      stop
    endif
    allocate(pisip(ipicount),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : pisip '')')
      stop
    endif
    allocate(pispop(ipicount),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : pispop '')')
      stop
    endif
    allocate(w1(3_i4*numat),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : w1 '')')
      stop
    endif
!
    itmp = 0 ! save pi atom info
    pisip = 0.0_dp
    pisea = 0.0_dp

    do pis = 1,ipicount ! loop over pi systems
      qpi    = 0.0_dp
      npi    = 0
      nelpi  = 0
      ipiadr3 = 0
      ipiadr4 = 0
      ipiel   = 0

      do k = 1,npiall
        if (ipimvec(k).eq.pis) then
          npi = npi + 1
          nati = nat(ipiadr(k))
          hybi = hybrid(ipiadr(k))
          qpi  = qpi - qf(ipiadr(k))
          ii   = nelpi
          if (nati.eq.5.and.hybi.eq.1)            nelpi = nelpi + 1  ! B in borine
          if (nati.eq.6.and.itag(ipiadr(k)).ne.1) nelpi = nelpi + 1  ! skip if its a carbene (tag itag=1)
          if (nati.eq.7.and.hybi.eq.2.and.itag(ipiadr(k)).eq.1) &
     &                                            nelpi = nelpi + 1  ! the itag=1 avoids an odd el number for the nitro group (its 4)
          if (nati.eq.7.and.hybi.le.2)            nelpi = nelpi + 1
          if (nati.eq.7.and.hybi.eq.3)            nelpi = nelpi + 2
          if (nati.eq.8.and.hybi.eq.1)            nelpi = nelpi + 1
          if (nati.eq.8.and.hybi.eq.2)            nelpi = nelpi + 1
          if (nati.eq.8.and.hybi.eq.3)            nelpi = nelpi + 2
          if (nati.eq.9.and.hybi.ne.1)            nelpi = nelpi + 2
          if (nati.eq.9.and.hybi.eq.1)            nelpi = nelpi + 3 !??? otherwise fluor-furan+ is wrong
          if (nati.eq.16.and.hybi.eq.1)           nelpi = nelpi + 1
          if (nati.eq.16.and.hybi.eq.2)           nelpi = nelpi + 1
          if (nati.eq.16.and.hybi.eq.3)           nelpi = nelpi + 2
          if (nati.eq.17.and.hybi.eq.0)           nelpi = nelpi + 2
          if (nati.eq.17.and.hybi.eq.1)           nelpi = nelpi + 3
          ipiadr3(npi) = ipiadr(k) ! map to original, full atom set
          ipiadr4(ipiadr(k)) = npi
          ipiel(ipiadr(k)) = nelpi - ii
          if (ipiel(ipiadr(k)).gt.2) ipiel(ipiadr(k)) = 2
        endif
      enddo
      nelpi = nelpi - ipis(pis)
      if (npi.lt.2.or.nelpi.lt.1) cycle
!
!  Set k points if PBC
!
      if (lpimpbc(pis)) then
        lgammaonly = .false.
        nkhloc = nkh
        wkhloc(1:nkh) = wkh(1:nkh)
      else
        lgammaonly = .true.
        nkhloc = 1_i4
        wkhloc(1) = 1.0_dp
      endif
!
!  Allocate workspace for Hueckel solve
!
      allocate(Api(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Api '')')
        stop
      endif
      allocate(Pnew(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Pnew '')')
        stop
      endif
      allocate(Pold(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Pold '')')
        stop
      endif
      allocate(Ptmp(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Ptmp '')')
        stop
      endif
      allocate(Spi(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Spi '')')
        stop
      endif
      allocate(eps(npi,nkh),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : eps '')')
        stop
      endif
      allocate(occ(npi,nkh),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : occ '')')
        stop
      endif
      allocate(occa(npi,nkh),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : occa '')')
        stop
      endif
      allocate(occb(npi,nkh),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : occb '')')
        stop
      endif
!
      nb = ilaenv( 1_i4, 'ZHETRD', 'U', npi, -1_i4, -1_i4, -1_i4 )
      lcwrk = (nb + 1)*npi
      allocate(cwrk(lcwrk),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : cwrk '')')
        stop
      endif
!
      eold = 0.0_dp
      Pold = 2.0_dp/3.0_dp
!****************************************************************************
!  Start of process that can be repeated if number of pi electrons changes  *
!****************************************************************************
      picfgloop: do npicfg = 1,2
!**************************************************************************************************
!  Iterative Hueckel loop, off-diag terms are reduced depending on P to avoid overdelocalization  *
!**************************************************************************************************
        iterloop: do nn = 1,maxhiter      ! just some iterations
!--------------------------
!  Find orbital energies  |
!--------------------------
!
!  Loop over k points
!
          do nk = 1,nkhloc
            if (lgammaonly) then
              xk = 0.0_dp
              yk = 0.0_dp
              zk = 0.0_dp
            else
              xk = xkv(nk)
              yk = ykv(nk)
              zk = zkv(nk)
            endif
!
!  Initialise matrix
!
            Api(1:npi,1:npi) = 0.0_dpc
!
!  Diagonal elements
!
            do i = 1,npi
              ii = ipiadr3(i)
              Api(i,i) = dcmplx(hdiag(nat(ii)) + qf(ii)*hueckelp3 - dble(ipiel(ii)-1)*pilpf,0.0_dp)
            enddo
!
!  Loop over bonds for pair interactions
!
            do i = 1,numat
              do ni = 1,nnbr_bond(i)
                j = nbrno_bond(ni,i)
                ia = ipiadr4(i)
                ja = ipiadr4(j)
                if (ia.gt.0.and.ja.gt.0) then
                  dum = 1.d-9*rbnbr(ni,i)                               ! distort so that Huckel for e.g. COT localizes to right bonds
                  dum = sqrt(hoffdiag(nat(i))*hoffdiag(nat(j))) - dum   ! better than arithmetic
                  dum2 = hiter
                  if (hybrid(i).eq.1) dum2 = dum2*htriple               ! triple bond is different
                  if (hybrid(j).eq.1) dum2 = dum2*htriple               ! triple bond is different
                  Atrm = dcmplx(dum*(1.0_dp-dum2*(2.0_dp/3.0_dp-Pold(ja,ia))),0.0_dp) ! Pmat scaling with benzene as reference
!
!  Compute phase factor
!
                  kr = xk*xbnbr(ni,i) + yk*ybnbr(ni,i) + zk*zbnbr(ni,i)
                  phase = dcmplx(cos(kr),-sin(kr))
!
!  Add phased term
!
                  Api(ja,ia) = Api(ja,ia) - Atrm*phase
                endif
              enddo
            enddo
!
!  Diagonalise matrix
!
            call zheev('N','U',npi,Api,npi,eps(1,nk),cwrk,lcwrk,w1,ifail)
            if (ifail.ne.0) then
              call pgfnff_error('diagonalisation of Huckel matrix has failed',0_i4)
              stop
            endif
!
!  Scale eigenvalues so that gaps are roughly ok (in eV)
!
            eps(1:npi,nk) = eps(1:npi,nk)*0.1_dp*27.2113957_dp
!
!  End of loop over k points
!
          enddo

!--------------------------
!  Solve for Fermi level  |
!--------------------------
!
!  Original implementation used a different temperature for the second configuration
!
          if (npicfg.eq.1) then
            kbt = gfnff_pi_temp1*kbinev
          else
            kbt = gfnff_pi_temp2*kbinev
          endif
!
!  Separate alpha and beta problem 
!
          nalpha = nelpi/2
          nbeta = nalpha
          if (2*nalpha.ne.nelpi) then
            nalpha = nalpha + 1
          endif
!
!  Solve for alpha occupancies
!
          qtot = dble(nalpha)
          call fermi(nkhloc, wkhloc, npi, npi, eps, kbt, qtot, occa, eafermi, eahomo, ealumo)
!
!  Solve for beta occupancies
!
          qtot = dble(nbeta)
          call fermi(nkhloc, wkhloc, npi, npi, eps, kbt, qtot, occb, ebfermi, ebhomo, eblumo)
!
!  Combine occupancies
!
          occ = occa + occb
!
!  Check for biradical for case with a single k point
!
          if (nkhloc.eq.1) then
            do i = npi,1,-1
              if (occ(i,1).gt.1.0d-6) then
                if (i.gt.1) then
                  lbiradical = (abs(occ(i,1)-occ(i-1,1)).lt.1.0d-4)
                  if (lbiradical) then
                    occ(1:npi,1) = 0.0_dp
                    do j = 1,nelpi/2
                      occ(j,1) = 2.0_dp
                    enddo
                    if (.not.lwarned) then
                      call pgfnff_warn('Breaking symmetry due to perfect biradical',0_i4)
                      lwarned = .true.
                    endif
                  endif
                endif
                exit
              endif         
            enddo
          endif
!------------------------------
!  Compute electronic energy  !
!------------------------------
          eel = 0.0_dp
          do nk = 1,nkhloc
            do j = 1,npi
              eel = eel + eps(j,nk)*occ(j,nk)
            enddo
          enddo
!-----------------------------
!  Construct density matrix  |
!-----------------------------
          Pnew(1:npi,1:npi) = 0.0_dp
!
!  Loop over k points
!
          do nk = 1,nkhloc
            if (lgammaonly) then
              xk = 0.0_dp
              yk = 0.0_dp
              zk = 0.0_dp
            else
              xk = xkv(nk)
              yk = ykv(nk)
              zk = zkv(nk)
            endif
!
!  Initialise matrix
!
            Spi(1:npi,1:npi) = 0.0_dpc
!
!  Diagonal elements
!
            do i = 1,npi
              ii = ipiadr3(i)
              Spi(i,i) = dcmplx(hdiag(nat(ii)) + qf(ii)*hueckelp3 - dble(ipiel(ii)-1)*pilpf,0.0_dp)
            enddo
!
!  Loop over bonds for pair interactions
!
            do i = 1,numat
              do ni = 1,nnbr_bond(i)
                j = nbrno_bond(ni,i)
                ia = ipiadr4(i)
                ja = ipiadr4(j)
                if (ia.gt.0.and.ja.gt.0) then
                  dum = 1.d-9*rbnbr(ni,i)                               ! distort so that Huckel for e.g. COT localizes to right bonds
                  dum = sqrt(hoffdiag(nat(i))*hoffdiag(nat(j))) - dum   ! better than arithmetic
                  dum2 = hiter
                  if (hybrid(i).eq.1) dum2 = dum2*htriple               ! triple bond is different
                  if (hybrid(j).eq.1) dum2 = dum2*htriple               ! triple bond is different
                  Atrm = dcmplx(dum*(1.0_dp-dum2*(2.0_dp/3.0_dp-Pold(ja,ia))),0.0_dp) ! Pmat scaling with benzene as reference
!
!  Compute phase factor
!
                  kr = xk*xbnbr(ni,i) + yk*ybnbr(ni,i) + zk*zbnbr(ni,i)
                  phase = dcmplx(cos(kr),-sin(kr))
!
!  Add phased term
!
                  Spi(ja,ia) = Spi(ja,ia) - Atrm*phase
                endif
              enddo
            enddo
!
!  Diagonalise matrix
!
            call zheev('V','U',npi,Spi,npi,eps(1,nk),cwrk,lcwrk,w1,ifail)
            if (ifail.ne.0) then
              call pgfnff_error('diagonalisation of Huckel matrix has failed',0_i4)
              stop
            endif
!
!  Scale eigenvalues so that gaps are roughly ok (in eV)
!
            eps(1:npi,nk) = eps(1:npi,nk)*0.1_dp*27.2113957_dp
!
!  Build unphased density matrix
!
            do i = 1,npi
              do j = 1,npi
                Api(j,i) = 0.0_dpc
                do n = 1,npi
                  Api(j,i) = Api(j,i) + dconjg(Spi(j,n))*occ(n,nk)*Spi(i,n)
                enddo
              enddo
            enddo
!
            do i = 1,numat
              ia = ipiadr4(i)
              if (ia.gt.0) then
                do j = 1,numat
                  ja = ipiadr4(j)
                  if (ja.gt.0) then
                    xd = xclat(j) - xclat(i)
                    yd = yclat(j) - yclat(i)
                    zd = zclat(j) - zclat(i)
!
!  Compute phase factor
!
                    kr = xk*xd + yk*yd + zk*zd
                    phase = dcmplx(cos(kr),-sin(kr))
!
                    Pnew(ja,ia) = Pnew(ja,ia) + dble(Api(ja,ia)*phase)
                  endif
                enddo
              endif
            enddo
!
!  End of loop over k points
!
          enddo
!
!  Set overall homo/lumo
!
          do i = 1,npi  ! save IP/EA
            bandocc = 0.0_dp
            epsmax = eps(i,1)
            epsmin = eps(i,1)
            if (i+1.lt.npi) epsmin = eps(i+1,1)
            do nk = 1,nkhloc
              bandocc = bandocc + occ(i,nk)
              epsmax = max(epsmax,eps(i,nk))
              if (i+1.lt.npi) then
                epsmin = min(epsmin,eps(i+1,nk))
              endif
            enddo
            if (bandocc.gt.0.5_dp) then
              pisip(pis) = epsmax   ! IP
              if (i+1.lt.npi) pisea(pis) = epsmin ! EA
            endif
          enddo
!
!  Copy Pnew to Pold
!
          Pold = Pnew
!
          if (abs(eel-eold).lt.1.d-4) exit  ! end of iterations
          eold = eel
        enddo iterloop
!**************************
!  End of iterative loop  *
!**************************
        if (npicfg.eq.1) then
          if (lverbose) then
            if (lpimpbc(pis)) then
              write(ioout,'(''Hueckel system :'',i3,'' charge : '',i3,'' ndim/Nel :'',2i5,3x,''eps(HOMO/LUMO)'',2f12.6,'' PBC'')') &
                pis,ipis(pis),npi,nelpi,pisip(pis),pisea(pis)
            else
              write(ioout,'(''Hueckel system :'',i3,'' charge : '',i3,'' ndim/Nel :'',2i5,3x,''eps(HOMO/LUMO)'',2f12.6)') &
                pis,ipis(pis),npi,nelpi,pisip(pis),pisea(pis)
            endif
          endif
          if (pisip(pis).gt.0.40_dp) then
            if (lverbose) then
              if (gfnff_pi_change.eq.0) then
                if (qpi.lt.0.0_dp) then
                  call pgfnff_warn('Probably wrong pi occupation in Huckel. 2nd attempt with Nel=Nel-1!',0_i4)
                else
                  call pgfnff_warn('Probably wrong pi occupation in Huckel. 2nd attempt with Nel=Nel+1!',0_i4)
                endif
              else
                if (gfnff_pi_change.lt.0) then
                  call pgfnff_warn('Probably wrong pi occupation in Huckel. 2nd attempt with Nel=Nel-1!',0_i4)
                else
                  call pgfnff_warn('Probably wrong pi occupation in Huckel. 2nd attempt with Nel=Nel+1!',0_i4)
                endif
              endif
              do i = 1,numat
                if (ipiadr4(i).ne.0) write(ioout,*) 'at,nb,hyb,Npiel:', i,nnbr_bond(i),hybrid(i),ipiel(i)
              enddo
            endif
            if (gfnff_pi_change.eq.0) then
              if (qpi.lt.0.0_dp) then
!
!  Reduce number of pi electrons by 1
!
                nelpi = nelpi - 1
              else
!
!  Increment number of pi electrons by 1
!
                nelpi = nelpi + 1
              endif
            else
              if (gfnff_pi_change.lt.0) then
!
!  Reduce number of pi electrons by 1
!
                nelpi = nelpi - 1
              else
!
!  Increment number of pi electrons by 1
!
                nelpi = nelpi + 1
              endif
            endif
!
!  Repeat procedure of building pibond with new number of electrons
!
            cycle picfgloop
          endif
        endif
!
!  Save BO or set to zero if not part of the pi system
!
        do i = 1,numat
          do ni = 1,nnbr_bond(i)
            j = nbrno_bond(ni,i)
            ia = ipiadr4(i)
            ja = ipiadr4(j)
            if (ia.gt.0.and.ja.gt.0) then
              pibond(ni,i) = Pold(ja,ia)
              itmp(i) = 1
              itmp(j) = 1
            endif
          enddo
        enddo
        exit
      enddo picfgloop
!
!  Free Hueckel workspace
!
      deallocate(cwrk,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : cwrk '')')
        stop
      endif
      deallocate(occb,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : occb '')')
        stop
      endif
      deallocate(occa,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : occa '')')
        stop
      endif
      deallocate(occ,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : occ '')')
        stop
      endif
      deallocate(eps,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : eps '')')
        stop
      endif
      deallocate(Ptmp,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Ptmp '')')
        stop
      endif
      deallocate(Pold,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Pold '')')
        stop
      endif
      deallocate(Pnew,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Pnew '')')
        stop
      endif
      deallocate(Spi,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Spi '')')
        stop
      endif
      deallocate(Api,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Api '')')
        stop
      endif
    enddo
!**************************
!  End of pi system loop  *
!**************************
    ipiadr = itmp  ! array used for identifying pi atoms in following codes
!
!  Free memory that is only for this section
!
    deallocate(w1,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : w1 '')')
      stop
    endif
    deallocate(pispop,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : pispop '')')
      stop
    endif
    deallocate(pisip,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : pisip '')')
      stop
    endif
    deallocate(pisea,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : pisea '')')
      stop
    endif
    deallocate(ipiel,stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Deallocation : ipiel '')')
      stop
    endif
  endif
!
!  Deallocate local memory
!
  deallocate(zkv,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : zkv '')')
    stop
  endif
  deallocate(ykv,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : ykv '')')
    stop
  endif
  deallocate(xkv,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : xkv '')')
    stop
  endif
  deallocate(wkhloc,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : wkhloc '')')
    stop
  endif
  deallocate(itmp,stat=ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Deallocation : itmp '')')
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
!
  end subroutine pgfnff_huckel_pbc
