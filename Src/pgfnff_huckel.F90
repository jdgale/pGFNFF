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
  subroutine pgfnff_huckel(numat,nat,qf,ipicount,npiall,ipimvec,ipiadr,ipis,itag,hybrid,pibond, &
                           maxnbr,nnbr_bond,nbrno_bond,rbnbr,lverbose)
!
!  Solve the Hueckel hamilton for a non-periodic system
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  implicit none
!
!  Passed variables
!
  integer(i4),                         intent(in)    :: numat                     ! Number of atoms
  integer(i4),                         intent(in)    :: nat(numat)                ! Atomic numbers
  integer(i4),                         intent(in)    :: hybrid(numat)             ! Hybridisation state
  integer(i4),                         intent(in)    :: ipicount
  integer(i4),                         intent(inout) :: ipiadr(numat)
  integer(i4),                         intent(in)    :: ipimvec(numat)
  integer(i4),                         intent(in)    :: ipis(numat)
  integer(i4),                         intent(in)    :: itag(numat)
  integer(i4),                         intent(in)    :: maxnbr                    ! Maximum number of neighbours
  integer(i4),                         intent(in)    :: nnbr_bond(numat)          ! Number of neighbours
  integer(i4),                         intent(in)    :: nbrno_bond(maxnbr,numat)  ! Pointers to neighbours
  integer(i4),                         intent(in)    :: npiall
  logical,                             intent(in)    :: lverbose                  ! If true then provide verbose output
  real(dp),                            intent(out)   :: pibond(maxnbr,numat)      ! pi bond orders
  real(dp),                            intent(in)    :: qf(numat)                 ! Charges
  real(dp),                            intent(in)    :: rbnbr(maxnbr,numat)       ! Distances for neighbours
!
!  Local variables
!
  integer(i4)                                        :: hybi
  integer(i4)                                        :: i
  integer(i4)                                        :: ia
  integer(i4)                                        :: ii
  integer(i4), dimension(:),       allocatable, save :: ipiadr3
  integer(i4), dimension(:),       allocatable, save :: ipiadr4
  integer(i4), dimension(:),       allocatable, save :: ipiel
  integer(i4), dimension(:),       allocatable, save :: itmp
  integer(i4)                                        :: j
  integer(i4)                                        :: ja
  integer(i4)                                        :: k
  integer(i4)                                        :: nati
  integer(i4)                                        :: nelpi
  integer(i4)                                        :: ni
  integer(i4)                                        :: nn
  integer(i4)                                        :: npi
  integer(i4)                                        :: pis
  integer(i4)                                        :: ierror
  real(dp),    dimension(:,:),     allocatable, save :: Api
  real(dp),    dimension(:,:),     allocatable, save :: Apisave
  real(dp),    dimension(:),       allocatable, save :: eps
  real(dp),    dimension(:),       allocatable, save :: occ
  real(dp),    dimension(:),       allocatable, save :: pisea
  real(dp),    dimension(:),       allocatable, save :: pisip
  real(dp),    dimension(:),       allocatable, save :: pispop
  real(dp),    dimension(:,:),     allocatable, save :: Pold
  real(dp),    dimension(:,:),     allocatable, save :: Spi
  real(dp)                                           :: dum
  real(dp)                                           :: dum2
  real(dp)                                           :: eold
  real(dp)                                           :: qpi
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
!***************
!  Do Hueckel  *
!***************
  if (ipicount.gt.0) then
    allocate(ipiel(numat),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : ipiel '')')
      stop
    endif
    allocate(pisea(ipicount),stat=ierror)
    if (ierror.gt.0) then
      write(ioout,'('' Error in Memory Allocation : piseap '')')
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
!  Allocate workspace for Hueckel solve
!
      allocate(Api(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Api '')')
        stop
      endif
      allocate(Apisave(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Apisave '')')
        stop
      endif
      allocate(Pold(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Pold '')')
        stop
      endif
      allocate(Spi(npi,npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : Spi '')')
        stop
      endif
      allocate(eps(npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : eps '')')
        stop
      endif
      allocate(occ(npi),stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Allocation : occ '')')
        stop
      endif
!
      eold = 0.0_dp
      Pold = 2.0_dp/3.0_dp
!
!  Iterative Hueckel loop, off-diag terms are reduced depending on P to avoid overdelocalization
!
      do nn = 1,maxhiter      ! just some iterations
        Api(1:npi,1:npi) = 0.0_dp
        do i = 1,npi
          ii = ipiadr3(i)
          Api(i,i) = hdiag(nat(ii)) + qf(ii)*hueckelp3 - dble(ipiel(ii)-1)*pilpf
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
              Api(ja,ia) = Api(ja,ia) - dum*(1.0_dp-dum2*(2.0_dp/3.0_dp-Pold(ja,ia))) ! Pmat scaling with benzene as reference
            endif
          enddo
        enddo
!
!  Save Api
!
        Apisave(1:npi,1:npi) = Api(1:npi,1:npi)
!
!  Solve
!
        call pgfnff_solve(Api,Spi,4000.0_dp,npi,nelpi,dum,occ,eps)  !diagonalize, 4000 better than 300
!
        do i = 1,npi  ! save IP/EA
          if (occ(i).gt.0.5_dp) then
            pisip(pis) = eps(i)   ! IP
            if (i+1.lt.npi) pisea(pis) = eps(i+1) ! EA
          endif
        enddo
        if (abs(dum-eold).lt.1.d-4) exit  ! end of iterations
        Pold(1:npi,1:npi) = Api(1:npi,1:npi)
        eold = dum
      enddo
!
!  End of iterative loop
!
      if (lverbose) then
        write(ioout,'(''Hueckel system :'',i3,'' charge : '',i3,'' ndim/Nel :'',2i5,3x,''eps(HOMO/LUMO)'',2f12.6)') &
          pis,ipis(pis),npi,nelpi,pisip(pis),pisea(pis)
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
        Api(1:npi,1:npi) = Apisave(1:npi,1:npi)
!
!  Solve
!
        call pgfnff_solve(Api,Spi,300.0_dp,npi,nelpi,dum,occ,eps)  !diagonalize
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
            pibond(ni,i) = Api(ja,ia)
            itmp(i) = 1
            itmp(j) = 1
          endif
        enddo
      enddo
!
!  Free Hueckel workspace
!
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
      deallocate(Pold,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Pold '')')
        stop
      endif
      deallocate(Spi,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Spi '')')
        stop
      endif
      deallocate(Apisave,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Apisave '')')
        stop
      endif
      deallocate(Api,stat=ierror)
      if (ierror.gt.0) then
        write(ioout,'('' Error in Memory Deallocation : Api '')')
        stop
      endif
    enddo
!
!  End of pi system loop
!
    ipiadr = itmp  ! array used for identifying pi atoms in following codes
!
!  Free memory that is only for this section
!
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
  end subroutine pgfnff_huckel
