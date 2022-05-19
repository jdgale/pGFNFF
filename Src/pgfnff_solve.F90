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
  subroutine pgfnff_solve(A,S,et,nbasis,nel,eel,focc,e)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)    :: nbasis              ! Number basis of functions
  integer(i4), intent(in)    :: nel                 ! Number of electrons
  real(dp),    intent(in)    :: et                  ! Electronic temp Fermi smear
  real(dp),    intent(out)   :: eel                 ! Electronic energy = sum_occ nocc*eps
  real(dp),    intent(out)   :: focc(nbasis)        ! Occupations
  real(dp),    intent(out)   :: e(nbasis)           ! Eigenvalues
  real(dp),    intent(inout) :: A(nbasis,nbasis)
  real(dp),    intent(out)   :: S(nbasis,nbasis)
!
!  Local variables
!
  integer(i4)                :: ihomoa
  integer(i4)                :: ihomob
  integer(i4)                :: i
  integer(i4)                :: info
  integer(i4)                :: lwork
  integer(i4)                :: m
  logical,              save :: lwarned = .false.
  real(dp)                   :: efa
  real(dp)                   :: efb
  real(dp)                   :: foda
  real(dp)                   :: fodb
  real(dp),      allocatable :: focca(:)
  real(dp),      allocatable :: foccb(:)
  real(dp),      allocatable :: aux(:)
  real(dp),      allocatable :: Ptmp(:,:)
!
  allocate(focca(nbasis),foccb(nbasis))
!
!  ZDO case
!
  lwork  = 1 + 6*nbasis + 2*nbasis**2
  allocate(aux(lwork))
  call dsyev('V','U',nbasis,A,nbasis,e,aux,lwork,info)
!
!  Scale energy so that gaps are roughly ok (in eV)
!
  e = e*0.1_dp*27.2113957_dp

  if (info.ne.0) then
    call pgfnff_error('Error in diagonalisation in pgfnff_solve',0_i4)
    stop
  endif
!
  if (et.gt.1.d-3) then
!
!  Fermi smearing, convert restricted occ first to alpha/beta
!
    call occu(nbasis,nel,ihomoa,ihomob,focc,focca,foccb)
    if (ihomoa.le.nbasis.and.ihomoa.gt.0) then
      call fermismear(nbasis,ihomoa,et,e,focca,foda,efa)
    else
      focca = 0
    endif
    if (ihomob.le.nbasis.and.ihomob.gt.0) then
      call fermismear(nbasis,ihomob,et,e,foccb,fodb,efb)
    else
      foccb = 0
    endif
    focc = focca + foccb
    if (ihomoa+1.le.nbasis) then
      if (abs(focc(ihomoa)-focc(ihomoa+1)).lt.1.d-4) then ! a perfect birad is anti-aromatic
        focc = 0                                          ! and hence we break the sym
        do i = 1,nel/2
          focc(i) = 2.0_dp
        enddo
        if (.not.lwarned) then
          call pgfnff_warn('Breaking symmetry due to perfect biradical',0_i4)
          lwarned = .true.
        endif
      endif
    endif
  else
    focc = 0
    do i = 1,nel/2
      focc(i) = 2.0_dp
    enddo
    if (2*(nel/2).ne.nel) focc(nel/2+1) = 1.0_dp
  endif

  focca = focc*e
  eel = sum(focca)
!
!  Density matrix
!
  S = A
  allocate(Ptmp(nbasis,nbasis))
  do m = 1,nbasis
    do i = 1,nbasis
      Ptmp(i,m) = S(i,m)*focc(m)
    enddo
  enddo
  call dgemm('n','t',nbasis,nbasis,nbasis,1.0_dp,S,nbasis,Ptmp,nbasis,0.0_dp,A,nbasis)

  deallocate(Ptmp)
  deallocate(focca,foccb)
      
  end subroutine pgfnff_solve
!
  subroutine fermismear(norbs,nel,t,eig,occ,fod,e_fermi)
  use m_pgfnff_types
!
!  Passed variables
!
  integer,      intent(in)  :: norbs
  integer,      intent(in)  :: nel
  real(dp),     intent(in)  :: eig(norbs)
  real(dp),     intent(out) :: occ(norbs)
  real(dp),     intent(in)  :: t
  real(dp),     intent(out) :: fod
  real(dp),     intent(out) :: e_fermi
!
!  Local variables
!
  integer(i4)               :: i
  integer(i4)               :: ncycle
  real(dp),       parameter :: kB = 3.166808578545117e-06_dp
  real(dp),       parameter :: autoev = 27.21138505_dp
  real(dp)                  :: boltz
  real(dp)                  :: bkt
  real(dp)                  :: change_fermi
  real(dp)                  :: dfermifunct
  real(dp)                  :: fermifunct
  real(dp)                  :: occt
  real(dp)                  :: thr
  real(dp)                  :: total_dfermi
  real(dp)                  :: total_number
!
  parameter (boltz = kB*autoev)
  parameter (thr   = 1.d-9)
!
  bkt = boltz*t
!
  if (nel+1.gt.norbs) then
    e_fermi = eig(nel)
    occ(1:norbs) = 1.0_dp
  else
    e_fermi = 0.5_dp*(eig(nel)+eig(nel+1))
    occt = nel
    do ncycle = 1,200  ! this loop would be possible instead of gotos
      total_number = 0.0_dp
      total_dfermi = 0.0_dp
      do i = 1,norbs
        fermifunct = 0.0_dp
        if ((eig(i)-e_fermi)/bkt.lt.50.0_dp) then
          fermifunct = 1.0_dp/(exp((eig(i)-e_fermi)/bkt)+1.0_dp)
          dfermifunct = exp((eig(i)-e_fermi)/bkt) / &
              &       (bkt*(exp((eig(i)-e_fermi)/bkt)+1.0_dp)**2)
        else
          dfermifunct = 0.0_dp
        endif
        occ(i) = fermifunct
        total_number = total_number + fermifunct
        total_dfermi = total_dfermi + dfermifunct
      enddo
      change_fermi = (occt-total_number)/total_dfermi
      e_fermi = e_fermi+change_fermi
      if (abs(occt-total_number).le.thr) exit
    enddo
  endif

  fod = 0.0_dp
  do i = 1,norbs
    if (eig(i).lt.e_fermi) then
      fod = fod + 1.0_dp - occ(i)
    else
      fod = fod + occ(i)
    endif
  enddo

  end subroutine fermismear
!
  subroutine occu(nbasis,nel,ihomoa,ihomob,focc,focca,foccb)
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4),    intent(in)  :: nel
  integer(i4),    intent(in)  :: nbasis
  integer(i4),    intent(out) :: ihomoa
  integer(i4),    intent(out) :: ihomob
  real(dp),       intent(out) :: focc(nbasis)
  real(dp),       intent(out) :: focca(nbasis)
  real(dp),       intent(out) :: foccb(nbasis)
!
!  Local variables
!
  integer(i4)                 :: i
  integer(i4)                 :: ihomo
  integer(i4)                 :: na
  integer(i4)                 :: nb
!
  focc  = 0.0_dp
  focca = 0.0_dp
  foccb = 0.0_dp
!
!  Even Nel
!
  if (mod(nel,2).eq.0) then
    ihomo = nel/2
    do i = 1,ihomo
      focc(i) = 2.0_dp
    enddo
    if (2*ihomo.ne.nel) then
      ihomo = ihomo + 1
      focc(ihomo) = 1.0_dp
    endif
!
!  Odd nel
!
  else
    na = nel/2 + 1
    nb = nel/2 
    do i = 1,na
      focc(i) = focc(i) + 1.0_dp
    enddo
    do i = 1,nb
      focc(i) = focc(i) + 1.0_dp
    enddo
  endif
!
  do i = 1,nbasis
    if (focc(i).eq.2) then
      focca(i) = 1.0_dp
      foccb(i) = 1.0_dp
    endif
    if (focc(i).eq.1) focca(i) = 1.0_dp
  enddo
!
  ihomoa = 0
  ihomob = 0
  do i = 1,nbasis
    if (focca(i).gt.0.99_dp) ihomoa = i
    if (foccb(i).gt.0.99_dp) ihomob = i
  enddo

  end subroutine occu
