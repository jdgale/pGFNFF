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
  subroutine pgfnff_outpot(numat,maxnbr,nnbr_bond,nbrno_bond)
!
!  Outputs GFNFF interatomic potential information 
!  NB: Should only be called if GFNFF is to be used and this is the I/O processor
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use m_io
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp
  use m_pgfnff_nbr_lib
  use m_pgfnff_topo
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)     :: numat                    ! Number of atoms
  integer(i4),                  intent(in)     :: maxnbr                   ! Maximum number of neighbours
  integer(i4),                  intent(in)     :: nnbr_bond(numat)         ! Number of neighbours
  integer(i4),                  intent(in)     :: nbrno_bond(maxnbr,numat) ! Pointer to neighbours
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: j
  integer(i4)                                  :: ni
  real(dp)                                     :: r0
  real(dp)                                     :: v2
  real(dp)                                     :: v3
!
!***************************************
!  Output pGFNFF potential parameters  *
!***************************************
  write(ioout,'(/,''  pGFNFF potential parameters: '',/)')
!
!  Bond energy
!
  write(ioout,'(''  Bond energy: '',/)')
  write(ioout,'(''----------------------------------------------------------------------'')')
  write(ioout,'(''  I     J            r0 (initial)      Exponent        Pre-exponential'')')
  write(ioout,'(''                        (Ang)           (Ang^-2)            (eV) '')')
  write(ioout,'(''----------------------------------------------------------------------'')')
  do i = 1,numat
    do ni = 1,nnbr_bond(i)
      j = nbrno_bond(ni,i)
      r0 = vbnbr(1,ni,i)
      v2 = vbnbr(2,ni,i)
      v3 = vbnbr(3,ni,i)
      write(ioout,'(i4,1x,i4,5x,f8.5,1x,f12.6,1x,f12.6)')  i,j,r0,v2,v3
    enddo
  enddo
!
  return
  end
