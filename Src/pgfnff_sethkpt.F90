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
  subroutine pgfnff_sethkpt(ndim,kv,lverbose)
!
!  Sets up k points according to shrinking factors for Huckel solve
!
!  nkh = total number of k points 
!  xkh = fractional x component of k point
!  ykh = fractional y component of k point
!  zkh = fractional z component of k point
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use m_pgfnff_types
  use m_io
  use m_pgfnff,        only : gfnff_ks
  use m_pgfnff_cfg,    only : maxnkh, changemaxnkh
  use m_pgfnff_cfg,    only : nkh, xkh, ykh, zkh, wkh
  implicit none
!
!  Passed variables
!
  integer(i4),   intent(in)   :: ndim      ! Dimensionality of system
  logical,       intent(in)   :: lverbose  ! If true then verbose output will be provided
  real(dp),      intent(in)   :: kv(3,3)   ! Reciprocal lattice vectors
!
!  Local variables
!
  integer(i4)                 :: i
  integer(i4)                 :: j
  integer(i4)                 :: k
  integer(i4)                 :: nx
  integer(i4)                 :: nx2
  integer(i4)                 :: ny
  integer(i4)                 :: nz
  real(dp)                    :: delx
  real(dp)                    :: dely
  real(dp)                    :: delz
  real(dp)                    :: ka
  real(dp)                    :: kb
  real(dp)                    :: kc
  real(dp)                    :: pi
  real(dp)                    :: rks
  real(dp)                    :: rx
  real(dp)                    :: ry
  real(dp)                    :: rz
  real(dp)                    :: sumwkh
  real(dp)                    :: wx
  real(dp)                    :: xadd
  real(dp)                    :: yadd
  real(dp)                    :: zadd
!
!  Initial call to ensure arrays are dimensioned as per initial value of maxnkh
!
  call changemaxnkh
!
  pi = 4.0_dp*atan(1.0_dp)
!
!  Compute reciprocal lattice vector lengths
!
  ka = sqrt(kv(1,1)**2 + kv(2,1)**2 + kv(3,1)**2)
  kb = sqrt(kv(1,2)**2 + kv(2,2)**2 + kv(3,2)**2)
  kc = sqrt(kv(1,3)**2 + kv(2,3)**2 + kv(3,3)**2)
!
!  Multiply gfnff_ks by 2*pi to match the reciprocal lattice vectors
!
  rks = 0.5_dp/(gfnff_ks*pi)
!
  nx = ka*rks + 1
  ny = kb*rks + 1
  nz = kc*rks + 1
!
  if (ndim.eq.2) then
    nz = 1
  elseif (ndim.eq.1) then
    ny = 1
    nz = 1
  endif
!
  delx = 1.0_dp/float(nx)
  dely = 1.0_dp/float(ny)
  delz = 1.0_dp/float(nz)
!
  xadd = 0.5_dp*delx
  yadd = 0.5_dp*dely
  zadd = 0.5_dp*delz
!
  if (ndim.le.2) then
    delz = 0.0_dp
    zadd = 0.0_dp
  endif
  if (ndim.eq.1) then
    dely = 0.0_dp
    yadd = 0.0_dp
  endif
!
!  Apply time reversal symmetry in x
!
  nx2 = nx/2
  if (2*nx2.ne.nx) nx2 = nx2 + 1
!*****************************************************
!  Generate k points across the full Brillouin zone  *
!*****************************************************
!
!  Check there is enough space to insert new K points
!
  if (nx*ny*nz.gt.maxnkh) then
    maxnkh = nx*ny*nz
    call changemaxnkh
  endif
  nkh = 0
  rx = - delx + xadd
  sumwkh = 0.0_dp
  do i = 1,nx2
    rx = rx + delx
    if (abs(rx).lt.1.0d-8.or.abs(rx-0.5_dp).lt.1.0d-8) then
      wx = 0.5_dp
    else
      wx = 1.0_dp
    endif
    ry = - dely + yadd
    do j = 1,ny
      ry = ry + dely
      rz = - delz + zadd
      do k = 1,nz
        rz = rz + delz
        nkh = nkh + 1
!
!  Add in new point
!
        xkh(nkh) = rx
        ykh(nkh) = ry
        zkh(nkh) = rz
        wkh(nkh) = wx
        sumwkh = sumwkh + wx
      enddo
    enddo
  enddo
!
!  Normalise weights
!
  do i = 1,nkh
    wkh(i) = wkh(i)/sumwkh
  enddo
!
  if (lverbose) then
    write(ioout,'(/,''  GFNFF: Huckel K points : '',/)')
    do i = 1,nkh
      write(ioout,'(i6,1x,3f8.5,1x,f6.4)') i,xkh(i),ykh(i),zkh(i),wkh(i)
    enddo
    write(ioout,'(/)')
  endif
!
  return
  end
