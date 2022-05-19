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
module m_nbr
! 
!  This module contains the variables that hold the configuration specific
!  neighbour list information for a system 
!
!  Julian Gale, Curtin University, July 2021
!
  use m_pgfnff_types
!
  implicit none
!
  integer(i4),                              save :: maxnbr = 6          ! Maximum number of neighbours
  real(dp),                                 save :: cutnbr = 0.0_dp     ! Cutoff used to generate list
!
!  General neighbour list
!
  integer(i4), dimension(:),       pointer, save :: nnbr => null()      ! Number of neighbours
  integer(i4), dimension(:,:),     pointer, save :: nbrno => null()     ! Pointer to atom that is neighbour
  integer(i4), dimension(:,:),     pointer, save :: ncnbr => null()     ! Pointer to cell of neighbour 
  real(dp),    dimension(:,:),     pointer, save :: rnbr => null()      ! Distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: xnbr => null()      ! x component of distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: ynbr => null()      ! y component of distance to neighbour
  real(dp),    dimension(:,:),     pointer, save :: znbr => null()      ! z component of distance to neighbour
!
!  Coordination number pointers to data within nnbr / nbrno and derivatives
!
  integer(i4), dimension(:),       pointer, save :: nnbr_cn   => null()      ! Number of neighbours with non-zero coordination number
  integer(i4), dimension(:,:),     pointer, save :: nbrno_cn  => null()      ! Pointer to neighbour that has non-zero coordination number
!
!  Bonding lists
!
  integer(i4), dimension(:),       pointer, save :: nnbr_bond => null()      ! Number of bonded neighbours
  integer(i4), dimension(:,:),     pointer, save :: nbrno_bond => null()     ! Pointer to atom that is a bonded neighbour
  integer(i4), dimension(:,:),     pointer, save :: ncnbr_bond => null()     ! Pointer to cell of bonded neighbour 
  real(dp),    dimension(:,:),     pointer, save :: rbnbr => null()          ! Distance to bonded neighbour
  real(dp),    dimension(:,:,:),   pointer, save :: vbnbr => null()          ! Parameters for the interaction with a bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: xbnbr => null()          ! x component of distance to bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: ybnbr => null()          ! y component of distance to bonded neighbour
  real(dp),    dimension(:,:),     pointer, save :: zbnbr => null()          ! z component of distance to bonded neighbour
!
!  Cell parameters lists
!
  integer(i4),                              save :: iimax2
  integer(i4),                              save :: iimid2
  integer(i4),                              save :: ivec2cell(3,125)
  real(dp),                                 save :: xvec2cell(125)
  real(dp),                                 save :: yvec2cell(125)
  real(dp),                                 save :: zvec2cell(125)

CONTAINS

  subroutine changemaxnbr(maxat)
!
!  Changes the size of arrays that hold the neighbour list for pGFNFF
!
!  Julian Gale, CIC, Curtin University, August 2020
!
  use m_io
  use m_pgfnff_reallocate
  implicit none
!
!  Passed variables
!
  integer(i4),  intent(in) :: maxat      ! Maximum number of atoms
!
!  Local variables
!
  integer(i4)              :: ierror
!
  call pgfnff_realloc(nnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr '')')
    stop
  endif
  call pgfnff_realloc(nbrno,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno '')')
    stop
  endif
  call pgfnff_realloc(ncnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ncnbr '')')
    stop
  endif
  call pgfnff_realloc(rnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : rnbr '')')
    stop
  endif
  call pgfnff_realloc(xnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xnbr '')')
    stop
  endif
  call pgfnff_realloc(ynbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ynbr '')')
    stop
  endif
  call pgfnff_realloc(znbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : znbr '')')
    stop
  endif
!
  call pgfnff_realloc(nnbr_cn,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_cn '')')
    stop
  endif
  call pgfnff_realloc(nbrno_cn,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_cn '')')
    stop
  endif
!
  call pgfnff_realloc(nnbr_bond,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nnbr_bond '')')
    stop
  endif
  call pgfnff_realloc(nbrno_bond,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : nbrno_bond '')')
    stop
  endif
  call pgfnff_realloc(ncnbr_bond,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ncnbr_bond '')')
    stop
  endif
  call pgfnff_realloc(rbnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : rbnbr '')')
    stop
  endif
  call pgfnff_realloc(vbnbr,3_i4,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : vbnbr '')')
    stop
  endif
  call pgfnff_realloc(xbnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : xbnbr '')')
    stop
  endif
  call pgfnff_realloc(ybnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : ybnbr '')')
    stop
  endif
  call pgfnff_realloc(zbnbr,maxnbr,maxat,ierror)
  if (ierror.gt.0) then
    write(ioout,'('' Error in Memory Allocation : zbnbr '')')
    stop
  endif
!
  end subroutine changemaxnbr

  subroutine getnbr(ndim,rv,kv,numat,nat,xclat,yclat,zclat,cut2,ldebug)
!
!  Computes the neighbour list for a system
!
!  On exit :
!
!  Neighbour list is set
!
!  Julian Gale, CIC, Curtin University, August 2021
!
  use m_pgfnff_types
  use m_io
  implicit none
!
!  Passed variables
!
  integer(i4), intent(in)                        :: ndim           ! Dimensionality
  integer(i4), intent(in)                        :: numat          ! Number of atoms
  integer(i4), intent(in)                        :: nat(numat)     ! Atomic numbers of atoms
  logical,     intent(in)                        :: ldebug         ! If true then output debug printing
  real(dp),    intent(in)                        :: cut2           ! Cutoff squared
  real(dp),    intent(in)                        :: kv(3,3)        ! Reciprocal space lattice vectors
  real(dp),    intent(in)                        :: rv(3,3)        ! Real space lattice vectors
  real(dp),    intent(in)                        :: xclat(numat)   ! Cartesian x coordinates
  real(dp),    intent(in)                        :: yclat(numat)   ! Cartesian y coordinates
  real(dp),    intent(in)                        :: zclat(numat)   ! Cartesian z coordinates
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ii
  integer(i4)                                    :: j
  integer(i4)                                    :: nn
  logical                                        :: lzeroR
  real(dp)                                       :: r2
  real(dp)                                       :: rij
  real(dp),                                 save :: rzerotol = 1.0d-12
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
!
  cutnbr = sqrt(cut2)
!
!  Set up the cell lists
!
  call set_cell_list(ndim,rv)
!
!  Check that arrays are initialised to the correct size
!
  call changemaxnbr(numat)
!
!  Set number of neighbours to zero
!
  nnbr(1:numat) = 0
  lzeroR = .false.
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  do i = 1,numat
!
!  Loop over atoms
!
    do j = 1,numat
!
!  Set centre cell coordinate differences
!
      xji0 = xclat(j) - xclat(i)
      yji0 = yclat(j) - yclat(i)
      zji0 = zclat(j) - zclat(i)
!
!  Loop over unit cells
!
      do ii = 1,iimax2
!
!  Exclude self term
!
        if (i.ne.j.or.ii.ne.iimid2) then
          xji = xji0 + xvec2cell(ii)
          yji = yji0 + yvec2cell(ii)
          zji = zji0 + zvec2cell(ii)
          r2 = xji*xji + yji*yji + zji*zji
          if (r2.lt.cut2) then
            rij = sqrt(r2)
            if (rij.lt.rzerotol) lzeroR = .true.
            nnbr(i) = nnbr(i) + 1
            if (nnbr(i).ge.maxnbr) then
              maxnbr = maxnbr + 2
              call changemaxnbr(numat)
            endif
            nbrno(nnbr(i),i) = j
            ncnbr(nnbr(i),i) = ii
            rnbr(nnbr(i),i) = rij
            xnbr(nnbr(i),i) = xji
            ynbr(nnbr(i),i) = yji
            znbr(nnbr(i),i) = zji
          endif
        endif
      enddo
    enddo
  enddo
!*******************
!  Debug printing  *
!*******************
  if (ldebug) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do i = 1,numat
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nnbr(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      write(ioout,'(i4,8(1x,i4))') i,(nbrno(nn,i),nn=1,nnbr(i))
    enddo
  endif
!
!  Trap zero distances
!
  if (lzeroR) then
    call pgfnff_error('atoms are too close together',0_i4)
    stop
  endif
!
  end subroutine getnbr

  subroutine set_cell_list(ndim,rv)
!
!  Store linear array of lattice vectors for ii/jj/kk=-2,2 
!  => 125 lattice vectors. This tidies up and speeds up the
!  loops over these indices.
!
!  Julian Gale, CIC, Curtin University, October 2021
!
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4),    intent(in)  :: ndim
  real(dp),       intent(in)  :: rv(3,3)
!
!  Local variables
!
  integer(i4)                 :: ii
  integer(i4)                 :: imaxl2
  integer(i4)                 :: jj
  integer(i4)                 :: jmaxl2
  integer(i4)                 :: kk
  integer(i4)                 :: kmaxl2
  real(dp)                    :: xcdi
  real(dp)                    :: ycdi
  real(dp)                    :: zcdi
  real(dp)                    :: xcdj
  real(dp)                    :: ycdj
  real(dp)                    :: zcdj
  real(dp)                    :: xcrd
  real(dp)                    :: ycrd
  real(dp)                    :: zcrd
!
  if (ndim.eq.3) then
!*************
!  3-D case  *
!*************
    imaxl2 = 2_i4
    jmaxl2 = 2_i4
    kmaxl2 = 2_i4
    xcdi = -3.0_dp*rv(1,1)
    ycdi = -3.0_dp*rv(2,1)
    zcdi = -3.0_dp*rv(3,1)
    iimax2 = 0
!
!  Loop over unit cells
!
    do ii = -2,2
      xcdi = xcdi + rv(1,1)
      ycdi = ycdi + rv(2,1)
      zcdi = zcdi + rv(3,1)
      xcdj = xcdi - 3.0_dp*rv(1,2)
      ycdj = ycdi - 3.0_dp*rv(2,2)
      zcdj = zcdi - 3.0_dp*rv(3,2)
      do jj = -2,2
        xcdj = xcdj + rv(1,2)
        ycdj = ycdj + rv(2,2)
        zcdj = zcdj + rv(3,2)
        xcrd = xcdj - 3.0_dp*rv(1,3)
        ycrd = ycdj - 3.0_dp*rv(2,3)
        zcrd = zcdj - 3.0_dp*rv(3,3)
        do kk = -2,2
          iimax2 = iimax2 + 1
          xcrd = xcrd + rv(1,3)
          ycrd = ycrd + rv(2,3)
          zcrd = zcrd + rv(3,3)
          ivec2cell(1,iimax2) = ii
          ivec2cell(2,iimax2) = jj
          ivec2cell(3,iimax2) = kk
          xvec2cell(iimax2) = xcrd
          yvec2cell(iimax2) = ycrd
          zvec2cell(iimax2) = zcrd
        enddo
      enddo
    enddo
    iimid2 = 63
  elseif (ndim.eq.2) then
!*************
!  2-D case  *
!*************
    imaxl2 = 2_i4
    jmaxl2 = 2_i4
    kmaxl2 = 0_i4
    xcdi = -3.0_dp*rv(1,1)
    ycdi = -3.0_dp*rv(2,1)
    iimax2 = 0
!
!  Loop over unit cells
!
    do ii = -2,2
      xcdi = xcdi + rv(1,1)
      ycdi = ycdi + rv(2,1)
      xcrd = xcdi - 3.0_dp*rv(1,2)
      ycrd = ycdi - 3.0_dp*rv(2,2)
      do jj = -2,2
        xcrd = xcrd + rv(1,2)
        ycrd = ycrd + rv(2,2)
        iimax2 = iimax2 + 1
        ivec2cell(1,iimax2) = ii
        ivec2cell(2,iimax2) = jj
        ivec2cell(3,iimax2) = 0
        xvec2cell(iimax2) = xcrd
        yvec2cell(iimax2) = ycrd
        zvec2cell(iimax2) = 0.0_dp
      enddo
    enddo
    iimid2 = 13
  elseif (ndim.eq.1) then
!*************
!  1-D case  *
!*************
    imaxl2 = 2_i4
    jmaxl2 = 0_i4
    kmaxl2 = 0_i4
    xcrd = -3.0_dp*rv(1,1)
    iimax2 = 0
!
!  Loop over unit cells
!
    do ii = -2,2
      xcrd = xcrd + rv(1,1)
      iimax2 = iimax2 + 1
      ivec2cell(1,iimax2) = ii
      ivec2cell(2,iimax2) = 0
      ivec2cell(3,iimax2) = 0
      xvec2cell(iimax2) = xcrd
      yvec2cell(iimax2) = 0.0_dp
      zvec2cell(iimax2) = 0.0_dp
    enddo
    iimid2 = 3
  else
!*************
!  0-D case  *
!*************
    imaxl2 = 0_i4
    jmaxl2 = 0_i4
    kmaxl2 = 0_i4
    iimax2 = 1
    iimid2 = 1
    ivec2cell(1,iimax2) = 0
    ivec2cell(2,iimax2) = 0
    ivec2cell(3,iimax2) = 0
    xvec2cell(1) = 0.0_dp
    yvec2cell(1) = 0.0_dp
    zvec2cell(1) = 0.0_dp
  endif
!
  return
  end

end module m_nbr
