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
  module m_pgfnff_reallocate
!
!  Set of routines and interface for handling memory
!  allocation and reallocation while preserving the
!  contents.
!
!  Julian Gale, CIC, Curtin University, October 2020
!
    use m_pgfnff_types
!
!  Local variables for storing memory usage information
!
    integer(i4), save :: losize  = 4
    integer(i4), save :: i4size  = 4
    integer(i4), save :: r8size  = 8
    integer(i4), save :: c16size = 16

    interface pgfnff_realloc

      module procedure &
        realloc_r8_1, &
        realloc_i4_1, &
        realloc_c16_1, &
        realloc_l_1, &
        realloc_r8_2, &
        realloc_i4_2, &
        realloc_c16_2, &
        realloc_l_2, &
        realloc_r8_3, &
        realloc_i4_3, &
        realloc_c16_3, &
        realloc_l_3, &
        realloc_r8_4, &
        realloc_i4_4, &
        realloc_c16_4, &
        realloc_l_4

    end interface 

    private :: realloc_l_1,   realloc_l_2,   realloc_l_3,   realloc_l_4
    private :: realloc_i4_1,  realloc_i4_2,  realloc_i4_3,  realloc_i4_4
    private :: realloc_r8_1,  realloc_r8_2,  realloc_r8_3,  realloc_r8_4
    private :: realloc_c16_1, realloc_c16_2, realloc_c16_3, realloc_c16_4

  contains

    subroutine realloc_r8_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      real(dp), dimension(:), pointer :: p
      integer(i4), intent(in)         :: n
      integer(i4), intent(out)        :: ierror
!
!  Local arguments
!
      real(dp), dimension(:), pointer :: plocal
      integer(i4)                     :: nold
      integer(i4)                     :: i
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0.0_dp
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          plocal(i) = p(i)
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_r8_1

    subroutine realloc_i4_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), dimension(:), pointer :: p
      integer(i4), intent(in)            :: n
      integer(i4), intent(out)           :: ierror
!
!  Local arguments
!
      integer(i4), dimension(:), pointer :: plocal
      integer(i4)                        :: nold
      integer(i4)                        :: i
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0_i4
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          plocal(i) = p(i)
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_i4_1

    subroutine realloc_c16_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      complex(dpc), dimension(:), pointer :: p
      integer(i4), intent(in)             :: n
      integer(i4), intent(out)            :: ierror
!
!  Local arguments
!
      complex(dpc), dimension(:), pointer :: plocal
      integer(i4)                         :: nold
      integer(i4)                         :: i
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = (0.0_dp,0.0_dp)
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          plocal(i) = p(i)
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_c16_1

    subroutine realloc_l_1(p,n,ierror)

      implicit none
!
!  Passed arguments
!
      logical, dimension(:), pointer :: p
      integer(i4), intent(in)        :: n
      integer(i4), intent(out)       :: ierror
!
!  Local arguments
!
      logical, dimension(:), pointer :: plocal
      integer(i4)                    :: nold
      integer(i4)                    :: i
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = .false.
!
      if (associated(p)) then
!
!  Find size of old array
!
        nold = size(p)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          plocal(i) = p(i)
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_l_1

    subroutine realloc_r8_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      real(dp), dimension(:,:), pointer :: p
      integer(i4), intent(in)           :: m,n
      integer(i4), intent(out)          :: ierror
!
!  Local arguments
!
      real(dp), dimension(:,:), pointer :: plocal
      integer(i4)                       :: mold,nold
      integer(i4)                       :: i, j
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0.0_dp
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            plocal(j,i) = p(j,i)
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_r8_2

    subroutine realloc_i4_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), dimension(:,:), pointer :: p
      integer(i4), intent(in)              :: m,n
      integer(i4), intent(out)             :: ierror
!
!  Local arguments
!
      integer(i4), dimension(:,:), pointer :: plocal
      integer(i4)                          :: mold,nold
      integer(i4)                          :: i, j
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0_i4
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            plocal(j,i) = p(j,i)
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_i4_2

    subroutine realloc_c16_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      complex(dpc), dimension(:,:), pointer :: p
      integer(i4), intent(in)               :: m,n
      integer(i4), intent(out)              :: ierror
!
!  Local arguments
!
      complex(dpc), dimension(:,:), pointer :: plocal
      integer(i4)                           :: mold,nold
      integer(i4)                           :: i, j
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = (0.0_dp,0.0_dp)
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            plocal(j,i) = p(j,i)
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_c16_2

    subroutine realloc_l_2(p,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      logical, dimension(:,:), pointer :: p
      integer(i4), intent(in)          :: m,n
      integer(i4), intent(out)         :: ierror
!
!  Local arguments
!
      logical, dimension(:,:), pointer :: plocal
      integer(i4)                      :: mold,nold
      integer(i4)                      :: i, j
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = .false.
!
      if (associated(p)) then
!
!  Find size of old array
!
        mold = size(p,1)
        nold = size(p,2)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            plocal(j,i) = p(j,i)
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_l_2

    subroutine realloc_r8_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      real(dp), dimension(:,:,:), pointer :: p
      integer(i4), intent(in)             :: k,m,n
      integer(i4), intent(out)            :: ierror
!
!  Local arguments
!
      real(dp), dimension(:,:,:), pointer :: plocal
      integer(i4)                         :: kold,mold,nold
      integer(i4)                         :: i, j, l
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0.0_dp
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            do l = 1,min(kold,k)
              plocal(l,j,i) = p(l,j,i)
            enddo
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_r8_3

    subroutine realloc_i4_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), dimension(:,:,:), pointer :: p
      integer(i4), intent(in)                :: k,m,n
      integer(i4), intent(out)               :: ierror
!
!  Local arguments
!
      integer(i4), dimension(:,:,:), pointer :: plocal
      integer(i4)                            :: kold,mold,nold
      integer(i4)                            :: i, j, l
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0_i4
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            do l = 1,min(kold,k)
              plocal(l,j,i) = p(l,j,i)
            enddo
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_i4_3

    subroutine realloc_c16_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      complex(dpc), dimension(:,:,:), pointer :: p
      integer(i4), intent(in)                 :: k,m,n
      integer(i4), intent(out)                :: ierror
!
!  Local arguments
!
      complex(dpc), dimension(:,:,:), pointer :: plocal
      integer(i4)                             :: kold,mold,nold
      integer(i4)                             :: i, j, l
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = (0.0_dp,0.0_dp)
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            do l = 1,min(kold,k)
              plocal(l,j,i) = p(l,j,i)
            enddo
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_c16_3

    subroutine realloc_l_3(p,k,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      logical, dimension(:,:,:), pointer :: p
      integer(i4), intent(in)            :: k,m,n
      integer(i4), intent(out)           :: ierror
!
!  Local arguments
!
      logical, dimension(:,:,:), pointer :: plocal
      integer(i4)                        :: kold,mold,nold
      integer(i4)                        :: i, j, l
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = .false.
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        mold = size(p,2)
        nold = size(p,3)
!
!  Preserve data
!
        do i = 1,min(nold,n)
          do j = 1,min(mold,m)
            do l = 1,min(kold,k)
              plocal(l,j,i) = p(l,j,i)
            enddo
          enddo
        enddo
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_l_3

    subroutine realloc_r8_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      real(dp), dimension(:,:,:,:), pointer :: p
      integer(i4), intent(in)               :: k,l,m,n
      integer(i4), intent(out)              :: ierror
!
!  Local arguments
!
      real(dp), dimension(:,:,:,:), pointer :: plocal
      integer(i4)                           :: kold,lold
      integer(i4)                           :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0.0_dp
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_r8_4

    subroutine realloc_i4_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      integer(i4), dimension(:,:,:,:), pointer :: p
      integer(i4), intent(in)                  :: k,l,m,n
      integer(i4), intent(out)                 :: ierror
!
!  Local arguments
!
      integer(i4), dimension(:,:,:,:), pointer :: plocal
      integer(i4)                              :: kold,lold
      integer(i4)                              :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = 0_i4
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_i4_4

    subroutine realloc_c16_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      complex(dpc), dimension(:,:,:,:), pointer :: p
      integer(i4), intent(in)                   :: k,l,m,n
      integer(i4), intent(out)                  :: ierror
!
!  Local arguments
!
      complex(dpc), dimension(:,:,:,:), pointer :: plocal
      integer(i4)                               :: kold,lold
      integer(i4)                               :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = (0.0_dp,0.0_dp)
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_c16_4

    subroutine realloc_l_4(p,k,l,m,n,ierror)

      implicit none
!
!  Passed arguments
!
      logical, dimension(:,:,:,:), pointer :: p
      integer(i4), intent(in)              :: k,l,m,n
      integer(i4), intent(out)             :: ierror
!
!  Local arguments
!
      logical, dimension(:,:,:,:), pointer :: plocal
      integer(i4)                          :: kold,lold
      integer(i4)                          :: mold,nold
!
!  Allocate new size of memory
!
      nullify(plocal)
      allocate(plocal(k,l,m,n),stat=ierror)
!
!  If there is an error then return
!
      if (ierror.ne.0) then
        return
      endif
      plocal = .false.
!
      if (associated(p)) then
!
!  Find size of old array
!
        kold = size(p,1)
        lold = size(p,2)
        mold = size(p,3)
        nold = size(p,4)
!
!  Preserve data
!
        plocal(1:min(kold,k),1:min(lold,l),1:min(mold,m), &
          1:min(nold,n)) = p(1:min(kold,k),1:min(lold,l), &
          1:min(mold,m),1:min(nold,n))
!
!  Free old memory
!
        deallocate(p)
      endif
!
!  Re-assign pointers so that p is returned as the new array
!
      p=>plocal

    end subroutine realloc_l_4

  end module m_pgfnff_reallocate
