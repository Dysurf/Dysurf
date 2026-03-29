! Dysurf, a program for simulating four-dimensional dynamical structure factors
! Copyright (C) 2023-2025 Yongheng Li <davy_li96@163.com>
! Copyright (C) 2020-2021 Changpeng Lin <changpeng.lin@epfl.ch>
! Copyright (C) 2020-2021 Jiawang Hong <hongjw@bit.edu.cn>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
module mpi_helper

  use mpi
  implicit none

  integer(kind=4) :: mpi_rank = 0
  integer(kind=4) :: mpi_size = 1
  logical :: mpi_is_initialized = .false.

contains

  subroutine mpi_helper_init()

    implicit none

    integer(kind=4) :: ierr
    logical :: flag

    call MPI_Initialized(flag, ierr)
    if (.not. flag) then
       call MPI_Init(ierr)
       mpi_is_initialized = .true.
    else
       mpi_is_initialized = .false.
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr)

  end subroutine mpi_helper_init

  subroutine mpi_helper_finalize()

    implicit none

    integer(kind=4) :: ierr
    logical :: flag

    call MPI_Finalized(flag, ierr)
    if (.not. flag .and. mpi_is_initialized) then
       call MPI_Finalize(ierr)
    end if

  end subroutine mpi_helper_finalize

  subroutine mpi_allreduce_sum_real8(arr)

    implicit none

    real(kind=8), intent(inout) :: arr(:)
    real(kind=8), allocatable :: work(:)
    integer(kind=4) :: ierr

    if (mpi_size <= 1) return

    allocate(work(size(arr)))
    call MPI_Allreduce(arr, work, size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    arr = work
    deallocate(work)

  end subroutine mpi_allreduce_sum_real8

  subroutine mpi_allreduce_sum_real8_5d(arr)

    implicit none

    real(kind=8), intent(inout) :: arr(:,:,:,:,:)
    real(kind=8), allocatable :: work(:,:,:,:,:)
    integer(kind=4) :: ierr

    if (mpi_size <= 1) return

    allocate(work(size(arr,1), size(arr,2), size(arr,3), size(arr,4), size(arr,5)))
    call MPI_Allreduce(arr, work, size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    arr = work
    deallocate(work)

  end subroutine mpi_allreduce_sum_real8_5d

  subroutine mpi_allreduce_sum_real8_3d(arr)

    implicit none

    real(kind=8), intent(inout) :: arr(:,:,:)
    real(kind=8), allocatable :: work(:,:,:)
    integer(kind=4) :: ierr

    if (mpi_size <= 1) return

    allocate(work(size(arr,1), size(arr,2), size(arr,3)))
    call MPI_Allreduce(arr, work, size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    arr = work
    deallocate(work)

  end subroutine mpi_allreduce_sum_real8_3d

end module mpi_helper
