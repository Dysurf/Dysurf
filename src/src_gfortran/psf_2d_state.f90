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

module psf_2d_state

  implicit none

  real(kind=8), allocatable :: q_slice_grid(:)
  real(kind=8), allocatable :: e_slice_grid(:)
  real(kind=8), allocatable :: intrinsic_slice(:,:,:)
  real(kind=8), allocatable :: psf2d_slice(:,:,:)
  real(kind=8), allocatable :: sigma_q_field(:,:,:)
  real(kind=8), allocatable :: sigma_e_left_field(:,:,:)
  real(kind=8), allocatable :: sigma_e_right_field(:,:,:)
  real(kind=8), allocatable :: shear_field(:,:,:)

contains

  subroutine psf2d_setup_state(ne, nqh, ntemps, qgrid, egrid)

    implicit none

    integer(kind=4), intent(in) :: ne, nqh, ntemps
    real(kind=8), intent(in) :: qgrid(:), egrid(:)

    call psf2d_free_state()
    allocate(q_slice_grid(nqh), e_slice_grid(ne))
    allocate(intrinsic_slice(ne, nqh, ntemps), psf2d_slice(ne, nqh, ntemps))
    allocate(sigma_q_field(ne, nqh, ntemps), sigma_e_left_field(ne, nqh, ntemps))
    allocate(sigma_e_right_field(ne, nqh, ntemps), shear_field(ne, nqh, ntemps))

    q_slice_grid = qgrid
    e_slice_grid = egrid
    intrinsic_slice = 0.d0
    psf2d_slice = 0.d0
    sigma_q_field = 0.d0
    sigma_e_left_field = 0.d0
    sigma_e_right_field = 0.d0
    shear_field = 0.d0

  end subroutine psf2d_setup_state

  subroutine psf2d_free_state()

    implicit none

    if (allocated(q_slice_grid)) deallocate(q_slice_grid)
    if (allocated(e_slice_grid)) deallocate(e_slice_grid)
    if (allocated(intrinsic_slice)) deallocate(intrinsic_slice)
    if (allocated(psf2d_slice)) deallocate(psf2d_slice)
    if (allocated(sigma_q_field)) deallocate(sigma_q_field)
    if (allocated(sigma_e_left_field)) deallocate(sigma_e_left_field)
    if (allocated(sigma_e_right_field)) deallocate(sigma_e_right_field)
    if (allocated(shear_field)) deallocate(shear_field)

  end subroutine psf2d_free_state

end module psf_2d_state
