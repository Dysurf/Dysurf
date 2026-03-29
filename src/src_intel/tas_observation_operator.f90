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

! Unified observed-spectrum operator for TAS resolution treatments
module tas_observation_operator

  use func, only: gauss
  use tas_resolution_cn, only: invert3

  implicit none

  real(kind=8), parameter :: tpi_local = acos(-1.d0) * 2.d0

contains

  ! Unified observed-spectrum entry point.
  ! linearized_4d: local linearized branch + full RM(4,4)
  subroutine calc_sobs(q0_cart, omega0, intensity0, e0, R0, RM, sobs, mode, ierr, grad_sample, &
                       q0_hkl, sample_x, sample_y, sample_z, rmsd_local, temp_local, branch_id)

    implicit none

    real(kind=8), intent(in) :: q0_cart(:)
    real(kind=8), intent(in) :: omega0, intensity0, e0, R0, RM(4,4)
    real(kind=8), intent(out) :: sobs
    character(len=*), intent(in) :: mode
    integer(kind=4), intent(out) :: ierr
    real(kind=8), intent(in), optional :: grad_sample(3)
    real(kind=8), intent(in), optional :: q0_hkl(3), sample_x(3), sample_y(3), sample_z(3)
    real(kind=8), intent(in), optional :: rmsd_local(:,:), temp_local
    integer(kind=4), intent(in), optional :: branch_id

    ierr = 0
    sobs = 0.d0

    if (size(q0_cart) /= 3) then
       ierr = 1
       return
    end if

    if (trim(mode) == 'linearized_4d') then
       if (.not. present(grad_sample)) then
          ierr = 2
          return
       end if
       call calc_sobs_linearized_4d(omega0, intensity0, e0, R0, RM, grad_sample, sobs, ierr)
    elseif (trim(mode) == 'psf2d') then
       ierr = 6
    else
       ierr = 4
    end if

  end subroutine calc_sobs

  subroutine calc_sobs_linearized_4d(omega0, intensity0, e0, R0, RM, grad_sample, sobs, ierr)

    implicit none

    real(kind=8), intent(in) :: omega0, intensity0, e0, R0, RM(4,4), grad_sample(3)
    real(kind=8), intent(out) :: sobs
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: info
    real(kind=8) :: A(3,3), Ainv(3,3), Mlin(3,3), Minv(3,3), b(3), c(3)
    real(kind=8) :: detA, detM, d, schur0, schur_lin, sigma_lin, scale

    ierr = 0
    sobs = 0.d0

    A = 0.5d0 * (RM(1:3, 1:3) + transpose(RM(1:3, 1:3)))
    b = RM(1:3, 4)
    d = RM(4, 4)

    call invert3(A, Ainv, detA, info)
    if (info /= 0 .or. detA <= 0.d0) then
      ierr = 11
      return
    end if

    schur0 = d - dot_product(b, matmul(Ainv, b))
    if (schur0 <= 0.d0) then
      ierr = 12
      return
    end if

    Mlin = A - spread(b, 2, 3) * spread(grad_sample, 1, 3) - spread(grad_sample, 2, 3) * spread(b, 1, 3) + &
           d * matmul(reshape(grad_sample, (/3, 1/)), reshape(grad_sample, (/1, 3/)))
    Mlin = 0.5d0 * (Mlin + transpose(Mlin))
    c = b - d * grad_sample

    call invert3(Mlin, Minv, detM, info)
    if (info /= 0 .or. detM <= 0.d0) then
      ierr = 13
      return
    end if

    schur_lin = d - dot_product(c, matmul(Minv, c))
    if (schur_lin <= 0.d0) then
      ierr = 14
      return
    end if

    sigma_lin = 1.d0 / sqrt(schur_lin)
    if (R0 > 0.d0) then
      scale = (R0 * sqrt(detA / detM)) / R0
    else
      scale = sqrt(detA / detM)
    end if

    sobs = intensity0 * scale * gauss(omega0, e0, sigma_lin)

  end subroutine calc_sobs_linearized_4d

end module tas_observation_operator
