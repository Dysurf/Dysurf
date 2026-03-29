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

! Local analytic 2D PSF kernel utilities.
module psf_2d_kernel

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

contains

  ! Return an asymmetric Gaussian value evaluated at dE.
  !
  ! Input:
  !   de                  : energy offset after shear correction
  !   sigma_e_left        : width used when de < 0
  !   sigma_e_right       : width used when de >= 0
  !
  ! Output:
  !   asymmetric_gaussian : unnormalized Gaussian value
  real(kind=8) function asymmetric_gaussian(de, sigma_e_left, sigma_e_right)

    implicit none

    real(kind=8), intent(in) :: de, sigma_e_left, sigma_e_right
    real(kind=8) :: sigma_use

    if (de < 0.d0) then
       sigma_use = sigma_e_left
    else
       sigma_use = sigma_e_right
    end if

    if (sigma_use <= 0.d0 .or. .not. ieee_is_finite(sigma_use)) then
       asymmetric_gaussian = 0.d0
       return
    end if

    asymmetric_gaussian = exp(-0.5d0 * (de / sigma_use) ** 2)

  end function asymmetric_gaussian

  ! Normalize a local 2D kernel so that sum(kernel) = 1.
  !
  ! Input/Output:
  !   kernel(:,:)         : local 2D PSF kernel
  !
  ! Output:
  !   ierr                : 0 on success
  subroutine normalize_kernel(kernel, ierr)

    implicit none

    real(kind=8), intent(inout) :: kernel(:,:)
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: sumk

    ierr = 0
    sumk = sum(kernel)
    if (.not. ieee_is_finite(sumk) .or. sumk <= 0.d0) then
       ierr = 1
       return
    end if

    kernel = kernel / sumk

  end subroutine normalize_kernel

  ! Placeholder parameter interpolation interface.
  ! Current version supports constant parameters only and simply returns them.
  !
  ! Input:
  !   sigma_q_in, sigma_e_left_in, sigma_e_right_in, shear_in : local PSF parameters
  !
  ! Output:
  !   sigma_q_out, sigma_e_left_out, sigma_e_right_out, shear_out
  subroutine interpolate_psf_params(sigma_q_in, sigma_e_left_in, sigma_e_right_in, shear_in, &
                                    sigma_q_out, sigma_e_left_out, sigma_e_right_out, shear_out, ierr)

    implicit none

    real(kind=8), intent(in) :: sigma_q_in, sigma_e_left_in, sigma_e_right_in, shear_in
    real(kind=8), intent(out) :: sigma_q_out, sigma_e_left_out, sigma_e_right_out, shear_out
    integer(kind=4), intent(out) :: ierr

    ierr = 0
    sigma_q_out = sigma_q_in
    sigma_e_left_out = sigma_e_left_in
    sigma_e_right_out = sigma_e_right_in
    shear_out = shear_in

  end subroutine interpolate_psf_params

  ! Build a local analytic 2D PSF kernel on a discrete (dq,dE) window.
  !
  ! Input:
  !   q_offsets(:)        : local q-window offsets in bin units or physical units
  !   e_offsets(:)        : local energy-window offsets in bin units or physical units
  !   sigma_q             : q Gaussian width, same unit as q_offsets
  !   sigma_e_left        : left-side energy width, same unit as e_offsets
  !   sigma_e_right       : right-side energy width, same unit as e_offsets
  !   shear               : q-E tilt, satisfying dE_eff = dE - shear * dQ
  !
  ! Output:
  !   kernel(ne,nq)       : normalized local 2D PSF kernel
  !   ierr                : 0 on success
  subroutine build_local_psf_kernel(q_offsets, e_offsets, sigma_q, sigma_e_left, sigma_e_right, shear, kernel, ierr)

    implicit none

    real(kind=8), intent(in) :: q_offsets(:), e_offsets(:)
    real(kind=8), intent(in) :: sigma_q, sigma_e_left, sigma_e_right, shear
    real(kind=8), allocatable, intent(out) :: kernel(:,:)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: iq, ie, nq, ne
    real(kind=8) :: q_part, de_eff

    ierr = 0
    nq = size(q_offsets)
    ne = size(e_offsets)

    if (nq <= 0 .or. ne <= 0) then
       ierr = 1
       return
    end if
    if (sigma_q <= 0.d0 .or. sigma_e_left <= 0.d0 .or. sigma_e_right <= 0.d0) then
       ierr = 2
       return
    end if
    if (.not. ieee_is_finite(sigma_q) .or. .not. ieee_is_finite(sigma_e_left) .or. &
        .not. ieee_is_finite(sigma_e_right) .or. .not. ieee_is_finite(shear)) then
       ierr = 3
       return
    end if

    allocate(kernel(ne, nq))
    kernel = 0.d0

    do iq = 1, nq
       q_part = exp(-0.5d0 * (q_offsets(iq) / sigma_q) ** 2)
       do ie = 1, ne
          de_eff = e_offsets(ie) - shear * q_offsets(iq)
          kernel(ie, iq) = q_part * asymmetric_gaussian(de_eff, sigma_e_left, sigma_e_right)
       end do
    end do

    call normalize_kernel(kernel, ierr)
    if (ierr /= 0) then
       if (allocated(kernel)) deallocate(kernel)
       ierr = 10 + ierr
    end if

  end subroutine build_local_psf_kernel

  ! Minimal selfcheck:
  ! 1. kernel sum should be approximately 1
  ! 2. asymmetric widths should produce a broader positive-energy side
  ! 3. positive shear should shift the peak toward larger dE at positive dQ
  subroutine selfcheck_psf_2d_kernel(ierr)

    implicit none

    integer(kind=4), intent(out) :: ierr

    real(kind=8), parameter :: q_offsets(5) = (/-2.d0, -1.d0, 0.d0, 1.d0, 2.d0/)
    real(kind=8), parameter :: e_offsets(9) = (/-4.d0, -3.d0, -2.d0, -1.d0, 0.d0, 1.d0, 2.d0, 3.d0, 4.d0/)
    real(kind=8), allocatable :: kernel(:,:)
    integer(kind=4) :: iq_left, iq_right, ie_peak_left, ie_peak_right, ie0, i
    real(kind=8) :: sumk, left_width, right_width

    ierr = 0

    call build_local_psf_kernel(q_offsets, e_offsets, 1.d0, 0.8d0, 1.6d0, 0.5d0, kernel, ierr)
    if (ierr /= 0) return

    sumk = sum(kernel)
    if (abs(sumk - 1.d0) > 1.d-12) then
       ierr = 11
       return
    end if

    ie0 = 5
    left_width = 0.d0
    right_width = 0.d0
    do i = 1, size(e_offsets)
       if (e_offsets(i) < 0.d0) left_width = left_width + kernel(i, 3)
       if (e_offsets(i) > 0.d0) right_width = right_width + kernel(i, 3)
    end do
    if (right_width <= left_width) then
       ierr = 12
       return
    end if

    iq_left = 2
    iq_right = 4
    ie_peak_left = maxloc(kernel(:, iq_left), dim=1)
    ie_peak_right = maxloc(kernel(:, iq_right), dim=1)
    if (e_offsets(ie_peak_right) <= e_offsets(ie_peak_left)) then
       ierr = 13
       return
    end if

    deallocate(kernel)

  end subroutine selfcheck_psf_2d_kernel

end module psf_2d_kernel
