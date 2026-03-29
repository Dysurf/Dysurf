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

! Position-dependent local 2D convolution on an intrinsic/model slice.
! This module is standalone and not connected to the main flow.
module psf_2d_convolution

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use psf_2d_kernel, only: build_local_psf_kernel

  implicit none

contains

  ! Extract the local patch of the input slice and the corresponding kernel axes.
  !
  ! Input:
  !   input_slice(ne, nq) : intrinsic/model slice to be convolved
  !   iq_center           : output q-pixel index
  !   ie_center           : output energy-pixel index
  !   window_q            : half-window in q pixels
  !   window_e            : half-window in energy pixels
  !
  ! Output:
  !   patch(:,:)          : extracted neighborhood patch
  !   q_offsets(:)        : local q offsets relative to (iq_center, ie_center)
  !   e_offsets(:)        : local energy offsets relative to (iq_center, ie_center)
  !   iq_lo, iq_hi        : clipped q-index bounds in the input slice
  !   ie_lo, ie_hi        : clipped energy-index bounds in the input slice
  !   ierr                : 0 on success
  subroutine extract_local_patch(input_slice, iq_center, ie_center, window_q, window_e, &
                                 patch, q_offsets, e_offsets, iq_lo, iq_hi, ie_lo, ie_hi, ierr)

    implicit none

    real(kind=8), intent(in) :: input_slice(:,:)
    integer(kind=4), intent(in) :: iq_center, ie_center, window_q, window_e
    real(kind=8), allocatable, intent(out) :: patch(:,:), q_offsets(:), e_offsets(:)
    integer(kind=4), intent(out) :: iq_lo, iq_hi, ie_lo, ie_hi, ierr

    integer(kind=4) :: nq, ne, iq, ie

    ierr = 0
    ne = size(input_slice, 1)
    nq = size(input_slice, 2)

    if (iq_center < 1 .or. iq_center > nq .or. ie_center < 1 .or. ie_center > ne) then
       ierr = 1
       return
    end if
    if (window_q < 0 .or. window_e < 0) then
       ierr = 2
       return
    end if

    iq_lo = max(1, iq_center - window_q)
    iq_hi = min(nq, iq_center + window_q)
    ie_lo = max(1, ie_center - window_e)
    ie_hi = min(ne, ie_center + window_e)

    allocate(patch(ie_hi - ie_lo + 1, iq_hi - iq_lo + 1))
    allocate(q_offsets(iq_hi - iq_lo + 1))
    allocate(e_offsets(ie_hi - ie_lo + 1))

    patch = input_slice(ie_lo:ie_hi, iq_lo:iq_hi)

    do iq = iq_lo, iq_hi
       q_offsets(iq - iq_lo + 1) = dble(iq - iq_center)
    end do
    do ie = ie_lo, ie_hi
       e_offsets(ie - ie_lo + 1) = dble(ie - ie_center)
    end do

  end subroutine extract_local_patch

  ! Apply a normalized local kernel to a local patch and return the convolved scalar.
  !
  ! Input:
  !   patch(:,:)          : local intrinsic/model patch
  !   kernel(:,:)         : local normalized PSF kernel, same shape as patch
  !
  ! Output:
  !   value               : convolved output at the center pixel
  !   ierr                : 0 on success
  subroutine apply_local_kernel(patch, kernel, value, ierr)

    implicit none

    real(kind=8), intent(in) :: patch(:,:), kernel(:,:)
    real(kind=8), intent(out) :: value
    integer(kind=4), intent(out) :: ierr

    ierr = 0
    value = 0.d0

    if (size(patch, 1) /= size(kernel, 1) .or. size(patch, 2) /= size(kernel, 2)) then
       ierr = 1
       return
    end if

    value = sum(patch * kernel)
    if (.not. ieee_is_finite(value)) then
       ierr = 2
       value = 0.d0
    end if

  end subroutine apply_local_kernel

  ! Convolve a 2D slice with a position-dependent local PSF kernel.
  !
  ! Input:
  !   input_slice(ne, nq)             : intrinsic/model slice
  !   sigma_q_field(ne, nq)           : q-width field
  !   sigma_e_left_field(ne, nq)      : left energy-width field
  !   sigma_e_right_field(ne, nq)     : right energy-width field
  !   shear_field(ne, nq)             : shear field
  !   window_q                        : half-window in q pixels
  !   window_e                        : half-window in energy pixels
  !
  ! Output:
  !   output_slice(ne, nq)            : PSF-convolved slice
  !   ierr                            : 0 on success
  subroutine convolve_slice_with_local_psf(input_slice, sigma_q_field, sigma_e_left_field, sigma_e_right_field, &
                                           shear_field, window_q, window_e, output_slice, ierr)

    implicit none

    real(kind=8), intent(in) :: input_slice(:,:), sigma_q_field(:,:), sigma_e_left_field(:,:), &
                                sigma_e_right_field(:,:), shear_field(:,:)
    integer(kind=4), intent(in) :: window_q, window_e
    real(kind=8), allocatable, intent(out) :: output_slice(:,:)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: nq, ne, iq, ie, iq_lo, iq_hi, ie_lo, ie_hi, info
    real(kind=8), allocatable :: patch(:,:), kernel(:,:), q_offsets(:), e_offsets(:)

    ierr = 0
    ne = size(input_slice, 1)
    nq = size(input_slice, 2)

    if (size(sigma_q_field, 1) /= ne .or. size(sigma_q_field, 2) /= nq) then
       ierr = 1
       return
    end if
    if (size(sigma_e_left_field, 1) /= ne .or. size(sigma_e_left_field, 2) /= nq) then
       ierr = 2
       return
    end if
    if (size(sigma_e_right_field, 1) /= ne .or. size(sigma_e_right_field, 2) /= nq) then
       ierr = 3
       return
    end if
    if (size(shear_field, 1) /= ne .or. size(shear_field, 2) /= nq) then
       ierr = 4
       return
    end if

    allocate(output_slice(ne, nq))
    output_slice = 0.d0

    do iq = 1, nq
       do ie = 1, ne
          call extract_local_patch(input_slice, iq, ie, window_q, window_e, patch, q_offsets, e_offsets, &
                                   iq_lo, iq_hi, ie_lo, ie_hi, info)
          if (info /= 0) then
             ierr = 10 + info
             deallocate(output_slice)
             return
          end if

          call build_local_psf_kernel(q_offsets, e_offsets, sigma_q_field(ie, iq), sigma_e_left_field(ie, iq), &
                                      sigma_e_right_field(ie, iq), shear_field(ie, iq), kernel, info)
          if (info /= 0) then
             ierr = 20 + info
             deallocate(patch, q_offsets, e_offsets)
             deallocate(output_slice)
             return
          end if

          call apply_local_kernel(patch, kernel, output_slice(ie, iq), info)
          if (info /= 0) then
             ierr = 30 + info
             deallocate(patch, q_offsets, e_offsets, kernel)
             deallocate(output_slice)
             return
          end if

          deallocate(patch, q_offsets, e_offsets, kernel)
       end do
    end do

  end subroutine convolve_slice_with_local_psf

  ! Minimal selfcheck:
  ! 1. delta-like input should reproduce a local kernel shape around the center
  ! 2. constant parameter field should produce a smooth finite output
  ! 3. position-dependent parameter field should broaden one side more than the other
  subroutine selfcheck_psf_2d_convolution(ierr)

    implicit none

    integer(kind=4), intent(out) :: ierr

    integer(kind=4), parameter :: ne = 9, nq = 7
    integer(kind=4) :: ie_peak_left, ie_peak_right
    real(kind=8) :: input_slice(ne, nq), sq(ne, nq), sel(ne, nq), ser(ne, nq), sh(ne, nq)
    real(kind=8), allocatable :: output_slice(:,:)

    ierr = 0

    input_slice = 0.d0
    input_slice(5, 4) = 1.d0

    sq = 1.d0
    sel = 0.8d0
    ser = 1.6d0
    sh = 0.5d0

    call convolve_slice_with_local_psf(input_slice, sq, sel, ser, sh, 2, 3, output_slice, ierr)
    if (ierr /= 0) return

    if (.not. ieee_is_finite(sum(output_slice))) then
       ierr = 11
       return
    end if
    if (abs(sum(output_slice) - 1.d0) > 1.d-12) then
       ierr = 12
       return
    end if
    if (output_slice(5, 4) <= 0.d0) then
       ierr = 13
       return
    end if

    sq = 1.d0
    sel = 1.d0
    ser = 1.d0
    sh = 0.d0
    input_slice = 0.d0
    input_slice(5, 4) = 1.d0
    call convolve_slice_with_local_psf(input_slice, sq, sel, ser, sh, 2, 2, output_slice, ierr)
    if (ierr /= 0) return
    if (.not. all(output_slice >= 0.d0)) then
      ierr = 14
      return
    end if

    input_slice = 0.d0
    input_slice(5, 4) = 1.d0
    sq = 1.d0
    sel = 0.8d0
    ser = 0.8d0
    sh = 0.d0
    sel(:, 1:nq/2) = 0.8d0
    ser(:, 1:nq/2) = 0.8d0
    sel(:, nq/2+1:nq) = 2.0d0
    ser(:, nq/2+1:nq) = 2.0d0
    call convolve_slice_with_local_psf(input_slice, sq, sel, ser, sh, 2, 3, output_slice, ierr)
    if (ierr /= 0) return

    ie_peak_left = maxloc(output_slice(:, 3), dim=1)
    ie_peak_right = maxloc(output_slice(:, 5), dim=1)
    if (output_slice(ie_peak_right, 5) >= output_slice(ie_peak_left, 3)) then
       ierr = 15
       return
    end if

  end subroutine selfcheck_psf_2d_convolution

end module psf_2d_convolution
