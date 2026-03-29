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
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Estimate analytic 2D PSF parameters from local 2D samples (dq, dE).
module psf_2d_param_estimation

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

contains

  ! Estimate analytic PSF parameters directly from a raw 2D PSF target using
  ! weighted moments on the (dq,dE) grid.
  !
  ! Input:
  !   dq_grid(:), dE_grid(:)          : target 2D grid coordinates
  !   raw_kernel(:,:)                 : normalized or non-normalized nonnegative 2D target
  !   sigma_q_default                 : fallback sigma_q
  !   sigma_e_left_default            : fallback left energy width
  !   sigma_e_right_default           : fallback right energy width
  !   shear_default                   : fallback shear
  !
  ! Output:
  !   sigma_q                         : weighted q width
  !   sigma_e_left                    : weighted left width after shear removal
  !   sigma_e_right                   : weighted right width after shear removal
  !   shear                           : weighted q-E tilt from covariance / var_q
  !   used_fallback                   : .TRUE. if defaults were used fully or partially
  !   ierr                            : 0 on success, nonzero only for gross input issues
  subroutine estimate_psf_params_from_2d_kernel(dq_grid, dE_grid, raw_kernel, sigma_q_default, sigma_e_left_default, &
                                                sigma_e_right_default, shear_default, sigma_q, sigma_e_left, sigma_e_right, &
                                                shear, used_fallback, ierr)

    implicit none

    real(kind=8), intent(in) :: dq_grid(:), dE_grid(:), raw_kernel(:,:)
    real(kind=8), intent(in) :: sigma_q_default, sigma_e_left_default, sigma_e_right_default, shear_default
    real(kind=8), intent(out) :: sigma_q, sigma_e_left, sigma_e_right, shear
    logical, intent(out) :: used_fallback
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: iq, ie, nleft, nright
    real(kind=8) :: wsum, mean_q, mean_e, var_q, cov_qe, sumsq_left, sumsq_right
    real(kind=8) :: qv, ev, w, de_eff
    real(kind=8), parameter :: tiny_var = 1.d-20

    ierr = 0
    used_fallback = .FALSE.
    sigma_q = sigma_q_default
    sigma_e_left = sigma_e_left_default
    sigma_e_right = sigma_e_right_default
    shear = shear_default

    if (size(raw_kernel,1) /= size(dE_grid) .or. size(raw_kernel,2) /= size(dq_grid)) then
       ierr = 1
       used_fallback = .TRUE.
       return
    end if

    wsum = sum(raw_kernel)
    if (.not. ieee_is_finite(wsum) .or. wsum <= 0.d0) then
       ierr = 2
       used_fallback = .TRUE.
       return
    end if

    mean_q = 0.d0
    mean_e = 0.d0
    do iq = 1, size(dq_grid)
       do ie = 1, size(dE_grid)
          w = raw_kernel(ie, iq)
          if (.not. ieee_is_finite(w) .or. w <= 0.d0) cycle
          mean_q = mean_q + w * dq_grid(iq)
          mean_e = mean_e + w * dE_grid(ie)
       end do
    end do
    mean_q = mean_q / wsum
    mean_e = mean_e / wsum

    var_q = 0.d0
    cov_qe = 0.d0
    do iq = 1, size(dq_grid)
       do ie = 1, size(dE_grid)
          w = raw_kernel(ie, iq)
          if (.not. ieee_is_finite(w) .or. w <= 0.d0) cycle
          qv = dq_grid(iq) - mean_q
          ev = dE_grid(ie) - mean_e
          var_q = var_q + w * qv * qv
          cov_qe = cov_qe + w * qv * ev
       end do
    end do
    var_q = var_q / wsum
    cov_qe = cov_qe / wsum

    if (.not. ieee_is_finite(var_q) .or. var_q <= tiny_var) then
       used_fallback = .TRUE.
       return
    end if

    sigma_q = sqrt(max(var_q, 0.d0))
    if (.not. ieee_is_finite(sigma_q) .or. sigma_q <= 0.d0) then
       sigma_q = sigma_q_default
       used_fallback = .TRUE.
    end if

    shear = cov_qe / var_q
    if (.not. ieee_is_finite(shear)) then
       shear = shear_default
       used_fallback = .TRUE.
    end if

    nleft = 0
    nright = 0
    sumsq_left = 0.d0
    sumsq_right = 0.d0
    do iq = 1, size(dq_grid)
       do ie = 1, size(dE_grid)
          w = raw_kernel(ie, iq)
          if (.not. ieee_is_finite(w) .or. w <= 0.d0) cycle
          de_eff = (dE_grid(ie) - mean_e) - shear * (dq_grid(iq) - mean_q)
          if (de_eff < 0.d0) then
             nleft = nleft + 1
             sumsq_left = sumsq_left + w * de_eff * de_eff
          elseif (de_eff > 0.d0) then
             nright = nright + 1
             sumsq_right = sumsq_right + w * de_eff * de_eff
          end if
       end do
    end do

    if (nleft > 0) then
       sigma_e_left = sqrt(sumsq_left / max(sum(raw_kernel, mask=((spread(dE_grid,2,size(dq_grid)) - mean_e) - &
                             shear * transpose(spread(dq_grid,2,size(dE_grid))) < 0.d0)), 1.d-30))
       if (.not. ieee_is_finite(sigma_e_left) .or. sigma_e_left <= 0.d0) then
          sigma_e_left = sigma_e_left_default
          used_fallback = .TRUE.
       end if
    else
       sigma_e_left = sigma_e_left_default
       used_fallback = .TRUE.
    end if

    if (nright > 0) then
       sigma_e_right = sqrt(sumsq_right / max(sum(raw_kernel, mask=((spread(dE_grid,2,size(dq_grid)) - mean_e) - &
                              shear * transpose(spread(dq_grid,2,size(dE_grid))) > 0.d0)), 1.d-30))
       if (.not. ieee_is_finite(sigma_e_right) .or. sigma_e_right <= 0.d0) then
          sigma_e_right = sigma_e_right_default
          used_fallback = .TRUE.
       end if
    else
       sigma_e_right = sigma_e_right_default
       used_fallback = .TRUE.
    end if

  end subroutine estimate_psf_params_from_2d_kernel

  ! Estimate analytic PSF parameters from local 2D samples.
  !
  ! Input:
  !   dq(:), dE(:)              : local 2D samples on the target slice
  !   sigma_q_default           : fallback sigma_q
  !   sigma_e_left_default      : fallback left energy width
  !   sigma_e_right_default     : fallback right energy width
  !   shear_default             : fallback shear
  !
  ! Output:
  !   sigma_q                   : estimated q-direction Gaussian width
  !   sigma_e_left              : estimated left-side dE_eff width
  !   sigma_e_right             : estimated right-side dE_eff width
  !   shear                     : estimated dE = shear * dq + residual slope
  !   used_fallback             : .TRUE. if defaults were used fully or partially
  !   ierr                      : 0 on success, nonzero only for gross input issues
  subroutine estimate_psf_params_from_4d_sample(dq, dE, sigma_q_default, sigma_e_left_default, &
                                                sigma_e_right_default, shear_default, &
                                                sigma_q, sigma_e_left, sigma_e_right, shear, &
                                                used_fallback, ierr)

    implicit none

    real(kind=8), intent(in) :: dq(:), dE(:)
    real(kind=8), intent(in) :: sigma_q_default, sigma_e_left_default, sigma_e_right_default, shear_default
    real(kind=8), intent(out) :: sigma_q, sigma_e_left, sigma_e_right, shear
    logical, intent(out) :: used_fallback
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: n, i, nvalid, nleft, nright
    real(kind=8), allocatable :: qv(:), ev(:), de_eff(:)
    real(kind=8) :: mean_q, mean_e, var_q, cov_qe, sumsq_left, sumsq_right
    real(kind=8), parameter :: tiny_var = 1.d-20

    ierr = 0
    used_fallback = .FALSE.
    sigma_q = sigma_q_default
    sigma_e_left = sigma_e_left_default
    sigma_e_right = sigma_e_right_default
    shear = shear_default

    n = size(dq)
    if (size(dE) /= n) then
       ierr = 1
       used_fallback = .TRUE.
       return
    end if
    if (n <= 0) then
       ierr = 2
       used_fallback = .TRUE.
       return
    end if

    allocate(qv(n), ev(n))
    nvalid = 0
    do i = 1, n
       if (ieee_is_finite(dq(i)) .and. ieee_is_finite(dE(i))) then
          nvalid = nvalid + 1
          qv(nvalid) = dq(i)
          ev(nvalid) = dE(i)
       end if
    end do

    if (nvalid < 4) then
       used_fallback = .TRUE.
       deallocate(qv, ev)
       return
    end if

    mean_q = sum(qv(1:nvalid)) / dble(nvalid)
    mean_e = sum(ev(1:nvalid)) / dble(nvalid)
    var_q = sum((qv(1:nvalid) - mean_q) ** 2) / dble(nvalid)

    if (.not. ieee_is_finite(var_q) .or. var_q <= tiny_var) then
       used_fallback = .TRUE.
       deallocate(qv, ev)
       return
    end if

    sigma_q = sqrt(max(var_q, 0.d0))
    if (.not. ieee_is_finite(sigma_q) .or. sigma_q <= 0.d0) then
       sigma_q = sigma_q_default
       used_fallback = .TRUE.
    end if

    cov_qe = sum((qv(1:nvalid) - mean_q) * (ev(1:nvalid) - mean_e)) / dble(nvalid)
    if (ieee_is_finite(cov_qe)) then
       shear = cov_qe / var_q
       if (.not. ieee_is_finite(shear)) then
          shear = shear_default
          used_fallback = .TRUE.
       end if
    else
       shear = shear_default
       used_fallback = .TRUE.
    end if

    allocate(de_eff(nvalid))
    de_eff = (ev(1:nvalid) - mean_e) - shear * (qv(1:nvalid) - mean_q)

    nleft = 0
    nright = 0
    sumsq_left = 0.d0
    sumsq_right = 0.d0
    do i = 1, nvalid
       if (.not. ieee_is_finite(de_eff(i))) cycle
       if (de_eff(i) < 0.d0) then
          nleft = nleft + 1
          sumsq_left = sumsq_left + de_eff(i) * de_eff(i)
       elseif (de_eff(i) > 0.d0) then
          nright = nright + 1
          sumsq_right = sumsq_right + de_eff(i) * de_eff(i)
       end if
    end do

    if (nleft >= 2) then
       sigma_e_left = sqrt(sumsq_left / dble(nleft))
       if (.not. ieee_is_finite(sigma_e_left) .or. sigma_e_left <= 0.d0) then
          sigma_e_left = sigma_e_left_default
          used_fallback = .TRUE.
       end if
    else
       sigma_e_left = sigma_e_left_default
       used_fallback = .TRUE.
    end if

    if (nright >= 2) then
       sigma_e_right = sqrt(sumsq_right / dble(nright))
       if (.not. ieee_is_finite(sigma_e_right) .or. sigma_e_right <= 0.d0) then
          sigma_e_right = sigma_e_right_default
          used_fallback = .TRUE.
       end if
    else
       sigma_e_right = sigma_e_right_default
       used_fallback = .TRUE.
    end if

    deallocate(qv, ev, de_eff)

  end subroutine estimate_psf_params_from_4d_sample

  ! Minimal selfcheck:
  ! construct deterministic samples with positive shear and asymmetric left/right widths.
  ! Check only finiteness and trend correctness, not exact recovery.
  subroutine selfcheck_estimate_psf_params_from_4d_sample(ierr)

    implicit none

    integer(kind=4), intent(out) :: ierr

    integer(kind=4), parameter :: n = 14
    integer(kind=4) :: i
    real(kind=8) :: dq(n), dE(n)
    real(kind=8) :: sigma_q, sigma_e_left, sigma_e_right, shear
    logical :: used_fallback

    ierr = 0

    do i = 1, 7
       dq(i) = -0.6d0 + 0.1d0 * dble(i - 1)
       dE(i) = 0.4d0 * dq(i) - 0.30d0 - 0.02d0 * dble(i - 1)
    end do
    do i = 8, 14
       dq(i) = -0.1d0 + 0.1d0 * dble(i - 8)
       dE(i) = 0.4d0 * dq(i) + 0.70d0 + 0.08d0 * dble(i - 8)
    end do

    call estimate_psf_params_from_4d_sample(dq, dE, 1.d0, 1.d0, 1.d0, 0.d0, &
                                            sigma_q, sigma_e_left, sigma_e_right, shear, &
                                            used_fallback, ierr)
    if (ierr /= 0) return

    if (.not. ieee_is_finite(sigma_q) .or. sigma_q <= 0.d0) then
       ierr = 11
       return
    end if
    if (.not. ieee_is_finite(shear)) then
       ierr = 12
       return
    end if
    if (.not. ieee_is_finite(sigma_e_left) .or. .not. ieee_is_finite(sigma_e_right)) then
       ierr = 13
       return
    end if
    if (shear <= 0.d0) then
       ierr = 14
       return
    end if
    if (sigma_e_right <= sigma_e_left) then
       ierr = 15
       return
    end if

  end subroutine selfcheck_estimate_psf_params_from_4d_sample

end module psf_2d_param_estimation
