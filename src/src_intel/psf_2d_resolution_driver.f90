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

! Driver utilities for the position-dependent 2D PSF resolution path.
module psf_2d_resolution_driver

  use iso_fortran_env, only: error_unit
  use tas_resolution_cn, only: tas_cn_resolution_hkle
  use psf_2d_projection, only: build_raw_2d_psf_from_rm
  use psf_2d_param_estimation, only: estimate_psf_params_from_2d_kernel
  use psf_2d_convolution, only: convolve_slice_with_local_psf

  implicit none

contains

  ! Deposit one intrinsic mode onto the current energy grid before any resolution.
  subroutine deposit_intrinsic_mode(slice_line, ecenter, intensity, delta_e)

    implicit none

    real(kind=8), intent(inout) :: slice_line(:)
    real(kind=8), intent(in) :: ecenter, intensity, delta_e

    integer(kind=4) :: ne_local, idx0
    real(kind=8) :: pos, frac

    ne_local = size(slice_line)
    if (ne_local <= 0) return
    if (delta_e <= 0.d0) return
    if (ecenter < 0.d0) return

    pos = ecenter / delta_e + 1.d0
    idx0 = floor(pos)
    frac = pos - dble(idx0)

    if (idx0 < 1) then
       slice_line(1) = slice_line(1) + intensity
    elseif (idx0 >= ne_local) then
       slice_line(ne_local) = slice_line(ne_local) + intensity
    else
       slice_line(idx0) = slice_line(idx0) + intensity * (1.d0 - frac)
       slice_line(idx0 + 1) = slice_line(idx0 + 1) + intensity * frac
    end if

  end subroutine deposit_intrinsic_mode

  ! Build a PSF2D-observed slice from an intrinsic slice and local TAS CN RM(4,4).
  subroutine build_psf2d_observed_slice(qlist, elist, ee, intrinsic_temp, &
                                        alat, blat, clat, alpha_deg, beta_deg, gamma_deg, &
                                        orient1, orient2, psf_mode, &
                                        psf_sigma_q_default, psf_sigma_e_left_default, &
                                        psf_sigma_e_right_default, psf_shear_default, &
                                        psf_window_q, psf_window_e, &
                                        psf_temp, sigma_q_temp, sigma_e_left_temp, sigma_e_right_temp, shear_temp, &
                                        rel_l2_psf, print_intrinsic_built, print_rm_built, print_params_estimated, print_conv_applied, ierr)

    implicit none

    real(kind=8), intent(in) :: qlist(:,:), ee(:), intrinsic_temp(:,:)
    integer(kind=4), intent(in) :: elist(:)
    real(kind=8), intent(in) :: alat, blat, clat, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    character(len=*), intent(in) :: psf_mode
    real(kind=8), intent(in) :: psf_sigma_q_default, psf_sigma_e_left_default, psf_sigma_e_right_default, psf_shear_default
    integer(kind=4), intent(in) :: psf_window_q, psf_window_e
    real(kind=8), intent(out) :: psf_temp(:,:), sigma_q_temp(:,:), sigma_e_left_temp(:,:), sigma_e_right_temp(:,:), shear_temp(:,:)
    real(kind=8), intent(out) :: rel_l2_psf
    logical, intent(inout) :: print_intrinsic_built, print_rm_built, print_params_estimated, print_conv_applied
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: ne_local, nqh_local, ih, kk, ierr_psf
    real(kind=8) :: R0cn, RMcn(4,4), sigma_q_psf, sigma_el_psf, sigma_er_psf, shear_psf
    real(kind=8) :: intrinsic_norm, diff_norm
    real(kind=8), allocatable :: q_grid(:), e_grid(:), raw_kernel(:,:), psf_slice_local(:,:)
    logical :: used_fallback
    real(kind=8) :: dq_step, de_step

    ierr = 0
    rel_l2_psf = 0.d0
    ne_local = size(intrinsic_temp, 1)
    nqh_local = size(intrinsic_temp, 2)

    if (size(psf_temp, 1) /= ne_local .or. size(psf_temp, 2) /= nqh_local) then
       ierr = 1
       return
    end if
    if (size(sigma_q_temp, 1) /= ne_local .or. size(sigma_q_temp, 2) /= nqh_local) then
       ierr = 2
       return
    end if
    if (trim(psf_mode) /= 'analytic') then
       ierr = 3
       return
    end if

    if (.not. print_intrinsic_built) then
      write(*,*) '2D PSF resolution: intrinsic_slice built from pre-resolution mode list.'
      write(*,*)
      print_intrinsic_built = .TRUE.
    end if

    de_step = 1.d0
    if (size(ee) >= 2) de_step = max(abs(ee(2) - ee(1)), 1.d-12)

    do ih = 1, nqh_local
       dq_step = local_q_step(qlist, elist, ih)
       allocate(q_grid(2 * psf_window_q + 1), e_grid(2 * psf_window_e + 1))
       q_grid = dq_step * (/(dble(kk), kk=-psf_window_q, psf_window_q)/)
       e_grid = de_step * (/(dble(kk), kk=-psf_window_e, psf_window_e)/)
       do kk = 1, ne_local
          call tas_cn_resolution_hkle(qlist(:,elist(ih)), ee(kk), alat, blat, clat, alpha_deg, beta_deg, gamma_deg, &
                                      orient1, orient2, R0cn, RMcn, ierr_psf)
          if (ierr_psf == 0 .and. .not. print_rm_built) then
             write(*,*) '2D PSF resolution: local 4D RM/RMM built from TAS CN resolution at nominal points.'
             write(*,*)
             print_rm_built = .TRUE.
          end if
          if (ierr_psf == 0) then
             call build_raw_2d_psf_from_rm(RMcn, q_grid, e_grid, raw_kernel, ierr_psf)
          end if
          if (ierr_psf == 0) then
             call estimate_psf_params_from_2d_kernel(q_grid, e_grid, raw_kernel, &
                                                    psf_sigma_q_default * dq_step, psf_sigma_e_left_default * de_step, &
                                                    psf_sigma_e_right_default * de_step, psf_shear_default * dq_step / de_step, &
                                                    sigma_q_psf, sigma_el_psf, sigma_er_psf, shear_psf, &
                                                    used_fallback, ierr_psf)
             if (.not. print_params_estimated) then
                write(*,*) '2D PSF resolution: local PSF parameters estimated from 4D RM/RMM-derived samples.'
                write(*,*)
                print_params_estimated = .TRUE.
             end if
          end if
          if (ierr_psf /= 0) then
             sigma_q_psf = psf_sigma_q_default
             sigma_el_psf = psf_sigma_e_left_default
             sigma_er_psf = psf_sigma_e_right_default
             shear_psf = psf_shear_default
          else
             sigma_q_psf = sigma_q_psf / dq_step
             sigma_el_psf = sigma_el_psf / de_step
             sigma_er_psf = sigma_er_psf / de_step
             shear_psf = shear_psf * dq_step / de_step
          end if
          sigma_q_temp(kk, ih) = sigma_q_psf
          sigma_e_left_temp(kk, ih) = sigma_el_psf
          sigma_e_right_temp(kk, ih) = sigma_er_psf
          shear_temp(kk, ih) = shear_psf
          if (allocated(raw_kernel)) deallocate(raw_kernel)
       end do
       deallocate(q_grid, e_grid)
    end do

    call convolve_slice_with_local_psf(intrinsic_temp, sigma_q_temp, sigma_e_left_temp, sigma_e_right_temp, shear_temp, &
                                       psf_window_q, psf_window_e, psf_slice_local, ierr_psf)
    if (ierr_psf /= 0) then
       write(error_unit, *) 'Failed to build 2D PSF-convolved slice.'
       ierr = 20 + ierr_psf
       return
    end if
    psf_temp = psf_slice_local
    intrinsic_norm = sqrt(sum(intrinsic_temp * intrinsic_temp))
    diff_norm = sqrt(sum((psf_temp - intrinsic_temp) * (psf_temp - intrinsic_temp)))
    rel_l2_psf = diff_norm / max(intrinsic_norm, 1.d-30)

    if (.not. print_conv_applied) then
       write(*,*) '2D PSF resolution: local 2D convolution applied to intrinsic_slice.'
       write(*,*) '2D PSF resolution: rel_l2(psf2d,intrinsic)=', rel_l2_psf
       write(*,*)
       print_conv_applied = .TRUE.
    end if

    deallocate(psf_slice_local)

  end subroutine build_psf2d_observed_slice

  real(kind=8) function local_q_step(qlist, elist, ih)

    implicit none

    real(kind=8), intent(in) :: qlist(:,:)
    integer(kind=4), intent(in) :: elist(:), ih
    real(kind=8) :: dqv(3)
    integer(kind=4) :: ih_prev, ih_next

    ih_prev = max(1, ih - 1)
    ih_next = min(size(elist), ih + 1)
    if (ih_next == ih_prev) then
       local_q_step = 1.d0
       return
    end if
    dqv = qlist(:, elist(ih_next)) - qlist(:, elist(ih_prev))
    local_q_step = 0.5d0 * sqrt(dot_product(dqv, dqv))
    if (local_q_step <= 0.d0) local_q_step = 1.d0

  end function local_q_step

end module psf_2d_resolution_driver
