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
! Local fixed-grid 4D resolution sampling utilities.
! This module provides internal sampled-Q helpers derived from the local

! RM(4,4) description. It is not a user-facing observed-spectrum mode.
module tas_resolution_sampling

  use material_sqw, only: material_mode_info_at_q
  use variables, only: nbands, crlatvec

  implicit none

  integer(kind=4), parameter :: tas_sampling_accuracy_default = 5

  type :: tas_resolution_sample_cache
     integer(kind=4) :: nsamp = 0
     real(kind=8) :: r0 = 0.d0
     real(kind=8) :: gamma_factor = 0.d0
     real(kind=8) :: omega_factor_x = 0.d0
     real(kind=8) :: omega_factor_y = 0.d0
     real(kind=8) :: prefactor = 0.d0
     real(kind=8), allocatable :: dq1(:)
     real(kind=8), allocatable :: dq2(:)
     real(kind=8), allocatable :: weight(:)
     real(kind=8), allocatable :: disp(:)
     real(kind=8), allocatable :: inte(:)
  end type tas_resolution_sample_cache

contains

  subroutine tas_resolution_sampling_prepare(q0_hkl, xvec, yvec, zvec, R0, RM, branch_id, rmsd_local, temp_local, cache, ierr)

    implicit none

    real(kind=8), intent(in) :: q0_hkl(3), xvec(3), yvec(3), zvec(3), R0, RM(4,4)
    integer(kind=4), intent(in) :: branch_id
    real(kind=8), intent(in) :: rmsd_local(:,:), temp_local
    type(tas_resolution_sample_cache), intent(out) :: cache
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: ix, iy, iz, isamp, nxy, nz, accuracy, info
    real(kind=8) :: Mxx, Mxy, Mxw, Myy, Myw, Mzz, Mww
    real(kind=8) :: detxy, detz
    real(kind=8) :: tqxx, tqxy, tqy, tqz
    real(kind=8) :: step1, step2, phi_x, phi_y, phi_z
    real(kind=8) :: tx, ty, tz, dQ1, dQ2, dQ3
    real(kind=8) :: norm, normz, h1, k1, l1, qcart(3)
    real(kind=8) :: omega_modes(nbands), intensities(nbands)
    logical :: mode_mask(nbands)

    ierr = 0
    cache%nsamp = 0
    cache%r0 = 0.d0
    cache%gamma_factor = 0.d0
    cache%omega_factor_x = 0.d0
    cache%omega_factor_y = 0.d0
    cache%prefactor = 0.d0
    if (allocated(cache%dq1)) deallocate(cache%dq1, cache%dq2, cache%weight, cache%disp, cache%inte)

    accuracy = tas_sampling_accuracy_default
    mode_mask = .false.
    if (branch_id < 1 .or. branch_id > nbands) then
       ierr = 10
       return
    end if
    mode_mask(branch_id) = .true.

    Mww = RM(4, 4)
    Mxw = RM(1, 4)
    Myw = RM(2, 4)
    Mzz = RM(3, 3)
    Mxx = RM(1, 1) - Mxw * Mxw / Mww
    Myy = RM(2, 2) - Myw * Myw / Mww
    Mxy = RM(1, 2) - Mxw * Myw / Mww

    if (Mww <= 0.d0 .or. Myy <= 0.d0 .or. Mzz <= 0.d0) then
       ierr = 11
       return
    end if

    detxy = sqrt(Mxx * Myy - Mxy * Mxy)
    detz = sqrt(Mzz)
    if (Mxx <= 0.d0 .or. detxy <= 0.d0 .or. detz <= 0.d0) then
       ierr = 12
       return
    end if

    cache%gamma_factor = sqrt(Mww / 2.d0)
    cache%omega_factor_x = Mxw / sqrt(2.d0 * Mww)
    cache%omega_factor_y = Myw / sqrt(2.d0 * Mww)
    cache%r0 = R0

    tqz = 1.d0 / detz
    tqy = sqrt(Mxx) / detxy
    tqxx = 1.d0 / sqrt(Mxx)
    tqxy = Mxy / sqrt(Mxx) / detxy

    nxy = 2 * accuracy + 1
    nz = 2 * accuracy + 1
    step1 = acos(-1.d0) / dble(nxy)
    step2 = acos(-1.d0) / dble(nz)
    cache%prefactor = step1 ** 2 * step2
    if (accuracy == 0) cache%prefactor = cache%prefactor * 0.79788d0 ** 3

    cache%nsamp = nxy * nxy * nz
    allocate(cache%dq1(cache%nsamp), cache%dq2(cache%nsamp), cache%weight(cache%nsamp), &
             cache%disp(cache%nsamp), cache%inte(cache%nsamp))

    isamp = 0
    do iz = 1, nz
       phi_z = -0.5d0 * acos(-1.d0) + 0.5d0 * step2 + dble(iz - 1) * step2
       tz = tan(phi_z)
       normz = exp(-0.5d0 * tz * tz) * (1.d0 + tz * tz)
       do iy = 1, nxy
          phi_y = -0.5d0 * acos(-1.d0) + 0.5d0 * step1 + dble(iy - 1) * step1
          ty = tan(phi_y)
          do ix = 1, nxy
             phi_x = -0.5d0 * acos(-1.d0) + 0.5d0 * step1 + dble(ix - 1) * step1
             tx = tan(phi_x)
             norm = exp(-0.5d0 * (tx * tx + ty * ty)) * (1.d0 + tx * tx) * (1.d0 + ty * ty)

             dQ1 = tqxx * tx - tqxy * ty
             dQ2 = tqy * ty
             dQ3 = tqz * tz
             h1 = q0_hkl(1) + dQ1 * xvec(1) + dQ2 * yvec(1) + dQ3 * zvec(1)
             k1 = q0_hkl(2) + dQ1 * xvec(2) + dQ2 * yvec(2) + dQ3 * zvec(2)
             l1 = q0_hkl(3) + dQ1 * xvec(3) + dQ2 * yvec(3) + dQ3 * zvec(3)
             qcart = matmul(crlatvec, (/h1, k1, l1/))

             call material_mode_info_at_q(qcart, rmsd_local, temp_local, omega_modes, intensities, info, mode_mask)
             if (info /= 0) then
                ierr = 40 + info
                return
             end if

             isamp = isamp + 1
             cache%dq1(isamp) = dQ1
             cache%dq2(isamp) = dQ2
             cache%weight(isamp) = norm * normz / detxy / detz
             cache%disp(isamp) = omega_modes(branch_id)
             cache%inte(isamp) = max(intensities(branch_id), 0.d0)
          end do
       end do
    end do

  end subroutine tas_resolution_sampling_prepare

  subroutine tas_resolution_sampling_evaluate(cache, omega0, conv)

    implicit none

    type(tas_resolution_sample_cache), intent(in) :: cache
    real(kind=8), intent(in) :: omega0
    real(kind=8), intent(out) :: conv

    integer(kind=4) :: isamp
    real(kind=8) :: omega_red

    conv = 0.d0
    do isamp = 1, cache%nsamp
       omega_red = cache%gamma_factor * (cache%disp(isamp) - omega0) + &
                   cache%omega_factor_x * cache%dq1(isamp) + &
                   cache%omega_factor_y * cache%dq2(isamp)
       conv = conv + cache%inte(isamp) * sampled_resolution_kernel_zero(omega_red) * cache%weight(isamp)
    end do
    conv = conv * cache%prefactor * cache%r0

  end subroutine tas_resolution_sampling_evaluate

  subroutine tas_resolution_sampling_free(cache)

    implicit none

    type(tas_resolution_sample_cache), intent(inout) :: cache

    if (allocated(cache%dq1)) deallocate(cache%dq1, cache%dq2, cache%weight, cache%disp, cache%inte)
    cache%nsamp = 0

  end subroutine tas_resolution_sampling_free

  real(kind=8) function sampled_resolution_kernel_zero(omega)

    implicit none

    real(kind=8), intent(in) :: omega

    sampled_resolution_kernel_zero = exp(-omega * omega) / sqrt(acos(-1.d0))

  end function sampled_resolution_kernel_zero

end module tas_resolution_sampling
