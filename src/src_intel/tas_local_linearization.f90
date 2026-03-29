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

! Local linearization of phonon branches around a TAS observation point
module tas_local_linearization

  use tas_resolution_cn, only: tas_cn_standard_system, scalar_metric, invert3

  implicit none

contains

  ! Build a local linear model for one phonon branch around a grid point.
  ! The output gradient is expressed in the TAS sample coordinates (Qx,Qy,Qz).
  subroutine tas_branch_local_model(ih, ik, il, ib, qlist, omega_grid, nqh, nqk, nql, &
                                    a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, &
                                    e0, grad_hkl, grad_sample, ierr)

    implicit none

    integer(kind=4), intent(in) :: ih, ik, il, ib
    integer(kind=4), intent(in) :: nqh, nqk, nql
    real(kind=8), intent(in) :: qlist(:, :)
    real(kind=8), intent(in) :: omega_grid(:, :)
    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    real(kind=8), intent(out) :: e0
    real(kind=8), intent(out) :: grad_hkl(3)
    real(kind=8), intent(out) :: grad_sample(3)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: ic, info
    real(kind=8) :: x(3), y(3), z(3), rlattice(6)
    real(kind=8) :: dq_hkl(3, 3), dq_sample(3, 3), rhs(3), dqt(3, 3), dqt_inv(3, 3), detv

    ierr = 0
    e0 = 0.d0
    grad_hkl = 0.d0
    grad_sample = 0.d0

    if (size(qlist, 1) /= 3 .or. size(omega_grid, 1) /= size(qlist, 2)) then
       ierr = 5
       return
    end if

    if (ih < 1 .or. ih > nqh .or. ik < 1 .or. ik > nqk .or. il < 1 .or. il > nql) then
       ierr = 1
       return
    end if

    ic = grid_index(ih, ik, il, nqk, nql)
    if (ib < 1 .or. ib > size(omega_grid, 2)) then
       ierr = 2
       return
    end if

    e0 = omega_grid(ic, ib)

    call tas_cn_standard_system(a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, x, y, z, rlattice, ierr)
    if (ierr /= 0) then
       ierr = 100 + ierr
       return
    end if

    call finite_difference_axis(1, ih, ik, il, ib, qlist, omega_grid, nqh, nqk, nql, dq_hkl(:, 1), rhs(1), ierr)
    if (ierr /= 0) return
    call finite_difference_axis(2, ih, ik, il, ib, qlist, omega_grid, nqh, nqk, nql, dq_hkl(:, 2), rhs(2), ierr)
    if (ierr /= 0) return
    call finite_difference_axis(3, ih, ik, il, ib, qlist, omega_grid, nqh, nqk, nql, dq_hkl(:, 3), rhs(3), ierr)
    if (ierr /= 0) return

    dq_sample(1, :) = (/scalar_metric(x, dq_hkl(:, 1), rlattice), scalar_metric(x, dq_hkl(:, 2), rlattice), scalar_metric(x, dq_hkl(:, 3), rlattice)/)
    dq_sample(2, :) = (/scalar_metric(y, dq_hkl(:, 1), rlattice), scalar_metric(y, dq_hkl(:, 2), rlattice), scalar_metric(y, dq_hkl(:, 3), rlattice)/)
    dq_sample(3, :) = (/scalar_metric(z, dq_hkl(:, 1), rlattice), scalar_metric(z, dq_hkl(:, 2), rlattice), scalar_metric(z, dq_hkl(:, 3), rlattice)/)

    dqt = transpose(dq_hkl)
    call invert3(dqt, dqt_inv, detv, info)
    if (info /= 0 .or. abs(detv) <= 1.d-14) then
       ierr = 3
       return
    end if
    grad_hkl = matmul(dqt_inv, rhs)

    dqt = transpose(dq_sample)
    call invert3(dqt, dqt_inv, detv, info)
    if (info /= 0 .or. abs(detv) <= 1.d-14) then
       ierr = 4
       return
    end if
    grad_sample = matmul(dqt_inv, rhs)

  end subroutine tas_branch_local_model

  ! Evaluate the local linearized phonon energy at a nearby point in sample coordinates.
  subroutine tas_linearized_energy(e0, grad_sample, dq_sample, elin)

    implicit none

    real(kind=8), intent(in) :: e0, grad_sample(3), dq_sample(3)
    real(kind=8), intent(out) :: elin

    elin = e0 + dot_product(grad_sample, dq_sample)

  end subroutine tas_linearized_energy

  subroutine finite_difference_axis(axis_id, ih, ik, il, ib, qlist, omega_grid, nqh, nqk, nql, dq_vec, de, ierr)

    implicit none

    integer(kind=4), intent(in) :: axis_id, ih, ik, il, ib
    integer(kind=4), intent(in) :: nqh, nqk, nql
    real(kind=8), intent(in) :: qlist(:, :)
    real(kind=8), intent(in) :: omega_grid(:, :)
    real(kind=8), intent(out) :: dq_vec(3), de
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: i1, i2
    integer(kind=4) :: ih1, ih2, ik1, ik2, il1, il2

    ierr = 0
    dq_vec = 0.d0
    de = 0.d0

    ih1 = ih
    ih2 = ih
    ik1 = ik
    ik2 = ik
    il1 = il
    il2 = il

    if (axis_id == 1) then
       if (nqh < 2) then
          ierr = 11
          return
       end if
       if (ih == 1) then
          ih2 = ih + 1
       elseif (ih == nqh) then
          ih1 = ih - 1
       else
          ih1 = ih - 1
          ih2 = ih + 1
       end if
    elseif (axis_id == 2) then
       if (nqk < 2) then
          ierr = 12
          return
       end if
       if (ik == 1) then
          ik2 = ik + 1
       elseif (ik == nqk) then
          ik1 = ik - 1
       else
          ik1 = ik - 1
          ik2 = ik + 1
       end if
    elseif (axis_id == 3) then
       if (nql < 2) then
          ierr = 13
          return
       end if
       if (il == 1) then
          il2 = il + 1
       elseif (il == nql) then
          il1 = il - 1
       else
          il1 = il - 1
          il2 = il + 1
       end if
    else
       ierr = 14
       return
    end if

    i1 = grid_index(ih1, ik1, il1, nqk, nql)
    i2 = grid_index(ih2, ik2, il2, nqk, nql)

    dq_vec = qlist(:, i2) - qlist(:, i1)
    de = omega_grid(i2, ib) - omega_grid(i1, ib)

    if (sqrt(dot_product(dq_vec, dq_vec)) <= 1.d-14) then
       ierr = 15
    end if

  end subroutine finite_difference_axis

  integer(kind=4) function grid_index(ih, ik, il, nqk, nql)

    implicit none

    integer(kind=4), intent(in) :: ih, ik, il, nqk, nql

    grid_index = (ih - 1) * nqk * nql + (ik - 1) * nql + il

  end function grid_index

end module tas_local_linearization
