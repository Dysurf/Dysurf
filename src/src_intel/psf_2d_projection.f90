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

! Local 4D -> 2D sample projection utilities for future PSF work.
! This module does not build kernels and does not modify the main flow.
module psf_2d_projection

  use tas_resolution_cn, only: invert4

  implicit none

contains

  ! Project local 4D response samples (dQx,dQy,dQz,dE) to a target 2D slice (dq,dE).
  !
  ! Input:
  !   q_offsets(3, nsamp) : local Q offsets in the same sample-frame basis used by RM
  !   e_offsets(nsamp)    : local energy offsets corresponding to q_offsets
  !   qdir(3)             : target slice q direction in the same basis as q_offsets
  !
  ! Output:
  !   dq(nsamp)           : projected offsets along the target q direction
  !   dE(nsamp)           : energy offsets copied from e_offsets
  !   ierr                : 0 on success
  subroutine project_local_4d_response_to_2d(q_offsets, e_offsets, qdir, dq, dE, ierr)

    implicit none

    real(kind=8), intent(in) :: q_offsets(:,:), e_offsets(:), qdir(3)
    real(kind=8), allocatable, intent(out) :: dq(:), dE(:)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: nsamp
    real(kind=8) :: qhat(3), normq

    ierr = 0
    if (size(q_offsets, 1) /= 3) then
       ierr = 1
       return
    end if

    nsamp = size(q_offsets, 2)
    if (size(e_offsets) /= nsamp) then
      ierr = 2
      return
    end if

    normq = sqrt(dot_product(qdir, qdir))
    if (normq <= 0.d0) then
       ierr = 3
       return
    end if
    qhat = qdir / normq

    allocate(dq(nsamp), dE(nsamp))
    dq = matmul(transpose(q_offsets), qhat)
    dE = e_offsets

  end subroutine project_local_4d_response_to_2d

  ! Minimal adapter/stub: convert an analytic local 4D Gaussian precision matrix RM
  ! into a deterministic set of local 4D samples. The returned samples are the center
  ! point plus +/- nsigma times the Cholesky columns of the covariance matrix inv(RM).
  !
  ! Input:
  !   RM(4,4)             : local 4D precision matrix in sample-frame coordinates
  !   nsigma              : distance from center in covariance-scaled coordinates
  !
  ! Output:
  !   q_offsets(3, 9)     : deterministic local Q offsets
  !   e_offsets(9)        : deterministic local energy offsets
  !   ierr                : 0 on success
  subroutine adapt_rm_to_local_4d_samples(RM, nsigma, q_offsets, e_offsets, ierr)

    implicit none

    real(kind=8), intent(in) :: RM(4,4)
    real(kind=8), intent(in) :: nsigma
    real(kind=8), allocatable, intent(out) :: q_offsets(:,:), e_offsets(:)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: info, iax, isamp
    real(kind=8) :: cov4(4,4), detv, L(4,4), vec4(4)

    ierr = 0
    if (nsigma <= 0.d0) then
       ierr = 1
       return
    end if

    call invert4(RM, cov4, detv, info)
    if (info /= 0 .or. detv <= 0.d0) then
       ierr = 2
       return
    end if

    call cholesky4(cov4, L, info)
    if (info /= 0) then
       ierr = 3
       return
    end if

    allocate(q_offsets(3, 9), e_offsets(9))
    q_offsets = 0.d0
    e_offsets = 0.d0

    isamp = 1
    q_offsets(:, isamp) = 0.d0
    e_offsets(isamp) = 0.d0

    do iax = 1, 4
       vec4 = nsigma * L(:, iax)
       isamp = isamp + 1
       q_offsets(:, isamp) = vec4(1:3)
       e_offsets(isamp) = vec4(4)
       isamp = isamp + 1
       q_offsets(:, isamp) = -vec4(1:3)
       e_offsets(isamp) = -vec4(4)
    end do

  end subroutine adapt_rm_to_local_4d_samples

  ! Build a raw local 2D PSF target on a (dq,dE) grid directly from RM(4,4).
  ! The returned kernel is the projected/marginalized 2D Gaussian in the chosen
  ! q direction and the energy axis.
  !
  ! Input:
  !   RM(4,4)             : local 4D precision matrix in sample-frame coordinates
  !   q_grid(:)           : target dq grid, in the same q unit used to define qdir
  !   e_grid(:)           : target dE grid, in the same energy unit as RM
  !
  ! Output:
  !   kernel(:,:)         : normalized raw 2D PSF target on (e_grid,q_grid)
  !   ierr                : 0 on success
  subroutine build_raw_2d_psf_from_rm(RM, q_grid, e_grid, kernel, ierr)

    implicit none

    real(kind=8), intent(in) :: RM(4,4), q_grid(:), e_grid(:)
    real(kind=8), allocatable, intent(out) :: kernel(:,:)
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: A(2,2), B(2,2), C(2,2), Cinv(2,2), P2(2,2), cov2(2,2), detc, detp
    real(kind=8) :: vec(2), expo
    integer(kind=4) :: iq, ie, info

    ierr = 0
    if (size(q_grid) <= 0 .or. size(e_grid) <= 0) then
       ierr = 1
       return
    end if

    A = reshape((/RM(1,1), RM(1,4), RM(4,1), RM(4,4)/), (/2,2/))
    B = reshape((/RM(1,2), RM(1,3), RM(4,2), RM(4,3)/), (/2,2/))
    C = reshape((/RM(2,2), RM(2,3), RM(3,2), RM(3,3)/), (/2,2/))

    call invert2x2(C, Cinv, detc, info)
    if (info /= 0 .or. detc <= 0.d0) then
       ierr = 2
       return
    end if

    P2 = A - matmul(B, matmul(Cinv, transpose(B)))
    call invert2x2(P2, cov2, detp, info)
    if (info /= 0 .or. detp <= 0.d0) then
       ierr = 3
       return
    end if

    allocate(kernel(size(e_grid), size(q_grid)))
    do iq = 1, size(q_grid)
       do ie = 1, size(e_grid)
          vec = (/q_grid(iq), e_grid(ie)/)
          expo = dot_product(vec, matmul(P2, vec))
          kernel(ie, iq) = exp(-0.5d0 * expo)
       end do
    end do
    kernel = kernel / max(sum(kernel), 1.d-30)

  end subroutine build_raw_2d_psf_from_rm

  ! Minimal selfcheck:
  ! 1. project a hand-built sample and verify dq
  ! 2. verify dE is passed through unchanged
  subroutine selfcheck_project_local_4d_response_to_2d(ierr)

    implicit none

    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: q_offsets_in(3,2), e_offsets_in(2), qdir(3)
    real(kind=8), allocatable :: dq(:), dE(:)

    ierr = 0
    q_offsets_in(:,1) = (/1.d0, 0.d0, 0.d0/)
    q_offsets_in(:,2) = (/1.d0, 1.d0, 0.d0/)
    e_offsets_in = (/0.25d0, -0.50d0/)
    qdir = (/2.d0, 0.d0, 0.d0/)

    call project_local_4d_response_to_2d(q_offsets_in, e_offsets_in, qdir, dq, dE, ierr)
    if (ierr /= 0) return

    if (abs(dq(1) - 1.d0) > 1.d-12) then
       ierr = 11
       return
    end if
    if (abs(dq(2) - 1.d0) > 1.d-12) then
       ierr = 12
       return
    end if
    if (abs(dE(1) - 0.25d0) > 1.d-12) then
       ierr = 13
       return
    end if
    if (abs(dE(2) + 0.50d0) > 1.d-12) then
       ierr = 14
       return
    end if

    deallocate(dq, dE)

  end subroutine selfcheck_project_local_4d_response_to_2d

  subroutine cholesky4(A, L, ierr)

    implicit none

    real(kind=8), intent(in) :: A(4,4)
    real(kind=8), intent(out) :: L(4,4)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: i, j, k
    real(kind=8) :: sumv

    ierr = 0
    L = 0.d0

    do i = 1, 4
       do j = 1, i
          sumv = A(i, j)
          do k = 1, j - 1
             sumv = sumv - L(i, k) * L(j, k)
          end do
          if (i == j) then
             if (sumv <= 0.d0) then
                ierr = 1
                L = 0.d0
                return
             end if
             L(i, j) = sqrt(sumv)
          else
             if (abs(L(j, j)) <= 1.d-20) then
                ierr = 2
                L = 0.d0
                return
             end if
             L(i, j) = sumv / L(j, j)
          end if
       end do
    end do

  end subroutine cholesky4

  subroutine invert2x2(A, Ainv, detv, ierr)

    implicit none

    real(kind=8), intent(in) :: A(2,2)
    real(kind=8), intent(out) :: Ainv(2,2), detv
    integer(kind=4), intent(out) :: ierr

    ierr = 0
    detv = A(1,1) * A(2,2) - A(1,2) * A(2,1)
    if (abs(detv) <= 1.d-30) then
       ierr = 1
       Ainv = 0.d0
       return
    end if
    Ainv(1,1) = A(2,2) / detv
    Ainv(1,2) = -A(1,2) / detv
    Ainv(2,1) = -A(2,1) / detv
    Ainv(2,2) = A(1,1) / detv

  end subroutine invert2x2

end module psf_2d_projection
