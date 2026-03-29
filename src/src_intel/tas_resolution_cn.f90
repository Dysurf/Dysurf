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

! Cooper-Nathans TAS resolution interface
module tas_resolution_cn

  use variables, only: tas_fix_mode, tas_e_fixed, tas_dm, tas_da, &
                       tas_coll_h_pre_mono, tas_coll_h_pre_samp, &
                       tas_coll_h_post_samp, tas_coll_h_post_ana, &
                       tas_mosaic_mono_h, tas_mosaic_ana_h, tas_mosaic_samp_h, &
                       tas_use_mosaic, tas_mono_dir, tas_sample_dir, tas_ana_dir

  implicit none

  real(kind=8), parameter :: mev_to_k2 = 2.072d0
  real(kind=8), parameter :: pi = acos(-1.d0)
  real(kind=8), parameter :: fwhm_to_sigma = 1.d0/2.3548200450309493d0
  real(kind=8), parameter :: default_vcol_arcmin = 120.d0
  real(kind=8), parameter :: tiny_mosaic_arcmin = 1.d-6

contains

  ! unified cartesian TAS resolution output for future 4D convolution
  subroutine tas_cn_resolution(qvec_cart, omega, R0, RM, ierr)

    implicit none

    real(kind=8), intent(in) :: qvec_cart(:)
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(out) :: R0
    real(kind=8), intent(out) :: RM(4,4)
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: info
    real(kind=8) :: qmag
    real(kind=8) :: qhat_cart(3), that_cart(3), vhat_cart(3), ref_axis(3)
    real(kind=8) :: rm_q(4,4), basis(3,3), evec(3), rm_cart(4,4)

    ierr = 0
    RM = 0.d0
    R0 = 0.d0

    if (size(qvec_cart) /= 3) then
       ierr = 10
       return
    end if

    qmag = sqrt(dot_product(qvec_cart, qvec_cart))
    if (qmag <= 0.d0) then
       ierr = 11
       return
    end if

    call tas_cn_rm_qframe(qmag, omega, rm_q, R0, ierr)
    if (ierr /= 0) return

    qhat_cart = qvec_cart/qmag
    ref_axis = (/0.d0, 0.d0, 1.d0/)
    call cross3(ref_axis, qhat_cart, that_cart)
    if (sqrt(dot_product(that_cart, that_cart)) <= 1.d-12) then
       ref_axis = (/1.d0, 0.d0, 0.d0/)
       call cross3(ref_axis, qhat_cart, that_cart)
    end if
    call normalize3(that_cart, info)
    if (info /= 0) then
       ierr = 14
       return
    end if
    call cross3(qhat_cart, that_cart, vhat_cart)
    call normalize3(vhat_cart, info)
    if (info /= 0) then
       ierr = 15
       return
    end if

    basis(:,1) = qhat_cart
    basis(:,2) = that_cart
    basis(:,3) = vhat_cart

    rm_cart = 0.d0
    rm_cart(1:3,1:3) = matmul(basis, matmul(rm_q(1:3,1:3), transpose(basis)))
    evec = matmul(basis, rm_q(1:3,4))
    rm_cart(1:3,4) = evec
    rm_cart(4,1:3) = evec
    rm_cart(4,4) = rm_q(4,4)

    RM = 0.5d0*(rm_cart+transpose(rm_cart))

  end subroutine tas_cn_resolution

  ! return an effective energy resolution width for TAS-CN mode
  subroutine tas_cn_sigmae(qvec_cart, omega, sigmae_eff, ierr)

    implicit none

    real(kind=8), intent(in) :: qvec_cart(:)
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(out) :: sigmae_eff
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: RM(4,4), R0

    ierr = 0
    sigmae_eff = 0.d0

    call tas_cn_resolution(qvec_cart, omega, R0, RM, ierr)
    if (ierr /= 0) return

    call tas_cn_project_sigmae(RM, sigmae_eff, ierr)
    if (ierr /= 0) then
       sigmae_eff = 0.d0
       ierr = 100+ierr
    end if

  end subroutine tas_cn_sigmae

  ! compute the basic TAS geometry from Q, energy transfer and fixed-energy mode
  subroutine tas_cn_geometry(qvec_cart, omega, ki, kf, theta_m, theta_a, two_theta_s, ierr)

    implicit none

    real(kind=8), intent(in) :: qvec_cart(:)
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(out) :: ki
    real(kind=8), intent(out) :: kf
    real(kind=8), intent(out) :: theta_m
    real(kind=8), intent(out) :: theta_a
    real(kind=8), intent(out) :: two_theta_s
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: ei, ef, qmag, cos_2theta, sin_theta_m, sin_theta_a

    ierr = 0
    ki = 0.d0
    kf = 0.d0
    theta_m = 0.d0
    theta_a = 0.d0
    two_theta_s = 0.d0

    if (size(qvec_cart) /= 3) then
       ierr = 1
       return
    end if

    if (trim(tas_fix_mode) == 'Ef') then
       ef = tas_e_fixed
       ei = ef+omega
    elseif (trim(tas_fix_mode) == 'Ei') then
       ei = tas_e_fixed
       ef = ei-omega
    else
       ierr = 2
       return
    end if

    if (ei <= 0.d0 .or. ef <= 0.d0) then
       ierr = 3
       return
    end if

    if (tas_dm <= 0.d0 .or. tas_da <= 0.d0) then
       ierr = 4
       return
    end if

    ki = sqrt(ei/mev_to_k2)
    kf = sqrt(ef/mev_to_k2)
    qmag = sqrt(dot_product(qvec_cart, qvec_cart))

    if (ki <= 0.d0 .or. kf <= 0.d0) then
       ierr = 5
       return
    end if

    sin_theta_m = acos(-1.d0)/(ki*tas_dm)
    sin_theta_a = acos(-1.d0)/(kf*tas_da)
    if (abs(sin_theta_m) > 1.d0 .or. abs(sin_theta_a) > 1.d0) then
       ierr = 6
       return
    end if

    theta_m = asin(sin_theta_m)
    theta_a = asin(sin_theta_a)

    cos_2theta = (ki**2+kf**2-qmag**2)/(2.d0*ki*kf)
    if (cos_2theta < -1.d0 .or. cos_2theta > 1.d0) then
       ierr = 7
       return
    end if

    two_theta_s = acos(cos_2theta)

  end subroutine tas_cn_geometry

  ! compatibility wrapper: keep the old cartesian RM interface
  subroutine tas_cn_rm(qvec_cart, omega, RM, R0, ierr)

    implicit none

    real(kind=8), intent(in) :: qvec_cart(:)
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(out) :: RM(4,4)
    real(kind=8), intent(out) :: R0
    integer(kind=4), intent(out) :: ierr

    call tas_cn_resolution(qvec_cart, omega, R0, RM, ierr)

  end subroutine tas_cn_rm

  subroutine tas_cn_rm_qframe(qmag, omega, RM, R0, ierr)

    implicit none

    real(kind=8), intent(in) :: qmag
    real(kind=8), intent(in) :: omega
    real(kind=8), intent(out) :: RM(4,4)
    real(kind=8), intent(out) :: R0
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: info, infin
    real(kind=8) :: ki, kf, theta_m, theta_a, two_theta_s
    real(kind=8) :: det_h, det_m, det_rm, etam, etamv, etaa, etaav, etas, etasv
    real(kind=8) :: taum, taua, ei, ef, sm, ss, sa, thetas, phi
    real(kind=8) :: rm_tmp, ra_tmp, r0_tmp, tempv, rmon, det_mon
    real(kind=8) :: alpha(4), beta(4)
    real(kind=8) :: minv(4,4)
    real(kind=8) :: G(8,8), F(4,4), A(6,8), C(4,8), B(4,6), H(8,8), Hinv(8,8)
    real(kind=8) :: Ninv(6,6)
    real(kind=8) :: gmon(4,4), fmon(2,2), cmon(2,4), hmon(4,4), tmon(2,7), dmon(4,7)

    ierr = 0
    RM = 0.d0
    R0 = 0.d0

    if (qmag <= 0.d0) then
       ierr = 11
       return
    end if

    taum = 2.d0*pi/tas_dm
    taua = 2.d0*pi/tas_da
    if (taum <= 0.d0 .or. taua <= 0.d0) then
       ierr = 12
       return
    end if

    alpha = (/angle_sigma_arcmin(tas_coll_h_pre_mono), &
              angle_sigma_arcmin(tas_coll_h_pre_samp), &
              angle_sigma_arcmin(tas_coll_h_post_samp), &
              angle_sigma_arcmin(tas_coll_h_post_ana)/)
    beta = angle_sigma_arcmin(default_vcol_arcmin)

    if (trim(tas_fix_mode) == 'Ei') then
       infin = 1
    else
       infin = -1
    end if

    ei = tas_e_fixed
    ef = tas_e_fixed
    if (infin > 0) then
       ef = tas_e_fixed-omega
    else
       ei = tas_e_fixed+omega
    end if
    if (ei <= 0.d0 .or. ef <= 0.d0) then
       ierr = 13
       return
    end if

    ki = sqrt(ei/mev_to_k2)
    kf = sqrt(ef/mev_to_k2)

    sm = dble(tas_mono_dir)
    ss = dble(tas_sample_dir)
    sa = dble(tas_ana_dir)
    tempv = taum/(2.d0*ki)
    if (abs(tempv) > 1.d0) then
       ierr = 14
       return
    end if
    theta_m = asin(tempv)*sm
    tempv = taua/(2.d0*kf)
    if (abs(tempv) > 1.d0) then
       ierr = 15
       return
    end if
    theta_a = asin(tempv)*sa
    tempv = (ki**2+kf**2-qmag**2)/(2.d0*ki*kf)
    if (abs(tempv) > 1.d0) then
       ierr = 16
       return
    end if
    two_theta_s = acos(tempv)*ss
    thetas = two_theta_s/2.d0
    phi = atan2(-kf*sin(two_theta_s), ki-kf*cos(two_theta_s))

    if (tas_use_mosaic) then
       etam = angle_sigma_arcmin(max(tas_mosaic_mono_h, tiny_mosaic_arcmin))
       etaa = angle_sigma_arcmin(max(tas_mosaic_ana_h, tiny_mosaic_arcmin))
    else
       etam = angle_sigma_arcmin(tiny_mosaic_arcmin)
       etaa = angle_sigma_arcmin(tiny_mosaic_arcmin)
    end if
    etamv = etam
    etaav = etaa

    G = 0.d0
    G(1,1) = 1.d0/alpha(1)**2
    G(2,2) = 1.d0/alpha(2)**2
    G(3,3) = 1.d0/beta(1)**2
    G(4,4) = 1.d0/beta(2)**2
    G(5,5) = 1.d0/alpha(3)**2
    G(6,6) = 1.d0/alpha(4)**2
    G(7,7) = 1.d0/beta(3)**2
    G(8,8) = 1.d0/beta(4)**2

    F = 0.d0
    F(1,1) = 1.d0/etam**2
    F(2,2) = 1.d0/etamv**2
    F(3,3) = 1.d0/etaa**2
    F(4,4) = 1.d0/etaav**2

    A = 0.d0
    A(1,:) = (/ki/(2.d0*tan(theta_m)), -ki/(2.d0*tan(theta_m)), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    A(2,:) = (/0.d0, ki, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    A(3,:) = (/0.d0, 0.d0, 0.d0, ki, 0.d0, 0.d0, 0.d0, 0.d0/)
    A(4,:) = (/0.d0, 0.d0, 0.d0, 0.d0, kf/(2.d0*tan(theta_a)), -kf/(2.d0*tan(theta_a)), 0.d0, 0.d0/)
    A(5,:) = (/0.d0, 0.d0, 0.d0, 0.d0, kf, 0.d0, 0.d0, 0.d0/)
    A(6,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, kf, 0.d0/)

    C = 0.d0
    C(1,:) = (/0.5d0, 0.5d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    C(2,:) = (/0.d0, 0.d0, 1.d0/(2.d0*sin(theta_m)), -1.d0/(2.d0*sin(theta_m)), 0.d0, 0.d0, 0.d0, 0.d0/)
    C(3,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.5d0, 0.5d0, 0.d0, 0.d0/)
    C(4,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/(2.d0*sin(theta_a)), -1.d0/(2.d0*sin(theta_a))/)

    B = 0.d0
    B(1,:) = (/cos(phi), sin(phi), 0.d0, -cos(phi-two_theta_s), -sin(phi-two_theta_s), 0.d0/)
    B(2,:) = (/-sin(phi), cos(phi), 0.d0, sin(phi-two_theta_s), -cos(phi-two_theta_s), 0.d0/)
    B(3,:) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0, -1.d0/)
    B(4,:) = (/2.d0*mev_to_k2*ki, 0.d0, 0.d0, -2.d0*mev_to_k2*kf, 0.d0, 0.d0/)

    H = G+matmul(transpose(C), matmul(F, C))
    call invert_matrix(H, Hinv, det_h, info)
    if (info /= 0 .or. det_h <= 0.d0) then
       ierr = 17
       return
    end if

    Ninv = matmul(A, matmul(Hinv, transpose(A)))
    minv = matmul(B, matmul(Ninv, transpose(B)))
    call invert_matrix(minv, RM, det_m, info)
    if (info /= 0 .or. det_m <= 0.d0) then
       ierr = 18
       return
    end if

    rm_tmp = ki**3/tan(theta_m)
    ra_tmp = kf**3/tan(theta_a)
    r0_tmp = rm_tmp*ra_tmp*(2.d0*pi)**4/(64.d0*pi**2*sin(theta_m)*sin(theta_a))
    r0_tmp = r0_tmp*sqrt(det4(F)/det_h)

    ! default NeutronPy/ResLib monitor normalization
    gmon = 0.d0
    gmon(1:4,1:4) = G(1:4,1:4)
    fmon = 0.d0
    fmon(1:2,1:2) = F(1:2,1:2)
    cmon = 0.d0
    cmon(1:2,1:4) = C(1:2,1:4)
    hmon = gmon+matmul(transpose(cmon), matmul(fmon, cmon))
    det_mon = det4(hmon)
    if (det_mon <= 0.d0) then
       ierr = 20
       return
    end if
    rmon = rm_tmp*(2.d0*pi)**2/(8.d0*pi*sin(theta_m))*sqrt(det2(fmon)/det_mon)
    r0_tmp = r0_tmp/rmon
    r0_tmp = r0_tmp*ki

    det_rm = 1.d0/det_m
    if (det_rm <= 0.d0) then
       ierr = 21
       return
    end if
    r0_tmp = r0_tmp/(2.d0*pi)**2*sqrt(det_rm)
    r0_tmp = r0_tmp*kf/ki

    if (tas_use_mosaic .and. tas_mosaic_samp_h > 0.d0) then
       etas = angle_sigma_arcmin(tas_mosaic_samp_h)
       etasv = etas
       r0_tmp = r0_tmp/sqrt((1.d0+(qmag*etas)**2*RM(3,3))*(1.d0+(qmag*etasv)**2*RM(2,2)))
       minv(2,2) = minv(2,2)+(qmag*etas)**2
       minv(3,3) = minv(3,3)+(qmag*etasv)**2
       call invert_matrix(minv, RM, det_m, info)
       if (info /= 0 .or. det_m <= 0.d0) then
          ierr = 19
          return
       end if
       det_rm = 1.d0/det_m
       if (det_rm <= 0.d0) then
          ierr = 22
          return
       end if
    end if

    RM = 0.5d0*(RM+transpose(RM))
    R0 = r0_tmp

  end subroutine tas_cn_rm_qframe

  subroutine tas_cn_debug_qframe(qmag, omega, ki, kf, theta_m, theta_a, two_theta_s, phi, &
                                 G, F, H, Ninv, Minv, RM, R0, ierr)

    implicit none

    real(kind=8), intent(in) :: qmag, omega
    real(kind=8), intent(out) :: ki, kf, theta_m, theta_a, two_theta_s, phi
    real(kind=8), intent(out) :: G(8,8), F(4,4), H(8,8), Ninv(6,6), Minv(4,4), RM(4,4), R0
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: info, infin
    real(kind=8) :: det_h, det_m, det_rm, etam, etamv, etaa, etaav, etas, etasv
    real(kind=8) :: taum, taua, ei, ef, sm, ss, sa, thetas
    real(kind=8) :: rm_tmp, ra_tmp, r0_tmp, tempv, rmon, det_mon
    real(kind=8) :: alpha(4), beta(4), A(6,8), C(4,8), B(4,6), Hinv(8,8), gmon(4,4), fmon(2,2), cmon(2,4), hmon(4,4)

    ierr = 0
    ki = 0.d0
    kf = 0.d0
    theta_m = 0.d0
    theta_a = 0.d0
    two_theta_s = 0.d0
    phi = 0.d0
    G = 0.d0
    F = 0.d0
    H = 0.d0
    Ninv = 0.d0
    Minv = 0.d0
    RM = 0.d0
    R0 = 0.d0

    if (qmag <= 0.d0) then
       ierr = 11
       return
    end if

    taum = 2.d0*pi/tas_dm
    taua = 2.d0*pi/tas_da
    if (taum <= 0.d0 .or. taua <= 0.d0) then
       ierr = 12
       return
    end if

    alpha = (/angle_sigma_arcmin(tas_coll_h_pre_mono), &
              angle_sigma_arcmin(tas_coll_h_pre_samp), &
              angle_sigma_arcmin(tas_coll_h_post_samp), &
              angle_sigma_arcmin(tas_coll_h_post_ana)/)
    beta = angle_sigma_arcmin(default_vcol_arcmin)

    if (trim(tas_fix_mode) == 'Ei') then
       infin = 1
    else
       infin = -1
    end if

    ei = tas_e_fixed
    ef = tas_e_fixed
    if (infin > 0) then
       ef = tas_e_fixed-omega
    else
       ei = tas_e_fixed+omega
    end if
    if (ei <= 0.d0 .or. ef <= 0.d0) then
       ierr = 13
       return
    end if

    ki = sqrt(ei/mev_to_k2)
    kf = sqrt(ef/mev_to_k2)
    sm = dble(tas_mono_dir)
    ss = dble(tas_sample_dir)
    sa = dble(tas_ana_dir)

    tempv = taum/(2.d0*ki)
    if (abs(tempv) > 1.d0) then
       ierr = 14
       return
    end if
    theta_m = asin(tempv)*sm

    tempv = taua/(2.d0*kf)
    if (abs(tempv) > 1.d0) then
       ierr = 15
       return
    end if
    theta_a = asin(tempv)*sa

    tempv = (ki**2+kf**2-qmag**2)/(2.d0*ki*kf)
    if (abs(tempv) > 1.d0) then
       ierr = 16
       return
    end if
    two_theta_s = acos(tempv)*ss
    thetas = two_theta_s/2.d0
    phi = atan2(-kf*sin(two_theta_s), ki-kf*cos(two_theta_s))

    if (tas_use_mosaic) then
       etam = angle_sigma_arcmin(max(tas_mosaic_mono_h, tiny_mosaic_arcmin))
       etaa = angle_sigma_arcmin(max(tas_mosaic_ana_h, tiny_mosaic_arcmin))
    else
       etam = angle_sigma_arcmin(tiny_mosaic_arcmin)
       etaa = angle_sigma_arcmin(tiny_mosaic_arcmin)
    end if
    etamv = etam
    etaav = etaa

    G(1,1) = 1.d0/alpha(1)**2
    G(2,2) = 1.d0/alpha(2)**2
    G(3,3) = 1.d0/beta(1)**2
    G(4,4) = 1.d0/beta(2)**2
    G(5,5) = 1.d0/alpha(3)**2
    G(6,6) = 1.d0/alpha(4)**2
    G(7,7) = 1.d0/beta(3)**2
    G(8,8) = 1.d0/beta(4)**2

    F(1,1) = 1.d0/etam**2
    F(2,2) = 1.d0/etamv**2
    F(3,3) = 1.d0/etaa**2
    F(4,4) = 1.d0/etaav**2

    A = 0.d0
    A(1,:) = (/ki/(2.d0*tan(theta_m)), -ki/(2.d0*tan(theta_m)), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    A(2,:) = (/0.d0, ki, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    A(3,:) = (/0.d0, 0.d0, 0.d0, ki, 0.d0, 0.d0, 0.d0, 0.d0/)
    A(4,:) = (/0.d0, 0.d0, 0.d0, 0.d0, kf/(2.d0*tan(theta_a)), -kf/(2.d0*tan(theta_a)), 0.d0, 0.d0/)
    A(5,:) = (/0.d0, 0.d0, 0.d0, 0.d0, kf, 0.d0, 0.d0, 0.d0/)
    A(6,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, kf, 0.d0/)

    C = 0.d0
    C(1,:) = (/0.5d0, 0.5d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    C(2,:) = (/0.d0, 0.d0, 1.d0/(2.d0*sin(theta_m)), -1.d0/(2.d0*sin(theta_m)), 0.d0, 0.d0, 0.d0, 0.d0/)
    C(3,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.5d0, 0.5d0, 0.d0, 0.d0/)
    C(4,:) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/(2.d0*sin(theta_a)), -1.d0/(2.d0*sin(theta_a))/)

    B = 0.d0
    B(1,:) = (/cos(phi), sin(phi), 0.d0, -cos(phi-two_theta_s), -sin(phi-two_theta_s), 0.d0/)
    B(2,:) = (/-sin(phi), cos(phi), 0.d0, sin(phi-two_theta_s), -cos(phi-two_theta_s), 0.d0/)
    B(3,:) = (/0.d0, 0.d0, 1.d0, 0.d0, 0.d0, -1.d0/)
    B(4,:) = (/2.d0*mev_to_k2*ki, 0.d0, 0.d0, -2.d0*mev_to_k2*kf, 0.d0, 0.d0/)

    H = G+matmul(transpose(C), matmul(F, C))
    call invert_matrix(H, Hinv, det_h, info)
    if (info /= 0 .or. det_h <= 0.d0) then
       ierr = 17
       return
    end if

    Ninv = matmul(A, matmul(Hinv, transpose(A)))
    Minv = matmul(B, matmul(Ninv, transpose(B)))
    call invert_matrix(Minv, RM, det_m, info)
    if (info /= 0 .or. det_m <= 0.d0) then
       ierr = 18
       return
    end if

    rm_tmp = ki**3/tan(theta_m)
    ra_tmp = kf**3/tan(theta_a)
    r0_tmp = rm_tmp*ra_tmp*(2.d0*pi)**4/(64.d0*pi**2*sin(theta_m)*sin(theta_a))
    r0_tmp = r0_tmp*sqrt(det4(F)/det_h)

    gmon = 0.d0
    gmon(1:4,1:4) = G(1:4,1:4)
    fmon = 0.d0
    fmon(1:2,1:2) = F(1:2,1:2)
    cmon = 0.d0
    cmon(1:2,1:4) = C(1:2,1:4)
    hmon = gmon+matmul(transpose(cmon), matmul(fmon, cmon))
    det_mon = det4(hmon)
    if (det_mon <= 0.d0) then
       ierr = 20
       return
    end if
    rmon = rm_tmp*(2.d0*pi)**2/(8.d0*pi*sin(theta_m))*sqrt(det2(fmon)/det_mon)
    r0_tmp = r0_tmp/rmon
    r0_tmp = r0_tmp*ki

    det_rm = 1.d0/det_m
    if (det_rm <= 0.d0) then
       ierr = 21
       return
    end if
    r0_tmp = r0_tmp/(2.d0*pi)**2*sqrt(det_rm)
    r0_tmp = r0_tmp*kf/ki

    if (tas_use_mosaic .and. tas_mosaic_samp_h > 0.d0) then
       etas = angle_sigma_arcmin(tas_mosaic_samp_h)
       etasv = etas
       r0_tmp = r0_tmp/sqrt((1.d0+(qmag*etas)**2*RM(3,3))*(1.d0+(qmag*etasv)**2*RM(2,2)))
       Minv(2,2) = Minv(2,2)+(qmag*etas)**2
       Minv(3,3) = Minv(3,3)+(qmag*etasv)**2
       call invert_matrix(Minv, RM, det_m, info)
       if (info /= 0 .or. det_m <= 0.d0) then
          ierr = 19
          return
       end if
       det_rm = 1.d0/det_m
       if (det_rm <= 0.d0) then
          ierr = 22
          return
       end if
    end if

    RM = 0.5d0*(RM+transpose(RM))
    R0 = r0_tmp

  end subroutine tas_cn_debug_qframe

  subroutine tas_cn_standard_system(a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, x, y, z, rlattice, ierr)

    implicit none

    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    real(kind=8), intent(out) :: x(3), y(3), z(3)
    real(kind=8), intent(out) :: rlattice(6)
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: lattice(6), proj, modx, mody, modz

    ierr = 0
    call reciprocal_lattice_params(a, b, c, alpha_deg, beta_deg, gamma_deg, lattice, rlattice, ierr)
    if (ierr /= 0) return

    modx = modvec_metric(orient1, rlattice)
    if (modx <= 0.d0) then
       ierr = 31
       return
    end if
    x = orient1/modx

    proj = scalar_metric(orient2, x, rlattice)
    y = orient2-x*proj
    mody = modvec_metric(y, rlattice)
    if (mody <= 0.d0) then
       ierr = 32
       return
    end if
    y = y/mody

    z = (/x(2)*y(3)-y(2)*x(3), x(3)*y(1)-y(3)*x(1), -x(2)*y(1)+y(2)*x(1)/)
    proj = scalar_metric(z, x, rlattice)
    z = z-x*proj
    proj = scalar_metric(z, y, rlattice)
    z = z-y*proj
    modz = modvec_metric(z, rlattice)
    if (modz <= 0.d0) then
       ierr = 33
       return
    end if
    z = z/modz

  end subroutine tas_cn_standard_system

  ! unified HKL/sample-frame TAS resolution output for future 4D convolution
  subroutine tas_cn_resolution_hkle(hkl, omega, a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, R0, RM, ierr)

    implicit none

    real(kind=8), intent(in) :: hkl(3), omega
    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    real(kind=8), intent(out) :: R0, RM(4,4)
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: x(3), y(3), z(3), rlattice(6), qmag, uq(3), xq, yq
    real(kind=8) :: tmat(4,4), RMq(4,4)

    ierr = 0
    RM = 0.d0
    R0 = 0.d0

    call tas_cn_standard_system(a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, x, y, z, rlattice, ierr)
    if (ierr /= 0) return

    qmag = modvec_metric(hkl, rlattice)
    if (qmag <= 0.d0) then
       ierr = 34
       return
    end if
    uq = hkl/qmag
    xq = scalar_metric(x, uq, rlattice)
    yq = scalar_metric(y, uq, rlattice)

    tmat = 0.d0
    tmat(1,:) = (/xq, yq, 0.d0, 0.d0/)
    tmat(2,:) = (/-yq, xq, 0.d0, 0.d0/)
    tmat(3,:) = (/0.d0, 0.d0, 1.d0, 0.d0/)
    tmat(4,:) = (/0.d0, 0.d0, 0.d0, 1.d0/)

    call tas_cn_rm_qframe(qmag, omega, RMq, R0, ierr)
    if (ierr /= 0) return

    RM = matmul(transpose(tmat), matmul(RMq, tmat))
    RM = 0.5d0*(RM+transpose(RM))

  end subroutine tas_cn_resolution_hkle

  subroutine tas_cn_rm_hkle(hkl, omega, a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, RMS, R0, ierr)

    implicit none

    real(kind=8), intent(in) :: hkl(3), omega
    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    real(kind=8), intent(out) :: RMS(4,4), R0
    integer(kind=4), intent(out) :: ierr

    call tas_cn_resolution_hkle(hkl, omega, a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, R0, RMS, ierr)

  end subroutine tas_cn_rm_hkle

  subroutine tas_cn_debug_hkle(hkl, omega, a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, &
                               qmag, x, y, z, rlattice, uq, xq, yq, tmat, RMq, RMS, R0, ierr)

    implicit none

    real(kind=8), intent(in) :: hkl(3), omega
    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    real(kind=8), intent(out) :: qmag, x(3), y(3), z(3), rlattice(6), uq(3), xq, yq
    real(kind=8), intent(out) :: tmat(4,4), RMq(4,4), RMS(4,4), R0
    integer(kind=4), intent(out) :: ierr

    ierr = 0
    qmag = 0.d0
    x = 0.d0
    y = 0.d0
    z = 0.d0
    rlattice = 0.d0
    uq = 0.d0
    xq = 0.d0
    yq = 0.d0
    tmat = 0.d0
    RMq = 0.d0
    RMS = 0.d0
    R0 = 0.d0

    call tas_cn_standard_system(a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, x, y, z, rlattice, ierr)
    if (ierr /= 0) return

    qmag = modvec_metric(hkl, rlattice)
    if (qmag <= 0.d0) then
       ierr = 34
       return
    end if
    uq = hkl/qmag
    xq = scalar_metric(x, uq, rlattice)
    yq = scalar_metric(y, uq, rlattice)

    tmat = 0.d0
    tmat(1,:) = (/xq, yq, 0.d0, 0.d0/)
    tmat(2,:) = (/-yq, xq, 0.d0, 0.d0/)
    tmat(3,:) = (/0.d0, 0.d0, 1.d0, 0.d0/)
    tmat(4,:) = (/0.d0, 0.d0, 0.d0, 1.d0/)

    call tas_cn_rm_qframe(qmag, omega, RMq, R0, ierr)
    if (ierr /= 0) return

    RMS = matmul(transpose(tmat), matmul(RMq, tmat))
    RMS = 0.5d0*(RMS+transpose(RMS))

  end subroutine tas_cn_debug_hkle

  subroutine tas_cn_sigmae_hkle(hkl, omega, a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, sigmae_eff, ierr)

    implicit none

    real(kind=8), intent(in) :: hkl(3), omega
    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(in) :: orient1(3), orient2(3)
    real(kind=8), intent(out) :: sigmae_eff
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: RMS(4,4), R0

    sigmae_eff = 0.d0
    call tas_cn_resolution_hkle(hkl, omega, a, b, c, alpha_deg, beta_deg, gamma_deg, orient1, orient2, R0, RMS, ierr)
    if (ierr /= 0) return
    call tas_cn_project_sigmae(RMS, sigmae_eff, ierr)
    if (ierr /= 0) ierr = 200+ierr

  end subroutine tas_cn_sigmae_hkle

  ! project the 4D precision matrix to an effective energy width
  ! after integrating out the three Q directions
  subroutine tas_cn_project_sigmae(RM, sigmae_eff, ierr)

    implicit none

    real(kind=8), intent(in) :: RM(4,4)
    real(kind=8), intent(out) :: sigmae_eff
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: info
    real(kind=8) :: A(3,3), Ainv(3,3), b(3), schur, detA

    ierr = 0
    sigmae_eff = 0.d0

    A = 0.5d0*(RM(1:3,1:3)+transpose(RM(1:3,1:3)))
    b = RM(1:3,4)

    call invert3(A, Ainv, detA, info)
    if (info /= 0 .or. detA <= 0.d0) then
       ierr = 1
       return
    end if

    schur = RM(4,4)-dot_product(b, matmul(Ainv, b))
    if (schur <= 0.d0) then
       ierr = 2
       return
    end if

    sigmae_eff = 1.d0/sqrt(schur)
    if (sigmae_eff <= 0.d0) then
       ierr = 3
       sigmae_eff = 0.d0
    end if

  end subroutine tas_cn_project_sigmae

  function angle_sigma_arcmin(value_arcmin)

    implicit none

    real(kind=8), intent(in) :: value_arcmin
    real(kind=8) :: angle_sigma_arcmin

    angle_sigma_arcmin = abs(value_arcmin)*pi/180.d0/60.d0*fwhm_to_sigma

  end function angle_sigma_arcmin

  subroutine reciprocal_lattice_params(a, b, c, alpha_deg, beta_deg, gamma_deg, lattice, rlattice, ierr)

    implicit none

    real(kind=8), intent(in) :: a, b, c, alpha_deg, beta_deg, gamma_deg
    real(kind=8), intent(out) :: lattice(6), rlattice(6)
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: alpha, beta, gamma, vol

    ierr = 0
    lattice = 0.d0
    rlattice = 0.d0

    if (a <= 0.d0 .or. b <= 0.d0 .or. c <= 0.d0) then
       ierr = 1
       return
    end if

    alpha = alpha_deg*pi/180.d0
    beta = beta_deg*pi/180.d0
    gamma = gamma_deg*pi/180.d0

    vol = 2.d0*a*b*c*sqrt(sin((alpha+beta+gamma)/2.d0)* &
                          sin((-alpha+beta+gamma)/2.d0)* &
                          sin((alpha-beta+gamma)/2.d0)* &
                          sin((alpha+beta-gamma)/2.d0))
    if (vol <= 0.d0) then
       ierr = 2
       return
    end if

    lattice = (/a, b, c, alpha, beta, gamma/)
    rlattice(1) = 2.d0*pi*b*c*sin(alpha)/vol
    rlattice(2) = 2.d0*pi*a*c*sin(beta)/vol
    rlattice(3) = 2.d0*pi*b*a*sin(gamma)/vol
    rlattice(4) = acos((cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma)))
    rlattice(5) = acos((cos(alpha)*cos(gamma)-cos(beta))/(sin(alpha)*sin(gamma)))
    rlattice(6) = acos((cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta)))

  end subroutine reciprocal_lattice_params

  function scalar_metric(v1, v2, lattice)

    implicit none

    real(kind=8), intent(in) :: v1(3), v2(3), lattice(6)
    real(kind=8) :: scalar_metric
    real(kind=8) :: a, b, c, alpha, beta, gamma

    a = lattice(1)
    b = lattice(2)
    c = lattice(3)
    alpha = lattice(4)
    beta = lattice(5)
    gamma = lattice(6)

    scalar_metric = v1(1)*v2(1)*a**2+v1(2)*v2(2)*b**2+v1(3)*v2(3)*c**2+ &
                    (v1(1)*v2(2)+v2(1)*v1(2))*a*b*cos(gamma)+ &
                    (v1(1)*v2(3)+v2(1)*v1(3))*a*c*cos(beta)+ &
                    (v1(3)*v2(2)+v2(3)*v1(2))*c*b*cos(alpha)

  end function scalar_metric

  function modvec_metric(v, lattice)

    implicit none

    real(kind=8), intent(in) :: v(3), lattice(6)
    real(kind=8) :: modvec_metric

    modvec_metric = sqrt(max(scalar_metric(v, v, lattice), 0.d0))

  end function modvec_metric

  subroutine cross3(a, b, c)

    implicit none

    real(kind=8), intent(in) :: a(3), b(3)
    real(kind=8), intent(out) :: c(3)

    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)

  end subroutine cross3

  subroutine normalize3(vec, ierr)

    implicit none

    real(kind=8), intent(inout) :: vec(3)
    integer(kind=4), intent(out) :: ierr

    real(kind=8) :: normv

    ierr = 0
    normv = sqrt(dot_product(vec, vec))
    if (normv <= 0.d0) then
      ierr = 1
      return
    end if
    vec = vec/normv

  end subroutine normalize3

  subroutine invert4(amat, ainv, detv, ierr)

    implicit none

    real(kind=8), intent(in) :: amat(4,4)
    real(kind=8), intent(out) :: ainv(4,4)
    real(kind=8), intent(out) :: detv
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: ii, jj, ipiv
    real(kind=8) :: aug(4,8), pivot, factor, maxv, tmp(8)

    ierr = 0
    detv = 1.d0
    aug = 0.d0
    aug(:,1:4) = amat
    do ii = 1, 4
       aug(ii,4+ii) = 1.d0
    end do

    do ii = 1, 4
       ipiv = ii
       maxv = abs(aug(ii,ii))
       do jj = ii+1, 4
          if (abs(aug(jj,ii)) > maxv) then
             maxv = abs(aug(jj,ii))
             ipiv = jj
          end if
       end do
       if (maxv <= 1.d-20) then
          ierr = 1
          ainv = 0.d0
          detv = 0.d0
          return
       end if
       if (ipiv /= ii) then
          tmp = aug(ii,:)
          aug(ii,:) = aug(ipiv,:)
          aug(ipiv,:) = tmp
          detv = -detv
       end if
       pivot = aug(ii,ii)
       detv = detv*pivot
       aug(ii,:) = aug(ii,:)/pivot
       do jj = 1, 4
          if (jj == ii) cycle
          factor = aug(jj,ii)
          aug(jj,:) = aug(jj,:)-factor*aug(ii,:)
       end do
    end do

    ainv = aug(:,5:8)
  end subroutine invert4

  subroutine invert3(amat, ainv, detv, ierr)

    implicit none

    real(kind=8), intent(in) :: amat(3,3)
    real(kind=8), intent(out) :: ainv(3,3)
    real(kind=8), intent(out) :: detv
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: ii, jj, ipiv
    real(kind=8) :: aug(3,6), pivot, factor, maxv, tmp(6)

    ierr = 0
    detv = 1.d0
    aug = 0.d0
    aug(:,1:3) = amat
    do ii = 1, 3
       aug(ii,3+ii) = 1.d0
    end do

    do ii = 1, 3
       ipiv = ii
       maxv = abs(aug(ii,ii))
       do jj = ii+1, 3
          if (abs(aug(jj,ii)) > maxv) then
             maxv = abs(aug(jj,ii))
             ipiv = jj
          end if
       end do
       if (maxv <= 1.d-20) then
          ierr = 1
          ainv = 0.d0
          detv = 0.d0
          return
       end if
       if (ipiv /= ii) then
          tmp = aug(ii,:)
          aug(ii,:) = aug(ipiv,:)
          aug(ipiv,:) = tmp
          detv = -detv
       end if
       pivot = aug(ii,ii)
       detv = detv*pivot
       aug(ii,:) = aug(ii,:)/pivot
       do jj = 1, 3
          if (jj == ii) cycle
          factor = aug(jj,ii)
          aug(jj,:) = aug(jj,:)-factor*aug(ii,:)
       end do
    end do

    ainv = aug(:,4:6)

  end subroutine invert3

  subroutine invert_matrix(amat, ainv, detv, ierr)

    implicit none

    real(kind=8), intent(in) :: amat(:,:)
    real(kind=8), intent(out) :: ainv(size(amat,1), size(amat,2))
    real(kind=8), intent(out) :: detv
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: n, ii, jj, ipiv
    real(kind=8), allocatable :: aug(:,:), tmp(:)
    real(kind=8) :: pivot, factor, maxv

    ierr = 0
    detv = 1.d0
    n = size(amat,1)
    if (size(amat,2) /= n) then
       ierr = 1
       ainv = 0.d0
       detv = 0.d0
       return
    end if

    allocate(aug(n, 2*n), tmp(2*n))
    aug = 0.d0
    aug(:,1:n) = amat
    do ii = 1, n
       aug(ii,n+ii) = 1.d0
    end do

    do ii = 1, n
       ipiv = ii
       maxv = abs(aug(ii,ii))
       do jj = ii+1, n
          if (abs(aug(jj,ii)) > maxv) then
             maxv = abs(aug(jj,ii))
             ipiv = jj
          end if
       end do
       if (maxv <= 1.d-20) then
          ierr = 2
          ainv = 0.d0
          detv = 0.d0
          deallocate(aug, tmp)
          return
       end if
       if (ipiv /= ii) then
          tmp = aug(ii,:)
          aug(ii,:) = aug(ipiv,:)
          aug(ipiv,:) = tmp
          detv = -detv
       end if
       pivot = aug(ii,ii)
       detv = detv*pivot
       aug(ii,:) = aug(ii,:)/pivot
       do jj = 1, n
          if (jj == ii) cycle
          factor = aug(jj,ii)
          aug(jj,:) = aug(jj,:)-factor*aug(ii,:)
       end do
    end do

    ainv = aug(:,n+1:2*n)
    deallocate(aug, tmp)

  end subroutine invert_matrix

  function det4(amat)

    implicit none

    real(kind=8), intent(in) :: amat(4,4)
    real(kind=8) :: det4, ainv(4,4)
    integer(kind=4) :: ierr

    call invert4(amat, ainv, det4, ierr)
    if (ierr /= 0) det4 = 0.d0

  end function det4

  function det2(amat)

    implicit none

    real(kind=8), intent(in) :: amat(2,2)
    real(kind=8) :: det2

    det2 = amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)

  end function det2

end module tas_resolution_cn
