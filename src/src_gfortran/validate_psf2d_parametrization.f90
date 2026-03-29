program validate_psf2d_parametrization

  use variables, only: tas_mode, tas_fix_mode, tas_use_cn, tas_use_mosaic, tas_debug_res, &
                       tas_e_fixed, tas_coll_h_pre_mono, tas_coll_h_pre_samp, tas_coll_h_post_samp, &
                       tas_coll_h_post_ana, tas_mosaic_mono_h, tas_mosaic_ana_h, tas_mosaic_samp_h, &
                       tas_dm, tas_da, tas_mono_dir, tas_sample_dir, tas_ana_dir
  use tas_resolution_cn, only: tas_cn_resolution_hkle
  use psf_2d_projection, only: build_raw_2d_psf_from_rm
  use psf_2d_param_estimation, only: estimate_psf_params_from_2d_kernel
  use psf_2d_kernel, only: build_local_psf_kernel

  implicit none

  integer(kind=4) :: iq, ie, ios, ierr, argc
  real(kind=8) :: h0, e0
  real(kind=8) :: hkl(3), orient1(3), orient2(3), RM(4,4), R0
  real(kind=8) :: a, b, c, alpha, beta, gamma
  real(kind=8) :: sigma_q, sigma_el, sigma_er, shear
  real(kind=8) :: cov22(2,2), dq_std, de_std, dq_max, de_max, rel_l2, rmse
  real(kind=8), allocatable :: q_grid(:), e_grid(:), raw_kernel(:,:), fit_kernel(:,:)
  real(kind=8), allocatable :: raw_dq_cut(:), fit_dq_cut(:), raw_de_cut(:), fit_de_cut(:)
  logical :: used_fallback
  integer(kind=4) :: nq, ne, jq, je
  character(len=256) :: arg, preset, prefix

  argc = command_argument_count()
  if (argc < 2) then
     write(*,*) 'Usage: ./validate_psf2d_parametrization <iq> <iE> [preset] [prefix] [u1 u2 u3 v1 v2 v3]'
     stop 1
  end if

  call get_command_argument(1, arg)
  read(arg, *, iostat=ios) iq
  if (ios /= 0) stop 2
  call get_command_argument(2, arg)
  read(arg, *, iostat=ios) ie
  if (ios /= 0) stop 3
  preset = 'base'
  prefix = 'psf2d_example'
  if (argc >= 3) then
     call get_command_argument(3, preset)
  end if
  if (argc >= 4) then
     call get_command_argument(4, prefix)
  end if

  a = 5.7409999999999997d0
  b = 5.8769999999999998d0
  c = 7.7050000000000001d0
  alpha = 90.d0
  beta = 90.d0
  gamma = 90.d0
  orient1 = (/1.d0, 0.d0, 0.d0/)
  orient2 = (/0.d0, 1.d0, 0.d0/)
  if (argc >= 10) then
     call get_command_argument(5, arg)
     read(arg, *, iostat=ios) orient1(1)
     if (ios /= 0) stop 10
     call get_command_argument(6, arg)
     read(arg, *, iostat=ios) orient1(2)
     if (ios /= 0) stop 11
     call get_command_argument(7, arg)
     read(arg, *, iostat=ios) orient1(3)
     if (ios /= 0) stop 12
     call get_command_argument(8, arg)
     read(arg, *, iostat=ios) orient2(1)
     if (ios /= 0) stop 13
     call get_command_argument(9, arg)
     read(arg, *, iostat=ios) orient2(2)
     if (ios /= 0) stop 14
     call get_command_argument(10, arg)
     read(arg, *, iostat=ios) orient2(3)
     if (ios /= 0) stop 15
  end if

  h0 = 2.d0 + 0.005d0 * dble(iq - 1)
  e0 = 0.05d0 * dble(ie - 1)
  hkl = (/h0, 0.d0, 0.d0/)

  tas_mode = 'CN'
  tas_fix_mode = 'Ef'
  tas_use_cn = .true.
  tas_use_mosaic = .true.
  tas_debug_res = .false.
  tas_mono_dir = 1
  tas_sample_dir = -1
  tas_ana_dir = 1
  tas_e_fixed = 14.7d0
  if (trim(preset) == 'tight') then
     tas_coll_h_pre_mono = 20.d0
     tas_coll_h_pre_samp = 20.d0
     tas_coll_h_post_samp = 20.d0
     tas_coll_h_post_ana = 20.d0
     tas_mosaic_mono_h = 20.d0
     tas_mosaic_ana_h = 20.d0
     tas_mosaic_samp_h = 30.d0
  elseif (trim(preset) == 'hcol_narrow') then
     tas_coll_h_pre_mono = 20.d0
     tas_coll_h_pre_samp = 20.d0
     tas_coll_h_post_samp = 20.d0
     tas_coll_h_post_ana = 20.d0
     tas_mosaic_mono_h = 35.d0
     tas_mosaic_ana_h = 35.d0
     tas_mosaic_samp_h = 60.d0
  elseif (trim(preset) == 'hcol_base') then
     tas_coll_h_pre_mono = 50.d0
     tas_coll_h_pre_samp = 80.d0
     tas_coll_h_post_samp = 50.d0
     tas_coll_h_post_ana = 120.d0
     tas_mosaic_mono_h = 35.d0
     tas_mosaic_ana_h = 35.d0
     tas_mosaic_samp_h = 60.d0
  elseif (trim(preset) == 'hcol_wide') then
     tas_coll_h_pre_mono = 80.d0
     tas_coll_h_pre_samp = 120.d0
     tas_coll_h_post_samp = 80.d0
     tas_coll_h_post_ana = 120.d0
     tas_mosaic_mono_h = 35.d0
     tas_mosaic_ana_h = 35.d0
     tas_mosaic_samp_h = 60.d0
  elseif (trim(preset) == 'loose') then
     tas_coll_h_pre_mono = 80.d0
     tas_coll_h_pre_samp = 120.d0
     tas_coll_h_post_samp = 80.d0
     tas_coll_h_post_ana = 120.d0
     tas_mosaic_mono_h = 50.d0
     tas_mosaic_ana_h = 50.d0
     tas_mosaic_samp_h = 80.d0
  else
     tas_coll_h_pre_mono = 50.d0
     tas_coll_h_pre_samp = 80.d0
     tas_coll_h_post_samp = 50.d0
     tas_coll_h_post_ana = 120.d0
     tas_mosaic_mono_h = 35.d0
     tas_mosaic_ana_h = 35.d0
     tas_mosaic_samp_h = 60.d0
  end if
  tas_dm = 3.3541627156970963d0
  tas_da = 3.3541627156970963d0

  call tas_cn_resolution_hkle(hkl, e0, a, b, c, alpha, beta, gamma, orient1, orient2, R0, RM, ierr)
  if (ierr /= 0) then
     write(*,*) 'Failed to build TAS CN local RM, ierr=', ierr
     stop 4
  end if

  call build_projected_raw_precision(RM, cov22, ierr)
  if (ierr /= 0) then
     write(*,*) 'Failed to build projected 2D precision, ierr=', ierr
     stop 5
  end if

  dq_std = sqrt(max(cov22(1,1), 1.d-12))
  de_std = sqrt(max(cov22(2,2), 1.d-12))
  dq_max = 4.d0 * dq_std
  de_max = 4.d0 * de_std
  nq = 121
  ne = 161

  allocate(q_grid(nq), e_grid(ne))
  call linspace(-dq_max, dq_max, q_grid)
  call linspace(-de_max, de_max, e_grid)

  call build_raw_2d_psf_from_rm(RM, q_grid, e_grid, raw_kernel, ierr)
  if (ierr /= 0) then
     write(*,*) 'Failed to build raw kernel, ierr=', ierr
     stop 6
  end if

  call estimate_psf_params_from_2d_kernel(q_grid, e_grid, raw_kernel, 1.d0, 1.5d0, 1.5d0, 0.d0, &
                                          sigma_q, sigma_el, sigma_er, shear, used_fallback, ierr)
  if (ierr /= 0) then
     write(*,*) 'Failed to estimate PSF parameters from raw kernel, ierr=', ierr
     stop 7
  end if

  call build_local_psf_kernel(q_grid, e_grid, sigma_q, sigma_el, sigma_er, shear, fit_kernel, ierr)
  if (ierr /= 0) then
     write(*,*) 'Failed to build fitted kernel, ierr=', ierr
     stop 8
  end if

  rel_l2 = sqrt(sum((fit_kernel - raw_kernel) ** 2)) / max(sqrt(sum(raw_kernel ** 2)), 1.d-30)
  rmse = sqrt(sum((fit_kernel - raw_kernel) ** 2) / dble(size(raw_kernel)))

  allocate(raw_dq_cut(nq), fit_dq_cut(nq), raw_de_cut(ne), fit_de_cut(ne))
  raw_dq_cut = raw_kernel(ne / 2 + 1, :)
  fit_dq_cut = fit_kernel(ne / 2 + 1, :)
  raw_de_cut = raw_kernel(:, nq / 2 + 1)
  fit_de_cut = fit_kernel(:, nq / 2 + 1)

  open(unit=10, file=trim(prefix)//'_raw.dat', status='replace', action='write')
  do je = 1, ne
     do jq = 1, nq
        write(10, '(3(E24.16,1x))') q_grid(jq), e_grid(je), raw_kernel(je, jq)
     end do
  end do
  close(10)

  open(unit=11, file=trim(prefix)//'_fit.dat', status='replace', action='write')
  do je = 1, ne
     do jq = 1, nq
        write(11, '(3(E24.16,1x))') q_grid(jq), e_grid(je), fit_kernel(je, jq)
     end do
  end do
  close(11)

  open(unit=12, file=trim(prefix)//'_dq_cut.dat', status='replace', action='write')
  do jq = 1, nq
     write(12, '(3(E24.16,1x))') q_grid(jq), raw_dq_cut(jq), fit_dq_cut(jq)
  end do
  close(12)

  open(unit=13, file=trim(prefix)//'_dE_cut.dat', status='replace', action='write')
  do je = 1, ne
     write(13, '(3(E24.16,1x))') e_grid(je), raw_de_cut(je), fit_de_cut(je)
  end do
  close(13)

  open(unit=14, file=trim(prefix)//'_params.txt', status='replace', action='write')
  write(14, '(A,3(E24.16,1x))') 'q0_hkl = ', hkl
  write(14, '(A,E24.16)') 'E0 = ', e0
  write(14, '(A,I0)') 'iq = ', iq
  write(14, '(A,I0)') 'iE = ', ie
  write(14, '(A,A)') 'preset = ', trim(preset)
  write(14, '(A,3(E24.16,1x))') 'u = ', orient1
  write(14, '(A,3(E24.16,1x))') 'v = ', orient2
  write(14, '(A,4(E24.16,1x))') 'hcol_arcmin = ', tas_coll_h_pre_mono, tas_coll_h_pre_samp, &
       tas_coll_h_post_samp, tas_coll_h_post_ana
  write(14, '(A,3(E24.16,1x))') 'mosaic_arcmin = ', tas_mosaic_mono_h, tas_mosaic_ana_h, tas_mosaic_samp_h
  write(14, '(A,E24.16)') 'sigma_q = ', sigma_q
  write(14, '(A,E24.16)') 'sigma_e_left = ', sigma_el
  write(14, '(A,E24.16)') 'sigma_e_right = ', sigma_er
  write(14, '(A,E24.16)') 'shear = ', shear
  write(14, '(A,L1)') 'used_fallback = ', used_fallback
  write(14, '(A,E24.16)') 'relative_l2 = ', rel_l2
  write(14, '(A,E24.16)') 'rmse = ', rmse
  write(14, '(A)') ''
  write(14, '(A)') 'raw 2D PSF: projected/marginalized Gaussian derived directly from local TAS CN RM'
  write(14, '(A)') 'fit 2D PSF: sheared 2D analytic PSF reconstructed from sigma_q, sigma_e_left, sigma_e_right, shear'
  close(14)

contains

  subroutine linspace(xmin, xmax, arr)
    real(kind=8), intent(in) :: xmin, xmax
    real(kind=8), intent(out) :: arr(:)
    integer(kind=4) :: i, n
    n = size(arr)
    if (n == 1) then
       arr(1) = xmin
       return
    end if
    do i = 1, n
       arr(i) = xmin + (xmax - xmin) * dble(i - 1) / dble(n - 1)
    end do
  end subroutine linspace

  subroutine build_projected_raw_precision(RM, cov22, ierr)
    real(kind=8), intent(in) :: RM(4,4)
    real(kind=8), intent(out) :: cov22(2,2)
    integer(kind=4), intent(out) :: ierr
    real(kind=8) :: A(2,2), B(2,2), C(2,2), Cinv(2,2), P2(2,2), detc, detp
    ierr = 0
    A = reshape((/RM(1,1), RM(1,4), RM(4,1), RM(4,4)/), (/2,2/))
    B = reshape((/RM(1,2), RM(1,3), RM(4,2), RM(4,3)/), (/2,2/))
    C = reshape((/RM(2,2), RM(2,3), RM(3,2), RM(3,3)/), (/2,2/))
    call invert2(C, Cinv, detc, ierr)
    if (ierr /= 0 .or. detc <= 0.d0) then
       ierr = 1
       return
    end if
    P2 = A - matmul(B, matmul(Cinv, transpose(B)))
    call invert2(P2, cov22, detp, ierr)
    if (ierr /= 0 .or. detp <= 0.d0) then
       ierr = 2
       return
    end if
  end subroutine build_projected_raw_precision

  subroutine invert2(A, Ainv, detv, ierr)
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
  end subroutine invert2

end program validate_psf2d_parametrization
