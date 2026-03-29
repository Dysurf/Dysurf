! Generic fixed-grid 4D TAS resolution convolution following
! neutronpy.instrument.tas_instrument::resolution_convolution(..., METHOD='fix')
module tas_numeric4d_convolution

  implicit none

  abstract interface
     subroutine tas_sqw_fix_interface(h, k, l, w, nmodes, sqw, ierr)
       integer(kind=4), intent(in) :: nmodes
       real(kind=8), intent(in) :: h, k, l, w
       real(kind=8), intent(out) :: sqw(nmodes)
       integer(kind=4), intent(out) :: ierr
     end subroutine tas_sqw_fix_interface

     subroutine tas_pref_fix_interface(h, k, l, w, nmodes, prefactor, bgr, ierr)
       integer(kind=4), intent(in) :: nmodes
       real(kind=8), intent(in) :: h, k, l, w
       real(kind=8), intent(out) :: prefactor(nmodes), bgr
       integer(kind=4), intent(out) :: ierr
     end subroutine tas_pref_fix_interface
  end interface

contains

  subroutine tas_resolution_convolution_fix(q0_hkl, omega0, xvec, yvec, zvec, R0, RM, &
                                            sqw_cb, pref_cb, nmodes, accuracy_xy, accuracy_z, conv, ierr)

    implicit none

    integer(kind=4), intent(in) :: nmodes
    integer(kind=4), intent(in) :: accuracy_xy, accuracy_z
    real(kind=8), intent(in) :: q0_hkl(3), omega0, xvec(3), yvec(3), zvec(3), R0, RM(4,4)
    procedure(tas_sqw_fix_interface) :: sqw_cb
    procedure(tas_pref_fix_interface) :: pref_cb
    real(kind=8), intent(out) :: conv
    integer(kind=4), intent(out) :: ierr

    integer(kind=4) :: ix, iy, iw, iz, nxy, nz, info
    real(kind=8) :: Mxx, Mxy, Mxw, Myy, Myw, Mzz, Mww, MMxx, detM
    real(kind=8) :: tqx, tqyy, tqyx, tqww, tqwy, tqwx, tqz
    real(kind=8) :: step1, step2, phi_x, phi_y, phi_w, phi_z
    real(kind=8) :: tx, ty, tw, tz, dQ1, dQ2, dW, dQ4
    real(kind=8) :: h1, k1, l1, w1, norm, normz, bgr
    real(kind=8) :: prefactor(nmodes), sqw(nmodes)

    ierr = 0
    conv = 0.d0

    Mxx = RM(1, 1)
    Mxy = RM(1, 2)
    Mxw = RM(1, 4)
    Myy = RM(2, 2)
    Myw = RM(2, 4)
    Mzz = RM(3, 3)
    Mww = RM(4, 4)

    if (Mww <= 0.d0 .or. Myy <= 0.d0 .or. Mzz <= 0.d0) then
       ierr = 11
       return
    end if

    Mxx = Mxx - Mxw * Mxw / Mww
    Mxy = Mxy - Mxw * Myw / Mww
    Myy = Myy - Myw * Myw / Mww
    MMxx = Mxx - Mxy * Mxy / Myy
    detM = MMxx * Myy * Mzz * Mww

    if (MMxx <= 0.d0 .or. detM <= 0.d0) then
       ierr = 12
       return
    end if

    tqz = 1.d0 / sqrt(Mzz)
    tqx = 1.d0 / sqrt(MMxx)
    tqyy = 1.d0 / sqrt(Myy)
    tqyx = -Mxy / Myy / sqrt(MMxx)
    tqww = 1.d0 / sqrt(Mww)
    tqwy = -Myw / Mww / sqrt(Myy)
    tqwx = -(Mxw / Mww - Myw / Mww * Mxy / Myy) / sqrt(MMxx)

    call pref_cb(q0_hkl(1), q0_hkl(2), q0_hkl(3), omega0, nmodes, prefactor, bgr, info)
    if (info /= 0) then
       ierr = 20 + info
       return
    end if

    nxy = 2 * accuracy_xy + 1
    nz = 2 * accuracy_z + 1
    step1 = acos(-1.d0) / dble(nxy)
    step2 = acos(-1.d0) / dble(nz)

    do iz = 1, nz
       phi_z = -0.5d0 * acos(-1.d0) + 0.5d0 * step2 + dble(iz - 1) * step2
       tz = tan(phi_z)
       normz = exp(-0.5d0 * tz * tz) * (1.d0 + tz * tz)
       do iw = 1, nxy
          phi_w = -0.5d0 * acos(-1.d0) + 0.5d0 * step1 + dble(iw - 1) * step1
          tw = tan(phi_w)
          do iy = 1, nxy
             phi_y = -0.5d0 * acos(-1.d0) + 0.5d0 * step1 + dble(iy - 1) * step1
             ty = tan(phi_y)
             do ix = 1, nxy
                phi_x = -0.5d0 * acos(-1.d0) + 0.5d0 * step1 + dble(ix - 1) * step1
                tx = tan(phi_x)
                norm = exp(-0.5d0 * (tx * tx + ty * ty + tw * tw)) * &
                       (1.d0 + tx * tx) * (1.d0 + ty * ty) * (1.d0 + tw * tw)

                dQ1 = tqx * tx
                dQ2 = tqyy * ty + tqyx * tx
                dW = tqwx * tx + tqwy * ty + tqww * tw
                dQ4 = tqz * tz

                h1 = q0_hkl(1) + dQ1 * xvec(1) + dQ2 * yvec(1) + dQ4 * zvec(1)
                k1 = q0_hkl(2) + dQ1 * xvec(2) + dQ2 * yvec(2) + dQ4 * zvec(2)
                l1 = q0_hkl(3) + dQ1 * xvec(3) + dQ2 * yvec(3) + dQ4 * zvec(3)
                w1 = omega0 + dW

                call sqw_cb(h1, k1, l1, w1, nmodes, sqw, info)
                if (info /= 0) then
                   ierr = 40 + info
                   return
                end if

                conv = conv + sum(sqw * prefactor) * norm * normz
             end do
          end do
       end do
    end do

    conv = conv * step1 ** 3 * step2 / sqrt(detM)
    if (accuracy_z == 0) conv = conv * 0.79788d0
    if (accuracy_xy == 0) conv = conv * 0.79788d0 ** 3
    conv = conv * R0
    conv = conv + bgr

  end subroutine tas_resolution_convolution_fix

end module tas_numeric4d_convolution
