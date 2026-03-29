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
! Material cross-section helpers for sampling phonon mode information at arbitrary Q points.
module material_sqw

  implicit none

contains

  subroutine sqeForGivenQ_material(qpt, omegas, rmsd, eigenvecs, sqetemp, T, mode_mask)

    use func, only: fBEmeV, faff_wk
    use constants, only: ci, eps3
    use variables, only: natoms, nbands, lneutron, lxray, &
                         masses2, coh_b2, xray_b2, aff_wk, lphase, &
                         aff_a, aff_b, aff_c, ntypes, nat, positions
    implicit none

    real(kind=8), intent(in) :: T, qpt(:), omegas(:), rmsd(:,:)
    real(kind=8), intent(out) :: sqetemp(:)
    complex(kind=8), intent(in) :: eigenvecs(:,:)
    logical, intent(in), optional :: mode_mask(:)

    integer(kind=4) :: ii, jj, kk
    real(kind=8) :: DWfac, aff(natoms)
    complex(kind=8) :: temp, tempsum
    sqetemp = 0.d0

    if (lneutron) then
       do ii = 1, nbands
          if (omegas(ii) .lt. eps3) cycle
          if (present(mode_mask)) then
             if (.not. mode_mask(ii)) cycle
          end if
          tempsum = 0.d0
          do jj = 1, natoms
             DWfac = 0.5d0*dot_product(qpt, rmsd(:,jj))**2
             if (lphase) then
                temp = coh_b2(jj)/sqrt(masses2(jj))* &
                       dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3))) * &
                       exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
             else
                temp = coh_b2(jj)/sqrt(masses2(jj))*exp(-DWfac)* &
                       dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))
             end if
             tempsum = tempsum+temp
          end do
          sqetemp(ii) = abs(tempsum)**2*(fBEmeV(omegas(ii), T)+1.d0)/omegas(ii)
       end do
    elseif (lxray) then
       if (aff_wk) then
          kk = 0
          do ii = 1, ntypes
             do jj = 1, nat(ii)
                kk = kk+1
                aff(kk) = faff_wk(aff_a(:,ii), aff_b(:,ii), aff_c(ii), qpt)
             end do
          end do
          do ii = 1, nbands
             if (omegas(ii) .lt. eps3) cycle
             if (present(mode_mask)) then
                if (.not. mode_mask(ii)) cycle
             end if
             tempsum = 0.d0
             do jj = 1, natoms
                DWfac = 0.5d0*dot_product(qpt, rmsd(:,jj))**2
                if (lphase) then
                   temp = aff(jj)/sqrt(masses2(jj))* &
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3))) * &
                          exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
                else
                   temp = aff(jj)/sqrt(masses2(jj))*exp(-DWfac)* &
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))
                end if
                tempsum = tempsum+temp
             end do
             sqetemp(ii) = abs(tempsum)**2*(fBEmeV(omegas(ii), T)+1.d0)/omegas(ii)
          end do
       else
          do ii = 1, nbands
             if (omegas(ii) .lt. eps3) cycle
             if (present(mode_mask)) then
                if (.not. mode_mask(ii)) cycle
             end if
             tempsum = 0.d0
             do jj = 1, natoms
                DWfac = 0.5d0*dot_product(qpt, rmsd(:,jj))**2
                if (lphase) then
                   temp = xray_b2(jj)/sqrt(masses2(jj))* &
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3))) * &
                          exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
                else
                   temp = xray_b2(jj)/sqrt(masses2(jj))*exp(-DWfac)* &
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))
                end if
                tempsum = tempsum+temp
             end do
             sqetemp(ii) = abs(tempsum)**2*(fBEmeV(omegas(ii), T)+1.d0)/omegas(ii)
          end do
       end if
    end if

  end subroutine sqeForGivenQ_material

  subroutine material_mode_info_at_q(qpt_cart, rmsd, temp, omegas_out, intensities_out, ierr, mode_mask)

    use phonon_spectra, only: DMsolver
    use constants, only: tpi, thz2mev
    use variables, only: nbands
    implicit none

    real(kind=8), intent(in) :: qpt_cart(3), rmsd(:,:), temp
    real(kind=8), intent(out) :: omegas_out(nbands), intensities_out(nbands)
    integer(kind=4), intent(out) :: ierr
    logical, intent(in), optional :: mode_mask(:)

    real(kind=8) :: qspace(3,1), omegas(1,nbands), sqetemp(nbands)
    complex(kind=8) :: evecs(1,nbands,nbands)

    ierr = 0
    omegas_out = 0.d0
    intensities_out = 0.d0

    qspace(:,1) = qpt_cart
    call DMsolver(qspace, omegas, evecs)
    omegas(1,:) = omegas(1,:)/tpi*thz2mev
    call sqeForGivenQ_material(qpt_cart, omegas(1,:), rmsd, evecs(1,:,:), sqetemp, temp, mode_mask)
    omegas_out = omegas(1,:)
    intensities_out = sqetemp

  end subroutine material_mode_info_at_q

end module material_sqw
