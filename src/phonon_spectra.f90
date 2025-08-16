! Dysurf, a program for simulating four-dimensional dynamical structure factors
!  Copyright (C) 2023-2025 Yongheng Li <davy_li96@163.com>
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

! Subroutines used to determine phonon band structure
! adapted from phonon_routines.f90 of ShengBTE ( 
! http://www.shengbte.org/ ), licensed under the GPL.

module phonon_spectra

  implicit none

contains
  
  ! phonon spectra driver
  subroutine phoneigen()
    
    use variables, only: natoms, nqtot, nbands, omega, eigenvec, cqlist
    use constants, only: eps1, eps3, tpi, thz2mev
    implicit none

    nbands = 3*natoms
    allocate(omega(nqtot,nbands), eigenvec(nqtot,nbands,nbands))
    omega = 0.d0
    eigenvec = 0.d0
    
    call DMsolver(cqlist(:,:), omega(:,:), eigenvec(:,:,:))

    omega = omega/tpi*thz2mev

    ! set the all imaginary frequency to zero
    if (any(omega .lt. -eps1)) then
       write(*,*) "Warning: the calculated system exists large imaginary frequency."
       write(*,*)
       !write(*,*) "Setting all imaginary frequencies to zero ..."
       !write(*,*)
       where (omega .lt. 0.d0) omega = 0.d0
       where(omega .lt. eps3) omega = omega+eps3
    end if
    
  end subroutine phoneigen

  ! dynamical matrix solver
  subroutine DMsolver(qspace, omegas, eigenvecs)

    use constants, only: thz2amua3, eps5
    use func, only: phexp
    use variables, only: natoms, nbands, nsize, masses2, vol, &
                         q0, nsize, latvec, positions, rlatvec, &
                         fc_short, nonanalytic, epsilon, born
    implicit none

    real(kind=8), intent(in) :: qspace(:,:)
    real(kind=8), intent(out) :: omegas(:,:)
    complex(kind=8), intent(out), optional :: eigenvecs(:,:,:)
 
    integer(kind=4) :: nqpt
    real(kind=8), allocatable :: mm(:,:)
    complex(kind=8), allocatable :: dyn_total(:,:), dyn_nac(:,:)
    real(kind=8), allocatable :: fc_diel(:,:,:,:,:,:,:)
    real(kind=8), allocatable :: fc_total(:,:,:,:,:,:,:)

    integer(kind=4) :: i , j, iq, ip, neq
    integer(kind=4) :: ix1, iy1, iz1, iatom1, ix2, iy2, iz2, iatom2
    real(kind=8) :: tmp1, tmp2, tmp3, dmin, Rnorm
    real(kind=8) :: rcell(3), r(3), rl(3), qr(27)
    complex(kind=8) :: ztmp
    
    real(kind=8), allocatable :: shortest(:,:)
    real(kind=8), allocatable :: omega2(:), rwork(:)
    complex(kind=8), allocatable :: work(:)
    integer(kind=4) :: nwork = 1
    integer(kind=4) :: ct(3)

    real(kind=8), external :: dnrm2

    nqpt = size(qspace,2)

    allocate(mm(natoms,natoms))
    allocate(omega2(nbands))
    allocate(rwork(max(1, 9*natoms-2)))

    allocate(fc_diel(natoms,3,nsize(1),nsize(2),nsize(3),natoms,3))
    allocate(fc_total(natoms,3,nsize(1),nsize(2),nsize(3),natoms,3))

    do i = 1, natoms
       mm(i,i) = masses2(i)
       do j = i+1, natoms
          mm(i,j) = sqrt(masses2(i)*masses2(j))
          mm(j,i) = mm(i,j)
       end do
    end do

    allocate(dyn_total(nbands,nbands))
    allocate(dyn_nac(nbands,nbands))
    allocate(work(nwork))
    allocate(shortest(3,nqpt))

    ! Use the 1st BZ image of each q point to improve the behavior of
    ! the non-analytic correction.
    ! how to deal with LO_TO splitting in non-1st BZ
    ct = floor(q0)
    do iq = 1, nqpt
       shortest(:,iq) = qspace(:,iq)
       tmp1 = dnrm2(3, shortest(:,iq), 1)
       do ix1 = -2-ct(1), 2-ct(1)
          do iy1 = -2-ct(2), 2-ct(2)
             do iz1 = -2-ct(3), 2-ct(3)
                r = qspace(:,iq)+ix1*rlatvec(:,1)+iy1*rlatvec(:,2)+iz1*rlatvec(:,3)
                tmp2 = dnrm2(3, r, 1)
                if (tmp2 .lt. tmp1) then
                   tmp1 = tmp2
                   shortest(:,iq) = r
                end if
             end do
          end do
       end do
    end do

    do iq = 1, nqpt
       dyn_total = 0.d0
       dyn_nac = 0.d0
       fc_diel = 0.d0
       ! If the nonanalytic flag is set to TRUE, add the electrostatic
       ! correction. No correction is applied exactly at \Gamma in
       ! order not to rely on guesses about directions.
       if (nonanalytic .and. (.not. all(shortest(:,iq) .eq. 0.))) then
          tmp3 = dot_product(shortest(:,iq), matmul(epsilon, shortest(:,iq)))
          do iatom1 = 1, natoms
             do iatom2 = 1, natoms
                do i = 1, 3
                   do j = 1, 3
                      tmp1 = dot_product(shortest(:,iq), born(:,i,iatom1))
                      tmp2 = dot_product(shortest(:,iq), born(:,j,iatom2))
                      dyn_nac(3*(iatom1-1)+i,3*(iatom2-1)+j) = tmp1*tmp2/&
                           mm(iatom1,iatom2)
                   end do
                end do
             end do
          end do
          dyn_nac = thz2amua3*dyn_nac/tmp3/vol
          ! Transform back to real space to obtain a correction to the
          ! short-range force constants.
          do iatom1 = 1, natoms
             do iatom2 = 1, natoms
                do i = 1, 3
                   do j = 1, 3
                      fc_diel(iatom1,i,:,:,:,iatom2,j) = real(dyn_nac(3*(iatom1-1)+i, &
                           3*(iatom2-1)+j))
                   end do
                end do
             end do
          end do
          fc_diel=fc_diel/(nsize(1)*nsize(2)*nsize(3))
       end if

       ! Force constants with long-range correction.
       fc_total=fc_short+fc_diel
       
       ! Build the dynamical matrix
       do iatom1 = 1, natoms
          do iatom2 = 1,natoms
             do ix1 = 1, nsize(1)
                do iy1 = 1, nsize(2)
                   do iz1 = 1, nsize(3)
                      rcell = matmul(latvec, (/ix1,iy1,iz1/)-(/1,1,1/))
                      r = positions(:,iatom1)-positions(:,iatom2)+rcell
                      dmin = huge(dmin)
                      do ix2 = -2, 2
                         do iy2 = -2, 2
                            do iz2 = -2, 2
                               rl = ix2*nsize(1)*latvec(:,1)+iy2*nsize(2)*latvec(:,2)+&
                                    iz2*nsize(3)*latvec(:,3)
                               Rnorm = dnrm2(3, rl+r, 1)
                               if (abs(Rnorm-dmin) .gt. eps5) then
                                  if (Rnorm .lt. dmin) then
                                     neq = 1
                                     dmin = Rnorm
                                     qr(neq) = dot_product(qspace(:,iq), rl+r)
                                  end if
                               else
                                  neq = neq+1
                                  qr(neq) = dot_product(qspace(:,iq), rl+r)
                               end if
                            end do
                         end do
                      end do
                      do ip = 1, neq
                         ztmp = phexp(-qr(ip))/neq
                         do i = 1, 3
                            do j = 1, 3
                               dyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j) = &
                                    dyn_total(3*(iatom1-1)+i,3*(iatom2-1)+j)+&
                                    ztmp*fc_total(iatom2,j,ix1,iy1,iz1,iatom1,i)
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do

       ! Frequencies squared result from a diagonalization of the
       ! dynamical matrix. The first call to zheev serves to ensure that
       ! enough space has been allocated for this.
       call zheev("V", "U", nbands, dyn_total, nbands, omega2, work, -1, rwork, i)

       if (real(work(1)) .gt. nwork) then
          nwork = nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V", "U", nbands, dyn_total, nbands, omega2, work, nwork, rwork, i)

       ! Eigenvectors are also returned if required.
       if (present(eigenvecs)) then
          eigenvecs(iq,:,:) = transpose(dyn_total)
       end if

       ! As is conventional, imaginary frequencies are returned as negative.
       omegas(iq,:) = sign(sqrt(abs(omega2)), omega2)

    end do

    deallocate(mm, omega2, rwork, fc_diel, fc_total, &
               dyn_total,dyn_nac, work, shortest)

  end subroutine DMsolver

end module phonon_spectra
