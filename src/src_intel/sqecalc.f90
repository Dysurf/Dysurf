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

! SQE calculator

module sqe_calculator

  implicit none

contains
  
  ! SQE driver
  subroutine sqecalc()
  
    use func, only: resfunc, gauss
    use tas_resolution_cn, only: tas_cn_sigmae_hkle, tas_cn_resolution_hkle, tas_cn_standard_system
    use tas_local_linearization, only: tas_branch_local_model
    use tas_observation_operator, only: calc_sobs
    use psf_2d_state, only: psf2d_setup_state, intrinsic_slice, psf2d_slice, sigma_q_field, &
                            sigma_e_left_field, sigma_e_right_field, shear_field
    use psf_2d_resolution_driver, only: deposit_intrinsic_mode, build_psf2d_observed_slice
    use mpi_helper, only: mpi_rank, mpi_size, mpi_allreduce_sum_real8_5d, mpi_allreduce_sum_real8_3d
    use constants, only: tpi, eps5
    use iso_fortran_env, only: error_unit, output_unit
    use variables, only: ne, deltaE, nqh, nqk, nql, ntemps, &
                         temps, sqebin, sqesum, omega, cqlist, qlist, path, u, v, clatvec, qdist, elist, &
                         nbands, lresfunc, degauss, order, paras, &
                         functype, xm, nqtot, eigenvec, armsd, ltds, l4d, &
                         tas_mode, tas_obs_mode, tas_use_cn, tas_debug_res, &
                         lresfunc2d_psf, lsave_intrinsic_slice, lsave_psf_slice, lsave_psf_params, &
                         psf_mode, psf_window_q, psf_window_e, psf_sigma_q_default, &
                         psf_sigma_e_left_default, psf_sigma_e_right_default, psf_shear_default
    implicit none
    
    integer(kind=4) :: ih, ik, il, ii, jj, kk, nq
    integer(kind=4) :: ierr_cn
    integer(kind=4) :: ierr_psf
    real(kind=8) :: smear, ee(ne), sigmae, ecenter
    real(kind=8) :: alat, blat, clat, alpha_deg, beta_deg, gamma_deg
    real(kind=8) :: orient1(3), orient2(3), psf_orient1(3), psf_orient2(3), sample_x(3), sample_y(3), sample_z(3), rlattice(6)
    real(kind=8) :: grad_hkl(3), grad_sample(3), e0_local, RMcn(4,4), R0cn, sobs
    real(kind=8) :: sqetemp(nbands)
    logical :: use_tas_cn_branch, cn_branch_reported
    logical :: use_psf2d_mode, need_psf2d_outputs
    logical :: intrinsic_reported, rm_reported, params_reported, conv_reported
    real(kind=8) :: rel_l2_psf
    
    allocate(sqebin(ne,nql,nqk,nqh,ntemps))
    if (ltds) then
       allocate(sqesum(nqk,nqh,ntemps))
    else
       allocate(sqesum(ne,nqh,ntemps))
    end if
    sqebin = 0.d0
    sqesum = 0.d0
    cn_branch_reported = .FALSE.
    intrinsic_reported = .FALSE.
    rm_reported = .FALSE.
    params_reported = .FALSE.
    conv_reported = .FALSE.
    use_psf2d_mode = (lresfunc .and. trim(tas_mode) .eq. 'CN' .and. tas_use_cn .and. trim(tas_obs_mode) .eq. 'psf2d')
    need_psf2d_outputs = use_psf2d_mode .or. lresfunc2d_psf

    alat = sqrt(dot_product(clatvec(:,1), clatvec(:,1)))
    blat = sqrt(dot_product(clatvec(:,2), clatvec(:,2)))
    clat = sqrt(dot_product(clatvec(:,3), clatvec(:,3)))
    alpha_deg = acos(dot_product(clatvec(:,2), clatvec(:,3))/blat/clat)*180.d0/tpi*2.d0
    beta_deg = acos(dot_product(clatvec(:,1), clatvec(:,3))/alat/clat)*180.d0/tpi*2.d0
    gamma_deg = acos(dot_product(clatvec(:,1), clatvec(:,2))/alat/blat)*180.d0/tpi*2.d0
    orient1 = path(:,1)
    orient2 = path(:,2)
    psf_orient1 = orient1
    psf_orient2 = orient2
    if (any(abs(u) > eps5)) psf_orient1 = u
    if (any(abs(v) > eps5)) psf_orient2 = v
    call tas_cn_standard_system(alat, blat, clat, alpha_deg, beta_deg, gamma_deg, orient1, orient2, &
                                sample_x, sample_y, sample_z, rlattice, ierr_cn)
    if (ierr_cn /= 0) then
       write(*,*) 'Failed to build TAS sample standard system in sqecalc.'
       stop
    end if

    ee = (/(deltaE*(ii-1), ii = 1, ne)/)
    if (need_psf2d_outputs) call psf2d_setup_state(ne, nqh, ntemps, qdist, ee)
    do ii = 1, ntemps
       nq = 0
       do ih = 1, nqh
          do ik = 1, nqk
             do il = 1, nql
                nq = nq+1
                if (mod(nq-1, mpi_size) == mpi_rank) then
                   call sqeForGivenQ(cqlist(:,nq), omega(nq,:), armsd(:,:,ii), &
                                     eigenvec(nq,:,:), sqetemp, temps(ii))
                   do jj = 1, nbands
                      ecenter = omega(nq,jj)
                      if (need_psf2d_outputs) then
                         call deposit_intrinsic_mode(intrinsic_slice(:, ih, ii), ecenter, sqetemp(jj), deltaE)
                      end if
                      if (use_psf2d_mode) cycle
                      use_tas_cn_branch = (lresfunc .and. trim(tas_mode) .eq. 'CN' .and. tas_use_cn)
                      if (use_tas_cn_branch) then
                         call tas_cn_resolution_hkle(qlist(:,nq), ecenter, alat, blat, clat, alpha_deg, beta_deg, gamma_deg, &
                                                     orient1, orient2, R0cn, RMcn, ierr_cn)
                         use_tas_cn_branch = (ierr_cn .eq. 0)
                         if (use_tas_cn_branch .and. trim(tas_obs_mode) .eq. 'linearized_4d') then
                            call tas_branch_local_model(ih, ik, il, jj, qlist, omega, nqh, nqk, nql, &
                                                       alat, blat, clat, alpha_deg, beta_deg, gamma_deg, orient1, orient2, &
                                                       e0_local, grad_hkl, grad_sample, ierr_cn)
                            use_tas_cn_branch = (ierr_cn .eq. 0)
                         else
                            e0_local = ecenter
                            grad_hkl = 0.d0
                            grad_sample = 0.d0
                         end if
                         if (use_tas_cn_branch) then
                            call tas_cn_sigmae_hkle(qlist(:,nq), ecenter, alat, blat, clat, alpha_deg, beta_deg, gamma_deg, &
                                                    orient1, orient2, sigmae, ierr_cn)
                            use_tas_cn_branch = (ierr_cn .eq. 0)
                         end if
                         if (use_tas_cn_branch .and. tas_debug_res .and. (.not. cn_branch_reported) .and. mpi_rank == 0) then
                            write(*,*) 'TAS observation branch enabled in sqecalc.'
                            write(*,*) 'mode=', trim(tas_obs_mode), ' Q(hkl)=', qlist(:,nq), ' E=', ecenter, ' sigmae=', sigmae
                            if (trim(tas_obs_mode) .eq. 'linearized_4d') then
                               write(*,*) 'local grad_sample=', grad_sample
                            end if
                            write(*,*)
                            cn_branch_reported = .TRUE.
                         end if
                      end if
                      do kk = 1, ne
                         if (use_tas_cn_branch) then
                            call calc_sobs(cqlist(:,nq), ee(kk), sqetemp(jj), e0_local, R0cn, RMcn, sobs, trim(tas_obs_mode), ierr_cn, &
                                           grad_sample, qlist(:,nq), sample_x, sample_y, sample_z, armsd(:,:,ii), temps(ii), jj)
                            if (ierr_cn /= 0) then
                               sobs = sqetemp(jj) * gauss(ee(kk), ecenter, sigmae)
                            end if
                            sqebin(kk,il,ik,ih,ii) = sqebin(kk,il,ik,ih,ii) + sobs
                         elseif (lresfunc) then
                            smear = resfunc(ee(kk), ecenter, functype, xm, order, paras)
                            sqebin(kk,il,ik,ih,ii) = sqebin(kk,il,ik,ih,ii)+sqetemp(jj)*smear
                         else
                            smear = gauss(ee(kk), ecenter, degauss)
                            sqebin(kk,il,ik,ih,ii) = sqebin(kk,il,ik,ih,ii)+sqetemp(jj)*smear
                         end if
                      end do
                   end do
                end if
                if (mpi_rank == 0) then
                   if (mod(nq, 1000) .eq. 0) then
                      write(*,*) "SQE info:", int(nq*100./nqtot), "percent of", nqtot, "q-points"
                      write(*,*)
                      call flush(output_unit)
                   end if
                end if
             end do
         end do
       end do
    end do

    call mpi_allreduce_sum_real8_5d(sqebin)
    if (need_psf2d_outputs) call mpi_allreduce_sum_real8_3d(intrinsic_slice)
    
    if (.not. use_psf2d_mode) then
       where(sqebin < eps5) sqebin = sqebin+eps5
    end if

    if (need_psf2d_outputs .and. mpi_rank == 0) then
       do ii = 1, ntemps
          call build_psf2d_observed_slice(qlist, elist, ee, intrinsic_slice(:,:,ii), &
                                          alat, blat, clat, alpha_deg, beta_deg, gamma_deg, &
                                          psf_orient1, psf_orient2, psf_mode, &
                                          psf_sigma_q_default, psf_sigma_e_left_default, &
                                          psf_sigma_e_right_default, psf_shear_default, &
                                          psf_window_q, psf_window_e, &
                                          psf2d_slice(:,:,ii), sigma_q_field(:,:,ii), sigma_e_left_field(:,:,ii), &
                                          sigma_e_right_field(:,:,ii), shear_field(:,:,ii), &
                                          rel_l2_psf, intrinsic_reported, rm_reported, params_reported, conv_reported, ierr_psf)
          if (ierr_psf /= 0) then
             write(error_unit, *) 'Failed to build 2D PSF-convolved slice in sqecalc.'
             stop
          end if
       end do
    end if
    
    if (use_psf2d_mode) then
       if (mpi_rank == 0) then
          sqesum = psf2d_slice
       else
          sqesum = 0.d0
       end if
    elseif (ltds) then
       forall(ii = 1:ntemps, jj = 1:nqh, kk = 1:nqk)
          sqesum(kk,jj,ii) = sum(reshape(sqebin(:,:,kk,jj,ii),(/nql*ne/)))
       end forall
    elseif (.not.l4d) then
       forall(ii = 1:ntemps, jj = 1:nqh, kk = 1:ne)
          sqesum(kk,jj,ii) = sum(reshape(sqebin(kk,:,:,jj,ii),(/nqk*nql/)))
       end forall
    end if
    
    deallocate(armsd)

  end subroutine sqecalc
  ! calculate the SQE at a given Q point
  ! Only one phonon emission is calculated
  subroutine sqeForGivenQ(qpt, omegas, rmsd, eigenvecs, sqetemp, T)
  
    use func, only: fBEmeV, faff_wk
    use constants, only: ci, eps3
    use variables, only: natoms, nbands, lneutron, lxray, &
                         masses2, coh_b2, xray_b2, aff_wk, lphase, &
                         aff_a, aff_b, aff_c, ntypes, nat, positions
    implicit none
    
    real(kind=8), intent(in) :: T, qpt(:), omegas(:), rmsd(:,:)
    real(kind=8), intent(out) :: sqetemp(:)
    complex(kind=8), intent(in) :: eigenvecs(:,:)

    integer(kind=4) :: ii, jj, kk
    real(kind=8) :: DWfac, aff(natoms)
    complex(kind=8) :: temp, tempsum

    if (lneutron) then
       do ii = 1, nbands
          if (omegas(ii) .lt. eps3) cycle
          tempsum = 0.d0
          do jj = 1, natoms
             DWfac = 0.5d0*dot_product(qpt, rmsd(:,jj))**2
             if (lphase) then
                temp = coh_b2(jj)/sqrt(masses2(jj))*&
                       dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))&
                       *exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
             else
                temp = coh_b2(jj)/sqrt(masses2(jj))*exp(-DWfac)*&
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
             tempsum = 0.d0
             do jj = 1, natoms
                DWfac = 0.5d0*dot_product(qpt, rmsd(:,jj))**2
                if (lphase) then
                   temp = aff(jj)/sqrt(masses2(jj))*&
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))&
                          *exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
                else
                   temp = aff(jj)/sqrt(masses2(jj))*exp(-DWfac)*&
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))
                          !exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
                end if
                tempsum = tempsum+temp
             end do
             sqetemp(ii) = abs(tempsum)**2*(fBEmeV(omegas(ii), T)+1.d0)/omegas(ii)
          end do
       else
          do ii = 1, nbands
             if (omegas(ii) .lt. eps3) cycle
             tempsum = 0.d0
             do jj = 1, natoms
                DWfac = 0.5d0*dot_product(qpt, rmsd(:,jj))**2
                if (lphase) then
                   temp = xray_b2(jj)/sqrt(masses2(jj))*&
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))&
                          *exp(ci*dot_product(qpt, positions(:,jj))-DWfac)                  
                else
                   temp = xray_b2(jj)/sqrt(masses2(jj))*exp(-DWfac)*&
                          dot_product(qpt, eigenvecs(ii,((jj-1)*3+1):((jj-1)*3+3)))
                       !exp(ci*dot_product(qpt, positions(:,jj))-DWfac)
                end if
                tempsum = tempsum+temp
             end do
             sqetemp(ii) = abs(tempsum)**2*(fBEmeV(omegas(ii), T)+1.d0)/omegas(ii)
          end do
       end if
    end if
  
  end subroutine sqeForGivenQ

  ! This subroutine give Mass for experiment(advice)
  ! We calculate as 1cm^2,so the total cross section should be 0.1cm^2
  subroutine given_Mass_reference(mass_ref,sigma_abs_total)
          use variables, only:ntypes,natoms,masses2,scatt_xs2,abs_xs2
      use constants, only: afj,barn
      implicit none
          real(kind=8), intent(out) ::mass_ref,sigma_abs_total
          real(kind=8) ::mass_molecule, tcro_molecule,abscro_molecule
          real(kind=8) :: ii

          mass_molecule=0.d0
          tcro_molecule=0.d0
          tcro_molecule=0.d0
      ! write(*,*) scatt_xs2
      ! write(*,*) abs_xs2
      ! first of all,get mass, etc
        do ii = 1, natoms
            mass_molecule = masses2(ii)+mass_molecule
            tcro_molecule = scatt_xs2(ii)+tcro_molecule
            abscro_molecule= abs_xs2(ii)+abscro_molecule
        end do
         write(*,*) 'Relative molecular mass(MW) is            ', mass_molecule
         write(*,*) 'total atoms_scatt_totoal_cross_section is ', tcro_molecule
                write(*,*) 'abscro_molecule is                        ', abscro_molecule
         ! Now, we will get mass_ref (0.1cm^2)
            mass_ref= 0.1d0/(tcro_molecule*barn*(1.d0/mass_molecule)*(afj))
         ! Now, we will get the total absorption
            sigma_abs_total=(abscro_molecule)*barn*(1.d0/mass_molecule)*(afj)*mass_ref
          write(*,*) 'Recommend mass of every 1 cm^2            ', mass_ref
          write(*,*) 'Total absobtion of every 1cm^2            ', sigma_abs_total
  end subroutine

end module sqe_calculator
