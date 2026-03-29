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

! Dysurf main program
! Only one phonon emission process is calculated

program DYSURF
  
  use readin
  use variables
  use constants
  use qpoints
  use phonon_spectra
  use rmsd
  use mpi_helper
  use mpi
  use func
  use sqe_calculator
  use psf_2d_state
  implicit none

  ! if any INCLUDE

  ! record cpu running time
  integer(kind=4) :: ii, jj, nargs
  integer(kind=4) :: ierr_mpi
  integer(kind=4) :: minute, second
  integer(kind=4) :: time_begin, time_end, time_last
  character(len=64) :: sqe_file, string, temp_string
  character(len=50) :: fmtstring

  ! record mass_ref and abstot_with_this_mass for experiment
  real(kind=8) ::mass_ref,sigma_abs_total
  real(kind=8) :: wall_begin, wall_end, wall_elapsed, wall_elapsed_max

  ! start
  call mpi_helper_init()
  nargs = command_argument_count()
  if (nargs .eq. 2) call cohb_print()
  call system_clock(time_begin)
  wall_begin = MPI_Wtime()
  if (mpi_rank == 0) then
     write(*,*) "----------------------------------"
     write(*,*) "     Dysurf Program Version 1.2   "
     write(*,*) "----------------------------------"
     write(*,*)
     write(*,*) "Dysurf program begins ..."
     write(*,*)
     if (mpi_size > 1) then
        write(*,*) "MPI enabled with", mpi_size, "ranks"
        write(*,*)
     end if
  end if

  ! read the input
  call input_parser()
  if (mpi_rank == 0) then
     write(*,*) "Read the input successfully."
     write(*,*)
  end if

  ! read second-order force constants and weight with atomic masses
  if (espresso) then
     call read2fc_espresso()
     if (mpi_rank == 0) then
        write(*,*) "Read the second-order force constants from espresso.fc successfully."
        write(*,*)
     end if
  else
     call read2fc()
     if (mpi_rank == 0) then
        write(*,*) "Read the second-order force constants from FORCE_CONSTANTS successfully."
        write(*,*)
     end if
  end if

  ! Generate q-point list
  call load_qpoints()

  ! write full q-points of bin box into file
  if (mpi_rank == 0) then
     write(*,*) "Write full q-points into file ..."
     write(*,*)
     open(1, file="qpoints_full.dat", status="replace")
     do ii = 1, nqtot
        write(1, "(I6,x,3(F16.8,x))") ii, pqlist(:,ii)
     end do
     close(1)
  end if

  ! phonon spectra
  if (mpi_rank == 0) then
     write(*,*) "Start to compute phonon spectra ..."
     write(*,*)
  end if
  call phoneigen()
  if (mpi_rank == 0) then
     write(*,*) "Phonon spectra calculation finished."
     write(*,*)
  end if

  ! write q-points of single dispersion path into file
  if (.not.ltds .and. .not.l4d) then
     if (mpi_rank == 0) then
        write(*,*) "Write q-points into file ..."
        write(*,*)
        open(1, file="qpoints.dat", status="replace")
        do ii = 1, nqh
           write(1, "(I6,x,3(F16.8,x))") ii, pqlist(:,elist(ii))
        end do
        close(1)

        ! write phonon dispersion into file
        write(*,*) "Write phonon energy into file ..."
        write(*,*)
        open(1, file=filename_omega, status="replace")
        write(string, *) nbands+1
        fmtstring = "("//trim(adjustl(string))//"E20.10)"
        do ii = 1, nqh
           write(1, fmtstring) qdist(ii), omega(elist(ii),:)
        end do
        close(1)
     end if
  end if

  ! RMSD
  call rmsdcalc()

  ! Dynamical structure factor
  if (mpi_rank == 0) then
     write(*,*) "Start to compute S(Q,E) ..."
     write(*,*)
  end if
  if (lneutron) then
     if (mpi_rank == 0) then
        write(*,*) "Enter inelastic neutron scattering ..."
        write(*,*)
     end if
     call sqecalc()
  elseif (lxray) then
     if (mpi_rank == 0) then
        write(*,*) "Enter inelastic X-ray scattering ..."
        write(*,*)
     end if
     call sqecalc()
  end if
  if (mpi_rank == 0) then
     write(*,*) "S(Q,E) calculation finished."
     write(*,*)
  end if

  ! write the output
  if (mpi_rank == 0) then
     write(*,*) "Write S(Q,E) into file ..."
     write(*,*)
     if (l4d) then
        do ii = 1, ntemps
           write(string, *) int(temps(ii))
           sqe_file = "SQEBIN_"//trim(adjustl(string))//"K"
           sqe_file = trim(sqe_file)//".dat"
           write(string, *) nqk*nql*ne
           fmtstring = "("//trim(adjustl(string))//"*(E24.16,x))"
           open(1, file=sqe_file, status="replace")
           do jj = 1, nqh
              write(1, fmtstring) reshape(sqebin(:,:,:,jj,ii),(/nqk*nql*ne/))
           end do
           close(1)
        end do
     else
        do ii = 1, ntemps
           write(string, *) int(temps(ii))
           sqe_file = "SQE_"//trim(adjustl(string))//"K"
           sqe_file = trim(sqe_file)//".dat"
           if (ltds) then
              write(string, *) nqk
           else
              write(string, *) ne
           end if
           fmtstring = "("//trim(adjustl(string))//"(E24.16,1x))"
           open(1, file=sqe_file, status="replace")
           do jj = 1, nqh
              write(1, fmtstring) sqesum(:,jj,ii)
           end do
           close(1)
        end do
     end if
     if (lresfunc2d_psf) then
        do ii = 1, ntemps
           write(string, *) int(temps(ii))
           if (lsave_intrinsic_slice .and. allocated(intrinsic_slice)) then
              temp_string = string
              sqe_file = "SQE_intrinsic_"//trim(adjustl(temp_string))//"K.dat"
              write(string, *) ne
              fmtstring = "("//trim(adjustl(string))//"(E24.16,1x))"
              open(1, file=sqe_file, status="replace")
              do jj = 1, nqh
                 write(1, fmtstring) intrinsic_slice(:,jj,ii)
              end do
              close(1)
              sqe_file = "line_intrinsic_qmid_"//trim(adjustl(temp_string))//"K.dat"
              open(1, file=sqe_file, status="replace")
              do jj = 1, ne
                 write(1, "(2(E24.16,1x))") deltaE * dble(jj - 1), intrinsic_slice(jj, max(1, (nqh + 1) / 2), ii)
              end do
              close(1)
           end if
           if (lsave_psf_slice .and. allocated(psf2d_slice)) then
              write(temp_string, *) int(temps(ii))
              sqe_file = "SQE_psf2d_"//trim(adjustl(temp_string))//"K.dat"
              write(string, *) ne
              fmtstring = "("//trim(adjustl(string))//"(E24.16,1x))"
              open(1, file=sqe_file, status="replace")
              do jj = 1, nqh
                 write(1, fmtstring) psf2d_slice(:,jj,ii)
              end do
              close(1)
              sqe_file = "line_psf2d_qmid_"//trim(adjustl(temp_string))//"K.dat"
              open(1, file=sqe_file, status="replace")
              do jj = 1, ne
                 write(1, "(2(E24.16,1x))") deltaE * dble(jj - 1), psf2d_slice(jj, max(1, (nqh + 1) / 2), ii)
              end do
              close(1)
           end if
           if (lsave_psf_params .and. allocated(sigma_q_field)) then
              write(string, *) ne
              fmtstring = "("//trim(adjustl(string))//"(E24.16,1x))"
              write(temp_string, *) int(temps(ii))

              sqe_file = "PSF_sigma_q_"//trim(adjustl(temp_string))//"K.dat"
              open(1, file=sqe_file, status="replace")
              do jj = 1, nqh
                 write(1, fmtstring) sigma_q_field(:,jj,ii)
              end do
              close(1)

              sqe_file = "PSF_sigma_e_left_"//trim(adjustl(temp_string))//"K.dat"
              open(1, file=sqe_file, status="replace")
              do jj = 1, nqh
                 write(1, fmtstring) sigma_e_left_field(:,jj,ii)
              end do
              close(1)

              sqe_file = "PSF_sigma_e_right_"//trim(adjustl(temp_string))//"K.dat"
              open(1, file=sqe_file, status="replace")
              do jj = 1, nqh
                 write(1, fmtstring) sigma_e_right_field(:,jj,ii)
              end do
              close(1)

              sqe_file = "PSF_shear_"//trim(adjustl(temp_string))//"K.dat"
              open(1, file=sqe_file, status="replace")
              do jj = 1, nqh
                 write(1, fmtstring) shear_field(:,jj,ii)
              end do
              close(1)
           end if
        end do
     end if
     if (trim(tas_obs_mode) .eq. 'psf2d') then
        write(*,*) '2D PSF resolution: main SQE output written from psf2d-convolved slice.'
        if (.not. lresfunc2d_psf) then
           write(*,*) '2D PSF resolution: debug side-output files disabled because lresfunc2d_psf=.FALSE..'
        end if
        write(*,*)
     end if
  end if

  ! give reference mass for experiment
  if (lneutron .and. mpi_rank == 0) then
      ! *********************** Neutron Experiment Recommendations ***********************
      write(*,*) '***** Recommended Sample Mass for Experiment (unit: grams) *****'
    
      ! Calculate reference mass and total absorption cross-section
      call given_Mass_reference(mass_ref, sigma_abs_total)
    
      write(*,*)  ! Blank line for better readability
    
      ! Provide specific recommendation for ARCS instrument
      write(*,*) 'Instrument example: ARCS with 5x5 cm sample container'
      write(*,*) '----------------------------------------------------'
      write(*,*) 'Recommended mass for ARCS:          ', mass_ref * 25.0d0, ' g'
      write(*,*) 'Total absorption for this mass:     ', sigma_abs_total * 25.0d0, ' barns'
      write(*,*) '****************************************************************'
  endif

  ! end
  call sqe_free()
  call system_clock(time_end)
  wall_end = MPI_Wtime()
  wall_elapsed = wall_end - wall_begin
  wall_elapsed_max = wall_elapsed
  if (mpi_size > 1) then
     call MPI_Reduce(wall_elapsed, wall_elapsed_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr_mpi)
  end if
  time_last = (time_end-time_begin)/10000+1
  minute = time_last/60
  second = time_last-minute*60
  if (mpi_rank == 0) then
     minute = int(wall_elapsed_max) / 60
     second = int(wall_elapsed_max) - minute * 60
     write(*,*) "The program lasts", minute, "min", second, "sec"
     write(*,*)
  end if
  call mpi_helper_finalize()

end program DYSURF
