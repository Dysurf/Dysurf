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
  use func
  use sqe_calculator
  implicit none

  ! if any INCLUDE

  ! record cpu running time
  integer(kind=4) :: ii, jj, nargs
  integer(kind=4) :: minute, second
  integer(kind=4) :: time_begin, time_end, time_last
  character(len=20) :: sqe_file, string
  character(len=50) :: fmtstring

  ! record mass_ref and abstot_with_this_mass for experiment
  real(kind=8) ::mass_ref,sigma_abs_total

  ! start
  nargs = command_argument_count()
  if (nargs .eq. 2) call cohb_print()
  call system_clock(time_begin)
  write(*,*) "----------------------------------"
  write(*,*) "     Dysurf Program Version 1.1   "
  write(*,*) "----------------------------------"
  write(*,*)
  write(*,*) "Dysurf program begins ..."
  write(*,*)

  ! read the input
  call input_parser()
  write(*,*) "Read the input successfully."
  write(*,*)

  ! read second-order force constants and weight with atomic masses
  if (espresso) then
     call read2fc_espresso()
     write(*,*) "Read the second-order force constants from espresso.fc successfully."
     write(*,*)
  else
     call read2fc()
     write(*,*) "Read the second-order force constants from FORCE_CONSTANTS successfully."
     write(*,*)
  end if

  ! Generate q-point list
  call load_qpoints()

  ! write full q-points of bin box into file
  write(*,*) "Write full q-points into file ..."
  write(*,*)
  open(1, file="qpoints_full.dat", status="replace")
  do ii = 1, nqtot
     write(1, "(I6,x,3(F16.8,x))") ii, pqlist(:,ii)
  end do
  close(1)

  ! phonon spectra
  write(*,*) "Start to compute phonon spectra ..."
  write(*,*)
  call phoneigen()
  write(*,*) "Phonon spectra calculation finished."
  write(*,*)

  ! write q-points of single dispersion path into file
  if (.not.ltds .and. .not.l4d) then
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

  ! RMSD
  call rmsdcalc()

  ! Dynamical structure factor
  write(*,*) "Start to compute S(Q,E) ..."
  write(*,*)
  if (lneutron) then
     write(*,*) "Enter inelastic neutron scattering ..."
     write(*,*)
     call sqecalc()
  elseif (lxray) then
     write(*,*) "Enter inelastic X-ray scattering ..."
     write(*,*)
     call sqecalc()
  end if
  write(*,*) "S(Q,E) calculation finished."
  write(*,*)

  ! write the output
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
        ! fmtstring = "("//trim(adjustl(string))//"*(E24.16,x))"
        fmtstring = "("//trim(adjustl(string))//"(E24.16,1x))"
        open(1, file=sqe_file, status="replace")
        do jj = 1, nqh
           write(1, fmtstring) sqesum(:,jj,ii)
        end do
        close(1)
     end do
  end if

  ! give reference mass for experiment
  if (lneutron) then
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
  time_last = (time_end-time_begin)/10000+1
  minute = time_last/60
  second = time_last-minute*60
  write(*,*) "The program lasts", minute, "min", second, "sec"
  write(*,*)

end program DYSURF
