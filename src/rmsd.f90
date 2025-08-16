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

! Subroutines used to calculate MSD and Debye-Waller factor

module rmsd

  implicit none

contains
  
  ! MSD calculator
  subroutine rmsdcalc()
    
    use phonon_spectra, only: DMsolver
    use func, only: fBE
    use constants, only: hbar, ps2s, amu, a2m, thz2mev, eps4
    use variables, only: qmesh, natoms, nbands, ntemps, ntypes, nat, &
                         temps, armsd, elements, rlatvec, masses2, &
                         write_rmsd, read_rmsd, filename_rmsd
    implicit none
    
    integer(kind=4) :: ntot, hh, ii, jj, kk, ll
    real(kind=8) :: T
    real(kind=8), allocatable :: qzone(:,:), omegas(:,:), amsd(:,:,:)
    complex(kind=8), allocatable :: eigenvecs(:,:,:)
    character(len=3) :: aux
    
    allocate(armsd(3,natoms,ntemps))
    
    if (read_rmsd) then
       write(*,*) "Read RMSD from file ..."
       write(*,*)
       write(*,*) "RMSD at each temperature is:"
       write(*,*)
       open(1, file=filename_rmsd, status="old")
       do ii = 1, ntemps
          read(1,*) T
          write(*,"(A4,x,F5.1,x,A1)") "T =",T, "K"
          do jj = 1, natoms
             read(1,*) aux, armsd(:,jj,ii)
             write(*,"(A4,x,3(F16.8,x))") aux, armsd(:,jj,ii)
          end do
          write(*,*)
       end do
       close(1)
       return
    end if
    
    write(*,*) "Start to compute RMSD ..."
    write(*,*)
  
    ntot = qmesh(1)*qmesh(2)*qmesh(3)
    allocate(qzone(3,ntot), omegas(ntot,nbands), eigenvecs(ntot,nbands,nbands))
    allocate(amsd(3,natoms,ntemps))
    omegas = 0.d0
    eigenvecs = 0.d0
    amsd = 0.d0
    
    
    ! discretize the q-points along Q1, Q2 and Q3
    ! Gamma centered mesh
    ll = 0
    do ii = 1, qmesh(1)
       do jj = 1, qmesh(2)
          do kk = 1, qmesh(3)
             ll = ll+1
             qzone(:,ll) = rlatvec(:,1)*(ii-1.d0)/qmesh(1)+&
                           rlatvec(:,2)*(jj-1.d0)/qmesh(2)+&
                           rlatvec(:,3)*(kk-1.d0)/qmesh(3)
          end do
       end do
    end do
    
    ! calculate phonon band structures
    call DMsolver(qzone, omegas, eigenvecs)

    ! set the all imaginary frequency to zero
    !where (omegas .lt. 0.d0) omegas = 0.d0
    !where(omegas < eps3) omegas = omegas+eps3

    ! MSD calculation
    do hh = 1, ntemps
       T = temps(hh)
       do ii = 1, natoms
          do jj = 1, 3
             do kk = 1, nbands
                do ll = 1, ntot
                   if (omegas(ll,kk) .lt. eps4) cycle
                   amsd(jj,ii,hh) = amsd(jj,ii,hh)+hbar/ntot/masses2(ii)/&
                                    omegas(ll,kk)*(0.5d0+fBE(omegas(ll,kk),T))*&
                                    abs(eigenvecs(ll,kk,((ii-1)*3+jj)))**2
                end do
             end do
          end do
       end do
    end do

    amsd = amsd*ps2s/amu/(a2m**2)
    armsd = sqrt(amsd)

    write(*,*) "RMSD calculation finished."
    write(*,*)

    ! write MSD into file
    write(*,*) "Write RMSD into file ..."
    write(*,*)
    write(*,*) "RMSD at each temperature is:"
    write(*,*)
    if (write_rmsd) then
       open(1, file=filename_rmsd, status="replace")
       do ii = 1, ntemps
          write(1,"(F7.2)") temps(ii)
          write(*,"(A4,x,F5.1,x,A1)") "T =", temps(ii), "K"
          ll = 0
          do jj = 1, ntypes
             do kk = 1, nat(jj)
                ll = ll+1
                write(1,"(A4,x,3(F16.8,x))") elements(jj), armsd(:,ll,ii)
                write(*,"(A4,x,3(F16.8,x))") elements(jj), armsd(:,ll,ii)
             end do
          end do
          write(*,*)
       end do
       close(1)
    end if

    deallocate(qzone, omegas, eigenvecs, amsd)
    
  end subroutine rmsdcalc

end module rmsd
