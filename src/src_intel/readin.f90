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

! Input parser and subrouines that read second-order force constants

module readin

  implicit none

contains
  
  ! read input file
  subroutine input_parser()
    
    use variables
    use constants, only: periodic_table, atomic_masses, scatt_b, tpi, eps5,scatt_totoal_cross_section,scatt_absob_cross_section
    use func, only: cross_prod
    use iso_fortran_env, only: error_unit
    implicit none
    
    integer(kind=4) :: ii, jj, kk, istat
    integer(kind=4), allocatable :: zid(:)
    character(len=10) :: coord_type
    character(len=20) :: inputfile
    character(len=100) :: buffer
    
    namelist / basic / ntypes, natoms, nsize
    namelist / inputsqe / &
        nat, ne, deltaE, nqh, nqk, nql, deltaH, deltaK, deltaL, q0, &
        masses, coh_b, xray_b, temp, temp_min, temp_max, temp_step, qmesh, &
        path, elements, nonanalytic, lneutron, lxray, write_rmsd, order, &
        xm, paras, functype, lresfunc, degauss, read_rmsd, filename_rmsd, l4d, &
        aff_wk, aff_a, aff_b, aff_c, espresso, clatvec, lphase, ltds, filename_omega, scatt_xs, abs_xs
    
    ! read the basic namelist
    ! ntypes, natoms, nsize, nat, path, elements, q0 must be specified in input file
    call get_command_argument(1, inputfile)
    open(1, file=inputfile, status="old")
    read(1, nml=basic)
    if (ntypes .lt. 1 .or. natoms .lt. 1 .or. natoms .lt. ntypes) then
       write(error_unit, *) "Error in input_parser: ntypes, natoms must be >= 1, natoms must be >= nelements"
       stop
    end if
    if(any(nsize .lt. 1)) then
       write(error_unit, *) "Error in input_parser: all components of nsize must be must be >= 1."
       stop
    end if
    allocate(elements(ntypes), nat(ntypes), positions(3,natoms), born(3,3,natoms), &
             masses(ntypes), coh_b(ntypes), xray_b(ntypes), masses2(natoms), &
             coh_b2(natoms), xray_b2(natoms), aff_a(5,ntypes), aff_b(5,ntypes), aff_c(ntypes),scatt_xs(ntypes),abs_xs(ntypes),scatt_xs2(natoms),abs_xs2(natoms))

    ! set the defaults
    ne = 1000
    deltaE = 0.01
    nqh = 100
    nqk = 6
    nql = 6
    deltaH = 0.01
    deltaK = 0.005
    deltaL = 0.005
    temp = 300.d0
    temp_min = 100.d0
    temp_max = 1000.d0
    temp_step = 0.d0
    qmesh = (/20, 20, 20/)
    path = 0.d0
    clatvec = 0.d0
    nonanalytic = .FALSE.
    espresso = .FALSE.
    !asr = .FALSE.
    !bornsym = .FALSE.
    lneutron = .TRUE.
    lxray = .TRUE.
    ltds = .FALSE.
    l4d = .FALSE.
    write_rmsd = .TRUE.
    read_rmsd = .FALSE.
    filename_rmsd = "rmsd.dat"
    filename_omega = "omega.dat"
    masses = 0.d0
    coh_b = 0.d0
    xray_b = 0.d0
    scatt_xs=0.d0
    abs_xs=0.d0
    epsilon = 0.d0
    epsilon(1,1) = 1.d0
    epsilon(2,2) = 1.d0
    epsilon(3,3) = 1.d0
    born = 0.d0
    aff_wk = .FALSE.
    aff_a = 0.d0
    aff_b = 0.d0
    aff_c = 0.d0
    degauss = 0.5
    order = 4
    xm = 1
    paras = 1.d0
    functype = "CNCS12"
    lresfunc = .FALSE.
    lphase = .FALSE.
    
    ! read the inputsqe namelist and crystal structure
    read(1, nml=inputsqe)
    ! to make the highest energy value reach ne*dE
    ne = ne+1
    
    if (temp_step .gt. 0.) then
       if (temp_min .gt. temp_max) then
          write(error_unit, *) "Error in input_parser: temp_max must be larger than temp_min."
          stop
       else
          ntemps = floor((temp_max-temp_min)/temp_step)+1
       end if
    else
       ntemps = 1
    end if
    allocate(temps(ntemps))
    if (ntemps .gt. 1) then
       do ii = 1, ntemps
          temps(ii) = temp_min+temp_step*(ii-1)
       end do
    else
       temps(1) = temp
    end if

    if (lresfunc) then
       if (.not. (functype .eq. "CNCS12" &
           .or. functype .eq. "CNCS20" &
           .or. functype .eq. "HNTAS" &
           .or. functype .eq. "CNTAS" &
           .or. functype .eq. "poly")) then
          write(error_unit, *) "Error in input_parser: unknown resolution function type."
       stop
       end if
    end if
    
    do while (.TRUE.)
       read(1, *, iostat=istat) buffer
       if (istat .ne. 0) exit
       if (index(buffer, "LATTICE_PARAMETERS") .ne. 0) then
          read(1, *, iostat=istat) lscale
          do ii = 1, 3
             read(1, *, iostat=istat) latvec(:,ii)
             latvec(:,ii) = latvec(:,ii)*lscale
          end do
          if (all(clatvec .lt. eps5)) clatvec = latvec
       end if
       if (index(buffer, "ATOMIC_POSITIONS") .ne. 0) then
          read(1, *, iostat=istat) coord_type
          coord_type = adjustl(coord_type)
          do ii = 1, natoms
             read(1, *, iostat=istat) positions(:,ii)
          end do
       end if
       if (index(buffer, "BORN") .ne. 0) then
          do ii = 1, 3
             read(1, *, iostat=istat) epsilon(:,ii)
          end do
          do ii = 1, natoms
             do jj = 1, 3
                read(1, *, iostat=istat) born(:,jj,ii)
             end do
          end do
       end if
    end do
    close(1)

    ! convert lattice basis to cartesian basis if coord_type(1) is 'D' or 'd'
    if (coord_type(1:1) == 'D' .or. coord_type(1:1) == 'd') then
       positions = matmul(latvec, positions)
    elseif (coord_type(1:1) .eq. 'C' .or. coord_type(1:1) .eq. 'c') then
       positions = positions
    else
       write(error_unit, *) "Error in input_parser: unknown atomic coordinate type."
       stop
    end if
    
    ! calculate the reciprocal lattice and volume of unit cell
    do ii = 1, 3
       jj = mod(ii,3)+1
       kk = mod(jj,3)+1
       call cross_prod(latvec(:,jj), latvec(:,kk), rlatvec(:,ii))
       call cross_prod(clatvec(:,jj), clatvec(:,kk), crlatvec(:,ii))
    end do
    vol = abs(dot_product(latvec(:,1), rlatvec(:,1)))
    cvol = abs(dot_product(clatvec(:,1), crlatvec(:,1)))
    rlatvec = tpi*rlatvec/vol
    crlatvec = tpi*crlatvec/cvol

    ! set atomic masses and scattering lengths
    allocate(zid(ntypes))
    do ii = 1, ntypes
       ! find the positions in periodic table
       jj = 1
       do while(jj .le. size(periodic_table))
          if (elements(ii) .eq. periodic_table(jj)) then
             zid(ii) = jj
             exit
          end if
          jj = jj+1
       end do
       if (jj .eq. (size(periodic_table)+1)) then
          write(error_unit, *) "Error in input_parser: specified elements not in current SQE program,&
                               & please manually set the corresponding masses, coh_b, xray_b,scatt_xs,abs_xs."
          stop
       end if
    end do
    if (any(masses .lt. eps5)) then
       do ii = 1, ntypes
          masses(ii) = atomic_masses(zid(ii))
       end do
    end if
    if (all(coh_b .lt. eps5) .and. lneutron) then
       do ii = 1, ntypes
          coh_b(ii) = scatt_b(zid(ii))
       end do
    end if
    if (any(xray_b .lt. eps5) .and. lxray .and. (.not. aff_wk)) then
       do ii = 1, ntypes
          xray_b(ii) = dble(zid(ii))
       end do
    end if
    !get total cross section
    if (all(scatt_xs .lt. eps5).and. lneutron) then
       do ii = 1, ntypes
          scatt_xs(ii) = scatt_totoal_cross_section(zid(ii))
       end do
    end if
    !get absorb cross section
    if (all(abs_xs .lt. eps5) .and. lneutron) then
       do ii = 1, ntypes
          abs_xs(ii) = scatt_absob_cross_section(zid(ii))
       end do
    end if
    deallocate(zid)
    
    kk = 0
    do ii = 1, ntypes
       do jj = 1, nat(ii)
          kk = kk+1
          masses2(kk) = masses(ii)
          coh_b2(kk) = coh_b(ii)
          xray_b2(kk) = xray_b(ii)
          scatt_xs2(kk)=scatt_xs(ii)
          abs_xs2(kk)=abs_xs(ii)
       end do
    end do

    deallocate(masses, coh_b, xray_b, scatt_xs,abs_xs)
  
  end subroutine input_parser

  ! scattering length printer
  subroutine cohb_print

    use constants, only: periodic_table, scatt_b
    use iso_fortran_env, only: error_unit
    implicit none

    character(len=3) :: xname
    integer(kind=4) :: ii, zid
    real(kind=8) :: coh_b, xray_b

    call get_command_argument(2, xname)

    ! find the positions in periodic table
    ii = 1
    do while(ii .le. size(periodic_table))
       if (xname .eq. periodic_table(ii)) then
          zid = ii
          exit
       end if
       ii = ii+1
    end do
    if (ii .eq. (size(periodic_table)+1)) then
       write(error_unit, *) "Specified elements not in current SQE program,&
                             & please manually set the corresponding masses, coh_b, xray_b."
       stop
    end if

    coh_b = scatt_b(zid)
    xray_b = dble(zid)
    write(*,"(A4,x,F8.3)") "INS:", coh_b
    write(*,"(A4,x,F8.3)") "IXS:", xray_b
    stop
  
  end subroutine cohb_print

  ! Free the memory used by input
  subroutine sqe_free()
    
    use variables, only: elements, nat, positions, born, masses2, sqesum, &
                         qlist, coh_b2, xray_b2, temps, cqlist, eigenvec, &
                         omega, fc_short, aff_a, aff_b, aff_c
    implicit none
    
    deallocate(elements, nat, positions, born, masses2, sqesum, &
               qlist, coh_b2, xray_b2, temps, cqlist, eigenvec, &
               omega, fc_short, aff_a, aff_b, aff_c)

  end subroutine sqe_free
  
  ! read FORCE_CONSTANTS file, Phonopy
  subroutine read2fc()
  
    use variables, only: natoms, nsize, masses2, fc_short
    use iso_fortran_env, only: error_unit
    use constants, only: dm2thz
    implicit none
    
    integer(kind=4) :: ntot, atom1, atom2, i, j, ip
    integer(kind=4) :: ix1, iy1, iz1, ix2, iy2, iz2, iatom1, iatom2
    real(kind=8) :: mm(natoms,natoms)

    do i = 1, natoms
       mm(i,i) = masses2(i)
       do j = i+1, natoms
          mm(i,j) = sqrt(masses2(i)*masses2(j))
          mm(j,i) = mm(i,j)
       end do
    end do

    allocate(fc_short(natoms, 3, nsize(1), nsize(2), nsize(3), natoms, 3))

    ! Phonopy's 2nd-order format is quite straightforward. Each file
    ! is essentially a sequence of 3x3 blocks, one for each pair of atoms
    ! in the supercell. A single header line contains the number of atoms.
    open(1, file="FORCE_CONSTANTS", status="old")
    read(1, *) ntot
    if (ntot .ne. nsize(1)*nsize(2)*nsize(3)*natoms) then
       write(error_unit,*) "Error in read2fc: wrong number of force constants for the specified supercell."
       stop
    end if
    do i=1, ntot
       do j=1, ntot
          read(1,*) atom1, atom2
          call split_index(atom1, nsize(1), nsize(2), nsize(3), ix1, iy1, iz1, iatom1)
          call split_index(atom2, nsize(1), nsize(2), nsize(3), ix2, iy2, iz2, iatom2)
          if(ix1.eq.1 .and. iy1.eq.1 .and. iz1.eq.1) then
             do ip = 1, 3
                read(1, *) fc_short(iatom1, ip, ix2, iy2, iz2, iatom2, :)
             end do
          else
             do ip = 1, 3
                read(1, *)
             end do
          end if
       end do
    end do
    close(1)

    ! After reading the force constants, they are reduced using the
    ! tensor product of the square root of the atomic masses. It is these
    ! reduced constants that enter the expression of the dynamical matrix.
    do iatom1 = 1, natoms
       do iatom2 = 1, natoms
          fc_short(iatom1, :, :, :, :, iatom2, :) = fc_short(iatom1, :, :, :, :, iatom2, :)/mm(iatom1, iatom2)
       end do
    end do
    fc_short = dm2thz*fc_short

  end subroutine read2fc

  ! read espresso.ifc2 file, Quantum Espresso
  subroutine read2fc_espresso()

    use variables, only: ntypes, natoms, nsize, masses2, fc_short
    use iso_fortran_env, only: error_unit
    use constants, only: dm2thz, ryd2ev, bohr2ang
    implicit none
    
    integer(kind=4) :: i, j, ipol, jpol, iat, jat, t1, t2, t3
    integer(kind=4) :: ntype, nat, ibrav, tipo(natoms), qscell(3)
    real(kind=8) :: mm(natoms,natoms), mass(ntypes), celldm(6), at(3,3), r(natoms,3)
    real(kind=8) :: eps(3,3), zeff(natoms,3,3)
    character(len=1) :: polar_key
    character(len=5) :: label(ntypes)

    do i = 1, natoms
       mm(i,i) = masses2(i)
       do j = i+1, natoms
          mm(i,j) = sqrt(masses2(i)*masses2(j))
          mm(j,i) = mm(i,j)
       end do
    end do

    allocate(fc_short(natoms, 3, nsize(1), nsize(2), nsize(3), natoms, 3))

    open(1, file="espresso.fc", status="old")
    read(1, *) ntype, nat, ibrav, celldm(1:6)
    if (nat .ne. natoms) then
       write(error_unit,*) "Error in read2fc_espresso: inconsistent number of atoms in unit cell."
       stop
    end if
    if (ibrav == 0) then
       read(1, *) ((at(i,j), i=1,3), j=1,3)
    end if

    do i = 1, ntypes
       read(1, *) j, label(i), mass(i)
    end do
    do i = 1, natoms
       read(1, *) j, tipo(i), r(i,1:3)
    end do

    read(1, *) polar_key
    if (polar_key .eq. "T") then
       do i = 1, 3
          read(1, *) eps(i,1:3)
       end do
       do i = 1, natoms
          read(1, *)
          do j = 1, 3
             read(1, *) zeff(i,j,1:3)
          end do
       end do
    end if
    read(1, *) qscell(1:3)
    if (qscell(1)*qscell(2)*qscell(3) .ne. nsize(1)*nsize(2)*nsize(3)) then
       write(error_unit,*) "Error in read2fc_espresso: inconsistent dimension of supercell."
       stop
    end if
    ! Read the force constants.
    do i = 1, 3*3*natoms*natoms
       read(1, *) ipol, jpol, iat, jat
       do j = 1, nsize(1)*nsize(2)*nsize(3)
          read(1, *) t1, t2, t3, fc_short(jat,jpol,t1,t2,t3,iat,ipol)
       end do
    end do
    close(1)

    ! convert from Ry/bohr^2 to eV/A^2
    fc_short = fc_short*ryd2ev/bohr2ang/bohr2ang

    ! After reading the force constants, they are reduced using the
    ! tensor product of the square root of the atomic masses. It is these
    ! reduced constants that enter the expression of the dynamical matrix.
    do iat = 1, natoms
       do jat = 1, natoms
          fc_short(iat, :, :, :, :, jat, :) = fc_short(iat, :, :, :, :, jat, :)/mm(iat, jat)
       end do
    end do
    fc_short = dm2thz*fc_short

  end subroutine read2fc_espresso
  
  ! Convert a supercell index of the kind used by Phonopy into
  ! a set of unit cell+atom indices.
  subroutine split_index(idx, nx, ny, nz, ix, iy, iz, iatom)
    
    use func, only: divmod
    implicit none

    integer(kind=4), intent(in) :: idx, nx, ny, nz
    integer(kind=4), intent(out) :: ix, iy, iz, iatom

    integer(kind=4) :: tmp1, tmp2

    call divmod(idx-1, nx, tmp1, ix)
    call divmod(tmp1, ny, tmp2, iy)
    call divmod(tmp2, nz, iatom, iz)

    ix = ix+1
    iy = iy+1
    iz = iz+1
    iatom = iatom+1
    
  end subroutine split_index

end module readin
