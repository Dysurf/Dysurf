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
! Variables in SQE program

module control

  ! common variables that control the SQE run

  integer(kind=4) :: ntypes
  ! no. of types of atoms in the system
  integer(kind=4) :: natoms
  ! totoal no. of atoms
  integer(kind=4), allocatable :: nat(:)
  ! no. of atoms in each type
  integer(kind=4) :: nsize(3)
  ! the size of supercell
  integer(kind=4) :: ne
  ! no. of points of phonon energy, Emax = ne*deltaE (meV)
  integer(kind=4) :: nqh
  ! no. of q-points along the calculated path
  integer(kind=4) :: nqk, nql
  ! no. of q-points perpendicular to the calculated path
  integer(kind=4) :: qmesh(3)
  ! q-mesh used in MSD calculation
  integer(kind=4) :: order, xm
  ! the order of fitting polynomial function for energy resolution
  ! multiple for standard deviation of energy
  
  real(kind=8) :: clatvec(3,3)
  ! lattice vector in angstrom for conventional unit cell
  real(kind=8), allocatable :: masses(:)
  ! atomic masses corresponding to each element
  real(kind=8), allocatable :: coh_b(:)
  ! coherent neutron scattering length corresponding to each element
  real(kind=8), allocatable :: scatt_xs(:)
  ! total neutron scattering cross section corresponding to each element
  real(kind=8), allocatable :: abs_xs(:)
  ! coherent neutron scattering length corresponding to each element
  real(kind=8), allocatable :: xray_b(:)
  ! scattering length of X-ray scattering corresponding to each element
  real(kind=8) :: temp, temp_min, temp_max, temp_step
  ! tempature to be used in calculations
  real(kind=8) :: deltaE
  ! phonon energy resolution in meV
  real(kind=8) :: deltaH, deltaK, deltaL
  ! wavevector resolution along and perpendicular to the calculated path in conventional cell basics
  real(kind=8) :: path(3,3)
  ! phonon dispersion path used in SQE calculation
  real(kind=8) :: u(3)
  ! sample orientation vector used only by the TAS PSF 4D-resolution/PSF path
  real(kind=8) :: v(3)
  ! second sample orientation vector used only by the TAS PSF 4D-resolution/PSF path
  real(kind=8) :: q0(3)
  ! orgins of q-path
  real(kind=8) :: degauss
  ! broadening in rad/ps of gaussian smearing for energy conservation
  real(kind=8) :: paras(-1:5)
  ! fitted parameters for stand deviation of energy
  ! polynomial function is implemented up to order 5
  real(kind=8), allocatable :: aff_a(:,:), aff_b(:,:), aff_c(:)
  ! parameters a1-a5, b1-b5, c for atomic form factor of X-ray
  ! real(kind=8) :: sigmaq
  ! stand deviation of q in resolution function
  ! e.g. 0.01, 0.025
  
  character(len=3), allocatable :: elements(:)
  ! name of each type of atom
  character(len=20) :: filename_rmsd
  ! filename of root mean square atomic displacements
  character(len=20) :: filename_omega
  ! filename of phonon energy
  character(len=10) :: functype
  ! type of resolution function
  ! it can be "CNCS12", "CNCS20" and "poly"
  character(len=16) :: tas_mode = 'none'
  ! TAS resolution mode, e.g. 'none' or 'CN'
  character(len=16) :: tas_fix_mode = 'Ef'
  ! TAS fixed-energy mode, e.g. 'Ef' or 'Ei'
  character(len=24) :: tas_obs_mode = 'linearized_4d'
  ! TAS observed-spectrum mode: linearized_4d or psf2d
  character(len=16) :: psf_mode = 'analytic'
  ! 2D PSF source mode, currently only analytic is implemented
  integer(kind=4) :: tas_mono_dir = -1
  ! TAS monochromator sign convention, matching NeutronPy mono.dir
  integer(kind=4) :: tas_sample_dir = 1
  ! TAS sample sign convention, matching NeutronPy sample.dir
  integer(kind=4) :: tas_ana_dir = -1
  ! TAS analyzer sign convention, matching NeutronPy ana.dir
  integer(kind=4) :: psf_window_q = 6
  ! local 2D PSF convolution half-window along q, in q bins
  integer(kind=4) :: psf_window_e = 12
  ! local 2D PSF convolution half-window along energy, in energy bins
  !character(len=20) :: filename_2fc = "FORCE_CONSTANTS"
  ! filename of second-order force constants generated in format of PHONOPY
  
  logical :: nonanalytic
  ! if .TRUE. compute and use the nonanalytic part of the dynamical matrix for LO-TO splitting
  logical :: lneutron
  ! if .TRUE. compute the dynamic structral factor of inelastic neutron scattering
  logical :: lxray
  ! if .TRUE. compute the dynamic structral factor of inelastic X-ray scattering
  logical ltds
  ! if .TRUE. compute thermal diffuse scattering
  logical l4d
  ! if .TRUE. compute 4D S(Q,E) and SQEBIN would be written into file
  logical :: aff_wk
  ! if .TRUE. use atomic form factor for inelastic X-ray scattering
  logical :: write_rmsd
  ! if .TRUE. write RMSD into file
  logical :: read_rmsd
  ! if .TRUE. read RMSD from file
  logical :: lresfunc
  ! if .TRUE. consider the Gaussian instrument resolution function
  logical :: espresso
  ! if .TRUE. read the second-order force constants in QE DFPT format
  logical :: lphase
  ! if .TRUE. consider phase from the dot product of qtransfer with atom position in cell
  logical :: tas_use_cn = .FALSE.
  ! explicit switch for Cooper-Nathans TAS resolution
  logical :: tas_use_mosaic = .TRUE.
  ! if .TRUE. include mosaic terms in TAS resolution
  logical :: tas_debug_res = .FALSE.
  ! if .TRUE. print extra information for TAS resolution debugging
  logical :: lresfunc2d_psf = .FALSE.
  ! if .TRUE. write additional debug outputs for the position-dependent 2D PSF path
  logical :: lsave_intrinsic_slice = .FALSE.
  ! if .TRUE. write the intrinsic/model slice used as 2D PSF input
  logical :: lsave_psf_slice = .FALSE.
  ! if .TRUE. write the 2D PSF-convolved slice
  logical :: lsave_psf_params = .FALSE.
  ! if .TRUE. write the local 2D PSF parameter fields
  !logical :: asr
  ! if .TRUE. perform acoustic sum rule for force constants
  !logical :: bornsym
  ! if .TRUE. csymmetrize the Born effective charges

  real(kind=8) :: tas_e_fixed = 14.7d0
  ! fixed incident or final energy in meV
  real(kind=8) :: tas_coll_h_pre_mono = 40.d0
  ! horizontal collimation before monochromator in arcmin
  real(kind=8) :: tas_coll_h_pre_samp = 40.d0
  ! horizontal collimation between monochromator and sample in arcmin
  real(kind=8) :: tas_coll_h_post_samp = 40.d0
  ! horizontal collimation between sample and analyzer in arcmin
  real(kind=8) :: tas_coll_h_post_ana = 40.d0
  ! horizontal collimation after analyzer in arcmin
  real(kind=8) :: tas_mosaic_mono_h = 30.d0
  ! monochromator horizontal mosaic in arcmin
  real(kind=8) :: tas_mosaic_ana_h = 30.d0
  ! analyzer horizontal mosaic in arcmin
  real(kind=8) :: tas_mosaic_samp_h = 20.d0
  ! sample horizontal mosaic in arcmin
  real(kind=8) :: tas_dm = 3.355d0
  ! monochromator d-spacing in angstrom
  real(kind=8) :: tas_da = 3.355d0
  ! analyzer d-spacing in angstrom
  real(kind=8) :: psf_sigma_q_default = 1.d0
  ! default q-direction width for 2D PSF fallback
  real(kind=8) :: psf_sigma_e_left_default = 1.5d0
  ! default left-side energy width for 2D PSF fallback
  real(kind=8) :: psf_sigma_e_right_default = 1.5d0
  ! default right-side energy width for 2D PSF fallback
  real(kind=8) :: psf_shear_default = 0.d0
  ! default shear for 2D PSF fallback
end module control

module phon

  ! Variables associated with SQE calculation

  integer(kind=4) :: nbands
  ! no. of bands in phonon spectra
  integer(kind=4) :: nqtot
  ! total no. of q-points
  integer(kind=4) :: ntemps
  ! no of temperatures in calculation
  integer(kind=4), allocatable :: elist(:)
  ! store qpoint id of phonon dispersion
  !integer(kind=4), allocatable :: zid(:)
  ! atomic number of each type

  real(kind=8) :: lscale
  ! scale factor of lattice parameter
  real(kind=8) :: latvec(3,3)
  ! lattice vectors in angstrom for primitive unit cell
  real(kind=8) :: vol
  ! volume of primitive unit cell
  real(kind=8) :: cvol
  ! volume of conventional unit cell
  real(kind=8) :: rlatvec(3,3)
  ! reciprocal lattice vectors for primitive unit cell
  real(kind=8) :: crlatvec(3,3)
  ! reciprocal lattice vectors for conventional unit cell
  real(kind=8), allocatable :: positions(:,:)
  ! atomic positions in Cartesian coordinates
  real(kind=8), allocatable :: temps(:)
  ! a list of temperatures to be used in calculations
  real(kind=8) :: epsilon(3,3)
  ! dielectric tensor of the system in the Cartesian basis
  real(kind=8), allocatable :: born(:,:,:)
  ! Born effective charge tensor of each atom in the Cartesian basis
  real(kind=8), allocatable :: masses2(:)
  ! atomic masses corresponding to each atom
  real(kind=8), allocatable :: coh_b2(:)
  ! coherent neutron scattering length corresponding to each atom
  real(kind=8), allocatable :: scatt_xs2(:)
  ! total neutron scattering cross section corresponding to each element
  real(kind=8), allocatable :: abs_xs2(:)
  ! coherent neutron scattering length corresponding to each atom
  real(kind=8), allocatable :: xray_b2(:)
  ! scattering length of X-ray scattering corresponding to each atom
  !real(kind=8), allocatable :: distances(:,:)
  ! measured distance from q0 center
  real(kind=8), allocatable :: qlist(:,:)
  ! a list of q-points in conventional cell basis
  ! dimension(3, nqtot)
  real(kind=8), allocatable :: pqlist(:,:)
  ! a list of q-points in primitive cell basis
  ! dimension(3, nqtot)
  real(kind=8), allocatable :: cqlist(:,:)
  ! a list of q-points in Cartesian basis
  ! dimension(3, nqtot)
  real(kind=8), allocatable :: qdist(:)
  ! a list of the distance along dispersion path
  ! dimension(nqtot)
  real(kind=8), allocatable :: omega(:,:)
  ! phonon frequency, dimension(nqtot,nbands)
  real(kind=8), allocatable :: fc_short(:,:,:,:,:,:,:)
  ! short-range force constants
  ! dimension(natoms, 3, nsize(1), nsize(2), nsize(3), natoms, 3)
  real(kind=8), allocatable :: armsd(:,:,:)
  ! root mean square displacement for each atom at each temperature
  !real(kind=8), allocatable :: DWfac(:,:)
  ! Debye_waller factor for each atom at each temperature
  ! norm square of q-points is not included
  real(kind=8), allocatable :: sqebin(:,:,:,:,:)
  ! SQE of each bin box, dimension(ne,nql,nqk,,nqh,ntemps)
  real(kind=8), allocatable :: sqesum(:,:,:)
  ! SQE, dimension(ne,nqh,ntemps)

  complex(kind=8), allocatable :: eigenvec(:,:,:)
  ! phonon eigenvector, dimension(nqtot,nbands,nbands)

end module phon

module variables
  use control
  use phon
end module variables
