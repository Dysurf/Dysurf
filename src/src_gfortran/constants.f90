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
! Module that contains physical constants, atomic masses and coherent neutron scattering length

module constants

  implicit none

  ! Mathematical constants
  real(kind=8), parameter :: pi = 3.141592653589793238d0
  real(kind=8), parameter :: tpi = 2.d0*pi
  real(kind=8), parameter :: fpi = 4.d0*pi
  real(kind=8), parameter :: one = 1.d0
  real(kind=8), parameter :: two = 2.d0
  real(kind=8), parameter :: zero = 0.d0
  complex(kind=8), parameter :: ci = (0.d0, 1.d0)
  complex(kind=8), parameter :: cone = (1.d0, 0.d0)
  complex(kind=8), parameter :: czero = (0.d0, 0.d0)

  ! Physical constants
  real(kind=8), parameter :: kbj = 1.380648813E-23 ! J/K
  real(kind=8), parameter :: kb = 8.61732835E-2 ! meV/K
  real(kind=8), parameter :: hj = 6.62607015E-34 ! J*s
  real(kind=8), parameter :: hh = 4.13566733E-12 ! meV*s
  !real(kind=8), parameter :: hbar = 6.5821189855E-13 ! meV*s
  real(kind=8), parameter :: hbar = 1.054571800E-34 ! J*s
  real(kind=8), parameter :: amu = 1.6605402E-27  ! kg
  real(kind=8), parameter :: ev = 1.60217733E-19 ! J
  real(kind=8), parameter :: vdc = 8.854187817E-12 ! F/m, vacuum dielectric constant
  real(kind=8), parameter :: afj = 6.02E23 !avogadro's number
  real(kind=8), parameter :: barn = 1E-24 !cm^2

  ! Unit conversion factors
  real(kind=8), parameter :: ang2m = 1.0d-10
  real(kind=8), parameter :: ryd2ev = 13.6056981d0
  real(kind=8), parameter :: bohr2ang = 0.52917721092d0
  real(kind=8), parameter :: dm2thz = 9648.53336213d0 ! from eV/(A^2*amu) to THz^2
  real(kind=8), parameter :: thz2amua3 = 1745913.10995753d0 ! THz^2*amu*A^3
  real(kind=8), parameter :: thz2mev = 4.13566733d0
  real(kind=8), parameter :: ps2s = 1.0E-12
  real(kind=8), parameter :: a2m = 1.0E-10
  real(kind=8), parameter :: vasp2thz = 15.633302d0
  real(kind=8), parameter :: qe2thz = 108.97077d0
  real(kind=8), parameter :: vasp2nac = 14.399652d0
  real(kind=8), parameter :: qe2nac = 2.d0
  
  ! ... zero up to a given accuracy
  real(kind=8), parameter :: eps1  = 1.0E-1
  real(kind=8), parameter :: eps2  = 1.0E-2
  real(kind=8), parameter :: eps3  = 1.0E-3
  real(kind=8), parameter :: eps4  = 1.0E-4
  real(kind=8), parameter :: eps5  = 1.0E-5
  real(kind=8), parameter :: eps6  = 1.0E-6
  real(kind=8), parameter :: eps7  = 1.0E-7
  real(kind=8), parameter :: eps8  = 1.0E-8
  
  ! List of elements, ordered by atomic number.
  character(len=3), parameter :: periodic_table(96) = [character(len=3) ::&
       "H","He",&
       "Li","Be","B","C","N","O","F","Ne",&
       "Na","Mg","Al","Si","P","S","Cl","Ar",&
       "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",&
       "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",&
       "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",&
       "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",&
       "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm"]
  
  ! Information about scattering length
  ! https://www.nist.gov/ncnr/planning-your-experiment/sld-periodic-table
  ! Note: the scattering length of some elements have the imaginary part,
  ! use the only real part isotopes instead.
  ! Take care of [108Cd,154Sm,156Gd,156Dy,Po,At,Rn,Fr,Ra,Ac,238Pu,244Cm]
  real(kind=8), parameter :: scatt_b(96) = [real(kind=8) ::&
  -3.7390,3.263,&   !H->He
  -1.90,7.79,6.65,6.6460,9.36,5.803,5.654,4.566,&   !Li->Ne
  3.63,5.375,3.449,4.1491,5.13,2.847,9.5770,1.909,&   !Na->Ar
  3.67,4.70,12.29,-3.438,-0.3824,3.635,-3.73,9.45,2.49,10.3,7.718,5.680,7.288,8.185,6.58,7.970,6.795,7.81,&   !K-Kr
  7.09,7.02,7.75,7.16,7.054,6.715,6.8,7.03,5.88,5.91,5.922,5.4,5.39,6.225,5.57,5.80,5.28,4.92,&   !Ru->Xe
  5.42,5.07,8.24,4.48,4.58,7.69,12.6,9.3,8.22,6.3,7.38,6.1,8.01,7.79,7.07,12.43,7.21,&   !Cs->Lu
  7.7,6.91,4.86,9.2,10.7,10.6,9.60,7.63,12.692,8.776,9.405,8.532,0.0,0.0,0.0,&   !Hf->Rn
  0.0,10.0,0.0,10.31,9.1,8.417,10.55,14.1,8.3,9.5]   !Fr->Cm
  
  real(kind=8), parameter :: atomic_masses(96) = [real(kind=8) ::&
  1.00794, 4.002602,&   !H->He
  6.941,9.012182,10.811,12.0107,14.0067,15.9994,18.9984032,20.1797,&   !Li->Ne
  22.98976928,24.3050,26.9815386,28.0855,30.973762,32.065,35.453,39.948,&   !Na->Ar
  39.0983,40.078,44.955912,47.867,50.9415,51.9961,54.938045,55.845,58.933195,&   !K->Co
  58.6934,63.546,65.38,69.723,72.64,74.92160,78.96,79.904,83.798,&   !Ni->Kr
  85.4678,87.62,88.90585,91.224,92.90638,95.96,98,101.07,102.90550,&   !Ru->Rh
  106.42,107.8682,112.411,114.818,118.710,121.760,127.60,126.90447,131.293,&   !Pd->Xe
  132.9054519,137.327,138.90547,140.116,140.90765,144.242,145,150.36,151.964,&   !Cs->Eu
  157.25,158.92535,162.500,164.93032,167.259,168.93421,173.054,174.9668,&   !Gd->Lu
  178.49,180.94788,183.84,186.207,190.23,192.217,195.084,196.966569,200.59,204.3833,207.2,208.98040,209,210,222,&   !Hf->Rn
  223,226,227,232.03806,231.0358,238.02891,237,244,243,247]   !Fr->Cm

  ! These data from https://www.ncnr.nist.gov/resources/n-lengths/list.html
  ! Note that； Some data donnot include, I replaced them zero(Data on some radioactive elements are inaccurate)
  real(kind=8), parameter :: scatt_totoal_cross_section(96) = [real(kind=8) ::&
  82.02, 1.34,&   !H->He
  1.37,7.63,5.24,5.551,11.51,4.232,4.018,2.628,&   !Li->Ne
  3.28,3.71,1.503,2.167,3.312,1.026,16.8,0.683,&   !Na->Ar
  1.96,2.83,23.5,4.35,5.1,3.49,2.15,11.62,5.6,18.5,8.03,4.131,6.83,8.6,5.5,8.3,5.9,7.68,&   !K-Kr
  6.8,6.25,7.7,6.46,6.255,5.71,6.3,6.6,4.6,4.48,4.99,6.5,2.62,4.892,3.9,4.32,3.81,2.96,&   !Ru->Xe
  3.9,3.38,9.66,2.94,2.66,16.6,21.3,39,9.2,180,6.84,90.3,8.42,8.7,6.38,23.4,7.2,&   !Cs->Lu
  10.2,6.01,4.6,11.5,14.7,14,11.71,7.75,26.8,9.89,11.118,9.156,0,0,12.6,&   !Hf->Rn
  0,13,0,13.36,10.5,8.908,14.5,0,9.0,0]  !Fr->Cm

  real(kind=8), parameter :: scatt_absob_cross_section(96) = [real(kind=8) ::&
  0.3326,0.00747,&   !H->He
  70.5,0.0076,767,0.0035,1.9,0.00019,0.0096,0.039,&   !Li->Ne
  0.53,0.063,0.231,0.171,0.172,0.53,33.5,0.675,&   !Na->Ar
  2.1,0.43,27.5,6.09,5.08,3.05,13.3,2.56,37.18,4.49,3.78,1.11,2.75,2.2,4.5,11.7,6.9,25,&   !K-Kr
  0.38,1.28,1.28,0.185,1.15,2.48,20,2.56,144.8,6.9,63.3,2520,193.8,0.626,4.91,4.7,6.15,23.9,& !Ru->Xe
  29.0,1.1,8.97,0.63,11.5,50.5,168.4,5922,4530,49700,23.4,994,64.7,159,100,34.8,74,&   !Cs->Lu
  104.1,20.6,18.3,89.7,16,425,10.3,98.65,372.30,3.43,0.171,0.0338,0,0,0,&   !Hf->Rn
  0,12.8,0,7.37,200.6,7.57,175.9,0,75.3,0]  !Fr->Cm

end module constants
