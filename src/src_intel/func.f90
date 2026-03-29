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
! Module contains some useful functions

module func
  
  use constants, only: tpi, ps2s, hbar, kbj, thz2mev
  implicit none
  
contains

  ! 3D cross product.
  subroutine cross_prod(a, b, res)
  
    real(kind=8), intent(in) :: a(3), b(3)
    real(kind=8), intent(out) :: res(3)

    integer(kind=4) :: i, j, k

    do i = 1, 3
       j = mod(i,3) + 1
       k = mod(j,3) + 1
       res(i) = a(j)*b(k)-a(k)*b(j)
    end do

  end subroutine cross_prod

  ! Quotient and remainder of integer division.
  subroutine divmod(a, b, q, r)
  
    implicit none

    integer(kind=4), intent(in) :: a, b
    integer(kind=4), intent(out) :: q, r

    q = a/b
    r = mod(a,b)

  end subroutine divmod

  ! exp(iunit*x). Sometimes it is much faster than the general
  ! complex exponential.
  function phexp(x)

    implicit none

    real(kind=8), intent(in) :: x
    complex(kind=8) :: phexp

    phexp = cmplx(cos(x), sin(x), kind=8)

  end function phexp
  
  ! Bose-Einstein distribution function
  function fBE(E, T)
  
    implicit none
    
    real(kind=8), intent(in) :: E, T
    real(kind=8) :: fBE
    
    real(kind=8) :: x
    
    x = E/ps2s
    fBE = 1.d0/(exp(hbar*x/kbj/T)-1.d0)
  
  end function fBE

  ! Bose-Einstein distribution function
  function fBEmeV(E, T)
  
    implicit none
    
    real(kind=8), intent(in) :: E, T
    real(kind=8) :: fBEmeV
    
    real(kind=8) :: x
    
    x = E/thz2mev*tpi/ps2s
    fBEmeV = 1.d0/(exp(hbar*x/kbj/T)-1.d0)
  
  end function fBEmeV
  
  ! Atomic form factor of X-ray for a given atom
  ! D. Waasmaier and A. Kirfel, Acta Cryst. A51, 416 (1995)
  function faff_wk(aff_a, aff_b, aff_c, qpt)
  
    implicit none
    
    real(kind=8), intent(in) :: aff_a(:), aff_b(:), aff_c, qpt(:)
    real(kind=8) :: faff_wk
    
    integer(kind=4) :: ii
    real(kind=8) :: q2
    
    q2 = dot_product(qpt, qpt)
    faff_wk = aff_c
    do ii = 1, 5
       faff_wk = faff_wk+aff_a(ii)*exp(-aff_b(ii)*q2)
    end do
  
  end function faff_wk
  
  ! Gaussian instrument resolution function
  function resfunc(thise, ecenter, ftype, xm, ord, a)
  
    implicit none
    
    integer(kind=4), intent(in) :: xm
    real(kind=8), intent(in) :: thise, ecenter
    character(len=*), intent(in) :: ftype
    integer(kind=4), intent(in), optional :: ord
    real(kind=8), intent(in), optional :: a(:)
    real(kind=8) :: resfunc
    
    real(kind=8) :: sigmae, pree
    
    if (ftype .eq. "CNCS12") then
       sigmae = sigma12(ecenter)
    elseif (ftype .eq. "CNCS20") then
       sigmae = sigma20(ecenter)
    elseif (ftype .eq. "Cuizhuhot") then
       sigmae = sigma_h(ecenter) 
    elseif (ftype .eq. "Cuizhucold") then
       sigmae = sigma_c(ecenter)
    elseif (ftype .eq. "poly") then
       sigmae = sigma_poly(ecenter, a, ord)
    end if
    pree = 1.d0/sqrt(tpi)/xm/sigmae
    resfunc = pree*exp(-(thise-ecenter)**2/2.d0/(xm*sigmae)**2)
    
  contains
    
    ! Energy resolution function for CNCS for Ei=12meV
    ! fitted to polynomial of order 4
    function sigma12(E)

      implicit none
      
      real(kind=8), intent(in) :: E
      real(kind=8) :: sigma12
      
      real(kind=8), parameter :: a(0:4) = (/0.4760462925081220, &
                                           -0.0461244825662230, &
                                           0.00100733850164346, &
                                           5.63726943406554E-5, &
                                           -5.91725715232495E-7/)
      real(kind=8), parameter :: b = 2.35
      
      sigma12 = (a(0)+a(1)*E+a(2)*E**2+a(3)*E**3+a(4)*E**4)/b
  
    end function sigma12
  
    ! Energy resolution function for CNCS for Ei=20meV
    ! fitted to polynomial of order 4
    function sigma20(E)
    
      implicit none
    
      real(kind=8), intent(in) :: E
      real(kind=8) :: sigma20
    
      real(kind=8), parameter :: a(0:4) = (/1.02188553396572, &
                                           -0.0594951491635698, &
                                           0.000791281838645456,&
                                           2.58258569481579E-5, &
                                           -1.64952077361938E-7/)
      real(kind=8), parameter :: b = 2.35
    
      sigma20 = (a(0)+a(1)*E+a(2)*E**2+a(3)*E**3+a(4)*E**4)/b
  
    end function sigma20

    ! Energy resolution function for Cuizhu Triple-Axis Spectrometer
    ! hot neutron, fitted to polynomial of order 4
    function sigma_h(E)
    
      implicit none
    
      real(kind=8), intent(in) :: E
      real(kind=8) :: sigma_h
    
      real(kind=8), parameter :: a(0:4) = (/1.56610382, &
                                           -0.395995079, &
                                           0.0416282809, &
                                           -1.55564158E-3, &
                                           2.15783291E-5/)
    
      sigma_h = a(0)+a(1)*E+a(2)*E**2+a(3)*E**3+a(4)*E**4
    
    end function sigma_h

    ! Energy resolution function for Cuizhu Triple-Axis Spectrometer
    ! cold neutron, fitted to polynomial of order 4
    function sigma_c(E)
    
      implicit none
    
      real(kind=8), intent(in) :: E
      real(kind=8) :: sigma_c
    
      real(kind=8), parameter :: a(0:4) = (/0.1, &
                                           -0.09271253,&
                                           0.03607527, &
                                           -0.00383385,&
                                           0.00015573/)
    
      sigma_c = a(0)+a(1)*E+a(2)*E**2+a(3)*E**3+a(4)*E**4
    
    end function sigma_c

    ! usr-defined polynomial function
    function sigma_poly(E, a, ord)
    
      implicit none
    
      integer(kind=4), intent(in) :: ord
      real(kind=8), intent(in) :: E
      real(kind=8), intent(in) :: a(:)
      real(kind=8) :: sigma_poly
    
      integer(kind=4) :: ii
    
      sigma_poly = 0.d0
      do ii = 0, ord
         sigma_poly = sigma_poly+a(ii+2)*E**ii
      end do
      sigma_poly = sigma_poly/a(1)
  
    end function sigma_poly
  
  end function resfunc
  
  function gauss(thise, ecenter, sig)
  
    implicit none
    
    real(kind=8), intent(in) :: thise, ecenter, sig
    real(kind=8) :: gauss
    
    real(kind=8) :: pree
    
    pree = 1.d0/sqrt(tpi)/sig
    gauss = pree*exp(-(thise-ecenter)**2/2.d0/sig**2)
    
  end function gauss
  
end module func
