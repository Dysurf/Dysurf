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
    use constants, only: tpi, thz2mev, eps5
    use variables, only: ne, deltaE, nqh, nqk, nql, ntemps, &
                         temps, sqebin, sqesum, omega, cqlist, &
                         nbands, lresfunc, degauss, order, paras, &
                         functype, xm, nqtot, eigenvec, armsd, ltds, l4d
    implicit none
    
    integer(kind=4) :: ih, ik, il, ii, jj, kk, nq
    real(kind=8) :: smear, ee(ne) 
    real(kind=8) :: sqetemp(nbands)
    
    allocate(sqebin(ne,nql,nqk,nqh,ntemps))
    if (ltds) then
       allocate(sqesum(nqk,nqh,ntemps))
    else
       allocate(sqesum(ne,nqh,ntemps))
    end if
    sqebin = 0.d0
    sqesum = 0.d0
    
    ee = (/(deltaE*(ii-1), ii = 1, ne)/)
    !ee = ee/thz2mev*tpi
  
    do ii = 1, ntemps
       nq = 0
       do ih = 1, nqh
          do ik = 1, nqk
             do il = 1, nql
                nq = nq+1
                call sqeForGivenQ(cqlist(:,nq), omega(nq,:), armsd(:,:,ii), &
                                  eigenvec(nq,:,:), sqetemp, temps(ii))
                do jj = 1, nbands
                   do kk = 1, ne
                      if (lresfunc) then
                         smear = resfunc(ee(kk), omega(nq,jj), functype, xm, order, paras)
                      else
                         smear = gauss(ee(kk), omega(nq,jj), degauss)
                      end if
                      sqebin(kk,il,ik,ih,ii) = sqebin(kk,il,ik,ih,ii)+sqetemp(jj)*smear
                   end do
                end do
                if (mod(nq, 1000) .eq. 0) then
                   write(*,*) "SQE info:", int(nq*100./nqtot), "percent of", nqtot, "q-points"
                   write(*,*)
                end if
             end do
          end do
       end do
    end do
    
    where(sqebin < eps5) sqebin = sqebin+eps5
    
    if (ltds) then
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
