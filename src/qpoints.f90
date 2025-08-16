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

! Bin box definition and q-points generation

module qpoints
  
  implicit none

contains
  
  ! generate a list of q-points in the defined bin box
  ! all in conventional lattice basis unless indicated
  subroutine load_qpoints()

    use variables, only: crlatvec, rlatvec, nqh, nqk, nql, path, nqtot, &
                         qlist, pqlist, cqlist, deltaH, deltaK, deltaL, &
                         q0, elist, qdist, ltds, l4d
    use constants, only: eps5
    implicit none
    
    integer(kind=4) :: ii, jj, kk, ll
    integer(kind=4) :: ih, ik, il
    integer(kind=4) :: info, int_arr(3)
    real(kind=8) :: qh, qk, ql, norm
    real(kind=8) :: qq(3), qt(3), dist(3), work(9)
    real(kind=8) :: rotmat(3,3), convert(3,3)
    
    real(kind=8), external :: dnrm2

    nqtot = nqh*nqk*nql
    allocate(qlist(3,nqtot), pqlist(3,nqtot), cqlist(3,nqtot), elist(nqh), qdist(nqh))

    ! calculate the conversion matrix from conventional basis to primitive basis
    convert = rlatvec
    call dgetrf(3, 3, convert, 3, int_arr, info)
    call dgetri(3, convert, 3, int_arr, work, 9, info)
    convert = matmul(convert, crlatvec)

    ! determine the two path-vector perpendicular to the long path
    ! Schmidt orthogonalization is used
    rotmat(:,1) = path(:,1)
    if (all(path(:,2) .lt. eps5)) then
       path(:,2) = (/2.d0, 1.d0, 2.d0/)  ! a trial vector 2 1 2
       path(:,2) = path(:,2)-&
                   path(:,1)*dot_product(path(:,1), path(:,2))/dot_product(path(:,1), path(:,1))
       norm = dnrm2(3, path(:,2), 1)
       path(:,2) = path(:,2)/norm
    end if
    rotmat(:,2) = path(:,2)
    if (all(path(:,3) .lt. eps5)) then
       path(:,3) = (/4.d0, 3.d0, 1.d0/)  ! a trial vector 4 3 1
       path(:,3) = path(:,3)-&
                   path(:,1)*dot_product(path(:,1), path(:,3))/dot_product(path(:,1), path(:,1))-&
                   path(:,2)*dot_product(path(:,2), path(:,3))/dot_product(path(:,2), path(:,2))
       norm = dnrm2(3, path(:,3), 1)
       path(:,3) = path(:,3)/norm
    end if
    rotmat(:,3) = path(:,3)

    ii = 0
    ll = 0
    do ih = 1, nqh
       qh = deltaH*(ih-1)
       do ik = 1, nqk
          if (ltds .or. l4d) then
             jj = ik-1
          else
             jj = ik-1-nqk/2
          end if
          qk = deltaK*jj
          do il = 1, nql
             if (l4d) then
                kk = il-1
             else
                kk = il-1-nql/2
             end if
             ql = deltaL*kk

             ii = ii+1
             qt = (/qh, qk, ql/)
             qq = q0+matmul(rotmat(:,:), qt)
             qlist(:,ii) = qq
             cqlist(:,ii) = matmul(crlatvec, qq)
             pqlist(:,ii) = matmul(convert, qq)
             if ((jj .eq. 0) .and. (kk .eq. 0)) then
                ll = ll+1
                elist(ll) = ii
                dist = cqlist(:,ii)-cqlist(:,elist(1))
                qdist(ll) = dnrm2(3, dist, 1)
             end if
          end do
       end do
    end do

  end subroutine load_qpoints
    
end module qpoints   
