module math_mod

! This module contains subroutines about math functions
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 11:34:28 -0500 (Fri, 16 Jan 2009) $
! $Revision: 58 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

implicit none

INTEGER, PARAMETER,private :: DP = KIND(1.0D0)
real(kind=DP),parameter,private ::  &
    PI=3.1415926535897932d0

#ifdef DataTypeDouble
INTEGER, PARAMETER,private :: SP = KIND(1.0D0)
#else
INTEGER, PARAMETER,private :: SP = KIND(1.0)
#endif

    interface invert
       module procedure invert_DP
       module procedure invert_real
    end interface

    interface times_product
       module procedure times_product_DP
       module procedure times_product_real
    end interface

contains

!-- determinants and inversion ---
function determinant (jac) result(det)
    real(SP),intent(in) :: jac(:,:)
    integer :: it
    real(SP) :: det
    it=ubound(jac,1)
    select case (it)
    case(1)
      det=1.0_SP
    case(2)
      det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
    case(3)
      det=jac(1,1)*( jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3) )
      det=det-jac(1,2)*( jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3) )
      det=det+jac(1,3)*( jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2) )
    case default
      print *, "wrong dimension for jacobian matrix"; stop 1
    end select
end function determinant

subroutine invert_real(matrix)
    real,intent(in out) :: matrix(:,:)
    integer :: i,k,n
    real :: con
    n=ubound(matrix,1)
    do k=1,n
       con=matrix(k,k); matrix(k,k)=1
       matrix(k,:)=matrix(k,:)/con
       do i=1,n
          if (i/=k) then
             con=matrix(i,k);matrix(i,k)=0.0_SP
             matrix(i,:)=matrix(i,:)-matrix(k,:)*con
          end if
       end do
    end do
end subroutine invert_real
subroutine invert_DP(matrix)
    real(DP),intent(in out) :: matrix(:,:)
    integer :: i,k,n
    real(DP) :: con
    n=ubound(matrix,1)
    do k=1,n
       con=matrix(k,k); matrix(k,k)=1
       matrix(k,:)=matrix(k,:)/con
       do i=1,n
          if (i/=k) then
             con=matrix(i,k);matrix(i,k)=0.0_SP
             matrix(i,:)=matrix(i,:)-matrix(k,:)*con
          end if
       end do
    end do
end subroutine invert_DP

subroutine times_product_real(A,B,C)
    real,intent(in) :: A(:),B(:)
    real,intent(out) :: C(:)
    integer n
    n=ubound(A,1)
    if (n==3) then
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
    else
      print *, "wrong dimension for times_product_real"; stop 1
    end if
end subroutine times_product_real

subroutine times_product_DP(A,B,C)
    real(DP),intent(in) :: A(:),B(:)
    real(DP),intent(out) :: C(:)
    integer n
    n=ubound(A,1)
    if (n==3) then
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
    else
      print *, "wrong dimension for times_product_DP"; stop 1
    end if
end subroutine times_product_DP

function norm(A) result(c)
   real(SP),dimension(:),intent(in) :: A
   real(SP) :: c
   c=sqrt(dot_product(A,A))
end function norm

!-- derivative and interpolate
function interp_2d(x,y,z,ni,nj,xi,yi) result(zi)
  integer,intent(in) :: ni,nj
  real(SP),intent(in) :: x(ni),y(nj),z(ni,nj)
  real(SP),intent(in) :: xi,yi
  real(SP) :: zi
  real(SP),dimension(ni) :: Lx,xt,xb
  real(SP),dimension(nj) :: Ly,yt,yb
  integer i,j

  Lx=1.0_SP
  do i=1,ni
     xb=x; xb=x(i)-xb; xb(i)=1.0_SP
     xt=x; xt=xi-xt; xt(i)=1.0_SP
     do j=1,ni
        Lx(i)=Lx(i)*xt(j)/xb(j)
     end do
  end do

  Ly=1.0_SP
  do i=1,nj
     yb=y; yb=y(i)-yb; yb(i)=1.0_SP
     yt=y; yt=yi-yt; yt(i)=1.0_SP
     do j=1,nj
        Ly(i)=Ly(i)*yt(j)/yb(j)
     end do
  end do

  zi=0.0_SP
  do j=1,nj
  do i=1,ni
     zi=zi+Lx(i)*Ly(j)*z(i,j)
  end do
  end do
end function interp_2d

function interp_3d(x,y,z,f,ni,nj,nk,xi,yi,zi) result(fi)
  integer,intent(in) :: ni,nj,nk
  real(SP),intent(in) :: x(ni),y(nj),z(nk),f(ni,nj,nk)
  real(SP),intent(in) :: xi,yi,zi
  real(SP) :: fi
  real(SP),dimension(ni) :: Lx,xt,xb
  real(SP),dimension(nj) :: Ly,yt,yb
  real(SP),dimension(nk) :: Lz,zt,zb
  integer i,j,k,n

  Lx=1.0_SP
  do i=1,ni
     xb=x; xb=x(i)-xb; xb(i)=1.0_SP
     xt=x; xt=xi-xt; xt(i)=1.0_SP
     do n=1,ni
        Lx(i)=Lx(i)*xt(n)/xb(n)
     end do
  end do

  Ly=1.0_SP
  do i=1,nj
     yb=y; yb=y(i)-yb; yb(i)=1.0_SP
     yt=y; yt=yi-yt; yt(i)=1.0_SP
     do n=1,nj
        Ly(i)=Ly(i)*yt(n)/yb(n)
     end do
  end do

  Lz=1.0_SP
  do i=1,nk
     zb=z; zb=z(i)-zb; zb(i)=1.0_SP
     zt=z; zt=zi-zt; zt(i)=1.0_SP
     do n=1,nk
        Lz(i)=Lz(i)*zt(n)/zb(n)
     end do
  end do

  fi=0.0_SP
  do k=1,nk
  do j=1,nj
  do i=1,ni
     fi=fi+Lx(i)*Ly(j)*Lz(k)*f(i,j,k)
  end do
  end do
  end do
end function interp_3d

end module math_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
