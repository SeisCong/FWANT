module abs_mod

! This module is used for absorbing outgoing waves based on
! Cerjan and Kosloff's nonreflecting boundary condition
! (Cerjan C., et al.(1985), Geophysics, 50(4):705-708)
! with improvement of taking into account time step value
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 11:21:43 -0500 (Fri, 16 Jan 2009) $
! $Revision: 54 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************


use constants_mod
use string_mod
use macdrp_mod, only : Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz
use para_mod
use grid_mod
use mpi_mod
implicit none

private

public ::            &
    abs_init,        &
    abs_destroy,     &
    abs_LxF_LyF_LzF, &
    abs_LxB_LyB_LzB, &
    abs_LxF_LyF_LzB, &
    abs_LxB_LyB_LzF, &
    abs_LxB_LyF_LzF, &
    abs_LxF_LyB_LzB, &
    abs_LxF_LyB_LzF, &
    abs_LxB_LyF_LzB, &
    abs_syn,         &
    abs_RK_beg,      &
    abs_RK_inn,      &
    abs_RK_fin

interface abs_LxF_LyF_LzF
  module procedure fun_noop
end interface
interface abs_LxB_LyB_LzB
  module procedure fun_noop
end interface
interface abs_LxF_LyF_LzB
  module procedure fun_noop
end interface
interface abs_LxB_LyB_LzF
  module procedure fun_noop
end interface
interface abs_LxB_LyF_LzF
  module procedure fun_noop
end interface
interface abs_LxF_LyB_LzB
  module procedure fun_noop
end interface
interface abs_LxF_LyB_LzF
  module procedure fun_noop
end interface
interface abs_LxB_LyF_LzB
  module procedure fun_noop
end interface
interface abs_destroy
  module procedure fun_noop
end interface
interface abs_syn
  module procedure fun_noop
end interface

interface abs_RK_beg
  module procedure fun_noop_para
end interface
interface abs_RK_inn
  module procedure fun_noop_para
end interface
interface abs_RK_fin
  module procedure abs_apply
end interface

integer,parameter :: num_blk=6

type ELEM
     logical :: isabs
     integer :: nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz
     integer :: ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk
end type ELEM

type (ELEM), dimension(num_blk) :: W

integer,dimension(SEIS_GEO,2) :: abs_number
real(SP),dimension(:),allocatable :: Ex,Ey,Ez

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine abs_init(fnm_conf)
use mpi_mod, only : absnode
character (len=*),intent(in) :: fnm_conf
integer fid,n,m,i,j,k
real(SP),dimension(SEIS_GEO,2) :: Vs
real(SP),dimension(SEIS_GEO) :: vec1
real(SP) :: L,stepx,stepy,stepz

fid=1001
abs_number=0
Vs=0.0_SP
open(fid,file=trim(fnm_conf),status="old")
do n=1,SEIS_GEO
do m=1,2
if (absnode(n,m)) then
   call string_conf(fid,1,'abs_number',(n-1)*2+m+1,abs_number(n,m))
   call string_conf(fid,1,'abs_velocity',(n-1)*2+m+1,Vs(n,m))
end if
end do
end do
close(fid)

do n=1,num_blk
   W(n)%isabs=.false.
   W(n)%nx1=nx1; W(n)%nx2=nx1-1; W(n)%nx=0
   W(n)%ny1=ny1; W(n)%ny2=ny1-1; W(n)%ny=0
   W(n)%nz1=nz1; W(n)%nz2=nz1-1; W(n)%nz=0
end do

allocate(Ex(nx1:nx2)); Ex=1.0_SP
allocate(Ey(ny1:ny2)); Ey=1.0_SP
allocate(Ez(nz1:nz2)); Ez=1.0_SP

!x1
n=1
W(n)%ni=abs_number(1,1);W(n)%ni1=ni1       ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj             ;W(n)%nj1=nj1       ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk             ;W(n)%nk1=nk1       ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
do i=W(n)%ni1,W(n)%ni2
   j=W(n)%nj1; k=W(n)%nk1
   stepx=(x(i+1)-x(i-1))/2.0_SP
   L=stepx*xsin(i)*z(k)
   Ex(i)=cal_e(W(n)%ni -(i-W(n)%ni1),Vs(1,1),L,W(n)%ni)
end do
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
end if
!x2
n=n+1
W(n)%ni=abs_number(1,2);W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj             ;W(n)%nj1=nj1                  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk             ;W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
do i=W(n)%ni1,W(n)%ni2
   j=W(n)%nj1; k=W(n)%nk1
   stepx=(x(i+1)-x(i-1))/2.0_SP
   L=stepx*xsin(i)*z(k)
   Ex(i)=cal_e(i-W(n)%ni1+1         ,Vs(1,2),L,W(n)%ni)
end do
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
end if

!y1
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1)        ;W(n)%nj1=nj1                ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk                     ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
do j=W(n)%nj1,W(n)%nj2
   i=W(n)%ni1; k=W(n)%nk1
   stepy=(y(j+1)-y(j-1))/2.0_SP
   L=stepy*z(k)
   Ey(j)=cal_e(W(n)%nj -(j-W(n)%nj1),Vs(2,1),L,W(n)%nj)
end do
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
end if
!y2
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2)        ;W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk                     ;W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
do j=W(n)%nj1,W(n)%nj2
   i=W(n)%ni1; k=W(n)%nk1
   stepy=(y(j+1)-y(j-1))/2.0_SP
   L=stepy*z(k)
   Ey(j)=cal_e(j-W(n)%nj1+1         ,Vs(2,2),L,W(n)%nj)
end do
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
end if

!z1
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1);W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1)        ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
do k=W(n)%nk1,W(n)%nk2
   i=W(n)%ni1; j=W(n)%nj1
   stepz=(z(k+1)-z(k-1))/2.0_SP
   L=stepz
   Ez(k)=cal_e(W(n)%nk -(k-W(n)%nk1),Vs(3,1),L,W(n)%nk)
end do
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
end if
!z2
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1)  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2)        ;W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
do k=W(n)%nk1,W(n)%nk2
   i=W(n)%ni1; j=W(n)%nj1
   stepz=(z(k+1)-z(k-1))/2.0_SP
   L=stepz
   Ez(k)=cal_e(k-W(n)%nk1+1         ,Vs(3,2),L,W(n)%nk)
end do
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
end if

do n=1,num_blk
if (W(n)%isabs) then
   W(n)%nx1=W(n)%ni1;W(n)%ny1=W(n)%nj1;W(n)%nz1=W(n)%nk1
   W(n)%nx2=W(n)%ni2;W(n)%ny2=W(n)%nj2;W(n)%nz2=W(n)%nk2
end if
end do
#ifdef DEBUG
do i=ni1,ni2
write(50+myid,*) Ex(i)
end do
do j=nj1,nj2
write(50+myid,*) Ey(j)
end do
do k=nk1,nk2
write(50+myid,*) Ez(k)
end do
#endif
end subroutine abs_init

subroutine fun_noop
   return
end subroutine fun_noop
subroutine fun_noop_para(rka,rkb)
real(SP),intent(in),optional :: rka,rkb
   return
end subroutine fun_noop_para

subroutine abs_apply(rka,rkb)
real(SP),intent(in),optional :: rka,rkb
integer :: n,i,j,k
real(SP) :: d

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d=min(Ex(i),Ey(j),Ez(k))

   Txx(i,j,k)=Txx(i,j,k)*d
   Tyy(i,j,k)=Tyy(i,j,k)*d
   Tzz(i,j,k)=Tzz(i,j,k)*d
   Txy(i,j,k)=Txy(i,j,k)*d
   Txz(i,j,k)=Txz(i,j,k)*d
   Tyz(i,j,k)=Tyz(i,j,k)*d
   Vx (i,j,k)=Vx (i,j,k)*d
   Vy (i,j,k)=Vy (i,j,k)*d
   Vz (i,j,k)=Vz (i,j,k)*d
end do
end do
end do
end do
end subroutine abs_apply


!function cal_e(i,fct,at,nb) result(d)
!integer i,nb
!real(SP) :: fct,d,at
!real(SP) :: ie
!!ie=i-3.0_SP
!ie=i
!if (nb==0) then
!   d=1.0_SP
!elseif (ie<0) then
!   d=1.0_SP
!else
!   d=exp(-(fct*ie)**2)
!end if
!end function cal_e
function cal_e(i,Vs,ah,nb) result(d)
integer,intent(in) :: i,nb
real(SP),intent(in) :: Vs,ah
real(SP) :: d
real(SP) :: ie
integer m,n
ie=i
!Vs=5000.0_SP
m=(nb*ah)/(Vs*stept)
d=0.0_SP
do n=1,m
   d=d+(n*stept*Vs)**2/(nb*ah)**2
end do
d=0.8_SP/d*1.1_SP
d=exp(-d*(ie/nb)**2)
end function cal_e

end module abs_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
