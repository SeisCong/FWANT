module abs_mod

! This module is used for absorbing outgoing waves based on
! unsplit-field ADE CFS-PML (auxiliary differential equation
! complex frequncy shifted PML)
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2008 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $
! $Revision: 510 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************
 
#include "mod_macdrp.h"
!#define AbsVzero
!#define DEBUG
!#define CorrAbs

use constants_mod
use string_mod
use para_mod
use nfseis_mod
use macdrp_mod, only :                        &
    Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,         &
    hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,hVx,hVy,hVz
use media_mod
use grid_mod
use mpi_mod
use mpi
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
  module procedure in_LxF_LyF_LzF
end interface
interface abs_LxB_LyB_LzB
  module procedure in_LxB_LyB_LzB
end interface
interface abs_LxF_LyF_LzB
  module procedure in_LxF_LyF_LzB
end interface
interface abs_LxB_LyB_LzF
  module procedure in_LxB_LyB_LzF
end interface
interface abs_LxB_LyF_LzF
  module procedure in_LxB_LyF_LzF
end interface
interface abs_LxF_LyB_LzB
  module procedure in_LxF_LyB_LzB
end interface
interface abs_LxF_LyB_LzF
  module procedure in_LxF_LyB_LzF
end interface
interface abs_LxB_LyF_LzB
  module procedure in_LxB_LyF_LzB
end interface

!-----------------------------------------------------------------------------
DEFFDWET 
DEFFDWET24
DEFFDWET22
DEFLDDRK4A
DEFLDDRK4B
HOCWETL
HOCWETR

!integer,parameter :: SP=kind(1.0)
integer,parameter :: CONSPD=2.0_SP
integer,parameter :: CONSPB=2.0_SP
integer,parameter :: CONSPA=1.0_SP

logical :: isx1,isx2,isy1,isy2,isz1,isz2

integer,dimension(SEIS_GEO,2) :: abs_number
real(SP),dimension(:),allocatable :: Ax,Ay,Az
real(SP),dimension(:),allocatable :: Bx,By,Bz
real(SP),dimension(:),allocatable :: Dx,Dy,Dz

real(SP),dimension(:,:,:),allocatable ::   &
         Vx1a, Vy1a, Vz1a, Txx1a, Txy1a, Txz1a, &  !== x1 direction ==
         Vx1b, Vy1b, Vz1b, Txx1b, Txy1b, Txz1b, &  !== x2 direction ==
        hVx1a,hVy1a,hVz1a,hTxx1a,hTxy1a,hTxz1a, &
        hVx1b,hVy1b,hVz1b,hTxx1b,hTxy1b,hTxz1b, &
        mVx1a,mVy1a,mVz1a,mTxx1a,mTxy1a,mTxz1a, &
        mVx1b,mVy1b,mVz1b,mTxx1b,mTxy1b,mTxz1b, &
        tVx1a,tVy1a,tVz1a,tTxx1a,tTxy1a,tTxz1a, &
        tVx1b,tVy1b,tVz1b,tTxx1b,tTxy1b,tTxz1b
real(SP),dimension(:,:,:),allocatable ::   &
         Vx2a, Vy2a, Vz2a, Tyy2a, Txy2a, Tyz2a, &  !== y1 direction ==
         Vx2b, Vy2b, Vz2b, Tyy2b, Txy2b, Tyz2b, &  !== y2 direction ==
        hVx2a,hVy2a,hVz2a,hTyy2a,hTxy2a,hTyz2a, &
        hVx2b,hVy2b,hVz2b,hTyy2b,hTxy2b,hTyz2b, &
        mVx2a,mVy2a,mVz2a,mTyy2a,mTxy2a,mTyz2a, &
        mVx2b,mVy2b,mVz2b,mTyy2b,mTxy2b,mTyz2b, &
        tVx2a,tVy2a,tVz2a,tTyy2a,tTxy2a,tTyz2a, &
        tVx2b,tVy2b,tVz2b,tTyy2b,tTxy2b,tTyz2b
real(SP),dimension(:,:,:),allocatable ::   &
         Vx3a, Vy3a, Vz3a, Tzz3a, Txz3a, Tyz3a, &  !== z1 direction ==
         Vx3b, Vy3b, Vz3b, Tzz3b, Txz3b, Tyz3b, &  !== z2 direction ==
        hVx3a,hVy3a,hVz3a,hTzz3a,hTxz3a,hTyz3a, &
        hVx3b,hVy3b,hVz3b,hTzz3b,hTxz3b,hTyz3b, &
        mVx3a,mVy3a,mVz3a,mTzz3a,mTxz3a,mTyz3a, &
        mVx3b,mVy3b,mVz3b,mTzz3b,mTxz3b,mTyz3b, &
        tVx3a,tVy3a,tVz3a,tTzz3a,tTxz3a,tTyz3a, &
        tVx3b,tVy3b,tVz3b,tTzz3b,tTxz3b,tTyz3b
#ifdef AnisGene
real(SP),dimension(:,:,:,:),allocatable :: &
     matCAa,matCAb,matCBa,matCBb
#endif

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine abs_init(fnm_conf)
use mpi_mod, only : absnode,freenode
character (len=*),intent(in) :: fnm_conf
integer fid,n,m,i,j,k,ierr,npt
real(SP),dimension(SEIS_GEO,2) :: Vs,fc,bmax,Rpp,dmax,amax
real(SP) :: x0,y0,z0,L0,Lx,Ly,Lz

fid=1001
abs_number=0
Vs=0.0_SP
Rpp=0.0_SP
dmax=0.0_SP
amax=0.0_SP

isx1=.false.
isx2=.false.
isy1=.false.
isy2=.false.
isz1=.false.
isz2=.false.

open(fid,file=trim(fnm_conf),status="old")
do n=1,SEIS_GEO
do m=1,2
if (absnode(n,m)) then
   call string_conf(fid,1,'abs_number',(n-1)*2+m+1,abs_number(n,m))
   ! reset absnode if layer <= 0
   if (abs_number(n,m)<=0) then
      absnode(n,m)=.false.
      cycle
   end if
   call string_conf(fid,1,'abs_velocity',(n-1)*2+m+1,Vs(n,m))
   call string_conf(fid,1,'CFS_bmax',(n-1)*2+m+1,bmax(n,m))
   !call string_conf(fid,1,'CFS_fc',(n-1)*2+m+1,fc(n,m))
   call string_conf(fid,1,'CFS_amax',(n-1)*2+m+1,amax(n,m))
   if (bmax(n,m)<1.0) then
      print *, "CFS_bmax should be large or equal to 1"
      print *, "n,m,CFS_bmax(n,m)=",n,m,bmax(n,m)
      stop 1
   end if
   Rpp(n,m)=cal_pml_R(abs_number(n,m))
end if
end do
end do
close(fid)

!=== if layer >0, deset freenode ===
if (abs_number(3,2)>0 .and. freenode) then
   freenode=.false.
end if

allocate(Ax(nx1:nx2)); Ax=0.0_SP
allocate(Bx(nx1:nx2)); Bx=1.0_SP
allocate(Dx(nx1:nx2)); Dx=0.0_SP
allocate(Ay(ny1:ny2)); Ay=0.0_SP
allocate(By(ny1:ny2)); By=1.0_SP
allocate(Dy(ny1:ny2)); Dy=0.0_SP
allocate(Az(nz1:nz2)); Az=0.0_SP
allocate(Bz(nz1:nz2)); Bz=1.0_SP
allocate(Dz(nz1:nz2)); Dz=0.0_SP

! -------------------  x layer --------------------------
! PML layer toward -x direction
if (abs_number(1,1)>0) then
   isx1=.true.
   x0=x(abs_number(1,1)+ni1)
   L0=z(nk2)*(x0-x(ni1))
   dmax(1,1)=cal_pml_dmax(L0,Vs(1,1),Rpp(1,1),bmax(1,1))
do i=ni1,abs_number(1,1)+ni1
   Lx=z(nk2)*(x0-x(i))
   Dx(i)=cal_pml_d(Lx,L0,dmax(1,1))
   Ax(i)=cal_pml_a(Lx,L0,amax(1,1))
   Bx(i)=cal_pml_b(Lx,L0,bmax(1,1))
end do
end if
! PML layer toward +x direction
if (abs_number(1,2)>0) then
   isx2=.true.
   x0=x(ni2-abs_number(1,2))
   L0=z(nk2)*(x(ni2)-x0)
   dmax(1,2)=cal_pml_dmax(L0,Vs(1,2),Rpp(1,2),bmax(1,2))
do i=ni2-abs_number(1,2),ni2
   Lx=z(nk2)*(x(i)-x0)
   Dx(i)=cal_pml_d(Lx,L0,dmax(1,2))
   Ax(i)=cal_pml_a(Lx,L0,amax(1,2))
   Bx(i)=cal_pml_b(Lx,L0,bmax(1,2))
end do
end if

! -------------------  y layer --------------------------
! PML layer toward -y direction
if (abs_number(2,1)>0) then
   isy1=.true.
   y0=y(abs_number(2,1)+nj1)
   L0=z(nk2)*(y0-y(nj1))
   dmax(2,1)=cal_pml_dmax(L0,Vs(2,1),Rpp(2,1),bmax(2,1))
do j=nj1,abs_number(2,1)+nj1
   Ly=z(nk2)*(y0-y(j))
   Dy(j)=cal_pml_d(Ly,L0,dmax(2,1))
   Ay(j)=cal_pml_a(Ly,L0,amax(2,1))
   By(j)=cal_pml_b(Ly,L0,bmax(2,1))
end do
end if
! PML layer toward +y direction
if (abs_number(2,2)>0) then
   isy2=.true.
   y0=y(nj2-abs_number(2,2))
   L0=z(nk2)*(y(nj2)-y0)
   dmax(2,2)=cal_pml_dmax(L0,Vs(2,2),Rpp(2,2),bmax(2,2))
do j=nj2-abs_number(2,2),nj2
   Ly=z(nk2)*(y(j)-y0)
   Dy(j)=cal_pml_d(Ly,L0,dmax(2,2))
   Ay(j)=cal_pml_a(Ly,L0,amax(2,2))
   By(j)=cal_pml_b(Ly,L0,bmax(2,2))
end do
end if

! -------------------  z layer --------------------------
! PML layer toward -z direction
if (abs_number(3,1)>0) then
   isz1=.true.
   z0=z(abs_number(3,1)+nk1)
   L0=z0-z(nk1)
   dmax(3,1)=cal_pml_dmax(L0,Vs(3,1),Rpp(3,1),bmax(3,1))
do k=nk1,nk1+abs_number(3,1)
   Lz=z0-z(k)
   Dz(k)=cal_pml_d(Lz,L0,dmax(3,1))
   Az(k)=cal_pml_a(Lz,L0,amax(3,1))
   Bz(k)=cal_pml_b(Lz,L0,bmax(3,1))
end do
end if

! PML layer toward +z direction
if (abs_number(3,2)>0) then
   isz2=.true.
   z0=z(nk2-abs_number(3,2))
   L0=z(nk2)-z0
   dmax(3,2)=cal_pml_dmax(L0,Vs(3,2),Rpp(3,2),bmax(3,2))
do k=nk2-abs_number(3,2),nk2
   Lz=z(k)-z0
   Dz(k)=cal_pml_d(Lz,L0,dmax(3,2))
   Az(k)=cal_pml_a(Lz,L0,amax(3,2))
   Bz(k)=cal_pml_b(Lz,L0,bmax(3,2))
end do
end if

! convert d_x to d_x/b_x
Dx=Dx/Bx
Dy=Dy/By
Dz=Dz/Bz

!=== allocate variables ===
if (isx1) then
   allocate(Vx1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(Vy1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(Vz1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(Txx1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(Txy1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(Txz1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(hVx1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(hVy1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(hVz1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(hTxx1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(hTxy1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(hTxz1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(mVx1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(mVy1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(mVz1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(mTxx1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(mTxy1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(mTxz1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(tVx1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(tVy1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(tVz1a (ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(tTxx1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(tTxy1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   allocate(tTxz1a(ni1:abs_number(1,1)+ni1,nj1:nj2,nk1:nk2))
   Vx1a =0.0; Vy1a =0.0; Vz1a =0.0
   Txx1a=0.0; Txy1a=0.0; Txz1a=0.0
   hVx1a =0.0; hVy1a =0.0; hVz1a =0.0
   hTxx1a=0.0; hTxy1a=0.0; hTxz1a=0.0
   mVx1a =0.0; mVy1a =0.0; mVz1a =0.0
   mTxx1a=0.0; mTxy1a=0.0; mTxz1a=0.0
   tVx1a =0.0; tVy1a =0.0; tVz1a =0.0
   tTxx1a=0.0; tTxy1a=0.0; tTxz1a=0.0
end if
if (isx2) then
   allocate(Vx1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(Vy1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(Vz1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(Txx1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(Txy1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(Txz1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(hVx1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(hVy1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(hVz1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(hTxx1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(hTxy1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(hTxz1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(mVx1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(mVy1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(mVz1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(mTxx1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(mTxy1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(mTxz1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(tVx1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(tVy1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(tVz1b (ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(tTxx1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(tTxy1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   allocate(tTxz1b(ni2-abs_number(1,2):ni2,nj1:nj2,nk1:nk2))
   Vx1b=0.0 ; Vy1b =0.0; Vz1b =0.0
   Txx1b=0.0; Txy1b=0.0; Txz1b=0.0
   hVx1b=0.0 ; hVy1b =0.0; hVz1b =0.0
   hTxx1b=0.0; hTxy1b=0.0; hTxz1b=0.0
   mVx1b=0.0 ; mVy1b =0.0; mVz1b =0.0
   mTxx1b=0.0; mTxy1b=0.0; mTxz1b=0.0
   tVx1b=0.0 ; tVy1b =0.0; tVz1b =0.0
   tTxx1b=0.0; tTxy1b=0.0; tTxz1b=0.0
end if

if (isy1) then
   allocate(Vx2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(Vy2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(Vz2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(Txy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(Tyy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(Tyz2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(hVx2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(hVy2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(hVz2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(hTxy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(hTyy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(hTyz2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(mVx2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(mVy2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(mVz2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(mTxy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(mTyy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(mTyz2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(tVx2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(tVy2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(tVz2a (ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(tTxy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(tTyy2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   allocate(tTyz2a(ni1:ni2,nj1:abs_number(2,1)+nj1,nk1:nk2))
   Vx2a =0.0; Vy2a =0.0; Vz2a =0.0
   Txy2a=0.0; Tyy2a=0.0; Tyz2a=0.0
   hVx2a =0.0; hVy2a =0.0; hVz2a =0.0
   hTxy2a=0.0; hTyy2a=0.0; hTyz2a=0.0
   mVx2a =0.0; mVy2a =0.0; mVz2a =0.0
   mTxy2a=0.0; mTyy2a=0.0; mTyz2a=0.0
   tVx2a =0.0; tVy2a =0.0; tVz2a =0.0
   tTxy2a=0.0; tTyy2a=0.0; tTyz2a=0.0
end if
if (isy2) then
   allocate(Vx2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(Vy2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(Vz2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(Txy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(Tyy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(Tyz2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(hVx2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(hVy2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(hVz2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(hTxy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(hTyy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(hTyz2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(mVx2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(mVy2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(mVz2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(mTxy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(mTyy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(mTyz2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(tVx2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(tVy2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(tVz2b (ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(tTxy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(tTyy2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   allocate(tTyz2b(ni1:ni2,nj2-abs_number(2,2):nj2,nk1:nk2))
   Vx2b =0.0; Vy2b =0.0; Vz2b =0.0
   Txy2b=0.0; Tyy2b=0.0; Tyz2b=0.0
   hVx2b =0.0; hVy2b =0.0; hVz2b =0.0
   hTxy2b=0.0; hTyy2b=0.0; hTyz2b=0.0
   mVx2b =0.0; mVy2b =0.0; mVz2b =0.0
   mTxy2b=0.0; mTyy2b=0.0; mTyz2b=0.0
   tVx2b =0.0; tVy2b =0.0; tVz2b =0.0
   tTxy2b=0.0; tTyy2b=0.0; tTyz2b=0.0
end if

if (isz1) then
   allocate(Vx3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(Vy3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(Vz3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(Txz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(Tyz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(Tzz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(hVx3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(hVy3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(hVz3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(hTxz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(hTyz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(hTzz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(mVx3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(mVy3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(mVz3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(mTxz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(mTyz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(mTzz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(tVx3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(tVy3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(tVz3a (ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(tTxz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(tTyz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   allocate(tTzz3a(ni1:ni2,nj1:nj2,nk1:abs_number(3,1)+nk1))
   Vx3a =0.0; Vy3a =0.0; Vz3a =0.0
   Txz3a=0.0; Tyz3a=0.0; Tzz3a=0.0
   hVx3a =0.0; hVy3a =0.0; hVz3a =0.0
   hTxz3a=0.0; hTyz3a=0.0; hTzz3a=0.0
   mVx3a =0.0; mVy3a =0.0; mVz3a =0.0
   mTxz3a=0.0; mTyz3a=0.0; mTzz3a=0.0
   tVx3a =0.0; tVy3a =0.0; tVz3a =0.0
   tTxz3a=0.0; tTyz3a=0.0; tTzz3a=0.0
end if
if (isz2) then
   allocate(Vx3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(Vy3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(Vz3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(Txz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(Tyz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(Tzz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(hVx3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(hVy3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(hVz3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(hTxz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(hTyz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(hTzz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(mVx3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(mVy3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(mVz3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(mTxz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(mTyz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(mTzz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(tVx3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(tVy3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(tVz3b (ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(tTxz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(tTyz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   allocate(tTzz3b(ni1:ni2,nj1:nj2,nk2-abs_number(3,2):nk2))
   Vx3b =0.0; Vy3b =0.0; Vz3b =0.0
   Txz3b=0.0; Tyz3b=0.0; Tzz3b=0.0
   hVx3b =0.0; hVy3b =0.0; hVz3b =0.0
   hTxz3b=0.0; hTyz3b=0.0; hTzz3b=0.0
   mVx3b =0.0; mVy3b =0.0; mVz3b =0.0
   mTxz3b=0.0; mTyz3b=0.0; mTzz3b=0.0
   tVx3b =0.0; tVy3b =0.0; tVz3b =0.0
   tTxz3b=0.0; tTyz3b=0.0; tTzz3b=0.0
end if

#ifdef AnisGene
!=== DxV[xyz] to DxV[xyz] ===
if (freenode) then
if (isx1) then
   allocate(matCAa(SEIS_GEO,SEIS_GEO,ni1:abs_number(1,1)+ni1,nj1:nj2))
   do j=nj1,nj2
   do i=ni1,ni1+abs_number(1,1)
      call coef_Dx2Dz(matCAa(:,:,i,j),i,j)
   end do
   end do
end if
if (isx2) then
   allocate(matCAb(SEIS_GEO,SEIS_GEO,ni2-abs_number(1,2):ni2,nj1:nj2))
   do j=nj1,nj2
   do i=ni2-abs_number(1,2),ni2
      call coef_Dx2Dz(matCAb(:,:,i,j),i,j)
   end do
   end do
end if
if (isy1) then
   allocate(matCBa(SEIS_GEO,SEIS_GEO,ni1:ni2,nj1:abs_number(2,1)+nj1))
   do j=nj1,nj1+abs_number(2,1)
   do i=ni1,ni2
      call coef_Dy2Dz(matCBa(:,:,i,j),i,j)
   end do
   end do
end if
if (isy2) then
   allocate(matCBb(SEIS_GEO,SEIS_GEO,ni1:ni2,nj2-abs_number(2,2):nj2))
   do j=nj2-abs_number(2,2),nj2
   do i=ni1,ni2
      call coef_Dy2Dz(matCBb(:,:,i,j),i,j)
   end do
   end do
end if
end if
#endif

#ifdef DEBUG
! -- is absorb --
write(200+myid,*) "-------- is absorb ----------"
write(200+myid,*) 'isx1: ', isx1
write(200+myid,*) 'isx2: ', isx2
write(200+myid,*) 'isy1: ', isy1
write(200+myid,*) 'isy2: ', isy2
write(200+myid,*) 'isz1: ', isz1
write(200+myid,*) 'isz2: ', isz2

! -- reflection coefficient --
write(200+myid,*) "-------- reflection coefficient ----------"
write(200+myid,"(a3,2es14.5)") 'x: ', Rpp(1,:)
write(200+myid,"(a3,2es14.5)") 'y: ', Rpp(2,:)
write(200+myid,"(a3,2es14.5)") 'z: ', Rpp(3,:)

write(200+myid,*) "-------- dmax ----------"
write(200+myid,"(a3,2es14.5)") 'x: ', dmax(1,:)
write(200+myid,"(a3,2es14.5)") 'y: ', dmax(2,:)
write(200+myid,"(a3,2es14.5)") 'z: ', dmax(3,:)

write(200+myid,*) "-------- amax ----------"
write(200+myid,"(a3,2es14.5)") 'x: ', amax(1,:)
write(200+myid,"(a3,2es14.5)") 'y: ', amax(2,:)
write(200+myid,"(a3,2es14.5)") 'z: ', amax(3,:)

! -- x damping --
if (isx1 .or. isx2) then
! Ax 
write(200+myid,*) "---- Ax(i) ----"
do i=nx1,nx2
write(200+myid,"(i6,1es14.5)") i,Ax(i)
end do
! Bx
write(200+myid,*) "---- Bx(i) ----"
do i=nx1,nx2
write(200+myid,"(i6,1es14.5)") i,Bx(i)
end do
! Dx
write(200+myid,*) "---- Dx(i) ----"
do i=nx1,nx2
write(200+myid,"(i6,1es14.5)") i,Dx(i)*Bx(i)
end do
end if

! -- y damping --
if (isy1 .or. isy2) then
! Ay
write(200+myid,*) "---- Ay(j) at half and regular points ----"
do j=ny1,ny2
write(200+myid,"(i6,1es14.5)") j,Ay(j)
end do
! By
write(200+myid,*) "---- By(j) at half and regular points ----"
do j=ny1,ny2
write(200+myid,"(i6,1es14.5)") j,By(j)
end do
! Dy
write(200+myid,*) "---- Dy(j) at half and regular points ----"
do j=ny1,ny2
write(200+myid,"(i6,1es14.5)") j,Dy(j)
end do
end if

! -- z damping --
if (isz1 .or. isz2) then
! Az
write(200+myid,*) "---- Az(k) at half and regular points ----"
do k=nz1,nz2
write(200+myid,"(i6,1es14.5)") k,Az(k)
end do
! Bz
write(200+myid,*) "---- Bz(k) at half and regular points ----"
do k=nz1,nz2
write(200+myid,"(i6,1es14.5)") k,Bz(k)
end do
! Dz
write(200+myid,*) "---- Dz(k) at half and regular points ----"
do k=nz1,nz2
write(200+myid,"(i6,1es14.5)") k,Dz(k)
end do
end if
#endif

end subroutine abs_init

subroutine abs_destroy
if (allocated(Txx1a)) deallocate(Txx1a)
if (allocated(Txy1a)) deallocate(Txy1a)
if (allocated(Txz1a)) deallocate(Txz1a)
if (allocated(Txy2a)) deallocate(Txy2a)
if (allocated(Tyy2a)) deallocate(Tyy2a)
if (allocated(Tyz2a)) deallocate(Tyz2a)
if (allocated(Txz3a)) deallocate(Txz3a)
if (allocated(Tyz3a)) deallocate(Tyz3a)
if (allocated(Tzz3a)) deallocate(Tzz3a)
if (allocated( Vx1a)) deallocate( Vx1a)
if (allocated( Vy1a)) deallocate( Vy1a)
if (allocated( Vz1a)) deallocate( Vz1a)
if (allocated( Vx2a)) deallocate( Vx2a)
if (allocated( Vy2a)) deallocate( Vy2a)
if (allocated( Vz2a)) deallocate( Vz2a)
if (allocated( Vx3a)) deallocate( Vx3a)
if (allocated( Vy3a)) deallocate( Vy3a)
if (allocated( Vz3a)) deallocate( Vz3a)
end subroutine abs_destroy

!-----------------------------------------------------------------------------

subroutine abs_syn
 integer n,i,j,k
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(n,i,j,k)

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)
    mVx1a(i,j,k)= Vx1a(i,j,k)
    mVy1a(i,j,k)= Vy1a(i,j,k)
    mVz1a(i,j,k)= Vz1a(i,j,k)
   mTxx1a(i,j,k)=Txx1a(i,j,k)
   mTxy1a(i,j,k)=Txy1a(i,j,k)
   mTxz1a(i,j,k)=Txz1a(i,j,k)
end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2
    mVx1b(i,j,k)= Vx1b(i,j,k)
    mVy1b(i,j,k)= Vy1b(i,j,k)
    mVz1b(i,j,k)= Vz1b(i,j,k)
   mTxx1b(i,j,k)=Txx1b(i,j,k)
   mTxy1b(i,j,k)=Txy1b(i,j,k)
   mTxz1b(i,j,k)=Txz1b(i,j,k)
end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2
    mVx2a(i,j,k)= Vx2a(i,j,k)
    mVy2a(i,j,k)= Vy2a(i,j,k)
    mVz2a(i,j,k)= Vz2a(i,j,k)
   mTyy2a(i,j,k)=Tyy2a(i,j,k)
   mTxy2a(i,j,k)=Txy2a(i,j,k)
   mTyz2a(i,j,k)=Tyz2a(i,j,k)
end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2
    mVx2b(i,j,k)= Vx2b(i,j,k)
    mVy2b(i,j,k)= Vy2b(i,j,k)
    mVz2b(i,j,k)= Vz2b(i,j,k)
   mTyy2b(i,j,k)=Tyy2b(i,j,k)
   mTxy2b(i,j,k)=Txy2b(i,j,k)
   mTyz2b(i,j,k)=Tyz2b(i,j,k)
end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2
    mVx3a(i,j,k)= Vx3a(i,j,k)
    mVy3a(i,j,k)= Vy3a(i,j,k)
    mVz3a(i,j,k)= Vz3a(i,j,k)
   mTzz3a(i,j,k)=Tzz3a(i,j,k)
   mTxz3a(i,j,k)=Txz3a(i,j,k)
   mTyz3a(i,j,k)=Tyz3a(i,j,k)
end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2
    mVx3b(i,j,k)= Vx3b(i,j,k)
    mVy3b(i,j,k)= Vy3b(i,j,k)
    mVz3b(i,j,k)= Vz3b(i,j,k)
   mTzz3b(i,j,k)=Tzz3b(i,j,k)
   mTxz3b(i,j,k)=Txz3b(i,j,k)
   mTyz3b(i,j,k)=Tyz3b(i,j,k)
end do
end do
end do
end if

!$OMP END PARALLEL DO
end subroutine abs_syn

subroutine abs_RK_beg(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) :: a,b
 integer :: i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)
     Vx1a(i,j,k)= mVx1a(i,j,k)+a* hVx1a(i,j,k)
     Vy1a(i,j,k)= mVy1a(i,j,k)+a* hVy1a(i,j,k)
     Vz1a(i,j,k)= mVz1a(i,j,k)+a* hVz1a(i,j,k)
    Txx1a(i,j,k)=mTxx1a(i,j,k)+a*hTxx1a(i,j,k)
    Txy1a(i,j,k)=mTxy1a(i,j,k)+a*hTxy1a(i,j,k)
    Txz1a(i,j,k)=mTxz1a(i,j,k)+a*hTxz1a(i,j,k)
    tVx1a(i,j,k)= mVx1a(i,j,k)+b* hVx1a(i,j,k)
    tVy1a(i,j,k)= mVy1a(i,j,k)+b* hVy1a(i,j,k)
    tVz1a(i,j,k)= mVz1a(i,j,k)+b* hVz1a(i,j,k)
   tTxx1a(i,j,k)=mTxx1a(i,j,k)+b*hTxx1a(i,j,k)
   tTxy1a(i,j,k)=mTxy1a(i,j,k)+b*hTxy1a(i,j,k)
   tTxz1a(i,j,k)=mTxz1a(i,j,k)+b*hTxz1a(i,j,k)
end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2
     Vx1b(i,j,k)= mVx1b(i,j,k)+a* hVx1b(i,j,k)
     Vy1b(i,j,k)= mVy1b(i,j,k)+a* hVy1b(i,j,k)
     Vz1b(i,j,k)= mVz1b(i,j,k)+a* hVz1b(i,j,k)
    Txx1b(i,j,k)=mTxx1b(i,j,k)+a*hTxx1b(i,j,k)
    Txy1b(i,j,k)=mTxy1b(i,j,k)+a*hTxy1b(i,j,k)
    Txz1b(i,j,k)=mTxz1b(i,j,k)+a*hTxz1b(i,j,k)
    tVx1b(i,j,k)= mVx1b(i,j,k)+b* hVx1b(i,j,k)
    tVy1b(i,j,k)= mVy1b(i,j,k)+b* hVy1b(i,j,k)
    tVz1b(i,j,k)= mVz1b(i,j,k)+b* hVz1b(i,j,k)
   tTxx1b(i,j,k)=mTxx1b(i,j,k)+b*hTxx1b(i,j,k)
   tTxy1b(i,j,k)=mTxy1b(i,j,k)+b*hTxy1b(i,j,k)
   tTxz1b(i,j,k)=mTxz1b(i,j,k)+b*hTxz1b(i,j,k)
end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2
     Vx2a(i,j,k)= mVx2a(i,j,k)+a* hVx2a(i,j,k)
     Vy2a(i,j,k)= mVy2a(i,j,k)+a* hVy2a(i,j,k)
     Vz2a(i,j,k)= mVz2a(i,j,k)+a* hVz2a(i,j,k)
    Tyy2a(i,j,k)=mTyy2a(i,j,k)+a*hTyy2a(i,j,k)
    Txy2a(i,j,k)=mTxy2a(i,j,k)+a*hTxy2a(i,j,k)
    Tyz2a(i,j,k)=mTyz2a(i,j,k)+a*hTyz2a(i,j,k)
    tVx2a(i,j,k)= mVx2a(i,j,k)+b* hVx2a(i,j,k)
    tVy2a(i,j,k)= mVy2a(i,j,k)+b* hVy2a(i,j,k)
    tVz2a(i,j,k)= mVz2a(i,j,k)+b* hVz2a(i,j,k)
   tTyy2a(i,j,k)=mTyy2a(i,j,k)+b*hTyy2a(i,j,k)
   tTxy2a(i,j,k)=mTxy2a(i,j,k)+b*hTxy2a(i,j,k)
   tTyz2a(i,j,k)=mTyz2a(i,j,k)+b*hTyz2a(i,j,k)
end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2
     Vx2b(i,j,k)= mVx2b(i,j,k)+a* hVx2b(i,j,k)
     Vy2b(i,j,k)= mVy2b(i,j,k)+a* hVy2b(i,j,k)
     Vz2b(i,j,k)= mVz2b(i,j,k)+a* hVz2b(i,j,k)
    Tyy2b(i,j,k)=mTyy2b(i,j,k)+a*hTyy2b(i,j,k)
    Txy2b(i,j,k)=mTxy2b(i,j,k)+a*hTxy2b(i,j,k)
    Tyz2b(i,j,k)=mTyz2b(i,j,k)+a*hTyz2b(i,j,k)
    tVx2b(i,j,k)= mVx2b(i,j,k)+b* hVx2b(i,j,k)
    tVy2b(i,j,k)= mVy2b(i,j,k)+b* hVy2b(i,j,k)
    tVz2b(i,j,k)= mVz2b(i,j,k)+b* hVz2b(i,j,k)
   tTyy2b(i,j,k)=mTyy2b(i,j,k)+b*hTyy2b(i,j,k)
   tTxy2b(i,j,k)=mTxy2b(i,j,k)+b*hTxy2b(i,j,k)
   tTyz2b(i,j,k)=mTyz2b(i,j,k)+b*hTyz2b(i,j,k)
end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2
     Vx3a(i,j,k)= mVx3a(i,j,k)+a* hVx3a(i,j,k)
     Vy3a(i,j,k)= mVy3a(i,j,k)+a* hVy3a(i,j,k)
     Vz3a(i,j,k)= mVz3a(i,j,k)+a* hVz3a(i,j,k)
    Tzz3a(i,j,k)=mTzz3a(i,j,k)+a*hTzz3a(i,j,k)
    Txz3a(i,j,k)=mTxz3a(i,j,k)+a*hTxz3a(i,j,k)
    Tyz3a(i,j,k)=mTyz3a(i,j,k)+a*hTyz3a(i,j,k)
    tVx3a(i,j,k)= mVx3a(i,j,k)+b* hVx3a(i,j,k)
    tVy3a(i,j,k)= mVy3a(i,j,k)+b* hVy3a(i,j,k)
    tVz3a(i,j,k)= mVz3a(i,j,k)+b* hVz3a(i,j,k)
   tTzz3a(i,j,k)=mTzz3a(i,j,k)+b*hTzz3a(i,j,k)
   tTxz3a(i,j,k)=mTxz3a(i,j,k)+b*hTxz3a(i,j,k)
   tTyz3a(i,j,k)=mTyz3a(i,j,k)+b*hTyz3a(i,j,k)
end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2
     Vx3b(i,j,k)= mVx3b(i,j,k)+a* hVx3b(i,j,k)
     Vy3b(i,j,k)= mVy3b(i,j,k)+a* hVy3b(i,j,k)
     Vz3b(i,j,k)= mVz3b(i,j,k)+a* hVz3b(i,j,k)
    Tzz3b(i,j,k)=mTzz3b(i,j,k)+a*hTzz3b(i,j,k)
    Txz3b(i,j,k)=mTxz3b(i,j,k)+a*hTxz3b(i,j,k)
    Tyz3b(i,j,k)=mTyz3b(i,j,k)+a*hTyz3b(i,j,k)
    tVx3b(i,j,k)= mVx3b(i,j,k)+b* hVx3b(i,j,k)
    tVy3b(i,j,k)= mVy3b(i,j,k)+b* hVy3b(i,j,k)
    tVz3b(i,j,k)= mVz3b(i,j,k)+b* hVz3b(i,j,k)
   tTzz3b(i,j,k)=mTzz3b(i,j,k)+b*hTzz3b(i,j,k)
   tTxz3b(i,j,k)=mTxz3b(i,j,k)+b*hTxz3b(i,j,k)
   tTyz3b(i,j,k)=mTyz3b(i,j,k)+b*hTyz3b(i,j,k)
end do
end do
end do
end if
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
   call swmpi_except('PML conflicts with charac boundary')
#endif
end subroutine abs_RK_beg

subroutine abs_RK_inn(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) :: a,b
 integer :: i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,k)
if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)
     Vx1a(i,j,k)= mVx1a(i,j,k)+a* hVx1a(i,j,k)
     Vy1a(i,j,k)= mVy1a(i,j,k)+a* hVy1a(i,j,k)
     Vz1a(i,j,k)= mVz1a(i,j,k)+a* hVz1a(i,j,k)
    Txx1a(i,j,k)=mTxx1a(i,j,k)+a*hTxx1a(i,j,k)
    Txy1a(i,j,k)=mTxy1a(i,j,k)+a*hTxy1a(i,j,k)
    Txz1a(i,j,k)=mTxz1a(i,j,k)+a*hTxz1a(i,j,k)
    tVx1a(i,j,k)= tVx1a(i,j,k)+b* hVx1a(i,j,k)
    tVy1a(i,j,k)= tVy1a(i,j,k)+b* hVy1a(i,j,k)
    tVz1a(i,j,k)= tVz1a(i,j,k)+b* hVz1a(i,j,k)
   tTxx1a(i,j,k)=tTxx1a(i,j,k)+b*hTxx1a(i,j,k)
   tTxy1a(i,j,k)=tTxy1a(i,j,k)+b*hTxy1a(i,j,k)
   tTxz1a(i,j,k)=tTxz1a(i,j,k)+b*hTxz1a(i,j,k)
end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2
     Vx1b(i,j,k)= mVx1b(i,j,k)+a* hVx1b(i,j,k)
     Vy1b(i,j,k)= mVy1b(i,j,k)+a* hVy1b(i,j,k)
     Vz1b(i,j,k)= mVz1b(i,j,k)+a* hVz1b(i,j,k)
    Txx1b(i,j,k)=mTxx1b(i,j,k)+a*hTxx1b(i,j,k)
    Txy1b(i,j,k)=mTxy1b(i,j,k)+a*hTxy1b(i,j,k)
    Txz1b(i,j,k)=mTxz1b(i,j,k)+a*hTxz1b(i,j,k)
    tVx1b(i,j,k)= tVx1b(i,j,k)+b* hVx1b(i,j,k)
    tVy1b(i,j,k)= tVy1b(i,j,k)+b* hVy1b(i,j,k)
    tVz1b(i,j,k)= tVz1b(i,j,k)+b* hVz1b(i,j,k)
   tTxx1b(i,j,k)=tTxx1b(i,j,k)+b*hTxx1b(i,j,k)
   tTxy1b(i,j,k)=tTxy1b(i,j,k)+b*hTxy1b(i,j,k)
   tTxz1b(i,j,k)=tTxz1b(i,j,k)+b*hTxz1b(i,j,k)
end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2
     Vx2a(i,j,k)= mVx2a(i,j,k)+a* hVx2a(i,j,k)
     Vy2a(i,j,k)= mVy2a(i,j,k)+a* hVy2a(i,j,k)
     Vz2a(i,j,k)= mVz2a(i,j,k)+a* hVz2a(i,j,k)
    Tyy2a(i,j,k)=mTyy2a(i,j,k)+a*hTyy2a(i,j,k)
    Txy2a(i,j,k)=mTxy2a(i,j,k)+a*hTxy2a(i,j,k)
    Tyz2a(i,j,k)=mTyz2a(i,j,k)+a*hTyz2a(i,j,k)
    tVx2a(i,j,k)= tVx2a(i,j,k)+b* hVx2a(i,j,k)
    tVy2a(i,j,k)= tVy2a(i,j,k)+b* hVy2a(i,j,k)
    tVz2a(i,j,k)= tVz2a(i,j,k)+b* hVz2a(i,j,k)
   tTyy2a(i,j,k)=tTyy2a(i,j,k)+b*hTyy2a(i,j,k)
   tTxy2a(i,j,k)=tTxy2a(i,j,k)+b*hTxy2a(i,j,k)
   tTyz2a(i,j,k)=tTyz2a(i,j,k)+b*hTyz2a(i,j,k)
end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2
     Vx2b(i,j,k)= mVx2b(i,j,k)+a* hVx2b(i,j,k)
     Vy2b(i,j,k)= mVy2b(i,j,k)+a* hVy2b(i,j,k)
     Vz2b(i,j,k)= mVz2b(i,j,k)+a* hVz2b(i,j,k)
    Tyy2b(i,j,k)=mTyy2b(i,j,k)+a*hTyy2b(i,j,k)
    Txy2b(i,j,k)=mTxy2b(i,j,k)+a*hTxy2b(i,j,k)
    Tyz2b(i,j,k)=mTyz2b(i,j,k)+a*hTyz2b(i,j,k)
    tVx2b(i,j,k)= tVx2b(i,j,k)+b* hVx2b(i,j,k)
    tVy2b(i,j,k)= tVy2b(i,j,k)+b* hVy2b(i,j,k)
    tVz2b(i,j,k)= tVz2b(i,j,k)+b* hVz2b(i,j,k)
   tTyy2b(i,j,k)=tTyy2b(i,j,k)+b*hTyy2b(i,j,k)
   tTxy2b(i,j,k)=tTxy2b(i,j,k)+b*hTxy2b(i,j,k)
   tTyz2b(i,j,k)=tTyz2b(i,j,k)+b*hTyz2b(i,j,k)
end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2
     Vx3a(i,j,k)= mVx3a(i,j,k)+a* hVx3a(i,j,k)
     Vy3a(i,j,k)= mVy3a(i,j,k)+a* hVy3a(i,j,k)
     Vz3a(i,j,k)= mVz3a(i,j,k)+a* hVz3a(i,j,k)
    Tzz3a(i,j,k)=mTzz3a(i,j,k)+a*hTzz3a(i,j,k)
    Txz3a(i,j,k)=mTxz3a(i,j,k)+a*hTxz3a(i,j,k)
    Tyz3a(i,j,k)=mTyz3a(i,j,k)+a*hTyz3a(i,j,k)
    tVx3a(i,j,k)= tVx3a(i,j,k)+b* hVx3a(i,j,k)
    tVy3a(i,j,k)= tVy3a(i,j,k)+b* hVy3a(i,j,k)
    tVz3a(i,j,k)= tVz3a(i,j,k)+b* hVz3a(i,j,k)
   tTzz3a(i,j,k)=tTzz3a(i,j,k)+b*hTzz3a(i,j,k)
   tTxz3a(i,j,k)=tTxz3a(i,j,k)+b*hTxz3a(i,j,k)
   tTyz3a(i,j,k)=tTyz3a(i,j,k)+b*hTyz3a(i,j,k)
end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2
     Vx3b(i,j,k)= mVx3b(i,j,k)+a* hVx3b(i,j,k)
     Vy3b(i,j,k)= mVy3b(i,j,k)+a* hVy3b(i,j,k)
     Vz3b(i,j,k)= mVz3b(i,j,k)+a* hVz3b(i,j,k)
    Tzz3b(i,j,k)=mTzz3b(i,j,k)+a*hTzz3b(i,j,k)
    Txz3b(i,j,k)=mTxz3b(i,j,k)+a*hTxz3b(i,j,k)
    Tyz3b(i,j,k)=mTyz3b(i,j,k)+a*hTyz3b(i,j,k)
    tVx3b(i,j,k)= tVx3b(i,j,k)+b* hVx3b(i,j,k)
    tVy3b(i,j,k)= tVy3b(i,j,k)+b* hVy3b(i,j,k)
    tVz3b(i,j,k)= tVz3b(i,j,k)+b* hVz3b(i,j,k)
   tTzz3b(i,j,k)=tTzz3b(i,j,k)+b*hTzz3b(i,j,k)
   tTxz3b(i,j,k)=tTxz3b(i,j,k)+b*hTxz3b(i,j,k)
   tTyz3b(i,j,k)=tTyz3b(i,j,k)+b*hTyz3b(i,j,k)
end do
end do
end do
end if
!$OMP END PARALLEL DO
end subroutine abs_RK_inn

subroutine abs_RK_fin(rkb)
 real(SP),intent(in) :: rkb
 real(SP) :: b
 integer :: i,j,k
 b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)
    Vx1a(i,j,k)= tVx1a(i,j,k)+b* hVx1a(i,j,k)
    Vy1a(i,j,k)= tVy1a(i,j,k)+b* hVy1a(i,j,k)
    Vz1a(i,j,k)= tVz1a(i,j,k)+b* hVz1a(i,j,k)
   Txx1a(i,j,k)=tTxx1a(i,j,k)+b*hTxx1a(i,j,k)
   Txy1a(i,j,k)=tTxy1a(i,j,k)+b*hTxy1a(i,j,k)
   Txz1a(i,j,k)=tTxz1a(i,j,k)+b*hTxz1a(i,j,k)
end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2
    Vx1b(i,j,k)= tVx1b(i,j,k)+b* hVx1b(i,j,k)
    Vy1b(i,j,k)= tVy1b(i,j,k)+b* hVy1b(i,j,k)
    Vz1b(i,j,k)= tVz1b(i,j,k)+b* hVz1b(i,j,k)
   Txx1b(i,j,k)=tTxx1b(i,j,k)+b*hTxx1b(i,j,k)
   Txy1b(i,j,k)=tTxy1b(i,j,k)+b*hTxy1b(i,j,k)
   Txz1b(i,j,k)=tTxz1b(i,j,k)+b*hTxz1b(i,j,k)
end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2
    Vx2a(i,j,k)= tVx2a(i,j,k)+b* hVx2a(i,j,k)
    Vy2a(i,j,k)= tVy2a(i,j,k)+b* hVy2a(i,j,k)
    Vz2a(i,j,k)= tVz2a(i,j,k)+b* hVz2a(i,j,k)
   Tyy2a(i,j,k)=tTyy2a(i,j,k)+b*hTyy2a(i,j,k)
   Txy2a(i,j,k)=tTxy2a(i,j,k)+b*hTxy2a(i,j,k)
   Tyz2a(i,j,k)=tTyz2a(i,j,k)+b*hTyz2a(i,j,k)
end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2
    Vx2b(i,j,k)= tVx2b(i,j,k)+b* hVx2b(i,j,k)
    Vy2b(i,j,k)= tVy2b(i,j,k)+b* hVy2b(i,j,k)
    Vz2b(i,j,k)= tVz2b(i,j,k)+b* hVz2b(i,j,k)
   Tyy2b(i,j,k)=tTyy2b(i,j,k)+b*hTyy2b(i,j,k)
   Txy2b(i,j,k)=tTxy2b(i,j,k)+b*hTxy2b(i,j,k)
   Tyz2b(i,j,k)=tTyz2b(i,j,k)+b*hTyz2b(i,j,k)
end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2
    Vx3a(i,j,k)= tVx3a(i,j,k)+b* hVx3a(i,j,k)
    Vy3a(i,j,k)= tVy3a(i,j,k)+b* hVy3a(i,j,k)
    Vz3a(i,j,k)= tVz3a(i,j,k)+b* hVz3a(i,j,k)
   Tzz3a(i,j,k)=tTzz3a(i,j,k)+b*hTzz3a(i,j,k)
   Txz3a(i,j,k)=tTxz3a(i,j,k)+b*hTxz3a(i,j,k)
   Tyz3a(i,j,k)=tTyz3a(i,j,k)+b*hTyz3a(i,j,k)
end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2
    Vx3b(i,j,k)= tVx3b(i,j,k)+b* hVx3b(i,j,k)
    Vy3b(i,j,k)= tVy3b(i,j,k)+b* hVy3b(i,j,k)
    Vz3b(i,j,k)= tVz3b(i,j,k)+b* hVz3b(i,j,k)
   Tzz3b(i,j,k)=tTzz3b(i,j,k)+b*hTzz3b(i,j,k)
   Txz3b(i,j,k)=tTxz3b(i,j,k)+b*hTxz3b(i,j,k)
   Tyz3b(i,j,k)=tTyz3b(i,j,k)+b*hTyz3b(i,j,k)
end do
end do
end do
end if
!$OMP END PARALLEL DO
end subroutine abs_RK_fin

!-----------------------------------------------------------------------------
subroutine in_LxF_LyF_LzF
  call LxF_LyF_LzF
end subroutine in_LxF_LyF_LzF
subroutine in_LxB_LyB_LzB
  call LxB_LyB_LzB
end subroutine in_LxB_LyB_LzB
subroutine in_LxF_LyF_LzB
  call LxF_LyF_LzB
end subroutine in_LxF_LyF_LzB
subroutine in_LxB_LyB_LzF
  call LxB_LyB_LzF
end subroutine in_LxB_LyB_LzF
subroutine in_LxB_LyF_LzF
  call LxB_LyF_LzF
end subroutine in_LxB_LyF_LzF
subroutine in_LxF_LyB_LzB
  call LxF_LyB_LzB
end subroutine in_LxF_LyB_LzB
subroutine in_LxF_LyB_LzF
  call LxF_LyB_LzF
end subroutine in_LxF_LyB_LzF
subroutine in_LxB_LyF_LzB
  call LxB_LyF_LzB
end subroutine in_LxB_LyF_LzB

!-----------------------------------------------------------------------------

   !xix = cos(y(j))**2 * (xsin(i)*xcot(i))**2 / z(k) &
   !     -sin(y(j))**2 * xsin(i)              / z(k) &
   !     -cos(y(j))    * xsin(i)**2           / z(k)
   !xiy = cos(y(j))*sin(y(j)) * (xsin(i)*xcot(i))**2 / z(k) &
   !     +cos(y(j))*sin(y(j)) * xsin(i)              / z(k) &
   !     -sin(y(j))           * xsin(i)**2           / z(k)
   !xiz = cos(y(j))           * xsin(i)**2 * xcot(i) / z(k) &
   !     -                      xsin(i)**2 * xcot(i) / z(k)
   !phx = -sin(y(j))*cos(y(j)) * xcot(i)       / z(k) &
   !      +sin(y(j))**2        * xcot(i)       / z(k) &
   !      +sin(y(j))                           / z(k)
   !rix =  cos(y(j))**2        * xsin(i)**2 * xcot(i) &
   !      +cos(y(j))*sin(y(j)) * 

subroutine LxF_LyF_LzF
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxF_LyF_LzF

subroutine LxB_LyB_LzB
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxB_LyB_LzB

subroutine LxF_LyF_LzB
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxF_LyF_LzB

subroutine LxB_LyB_LzF
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxB_LyB_LzF

subroutine LxB_LyF_LzF
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxB_LyF_LzF

subroutine LxF_LyB_LzB
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxF_LyB_LzB

subroutine LxF_LyB_LzF
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzF1(Vx,i,j,k) &
     m3d_FDzF2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzF1(Vy,i,j,k) &
     m3d_FDzF2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzF1(Vz,i,j,k) &
     m3d_FDzF2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxF_LyB_LzF

subroutine LxB_LyF_LzB
integer :: n,i,j,k
real(SP) :: DxVx,DxVy,DxVz,DxTxx,            DxTxy,      DxTxz
real(SP) :: DyVx,DyVy,DyVz,      DyTyy,      DyTxy,DyTyz
real(SP) :: DzVx,DzVy,DzVz,            DzTzz,      DzTyz,DzTxz

real(SP) :: rrho
real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
real(SP) :: tc45,tc46,tc56
real(SP) :: c2zvx,c2zvy,c2zvz
#endif

real(SP) :: b1,b1_1,ad

call abs_tsymm
call abs_vsymm

if (isx1) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1a(i,j,k))*rrho
   hTxx1a(i,j,k) = Dx(i)*DxTxx - ad * Txx1a(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1a(i,j,k))*rrho
   hTxy1a(i,j,k) = Dx(i)*DxTxy - ad * Txy1a(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1a(i,j,k))*rrho
   hTxz1a(i,j,k) = Dx(i)*DxTxz - ad * Txz1a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1a(i,j,k)
#endif
   hVx1a(i,j,k) = Dx(i)*DxVx - ad * Vx1a(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1a(i,j,k)
#endif
   hVy1a(i,j,k) = Dx(i)*DxVy - ad * Vy1a(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1a(i,j,k)
#endif
   hVz1a(i,j,k) = Dx(i)*DxVz - ad * Vz1a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1a(i,j,k)+c2zvy*Vy1a(i,j,k)+c2zvz*Vz1a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isx2) then
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bx(i)
   b1_1=b1-1.0_SP
   ad=Ax(i)+Dx(i)

   !=== correct Vx,t ===
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxx = ( DxTxx + 2.0_SP * Txz(i,j,k) ) / z(k) 
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DxTxx - b1*Txx1b(i,j,k))*rrho
   hTxx1b(i,j,k) = Dx(i)*DxTxx - ad * Txx1b(i,j,k)

   !=== correct Vy,t ===
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxy = ( DxTxy + Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DxTxy - b1*Txy1b(i,j,k))*rrho
   hTxy1b(i,j,k) = Dx(i)*DxTxy - ad * Txy1b(i,j,k)

   !=== correct Vz,t ===
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxTxz = ( DxTxz - Txx(i,j,k) + Tzz(i,j,k) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DxTxz - b1*Txz1b(i,j,k))*rrho
   hTxz1b(i,j,k) = Dx(i)*DxTxz - ad * Txz1b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,x related stress terms ===
   DxVx = (              &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVx =  ( DxVx + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc11*DxVx-b1*tc11*Vx1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc12*DxVx-b1*tc12*Vx1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc13*DxVx-b1*tc13*Vx1b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc14*DxVx-b1*tc14*Vx1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc15*DxVx-b1*tc15*Vx1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc16*DxVx-b1*tc16*Vx1b(i,j,k)
#endif
   hVx1b(i,j,k) = Dx(i)*DxVx - ad * Vx1b(i,j,k)

   !=== correct Vy,x related stress terms ===
   DxVy = (              &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVy =  DxVy / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DxVy-b1*tc66*Vy1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DxVy-b1*tc16*Vy1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DxVy-b1*tc26*Vy1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DxVy-b1*tc36*Vy1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DxVy-b1*tc46*Vy1b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DxVy-b1*tc56*Vy1b(i,j,k)
#endif
   hVy1b(i,j,k) = Dx(i)*DxVy - ad * Vy1b(i,j,k)

   !=== correct Vz,x related stress terms ===
   DxVz = (              &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)
   DxVz =  (DxVz-Vx(i,j,k)) / z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DxVz-b1*tc55*Vz1b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DxVz-b1*tc15*Vz1b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DxVz-b1*tc25*Vz1b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DxVz-b1*tc35*Vz1b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DxVz-b1*tc45*Vz1b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DxVz-b1*tc56*Vz1b(i,j,k)
#endif
   hVz1b(i,j,k) = Dx(i)*DxVz - ad * Vz1b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCAa(1,1,i,j)+tc14*matCAa(2,1,i,j)+tc13*matCAa(3,1,i,j)
      c2zvy=tc15*matCAa(1,2,i,j)+tc14*matCAa(2,2,i,j)+tc13*matCAa(3,2,i,j)
      c2zvz=tc15*matCAa(1,3,i,j)+tc14*matCAa(2,3,i,j)+tc13*matCAa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCAa(1,1,i,j)+tc24*matCAa(2,1,i,j)+tc23*matCAa(3,1,i,j)
      c2zvy=tc25*matCAa(1,2,i,j)+tc24*matCAa(2,2,i,j)+tc23*matCAa(3,2,i,j)
      c2zvz=tc25*matCAa(1,3,i,j)+tc24*matCAa(2,3,i,j)+tc23*matCAa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
      != Txy =
      c2zvx=tc56*matCAa(1,1,i,j)+tc46*matCAa(2,1,i,j)+tc36*matCAa(3,1,i,j)
      c2zvy=tc56*matCAa(1,2,i,j)+tc46*matCAa(2,2,i,j)+tc36*matCAa(3,2,i,j)
      c2zvz=tc56*matCAa(1,3,i,j)+tc46*matCAa(2,3,i,j)+tc36*matCAa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DxVx+c2zvy*DxVy+c2zvz*DxVz)                    &
        +b1*(c2zvx*Vx1b(i,j,k)+c2zvy*Vy1b(i,j,k)+c2zvz*Vz1b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc13/tc33*DxVx+b1*tc13*tc13/tc33*Vx1b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc13/tc33*DxVx+b1*tc23*tc13/tc33*Vx1b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy1) then
do k=nk1,nk2
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2a(i,j,k))*rrho
   hTxy2a(i,j,k) = Dy(j)*DyTxy - ad * Txy2a(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2a(i,j,k))*rrho
   hTyy2a(i,j,k) = Dy(j)*DyTyy - ad * Tyy2a(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2a(i,j,k))*rrho
   hTyz2a(i,j,k) = Dy(j)*DyTyz - ad * Tyz2a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2a(i,j,k)
#endif
   hVx2a(i,j,k) = Dy(j)*DyVx - ad * Vx2a(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2a(i,j,k)
#endif
   hVy2a(i,j,k) = Dy(j)*DyVy - ad * Vy2a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2a(i,j,k)
#endif
   hVz2a(i,j,k) = Dy(j)*DyVz - ad * Vz2a(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2a(i,j,k)+c2zvy*Vy2a(i,j,k)+c2zvz*Vz2a(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2a(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2a(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isy2) then
do k=nk1,nk2
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/By(j)
   b1_1=b1-1.0_SP
   ad=Ay(j)+Dy(j)

   !=== correct Vx,t ===
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTxy = ( DyTxy/xsin(i) +(Txx(i,j,k)- Tyy(i,j,k))*xcot(i)+Txz(i,j,k) ) / z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DyTxy - b1*Txy2b(i,j,k))*rrho
   hTxy2b(i,j,k) = Dy(j)*DyTxy - ad * Txy2b(i,j,k)

   !=== correct Vy,t ===
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTyy = ( DyTyy/xsin(i) + 2.0*Txy(i,j,k)*xcot(i) + 2.0*Tyz(i,j,k) ) / z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DyTyy - b1*Tyy2b(i,j,k))*rrho
   hTyy2b(i,j,k) = Dy(j)*DyTyy - ad * Tyy2b(i,j,k)

   !=== correct Vz,t ===
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyTyz = ( DyTyz/xsin(i) + Txz(i,j,k)*xcot(i) + (Tzz(i,j,k)-Tyy(i,j,k)) ) / z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DyTyz - b1*Tyz2b(i,j,k))*rrho
   hTyz2b(i,j,k) = Dy(j)*DyTyz - ad * Tyz2b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,y related stress terms ===
   DyVx = (              &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVx = ( DyVx/xsin(i) - Vy(i,j,k)*xcot(i) ) / z(k)

   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc66*DyVx-b1*tc66*Vx2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc16*DyVx-b1*tc16*Vx2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc26*DyVx-b1*tc26*Vx2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc36*DyVx-b1*tc36*Vx2b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc46*DyVx-b1*tc46*Vx2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc56*DyVx-b1*tc56*Vx2b(i,j,k)
#endif
   hVx2b(i,j,k) = Dy(j)*DyVx - ad * Vx2b(i,j,k)

   !=== correct Vy,y related stress terms ===
   DyVy = (              &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVy = ( DyVy/xsin(i) + Vx(i,j,k)*xcot(i) + Vz(i,j,k) ) / z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc12*DyVy-b1*tc12*Vy2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc22*DyVy-b1*tc22*Vy2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc23*DyVy-b1*tc23*Vy2b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc24*DyVy-b1*tc24*Vy2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc25*DyVy-b1*tc25*Vy2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc26*DyVy-b1*tc26*Vy2b(i,j,k)
#endif
   hVy2b(i,j,k) = Dy(j)*DyVy - ad * Vy2b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DyVz = (              &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)
   DyVz = (DyVz/xsin(i) - Vy(i,j,k)) / z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DyVz-b1*tc44*Vz2b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DyVz-b1*tc14*Vz2b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DyVz-b1*tc24*Vz2b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DyVz-b1*tc34*Vz2b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DyVz-b1*tc45*Vz2b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DyVz-b1*tc46*Vz2b(i,j,k)
#endif
   hVz2b(i,j,k) = Dy(j)*DyVz - ad * Vz2b(i,j,k)

   !=== correct free surface boundary condition if freenode ===
   if (freenode .and. k==nk2) then
#ifdef AnisGene
      != Txx =
      c2zvx=tc15*matCBa(1,1,i,j)+tc14*matCBa(2,1,i,j)+tc13*matCBa(3,1,i,j)
      c2zvy=tc15*matCBa(1,2,i,j)+tc14*matCBa(2,2,i,j)+tc13*matCBa(3,2,i,j)
      c2zvz=tc15*matCBa(1,3,i,j)+tc14*matCBa(2,3,i,j)+tc13*matCBa(3,3,i,j)
      hTxx(i,j,k)=hTxx(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Tyy =
      c2zvx=tc25*matCBa(1,1,i,j)+tc24*matCBa(2,1,i,j)+tc23*matCBa(3,1,i,j)
      c2zvy=tc25*matCBa(1,2,i,j)+tc24*matCBa(2,2,i,j)+tc23*matCBa(3,2,i,j)
      c2zvz=tc25*matCBa(1,3,i,j)+tc24*matCBa(2,3,i,j)+tc23*matCBa(3,3,i,j)
      hTyy(i,j,k)=hTyy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
      != Txy =
      c2zvx=tc56*matCBa(1,1,i,j)+tc46*matCBa(2,1,i,j)+tc36*matCBa(3,1,i,j)
      c2zvy=tc56*matCBa(1,2,i,j)+tc46*matCBa(2,2,i,j)+tc36*matCBa(3,2,i,j)
      c2zvz=tc56*matCBa(1,3,i,j)+tc46*matCBa(2,3,i,j)+tc36*matCBa(3,3,i,j)
      hTxy(i,j,k)=hTxy(i,j,k)                                       &
        -b1_1*(c2zvx*DyVx+c2zvy*DyVy+c2zvz*DyVz)                    &
        +b1*(c2zvx*Vx2b(i,j,k)+c2zvy*Vy2b(i,j,k)+c2zvz*Vz2b(i,j,k))
#else                                                                       
      != Txx =
      hTxx(i,j,k)=hTxx(i,j,k)-b1_1*tc13*tc23/tc33*DyVy+b1*tc13*tc23/tc33*Vy2b(i,j,k)
      != Tyy =
      hTyy(i,j,k)=hTyy(i,j,k)-b1_1*tc23*tc23/tc33*DyVy+b1*tc23*tc23/tc33*Vy2b(i,j,k)
#endif
   end if

end do
end do
end do
end if

if (isz1) then
do k=nk1,nk1+abs_number(3,1)
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3a(i,j,k))*rrho
   hTxz3a(i,j,k) = Dz(k)*DzTxz - ad * Txz3a(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3a(i,j,k))*rrho
   hTyz3a(i,j,k) = Dz(k)*DzTyz - ad * Tyz3a(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3a(i,j,k))*rrho
   hTzz3a(i,j,k) = Dz(k)*DzTzz - ad * Tzz3a(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3a(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3a(i,j,k)
#endif
   hVx3a(i,j,k) = Dz(k)*DzVx - ad * Vx3a(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3a(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3a(i,j,k)
#endif
   hVy3a(i,j,k) = Dz(k)*DzVy - ad * Vy3a(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3a(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3a(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3a(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3a(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3a(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3a(i,j,k)
#endif
   hVz3a(i,j,k) = Dz(k)*DzVz - ad * Vz3a(i,j,k)

end do
end do
end do
end if

if (isz2) then
do k=nk2-abs_number(3,2),nk2
do j=nj1,nj2
do i=ni1,ni2

   rrho=1.0/rho(i,j,k)
   b1=1.0_SP/Bz(k)
   b1_1=b1-1.0_SP
   ad=Az(k)+Dz(k)

   !=== correct Vx,t ===
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   hVx(i,j,k) =  hVx(i,j,k) + (b1_1*DzTxz - b1*Txz3b(i,j,k))*rrho
   hTxz3b(i,j,k) = Dz(k)*DzTxz - ad * Txz3b(i,j,k)

   !=== correct Vy,t ===
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
   hVy(i,j,k) =  hVy(i,j,k) + (b1_1*DzTyz - b1*Tyz3b(i,j,k))*rrho
   hTyz3b(i,j,k) = Dz(k)*DzTyz - ad * Tyz3b(i,j,k)

   !=== correct Vz,t ===
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   hVz(i,j,k) =  hVz(i,j,k) + (b1_1*DzTzz - b1*Tzz3b(i,j,k))*rrho
   hTzz3b(i,j,k) = Dz(k)*DzTzz - ad * Tzz3b(i,j,k)

!===== hook equation =====
#if defined AnisGene
   tc11=C11(i,j,k)
   tc12=C12(i,j,k)
   tc13=C13(i,j,k)
   tc14=C14(i,j,k)
   tc15=C15(i,j,k)
   tc16=C16(i,j,k)
   tc22=C22(i,j,k)
   tc23=C23(i,j,k)
   tc24=C24(i,j,k)
   tc25=C25(i,j,k)
   tc26=C26(i,j,k)
   tc33=C33(i,j,k)
   tc34=C34(i,j,k)
   tc35=C35(i,j,k)
   tc36=C36(i,j,k)
   tc44=C44(i,j,k)
   tc45=C45(i,j,k)
   tc46=C46(i,j,k)
   tc55=C55(i,j,k)
   tc56=C56(i,j,k)
   tc66=C66(i,j,k)
#elif defined AnisVTI
   tc11=C11(i,j,k)
   tc13=C13(i,j,k)
   tc33=C33(i,j,k)
   tc44=C44(i,j,k)
   tc66=C66(i,j,k)

   tc12=tc11-2.0_SP*tc66
   tc22=tc11
   tc23=tc13
   tc55=tc44
#else
   tc13=lambda(i,j,k)
   tc66=mu(i,j,k)

   tc11=tc13+2.0_SP*tc66
   tc12=tc13
   tc22=tc11
   tc23=tc13
   tc33=tc11
   tc44=tc66
   tc55=tc66
#endif

   !=== correct Vx,z related stress terms ===
   DzVx = (              &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )*zeta_z(k)

   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc55*DzVx-b1*tc55*Vx3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc15*DzVx-b1*tc15*Vx3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc25*DzVx-b1*tc25*Vx3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc35*DzVx-b1*tc35*Vx3b(i,j,k)
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc45*DzVx-b1*tc45*Vx3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc56*DzVx-b1*tc56*Vx3b(i,j,k)
#endif
   hVx3b(i,j,k) = Dz(k)*DzVx - ad * Vx3b(i,j,k)

   !=== correct Vy,z related stress terms ===
   DzVy = (              &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )*zeta_z(k)

   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc44*DzVy-b1*tc44*Vy3b(i,j,k)
#ifdef AnisGene
   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc14*DzVy-b1*tc14*Vy3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc24*DzVy-b1*tc24*Vy3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc34*DzVy-b1*tc34*Vy3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc45*DzVy-b1*tc45*Vy3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc46*DzVy-b1*tc46*Vy3b(i,j,k)
#endif
   hVy3b(i,j,k) = Dz(k)*DzVy - ad * Vy3b(i,j,k)

   !=== correct Vz,y related stress terms ===
   DzVz = (              &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )*zeta_z(k)

   hTxx(i,j,k)=hTxx(i,j,k)+b1_1*tc13*DzVz-b1*tc13*Vz3b(i,j,k)
   hTyy(i,j,k)=hTyy(i,j,k)+b1_1*tc23*DzVz-b1*tc23*Vz3b(i,j,k)
   hTzz(i,j,k)=hTzz(i,j,k)+b1_1*tc33*DzVz-b1*tc33*Vz3b(i,j,k)
#ifdef AnisGene
   hTyz(i,j,k)=hTyz(i,j,k)+b1_1*tc34*DzVz-b1*tc34*Vz3b(i,j,k)
   hTxz(i,j,k)=hTxz(i,j,k)+b1_1*tc35*DzVz-b1*tc35*Vz3b(i,j,k)
   hTxy(i,j,k)=hTxy(i,j,k)+b1_1*tc36*DzVz-b1*tc36*Vz3b(i,j,k)
#endif
   hVz3b(i,j,k) = Dz(k)*DzVz - ad * Vz3b(i,j,k)

end do
end do
end do
end if
end subroutine LxB_LyF_LzB

!=======================================================
subroutine abs_tsymm
#ifdef AnisGene
integer i,j,m
if (freenode) then
if (isx1) then
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)
   do m=1,LenFD
      Txx(i,j,nk2+m)=2.0*Txx(i,j,nk2)-Txx(i,j,nk2-m)
      Tyy(i,j,nk2+m)=2.0*Tyy(i,j,nk2)-Tyy(i,j,nk2-m)
      Txy(i,j,nk2+m)=2.0*Txy(i,j,nk2)-Txy(i,j,nk2-m)
   end do
end do
end do
end if

if (isx2) then
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2
   do m=1,LenFD
      Txx(i,j,nk2+m)=2.0*Txx(i,j,nk2)-Txx(i,j,nk2-m)
      Tyy(i,j,nk2+m)=2.0*Tyy(i,j,nk2)-Tyy(i,j,nk2-m)
      Txy(i,j,nk2+m)=2.0*Txy(i,j,nk2)-Txy(i,j,nk2-m)
   end do
end do
end do
end if

if (isy1) then
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2
   do m=1,LenFD
      Txx(i,j,nk2+m)=2.0*Txx(i,j,nk2)-Txx(i,j,nk2-m)
      Tyy(i,j,nk2+m)=2.0*Tyy(i,j,nk2)-Tyy(i,j,nk2-m)
      Txy(i,j,nk2+m)=2.0*Txy(i,j,nk2)-Txy(i,j,nk2-m)
   end do
end do
end do
end if

if (isy2) then
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2
   do m=1,LenFD
      Txx(i,j,nk2+m)=2.0*Txx(i,j,nk2)-Txx(i,j,nk2-m)
      Tyy(i,j,nk2+m)=2.0*Tyy(i,j,nk2)-Tyy(i,j,nk2-m)
      Txy(i,j,nk2+m)=2.0*Txy(i,j,nk2)-Txy(i,j,nk2-m)
   end do
end do
end do
end if

end if ! freenode
#endif
end subroutine abs_tsymm

subroutine abs_vsymm
integer i,j,n,m

if (freenode) then
if (isx1) then
do j=nj1,nj2
do i=ni1,ni1+abs_number(1,1)
   do m=1,LenFD
      Vx(i,j,nk2+m)=2.0*Vx(i,j,nk2)-Vx(i,j,nk2-m)
      Vy(i,j,nk2+m)=2.0*Vy(i,j,nk2)-Vy(i,j,nk2-m)
      Vz(i,j,nk2+m)=2.0*Vz(i,j,nk2)-Vz(i,j,nk2-m)
   end do
end do
end do
end if

if (isx2) then
do j=nj1,nj2
do i=ni2-abs_number(1,2),ni2
   do m=1,LenFD
      Vx(i,j,nk2+m)=2.0*Vx(i,j,nk2)-Vx(i,j,nk2-m)
      Vy(i,j,nk2+m)=2.0*Vy(i,j,nk2)-Vy(i,j,nk2-m)
      Vz(i,j,nk2+m)=2.0*Vz(i,j,nk2)-Vz(i,j,nk2-m)
   end do
end do
end do
end if

if (isy1) then
do j=nj1,nj1+abs_number(2,1)
do i=ni1,ni2
   do m=1,LenFD
      Vx(i,j,nk2+m)=2.0*Vx(i,j,nk2)-Vx(i,j,nk2-m)
      Vy(i,j,nk2+m)=2.0*Vy(i,j,nk2)-Vy(i,j,nk2-m)
      Vz(i,j,nk2+m)=2.0*Vz(i,j,nk2)-Vz(i,j,nk2-m)
   end do
end do
end do
end if

if (isy2) then
do j=nj2-abs_number(2,2),nj2
do i=ni1,ni2
   do m=1,LenFD
      Vx(i,j,nk2+m)=2.0*Vx(i,j,nk2)-Vx(i,j,nk2-m)
      Vy(i,j,nk2+m)=2.0*Vy(i,j,nk2)-Vy(i,j,nk2-m)
      Vz(i,j,nk2+m)=2.0*Vz(i,j,nk2)-Vz(i,j,nk2-m)
   end do
end do
end do
end if

end if ! freenode
end subroutine abs_vsymm

!--------------------------------------------------------------------}
function cal_pml_dmax(L,Vp,Rpp,bmax) result(dmax)
real(SP) :: L,Vp,Rpp,bmax
real(SP) :: dmax
!dmax=-(p+1.0_SP)*Vp/2.0_SP/L*log(Rpp)
dmax=-Vp/2.0_SP/L*log(Rpp)*(CONSPD+1.0_SP)
!dmax=dmax*(CONSPD+CONSPB+1)/( (CONSPD+1.0_SP)*bmax+CONSPB )
end function cal_pml_dmax

function cal_pml_d(x,L,dmax) result(d)
real(SP) :: x,L,dmax
real(SP) :: d
!real(SP),parameter :: p=2.0_SP

!dmax=-(p+1.0_SP)*Vp/2.0_SP/L*log(Rpp)

if (x<0) then
   d=0.0
else
   d=dmax*(x/L)**CONSPD
end if
end function cal_pml_d

function cal_pml_amax(fc) result(amax)
real(SP) :: fc
real(SP) :: amax
amax=PI*fc
end function cal_pml_amax

function cal_pml_a(x,L,amax) result(a)
real(SP) :: x,L,amax
real(SP) :: a
!amax=PI*fc
if (x<0) then
   a=0.0
else
   !a=amax*(L-x)/L
   a=amax*(1.0_SP-(x/L)**CONSPA)
end if
end function cal_pml_a

function cal_pml_b(x,L,bmax) result(b)
real(SP) :: x,L,bmax
real(SP) :: b
!real(SP),parameter :: p=2.0_SP
if (x<0) then
   b=1.0
else
   b=1.0+(bmax-1.0)*(x/L)**CONSPB
end if
end function cal_pml_b
function cal_pml_R(n) result(r)
integer :: n
real(SP) :: r
r = real(10**(-( log10(real(N,DP)) - 1.0_DP)/log10(2.0_DP) -3.0_DP),SP)
end function cal_pml_R

#ifdef AnisGene
subroutine coef_Dx2Dz(CA,i,j)
  use math_mod, only : invert
  real(SP),dimension(SEIS_GEO,SEIS_GEO),intent(out) :: CA
  integer,intent(in) :: i,j

  real(SP),dimension(SEIS_GEO,SEIS_GEO) :: A,C
  integer :: k

  k=nk2
  !----
  C(1,:)=(/ C35(i,j,k),C34(i,j,k),C23(i,j,k) /)
  C(2,:)=(/ C45(i,j,k),C44(i,j,k),C34(i,j,k) /)
  C(3,:)=(/ C55(i,j,k),C45(i,j,k),C35(i,j,k) /)
  call invert(C)
  !----
  A(1,:)=(/ C13(i,j,k),C36(i,j,k),C35(i,j,k) /)
  A(2,:)=(/ C14(i,j,k),C46(i,j,k),C45(i,j,k) /)
  A(3,:)=(/ C15(i,j,k),C56(i,j,k),C55(i,j,k) /)
  !----
  CA=matmul(C,A)
end subroutine coef_Dx2Dz

subroutine coef_Dy2Dz(CA,i,j)
  use math_mod, only : invert
  real(SP),dimension(SEIS_GEO,SEIS_GEO),intent(out) :: CB
  integer,intent(in) :: i,j

  real(SP),dimension(SEIS_GEO,SEIS_GEO) :: B,C
  integer :: k

  k=nk2
  !----
  C(1,:)=(/ C35(i,j,k),C34(i,j,k),C23(i,j,k) /)
  C(2,:)=(/ C45(i,j,k),C44(i,j,k),C34(i,j,k) /)
  C(3,:)=(/ C55(i,j,k),C45(i,j,k),C35(i,j,k) /)
  call invert(C)
  !----
  B(1,:)=(/ C36(i,j,k),C23(i,j,k),C34(i,j,k) /)
  B(2,:)=(/ C46(i,j,k),C24(i,j,k),C44(i,j,k) /)
  B(3,:)=(/ C56(i,j,k),C25(i,j,k),C45(i,j,k) /)
  !----
  CB=matmul(C,B)
end subroutine coef_Dx2Dz
#endif

end module abs_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
