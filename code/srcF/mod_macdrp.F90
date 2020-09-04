module macdrp_mod

! This module contains the variables and subroutines
! used in the DRP/opt MacCormack fd operator
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 11:57:10 -0500 (Fri, 16 Jan 2009) $
! $Revision: 59 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"
#define WATER

use constants_mod, only : SEIS_GEO
use math_mod
use para_mod
use mpi
use mpi_mod
use media_mod
use grid_mod

implicit none
private
public ::                                      &
    Txx, Tyy, Txy, Vx, Vy, Tzz, Txz, Tyz, Vz,  &
    hTxx,hTyy,hTxy,hVx,hVy,hTzz,hTxz,hTyz,hVz, &
    macdrp_init,                               &
    macdrp_syn,                                &
    macdrp_mesg_init,                          &
    macdrp_destroy,                            &
    macdrp_LxF_LyF_LzF,                        &
    macdrp_LxB_LyB_LzB,                        &
    macdrp_LxB_LyB_LzF,                        &
    macdrp_LxF_LyF_LzB,                        &
    macdrp_LxB_LyF_LzF,                        &
    macdrp_LxF_LyB_LzB,                        &
    macdrp_LxF_LyB_LzF,                        &
    macdrp_LxB_LyF_LzB,                        &
    macdrp_RK_beg,                             &
    macdrp_RK_inn,                             &
    macdrp_RK_fin,                             &
    macdrp_check,                              &
    atten_graves

interface macdrp_LxF_LyF_LzF
  module procedure in_LxF_LyF_LzF
end interface
interface macdrp_LxB_LyB_LzB
  module procedure in_LxB_LyB_LzB
end interface
interface macdrp_LxF_LyF_LzB
  module procedure in_LxF_LyF_LzB
end interface
interface macdrp_LxB_LyB_LzF
  module procedure in_LxB_LyB_LzF
end interface
interface macdrp_LxB_LyF_LzF
  module procedure in_LxB_LyF_LzF
end interface
interface macdrp_LxF_LyB_LzB
  module procedure in_LxF_LyB_LzB
end interface
interface macdrp_LxF_LyB_LzF
  module procedure in_LxF_LyB_LzF
end interface
interface macdrp_LxB_LyF_LzB
  module procedure in_LxB_LyF_LzB
end interface

DEFFDWET
DEFFDWET24
DEFFDWET22
DEFLDDRK2A
DEFLDDRK2B
DEFLDDRK4A
DEFLDDRK4B
HOCWETL
HOCWETR

#ifdef MPIBuffered
real(SP),dimension(:),allocatable ::        &
     BufX1,BufX2,BufY1,BufY2,BufZ1,BufZ2,   &
     RevX1,RevX2,RevY1,RevY2,RevZ1,RevZ2
integer,parameter,private :: NREQ=4
integer,private :: &
     NBufXL, NBufXS, &
     NBufYL, NBufYS, &
     NBufZL, NBufZS
#else
!integer,parameter,private :: NREQ=36
integer,parameter,private :: NREQ=24
#endif

real(SP),dimension(:,:,:),allocatable ::        &
      Txx, Tyy, Txy, Vx, Vy, Tzz, Txz, Tyz, Vz, &
     hTxx,hTyy,hTxy,hVx,hVy,hTzz,hTxz,hTyz,hVz, &
     mTxx,mTyy,mTxy,mVx,mVy,mTzz,mTxz,mTyz,mVz, &
     tTxx,tTyy,tTxy,tVx,tVy,tTzz,tTxz,tTyz,tVz
real(SP),dimension(:,:,:,:),allocatable,public :: &
     matVx2Vz,matVy2Vz
real(SP),dimension(:,:),allocatable,public :: &
     TxSrc,TySrc,TzSrc,                       &
     VxSrc,VySrc,VzSrc
real(SP),dimension(4),public :: firRKa,firRKb, secRKa,secRKb
integer,dimension(SEIS_GEO*2,SEIS_GEO*2+1),public :: indx
integer ierr
integer,dimension(MPI_STATUS_SIZE) :: istatus
integer,dimension(NREQ) :: reqXB, reqXF, reqYB, reqYF, reqZB, reqZF
integer,dimension(MPI_STATUS_SIZE,NREQ) :: reqstat
#ifdef VERBOSE
integer fid_out
#endif

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine macdrp_init
integer ierr
allocate( Txx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txx=0.0_SP
allocate( Tyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tyy=0.0_SP
allocate( Txy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txy=0.0_SP
allocate( Vx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vx =0.0_SP
allocate( Vy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vy =0.0_SP
allocate( Tzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tzz=0.0_SP
allocate( Txz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txz=0.0_SP
allocate( Tyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tyz=0.0_SP
allocate( Vz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vz =0.0_SP
allocate(hTxx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxx=0.0_SP
allocate(hTyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTyy=0.0_SP
allocate(hTxy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxy=0.0_SP
allocate(hVx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVx =0.0_SP
allocate(hVy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVy =0.0_SP
allocate(hTzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTzz=0.0_SP
allocate(hTxz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxz=0.0_SP
allocate(hTyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTyz=0.0_SP
allocate(hVz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVz =0.0_SP
allocate(mTxx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mTxx=0.0_SP
allocate(mTyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mTyy=0.0_SP
allocate(mTxy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mTxy=0.0_SP
allocate(mVx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mVx =0.0_SP
allocate(mVy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mVy =0.0_SP
allocate(mTzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mTzz=0.0_SP
allocate(mTxz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mTxz=0.0_SP
allocate(mTyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mTyz=0.0_SP
allocate(mVz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); mVz =0.0_SP
allocate(tTxx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tTxx=0.0_SP
allocate(tTyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tTyy=0.0_SP
allocate(tTxy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tTxy=0.0_SP
allocate(tVx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tVx =0.0_SP
allocate(tVy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tVy =0.0_SP
allocate(tTzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tTzz=0.0_SP
allocate(tTxz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tTxz=0.0_SP
allocate(tTyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tTyz=0.0_SP
allocate(tVz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); tVz =0.0_SP
if (ierr>0) then
   print *, "can't allocate variable in macdrp_init"
   stop 1
end if
allocate(matVx2Vz(SEIS_GEO,SEIS_GEO,nx,ny),stat=ierr); matVx2Vz=0.0_SP
allocate(matVy2Vz(SEIS_GEO,SEIS_GEO,nx,ny),stat=ierr); matVy2Vz=0.0_SP
allocate(TxSrc(nx1:nx2,ny1:ny2),stat=ierr); TxSrc=0.0_SP
allocate(TySrc(nx1:nx2,ny1:ny2),stat=ierr); TySrc=0.0_SP
allocate(TzSrc(nx1:nx2,ny1:ny2),stat=ierr); TzSrc=0.0_SP
allocate(VxSrc(nx1:nx2,ny1:ny2),stat=ierr); VxSrc=0.0_SP
allocate(VySrc(nx1:nx2,ny1:ny2),stat=ierr); VySrc=0.0_SP
allocate(VzSrc(nx1:nx2,ny1:ny2),stat=ierr); VzSrc=0.0_SP
#ifdef MPIBuffered
allocate(BufX1(LenFD*nj*nk*6),stat=ierr); BufX1=0.0_SP
allocate(BufX2(LenFD*nj*nk*6),stat=ierr); BufX2=0.0_SP
allocate(BufY1(LenFD*ni*nk*6),stat=ierr); BufY1=0.0_SP
allocate(BufY2(LenFD*ni*nk*6),stat=ierr); BufY2=0.0_SP
allocate(BufZ1(LenFD*ni*nj*6),stat=ierr); BufZ1=0.0_SP
allocate(BufZ2(LenFD*ni*nj*6),stat=ierr); BufZ2=0.0_SP
allocate(RevX1(LenFD*nj*nk*6),stat=ierr); RevX1=0.0_SP
allocate(RevX2(LenFD*nj*nk*6),stat=ierr); RevX2=0.0_SP
allocate(RevY1(LenFD*ni*nk*6),stat=ierr); RevY1=0.0_SP
allocate(RevY2(LenFD*ni*nk*6),stat=ierr); RevY2=0.0_SP
allocate(RevZ1(LenFD*ni*nj*6),stat=ierr); RevZ1=0.0_SP
allocate(RevZ2(LenFD*ni*nj*6),stat=ierr); RevZ2=0.0_SP
NBufXL=nj*nk*LenFDL*6
NBufXS=nj*nk*LenFDS*6
NBufYL=ni*nk*LenFDL*6
NBufYS=ni*nk*LenFDS*6
NBufZL=ni*nj*LenFDL*6
NBufZS=ni*nj*LenFDS*6
#endif
! main
indx(:,SEIS_GEO*2+1)=(/ ni1+LenFD,ni2-LenFD, &
                        nj1+LenFD,nj2-LenFD, &
                        nk1+LenFD,nk2-LenFD /)
indx(:,SEIS_GEO*2  )=(/ ni1,ni2,nj1,nj2,nk2-LenFD+1,nk2 /) ! z2
indx(:,SEIS_GEO*2-1)=(/ ni1,ni2,nj1,nj2,nk1,nk1+LenFD-1 /) ! z1
indx(:,1)=(/ ni1,ni2,nj1,nj1+LenFD-1,nk1+LenFD,nk2-LenFD /) ! y1
indx(:,2)=(/ ni1,ni2,nj2-LenFD+1,nj2,nk1+LenFD,nk2-LenFD /) ! y2
indx(:,3)=(/ ni1,ni1+LenFD-1,nj1+LenFD,nj2-LenFD,nk1+LenFD,nk2-LenFD /) ! x1
indx(:,4)=(/ ni2-LenFD+1,ni2,nj1+LenFD,nj2-LenFD,nk1+LenFD,nk2-LenFD /) ! x2
! rk coefficient
firRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
firRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
!secRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
!secRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
secRKa=(/ RK2a2, 0.0, 0.0, 0.0 /)
secRKb=(/ RK2b1, RK2b2, 0.0, 0.0 /)
! mat to convert V,z
#ifdef VERBOSE
  fid_out=9050
  open(fid_out,                                                                      &
       file='log_maxval_'//trim(set_mpi_subfix(thisid(1),thisid(2),thisid(3)))//'.dat', &
       status='unknown')
#endif
end subroutine macdrp_init
subroutine macdrp_destroy
deallocate( Txx, Tyy, Txy, Vx, Vy)
deallocate(hTxx,hTyy,hTxy,hVx,hVy)
deallocate(mTxx,mTyy,mTxy,mVx,mVy)
deallocate(tTxx,tTyy,tTxy,tVx,tVy)
#ifdef VERBOSE
  close(fid_out)
#endif
end subroutine macdrp_destroy
subroutine macdrp_check(ntime)
integer,intent(in) :: ntime
real(SP) :: V1,V2,V3,T11,T22,T33,T12,T13,T23,W
integer ierr
#ifndef CheckOverFlow
 return
#endif
if (mod(ntime,1)==0) then
    V1=maxval(abs(Vx))
    V2=maxval(abs(Vy))
    V3=maxval(abs(Vz))
   T11=maxval(abs(Txx))
   T22=maxval(abs(Tyy))
   T33=maxval(abs(Tzz))
   T12=maxval(abs(Txy))
   T13=maxval(abs(Txz))
   T23=maxval(abs(Tyz))
#ifdef VERBOSE
   write(fid_out,'(i5,9es12.5)') ntime, V1,V2,V3,T11,T22,T33,T12,T13,T23
#endif
   W=max(V1,V2,V3,T11,T22,T33,T12,T13,T23)
   if (W>=huge(1.0)) then
      print *, "Overflow error: "
      write(*,"(i5,i3.2,2(i2.2),9(es12.5,3i5))") ntime,thisid(1),thisid(2),thisid(3), &
         V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
        T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
        T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))
#ifdef VERBOSE
      write(fid_out,"(i5,9(es12.5,3i5))") ntime,                                      &
         V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
        T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
        T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))
#endif
      call MPI_ABORT(SWMPI_COMM,1,ierr)
   end if
end if
end subroutine macdrp_check

!{----------- 4-6 LDDRK stages ----------------

subroutine macdrp_syn

 integer i,j,k
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
  mTxx(i,j,k)=Txx(i,j,k)
  mTyy(i,j,k)=Tyy(i,j,k)
  mTxy(i,j,k)=Txy(i,j,k)
  mVx (i,j,k)=Vx (i,j,k)
  mVy (i,j,k)=Vy (i,j,k)
  mTzz(i,j,k)=Tzz (i,j,k)
  mTxz(i,j,k)=Txz (i,j,k)
  mTyz(i,j,k)=Tyz (i,j,k)
  mVz (i,j,k)=Vz  (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
end subroutine macdrp_syn

!-- Generic RK --
subroutine macdrp_RK_beg(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) :: a,b
 integer i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
    Txx(i,j,k)=mTxx(i,j,k)+a*hTxx(i,j,k)
    Tyy(i,j,k)=mTyy(i,j,k)+a*hTyy(i,j,k)
    Tzz(i,j,k)=mTzz(i,j,k)+a*hTzz(i,j,k)
    Txy(i,j,k)=mTxy(i,j,k)+a*hTxy(i,j,k)
    Txz(i,j,k)=mTxz(i,j,k)+a*hTxz(i,j,k)
    Tyz(i,j,k)=mTyz(i,j,k)+a*hTyz(i,j,k)
    Vx (i,j,k)=mVx (i,j,k)+a*hVx (i,j,k)
    Vy (i,j,k)=mVy (i,j,k)+a*hVy (i,j,k)
    Vz (i,j,k)=mVz (i,j,k)+a*hVz (i,j,k)

    tTxx(i,j,k)=mTxx(i,j,k)+b*hTxx(i,j,k)
    tTyy(i,j,k)=mTyy(i,j,k)+b*hTyy(i,j,k)
    tTzz(i,j,k)=mTzz(i,j,k)+b*hTzz(i,j,k)
    tTxy(i,j,k)=mTxy(i,j,k)+b*hTxy(i,j,k)
    tTxz(i,j,k)=mTxz(i,j,k)+b*hTxz(i,j,k)
    tTyz(i,j,k)=mTyz(i,j,k)+b*hTyz(i,j,k)
    tVx (i,j,k)=mVx (i,j,k)+b*hVx (i,j,k)
    tVy (i,j,k)=mVy (i,j,k)+b*hVy (i,j,k)
    tVz (i,j,k)=mVz (i,j,k)+b*hVz (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
call free_charac
call free_extrap
#endif
end subroutine macdrp_RK_beg

subroutine macdrp_RK_inn(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) :: a,b
 integer i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
    Txx(i,j,k)=mTxx(i,j,k)+a*hTxx(i,j,k)
    Tyy(i,j,k)=mTyy(i,j,k)+a*hTyy(i,j,k)
    Tzz(i,j,k)=mTzz(i,j,k)+a*hTzz(i,j,k)
    Txy(i,j,k)=mTxy(i,j,k)+a*hTxy(i,j,k)
    Txz(i,j,k)=mTxz(i,j,k)+a*hTxz(i,j,k)
    Tyz(i,j,k)=mTyz(i,j,k)+a*hTyz(i,j,k)
    Vx (i,j,k)=mVx (i,j,k)+a*hVx (i,j,k)
    Vy (i,j,k)=mVy (i,j,k)+a*hVy (i,j,k)
    Vz (i,j,k)=mVz (i,j,k)+a*hVz (i,j,k)

    tTxx(i,j,k)=tTxx(i,j,k)+b*hTxx(i,j,k)
    tTyy(i,j,k)=tTyy(i,j,k)+b*hTyy(i,j,k)
    tTzz(i,j,k)=tTzz(i,j,k)+b*hTzz(i,j,k)
    tTxy(i,j,k)=tTxy(i,j,k)+b*hTxy(i,j,k)
    tTxz(i,j,k)=tTxz(i,j,k)+b*hTxz(i,j,k)
    tTyz(i,j,k)=tTyz(i,j,k)+b*hTyz(i,j,k)
    tVx (i,j,k)=tVx (i,j,k)+b*hVx (i,j,k)
    tVy (i,j,k)=tVy (i,j,k)+b*hVy (i,j,k)
    tVz (i,j,k)=tVz (i,j,k)+b*hVz (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
call free_charac
call free_extrap
#endif
end subroutine macdrp_RK_inn

subroutine macdrp_RK_fin(rkb)
 real(SP),intent(in) :: rkb
 real(SP) :: b
 integer i,j,k
 b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
    Txx(i,j,k)=tTxx(i,j,k)+b*hTxx(i,j,k)
    Tyy(i,j,k)=tTyy(i,j,k)+b*hTyy(i,j,k)
    Tzz(i,j,k)=tTzz(i,j,k)+b*hTzz(i,j,k)
    Txy(i,j,k)=tTxy(i,j,k)+b*hTxy(i,j,k)
    Txz(i,j,k)=tTxz(i,j,k)+b*hTxz(i,j,k)
    Tyz(i,j,k)=tTyz(i,j,k)+b*hTyz(i,j,k)
    Vx (i,j,k)=tVx (i,j,k)+b*hVx (i,j,k)
    Vy (i,j,k)=tVy (i,j,k)+b*hVy (i,j,k)
    Vz (i,j,k)=tVz (i,j,k)+b*hVz (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
call free_charac
call free_extrap
#endif
end subroutine macdrp_RK_fin

subroutine atten_graves
 integer :: i,j,k
 real(SP) :: Qatt
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
 if (Qs(i,j,k)<QsINF) then
    ! Qs is Q value not attenuation coefficient now
    Qatt=exp((-PI*QsF0*stept)/Qs(i,j,k))
    Txx(i,j,k)=Txx(i,j,k)*Qatt
    Tyy(i,j,k)=Tyy(i,j,k)*Qatt
    Tzz(i,j,k)=Tzz(i,j,k)*Qatt
    Txy(i,j,k)=Txy(i,j,k)*Qatt
    Txz(i,j,k)=Txz(i,j,k)*Qatt
    Tyz(i,j,k)=Tyz(i,j,k)*Qatt
    Vx (i,j,k)=Vx (i,j,k)*Qatt
    Vy (i,j,k)=Vy (i,j,k)*Qatt
    Vz (i,j,k)=Vz (i,j,k)*Qatt
 end if
 end do
 end do
 end do
!$OMP END PARALLEL DO
end subroutine atten_graves
!---------------------------------------------------}

!{------------------------ DRP/opt macdrp --------------------------------

!{----- wrapper of LxF_LyF_LzF ---------
subroutine in_LxF_LyF_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1
#ifdef MPIBuffered
  call fill_buff_LxF
  call fill_buff_LyF
  call fill_buff_LzF
#endif
  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxF_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxF
  call recv_buff_LyF
  call recv_buff_LzF
#endif

  do n=1,SEIS_GEO*2
  call LxF_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyF_LzF_VHOC
#endif

end subroutine in_LxF_LyF_LzF
subroutine in_LxB_LyB_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif

#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  !Txx=real(myid)
  n=SEIS_GEO*2+1
  !if (myid==0) then
  !  Vx=1.0;Vy=2.0;Vz=3.0;Txx=4.0;Tyy=5.0;Tzz=6.0;Txy=7.0;Txz=8.0;Tyz=9.0
  !end if

#ifdef MPIBuffered
  call fill_buff_LxB
  call fill_buff_LyB
  call fill_buff_LzB
#endif
  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxB_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxB
  call recv_buff_LyB
  call recv_buff_LzB
#endif

  do n=1,SEIS_GEO*2
  call LxB_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyB_LzB_VHOC
#endif
end subroutine in_LxB_LyB_LzB
subroutine in_LxB_LyB_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif

#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1

#ifdef MPIBuffered
  call fill_buff_LxB
  call fill_buff_LyB
  call fill_buff_LzF
#endif
  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxB_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxB
  call recv_buff_LyB
  call recv_buff_LzF
#endif

  do n=1,SEIS_GEO*2
  call LxB_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyB_LzF_VHOC
#endif
end subroutine in_LxB_LyB_LzF
subroutine in_LxF_LyF_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif

#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1

#ifdef MPIBuffered
  call fill_buff_LxF
  call fill_buff_LyF
  call fill_buff_LzB
#endif
  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxF_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxF
  call recv_buff_LyF
  call recv_buff_LzB
#endif

  do n=1,SEIS_GEO*2
  call LxF_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyF_LzB_VHOC
#endif
end subroutine in_LxF_LyF_LzB
subroutine in_LxB_LyF_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1

#ifdef MPIBuffered
  call fill_buff_LxB
  call fill_buff_LyF
  call fill_buff_LzF
#endif
  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxB_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxB
  call recv_buff_LyF
  call recv_buff_LzF
#endif

  do n=1,SEIS_GEO*2
  call LxB_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyF_LzF_VHOC
#endif
end subroutine in_LxB_LyF_LzF
subroutine in_LxF_LyB_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1

#ifdef MPIBuffered
  call fill_buff_LxF
  call fill_buff_LyB
  call fill_buff_LzB
#endif
  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxF_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxF
  call recv_buff_LyB
  call recv_buff_LzB
#endif

  do n=1,SEIS_GEO*2
  call LxF_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyB_LzB_VHOC
#endif
end subroutine in_LxF_LyB_LzB
subroutine in_LxF_LyB_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1

#ifdef MPIBuffered
  call fill_buff_LxF
  call fill_buff_LyB
  call fill_buff_LzF
#endif
  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxF_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxF
  call recv_buff_LyB
  call recv_buff_LzF
#endif

  do n=1,SEIS_GEO*2
  call LxF_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyB_LzF_VHOC
#endif
end subroutine in_LxF_LyB_LzF
subroutine in_LxB_LyF_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
#ifdef CondFreeTIMG
  if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
  if (freenode) call free_vext
#endif
  n=SEIS_GEO*2+1

#ifdef MPIBuffered
  call fill_buff_LxB
  call fill_buff_LyF
  call fill_buff_LzB
#endif
  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxB_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)
#ifdef MPIBuffered
  call recv_buff_LxB
  call recv_buff_LyF
  call recv_buff_LzB
#endif

  do n=1,SEIS_GEO*2
  call LxB_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyF_LzB_VHOC
#endif
end subroutine in_LxB_LyF_LzB
!-- private subroutine --
subroutine LxF_LyF_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyF_LzF

subroutine LxB_LyB_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyB_LzB

subroutine LxF_LyF_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif


   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyF_LzB

subroutine LxB_LyB_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyB_LzF

subroutine LxB_LyF_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyF_LzF

subroutine LxF_LyB_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyB_LzB

subroutine LxF_LyB_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyB_LzF

subroutine LxB_LyF_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )*xi_x(i)
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )*xi_x(i)
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )*xi_x(i)
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )*xi_x(i)
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )*xi_x(i)
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )*xi_x(i)

   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )*eta_y(j)
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )*eta_y(j)
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )*eta_y(j)
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )*eta_y(j)
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )*eta_y(j)
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )*eta_y(j)

   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )*zeta_z(k)
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )*zeta_z(k)
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )*zeta_z(k)
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )*zeta_z(k)
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )*zeta_z(k)
#ifdef CondFreeVLOW
   end if
#endif

#ifdef WATER
   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if
#endif

   hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
        +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
   hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
        +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
   hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
        +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )

   E11=(DxVx+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
   E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP

   E33=DzVz
   E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))/2.0_SP
#ifndef CondFreeCharac
#ifndef CondFreeVHOC
   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if
#endif
#endif

   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyF_LzB

#ifdef CondFreeVHOC
!*************************************************************************
!* use Compact MacCormack scheme to calculate velocities fd with         *
!* respect to eta and assemble the right hand side to update stresses    *
!*************************************************************************
subroutine LxF_LyF_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyF_LzF_VHOC
subroutine LxB_LyB_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyB_LzB_VHOC

subroutine LxF_LyF_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyF_LzB_VHOC
subroutine LxB_LyB_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyB_LzF_VHOC

subroutine LxB_LyF_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyF_LzF_VHOC
subroutine LxF_LyB_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyB_LzB_VHOC

subroutine LxF_LyB_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyB_LzF_VHOC
subroutine LxB_LyF_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: E11,E22,E33,E12,E13,E23

if (.not. freenode) return

!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)

   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)

   DzVx(n) = (-DxVz(n)/z(k)+Vx(i,j,k)/z(k)+VxSrc(i,j)/miu)/zeta_z(k)
   DzVy(n) = (-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k)+VySrc(i,j)/miu)/zeta_z(k)
   DzVz(n) = (-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu)/zeta_z(k)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )*xi_x(i)
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )*xi_x(i)
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )*xi_x(i)

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )*eta_y(j)
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )*eta_y(j)
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )*eta_y(j)

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   DzVx(n)=DzVx(n)*zeta_z(k)
   DzVy(n)=DzVy(n)*zeta_z(k)
   DzVz(n)=DzVz(n)*zeta_z(k)
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   E11=(DxVx(n)+Vz(i,j,k))/z(k)
   E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
   E33=DzVz(n)
   E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)/2.0_SP
   E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))/2.0_SP
   E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))/2.0_SP
   hTxx(i,j,k)=lam2mu*E11+lam*(E22+E33)
   hTyy(i,j,k)=lam2mu*E22+lam*(E11+E33)
   hTzz(i,j,k)=lam2mu*E33+lam*(E11+E22)
   hTxy(i,j,k)=2.0_SP*miu*E12
   hTxz(i,j,k)=2.0_SP*miu*E13
   hTyz(i,j,k)=2.0_SP*miu*E23
   if (k==nk2) then
      hTxz(i,j,k)=VxSrc(i,j)
      hTyz(i,j,k)=VySrc(i,j)
      hTzz(i,j,k)=VzSrc(i,j)
   end if
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyF_LzB_VHOC
#endif

subroutine free_charac
integer i,j,k
real(SP) v1,v2,v3,t11,t22,t33,t12,t13,t23
real(SP) lam,miu,lam2mu,rrho,f1,f2,fct
real(SP) :: Tx,Ty,Tz

if (.not. freenode) return

   k=nk2
do j=nj1,nj2
do i=ni1,ni2

   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)

   Tx=-TxSrc(i,j); Ty=-TySrc(i,j); Tz= TzSrc(i,j)
   Vx (i,j,k)=v1-(t13+Tx)/f2
   Vy (i,j,k)=v2-(t23+Ty)/f2
   Vz (i,j,k)=v3-(t33-Tz)/f1
   Txx(i,j,k)=t11-fct*(t33-Tz)
   Tyy(i,j,k)=t22-fct*(t33-Tz)
   Tzz(i,j,k)=Tz
   Txz(i,j,k)=Tx
   Tyz(i,j,k)=Ty
end do
end do
end subroutine free_charac

subroutine free_extrap
integer i,j,k

if (freenode) then
do k=nk2+1,nz2
do j=nj1,nj2
do i=ni1,ni2
   Vx (i,j,k)=4.0*Vx (i,j,k-1)-6.0*Vx (i,j,k-2)+4.0*Vx (i,j,k-3)-Vx (i,j,k-4)
   Vy (i,j,k)=4.0*Vy (i,j,k-1)-6.0*Vy (i,j,k-2)+4.0*Vy (i,j,k-3)-Vy (i,j,k-4)
   Vz (i,j,k)=4.0*Vz (i,j,k-1)-6.0*Vz (i,j,k-2)+4.0*Vz (i,j,k-3)-Vz (i,j,k-4)
   Txx(i,j,k)=4.0*Txx(i,j,k-1)-6.0*Txx(i,j,k-2)+4.0*Txx(i,j,k-3)-Txx(i,j,k-4)
   Tzz(i,j,k)=4.0*Tzz(i,j,k-1)-6.0*Tzz(i,j,k-2)+4.0*Tzz(i,j,k-3)-Tzz(i,j,k-4)
   Tyy(i,j,k)=4.0*Tyy(i,j,k-1)-6.0*Tyy(i,j,k-2)+4.0*Tyy(i,j,k-3)-Tyy(i,j,k-4)
   Txy(i,j,k)=4.0*Txy(i,j,k-1)-6.0*Txy(i,j,k-2)+4.0*Txy(i,j,k-3)-Txy(i,j,k-4)
   Txz(i,j,k)=4.0*Txz(i,j,k-1)-6.0*Txz(i,j,k-2)+4.0*Txz(i,j,k-3)-Txz(i,j,k-4)
   Tyz(i,j,k)=4.0*Tyz(i,j,k-1)-6.0*Tyz(i,j,k-2)+4.0*Tyz(i,j,k-3)-Tyz(i,j,k-4)
end do
end do
end do
end if
end subroutine free_extrap
subroutine free_vext
integer i,j,k

if (freenode) then
do k=nk2+1,nz2
do j=nj1,nj2
do i=ni1,ni2
   Vx (i,j,k)=4.0*Vx (i,j,k-1)-6.0*Vx (i,j,k-2)+4.0*Vx (i,j,k-3)-Vx (i,j,k-4)
   Vy (i,j,k)=4.0*Vy (i,j,k-1)-6.0*Vy (i,j,k-2)+4.0*Vy (i,j,k-3)-Vy (i,j,k-4)
   Vz (i,j,k)=4.0*Vz (i,j,k-1)-6.0*Vz (i,j,k-2)+4.0*Vz (i,j,k-3)-Vz (i,j,k-4)
end do
end do
end do
end if
end subroutine free_vext
subroutine free_timg
integer i,j,k

if (freenode) then
do k=1,SEIS_GEO
do j=nj1,nj2
do i=ni1,ni2
   Tzz(i,j,nk2+k)=2.0_SP*TzSrc(i,j)-Tzz(i,j,nk2-k)
   Txz(i,j,nk2+k)=2.0_SP*TxSrc(i,j)-Txz(i,j,nk2-k)
   Tyz(i,j,nk2+k)=2.0_SP*TySrc(i,j)-Tyz(i,j,nk2-k)
end do
end do
end do

do j=nj1,nj2
do i=ni1,ni2
   Tzz(i,j,nk2)=TzSrc(i,j)
   Txz(i,j,nk2)=TxSrc(i,j)
   Tyz(i,j,nk2)=TySrc(i,j)
end do
end do
end if
end subroutine free_timg

!--------------------------------------------------------------------}
subroutine macdrp_mesg_init
  call mesg_init_LxF
  call mesg_init_LxB
  call mesg_init_LyF
  call mesg_init_LyB
  call mesg_init_LzF
  call mesg_init_LzB
end subroutine macdrp_mesg_init

subroutine mesg_init_LxF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufX1,NBufXL,SEISMPI_DATATYPE,neigid(1,1),1131,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevX2,NBufXL,SEISMPI_DATATYPE,neigid(1,2),1131,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufX2,NBufXS,SEISMPI_DATATYPE,neigid(1,2),1211,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevX1,NBufXS,SEISMPI_DATATYPE,neigid(1,1),1211,SWMPI_COMM,r2,ierr)
reqXF=(/ s1,r1,s2,r2 /)
#else
! --- LxF ------------------------------------------------------------------
! to X1
call MPI_SEND_INIT(Txx(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1131,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1132,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1135,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1139,SWMPI_COMM,s9,ierr)
! from X2
call MPI_RECV_INIT(Txx(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1131,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1132,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1135,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1139,SWMPI_COMM,r9,ierr)
! put into array
!reqXF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXF(1:12)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)

!to X2
call MPI_SEND_INIT(Txx(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1211,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1212,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1215,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1219,SWMPI_COMM,s9,ierr)
!from X1
call MPI_RECV_INIT(Txx(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1211,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1212,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1215,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1219,SWMPI_COMM,r9,ierr)
! put into array
!reqXF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXF(13:24)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
#endif
end subroutine mesg_init_LxF

subroutine mesg_init_LxB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufX1,NBufXS,SEISMPI_DATATYPE,neigid(1,1),1111,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevX2,NBufXS,SEISMPI_DATATYPE,neigid(1,2),1111,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufX2,NBufXL,SEISMPI_DATATYPE,neigid(1,2),1231,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevX1,NBufXL,SEISMPI_DATATYPE,neigid(1,1),1231,SWMPI_COMM,r2,ierr)
reqXB=(/ s1,r1,s2,r2 /)
#else
! --- LxB ------------------------------------------------------------------
! to X1
call MPI_SEND_INIT(Txx(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1111,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1112,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1115,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1119,SWMPI_COMM,s9,ierr)
! from X2
call MPI_RECV_INIT(Txx(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1111,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1112,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1115,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1119,SWMPI_COMM,r9,ierr)
! put into array
!reqXB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXB(1:12)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)

! to X2
call MPI_SEND_INIT(Txx(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1231,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1232,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1235,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1239,SWMPI_COMM,s9,ierr)
! from  X1
call MPI_RECV_INIT(Txx(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1231,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1232,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1235,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1239,SWMPI_COMM,r9,ierr)
! put into array
!reqXB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXB(13:24)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
#endif
end subroutine mesg_init_LxB

subroutine mesg_init_LyF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufY1,NBufYL,SEISMPI_DATATYPE,neigid(2,1),2131,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevY2,NBufYL,SEISMPI_DATATYPE,neigid(2,2),2131,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufY2,NBufYS,SEISMPI_DATATYPE,neigid(2,2),2211,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevY1,NBufYS,SEISMPI_DATATYPE,neigid(2,1),2211,SWMPI_COMM,r2,ierr)
reqYF=(/ s1,r1,s2,r2 /)
#else
! --- LyF ------------------------------------------------------------------
! to Y1
!call MPI_SEND_INIT(Txx(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2131,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2132,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2134,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2139,SWMPI_COMM,s9,ierr)
! from Y2
!call MPI_RECV_INIT(Txx(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2131,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2132,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2134,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2139,SWMPI_COMM,r9,ierr)
! put into array
!reqYF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYF(1:12)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
!-----

! to Y2
!call MPI_SEND_INIT(Txx(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2211,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2212,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2214,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2219,SWMPI_COMM,s9,ierr)
! from Y1
!call MPI_RECV_INIT(Txx(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2211,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2212,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2214,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2219,SWMPI_COMM,r9,ierr)
! put into array
!reqYF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYF(13:24)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
#endif
end subroutine mesg_init_LyF

subroutine mesg_init_LyB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufY1,NBufYS,SEISMPI_DATATYPE,neigid(2,1),2111,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevY2,NBufYS,SEISMPI_DATATYPE,neigid(2,2),2111,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufY2,NBufYL,SEISMPI_DATATYPE,neigid(2,2),2231,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevY1,NBufYL,SEISMPI_DATATYPE,neigid(2,1),2231,SWMPI_COMM,r2,ierr)
reqYB=(/ s1,r1,s2,r2 /)
#else
! --- LyB ------------------------------------------------------------------
! to Y1
!call MPI_SEND_INIT(Txx(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2111,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2112,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2114,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2119,SWMPI_COMM,s9,ierr)
! from Y2
!call MPI_RECV_INIT(Txx(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2111,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2112,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2114,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2119,SWMPI_COMM,r9,ierr)
! put into array
!reqYB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYB(1:12)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)

! to Y2
!call MPI_SEND_INIT(Txx(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2231,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2232,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2234,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2239,SWMPI_COMM,s9,ierr)
! from Y1
!call MPI_RECV_INIT(Txx(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2231,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2232,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2234,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2239,SWMPI_COMM,r9,ierr)
! put into array
!reqYB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYB(13:24)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
#endif
end subroutine mesg_init_LyB

subroutine mesg_init_LzF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufZ1,NBufZL,SEISMPI_DATATYPE,neigid(3,1),3131,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevZ2,NBufZL,SEISMPI_DATATYPE,neigid(3,2),3131,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufZ2,NBufZS,SEISMPI_DATATYPE,neigid(3,2),3211,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevZ1,NBufZS,SEISMPI_DATATYPE,neigid(3,1),3211,SWMPI_COMM,r2,ierr)
reqZF=(/ s1,r1,s2,r2 /)
#else
! --- LzF ------------------------------------------------------------------
!to Z1
!call MPI_SEND_INIT(Txx(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3131,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3132,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3133,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3139,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3131,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3132,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3133,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3139,SWMPI_COMM,r9,ierr)
! put into array
!reqZF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZF(1:12)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to Z2
!call MPI_SEND_INIT(Txx(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3211,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3212,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3213,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3219,SWMPI_COMM,s9,ierr)
!from Z1
!call MPI_RECV_INIT(Txx(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3211,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3212,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3213,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3219,SWMPI_COMM,r9,ierr)
! put into array
!reqZF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZF(13:24)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
#endif
end subroutine mesg_init_LzF

subroutine mesg_init_LzB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufZ1,NBufZS,SEISMPI_DATATYPE,neigid(3,1),3111,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevZ2,NBufZS,SEISMPI_DATATYPE,neigid(3,2),3111,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufZ2,NBufZL,SEISMPI_DATATYPE,neigid(3,2),3231,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevZ1,NBufZL,SEISMPI_DATATYPE,neigid(3,1),3231,SWMPI_COMM,r2,ierr)
reqZB=(/ s1,r1,s2,r2 /)
#else
! --- LzB ------------------------------------------------------------------
!to Z1
!call MPI_SEND_INIT(Txx(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3111,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3112,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3113,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3119,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3111,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3112,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3113,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3119,SWMPI_COMM,r9,ierr)
! put into array
!reqZB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZB(1:12)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to Z2
!call MPI_SEND_INIT(Txx(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3231,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3232,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3233,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3239,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3231,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3232,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3233,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3239,SWMPI_COMM,r9,ierr)
! put into array
!reqZB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZB(13:24)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
#endif
end subroutine mesg_init_LzB

#ifdef MPIBuffered
subroutine fill_buff_LxF
integer :: i,j,k,n
! to X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX1(n)=Tyy(i,j,k)
   !n=n+1; BufX1(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)=Txz(i,j,k)
   !n=n+1; BufX1(n)=Tyz(i,j,k)
end do
end do
end do
! to X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX2(n)=Tyy(i,j,k)
   !n=n+1; BufX2(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)=Txz(i,j,k)
   !n=n+1; BufX2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LxF
subroutine recv_buff_LxF
integer :: i,j,k,n
! from X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1;  Vx(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1;  Vy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1;  Vz(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1; Txx(i,j,k)=RevX1(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX1(n)
   !n=n+1; Tzz(i,j,k)=RevX1(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1; Txy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1; Txz(i,j,k)=RevX1(n)
   !n=n+1; Tyz(i,j,k)=RevX1(n)
end do
end do
end do
! from X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1;  Vx(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1;  Vy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1;  Vz(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1; Txx(i,j,k)=RevX2(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX2(n)
   !n=n+1; Tzz(i,j,k)=RevX2(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1; Txy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1; Txz(i,j,k)=RevX2(n)
   !n=n+1; Tyz(i,j,k)=RevX2(n)
end do
end do
end do
end subroutine recv_buff_LxF

subroutine fill_buff_LxB
integer :: i,j,k,n
! to X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX1(n)=Tyy(i,j,k)
   !n=n+1; BufX1(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)=Txz(i,j,k)
   !n=n+1; BufX1(n)=Tyz(i,j,k)
end do
end do
end do
! to X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX2(n)=Tyy(i,j,k)
   !n=n+1; BufX2(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)=Txz(i,j,k)
   !n=n+1; BufX2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LxB
subroutine recv_buff_LxB
integer :: i,j,k,n
! from X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1;  Vx(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1;  Vy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1;  Vz(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1; Txx(i,j,k)=RevX1(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX1(n)
   !n=n+1; Tzz(i,j,k)=RevX1(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1; Txy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1; Txz(i,j,k)=RevX1(n)
   !n=n+1; Tyz(i,j,k)=RevX1(n)
end do
end do
end do
! from X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1;  Vx(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1;  Vy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1;  Vz(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1; Txx(i,j,k)=RevX2(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX2(n)
   !n=n+1; Tzz(i,j,k)=RevX2(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1; Txy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1; Txz(i,j,k)=RevX2(n)
   !n=n+1; Tyz(i,j,k)=RevX2(n)
end do
end do
end do
end subroutine recv_buff_LxB
!---------------------------------------------------
subroutine fill_buff_LyF
integer :: i,j,k,n
! to Y1
n=0
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vz(i,j,k)
end do
end do
end do
   !n=n+1; BufY1(n)=Txx(i,j,k)
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)=Tyy(i,j,k)
end do
end do
end do
   !n=n+1; BufY1(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)=Txy(i,j,k)
end do
end do
end do
   !n=n+1; BufY1(n)=Txz(i,j,k)
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)=Tyz(i,j,k)
end do
end do
end do
! to Y2
n=0
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vz(i,j,k)
end do
end do
end do
   !n=n+1; BufY2(n)=Txx(i,j,k)
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)=Tyy(i,j,k)
end do
end do
end do
   !n=n+1; BufY2(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)=Txy(i,j,k)
end do
end do
end do
   !n=n+1; BufY2(n)=Txz(i,j,k)
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LyF
subroutine recv_buff_LyF
integer :: i,j,k,n
! from Y1
n=0
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY1(n)
   n=n+1; Tyy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY1(n)
   n=n+1; Txy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY1(n)
   n=n+1; Tyz(i,j,k)=RevY1(n)
end do
end do
end do
! from Y2
n=0
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY2(n)
   n=n+1; Tyy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY2(n)
   n=n+1; Txy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY2(n)
   n=n+1; Tyz(i,j,k)=RevY2(n)
end do
end do
end do
end subroutine recv_buff_LyF
subroutine fill_buff_LyB
integer :: i,j,k,n
! to Y1
n=0
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   !n=n+1; BufY1(n)=Txx(i,j,k)
   n=n+1; BufY1(n)=Tyy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   !n=n+1; BufY1(n)=Tzz(i,j,k)
   n=n+1; BufY1(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   !n=n+1; BufY1(n)=Txz(i,j,k)
   n=n+1; BufY1(n)=Tyz(i,j,k)
end do
end do
end do
! to Y2
n=0
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   !n=n+1; BufY2(n)=Txx(i,j,k)
   n=n+1; BufY2(n)=Tyy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   !n=n+1; BufY2(n)=Tzz(i,j,k)
   n=n+1; BufY2(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   !n=n+1; BufY2(n)=Txz(i,j,k)
   n=n+1; BufY2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LyB
subroutine recv_buff_LyB
integer :: i,j,k,n
! from Y1
n=0
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY1(n)
   n=n+1; Tyy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY1(n)
   n=n+1; Txy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY1(n)
   n=n+1; Tyz(i,j,k)=RevY1(n)
end do
end do
end do
! from Y2
n=0
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY2(n)
   n=n+1; Tyy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY2(n)
   n=n+1; Txy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY2(n)
   n=n+1; Tyz(i,j,k)=RevY2(n)
end do
end do
end do
end subroutine recv_buff_LyB
!---------------------------------------------------
subroutine fill_buff_LzF
integer :: i,j,k,n
! to Z1
n=0
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txx(i,j,k)
   !n=n+1; BufZ1(n)=Tyy(i,j,k)
   n=n+1; BufZ1(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txy(i,j,k)
   n=n+1; BufZ1(n)=Txz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)=Tyz(i,j,k)
end do
end do
end do
! to Z2
n=0
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txx(i,j,k)
   !n=n+1; BufZ2(n)=Tyy(i,j,k)
   n=n+1; BufZ2(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txy(i,j,k)
   n=n+1; BufZ2(n)=Txz(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LzF
subroutine recv_buff_LzF
integer :: i,j,k,n
! from Z1
n=0
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ1(n)
   !n=n+1; Tyy(i,j,k)=RevZ1(n)
   n=n+1; Tzz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ1(n)
   n=n+1; Txz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ1(n)
end do
end do
end do
! from Z2
n=0
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ2(n)
   !n=n+1; Tyy(i,j,k)=RevZ2(n)
   n=n+1; Tzz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ2(n)
   n=n+1; Txz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ2(n)
end do
end do
end do
end subroutine recv_buff_LzF

subroutine fill_buff_LzB
integer :: i,j,k,n
! to Z1
n=0
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txx(i,j,k)
   !n=n+1; BufZ1(n)=Tyy(i,j,k)
   n=n+1; BufZ1(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txy(i,j,k)
   n=n+1; BufZ1(n)=Txz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)=Tyz(i,j,k)
end do
end do
end do
! to Z2
n=0
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txx(i,j,k)
   !n=n+1; BufZ2(n)=Tyy(i,j,k)
   n=n+1; BufZ2(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txy(i,j,k)
   n=n+1; BufZ2(n)=Txz(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LzB
subroutine recv_buff_LzB
integer :: i,j,k,n
! from Z1
n=0
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ1(n)
   !n=n+1; Tyy(i,j,k)=RevZ1(n)
   n=n+1; Tzz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ1(n)
   n=n+1; Txz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ1(n)
end do
end do
end do
! from Z2
n=0
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ2(n)
   !n=n+1; Tyy(i,j,k)=RevZ2(n)
   n=n+1; Tzz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ2(n)
   n=n+1; Txz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ2(n)
end do
end do
end do
end subroutine recv_buff_LzB

#endif

end module macdrp_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
