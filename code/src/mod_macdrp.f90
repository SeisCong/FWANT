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

real(SP),parameter,private :: MCA_1=-0.30874_SP,MCA0=-0.6326_SP,MCA1=1.233_SP,MCA2=-0.3334_SP,MCA3=0.04168_SP
real(SP),parameter,private :: MAC24A0=-7.0_SP/6.0_SP,MAC24A1=8.0_SP/6.0_SP,MAC24A2=-1.0_SP/6.0_SP
real(SP),parameter,private :: MAC22A0=-1.0_SP,MAC22A1=1.0_SP
real(SP),parameter :: RK2a2=1.0_SP
real(SP),parameter :: RK2b1=1.0_SP/2.0_SP,RK2b2=1.0_SP/2.0_SP
real(SP),parameter :: RK4a2=0.5_SP,RK4a3=0.5_SP,RK4a4=1.0_SP
real(SP),parameter :: RK4b1=1.0_SP/6.0_SP,RK4b2=1.0_SP/3.0_SP,RK4b3=1.0_SP/3.0_SP,RK4b4=1.0_SP/6.0_SP
real(SP),parameter :: LA0=2.0_SP/3.0_SP, LA1=1.0_SP/3.0_SP
real(SP),parameter :: RA_1=-1.0_SP/6.0_SP, RA0=-2.0_SP/3.0_SP, RA1=5.0_SP/6.0_SP

!integer,parameter,private :: NREQ=36
integer,parameter,private :: NREQ=24


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
! main
indx(:,SEIS_GEO*2+1)=(/ ni1+3,ni2-3, &
                        nj1+3,nj2-3, &
                        nk1+3,nk2-3 /)
indx(:,SEIS_GEO*2  )=(/ ni1,ni2,nj1,nj2,nk2-3+1,nk2 /) ! z2
indx(:,SEIS_GEO*2-1)=(/ ni1,ni2,nj1,nj2,nk1,nk1+3-1 /) ! z1
indx(:,1)=(/ ni1,ni2,nj1,nj1+3-1,nk1+3,nk2-3 /) ! y1
indx(:,2)=(/ ni1,ni2,nj2-3+1,nj2,nk1+3,nk2-3 /) ! y2
indx(:,3)=(/ ni1,ni1+3-1,nj1+3,nj2-3,nk1+3,nk2-3 /) ! x1
indx(:,4)=(/ ni2-3+1,ni2,nj1+3,nj2-3,nk1+3,nk2-3 /) ! x2
! rk coefficient
firRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
firRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
!secRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
!secRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
secRKa=(/ RK2a2, 0.0, 0.0, 0.0 /)
secRKb=(/ RK2b1, RK2b2, 0.0, 0.0 /)
! mat to convert V,z






end subroutine macdrp_init
subroutine macdrp_destroy
deallocate( Txx, Tyy, Txy, Vx, Vy)
deallocate(hTxx,hTyy,hTxy,hVx,hVy)
deallocate(mTxx,mTyy,mTxy,mVx,mVy)
deallocate(tTxx,tTyy,tTxy,tVx,tVy)



end subroutine macdrp_destroy
subroutine macdrp_check(ntime)
integer,intent(in) :: ntime
real(SP) :: V1,V2,V3,T11,T22,T33,T12,T13,T23,W
integer ierr

 return

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



   W=max(V1,V2,V3,T11,T22,T33,T12,T13,T23)
   if (W>=huge(1.0)) then
      print *, "Overflow error: "
      write(*,"(i5,i3.2,2(i2.2),9(es12.5,3i5))") ntime,thisid(1),thisid(2),thisid(3), &
         V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
        T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
        T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))






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




  !Txx=real(myid)

  if (freenode) call free_timg




  n=SEIS_GEO*2+1





  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxF_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxF_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do





end subroutine in_LxF_LyF_LzF
subroutine in_LxB_LyB_LzB
  integer n






  if (freenode) call free_timg




  !Txx=real(myid)
  n=SEIS_GEO*2+1
  !if (myid==0) then
  !  Vx=1.0;Vy=2.0;Vz=3.0;Txx=4.0;Tyy=5.0;Tzz=6.0;Txy=7.0;Txz=8.0;Tyz=9.0
  !end if






  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxB_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxB_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




end subroutine in_LxB_LyB_LzB
subroutine in_LxB_LyB_LzF
  integer n






  if (freenode) call free_timg




  n=SEIS_GEO*2+1






  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxB_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxB_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




end subroutine in_LxB_LyB_LzF
subroutine in_LxF_LyF_LzB
  integer n






  if (freenode) call free_timg




  n=SEIS_GEO*2+1






  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxF_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxF_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




end subroutine in_LxF_LyF_LzB
subroutine in_LxB_LyF_LzF
  integer n




  !Txx=real(myid)

  if (freenode) call free_timg




  n=SEIS_GEO*2+1






  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxB_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxB_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




end subroutine in_LxB_LyF_LzF
subroutine in_LxF_LyB_LzB
  integer n




  !Txx=real(myid)

  if (freenode) call free_timg




  n=SEIS_GEO*2+1






  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxF_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxF_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




end subroutine in_LxF_LyB_LzB
subroutine in_LxF_LyB_LzF
  integer n




  !Txx=real(myid)

  if (freenode) call free_timg




  n=SEIS_GEO*2+1






  call MPI_STARTALL(NREQ,reqXF,ierr)
  call MPI_STARTALL(NREQ,reqYB,ierr)
  call MPI_STARTALL(NREQ,reqZF,ierr)
  call LxF_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZF,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxF_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




end subroutine in_LxF_LyB_LzF
subroutine in_LxB_LyF_LzB
  integer n




  !Txx=real(myid)

  if (freenode) call free_timg




  n=SEIS_GEO*2+1






  call MPI_STARTALL(NREQ,reqXB,ierr)
  call MPI_STARTALL(NREQ,reqYF,ierr)
  call MPI_STARTALL(NREQ,reqZB,ierr)
  call LxB_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(NREQ,reqXB,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqYF,reqstat,ierr)
  call MPI_WAITALL(NREQ,reqZB,reqstat,ierr)






  do n=1,SEIS_GEO*2
  call LxB_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do




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
     (MCA_1*Txx(i-1,j,k)+MCA0*Txx(i,j,k) &
     +MCA1*Txx(i+1,j,k)+MCA2*Txx(i+2,j,k)+MCA3*Txx(i+3,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (MCA_1*Txy(i-1,j,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i+1,j,k)+MCA2*Txy(i+2,j,k)+MCA3*Txy(i+3,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (MCA_1*Txz(i-1,j,k)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i+1,j,k)+MCA2*Txz(i+2,j,k)+MCA3*Txz(i+3,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (MCA_1*Vx(i-1,j,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i+1,j,k)+MCA2*Vx(i+2,j,k)+MCA3*Vx(i+3,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (MCA_1*Vy(i-1,j,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i+1,j,k)+MCA2*Vy(i+2,j,k)+MCA3*Vy(i+3,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (MCA_1*Vz(i-1,j,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i+1,j,k)+MCA2*Vz(i+2,j,k)+MCA3*Vz(i+3,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (MCA_1*Tyy(i,j-1,k)+MCA0*Tyy(i,j,k) &
     +MCA1*Tyy(i,j+1,k)+MCA2*Tyy(i,j+2,k)+MCA3*Tyy(i,j+3,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (MCA_1*Txy(i,j-1,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i,j+1,k)+MCA2*Txy(i,j+2,k)+MCA3*Txy(i,j+3,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (MCA_1*Tyz(i,j-1,k)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j+1,k)+MCA2*Tyz(i,j+2,k)+MCA3*Tyz(i,j+3,k)) &
     )*eta_y(j)
   DyVx = (               &
     (MCA_1*Vx(i,j-1,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j+1,k)+MCA2*Vx(i,j+2,k)+MCA3*Vx(i,j+3,k))  &
     )*eta_y(j)
   DyVy = (               &
     (MCA_1*Vy(i,j-1,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j+1,k)+MCA2*Vy(i,j+2,k)+MCA3*Vy(i,j+3,k))  &
     )*eta_y(j)
   DyVz = (               &
     (MCA_1*Vz(i,j-1,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j+1,k)+MCA2*Vz(i,j+2,k)+MCA3*Vz(i,j+3,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (MCA_1*Tzz(i,j,k-1)+MCA0*Tzz(i,j,k) &
     +MCA1*Tzz(i,j,k+1)+MCA2*Tzz(i,j,k+2)+MCA3*Tzz(i,j,k+3)) &
     )*zeta_z(k)
   DzTxz = (              &
     (MCA_1*Txz(i,j,k-1)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i,j,k+1)+MCA2*Txz(i,j,k+2)+MCA3*Txz(i,j,k+3)) &
     )*zeta_z(k)
   DzTyz = (              &
     (MCA_1*Tyz(i,j,k-1)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j,k+1)+MCA2*Tyz(i,j,k+2)+MCA3*Tyz(i,j,k+3)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (MAC22A0*Vx(i,j,k)  &
     +MAC22A1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC22A0*Vy(i,j,k)  &
     +MAC22A1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC22A0*Vz(i,j,k)  &
     +MAC22A1*Vz(i,j,k+1))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (MAC24A0*Vx(i,j,k)  &
     +MAC24A1*Vx(i,j,k+1)+MAC24A2*Vx(i,j,k+2))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC24A0*Vy(i,j,k)  &
     +MAC24A1*Vy(i,j,k+1)+MAC24A2*Vy(i,j,k+2))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC24A0*Vz(i,j,k)  &
     +MAC24A1*Vz(i,j,k+1)+MAC24A2*Vz(i,j,k+2))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (MCA_1*Vx(i,j,k-1)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j,k+1)+MCA2*Vx(i,j,k+2)+MCA3*Vx(i,j,k+3))  &
     )*zeta_z(k)
   DzVy = (               &
     (MCA_1*Vy(i,j,k-1)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j,k+1)+MCA2*Vy(i,j,k+2)+MCA3*Vy(i,j,k+3))  &
     )*zeta_z(k)
   DzVz = (               &
     (MCA_1*Vz(i,j,k-1)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j,k+1)+MCA2*Vz(i,j,k+2)+MCA3*Vz(i,j,k+3))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
     (-MCA3*Txx(i-3,j,k)-MCA2*Txx(i-2,j,k) &
     -MCA1*Txx(i-1,j,k)-MCA0*Txx(i,j,k)-MCA_1*Txx(i+1,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (-MCA3*Txy(i-3,j,k)-MCA2*Txy(i-2,j,k) &
     -MCA1*Txy(i-1,j,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i+1,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (-MCA3*Txz(i-3,j,k)-MCA2*Txz(i-2,j,k) &
     -MCA1*Txz(i-1,j,k)-MCA0*Txz(i,j,k)-MCA_1*Txz(i+1,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (-MCA3*Vx(i-3,j,k)-MCA2*Vx(i-2,j,k)  &
     -MCA1*Vx(i-1,j,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i+1,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (-MCA3*Vy(i-3,j,k)-MCA2*Vy(i-2,j,k)  &
     -MCA1*Vy(i-1,j,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i+1,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (-MCA3*Vz(i-3,j,k)-MCA2*Vz(i-2,j,k)  &
     -MCA1*Vz(i-1,j,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i+1,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (-MCA3*Tyy(i,j-3,k)-MCA2*Tyy(i,j-2,k) &
     -MCA1*Tyy(i,j-1,k)-MCA0*Tyy(i,j,k)-MCA_1*Tyy(i,j+1,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (-MCA3*Txy(i,j-3,k)-MCA2*Txy(i,j-2,k) &
     -MCA1*Txy(i,j-1,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i,j+1,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (-MCA3*Tyz(i,j-3,k)-MCA2*Tyz(i,j-2,k) &
     -MCA1*Tyz(i,j-1,k)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j+1,k)) &
     )*eta_y(j)
   DyVx = (               &
     (-MCA3*Vx(i,j-3,k)-MCA2*Vx(i,j-2,k)  &
     -MCA1*Vx(i,j-1,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j+1,k))  &
     )*eta_y(j)
   DyVy = (               &
     (-MCA3*Vy(i,j-3,k)-MCA2*Vy(i,j-2,k)  &
     -MCA1*Vy(i,j-1,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j+1,k))  &
     )*eta_y(j)
   DyVz = (               &
     (-MCA3*Vz(i,j-3,k)-MCA2*Vz(i,j-2,k)  &
     -MCA1*Vz(i,j-1,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j+1,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (-MCA3*Tzz(i,j,k-3)-MCA2*Tzz(i,j,k-2) &
     -MCA1*Tzz(i,j,k-1)-MCA0*Tzz(i,j,k)-MCA_1*Tzz(i,j,k+1)) &
     )*zeta_z(k)
   DzTxz = (              &
     (-MCA3*Txz(i,j,k-3)-MCA2*Txz(i,j,k-2) &
     -MCA1*Txz(i,j,k-1)-MCA0*Txz(i,j,k)-MCA_1*Txz(i,j,k+1)) &
     )*zeta_z(k)
   DzTyz = (              &
     (-MCA3*Tyz(i,j,k-3)-MCA2*Tyz(i,j,k-2) &
     -MCA1*Tyz(i,j,k-1)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j,k+1)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (-MAC22A1*Vx(i,j,k-1)  &
     -MAC22A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC22A1*Vy(i,j,k-1)  &
     -MAC22A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC22A1*Vz(i,j,k-1)  &
     -MAC22A0*Vz(i,j,k))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (-MAC24A2*Vx(i,j,k-2)  &
     -MAC24A1*Vx(i,j,k-1)-MAC24A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC24A2*Vy(i,j,k-2)  &
     -MAC24A1*Vy(i,j,k-1)-MAC24A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC24A2*Vz(i,j,k-2)  &
     -MAC24A1*Vz(i,j,k-1)-MAC24A0*Vz(i,j,k))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (-MCA3*Vx(i,j,k-3)-MCA2*Vx(i,j,k-2)  &
     -MCA1*Vx(i,j,k-1)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MCA3*Vy(i,j,k-3)-MCA2*Vy(i,j,k-2)  &
     -MCA1*Vy(i,j,k-1)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MCA3*Vz(i,j,k-3)-MCA2*Vz(i,j,k-2)  &
     -MCA1*Vz(i,j,k-1)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j,k+1))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
     (MCA_1*Txx(i-1,j,k)+MCA0*Txx(i,j,k) &
     +MCA1*Txx(i+1,j,k)+MCA2*Txx(i+2,j,k)+MCA3*Txx(i+3,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (MCA_1*Txy(i-1,j,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i+1,j,k)+MCA2*Txy(i+2,j,k)+MCA3*Txy(i+3,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (MCA_1*Txz(i-1,j,k)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i+1,j,k)+MCA2*Txz(i+2,j,k)+MCA3*Txz(i+3,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (MCA_1*Vx(i-1,j,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i+1,j,k)+MCA2*Vx(i+2,j,k)+MCA3*Vx(i+3,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (MCA_1*Vy(i-1,j,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i+1,j,k)+MCA2*Vy(i+2,j,k)+MCA3*Vy(i+3,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (MCA_1*Vz(i-1,j,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i+1,j,k)+MCA2*Vz(i+2,j,k)+MCA3*Vz(i+3,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (MCA_1*Tyy(i,j-1,k)+MCA0*Tyy(i,j,k) &
     +MCA1*Tyy(i,j+1,k)+MCA2*Tyy(i,j+2,k)+MCA3*Tyy(i,j+3,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (MCA_1*Txy(i,j-1,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i,j+1,k)+MCA2*Txy(i,j+2,k)+MCA3*Txy(i,j+3,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (MCA_1*Tyz(i,j-1,k)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j+1,k)+MCA2*Tyz(i,j+2,k)+MCA3*Tyz(i,j+3,k)) &
     )*eta_y(j)
   DyVx = (               &
     (MCA_1*Vx(i,j-1,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j+1,k)+MCA2*Vx(i,j+2,k)+MCA3*Vx(i,j+3,k))  &
     )*eta_y(j)
   DyVy = (               &
     (MCA_1*Vy(i,j-1,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j+1,k)+MCA2*Vy(i,j+2,k)+MCA3*Vy(i,j+3,k))  &
     )*eta_y(j)
   DyVz = (               &
     (MCA_1*Vz(i,j-1,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j+1,k)+MCA2*Vz(i,j+2,k)+MCA3*Vz(i,j+3,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (-MCA3*Tzz(i,j,k-3)-MCA2*Tzz(i,j,k-2) &
     -MCA1*Tzz(i,j,k-1)-MCA0*Tzz(i,j,k)-MCA_1*Tzz(i,j,k+1)) &
     )*zeta_z(k)
   DzTxz = (              &
     (-MCA3*Txz(i,j,k-3)-MCA2*Txz(i,j,k-2) &
     -MCA1*Txz(i,j,k-1)-MCA0*Txz(i,j,k)-MCA_1*Txz(i,j,k+1)) &
     )*zeta_z(k)
   DzTyz = (              &
     (-MCA3*Tyz(i,j,k-3)-MCA2*Tyz(i,j,k-2) &
     -MCA1*Tyz(i,j,k-1)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j,k+1)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (-MAC22A1*Vx(i,j,k-1)  &
     -MAC22A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC22A1*Vy(i,j,k-1)  &
     -MAC22A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC22A1*Vz(i,j,k-1)  &
     -MAC22A0*Vz(i,j,k))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (-MAC24A2*Vx(i,j,k-2)  &
     -MAC24A1*Vx(i,j,k-1)-MAC24A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC24A2*Vy(i,j,k-2)  &
     -MAC24A1*Vy(i,j,k-1)-MAC24A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC24A2*Vz(i,j,k-2)  &
     -MAC24A1*Vz(i,j,k-1)-MAC24A0*Vz(i,j,k))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (-MCA3*Vx(i,j,k-3)-MCA2*Vx(i,j,k-2)  &
     -MCA1*Vx(i,j,k-1)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MCA3*Vy(i,j,k-3)-MCA2*Vy(i,j,k-2)  &
     -MCA1*Vy(i,j,k-1)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MCA3*Vz(i,j,k-3)-MCA2*Vz(i,j,k-2)  &
     -MCA1*Vz(i,j,k-1)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j,k+1))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if




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
     (-MCA3*Txx(i-3,j,k)-MCA2*Txx(i-2,j,k) &
     -MCA1*Txx(i-1,j,k)-MCA0*Txx(i,j,k)-MCA_1*Txx(i+1,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (-MCA3*Txy(i-3,j,k)-MCA2*Txy(i-2,j,k) &
     -MCA1*Txy(i-1,j,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i+1,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (-MCA3*Txz(i-3,j,k)-MCA2*Txz(i-2,j,k) &
     -MCA1*Txz(i-1,j,k)-MCA0*Txz(i,j,k)-MCA_1*Txz(i+1,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (-MCA3*Vx(i-3,j,k)-MCA2*Vx(i-2,j,k)  &
     -MCA1*Vx(i-1,j,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i+1,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (-MCA3*Vy(i-3,j,k)-MCA2*Vy(i-2,j,k)  &
     -MCA1*Vy(i-1,j,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i+1,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (-MCA3*Vz(i-3,j,k)-MCA2*Vz(i-2,j,k)  &
     -MCA1*Vz(i-1,j,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i+1,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (-MCA3*Tyy(i,j-3,k)-MCA2*Tyy(i,j-2,k) &
     -MCA1*Tyy(i,j-1,k)-MCA0*Tyy(i,j,k)-MCA_1*Tyy(i,j+1,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (-MCA3*Txy(i,j-3,k)-MCA2*Txy(i,j-2,k) &
     -MCA1*Txy(i,j-1,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i,j+1,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (-MCA3*Tyz(i,j-3,k)-MCA2*Tyz(i,j-2,k) &
     -MCA1*Tyz(i,j-1,k)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j+1,k)) &
     )*eta_y(j)
   DyVx = (               &
     (-MCA3*Vx(i,j-3,k)-MCA2*Vx(i,j-2,k)  &
     -MCA1*Vx(i,j-1,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j+1,k))  &
     )*eta_y(j)
   DyVy = (               &
     (-MCA3*Vy(i,j-3,k)-MCA2*Vy(i,j-2,k)  &
     -MCA1*Vy(i,j-1,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j+1,k))  &
     )*eta_y(j)
   DyVz = (               &
     (-MCA3*Vz(i,j-3,k)-MCA2*Vz(i,j-2,k)  &
     -MCA1*Vz(i,j-1,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j+1,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (MCA_1*Tzz(i,j,k-1)+MCA0*Tzz(i,j,k) &
     +MCA1*Tzz(i,j,k+1)+MCA2*Tzz(i,j,k+2)+MCA3*Tzz(i,j,k+3)) &
     )*zeta_z(k)
   DzTxz = (              &
     (MCA_1*Txz(i,j,k-1)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i,j,k+1)+MCA2*Txz(i,j,k+2)+MCA3*Txz(i,j,k+3)) &
     )*zeta_z(k)
   DzTyz = (              &
     (MCA_1*Tyz(i,j,k-1)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j,k+1)+MCA2*Tyz(i,j,k+2)+MCA3*Tyz(i,j,k+3)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (MAC22A0*Vx(i,j,k)  &
     +MAC22A1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC22A0*Vy(i,j,k)  &
     +MAC22A1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC22A0*Vz(i,j,k)  &
     +MAC22A1*Vz(i,j,k+1))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (MAC24A0*Vx(i,j,k)  &
     +MAC24A1*Vx(i,j,k+1)+MAC24A2*Vx(i,j,k+2))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC24A0*Vy(i,j,k)  &
     +MAC24A1*Vy(i,j,k+1)+MAC24A2*Vy(i,j,k+2))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC24A0*Vz(i,j,k)  &
     +MAC24A1*Vz(i,j,k+1)+MAC24A2*Vz(i,j,k+2))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (MCA_1*Vx(i,j,k-1)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j,k+1)+MCA2*Vx(i,j,k+2)+MCA3*Vx(i,j,k+3))  &
     )*zeta_z(k)
   DzVy = (               &
     (MCA_1*Vy(i,j,k-1)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j,k+1)+MCA2*Vy(i,j,k+2)+MCA3*Vy(i,j,k+3))  &
     )*zeta_z(k)
   DzVz = (               &
     (MCA_1*Vz(i,j,k-1)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j,k+1)+MCA2*Vz(i,j,k+2)+MCA3*Vz(i,j,k+3))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
     (-MCA3*Txx(i-3,j,k)-MCA2*Txx(i-2,j,k) &
     -MCA1*Txx(i-1,j,k)-MCA0*Txx(i,j,k)-MCA_1*Txx(i+1,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (-MCA3*Txy(i-3,j,k)-MCA2*Txy(i-2,j,k) &
     -MCA1*Txy(i-1,j,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i+1,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (-MCA3*Txz(i-3,j,k)-MCA2*Txz(i-2,j,k) &
     -MCA1*Txz(i-1,j,k)-MCA0*Txz(i,j,k)-MCA_1*Txz(i+1,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (-MCA3*Vx(i-3,j,k)-MCA2*Vx(i-2,j,k)  &
     -MCA1*Vx(i-1,j,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i+1,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (-MCA3*Vy(i-3,j,k)-MCA2*Vy(i-2,j,k)  &
     -MCA1*Vy(i-1,j,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i+1,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (-MCA3*Vz(i-3,j,k)-MCA2*Vz(i-2,j,k)  &
     -MCA1*Vz(i-1,j,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i+1,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (MCA_1*Tyy(i,j-1,k)+MCA0*Tyy(i,j,k) &
     +MCA1*Tyy(i,j+1,k)+MCA2*Tyy(i,j+2,k)+MCA3*Tyy(i,j+3,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (MCA_1*Txy(i,j-1,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i,j+1,k)+MCA2*Txy(i,j+2,k)+MCA3*Txy(i,j+3,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (MCA_1*Tyz(i,j-1,k)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j+1,k)+MCA2*Tyz(i,j+2,k)+MCA3*Tyz(i,j+3,k)) &
     )*eta_y(j)
   DyVx = (               &
     (MCA_1*Vx(i,j-1,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j+1,k)+MCA2*Vx(i,j+2,k)+MCA3*Vx(i,j+3,k))  &
     )*eta_y(j)
   DyVy = (               &
     (MCA_1*Vy(i,j-1,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j+1,k)+MCA2*Vy(i,j+2,k)+MCA3*Vy(i,j+3,k))  &
     )*eta_y(j)
   DyVz = (               &
     (MCA_1*Vz(i,j-1,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j+1,k)+MCA2*Vz(i,j+2,k)+MCA3*Vz(i,j+3,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (MCA_1*Tzz(i,j,k-1)+MCA0*Tzz(i,j,k) &
     +MCA1*Tzz(i,j,k+1)+MCA2*Tzz(i,j,k+2)+MCA3*Tzz(i,j,k+3)) &
     )*zeta_z(k)
   DzTxz = (              &
     (MCA_1*Txz(i,j,k-1)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i,j,k+1)+MCA2*Txz(i,j,k+2)+MCA3*Txz(i,j,k+3)) &
     )*zeta_z(k)
   DzTyz = (              &
     (MCA_1*Tyz(i,j,k-1)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j,k+1)+MCA2*Tyz(i,j,k+2)+MCA3*Tyz(i,j,k+3)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (MAC22A0*Vx(i,j,k)  &
     +MAC22A1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC22A0*Vy(i,j,k)  &
     +MAC22A1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC22A0*Vz(i,j,k)  &
     +MAC22A1*Vz(i,j,k+1))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (MAC24A0*Vx(i,j,k)  &
     +MAC24A1*Vx(i,j,k+1)+MAC24A2*Vx(i,j,k+2))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC24A0*Vy(i,j,k)  &
     +MAC24A1*Vy(i,j,k+1)+MAC24A2*Vy(i,j,k+2))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC24A0*Vz(i,j,k)  &
     +MAC24A1*Vz(i,j,k+1)+MAC24A2*Vz(i,j,k+2))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (MCA_1*Vx(i,j,k-1)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j,k+1)+MCA2*Vx(i,j,k+2)+MCA3*Vx(i,j,k+3))  &
     )*zeta_z(k)
   DzVy = (               &
     (MCA_1*Vy(i,j,k-1)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j,k+1)+MCA2*Vy(i,j,k+2)+MCA3*Vy(i,j,k+3))  &
     )*zeta_z(k)
   DzVz = (               &
     (MCA_1*Vz(i,j,k-1)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j,k+1)+MCA2*Vz(i,j,k+2)+MCA3*Vz(i,j,k+3))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
     (MCA_1*Txx(i-1,j,k)+MCA0*Txx(i,j,k) &
     +MCA1*Txx(i+1,j,k)+MCA2*Txx(i+2,j,k)+MCA3*Txx(i+3,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (MCA_1*Txy(i-1,j,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i+1,j,k)+MCA2*Txy(i+2,j,k)+MCA3*Txy(i+3,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (MCA_1*Txz(i-1,j,k)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i+1,j,k)+MCA2*Txz(i+2,j,k)+MCA3*Txz(i+3,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (MCA_1*Vx(i-1,j,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i+1,j,k)+MCA2*Vx(i+2,j,k)+MCA3*Vx(i+3,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (MCA_1*Vy(i-1,j,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i+1,j,k)+MCA2*Vy(i+2,j,k)+MCA3*Vy(i+3,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (MCA_1*Vz(i-1,j,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i+1,j,k)+MCA2*Vz(i+2,j,k)+MCA3*Vz(i+3,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (-MCA3*Tyy(i,j-3,k)-MCA2*Tyy(i,j-2,k) &
     -MCA1*Tyy(i,j-1,k)-MCA0*Tyy(i,j,k)-MCA_1*Tyy(i,j+1,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (-MCA3*Txy(i,j-3,k)-MCA2*Txy(i,j-2,k) &
     -MCA1*Txy(i,j-1,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i,j+1,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (-MCA3*Tyz(i,j-3,k)-MCA2*Tyz(i,j-2,k) &
     -MCA1*Tyz(i,j-1,k)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j+1,k)) &
     )*eta_y(j)
   DyVx = (               &
     (-MCA3*Vx(i,j-3,k)-MCA2*Vx(i,j-2,k)  &
     -MCA1*Vx(i,j-1,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j+1,k))  &
     )*eta_y(j)
   DyVy = (               &
     (-MCA3*Vy(i,j-3,k)-MCA2*Vy(i,j-2,k)  &
     -MCA1*Vy(i,j-1,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j+1,k))  &
     )*eta_y(j)
   DyVz = (               &
     (-MCA3*Vz(i,j-3,k)-MCA2*Vz(i,j-2,k)  &
     -MCA1*Vz(i,j-1,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j+1,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (-MCA3*Tzz(i,j,k-3)-MCA2*Tzz(i,j,k-2) &
     -MCA1*Tzz(i,j,k-1)-MCA0*Tzz(i,j,k)-MCA_1*Tzz(i,j,k+1)) &
     )*zeta_z(k)
   DzTxz = (              &
     (-MCA3*Txz(i,j,k-3)-MCA2*Txz(i,j,k-2) &
     -MCA1*Txz(i,j,k-1)-MCA0*Txz(i,j,k)-MCA_1*Txz(i,j,k+1)) &
     )*zeta_z(k)
   DzTyz = (              &
     (-MCA3*Tyz(i,j,k-3)-MCA2*Tyz(i,j,k-2) &
     -MCA1*Tyz(i,j,k-1)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j,k+1)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (-MAC22A1*Vx(i,j,k-1)  &
     -MAC22A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC22A1*Vy(i,j,k-1)  &
     -MAC22A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC22A1*Vz(i,j,k-1)  &
     -MAC22A0*Vz(i,j,k))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (-MAC24A2*Vx(i,j,k-2)  &
     -MAC24A1*Vx(i,j,k-1)-MAC24A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC24A2*Vy(i,j,k-2)  &
     -MAC24A1*Vy(i,j,k-1)-MAC24A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC24A2*Vz(i,j,k-2)  &
     -MAC24A1*Vz(i,j,k-1)-MAC24A0*Vz(i,j,k))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (-MCA3*Vx(i,j,k-3)-MCA2*Vx(i,j,k-2)  &
     -MCA1*Vx(i,j,k-1)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MCA3*Vy(i,j,k-3)-MCA2*Vy(i,j,k-2)  &
     -MCA1*Vy(i,j,k-1)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MCA3*Vz(i,j,k-3)-MCA2*Vz(i,j,k-2)  &
     -MCA1*Vz(i,j,k-1)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j,k+1))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
     (MCA_1*Txx(i-1,j,k)+MCA0*Txx(i,j,k) &
     +MCA1*Txx(i+1,j,k)+MCA2*Txx(i+2,j,k)+MCA3*Txx(i+3,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (MCA_1*Txy(i-1,j,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i+1,j,k)+MCA2*Txy(i+2,j,k)+MCA3*Txy(i+3,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (MCA_1*Txz(i-1,j,k)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i+1,j,k)+MCA2*Txz(i+2,j,k)+MCA3*Txz(i+3,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (MCA_1*Vx(i-1,j,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i+1,j,k)+MCA2*Vx(i+2,j,k)+MCA3*Vx(i+3,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (MCA_1*Vy(i-1,j,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i+1,j,k)+MCA2*Vy(i+2,j,k)+MCA3*Vy(i+3,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (MCA_1*Vz(i-1,j,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i+1,j,k)+MCA2*Vz(i+2,j,k)+MCA3*Vz(i+3,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (-MCA3*Tyy(i,j-3,k)-MCA2*Tyy(i,j-2,k) &
     -MCA1*Tyy(i,j-1,k)-MCA0*Tyy(i,j,k)-MCA_1*Tyy(i,j+1,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (-MCA3*Txy(i,j-3,k)-MCA2*Txy(i,j-2,k) &
     -MCA1*Txy(i,j-1,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i,j+1,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (-MCA3*Tyz(i,j-3,k)-MCA2*Tyz(i,j-2,k) &
     -MCA1*Tyz(i,j-1,k)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j+1,k)) &
     )*eta_y(j)
   DyVx = (               &
     (-MCA3*Vx(i,j-3,k)-MCA2*Vx(i,j-2,k)  &
     -MCA1*Vx(i,j-1,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j+1,k))  &
     )*eta_y(j)
   DyVy = (               &
     (-MCA3*Vy(i,j-3,k)-MCA2*Vy(i,j-2,k)  &
     -MCA1*Vy(i,j-1,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j+1,k))  &
     )*eta_y(j)
   DyVz = (               &
     (-MCA3*Vz(i,j-3,k)-MCA2*Vz(i,j-2,k)  &
     -MCA1*Vz(i,j-1,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j+1,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (MCA_1*Tzz(i,j,k-1)+MCA0*Tzz(i,j,k) &
     +MCA1*Tzz(i,j,k+1)+MCA2*Tzz(i,j,k+2)+MCA3*Tzz(i,j,k+3)) &
     )*zeta_z(k)
   DzTxz = (              &
     (MCA_1*Txz(i,j,k-1)+MCA0*Txz(i,j,k) &
     +MCA1*Txz(i,j,k+1)+MCA2*Txz(i,j,k+2)+MCA3*Txz(i,j,k+3)) &
     )*zeta_z(k)
   DzTyz = (              &
     (MCA_1*Tyz(i,j,k-1)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j,k+1)+MCA2*Tyz(i,j,k+2)+MCA3*Tyz(i,j,k+3)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (MAC22A0*Vx(i,j,k)  &
     +MAC22A1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC22A0*Vy(i,j,k)  &
     +MAC22A1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC22A0*Vz(i,j,k)  &
     +MAC22A1*Vz(i,j,k+1))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (MAC24A0*Vx(i,j,k)  &
     +MAC24A1*Vx(i,j,k+1)+MAC24A2*Vx(i,j,k+2))  &
     )*zeta_z(k)
   DzVy = (               &
     (MAC24A0*Vy(i,j,k)  &
     +MAC24A1*Vy(i,j,k+1)+MAC24A2*Vy(i,j,k+2))  &
     )*zeta_z(k)
   DzVz = (               &
     (MAC24A0*Vz(i,j,k)  &
     +MAC24A1*Vz(i,j,k+1)+MAC24A2*Vz(i,j,k+2))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (MCA_1*Vx(i,j,k-1)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j,k+1)+MCA2*Vx(i,j,k+2)+MCA3*Vx(i,j,k+3))  &
     )*zeta_z(k)
   DzVy = (               &
     (MCA_1*Vy(i,j,k-1)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j,k+1)+MCA2*Vy(i,j,k+2)+MCA3*Vy(i,j,k+3))  &
     )*zeta_z(k)
   DzVz = (               &
     (MCA_1*Vz(i,j,k-1)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j,k+1)+MCA2*Vz(i,j,k+2)+MCA3*Vz(i,j,k+3))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
     (-MCA3*Txx(i-3,j,k)-MCA2*Txx(i-2,j,k) &
     -MCA1*Txx(i-1,j,k)-MCA0*Txx(i,j,k)-MCA_1*Txx(i+1,j,k)) &
     )*xi_x(i)
   DxTxy = (              &
     (-MCA3*Txy(i-3,j,k)-MCA2*Txy(i-2,j,k) &
     -MCA1*Txy(i-1,j,k)-MCA0*Txy(i,j,k)-MCA_1*Txy(i+1,j,k)) &
     )*xi_x(i)
   DxTxz = (              &
     (-MCA3*Txz(i-3,j,k)-MCA2*Txz(i-2,j,k) &
     -MCA1*Txz(i-1,j,k)-MCA0*Txz(i,j,k)-MCA_1*Txz(i+1,j,k)) &
     )*xi_x(i)
   DxVx =  (              &
     (-MCA3*Vx(i-3,j,k)-MCA2*Vx(i-2,j,k)  &
     -MCA1*Vx(i-1,j,k)-MCA0*Vx(i,j,k)-MCA_1*Vx(i+1,j,k))  &
     )*xi_x(i)
   DxVy = (               &
     (-MCA3*Vy(i-3,j,k)-MCA2*Vy(i-2,j,k)  &
     -MCA1*Vy(i-1,j,k)-MCA0*Vy(i,j,k)-MCA_1*Vy(i+1,j,k))  &
     )*xi_x(i)
   DxVz = (               &
     (-MCA3*Vz(i-3,j,k)-MCA2*Vz(i-2,j,k)  &
     -MCA1*Vz(i-1,j,k)-MCA0*Vz(i,j,k)-MCA_1*Vz(i+1,j,k))  &
     )*xi_x(i)

   DyTyy = (              &
     (MCA_1*Tyy(i,j-1,k)+MCA0*Tyy(i,j,k) &
     +MCA1*Tyy(i,j+1,k)+MCA2*Tyy(i,j+2,k)+MCA3*Tyy(i,j+3,k)) &
     )*eta_y(j)
   DyTxy = (              &
     (MCA_1*Txy(i,j-1,k)+MCA0*Txy(i,j,k) &
     +MCA1*Txy(i,j+1,k)+MCA2*Txy(i,j+2,k)+MCA3*Txy(i,j+3,k)) &
     )*eta_y(j)
   DyTyz = (              &
     (MCA_1*Tyz(i,j-1,k)+MCA0*Tyz(i,j,k) &
     +MCA1*Tyz(i,j+1,k)+MCA2*Tyz(i,j+2,k)+MCA3*Tyz(i,j+3,k)) &
     )*eta_y(j)
   DyVx = (               &
     (MCA_1*Vx(i,j-1,k)+MCA0*Vx(i,j,k)  &
     +MCA1*Vx(i,j+1,k)+MCA2*Vx(i,j+2,k)+MCA3*Vx(i,j+3,k))  &
     )*eta_y(j)
   DyVy = (               &
     (MCA_1*Vy(i,j-1,k)+MCA0*Vy(i,j,k)  &
     +MCA1*Vy(i,j+1,k)+MCA2*Vy(i,j+2,k)+MCA3*Vy(i,j+3,k))  &
     )*eta_y(j)
   DyVz = (               &
     (MCA_1*Vz(i,j-1,k)+MCA0*Vz(i,j,k)  &
     +MCA1*Vz(i,j+1,k)+MCA2*Vz(i,j+2,k)+MCA3*Vz(i,j+3,k))  &
     )*eta_y(j)

   DzTzz = (              &
     (-MCA3*Tzz(i,j,k-3)-MCA2*Tzz(i,j,k-2) &
     -MCA1*Tzz(i,j,k-1)-MCA0*Tzz(i,j,k)-MCA_1*Tzz(i,j,k+1)) &
     )*zeta_z(k)
   DzTxz = (              &
     (-MCA3*Txz(i,j,k-3)-MCA2*Txz(i,j,k-2) &
     -MCA1*Txz(i,j,k-1)-MCA0*Txz(i,j,k)-MCA_1*Txz(i,j,k+1)) &
     )*zeta_z(k)
   DzTyz = (              &
     (-MCA3*Tyz(i,j,k-3)-MCA2*Tyz(i,j,k-2) &
     -MCA1*Tyz(i,j,k-1)-MCA0*Tyz(i,j,k)-MCA_1*Tyz(i,j,k+1)) &
     )*zeta_z(k)

   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     (-MAC22A1*Vx(i,j,k-1)  &
     -MAC22A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC22A1*Vy(i,j,k-1)  &
     -MAC22A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC22A1*Vz(i,j,k-1)  &
     -MAC22A0*Vz(i,j,k))  &
     )*zeta_z(k)
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     (-MAC24A2*Vx(i,j,k-2)  &
     -MAC24A1*Vx(i,j,k-1)-MAC24A0*Vx(i,j,k))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MAC24A2*Vy(i,j,k-2)  &
     -MAC24A1*Vy(i,j,k-1)-MAC24A0*Vy(i,j,k))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MAC24A2*Vz(i,j,k-2)  &
     -MAC24A1*Vz(i,j,k-1)-MAC24A0*Vz(i,j,k))  &
     )*zeta_z(k)
   else

   DzVx = (               &
     (-MCA3*Vx(i,j,k-3)-MCA2*Vx(i,j,k-2)  &
     -MCA1*Vx(i,j,k-1)-MCA0*Vx(i,j,k)-MCA_1*Vx(i,j,k+1))  &
     )*zeta_z(k)
   DzVy = (               &
     (-MCA3*Vy(i,j,k-3)-MCA2*Vy(i,j,k-2)  &
     -MCA1*Vy(i,j,k-1)-MCA0*Vy(i,j,k)-MCA_1*Vy(i,j,k+1))  &
     )*zeta_z(k)
   DzVz = (               &
     (-MCA3*Vz(i,j,k-3)-MCA2*Vz(i,j,k-2)  &
     -MCA1*Vz(i,j,k-1)-MCA0*Vz(i,j,k)-MCA_1*Vz(i,j,k+1))  &
     )*zeta_z(k)

   end if



   if (miu<=SEIS_ZERO) then
      DxTxy=0.0
      DxTxz=0.0
      DyTxy=0.0
      DyTyz=0.0
      DzTxz=0.0
      DzTyz=0.0
   end if


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


   if (freenode .and. k==nk2) then
      E33=-lam/lam2mu*(E11+E22)+VzSrc(i,j)/lam2mu
      E13=VxSrc(i,j)/miu/2.0_SP
      E23=VySrc(i,j)/miu/2.0_SP
   end if



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
call MPI_SEND_INIT(Txx(ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1211,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1212,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1215,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-1+1,ny1,nz1),1,DTypeXS,neigid(1,2),1219,SWMPI_COMM,s9,ierr)
!from X1
call MPI_RECV_INIT(Txx(ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1211,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1212,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1215,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-1,ny1,nz1),1,DTypeXS,neigid(1,1),1219,SWMPI_COMM,r9,ierr)
! put into array
!reqXF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXF(13:24)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------

end subroutine mesg_init_LxF

subroutine mesg_init_LxB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9







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
call MPI_SEND_INIT(Txx(ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1231,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1232,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1235,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-3+1,ny1,nz1),1,DTypeXL,neigid(1,2),1239,SWMPI_COMM,s9,ierr)
! from  X1
call MPI_RECV_INIT(Txx(ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1231,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1232,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1235,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-3,ny1,nz1),1,DTypeXL,neigid(1,1),1239,SWMPI_COMM,r9,ierr)
! put into array
!reqXB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXB(13:24)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------

end subroutine mesg_init_LxB

subroutine mesg_init_LyF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9







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
!call MPI_SEND_INIT(Txx(nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2211,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2212,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2214,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-1+1,nz1),1,DTypeYS,neigid(2,2),2219,SWMPI_COMM,s9,ierr)
! from Y1
!call MPI_RECV_INIT(Txx(nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2211,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2212,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2214,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-1,nz1),1,DTypeYS,neigid(2,1),2219,SWMPI_COMM,r9,ierr)
! put into array
!reqYF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYF(13:24)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------

end subroutine mesg_init_LyF

subroutine mesg_init_LyB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9







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
!call MPI_SEND_INIT(Txx(nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2231,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2232,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2234,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-3+1,nz1),1,DTypeYL,neigid(2,2),2239,SWMPI_COMM,s9,ierr)
! from Y1
!call MPI_RECV_INIT(Txx(nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2231,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2232,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2234,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-3,nz1),1,DTypeYL,neigid(2,1),2239,SWMPI_COMM,r9,ierr)
! put into array
!reqYB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYB(13:24)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------

end subroutine mesg_init_LyB

subroutine mesg_init_LzF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9







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
!call MPI_SEND_INIT(Txx(nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3211,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3212,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3213,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-1+1),1,DTypeZS,neigid(3,2),3219,SWMPI_COMM,s9,ierr)
!from Z1
!call MPI_RECV_INIT(Txx(nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3211,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3212,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3213,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-1),1,DTypeZS,neigid(3,1),3219,SWMPI_COMM,r9,ierr)
! put into array
!reqZF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZF(13:24)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------

end subroutine mesg_init_LzF

subroutine mesg_init_LzB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9







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
!call MPI_SEND_INIT(Txx(nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3231,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3232,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3233,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-3+1),1,DTypeZL,neigid(3,2),3239,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3231,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3232,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3233,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-3),1,DTypeZL,neigid(3,1),3239,SWMPI_COMM,r9,ierr)
! put into array
!reqZB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZB(13:24)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

end subroutine mesg_init_LzB


end module macdrp_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
