module src_mod

! This module is used for incorporating seismic source
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:46:39 -0500 (Fri, 16 Jan 2009) $
! $Revision: 64 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod
use string_mod
use math_mod
use para_mod
use mpi_mod
use grid_mod
use media_mod
use nfseis_mod
implicit none
private

public ::           &
  src_fnm_init,     &
  src_alloc_moment, &
  src_alloc_force,  &
  src_destroy,      &
  src_import,       &
  src_fnm_get,      &
  src_stress,       &
  src_force,        &
  src_surface
public ::        &
  fun_bell,      &
  fun_ricker,    &
  fun_gauss,     &
  src_stf,       &
  stf_name2id,   &
  stf_id2name

!----------------------------------------------------------------------!

integer,parameter,public ::  &
     SIG_SVF_BELL     =100,  &
     SIG_SVF_TRIANGLE =110,  &
     SIG_SVF_GAUSS    =120,  &
     SIG_STF_RICKER   =200,  &
     SIG_STF_BSHIFT   =210,  &
     SIG_STF_GAUSS    =220,  &
     SIG_STF_BELL     =230,  &
     SIG_STF_STEP     =240,  &
     SIG_STF_DELTA    =250

integer,public :: num_force,num_moment,ntwin_force,ntwin_moment
integer,allocatable,public :: force_indx(:,:),moment_indx(:,:)
real(SP),allocatable,public ::                        &
    force_shift(:,:), force_t0(:), frcstf_time(:),frcstf_freq(:),   &
    moment_shift(:,:),moment_t0(:),momstf_time(:),momstf_freq(:)
real(SP),dimension(:,:),allocatable,public ::  &
    MomTxx,MomTyy,MomTzz,MomTxy,MomTxz,MomTyz, &
    ForceX,ForceY,ForceZ
character (len=SEIS_STRLEN),public :: fnm_src_conf, pnm_src
integer,public :: momstf_id,frcstf_id
#ifdef VERBOSE
integer fid_out
#endif

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

!*************************************************************************
!*                ------ src alloc and dealloc ------                    *
!*************************************************************************
subroutine src_alloc_moment(npt,ntw)
  integer,intent(in) :: npt,ntw
  allocate(moment_indx(SEIS_GEO,npt)); moment_indx=0
  allocate(moment_shift(SEIS_GEO,npt)); moment_shift=0.0_SP
  allocate(moment_t0(npt)); moment_t0=0.0_SP
  allocate(MomTxx(ntw,npt)); MomTxx=0.0_SP
  allocate(MomTyy(ntw,npt)); MomTyy=0.0_SP
  allocate(MomTzz(ntw,npt)); MomTzz=0.0_SP
  allocate(MomTxy(ntw,npt)); MomTxy=0.0_SP
  allocate(MomTxz(ntw,npt)); MomTxz=0.0_SP
  allocate(MomTyz(ntw,npt)); MomTyz=0.0_SP
  allocate(momstf_time(ntw)); momstf_time=0.0_SP
  allocate(momstf_freq(ntw)); momstf_freq=0.0_SP
end subroutine src_alloc_moment
subroutine src_alloc_force(npt,ntw)
  integer,intent(in) :: npt,ntw
  allocate(force_indx(SEIS_GEO,npt)); force_indx=0
  allocate(force_shift(SEIS_GEO,npt)); force_shift=0.0_SP
  allocate(force_t0(npt)); force_t0=0.0_SP
  allocate(ForceX(ntw,npt)); ForceX=0.0_SP
  allocate(ForceY(ntw,npt)); ForceY=0.0_SP
  allocate(ForceZ(ntw,npt)); ForceZ=0.0_SP
  allocate(frcstf_time(ntw)); frcstf_time=0.0_SP
  allocate(frcstf_freq(ntw)); frcstf_freq=0.0_SP
end subroutine src_alloc_force

subroutine src_destroy
  if (allocated(moment_indx)) deallocate(moment_indx)
  if (allocated(moment_shift)) deallocate(moment_shift)
  if (allocated(momstf_time)) deallocate(momstf_time)
  if (allocated(momstf_freq)) deallocate(momstf_freq)
  if (allocated(MomTxx)) deallocate(MomTxx)
  if (allocated(MomTyy)) deallocate(MomTyy)
  if (allocated(MomTzz)) deallocate(MomTzz)
  if (allocated(MomTxy)) deallocate(MomTxy)
  if (allocated(MomTxz)) deallocate(MomTxz)
  if (allocated(MomTyz)) deallocate(MomTyz)
  if (allocated(force_indx)) deallocate(force_indx)
  if (allocated(force_shift)) deallocate(force_shift)
  if (allocated(frcstf_time)) deallocate(frcstf_time)
  if (allocated(frcstf_freq)) deallocate(frcstf_freq)
  if (allocated(ForceX)) deallocate(ForceX)
  if (allocated(ForceY)) deallocate(ForceY)
  if (allocated(ForceZ)) deallocate(ForceZ)
#ifdef VERBOSE
  close(fid_out)
#endif
end subroutine src_destroy

!*************************************************************************
!*                     ------ source io ------                           *
!*************************************************************************
subroutine src_fnm_init(filenm)
character (len=*),intent(in) :: filenm
integer fid
fid=1001
open(fid,file=trim(filenm),status="old")
  call string_conf(fid,1,'SOURCE_CONF',2,fnm_src_conf)
  call string_conf(fid,1,'SOURCE_ROOT',2,pnm_src)
close(fid)
end subroutine src_fnm_init
function src_fnm_get(n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_src)//'/'//'source'//'_'//set_mpi_subfix(n_i,n_j,n_k)//'.nc'
end function src_fnm_get

subroutine src_import(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  integer n,i,j,k
  character (len=SEIS_STRLEN) :: filenm
  filenm=src_fnm_get(n_i,n_j,n_k)
  call nfseis_attget(filenm,'number_of_moment',num_moment)
  call nfseis_attget(filenm,'number_of_force',num_force)
if (num_moment>0) then
  call nfseis_diminfo(filenm,'moment_time_window',ntwin_moment)
  call src_alloc_moment(num_moment,ntwin_moment)
  call nfseis_varget(filenm,'Mxx',MomTxx,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'Myy',MomTyy,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'Mzz',MomTzz,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'Mxy',MomTxy,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'Mxz',MomTxz,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'Myz',MomTyz,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'moment_indx',moment_indx,(/1,1/),(/SEIS_GEO,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'moment_shift',moment_shift,(/1,1/),(/SEIS_GEO,num_moment/),(/1,1/))
  call nfseis_varget(filenm,'moment_start_time',moment_t0,(/1/),(/num_moment/),(/1/))
  call nfseis_varget(filenm,'moment_stf_time',momstf_time,(/1/),(/ntwin_moment/),(/1/))
  call nfseis_varget(filenm,'moment_stf_freq',momstf_freq,(/1/),(/ntwin_moment/),(/1/))
  !att
  call nfseis_attget(filenm,'moment_stf_id',momstf_id)
  !correct
  !do n=1,num_moment
  !   moment_indx(:,n)=inn_ijk(moment_indx(:,n))
  !end do
end if
if (num_force>0) then
  call nfseis_diminfo(filenm,'force_time_window',ntwin_force)
  call src_alloc_force(num_force,ntwin_force)
  call nfseis_varget(filenm,'Fx',ForceX,(/1,1/),(/ntwin_force,num_force/),(/1,1/))
  call nfseis_varget(filenm,'Fy',ForceY,(/1,1/),(/ntwin_force,num_force/),(/1,1/))
  call nfseis_varget(filenm,'Fz',ForceZ,(/1,1/),(/ntwin_force,num_force/),(/1,1/))
  call nfseis_varget(filenm,'force_indx',force_indx,(/1,1/),(/SEIS_GEO,num_force/),(/1,1/))
  call nfseis_varget(filenm,'force_shift',force_shift,(/1,1/),(/SEIS_GEO,num_force/),(/1,1/))
  call nfseis_varget(filenm,'force_start_time',force_t0,(/1/),(/num_force/),(/1/))
  call nfseis_varget(filenm,'force_stf_time',frcstf_time,(/1/),(/ntwin_force/),(/1/))
  call nfseis_varget(filenm,'force_stf_freq',frcstf_freq,(/1/),(/ntwin_force/),(/1/))
  !att
  call nfseis_attget(filenm,'force_stf_id',frcstf_id)
  !correct
  !do n=1,num_force
  !   force_indx(:,n)=inn_ijk(force_indx(:,n))
  !end do
end if
#ifdef VERBOSE
  fid_out=9010
  open(fid_out,                                                                      &
       file='log_stf_'//trim(set_mpi_subfix(thisid(1),thisid(2),thisid(3)))//'.dat', &
       status='unknown')
  write(fid_out,*) 'stf=',momstf_id,frcstf_id
  write(fid_out,*) 'ntime  time  stf'
#endif
end subroutine src_import

!*************************************************************************
!*                    ------ source coupling ------                      *
!*************************************************************************
subroutine src_stress(Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,tinc,stept)
real(SP),dimension(:,:,:),intent(inout) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer,intent(in) :: ntime
real(SP),intent(in) :: tinc,stept

integer :: i,j,k,n,m,si,sj,sk
real(SP) :: t,rate,Mxx,Mxy,Mxz,Myy,Myz,Mzz,d,V
#ifdef SrcSmooth
integer :: Li,Lj,Lk
real(SP),dimension(-LenFD:LenFD,-LenFD:LenFD,-LenFD:LenFD) :: normd
real(SP) :: x0,y0,z0
#endif

if ( num_moment<=0 ) return

t=(ntime+tinc)*stept

do n=1,num_moment
   si=moment_indx(1,n); sj=moment_indx(2,n); sk=moment_indx(3,n)
#ifdef SrcSmooth
   x0=moment_shift(1,n); y0=moment_shift(2,n); z0=moment_shift(3,n)
   call cal_norm_delt(normd, x0,y0,z0,                            &
        real(LenFD/2.0,SP), real(LenFD/2.0,SP), real(LenFD/2.0,SP) )
   if (freenode .and. sk+LenFD>nk2)  then
      normd=normd/sum(normd(:,:,-LenFD:nk2-sk))
   end if
#endif

do m=1,ntwin_moment
   rate=src_svf(t-moment_t0(n),momstf_time(m),momstf_freq(m),momstf_id,ntime)
#ifdef VERBOSE
   if (n==1) then
      write(fid_out,'(i5,2es12.5)') ntime, t, rate
   end if
#endif
   Mxx=MomTxx(m,n); Myy=MomTyy(m,n); Mzz=MomTzz(m,n)
   Mxy=MomTxy(m,n); Mxz=MomTxz(m,n); Myz=MomTyz(m,n)
#ifdef SrcSmooth
do Lk=-LenFD,LenFD
do Lj=-LenFD,LenFD
do Li=-LenFD,LenFD
   i=Li+si; j=Lj+sj; k=Lk+sk; d=normd(Li,Lj,Lk)
#else
   i=si; j=sj; k=sk; d=1.0_SP
#endif
   if (      i<=ni2 .and. i>=ni1  &
       .and. j<=nj2 .and. j>=nj1  &
       .and. k<=nk2 .and. k>=nk1 ) then
       !V=stepx*stepy*stepz*z(k)**2*xsin(i)
       V=(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))/8.0*z(k)**2*xsin(i)
       V=d*rate/V
       Txx(i,j,k)=Txx(i,j,k)-Mxx*V
       Tyy(i,j,k)=Tyy(i,j,k)-Myy*V
       Tzz(i,j,k)=Tzz(i,j,k)-Mzz*V
       Txy(i,j,k)=Txy(i,j,k)-Mxy*V
       Txz(i,j,k)=Txz(i,j,k)-Mxz*V
       Tyz(i,j,k)=Tyz(i,j,k)-Myz*V
   end if
#ifdef SrcSmooth
end do
end do
end do  ! loop_of Li,Lj,Lk
#endif

end do ! ntwin
end do ! loop_of_npt
end subroutine src_stress

subroutine src_force(Vx,Vy,Vz,ntime,tinc,stept)
real(SP),dimension(:,:,:),intent(inout) :: Vx,Vy,Vz
integer,intent(in) :: ntime
real(SP),intent(in) :: tinc,stept

integer :: i,j,k,n,m,si,sj,sk
real(SP) :: t,disp,d,fx0,fy0,fz0,fx,fy,fz,V
#ifdef SrcSmooth
integer :: Li,Lj,Lk
real(SP),dimension(-LenFD:LenFD,-LenFD:LenFD,-LenFD:LenFD) :: normd
real(SP) :: x0,y0,z0
#endif

if ( num_force<=0 ) return

t=(ntime+tinc)*stept

do n=1,num_force
   si=force_indx(1,n); sj=force_indx(2,n); sk=force_indx(3,n)
#ifdef SrcSmooth
   x0=force_shift(1,n); y0=force_shift(2,n); z0=force_shift(3,n)
   call cal_norm_delt(normd, x0,y0,z0,                            &
        real(LenFD/2.0,SP), real(LenFD/2.0,SP), real(LenFD/2.0,SP) )
   if (freenode .and. sk+LenFD>nk2)  then
      normd=normd/sum(normd(:,:,-LenFD:nk2-sk))
   end if
#endif

do m=1,ntwin_force
   disp=src_stf(t-force_t0(n),frcstf_time(m),frcstf_freq(m),frcstf_id,ntime)
#ifdef VERBOSE
if (n==1) then
   write(fid_out,'(i5,2es12.5)') ntime, t, disp
end if
#endif
   fx0=ForceX(m,n)*disp; fy0=ForceY(m,n)*disp; fz0=ForceZ(m,n)*disp
#ifdef SrcSmooth
do Lk=-LenFD,LenFD
do Lj=-LenFD,LenFD
do Li=-LenFD,LenFD
   i=Li+si; j=Lj+sj; k=Lk+sk; d=normd(Li,Lj,Lk)
#else
   i=si; j=sj; k=sk; d=1.0_SP
#endif
   if (      i<=ni2 .and. i>=ni1  &
       .and. j<=nj2 .and. j>=nj1  &
       .and. k<=nk2 .and. k>=nk1 ) then
       V=(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*(z(k+1)-z(k-1))/8.0*z(k)**2*xsin(i)
       V=1.0/(V*rho(i,j,k))*d
       fx=fx0*V; fy=fy0*V; fz=fz0*V
#ifdef SrcSurface
       if (freenode .and. k==nk2) cycle
#else
       if (freenode .and. k==nk2) then
          fx=2.0_SP*fx;fy=2.0_SP*fy;fz=2.0_SP*fz
       end if
#endif
       Vx(i,j,k)=Vx(i,j,k)+fx
       Vy(i,j,k)=Vy(i,j,k)+fy
       Vz(i,j,k)=Vz(i,j,k)+fz
   end if !i[jk]
#ifdef SrcSmooth
end do !Li
end do !Lj
end do !Lk
#endif

end do ! ntwin
end do ! loop_of_npt
!call src_charac(t,ntime)
end subroutine src_force

!subroutine src_charac(t,ntime)
!use macdrp_mod, only : TxSrc,TySrc,TzSrc
!real(SP),intent(in) :: t
!integer,intent(in) :: ntime
!
!real(SP) :: disp,fx,fy,fz,area
!integer :: n,m,si,sj,sk,nf
!
!if (.not. freenode) return
!if (num_force<=0) return
!
!m=1 ! for ntwin
!
!do nf=1,num_force
!   si=force_indx(1,nf); sj=force_indx(2,nf); sk=force_indx(3,nf)
!   if (sk/=nk2) cycle
!   disp=src_stf(t-force_t0(n),frcstf_time(m),frcstf_freq(m),frcstf_id,ntime)
!   !area=z(sk)**2*xsin(si)*stepx*stepy*stepz
!   area=z(sk)**2*xsin(si)*(x(si+1)-x(si-1))*(y(sj+1)-y(sj-1))/4.0
!   fx=ForceX(m,nf); fy=ForceY(m,nf); fz=ForceZ(m,nf)
!   TxSrc(si,sj)=disp*fx/area
!   TySrc(si,sj)=disp*fy/area
!   TzSrc(si,sj)=disp*fz/area
!end do ! loop_of_npt_fault
!end subroutine src_charac

subroutine src_surface(ntime,tinc,stept)
use macdrp_mod, only : &
    TxSrc,TySrc,TzSrc
use macdrp_mod, only : &
    VxSrc,VySrc,VzSrc
integer,intent(in) :: ntime
real(SP),intent(in) :: tinc,stept

integer :: n,m,i,j,k
real(SP) :: t,disp,rate,S,fx0,fy0,fz0

if (.not. freenode) return

if ( num_force<=0 ) return

t=(ntime+tinc)*stept

m=1 ! for ntwin

do n=1,num_force
   i=force_indx(1,n); j=force_indx(2,n); k=force_indx(3,n)
do m=1,ntwin_force
   disp=src_stf(t-force_t0(n),frcstf_time(m),frcstf_freq(m),frcstf_id,ntime)
   rate=src_svf(t-force_t0(n),frcstf_time(m),frcstf_freq(m),frcstf_id,ntime)
   fx0=ForceX(m,n); fy0=ForceY(m,n); fz0=ForceZ(m,n)
   if (      i<=ni2 .and. i>=ni1  &
       .and. j<=nj2 .and. j>=nj1  &
       .and. k==nk2  ) then

       S=(x(i+1)-x(i-1))*(y(j+1)-y(j-1))/4.0*z(k)**2*xsin(i)

       TxSrc(i,j)=fx0/S*disp
       TySrc(i,j)=fy0/S*disp
       TzSrc(i,j)=fz0/S*disp
       VxSrc(i,j)=fx0/S*rate
       VySrc(i,j)=fy0/S*rate
       VzSrc(i,j)=fz0/S*rate
    end if
end do
end do

end subroutine src_surface

!*************************************************************************
!*                --- source signal generator --                         *
!*************************************************************************
subroutine cal_norm_delt(delt,x0,y0,z0,rx0,ry0,rz0)
  real(SP),dimension(-LenFD:LenFD,-LenFD:LenFD,-LenFD:LenFD),intent(out) :: delt
  real(SP),intent(in) :: x0,y0,z0,rx0,ry0,rz0
  real(SP) :: D1,D2,D3
  integer i,j,k
  delt=0.0_SP
  do k=-LenFD,LenFD
  do j=-LenFD,LenFD
  do i=-LenFD,LenFD
     D1=fun_gauss(i-x0,rx0,0.0_SP)
     D2=fun_gauss(j-y0,ry0,0.0_SP)
     D3=fun_gauss(k-z0,rz0,0.0_SP)
     delt(i,j,k)=D1*D2*D3
  end do
  end do
  end do
  if (sum(delt)<SEIS_ZERO) call swmpi_except('cal_norm_delt is zero')
  delt=delt/sum(delt)
end subroutine cal_norm_delt
function src_svf(t,t0,f0,flag_stf_type,ntime) result(rate)
real(SP),intent(in) :: t,t0,f0
integer,intent(in) :: flag_stf_type
integer,intent(in) :: ntime
real(SP) :: rate
select case(flag_stf_type)
case (SIG_SVF_GAUSS)
   rate=fun_gauss(t,f0,t0)
case (SIG_STF_GAUSS)
   rate= fun_gauss_deriv(t,f0,t0)
case (SIG_STF_RICKER)
   rate=fun_ricker_deriv(t,f0,t0)
case (SIG_SVF_BELL)
   rate= fun_bell(t-t0,f0)
case (SIG_STF_BELL)
   rate= fun_bell_deriv(t-t0,f0)
case (SIG_SVF_TRIANGLE)
   rate= fun_triangle(t-t0,f0)
case (SIG_STF_STEP)
   rate= fun_delta(t,ntime)
case default
   print *, "Have you told me how to generate the SVF of ", flag_stf_type
   stop 1
end select
if (abs(rate)<SEIS_ZERO) rate=0.0
end function src_svf
function src_stf(t,t0,f0,flag_stf_type,ntime) result(disp)
real(SP),intent(in) :: t,t0,f0
integer,intent(in) :: flag_stf_type
integer,intent(in) :: ntime
real(SP) :: disp
select case(flag_stf_type)
case (SIG_STF_GAUSS)
   disp=fun_gauss(t,f0,t0)
case (SIG_STF_RICKER)
   disp=fun_ricker(t,f0,t0)
case (SIG_STF_BELL)
   disp= fun_bell(t-t0,f0)
case (SIG_SVF_BELL)
   disp= fun_bell_int(t-t0,f0)
case (SIG_STF_STEP)
   disp=fun_step(t)
case (SIG_STF_DELTA)
   disp=fun_delta(t,ntime)
case default
   print *, "Have you told me how to generate the STF of ", flag_stf_type
   stop 1
end select
!if (abs(disp)<SEIS_ZERO) disp=0.0
end function src_stf

! gauss
function fun_gauss(t,a,t0) result(f)
  real(SP),intent(in) :: t,t0,a
  real(SP) :: f
  if (abs(t0)>=SEIS_ZERO .and. (t<=0.0 .or. t>=2.0*t0) ) then
     f=0.0
  else
     f=exp(-(t-t0)**2/(a**2))/(sqrt(PI)*a)
  end if
end function fun_gauss
function fun_gauss_deriv(t,a,t0) result(f)
  real(SP),intent(in) :: t,t0,a
  real(SP) :: f
  if (abs(t0)>=SEIS_ZERO .and. (t<=0.0 .or. t>=2.0*t0) ) then
     f=0.0
  else
     f=exp(-(t-t0)**2/(a**2))/(sqrt(PI)*a)*(-2*(t-t0)/a**2)
  end if
end function fun_gauss_deriv
! ricker
function fun_ricker(t,fc,t0) result (v)
    real(SP),intent(in) :: t,t0,fc
    real(SP) :: u,f0,v
    if (t<=0.0) then
       v=0.0; return
    end if
    f0=sqrt(PI)/2.0
    u=(t-t0)*2.0*PI*fc
    v=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
    !v=u*(1.5-u**2/4.0)*exp(-u**2.0/4.0)*f0*PI*fc
end function fun_ricker
function fun_ricker_deriv(t,fc,t0) result (v)
    real(SP),intent(in) :: t,t0,fc
    real(SP) :: u,f0,v
    if (t<=0.0) then
       v=0.0; return
    end if
    f0=sqrt(PI)/2.0
    u=(t-t0)*2.0*PI*fc
    !f=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
    v=u*(1.5-u**2/4.0)*exp(-u**2.0/4.0)*f0*PI*fc
end function fun_ricker_deriv
!bell
function fun_bell(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: v
  if (t>0.0 .and. t<riset) then
     v=(1.0 - cos(2*PI*t/riset))/riset
  else
     v=0.0
  end if
end function fun_bell
function fun_bell_deriv(t,riset) result(v)
  real,intent(in) :: t,riset
  real :: v
  if (t>0.0 .and. t<riset) then
     v=2.0*PI*sin(2*PI*t/riset)/riset
  else
     v=0.0
  end if
end function fun_bell_deriv
function fun_bell_int(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: v
  if (t<=0.0) then
     v=0.0
  elseif (t<riset) then
     v=t/riset - sin(2*PI*t/riset)/(2.0*pi)
  else
     v=1.0
  end if
end function fun_bell_int
!triangle
function fun_triangle(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: t0,v
  t0=riset/2.0;
  if (t>riset) then
     v=0.0
  elseif (t>t0) then
     v=2.0/t0-t/(t0**2)
  elseif (t>0.0) then
     v=t/(t0**2)
  else
     v=0.0
  end if
end function fun_triangle
function fun_triangle_int(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: t0,v
  t0=riset/2.0;
  if (t>riset) then
     v=1.0
  elseif (t>t0) then
     v=-0.5*t**2/(t0**2)+2.0*t/t0-1.0
  elseif (t>0.0) then
     v=0.5*t**2/(t0**2)
  else
     v=0.0
  end if
end function fun_triangle_int
!bshift
function fun_bshift(t,riset,t0) result(v)
    real(SP),intent(in) :: t,t0,riset
    real(SP) :: bshift,v
    bshift=riset/2.0;
    if (t>0.0 .and. t<riset) then
       v=0.5/t0/(cosh((t-bshift)/t0)**2.0)
    else
       v=0.0
    end if
end function fun_bshift
!step
function fun_step(t) result(v)
    real(SP),intent(in) :: t
    real(SP) :: v
    !if (t>=0.0 ) then
       v=1.0
    !else
    !   v=0.0
    !end if
end function fun_step
!delta
function fun_delta(t,ntime) result(v)
    real(SP),intent(in) :: t
    integer,intent(in) :: ntime
    real(SP) :: v
    !if (abs(t)<= stept ) then
    if ( ntime == 0 ) then
       v=1.0_SP/stept
    else
       v=0.0
    end if
end function fun_delta

function stf_name2id(stfname) result(id)
character (len=*),intent(in) :: stfname
integer :: id
select case (trim(stfname))
case ('gauss')
     id=SIG_STF_GAUSS
case ('gauss_int')
     id=SIG_SVF_GAUSS
case ('ricker')
     id=SIG_STF_RICKER
case ('bell')
     id=SIG_STF_BELL
case ('bell_int')
     id=SIG_SVF_BELL
case ('triangle_int')
     id=SIG_SVF_TRIANGLE
case ('B(shift)')
     id=SIG_STF_BSHIFT
case ('step')
     id=SIG_STF_STEP
case ('delta')
     id=SIG_STF_DELTA
case default
     print *, "Have you told me how to generate the STF of ", trim(stfname)
     stop 1
end select
end function stf_name2id
function stf_id2name(id) result(stfname)
integer,intent(in) :: id
character (len=SEIS_STRLEN) :: stfname
select case (id)
case (SIG_STF_GAUSS)
     stfname='gauss'
case (SIG_SVF_GAUSS)
     stfname='gauss_int'
case (SIG_STF_RICKER)
     stfname='ricker'
case (SIG_STF_BELL)
     stfname='bell'
case (SIG_SVF_BELL)
     stfname='bell_int'
case (SIG_SVF_TRIANGLE)
     stfname='triangle_int'
case (SIG_STF_BSHIFT)
     stfname='B(shift)'
case (SIG_STF_STEP)
     stfname='step'
case (SIG_STF_DELTA)
     stfname='delta'
case default
     print *, "Cann't find the stf name for ", id
     stop 1
end select
end function stf_id2name

end module src_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
