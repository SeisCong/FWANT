module io_mod

! This module is used for input and output data
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 11:30:33 -0500 (Fri, 16 Jan 2009) $
! $Revision: 56 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

use constants_mod, only : SEIS_STRLEN,SEIS_GEO
use string_mod, only : string_conf
use para_mod
use nfseis_mod
use mpi_mod
implicit none

private
public :: io_init
public :: io_pt_read,io_pt_import,io_snap_read,io_snap_locate
public :: io_rest_export,io_rest_import
public :: io_seismo_init,io_seismo_put,io_seismo_close
public :: io_wave_export,io_wave_close
public :: io_enum,io_out_pattern
public :: io_delete
public :: get_fnm_seismo,get_fnm_station,   &
          get_fnm_snapnode,get_fnm_snapnode_n
public :: seismo_on_this,retrieve_recvline
public :: retrieve_snap_seis,retrieve_snap_time
public :: corr_subse,corr_indx
public :: read_nc_list

character (len=SEIS_STRLEN),public ::        &
    fnm_log

!---------------------------------------------
integer,public :: num_recv,num_line,num_snap,num_pt

integer,public :: pt_tinv
integer,allocatable,public :: &
     pt_indx(:,:),            &
     pt_id  (:,:),            &
     num_line_pt(:)
real(SP),allocatable,public :: &
     pt_xyz(:,:)
real(SP),public :: topo_hyper_height

integer,allocatable,public :: &
     snap_tinv(:),            &
     snap_tcnt(:),            &
     snap_gsubs(:,:),         &
     snap_gsubt(:,:),         &
     snap_gsube(:,:),         &
     snap_gsubc(:,:),         &
     snap_subs(:,:),          &
     snap_subt(:,:),          &
     snap_sube(:,:),          &
     snap_subc(:,:),          &
     sgt_ncid(:),             &
     sgt_vid(:,:),            &
     sgt_tid(:),              &
     vel_ncid(:),             &
     vel_vid(:,:),            &
     vel_tid(:)
logical,allocatable,public :: &
     snap_ishere(:),          &
     snap_oflag(:),           &
     sgt_out(:),              &
     vel_out(:)

character (len=SEIS_STRLEN),public :: &
    pnm_out,pnm_station,           &
    pnm_rest, fnm_rest, fnm_rest_point

!---------------------------------------------
integer pt_ncid,pt_tid,pt_vid(SEIS_GEO),pt_sid(6)
integer rest_tinv,run_from_rest
integer :: nt_dyn_rest, nt_dyn_sync, nt_dyn_new

!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
subroutine io_init(fnm_input)
character (len=*),intent(in) :: fnm_input
integer :: fid

fid=1001
open(fid,file=trim(fnm_input),status="old")

call string_conf(fid,1,'OUTPUT_ROOT',2,pnm_out)
call string_conf(fid,1,'STATION_ROOT',2,pnm_station)
!-- log file --
call string_conf(fid,1,'fnm_log',2,fnm_log)

!-- restart --
call string_conf(fid,1,'CHECKPOINT_ROOT',2,pnm_rest)
call string_conf(fid,1,'checkpoint_tinv',2,rest_tinv)
call string_conf(fid,1,'run_from_checkpoint',2,run_from_rest)
call string_conf(fid,1,'urgent_checkpoint',2,fnm_rest_point)
nt_dyn_rest=0; nt_dyn_sync=0;

close(fid)
end subroutine io_init

!---------------------------------------------------------------------------

subroutine io_pt_read(fnm_input)
character (len=*),intent(in) :: fnm_input
real(SP),dimension(SEIS_GEO) :: xyz0,dxyz
integer :: fid,n,m,npt

fid=1001
open(fid,file=trim(fnm_input),status="old")
call string_conf(fid,1,'number_of_recv',2,num_recv)
call string_conf(fid,1,'number_of_inline',2,num_line)
call string_conf(fid,1,'tinv_of_seismo',2,pt_tinv)
call string_conf(fid,1,'topo_hyper_height',2,topo_hyper_height)
num_pt=num_recv
do n=1,num_line
   call string_conf(fid,1,trim(io_enum('line_',n)),8,m)
   num_pt=num_pt+m
end do
call alloc_pt(num_pt,num_line)
!-- recv --
npt=0
do n=1,num_recv
   npt=npt+1
   do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('recv_',n)),m+1,pt_xyz(m,npt))
   end do
   pt_id(:,npt)=(/ n,0 /)
end do
!-- inline --
do n=1,num_line
   do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('line_',n)),m+1,xyz0(m))
      call string_conf(fid,1,trim(io_enum('line_',n)),m+1+SEIS_GEO,dxyz(m))
   end do
   call string_conf(fid,1,trim(io_enum('line_',n)),1+1+2*SEIS_GEO,num_line_pt(n))
   do m=1,num_line_pt(n)
      npt=npt+1
      pt_xyz(:,npt)=xyz0+dxyz*(m-1)
      pt_id(:,npt)=(/ m,n /)
   end do
end do
pt_xyz(1:2,:)=pt_xyz(1:2,:)*PI/180.0_SP
end subroutine io_pt_read
subroutine alloc_pt(npt,nline)
 integer,intent(in) :: npt,nline
 allocate(pt_xyz(SEIS_GEO,npt)); pt_xyz=0.0
 allocate(pt_id (       2,npt)); pt_id =0
 allocate(num_line_pt(nline));   num_line_pt=0
end subroutine alloc_pt

subroutine io_pt_import
integer :: n
character (len=SEIS_STRLEN) :: filenm
filenm= get_fnm_station(pnm_station,thisid(1),thisid(2),thisid(3))
call nfseis_diminfo(trim(filenm),'num_pt',num_pt)
call nfseis_attget(trim(filenm),'tinv',pt_tinv)
if (num_pt>0) then
   allocate(pt_indx(SEIS_GEO,num_pt))
   call nfseis_varget(trim(filenm),'indx', pt_indx, &
        (/1,1/),(/SEIS_GEO,num_pt/),(/1,1/))
!do n=1,num_pt
!   pt_indx(1,n)=inn_i(pt_indx(1,n))
!   pt_indx(2,n)=inn_j(pt_indx(2,n))
!   pt_indx(3,n)=inn_k(pt_indx(3,n))
!end do
end if
end subroutine io_pt_import

!---------------------------------------------------------------------------

subroutine io_snap_read(fnm_input)
character (len=*),intent(in) :: fnm_input
character (len=2) :: varnm
integer :: fid,n,m

fid=1001
open(fid,file=trim(fnm_input),status="old")
!-- output --
call string_conf(fid,1,'number_of_snap',2,num_snap)
call alloc_snap(num_snap)
do n=1,num_snap
do m=1,SEIS_GEO
 call string_conf(fid,1,trim(io_enum('snap_',n)),m+1,snap_subs(m,n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),m+1+SEIS_GEO,snap_subc(m,n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),m+1+2*SEIS_GEO,snap_subt(m,n))
end do
 call string_conf(fid,1,trim(io_enum('snap_',n)),1+1+3*SEIS_GEO,snap_tinv(n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),2+1+3*SEIS_GEO,snap_tcnt(n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),3+1+3*SEIS_GEO,varnm)
 sgt_out(n)=.false.; if (index(varnm,"T")/=0) sgt_out(n)=.true.
 vel_out(n)=.false.; if (index(varnm,"V")/=0) vel_out(n)=.true.
 call corr_indx(snap_subs(:,n),snap_subc(:,n),snap_subt(:,n),snap_sube(:,n))
end do
!snap_sube=snap_subs+snap_subt*(snap_subc-1)
close(fid)
end subroutine io_snap_read

subroutine alloc_snap(nsnap)
  integer,intent(in) :: nsnap
  allocate(snap_gsubs(SEIS_GEO,nsnap))
  allocate(snap_gsube(SEIS_GEO,nsnap))
  allocate(snap_gsubt(SEIS_GEO,nsnap))
  allocate(snap_gsubc(SEIS_GEO,nsnap))
  allocate(snap_subs(SEIS_GEO,nsnap))
  allocate(snap_sube(SEIS_GEO,nsnap))
  allocate(snap_subt(SEIS_GEO,nsnap))
  allocate(snap_subc(SEIS_GEO,nsnap))
  allocate(snap_tinv(nsnap))
  allocate(snap_tcnt(nsnap))
  allocate(snap_ishere(nsnap)); snap_ishere=.false.
  allocate(snap_oflag(nsnap)); snap_oflag=.false.
  allocate(vel_ncid(nsnap))
  allocate(vel_vid(SEIS_GEO,nsnap))
  allocate(vel_tid(nsnap))
  allocate(vel_out(nsnap)); vel_out=.false.
  allocate(sgt_ncid(nsnap))
  allocate(sgt_vid(6,nsnap))
  allocate(sgt_tid(nsnap))
  allocate(sgt_out(nsnap)); sgt_out=.false.
end subroutine alloc_snap

subroutine io_snap_locate(n_i,n_j,n_k)
integer,intent(in) :: n_i,n_j,n_K
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
integer,dimension(SEIS_GEO) :: gsubs,gsube
integer n

do n=1,num_snap
   subs=snap_subs(:,n);subc=snap_subc(:,n);subt=snap_subt(:,n)
   sube=snap_sube(:,n)
   ! convert into this thread
if (      subs(1)<=ngi2 .and. sube(1)>=ngi1 &
    .and. subs(2)<=ngj2 .and. sube(2)>=ngj1 &
    .and. subs(3)<=ngk2 .and. sube(3)>=ngk1 ) then
   snap_ishere(n)=.true.
   call corr_subse(subs,subc,subt,sube)
   gsubs(1)=out_i(subs(1)); gsubs(2)=out_j(subs(2)); gsubs(3)=out_k(subs(3))
   gsube(1)=out_i(sube(1)); gsube(2)=out_j(sube(2)); gsube(3)=out_k(sube(3))
   snap_gsubs(:,n)=gsubs
   snap_gsube(:,n)=gsube
   snap_gsubc(:,n)=subc
   snap_subs(1,n)=swmpi_locli(subs(1),n_i)
   snap_subs(2,n)=swmpi_loclj(subs(2),n_j)
   snap_subs(3,n)=swmpi_loclk(subs(3),n_k)
   snap_sube(1,n)=swmpi_locli(sube(1),n_i)
   snap_sube(2,n)=swmpi_loclj(sube(2),n_j)
   snap_sube(3,n)=swmpi_loclk(sube(3),n_k)
   snap_subc(:,n)=subc
end if
end do
end subroutine io_snap_locate
!---------------------------------------------------------------------------

subroutine corr_indx(subs,subc,subt,sube)
integer,dimension(SEIS_GEO),intent(inout) :: subs,subc,subt,sube
   ! -1 global first index
   ! 0 index of source center
   ! -2 global index of free surface
   if (subs(1)==-1) then
       subs(1)=ni1
   else
       subs(1)=inn_i(subs(1))
   end if

   if (subs(2)==-1) then
       subs(2)=nj1
   else
       subs(2)=inn_j(subs(2))
   end if

   if (subs(3)==-2) then
      subs(3)=swmpi_globk(nk2,dims(3)-1)
   elseif (subs(3)==-1) then
      subs(3)=nk1
   else
      subs(3)=inn_k(subs(3))
   end if

   sube=subs+subt*(subc-1)
   if (subc(1)==-1) then
      sube(1)=swmpi_globi(ni2,dims(1)-1)
      subc(1)=(sube(1)-subs(1))/subt(1)+1
   end if
   if (subc(2)==-1) then
      sube(2)=swmpi_globj(nj2,dims(2)-1)
      subc(2)=(sube(2)-subs(2))/subt(2)+1
   end if
   if (subc(3)==-1) then
      sube(3)=swmpi_globk(nk2,dims(3)-1)
      subc(3)=(sube(3)-subs(3))/subt(3)+1
   end if
end subroutine corr_indx
subroutine corr_subse(subs,subc,subt,sube)
! in global mode
integer,dimension(SEIS_GEO),intent(inout) :: subs,subc,subt,sube
integer n
 if (ngi1>subs(1)) then
    n=ngi1-subs(1)
    if (mod(n,subt(1))==0) then
       subs(1)=ngi1
    else
       n=(n+subt(1)-1)/subt(1)
       subs(1)=n*subt(1)+subs(1)
    end if
 end if
 if (ngj1>subs(2)) then
    n=ngj1-subs(2)
    if (mod(n,subt(2))==0) then
       subs(2)=ngj1
    else
       n=(n+subt(2)-1)/subt(2)
       subs(2)=n*subt(2)+subs(2)
    end if
 end if
 if (ngk1>subs(3)) then
    n=ngk1-subs(3)
    if (mod(n,subt(3))==0) then
       subs(3)=ngk1
   else
       n=(n+subt(3)-1)/subt(3)
       subs(3)=n*subt(3)+subs(3)
   end if
 end if
 subc(1)=(ngi2-subs(1))/subt(1)
 subc(2)=(ngj2-subs(2))/subt(2)
 subc(3)=(ngk2-subs(3))/subt(3)
 sube(1)=min(sube(1),subc(1)*subt(1)+subs(1))
 sube(2)=min(sube(2),subc(2)*subt(2)+subs(2))
 sube(3)=min(sube(3),subc(3)*subt(3)+subs(3))
 subc=(sube-subs)/subt+1
end subroutine corr_subse
!---------------------------------------------------------------------------
subroutine io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
use mpi
integer,intent(in) :: ntime
real(SP),dimension(:,:,:),intent(in) :: Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz
integer,dimension(SEIS_GEO) :: subs,subc,subt
character (len=SEIS_STRLEN) :: filenm
integer :: fid,n,ierr,ncid,vid(3),sid(6)
  fid=5001
if (mod(ntime,10)==0) then
if (masternode) then
  open(fid,file=trim(fnm_rest_point),status='old',iostat=ierr)
  if (ierr/=0)  &
    call swmpi_except("io_rest_export: error when open "//trim(fnm_rest_point))
  read(fid,*) nt_dyn_rest, nt_dyn_sync,nt_dyn_new
  close(fid)
end if
  call MPI_BCAST(nt_dyn_rest,1,MPI_INTEGER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(nt_dyn_sync,1,MPI_INTEGER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(nt_dyn_new ,1,MPI_INTEGER,0,SWMPI_COMM,ierr)
  if (nt_dyn_new>ntime) call reset_nt(nt_dyn_new)
end if

if (ntime==nt_dyn_rest .or. ntime==nt_dyn_sync .or. mod(ntime,rest_tinv)==0) then
  if (num_pt>0) ierr=nf90_sync(pt_ncid)
  do n=1,num_snap
     if (snap_oflag(n)) then
        if (vel_out(n)) ierr=nf90_sync(vel_ncid(n))
        if (sgt_out(n)) ierr=nf90_sync(sgt_ncid(n))
     end if
  end do
end if

if (ntime/=nt_dyn_rest .and. mod(ntime,rest_tinv)/=0) return

  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=get_fnm_rest(pnm_rest,ntime,thisid(1),thisid(2),thisid(3))
  call nfseis_grid3d_def(filenm,nx,ny,nz,ncid, &
       "Restart file generated by seis3d_wave" )
  call nfseis_grid3d_defvar(ncid,'Vx' ,vid(1))
  call nfseis_grid3d_defvar(ncid,'Vy' ,vid(2))
  call nfseis_grid3d_defvar(ncid,'Vz' ,vid(3))
  call nfseis_grid3d_defvar(ncid,'Txx',sid(1))
  call nfseis_grid3d_defvar(ncid,'Tyy',sid(2))
  call nfseis_grid3d_defvar(ncid,'Tzz',sid(3))
  call nfseis_grid3d_defvar(ncid,'Txy',sid(4))
  call nfseis_grid3d_defvar(ncid,'Txz',sid(5))
  call nfseis_grid3d_defvar(ncid,'Tyz',sid(6))
  call nfseis_grid3d_enddef(ncid)
  call nfseis_put(ncid,vid(1),Vx, subs,subc,subt)
  call nfseis_put(ncid,vid(2),Vy, subs,subc,subt)
  call nfseis_put(ncid,vid(3),Vz, subs,subc,subt)
  call nfseis_put(ncid,sid(1),Txx,subs,subc,subt)
  call nfseis_put(ncid,sid(2),Tyy,subs,subc,subt)
  call nfseis_put(ncid,sid(3),Tzz,subs,subc,subt)
  call nfseis_put(ncid,sid(4),Txy,subs,subc,subt)
  call nfseis_put(ncid,sid(5),Txz,subs,subc,subt)
  call nfseis_put(ncid,sid(6),Tyz,subs,subc,subt)
  call nfseis_close(ncid)

  if (ntime-rest_tinv>0) then
     filenm=get_fnm_rest(pnm_rest,ntime-rest_tinv,thisid(1),thisid(2),thisid(3))
     call io_delete(filenm)
  end if
end subroutine io_rest_export
subroutine io_rest_import(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
integer,intent(out) :: ntime
real(SP),dimension(:,:,:),intent(out) :: Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz
integer,dimension(SEIS_GEO) :: subs,subc,subt
character (len=SEIS_STRLEN) :: filenm
if (run_from_rest==0) then
  ntime=0
else
  ntime=run_from_rest
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=get_fnm_rest(pnm_rest,ntime,thisid(1),thisid(2),thisid(3))
  call nfseis_varget( filenm,'Vx', Vx, subs,subc,subt)
  call nfseis_varget( filenm,'Vy', Vy, subs,subc,subt)
  call nfseis_varget( filenm,'Vz', Vz, subs,subc,subt)
  call nfseis_varget( filenm,'Txx',Txx,subs,subc,subt)
  call nfseis_varget( filenm,'Tyy',Tyy,subs,subc,subt)
  call nfseis_varget( filenm,'Tzz',Tzz,subs,subc,subt)
  call nfseis_varget( filenm,'Txy',Txy,subs,subc,subt)
  call nfseis_varget( filenm,'Txz',Txz,subs,subc,subt)
  call nfseis_varget( filenm,'Tyz',Tyz,subs,subc,subt)
end if
end subroutine io_rest_import
!---------------------------------------------------------------------------

subroutine io_seismo_init
character (len=SEIS_STRLEN) :: filenm
if (num_pt<=0) return
filenm=get_fnm_seismo(pnm_out,thisid(1),thisid(2),thisid(3))
if (run_from_rest==0) then
   call nfseis_seismo_def(filenm,num_pt,pt_ncid,pt_tid,title="seismograms")
   call nfseis_seismo_defvar(pt_ncid,'Vx',pt_vid(1))
   call nfseis_seismo_defvar(pt_ncid,'Vy',pt_vid(2))
   call nfseis_seismo_defvar(pt_ncid,'Vz',pt_vid(3))
   call nfseis_seismo_defvar(pt_ncid,'Txx',pt_sid(1))
   call nfseis_seismo_defvar(pt_ncid,'Tyy',pt_sid(2))
   call nfseis_seismo_defvar(pt_ncid,'Tzz',pt_sid(3))
   call nfseis_seismo_defvar(pt_ncid,'Txy',pt_sid(4))
   call nfseis_seismo_defvar(pt_ncid,'Txz',pt_sid(5))
   call nfseis_seismo_defvar(pt_ncid,'Tyz',pt_sid(6))
   call nfseis_seismo_enddef(pt_ncid)
else
   call nfseis_open(filenm,pt_ncid)
   call nfseis_inq_varid(pt_ncid,'time',pt_tid)
   call nfseis_inq_varid(pt_ncid,'Vx',pt_vid(1))
   call nfseis_inq_varid(pt_ncid,'Vy',pt_vid(2))
   call nfseis_inq_varid(pt_ncid,'Vz',pt_vid(3))
   call nfseis_inq_varid(pt_ncid,'Txx',pt_sid(1))
   call nfseis_inq_varid(pt_ncid,'Tyy',pt_sid(2))
   call nfseis_inq_varid(pt_ncid,'Tzz',pt_sid(3))
   call nfseis_inq_varid(pt_ncid,'Txy',pt_sid(4))
   call nfseis_inq_varid(pt_ncid,'Txz',pt_sid(5))
   call nfseis_inq_varid(pt_ncid,'Tyz',pt_sid(6))
end if
end subroutine io_seismo_init

subroutine io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
real(SP),dimension(:,:,:),intent(in) :: Vx,Vy,Vz
real(SP),dimension(:,:,:),intent(in) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer,intent(in) :: ntime
real(SP) :: t
integer i,j,k,n,m,ierr

if (num_pt<=0) return

t=ntime*stept
if (mod(ntime,pt_tinv)==0) then
m=ntime/pt_tinv
do n=1,num_pt
   i=pt_indx(1,n);j=pt_indx(2,n);k=pt_indx(3,n)
   call nfseis_put(pt_ncid,pt_tid,t,(/m/),(/1/),(/1/))
   call nfseis_put(pt_ncid,pt_vid(1),Vx(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_vid(2),Vy(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_vid(3),Vz(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_sid(1),Txx(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_sid(2),Tyy(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_sid(3),Tzz(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_sid(4),Txy(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_sid(5),Txz(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_sid(6),Tyz(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
end do
end if
if (mod(ntime,pt_tinv*100)==0) ierr=nf90_sync(pt_ncid)
end subroutine io_seismo_put

subroutine io_seismo_close
if (num_pt<=0) return
call nfseis_close(pt_ncid)
end subroutine io_seismo_close

!---------------------------------------------------------------------------
subroutine io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
real(SP),dimension(:,:,:),intent(in) :: Vx,Vy,Vz
real(SP),dimension(:,:,:),intent(in) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer,intent(in) :: ntime
real(SP),intent(in) :: stept

integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
integer,dimension(SEIS_GEO) :: isubs,isube
integer,dimension(SEIS_GEO+1) :: startput,countput,strideput
character (len=SEIS_STRLEN) :: filenm
real(SP) :: t
integer :: n,n1
t=ntime*stept
do n=1,num_snap
   if (.not. snap_ishere(n)) cycle
   if (mod(ntime,snap_tinv(n))/=0) cycle

   n1=mod(ntime/snap_tinv(n)-1,snap_tcnt(n))
   subs=snap_subs(:,n); subt=snap_subt(:,n); sube=snap_sube(:,n)
   subc=snap_subc(:,n)
   isubs(1) =out_i(subs(1));isubs(2)=out_j(subs(2));isubs(3)=out_k(subs(3))
   isube(1) =out_i(sube(1));isube(2)=out_j(sube(2));isube(3)=out_k(sube(3))

if ( n1==0 ) then
   if (vel_out(n)) then
   filenm=get_fnm_snapnode(pnm_out,'vel_',n,ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_snap_def(filenm,vel_ncid(n),vel_tid(n),stept*snap_tinv(n),     &
          subc,"snap of velocity feilds")
   call nfseis_snap_attdef(vel_ncid(n), &
          snap_gsubs(:,n), snap_gsubc(:,n), snap_gsubt(:,n), snap_gsube(:,n), &
          isubs, subc, subt, isube)
   call nfseis_snap_defvar(vel_ncid(n),'Vx',vel_vid(1,n))
   call nfseis_snap_defvar(vel_ncid(n),'Vy',vel_vid(2,n))
   call nfseis_snap_defvar(vel_ncid(n),'Vz',vel_vid(3,n))
   call nfseis_snap_enddef(vel_ncid(n))
   snap_oflag(n)=.true.
   end if
   if (sgt_out(n)) then
   filenm=get_fnm_snapnode(pnm_out,'sgt_',n,ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_snap_def(filenm,sgt_ncid(n),sgt_tid(n),stept*snap_tinv(n),     &
          subc,"snap of stress feilds")
   call nfseis_snap_attdef(sgt_ncid(n), &
          snap_gsubs(:,n), snap_gsubc(:,n), snap_gsubt(:,n), snap_gsube(:,n), &
          isubs, subc, subt, isube)
   call nfseis_snap_defvar(sgt_ncid(n),'Txx',sgt_vid(1,n))
   call nfseis_snap_defvar(sgt_ncid(n),'Tyy',sgt_vid(2,n))
   call nfseis_snap_defvar(sgt_ncid(n),'Tzz',sgt_vid(3,n))
   call nfseis_snap_defvar(sgt_ncid(n),'Txy',sgt_vid(4,n))
   call nfseis_snap_defvar(sgt_ncid(n),'Txz',sgt_vid(5,n))
   call nfseis_snap_defvar(sgt_ncid(n),'Tyz',sgt_vid(6,n))
   call nfseis_snap_enddef(sgt_ncid(n))
   snap_oflag(n)=.true.
   end if
elseif ( run_from_rest>0 .and. ntime-run_from_rest<=snap_tinv(n) ) then
   if (vel_out(n)) then
   filenm=get_fnm_snapnode(pnm_out,'vel_',n,ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_open(filenm,vel_ncid(n))
   call nfseis_inq_varid(vel_ncid(n),'Vx',vel_vid(1,n))
   call nfseis_inq_varid(vel_ncid(n),'Vy',vel_vid(2,n))
   call nfseis_inq_varid(vel_ncid(n),'Vz',vel_vid(3,n))
   call nfseis_inq_varid(vel_ncid(n),'time',vel_tid(n))
   snap_oflag(n)=.true.
   end if
   if (sgt_out(n)) then
   filenm=get_fnm_snapnode(pnm_out,'sgt_',n,ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_open(filenm,sgt_ncid(n))
   call nfseis_inq_varid(sgt_ncid(n),'Txx',sgt_vid(1,n))
   call nfseis_inq_varid(sgt_ncid(n),'Tyy',sgt_vid(2,n))
   call nfseis_inq_varid(sgt_ncid(n),'Tzz',sgt_vid(3,n))
   call nfseis_inq_varid(sgt_ncid(n),'Txy',sgt_vid(4,n))
   call nfseis_inq_varid(sgt_ncid(n),'Txz',sgt_vid(5,n))
   call nfseis_inq_varid(sgt_ncid(n),'Tyz',sgt_vid(6,n))
   call nfseis_inq_varid(sgt_ncid(n),'time',sgt_tid(n))
   snap_oflag(n)=.true.
   end if
end if

   startput=(/1,1,1,n1+1/);countput=(/subc,1/);strideput=(/1,1,1,1/)

   if (vel_out(n)) then
   call nfseis_put(vel_ncid(n),vel_tid(n),t, &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(vel_ncid(n),vel_vid(1,n), &
        Vx(subs(1):sube(1):subt(1),          &
           subs(2):sube(2):subt(2),          &
           subs(3):sube(3):subt(3)),         &
        startput,countput,strideput )
   call nfseis_put(vel_ncid(n),vel_vid(2,n), &
        Vy(subs(1):sube(1):subt(1),          &
           subs(2):sube(2):subt(2),          &
           subs(3):sube(3):subt(3)),         &
        startput,countput,strideput )
   call nfseis_put(vel_ncid(n),vel_vid(3,n), &
        Vz(subs(1):sube(1):subt(1),          &
           subs(2):sube(2):subt(2),          &
           subs(3):sube(3):subt(3)),         &
        startput,countput,strideput )
   end if
   if (sgt_out(n)) then
   call nfseis_put(sgt_ncid(n),sgt_tid(n),t, &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(sgt_ncid(n),sgt_vid(1,n), &
        Txx(subs(1):sube(1):subt(1),         &
            subs(2):sube(2):subt(2),         &
            subs(3):sube(3):subt(3)),        &
        startput,countput,strideput )
   call nfseis_put(sgt_ncid(n),sgt_vid(2,n), &
        Tyy(subs(1):sube(1):subt(1),         &
            subs(2):sube(2):subt(2),         &
            subs(3):sube(3):subt(3)),        &
        startput,countput,strideput )
   call nfseis_put(sgt_ncid(n),sgt_vid(3,n), &
        Tzz(subs(1):sube(1):subt(1),         &
            subs(2):sube(2):subt(2),         &
            subs(3):sube(3):subt(3)),        &
        startput,countput,strideput )
   call nfseis_put(sgt_ncid(n),sgt_vid(4,n), &
        Txy(subs(1):sube(1):subt(1),         &
            subs(2):sube(2):subt(2),         &
            subs(3):sube(3):subt(3)),        &
        startput,countput,strideput )
   call nfseis_put(sgt_ncid(n),sgt_vid(5,n), &
        Txz(subs(1):sube(1):subt(1),         &
            subs(2):sube(2):subt(2),         &
            subs(3):sube(3):subt(3)),        &
        startput,countput,strideput )
   call nfseis_put(sgt_ncid(n),sgt_vid(6,n), &
        Tyz(subs(1):sube(1):subt(1),         &
            subs(2):sube(2):subt(2),         &
            subs(3):sube(3):subt(3)),        &
        startput,countput,strideput )
   end if

if ( n1==snap_tcnt(n)-1 ) then
   if (vel_out(n)) call nfseis_close(vel_ncid(n))
   if (sgt_out(n)) call nfseis_close(sgt_ncid(n))
   snap_oflag(n)=.false.
end if

end do
end subroutine io_wave_export

subroutine io_wave_close
integer n
do n=1,num_snap
if (snap_oflag(n)) then
   if (vel_out(n)) call nfseis_close(vel_ncid(n))
   if (sgt_out(n)) call nfseis_close(sgt_ncid(n))
end if
end do
end subroutine io_wave_close

subroutine read_nc_list(pnm_list,fnm_list,varnm,var,subs1,subc1)
character (len=*),intent(in) :: pnm_list,fnm_list,varnm
real(SP),dimension(:,:),intent(out) :: var
integer,dimension(SEIS_GEO),intent(in) :: subs1,subc1

character (len=SEIS_STRLEN) :: fnm_nc
integer,dimension(SEIS_GEO) :: subs,subc,subt
integer fid,ierr

integer i1,i2,j1,j2,k1,k2
integer gi1,gi2,gj1,gj2,gk1,gk2
integer x1,y1,z1

gi1=subs1(1); gi2=subs1(1)+subc1(1)-1
gj1=subs1(2); gj2=subs1(2)+subc1(2)-1
gk1=subs1(3); gk2=subs1(3)+subc1(3)-1

fid=4001
open(fid,file=trim(pnm_list)//'/'//trim(fnm_list),status='old')
do
   read(fid,'(a)',iostat=ierr) fnm_nc
   if (ierr<0) exit

   fnm_nc=trim(pnm_list)//'/'//trim(fnm_nc)
   call nfseis_attget(fnm_nc,'indx',subs)
   call nfseis_diminfo(fnm_nc,'I',i2)
   call nfseis_diminfo(fnm_nc,'J',j2)
   call nfseis_diminfo(fnm_nc,'K',k2)
   i1=subs(1);j1=subs(2);k1=subs(3)
   i2=i2+i1-1;j2=j2+j1-1;k2=k2+k1-1

if (      i1<=gi2 .and. i2>=gi1 &
    .and. j1<=gj2 .and. j2>=gj1 &
    .and. k1<=gk2 .and. k2>=gk1 ) then
   x1=max(gi1,i1); y1=max(gj1,j1); z1=max(gk1,k1)

   subs(1)=x1-i1+1;
   subs(2)=y1-j1+1;
   subs(3)=z1-k1+1;
   subc(1)=min(gi2,i2)-x1+1
   subc(2)=min(gj2,j2)-y1+1
   subc(3)=min(gk2,k2)-z1+1
   subt=1
   i1=x1-gi1+1; i2=i1+subc(1)-1
   j1=y1-gj1+1; j2=j1+subc(2)-1
   k1=z1-gk1+1; k2=k1+subc(3)-1
   call nfseis_varget(fnm_nc,trim(varnm),var(i1:i2,j1:j2), &
        subs,subc,subt)
end if

end do !read
close(fid)
end subroutine read_nc_list

function seismo_on_this(pnm_sta,id,indx,n_i,n_j,n_k,npt) result(isIn)
character (len=*) :: pnm_sta
integer id,indx,n_i,n_j,n_k,npt
logical isIn

integer n,num_pt
character (len=SEIS_STRLEN) :: filenm
integer,allocatable :: ids(:,:)

isIn=.false.
filenm= get_fnm_station(pnm_sta,n_i,n_j,n_k)
call nfseis_diminfo(trim(filenm),'num_pt',num_pt)
if (num_pt>0) then
   allocate(ids(2,num_pt))
   call nfseis_varget(trim(filenm),'id', ids, &
        (/1,1/),(/2,num_pt/),(/1,1/))
do n=1,num_pt
if (ids(2,n)==id .and. ids(1,n)==indx) then
   isIn=.true.
   npt=n
   exit
end if
end do
   deallocate(ids)
end if
end function seismo_on_this
subroutine retrieve_recvline(fnm_nc,id,vnm,V,nts,ntc,ntt,T)
character (len=*),intent(in) :: fnm_nc,vnm
integer,intent(in) :: id,nts,ntc,ntt
real,dimension(:),intent(out) :: V
real,dimension(:),optional,intent(out) :: T

integer numt

call nfseis_diminfo(fnm_nc,'time',numt)
if (numt<1) then
   print *, trim(fnm_nc), ' no data'
   stop 1
end if
if (numt-nts+1<ntc/ntt) then
   print *, trim(fnm_nc), ' time.length-nts+1<ntc/ntt'
   print *, numt,nts,ntc,ntt
   stop 1
end if

call nfseis_varget(fnm_nc,trim(vnm),V,(/id,nts/),(/1,ntc/),(/1,ntt/))
if (present(T)) then
call nfseis_varget(fnm_nc,'time',T,(/nts/),(/ntc/),(/ntt/))
end if
end subroutine retrieve_recvline

subroutine retrieve_snap_seis(fnm_prefix,i,j,k,varnm,var,num_nt)
character (len=*),intent(in) :: fnm_prefix,varnm
integer,intent(in) :: i,j,k,num_nt
real,dimension(:),intent(out) :: var
character (len=SEIS_STRLEN) :: filenm

integer m,n,mt
m=0; n=0
do 
   if (m>=num_nt) exit
   n=n+1
   filenm=trim(fnm_prefix)               &
       //'_n'//trim(io_out_pattern(n,5)) &
       //'.nc'
   call nfseis_diminfo(filenm,'time',mt)
   if (m+mt>num_nt) mt=num_nt-m
   call nfseis_varget(filenm,varnm,var(m+1:m+mt), &
        (/i,j,k,1/),(/1,1,1,mt/),(/1,1,1,1/)  &
        )
   m=m+mt
end do
end subroutine retrieve_snap_seis
subroutine retrieve_snap_time(fnm_prefix,T,num_nt)
character (len=*),intent(in) :: fnm_prefix
real,dimension(:),intent(out) :: T
integer,intent(in) :: num_nt
character (len=SEIS_STRLEN) :: filenm

integer m,n,mt
m=0; n=0
do 
   if (m>=num_nt) exit
   n=n+1
   filenm=trim(fnm_prefix)               &
       //'_n'//trim(io_out_pattern(n,5)) &
       //'.nc'
   call nfseis_diminfo(filenm,'time',mt)
   if (m+mt>num_nt) mt=num_nt-m
   call nfseis_varget(filenm,'time',T(m+1:m+mt), &
        (/1/),(/mt/),(/1/))
   m=m+mt
end do
end subroutine retrieve_snap_time

!---------------------------------------------------------------------------
function get_fnm_station(pnm_sta,n_i,n_j,n_k) result(filenm)
  character (len=*),intent(in) :: pnm_sta
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_sta)                    &
       //'/'//'station'                   &
       //'_'//set_mpi_subfix(n_i,n_j,n_k) &
       //'.nc'
end function get_fnm_station
function get_fnm_seismo(pnm,n_i,n_j,n_k) result(filenm)
  character (len=*),intent(in) :: pnm
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm    )                          &
       //'/'//'seismo'                          &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'.nc'
end function get_fnm_seismo
function get_fnm_snapnode(pnm,prefix,n,ntime,n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n,ntime,n_i,n_j,n_k
  character (len=*),intent(in) :: pnm,prefix
  character (len=SEIS_STRLEN) :: filenm
  integer n0
  n0=(ntime+snap_tinv(n)*snap_tcnt(n)-1)/(snap_tinv(n)*snap_tcnt(n))
  filenm=trim(pnm)                              &
       //'/'//trim(io_enum(prefix,n))           &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'_n'//trim(io_out_pattern(n0,5))       &
       //'.nc'
end function get_fnm_snapnode
function get_fnm_snapnode_n(pnm,prefix,id,n0,n_i,n_j,n_k) result(filenm)
  character (len=*),intent(in) :: pnm,prefix
  integer,intent(in) :: id,n0,n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm)                              &
       //'/'//trim(io_enum(prefix,id))          &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'_n'//trim(io_out_pattern(n0,5))       &
       //'.nc'
end function get_fnm_snapnode_n
function get_fnm_rest(pnm,ntime,n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: ntime,n_i,n_j,n_k
  character (len=*),intent(in) :: pnm
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm)//'/'                          &
       //'rest_t'//trim(io_out_pattern(ntime,5)) &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))  &
       //'.nc'
end function get_fnm_rest

!function set_mpi_subfix(i,j,k) result(filenm)
!    integer i,j,k
!    character (len=SEIS_STRLEN) :: filenm
!    character (len=SEIS_STRLEN) :: str1,str2,str3
!    write(str1,"(i2.2)") i
!    write(str2,"(i2.2)") j
!    write(str3,"(i2.2)") k
!    filenm  ='mpi'//trim(adjustl(str1))  &
!                    //trim(adjustl(str2))  &
!                    //trim(adjustl(str3))
!end function set_mpi_subfix

function io_enum(prefix,num) result(ioname)
character (len=*),intent(in) :: prefix
integer,intent(in) :: num
character (len=SEIS_STRLEN) :: ioname
character (len=SEIS_STRLEN) :: str
write(str,"(i3.3)") num
ioname=trim(prefix)//trim(str)
end function io_enum

function io_out_pattern(num,width) result(ioname)
integer,intent(in) :: num
integer,optional,intent(in) :: width
character (len=SEIS_STRLEN) :: ioname,str,fmt_str
if (present(width)) then
   write(str,"(i)") width
   fmt_str="(i"//trim(str)//"."//trim(str)//")"
else
   fmt_str="(i4.4)"
end if
write(ioname,fmt_str) num
end function io_out_pattern

subroutine io_delete(filenm)
character (len=*),intent(in) :: filenm
integer fid
fid=5001
open(fid,file=trim(filenm),status='unknown')
close(fid,status='delete')
end subroutine io_delete

end module io_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
