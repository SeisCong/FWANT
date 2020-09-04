program tool_expt_snap

! This program gathers one snap data from output node nc file
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 13:01:39 -0500 (Fri, 16 Jan 2009) $
! $Revision: 71 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

use constants_mod
use string_mod, only : string_conf
use math_mod
use para_mod
use mpi_mod
use nfseis_mod
use grid_mod
use media_mod
use io_mod

implicit none

type SNAPBD
     integer n_i,n_j,n_k
     integer,dimension(SEIS_GEO) :: indxs,indxe,indxc
     integer,dimension(SEIS_GEO) :: subs,sube,subt,subc
end type SNAPBD

integer nthd
type(SNAPBD),allocatable,dimension(:) :: info_list
integer id,nlayer
character (len=SEIS_STRLEN) :: filenm,pnm_expt,fnm_expt
real,dimension(:,:,:),allocatable :: Vx,Vy,Vz
integer n,i,j,k,si,sj,sk,n_i,n_j,n_k,i1,i2,j1,j2,k1,k2,nofnc,ninnc
integer,dimension(SEIS_GEO) :: subs,sube,subt,subc

!---------------------------------------------

pnm_expt='./export/'

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)

write(*,*) "Input snap id and time layer number(0 means coordinate)"
read(*,*) id, nlayer

allocate(info_list(dims(1)*dims(2)*dims(3)))
allocate(Vx(snap_subc(1,id),snap_subc(2,id),snap_subc(3,id)))
allocate(Vy(snap_subc(1,id),snap_subc(2,id),snap_subc(3,id)))
allocate(Vz(snap_subc(1,id),snap_subc(2,id),snap_subc(3,id)))

call locate_snap(id,nthd)

if (nlayer==0) then 

!-- get coord data
do n=1,nthd
   n_i=info_list(n)%n_i;n_j=info_list(n)%n_j;n_k=info_list(n)%n_k
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   i1=info_list(n)%indxs(1);j1=info_list(n)%indxs(2);k1=info_list(n)%indxs(3);
   i2=info_list(n)%indxe(1);j2=info_list(n)%indxe(2);k2=info_list(n)%indxe(3);
   subs=info_list(n)%subs;subc=info_list(n)%subc;subt=info_list(n)%subt;
   sube=info_list(n)%sube
   filenm=grid_coordfnm_get(n_i,n_j,n_k)

   call nfseis_varget(filenm,'x',Vx(i1:i2,j1:j2,k1:k2),subs,subc,subt) 
   call nfseis_varget(filenm,'y',Vy(i1:i2,j1:j2,k1:k2),subs,subc,subt) 
   call nfseis_varget(filenm,'z',Vz(i1:i2,j1:j2,k1:k2),subs,subc,subt) 
end do
filenm=trim(pnm_expt)//'/'//trim(io_enum('snap_',id))//'_coord.grd'
call export_grd(filenm,Vx,Vy,Vz)

else

nofnc=(nlayer-1)/snap_tcnt(id)+1
ninnc=mod( (nlayer-1),snap_tcnt(id) )+1
do n=1,nthd
   n_i=info_list(n)%n_i;n_j=info_list(n)%n_j;n_k=info_list(n)%n_k
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   i1=info_list(n)%indxs(1);j1=info_list(n)%indxs(2);k1=info_list(n)%indxs(3);
   i2=info_list(n)%indxe(1);j2=info_list(n)%indxe(2);k2=info_list(n)%indxe(3);
   subc=info_list(n)%indxc
   filenm=get_fnm_snapnode_n(pnm_out,'vel_',id,nofnc,n_i,n_j,n_k)
   call nfseis_varget(filenm,'Vx',Vx(i1:i2,j1:j2,k1:k2), &
        (/1,1,1,ninnc/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_varget(filenm,'Vy',Vy(i1:i2,j1:j2,k1:k2), &
        (/1,1,1,ninnc/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_varget(filenm,'Vz',Vz(i1:i2,j1:j2,k1:k2), &
        (/1,1,1,ninnc/),(/subc,1/),(/1,1,1,1/) )
end do
filenm=trim(pnm_expt)//'/'//trim(io_enum('vel_',id)) &
       //'_ndim'//io_out_pattern(nlayer,5)//'.grd'
call export_grd(filenm,Vx,Vy,Vz)

end if

deallocate(info_list)
deallocate(Vx,Vy,Vz)

print *, 'export snap finished'

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine locate_snap(id,nthd)
integer id,nthd
integer n_i1,n_i2,n_j1,n_j2,n_k1,n_k2
integer n_i,n_j,n_k
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube

nthd=0;
n_i1=(snap_subs(1,id)-ni1)/ni
n_i2=(snap_sube(1,id)-ni1)/ni
n_j1=(snap_subs(2,id)-ni1)/nj
n_j2=(snap_sube(2,id)-ni1)/nj
n_k1=(snap_subs(3,id)-ni1)/nk
n_k2=(snap_sube(3,id)-ni1)/nk

do n_i=n_i1,n_i2
do n_j=n_j1,n_j2
do n_k=n_k1,n_k2
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   nthd=nthd+1;
   subs=snap_subs(:,id); sube=snap_sube(:,id)
   subc=snap_subc(:,id); subt=snap_subt(:,id)
   call corr_subse(subs,subc,subt,sube)
   info_list(nthd)%n_i=n_i
   info_list(nthd)%n_j=n_j
   info_list(nthd)%n_k=n_k
   info_list(nthd)%indxs=(subs-snap_subs(:,id))/subt+1
   info_list(nthd)%indxe=(sube-snap_subs(:,id))/subt+1
   info_list(nthd)%indxc=subc
   info_list(nthd)%subs(1)=swmpi_locli(subs(1),n_i)
   info_list(nthd)%subs(2)=swmpi_loclj(subs(2),n_j)
   info_list(nthd)%subs(3)=swmpi_loclk(subs(3),n_k)
   info_list(nthd)%sube(1)=swmpi_locli(sube(1),n_i)
   info_list(nthd)%sube(2)=swmpi_loclj(sube(2),n_j)
   info_list(nthd)%sube(3)=swmpi_loclk(sube(3),n_k)
   info_list(nthd)%subc=subc
   info_list(nthd)%subt=subt
end do
end do
end do
end subroutine locate_snap

subroutine export_grd(filenm,x,y,z)
character (len=*) :: filenm
real,dimension(:,:,:) :: x,y,z
integer fid,i,j,k,ni,nj,nk

ni=size(x,1); nj=size(x,2); nk=size(x,3)
fid=1001
open(fid,file=trim(filenm),status='unknown')
write(fid,*) 1
write(fid,*) ni,nj,nk
write(fid,*) ( ((x(i,j,k), i=1,ni ), j=1,nj),k=1,nk ),  &
             ( ((y(i,j,k), i=1,ni ), j=1,nj),k=1,nk ),  &
             ( ((z(i,j,k), i=1,ni ), j=1,nj),k=1,nk )
close(fid)
end subroutine export_grd

end program tool_expt_snap

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
