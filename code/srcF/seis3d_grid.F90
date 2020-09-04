program seis3d_grid

! This program generates grid coordinate
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:49:59 -0500 (Fri, 16 Jan 2009) $
! $Revision: 65 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod
use string_mod
use para_mod
use math_mod
use mpi_mod
use nfseis_mod
use grid_mod
use io_mod

implicit none

!-----------------------------------------------------------------------------

real,allocatable :: gx(:),gy(:),gz(:)
character (len=SEIS_STRLEN) :: filenm

integer :: ierr
integer :: n,i,j,k,num_x,num_y,num_z,num_i,num_j,num_k
integer :: n_i,n_j,n_k

!-----------------------------------------------------------------------------

call get_conf_name(fnm_conf)

call swmpi_init(fnm_conf)
call para_init(fnm_conf)
call grid_fnm_init(fnm_conf)

!-----------------------------------------------------------------------------

! calculate coordinate
num_x=swmpi_globi(nx2,dims(1)-1)
num_y=swmpi_globj(ny2,dims(2)-1)
num_z=swmpi_globk(nz2,dims(3)-1)
num_i=swmpi_globi(ni2,dims(1)-1)
num_j=swmpi_globj(nj2,dims(2)-1)
num_k=swmpi_globk(nk2,dims(3)-1)

call alloc_local(num_x,num_y,num_z)

call read_grid_coord(fnm_grid_conf) 

print *, "output grid coordinate and locating ..."
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   !n_k=dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   !output grid
   filenm=grid_coordfnm_get(n_i,n_j,n_k)
   call coord_skel(filenm,nx,ny,nz)
   call nfseis_varput(filenm,'x',gx(ngx1:ngx2),(/1/),(/nx/),(/1/))
   call nfseis_varput(filenm,'y',gy(ngy1:ngy2),(/1/),(/ny/),(/1/))
   call nfseis_varput(filenm,'z',gz(ngz1:ngz2),(/1/),(/nz/),(/1/))
end do
end do
end do

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
subroutine alloc_local(n1,n2,n3)
integer,intent(in) :: n1,n2,n3
allocate(gx(nx1:n1)); gx=0.0
allocate(gy(ny1:n2)); gy=0.0
allocate(gz(nz1:n3)); gz=0.0
end subroutine

subroutine read_grid_coord(fnm_conf)
character (len=*),intent(in) :: fnm_conf

character (len=SEIS_STRLEN) :: str
real(SP),dimension(SEIS_GEO) :: dh,xyz0
integer :: fid,i,j,k,n,n1,n2,n3
real :: d2m
fid=1001
open(fid,file=trim(fnm_conf),status="old")

  call string_conf(fid,1,'distance2meter',2,d2m)
  do n=1,SEIS_GEO
     call string_conf(fid,1,'steph',n+1,dh(n))
     call string_conf(fid,1,'theta0_phi0_rmax',n+1,xyz0(n))
  end do
  !x
  if (abs(dh(1)>=SEIS_ZERO)) then
     do i=nx1,num_x
        gx(i)=(i-ni1)*dh(1)+xyz0(1)
     end do
  else
     n1=ni*dims(1)
     call string_conf(fid,1,'<x',2,str)
     read(fid,*) ( gx(i),i=ni1,n1+ni1-1 )
  end if
  !y
  if (abs(dh(2)>=SEIS_ZERO)) then
     do j=ny1,num_y
        gy(j)=(j-nj1)*dh(2)+xyz0(2)
     end do
  else
     n2=nj*dims(2)
     call string_conf(fid,1,'<y',2,str)
     read(fid,*) ( gy(j),j=nj1,n2+nj1-1 )
  end if
  !z
  if (abs(dh(3)>=SEIS_ZERO)) then
     do k=nz1,num_z
        gz(k)=(k-num_z+LenFD)*dh(3)+xyz0(3)
     end do
  else
     n3=nk*dims(3)
     call string_conf(fid,1,'<z',2,str)
     read(fid,*) ( gz(k),k=nk1,n3+nk1-1 )
  end if

  gx=gx/180.0_DP*PI; gy=gy/180.0_DP*PI; gz=gz*d2m

  n1=swmpi_globi(ni2,dims(1)-1)
  n2=swmpi_globj(nj2,dims(2)-1)
  n3=swmpi_globk(nk2,dims(3)-1)
  do n=1,LenFD
     gx(ni1-n)=gx(ni1-n+1)+gx(ni1)-gx(ni1+1)
     gx(n1+n)=gx(n1+n-1)+gx(n1)-gx(n1-1)
     gy(nj1-n)=gy(nj1-n+1)+gy(nj1)-gy(nj1+1)
     gy(n2+n)=gy(n2+n-1)+gy(n2)-gy(n2-1)
     gz(nk1-n)=gz(nk1-n+1)+gz(nk1)-gz(nk1+1)
     gz(n3+n)=gz(n3+n-1)+gz(n3)-gz(n3-1)
  end do
close(fid)
end subroutine read_grid_coord

subroutine coord_skel(filenm,nx,ny,nz)
character (len=*),intent(in) :: filenm
integer,intent(in) :: nx,ny,nz

integer ncid,ierr,oldMode
integer xid,yid,zid
integer vid

ierr=nf90_create( path=trim(filenm),cmode=nf90_clobber,ncid=ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_handle_err(ierr)

! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, xid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'J', ny, yid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'K', nz, zid);     call nfseis_handle_err(ierr)

! -- variable
ierr=nf90_def_var(ncid, 'x', nf90_float, (/ xid /), vid )
ierr=nf90_def_var(ncid, 'y', nf90_float, (/ yid /), vid )
ierr=nf90_def_var(ncid, 'z', nf90_float, (/ zid /), vid )
!--
ierr=nf90_enddef(ncid)
!--
ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine coord_skel

end program seis3d_grid

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
