program seis3d_metric

! This program calculates the metric coefficient
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2009 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:53:22 -0500 (Fri, 16 Jan 2009) $
! $Revision: 67 $
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

real,allocatable :: gx(:),gy(:),gz(:),gxsin(:),gxcot(:), &
     x_xi(:),y_eta(:),z_zeta(:)
character (len=SEIS_STRLEN) :: filenm

integer :: ierr
integer :: n,i,j,k,num_x,num_y,num_z,num_i,num_j,num_k
integer :: n_i,n_j,n_k

integer :: vecF(5),vecB(5)
real(SP) :: coeF(5),coeB(5)

vecF=(/ -1,0,1,2,3 /)
coeF=(/ -0.30874,-0.6326,1.2330,-0.3334,0.04168 /)
vecB=(/ -3,-2,-1,0,1 /)
coeB=(/ -0.04168,0.3334,-1.2330,0.6326,0.30874 /)

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
call grid_alloc

! read to global coordinate
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_coord_import(n_i,n_j,n_k)
   gx(ngx1:ngx2)=x(nx1:nx2)
   gy(ngy1:ngy2)=y(ny1:ny2)
   gz(ngz1:ngz2)=z(nz1:nz2)
end do
end do
end do

! calculate metric
do i=nx1,num_x
   gxsin(i)=sin(real(gx(i),DP))
   gxcot(i)=1.0_DP/tan(real(gx(i),DP))
end do
do i=ni1,num_i
   x_xi(i) = (                        &
      dot_product(gx(i+vecB),coeB) &
     +dot_product(gx(i+vecF),coeF) &
     )*0.5
   xi_x(i)=1.0_SP/x_xi(i)
end do
do j=nj1,num_j
   y_eta(j) = (                       &
      dot_product(gy(j+vecB),coeB) &
     +dot_product(gy(j+vecF),coeF) &
     )*0.5
   eta_y(j)=1.0_SP/y_eta(j)
end do
do k=nk1,num_k
   z_zeta(k) = (                      &
      dot_product(gz(k+vecB),coeB) &
     +dot_product(gz(k+vecF),coeF) &
     )*0.5
   zeta_z(k)=1.0_SP/z_zeta(k)
end do
do n=1,LenFD
   xi_x(ni1-n)=xi_x(ni1); xi_x(num_i+n)=xi_x(num_i)
   eta_y(nj1-n)=eta_y(nj1); eta_y(num_j+n)=eta_y(num_j)
   zeta_z(nk1-n)=zeta_z(nk1); zeta_z(num_k+n)=zeta_z(num_k)
end do

print *, "output metric ..."
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   !n_k=dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   !output grid
   filenm=grid_metricfnm_get(n_i,n_j,n_k)
   call metric_skel(filenm,nx,ny,nz)
   call nfseis_varput(filenm,'xsin',gxsin(ngx1:ngx2),(/1/),(/nx/),(/1/))
   call nfseis_varput(filenm,'xcot',gxcot(ngx1:ngx2),(/1/),(/nx/),(/1/))
   call nfseis_varput(filenm,'xi_x',xi_x(ngx1:ngx2),(/1/),(/nx/),(/1/))
   call nfseis_varput(filenm,'eta_y',eta_y(ngy1:ngy2),(/1/),(/ny/),(/1/))
   call nfseis_varput(filenm,'zeta_z',zeta_z(ngz1:ngz2),(/1/),(/nz/),(/1/))

end do
end do
end do

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
subroutine alloc_local(n1,n2,n3)
integer,intent(in) :: n1,n2,n3
allocate(gx(nx1:n1)); gx=0.0
allocate(gxsin(nx1:n1)); gxsin=0.0
allocate(gxcot(nx1:n1)); gxcot=0.0
allocate(gy(ny1:n2)); gy=0.0
allocate(gz(nz1:n3)); gz=0.0
allocate(x_xi(nx1:n1)); x_xi=0.0
allocate(y_eta(ny1:n2)); y_eta=0.0
allocate(z_zeta(nz1:n3)); z_zeta=0.0
allocate(xi_x(nx1:n1)); xi_x=0.0
allocate(eta_y(ny1:n2)); eta_y=0.0
allocate(zeta_z(nz1:n3)); zeta_z=0.0
end subroutine

subroutine metric_skel(filenm,nx,ny,nz)
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
ierr=nf90_def_var(ncid, 'xsin', nf90_float, (/ xid /), vid )
ierr=nf90_def_var(ncid, 'xcot', nf90_float, (/ xid /), vid )
ierr=nf90_def_var(ncid, 'xi_x', nf90_float, (/ xid /), vid )
ierr=nf90_def_var(ncid, 'eta_y', nf90_float, (/ yid /), vid )
ierr=nf90_def_var(ncid, 'zeta_z', nf90_float, (/ zid /), vid )
!--
ierr=nf90_enddef(ncid)
!--
ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine metric_skel

end program seis3d_metric

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
