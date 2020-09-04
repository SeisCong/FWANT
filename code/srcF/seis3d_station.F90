program seis3d_station

! This program set distribute station among mpi threads
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2007 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:56:20 -0500 (Fri, 16 Jan 2009) $
! $Revision: 69 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod
use string_mod, only : string_conf
use para_mod
use math_mod
use mpi_mod
use nfseis_mod
use grid_mod
use io_mod

implicit none

integer :: n_i,n_j,n_k

!-----------------------------------------------------------------------------

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)

call grid_fnm_init(fnm_conf)
call grid_alloc

call io_init(fnm_conf)
call io_pt_read(fnm_conf)

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1

  write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
  call swmpi_change_fnm(n_i,n_j,n_k)
  call swmpi_set_gindx(n_i,n_j,n_k)
  call grid_coord_import(n_i,n_j,n_k)

  call receiver_locate(n_i,n_j,n_k)

end do
end do
end do

call grid_dealloc

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine receiver_locate(n_i,n_j,n_k)
integer,intent(in) :: n_i,n_j,n_k
character (len=SEIS_STRLEN) :: filenm
integer,dimension(1) :: p
integer npt,n,i,j,k,gi,gj,gk
real(SP) x0,y0,z0
npt=0
filenm=get_fnm_station(pnm_station,n_i,n_j,n_k)
call station_skel(filenm,pt_tinv)
do n=1,num_pt
   x0=pt_xyz(1,n);y0=pt_xyz(2,n);z0=pt_xyz(3,n)
   if (z0<=topo_hyper_height) then
      p=minloc(abs(x-x0)); i=p(1)
      p=minloc(abs(y-y0)); j=p(1)
      p=minloc(abs(z-z0)); k=p(1)
   elseif (n_k==dims(3)-1) then
      p=minloc(abs(x-x0)); i=p(1)
      p=minloc(abs(y-y0)); j=p(1)
      k=nk2
   else
      i=nx1;j=ny1;k=nz1
   end if
   if ( i>=ni1 .and. i<=ni2      &
       .and. j>=nj1 .and. j<=nj2 &
       .and. k>=nk1 .and. k<=nk2 ) then
      npt=npt+1
      gi=swmpi_globi(i,n_i)
      gj=swmpi_globj(j,n_j)
      gk=swmpi_globk(k,n_k)
      call nfseis_varput(filenm,'indx',(/i,j,k/),                          &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'gindx',(/out_i(gi),out_j(gj),out_k(gk)/), &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'coord',(/x0,y0,z0/),                      &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'grid',(/x(i),y(j),z(k)/),                 &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'id',pt_id(:,n),                           &
           (/1,npt/),(/2,1/),(/1,1/))
   end if
end do
end subroutine receiver_locate

subroutine station_skel( filenm,tinv,title )
character (len=*),intent(in) :: filenm
integer,intent(in) :: tinv
character (len=*),optional,intent(in) :: title
integer :: ierr,ncid,gdimid,ndimid,oldMode
integer :: indxid,gindxid,coordid,loctid,idid,twoelemid
!--
ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'station_skel:'//trim(filenm))
!--
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in station_skel')
!-- define dim
ierr=nf90_def_dim(ncid, "num_pt", nf90_unlimited, ndimid)
     call nfseis_except(ierr,'num_pt dim in station_skel')
ierr=nf90_def_dim(ncid, "geo_dim", SEIS_GEO, gdimid)
     call nfseis_except(ierr,'geo_dim dim in station_skel')
ierr=nf90_def_dim(ncid, "twoelem", 2, twoelemid)
     call nfseis_except(ierr,'twoelem dim in station_skel')
! -- define variable
ierr=nf90_def_var(ncid, 'indx', nf90_int, (/ gdimid, ndimid /), indxid )
     call nfseis_except(ierr,'indx var in station_skel')
ierr=nf90_def_var(ncid, 'gindx', nf90_int, (/ gdimid, ndimid /), gindxid )
     call nfseis_except(ierr,'gindx var in station_skel')
ierr=nf90_def_var(ncid, 'coord', SEISNC_DATATYPE, (/ gdimid,ndimid /), coordid )
     call nfseis_except(ierr,'coord var in station_skel')
ierr=nf90_def_var(ncid, 'grid', SEISNC_DATATYPE, (/ gdimid,ndimid /), loctid )
     call nfseis_except(ierr,'grid var in station_skel')
ierr=nf90_def_var(ncid, 'id', nf90_int, (/ twoelemid,ndimid /), idid )
     call nfseis_except(ierr,'id var in station_skel')
!-- define global attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"tinv",tinv )
     call nfseis_except(ierr,'tinv att in station_skel')
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in station_skel')
end if
!--
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in station_skel')
!-- 
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in station_skel')
end subroutine station_skel

end program seis3d_station

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
