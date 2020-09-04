module nfseis_mod

! This module relates with netcdf subroutines
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:45:22 -0500 (Fri, 16 Jan 2009) $
! $Revision: 62 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use netcdf
use constants_mod
implicit none

interface nfseis_attget
   module procedure nfseis_attget_int
   module procedure nfseis_attget_int1d
   module procedure nfseis_attget_real
   module procedure nfseis_attget_real1d
   module procedure nfseis_attget_character
end interface

interface nfseis_attput
   module procedure nfseis_attput_int
   module procedure nfseis_attput_int1d
   module procedure nfseis_attput_real
   module procedure nfseis_attput_real1d
end interface

interface nfseis_varget
   module procedure nfseis_varget_int
   module procedure nfseis_varget_int1d
   module procedure nfseis_varget_int2d
   module procedure nfseis_varget_int3d
   module procedure nfseis_varget_int4d
   module procedure nfseis_varget_real
   module procedure nfseis_varget_real1d
   module procedure nfseis_varget_real2d
   module procedure nfseis_varget_real3d
   module procedure nfseis_varget_real4d
end interface

interface nfseis_varput
   module procedure nfseis_varput_int
   module procedure nfseis_varput_int1d
   module procedure nfseis_varput_int2d
   module procedure nfseis_varput_int3d
   module procedure nfseis_varput_real
   module procedure nfseis_varput_real1d
   module procedure nfseis_varput_real2d
   module procedure nfseis_varput_real3d
   module procedure nfseis_varput_real4d
end interface

interface nfseis_put
   module procedure nfseis_put_int
   module procedure nfseis_put_int1d
   module procedure nfseis_put_int2d
   module procedure nfseis_put_int3d
   module procedure nfseis_put_real
   module procedure nfseis_put_real1d
   module procedure nfseis_put_real2d
   module procedure nfseis_put_real3d
   module procedure nfseis_put_real4d
end interface

!-------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------!

!*************************************************************************
!*                    seismogram data file                               *
!*************************************************************************
subroutine nfseis_seismo_def(filenm,num_pt,ncid,tid,title)
character (len=*),intent(in) :: filenm
integer,intent(in) :: num_pt
integer,intent(out) :: ncid,tid
character (len=*),intent(in),optional :: title
integer :: ierr,ndimid,tdimid,oldMode
!--
ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'seismo_def:'//trim(filenm))
!--
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in seismo_def')
!-- define dim
ierr=nf90_def_dim(ncid, "num_pt", num_pt, ndimid)
     call nfseis_except(ierr,'num_pt dim in seismo_def')
ierr=nf90_def_dim(ncid, "time", nf90_unlimited, tdimid)
     call nfseis_except(ierr,'time dim in seismo_def')
! -- define variable
ierr=nf90_def_var(ncid, 'time', SEISNC_DATATYPE, (/ tdimid /), tid )
     call nfseis_except(ierr,'time var in seismo_def')
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in seismo_def')
end if
end subroutine nfseis_seismo_def
subroutine nfseis_seismo_defvar(ncid,varnm,vid)
integer,intent(in) :: ncid
character (len=*),intent(in) :: varnm
integer,intent(out) :: vid
integer :: ierr,ndimid,tdimid
!--
ierr=nf90_inq_dimid(ncid,'time',tdimid)
     call nfseis_except(ierr,'time dim in seismo_addvar')
ierr=nf90_inq_dimid(ncid,'num_pt',ndimid)
     call nfseis_except(ierr,'num_pt dim in seismo_addvar')
ierr=nf90_def_var(ncid,varnm,SEISNC_DATATYPE, &
     (/ ndimid,tdimid /), vid )
     call nfseis_except(ierr,'add var in sesimo_addvar')
end subroutine nfseis_seismo_defvar
subroutine nfseis_seismo_enddef(ncid)
integer ncid
integer ierr
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in seismo_enddef')
end subroutine nfseis_seismo_enddef

!*************************************************************************
!*                     3D snap wave field file                           *
!*************************************************************************
subroutine nfseis_snap_def(filenm,ncid,tid,stept,subc,title)
character (len=*),intent(in) :: filenm
integer,intent(out) :: ncid,tid
real(SP),intent(in) :: stept
integer,dimension(SEIS_GEO),intent(in) :: subc
character (len=*),intent(in),optional :: title

integer ierr,oldMode
integer tdimid,xdimid,ydimid,zdimid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'wave_init:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in wave_init')
! -- define dim
ierr=nf90_def_dim(ncid, 'I', subc(1), xdimid)
     call nfseis_except(ierr,'I dim in wave_init')
ierr=nf90_def_dim(ncid, 'J', subc(2), ydimid)
     call nfseis_except(ierr,'J dim in wave_init')
ierr=nf90_def_dim(ncid, 'K', subc(3), zdimid)
     call nfseis_except(ierr,'K dim in wave_init')
ierr=nf90_def_dim(ncid, 'time', nf90_unlimited, tdimid)
     call nfseis_except(ierr,'time dim in wave_init')
ierr=nf90_def_var(ncid,'time',SEISNC_DATATYPE,(/ tdimid /), tid )
     call nfseis_except(ierr,'time var in wave_init')
ierr=nf90_put_att(ncid,NF90_GLOBAL,"stept",stept )
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in wave_init')
end if
!--
!ierr=nf90_enddef(ncid)
!     call nfseis_except(ierr,'enddef in wave_init')
end subroutine nfseis_snap_def
subroutine nfseis_snap_defvar(ncid,varnm,vid)
character (len=*),intent(in) :: varnm
integer,intent(out) :: ncid,vid
integer ierr,tdimid,xdimid,ydimid,zdimid
ierr=nf90_inq_dimid(ncid,'time',tdimid)
     call nfseis_except(ierr,'time dim in snap_addvar')
ierr=nf90_inq_dimid(ncid,'I',xdimid)
     call nfseis_except(ierr,'I dim in snap_addvar')
ierr=nf90_inq_dimid(ncid,'J',ydimid)
     call nfseis_except(ierr,'J dim in snap_addvar')
ierr=nf90_inq_dimid(ncid,'K',zdimid)
     call nfseis_except(ierr,'K dim in snap_addvar')
ierr=nf90_def_var(ncid,varnm,SEISNC_DATATYPE, &
     (/ xdimid, ydimid, zdimid,tdimid /), vid )
     call nfseis_except(ierr,'add var in snap_addvar')
end subroutine nfseis_snap_defvar
subroutine nfseis_snap_attdef(ncid, &
  gsubs,gsubc,gsubt,gsube,subs,subc,subt,sube)
integer,intent(in) :: ncid
integer,dimension(SEIS_GEO),intent(in) :: gsubs,gsubc,gsubt,gsube
integer,dimension(SEIS_GEO),intent(in) :: subs,subc,subt,sube
integer ierr
!-- define variable attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsubs",gsubs)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsubc",gsubc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsubt",gsubt)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsube",gsube)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"subs",subs)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"subc",subc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"subt",subt)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"sube",sube)
end subroutine nfseis_snap_attdef
subroutine nfseis_snap_enddef(ncid)
integer ncid
integer ierr
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in snap_enddef')
end subroutine nfseis_snap_enddef

!*************************************************************************
!*                    3D grid point value file                           *
!*************************************************************************
subroutine nfseis_grid3d_def(filenm,nx,ny,nz,ncid,title)
character (len=*),intent(in) :: filenm
integer,intent(in) :: nx,ny,nz
character (len=*),intent(in),optional :: title
integer :: ncid

integer :: ierr,oldMode
integer :: Iid,Jid,Kid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'grid3d_def:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in grid3d_def')
! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, Iid)
     call nfseis_except(ierr,'I dim in grid3d_def')
ierr=nf90_def_dim(ncid, 'J', ny, Jid)
     call nfseis_except(ierr,'J dim in grid3d_def')
ierr=nf90_def_dim(ncid, 'K', nz, Kid)
     call nfseis_except(ierr,'K dim in grid3d_def')
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in grid3d_def')
end if
end subroutine nfseis_grid3d_def
subroutine nfseis_grid3d_defvar(ncid,varnm,vid)
integer,intent(in) :: ncid
character (len=*),intent(in) :: varnm
integer,intent(out) :: vid
integer :: ierr,Iid,Jid,Kid
!--
ierr=nf90_inq_dimid(ncid,'I',Iid)
     call nfseis_except(ierr,'I dim in grid3d_addvar')
ierr=nf90_inq_dimid(ncid,'J',Jid)
     call nfseis_except(ierr,'J dim in grid3d_addvar')
ierr=nf90_inq_dimid(ncid,'K',Kid)
     call nfseis_except(ierr,'K dim in grid3d_addvar')
ierr=nf90_def_var(ncid,varnm,SEISNC_DATATYPE, &
     (/ Iid,Jid,Kid /), vid )
     call nfseis_except(ierr,'add var in grid3d_addvar')
end subroutine nfseis_grid3d_defvar
subroutine nfseis_grid3d_enddef(ncid)
integer ncid
integer ierr
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in grid3d_enddef')
end subroutine nfseis_grid3d_enddef

subroutine nfseis_grid3d_skel(filenm,nx,ny,nz,title)
character (len=*),intent(in) :: filenm
integer,intent(in) :: nx,ny,nz
character (len=*),intent(in),optional :: title

integer :: ncid,ierr,oldMode
integer :: Iid,Jid,Kid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'grid3d_header:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in grid3d_header')
! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, Iid)
     call nfseis_except(ierr,'I dim in grid3d_header')
ierr=nf90_def_dim(ncid, 'J', ny, Jid)
     call nfseis_except(ierr,'J dim in grid3d_header')
ierr=nf90_def_dim(ncid, 'K', nz, Kid)
     call nfseis_except(ierr,'K dim in grid3d_header')
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in grid3d_header')
end if
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'endedf in grid3d_header')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in grid3d_header')
end subroutine nfseis_grid3d_skel

subroutine nfseis_grid3d_addvar(filenm,aVarNm)
character (len=*),intent(in) :: aVarNm,filenm

integer ncid,ierr,oldMode
integer Iid,Jid,Kid
integer varid

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_except(ierr,'grid3d_addvar:'//trim(filenm))
ierr=nf90_redef( ncid)
     call nfseis_except(ierr,'redef in grid3d_addvar')
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in grid3d_addvar')
! -- req dim
ierr=nf90_inq_dimid(ncid,'I',Iid)
     call nfseis_except(ierr,'I dim in grid3d_addvar')
ierr=nf90_inq_dimid(ncid,'J',Jid)
     call nfseis_except(ierr,'J dim in grid3d_addvar')
ierr=nf90_inq_dimid(ncid,'K',Kid)
     call nfseis_except(ierr,'K dim in grid3d_addvar')
! -- define variable
ierr=nf90_def_var(ncid, trim(aVarNm), SEISNC_DATATYPE,(/ Iid,Jid,Kid /),varid)
     call nfseis_except(ierr,'var def in grid3d_addvar')
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in grid3d_addvar')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in grid3d_addvar')
end subroutine nfseis_grid3d_addvar

subroutine nfseis_grid3d_attput(filenm,         &
           subs,subc,subt,                    &
           ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,  &
           nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz,  &
           ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,     &
           ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,     &
           point_in_this)
character (len=*),intent(in) :: filenm
integer,dimension(SEIS_GEO),intent(in) :: subs,subc,subt
integer,intent(in) ::                         &
        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,     &
        nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz,     &
        ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,        &
        ngx1,ngx2,ngy1,ngy2,ngz1,ngz2 !,      &
integer,dimension(SEIS_GEO*2) :: point_in_this

integer ierr,ncid

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
ierr=nf90_redef( ncid)
!-- define variable attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ni1",ni1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ni2",ni2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nj1",nj1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nj2",nj2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nk1",nk1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nk2",nk2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ni",ni)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nj",nj)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nk",nk)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nx1",nx1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nx2",nx2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ny1",ny1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ny2",ny2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nz1",nz1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nz2",nz2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nx",nx)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ny",ny)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nz",nz)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngi1",ngi1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngi2",ngi2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngj1",ngj1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngj2",ngj2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngk1",ngk1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngk2",ngk2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngx1",ngx1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngx2",ngx2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngy1",ngy1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngy2",ngy2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngz1",ngz1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngz2",ngz2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"point_in_this",point_in_this)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"start",subs)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"count",subc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"stride",subt)
!--
ierr=nf90_enddef(ncid)
!--
ierr=nf90_close(ncid)
end subroutine nfseis_grid3d_attput

!----------------------------------------------------------------

subroutine nfseis_diminfo(filenm,dimnm,dimlen)
character (len=*),intent(in) :: filenm,dimnm
integer,intent(out) :: dimlen

integer :: ierr,ncid,dimid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'diminfo:'//trim(filenm))
ierr=nf90_inq_dimid(ncid, dimnm, dimid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'dim name in diminfo:'//trim(dimnm))
ierr=nf90_inquire_dimension( ncid, dimid, len=dimlen)
     call nfseis_except(ierr,'inquire dimlen in diminfo')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in diminfo')
end subroutine nfseis_diminfo

! --- put ---
subroutine nfseis_put_int(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int')
     end if
end subroutine nfseis_put_int
subroutine nfseis_put_int1d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int1d')
     end if
end subroutine nfseis_put_int1d
subroutine nfseis_put_int2d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int2d')
     end if
end subroutine nfseis_put_int2d
subroutine nfseis_put_int3d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int3d')
     end if
end subroutine nfseis_put_int3d
subroutine nfseis_put_real(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real')
     end if
end subroutine nfseis_put_real
subroutine nfseis_put_real1d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real1d')
     end if
end subroutine nfseis_put_real1d
subroutine nfseis_put_real2d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real2d')
     end if
end subroutine nfseis_put_real2d
subroutine nfseis_put_real3d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real3d')
     end if
end subroutine nfseis_put_real3d
subroutine nfseis_put_real4d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:,:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real4d')
     end if
end subroutine nfseis_put_real4d

! --- varput ---
subroutine nfseis_varput_int(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int')
end subroutine nfseis_varput_int
subroutine nfseis_varput_int1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int1d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int1d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int1d')
end subroutine nfseis_varput_int1d
subroutine nfseis_varput_int2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int2d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int2d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int2d')
end subroutine nfseis_varput_int2d
subroutine nfseis_varput_int3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int3d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int3d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int3d')
end subroutine nfseis_varput_int3d
subroutine nfseis_varput_real(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real')
end subroutine nfseis_varput_real
subroutine nfseis_varput_real1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real1d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real1d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real1d')
end subroutine nfseis_varput_real1d
subroutine nfseis_varput_real2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real2d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real2d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real2d')
end subroutine nfseis_varput_real2d
subroutine nfseis_varput_real3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real3d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real3d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real3d')
end subroutine nfseis_varput_real3d
subroutine nfseis_varput_real4d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real4d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real4d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real4d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real4d')
end subroutine nfseis_varput_real4d

! --- get var ---
subroutine nfseis_varget_int(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
integer,dimension(1:1) :: a
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,a,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int error:'//trim(filenm))
     end if
var=a(1)
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int')
end subroutine nfseis_varget_int
subroutine nfseis_varget_int1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int1d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int1d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int1d')
end subroutine nfseis_varget_int1d
subroutine nfseis_varget_int2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int2d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int2d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int2d')
end subroutine nfseis_varget_int2d
subroutine nfseis_varget_int3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int3d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int3d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int3d')
end subroutine nfseis_varget_int3d
subroutine nfseis_varget_int4d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int4d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int4d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int4d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int4d')
end subroutine nfseis_varget_int4d
subroutine nfseis_varget_real(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
real(SP),dimension(1:1) :: a
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,a,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real error:'//trim(filenm))
     end if
var=a(1)
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real')
end subroutine nfseis_varget_real
subroutine nfseis_varget_real1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real1d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real1d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real1d')
end subroutine nfseis_varget_real1d
subroutine nfseis_varget_real2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real2d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real2d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real2d')
end subroutine nfseis_varget_real2d
subroutine nfseis_varget_real3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real3d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real3d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real3d')
end subroutine nfseis_varget_real3d
subroutine nfseis_varget_real4d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real4d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real4d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real4d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real4d')
end subroutine nfseis_varget_real4d

! --- get att ---
subroutine nfseis_attget_int(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_int
subroutine nfseis_attget_int1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,dimension(:),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_int1d
subroutine nfseis_attget_real(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_real
subroutine nfseis_attget_real1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),dimension(:),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_real1d
subroutine nfseis_attget_character(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
character (len=SEIS_STRLEN),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_character

! --- put att ---
subroutine nfseis_attput_int(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_int
subroutine nfseis_attput_int1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,dimension(:),intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_int1d
subroutine nfseis_attput_real(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_real
subroutine nfseis_attput_real1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),dimension(:),intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_real1d

!----------------------------------------------------------------

subroutine nfseis_open(filenm,ncid)
character (len=*),intent(in) :: filenm
integer,intent(out) :: ncid
integer :: ierr
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'nfseis_open:'//trim(filenm))
end subroutine nfseis_open
subroutine nfseis_inq_varid(ncid,vnm,vid)
integer,intent(in) :: ncid
integer,intent(out) :: vid
character (len=*),intent(in) :: vnm
integer ierr
ierr=nf90_inq_varid(ncid, vnm, vid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'nfseis_inq_varid:'//trim(vnm))
end subroutine nfseis_inq_varid

subroutine nfseis_close(ncid)
    integer,intent(in) :: ncid
    integer ierr
    ierr=nf90_close(ncid)
end subroutine nfseis_close

subroutine nfseis_except(ierr,msg)
    integer,intent(in) :: ierr
    character (len=*),intent(in) :: msg
    if (ierr /= nf90_noerr) then
       print *, trim(msg)
       print *, trim(nf90_strerror(ierr))
       stop 2
    end if
end subroutine nfseis_except
subroutine nfseis_handle_err(ierr)
    integer,intent(in) :: ierr
    if (ierr /= nf90_noerr) then
       print *, trim(nf90_strerror(ierr))
       stop 2
    end if
end subroutine nfseis_handle_err

end module nfseis_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
