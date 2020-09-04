program seis3d_media

! This program inits the medium of the 3D geology model
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:52:03 -0500 (Fri, 16 Jan 2009) $
! $Revision: 66 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod
use string_mod
use math_mod
use para_mod
use mpi_mod
use nfseis_mod
use grid_mod
use media_mod
use custom_mod
use io_mod
#ifdef MediaMPI
use mpi
#endif

implicit none

integer :: NGRIDX,NGRIDY,NGRIDZ

type STRUCT_1D
     logical :: yes
     integer :: nk
     character (len=SEIS_STRLEN) :: fnm
     character (len=SEIS_STRLEN) :: filetype
     real(SP),dimension(:),pointer :: z
     real(SP),dimension(:),pointer :: d
     real(SP),dimension(:),pointer :: Vp,Vs,Dp
     real(SP),dimension(:),pointer :: Qs
     real(SP) :: QsF0,QsINF
end type STRUCT_1D

type STRUCT_INTERFACE
     logical :: yes
     character (len=SEIS_STRLEN) :: fnm
     character (len=SEIS_STRLEN) :: filetype
     integer :: ni,nj,nk
     real(SP),dimension(:),pointer :: x,y
     real(SP),dimension(:,:,:),pointer :: z
     real(SP),dimension(:,:,:),pointer :: h
     real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
     real(SP),dimension(:,:),pointer :: Qs
     real(SP) :: QsF0,QsINF
end type STRUCT_INTERFACE

type STRUCT_LAYERED
     logical :: yes
     character (len=SEIS_STRLEN) :: fnm
     character (len=SEIS_STRLEN) :: filetype
     integer :: ni,nj,nk
     real(SP),dimension(:),pointer :: x,y
     real(SP),dimension(:,:,:),pointer :: d
     real(SP),dimension(:,:,:),pointer :: h
     real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
     real(SP),dimension(:,:),pointer :: Qs
     real(SP) :: QsF0,QsINF
end type STRUCT_LAYERED

type STRUCT_COMPOSITE
     logical :: yes
     character (len=SEIS_STRLEN) :: fnm
     character (len=SEIS_STRLEN) :: filetype
     integer :: ni,nj,nk
     real(SP),dimension(:),pointer :: x,y
     real(SP),dimension(:,:,:),pointer :: h
     real(SP),dimension(:,:,:,:),pointer :: Vp,Vs,Dp
     real(SP),dimension(:,:,:,:),pointer :: Qs
     real(SP),dimension(:),pointer :: Vp_poly_d,Vs_poly_d,Dp_poly_d,Qs_poly_d
     real(SP) :: QsF0,QsINF
end type STRUCT_COMPOSITE

type STRUCT_VOLUME
     logical :: yes
     character (len=SEIS_STRLEN) :: fnm
     character (len=SEIS_STRLEN) :: filetype
     integer :: imax,jmax,kmax
     integer :: ni,nj,nk
     real(SP) :: rz0
     real(SP),dimension(:),pointer :: gx,gy,grz
     real(SP),dimension(:),pointer :: x,y,rz
     real(SP),dimension(:,:,:),pointer :: Vp,Vs,Dp
     real(SP),dimension(:,:,:),pointer :: Qs
     real(SP) :: QsF0,QsINF
end type STRUCT_VOLUME

type STRUCT_VERPOLY
     logical :: yes
     character (len=SEIS_STRLEN) :: fnm
     character (len=SEIS_STRLEN) :: filetype
     integer :: imax,jmax,kmax,nmax
     integer :: ni,nj,nk
     real(SP) :: rz0
     real(SP),dimension(:),pointer :: gx,gy
     real(SP),dimension(:),pointer :: x,y
     real(SP),dimension(:,:,:),pointer :: rz
     real(SP),pointer :: Vp_poly_c(:,:,:,:),Vp_poly_d(:)
     real(SP),pointer :: Vs_poly_c(:,:,:,:),Vs_poly_d(:)
     real(SP),pointer :: Dp_poly_c(:,:,:,:),Dp_poly_d(:)
     real(SP),pointer :: Qs_poly_c(:,:,:,:),Qs_poly_d(:)
     real(SP) :: QsF0,QsINF
end type STRUCT_VERPOLY

type (STRUCT_1D) :: L1D
type (STRUCT_INTERFACE) :: LA
type (STRUCT_LAYERED) :: L3D
type (STRUCT_COMPOSITE) :: LC
type (STRUCT_VOLUME) :: BV,PV
type (STRUCT_VERPOLY) :: BP,PP

real(SP) :: ztopo

real(DP),dimension(:,:,:),allocatable :: gx,gy,gz
real(SP) :: dtmax,dtmaxVp,dtmaxL
integer,dimension(SEIS_GEO) :: dtindx,dtnode

integer :: n_i,n_j,n_k
#ifdef MediaMPI
integer :: ierr,ROOT
real(SP),dimension(2) :: dtin,dtout
#endif

!-----------------------------------------------------------------------------

#ifdef MediaMPI
call MPI_INIT(ierr)
#endif

call get_conf_name(fnm_conf)

call swmpi_init(fnm_conf)
call para_init(fnm_conf)

#ifdef MediaMPI
call swmpi_cart_creat
call swmpi_reinit_para
call swmpi_datatype
#else
call swmpi_set_gindx(0,0,0)
#endif

call grid_fnm_init(fnm_conf)
call media_fnm_init(fnm_conf)
call media_alloc
call grid_alloc

call alloc_local

!-----------------------------------------------------------------------------

#ifdef MediaMPI
if (masternode)  &
#endif
   print *, 'init media ...'
call init_media(fnm_media_conf)

dtmax=1.0e10
!-----------------------------------------------------------------------------

#ifdef MediaMPI
if (masternode)  &
#endif
   print *, 'calculate effective media ...'

#ifdef MediaMPI
   n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#else
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=dims(3)-1,0,-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
#endif

   call grid_coord_import(n_i,n_j,dims(3)-1); ztopo=z(nk2)
   call grid_coord_import(n_i,n_j,n_k)

  if (BV%yes) call volume_read(BV)
  if (PV%yes) call volume_read(PV)

  if (BP%yes) call verpoly_read(BP)
  if (PP%yes) call verpoly_read(PP)

  call effmedia_eval

  call media_extend

#ifdef MediaMPI
   if (masternode) write(*,*) "exchange media on boundary stencil ..."
   call media_exchange
#endif

#ifdef MediaMPI
   if (masternode)  &
#endif
      print *, '  export media...'
   call media_create(n_i,n_j,n_k)

   call media_stept_check(n_i,n_j,n_k)

#ifndef MediaMPI
end do
end do
end do

   print *, "exchange media on boundary stencil ..."
   call media_exchange
#endif

#ifdef MediaMPI
dtin=(/ dtmax,real(myid,SP) /)
ROOT=0
#ifdef DataTypeDouble
call MPI_REDUCE(dtin,dtout,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,ROOT,SWMPI_COMM,ierr)
#else
call MPI_REDUCE(dtin,dtout,1,MPI_2REAL,MPI_MINLOC,ROOT,SWMPI_COMM,ierr)
#endif
dtmax=dtout(1)
if (myid==ROOT) then
print *, "Maximum allowed time step is", dtmax, ' in thread', nint(dtout(2))
if (dtmax<stept) then
   print *, "Serious Error: stept>dtmax", stept,dtmax
   print *, " occurs on thread", nint(dtout(2))
   !call MPI_ABORT(SWMPI_COMM,1,ierr)
end if
end if
#else
print *, "Maximum allowed time step is", dtmax
write(*,"(a,3i5,a,3i5)") "located on", dtindx,' in thread', dtnode
print *, " Vp and dL are:", dtmaxVp,dtmaxL
if (dtmax<stept) then
   print *, "Serious Error: stept>dtmax", stept,dtmax
   !stop 1
end if
#endif

call media_destroy
call grid_dealloc

#ifdef MediaMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

subroutine init_media(fnm_conf)
character (len=*),intent(in) :: fnm_conf
character (len=SEIS_STRLEN) :: str
integer :: fid

fid=1001
open(fid,file=trim(fnm_conf),status="old")

call string_conf(fid,1,'half_sample_point',2,NGRIDX)
call string_conf(fid,1,'half_sample_point',3,NGRIDY)
call string_conf(fid,1,'half_sample_point',4,NGRIDZ)

L1D%yes=.false.; L3D%yes=.false.;
LA%yes=.false.; LC%yes=.false.;
BP%yes=.false.; BV%yes=.false.;
PP%yes=.false.; PV%yes=.false.;

! background model
call string_conf(fid,1,'background_type',2,str)
select case (trim(str))
case ('cart1d')
     L1D%yes=.true.
     call string_conf(fid,1,'background_format',2,L1D%filetype)
     call string_conf(fid,1,'background_filename',2,L1D%fnm)
     call layer1d_read(L1D)
case ('interface')
     LA%yes=.true.
     call string_conf(fid,1,'background_format',2,LA%filetype)
     call string_conf(fid,1,'background_filename',2,LA%fnm)
     call interface_read(LA)
case ('layered')
     L3D%yes=.true.
     call string_conf(fid,1,'background_format',2,L3D%filetype)
     call string_conf(fid,1,'background_filename',2,L3D%fnm)
     call layered_read(L3D)
case ('composite')
     LC%yes=.true.
     call string_conf(fid,1,'background_format',2,LC%filetype)
     call string_conf(fid,1,'background_filename',2,LC%fnm)
     call composite_read(LC)
case ('volume')
     BV%yes=.true.
     call string_conf(fid,1,'background_format',2,BV%filetype)
     call string_conf(fid,1,'background_filename',2,BV%fnm)
     call volume_init(BV)
case ('verpoly')
     BP%yes=.true.
     call string_conf(fid,1,'background_format',2,BP%filetype)
     call string_conf(fid,1,'background_filename',2,BP%fnm)
     call verpoly_init(BP)
end select

! perturbed model
call string_conf(fid,1,'perturbed_type',2,str)
select case (trim(str))
case ('volume')
     PV%yes=.true.
     call string_conf(fid,1,'perturbed_format',2,PV%filetype)
     call string_conf(fid,1,'perturbed_filename',2,PV%fnm)
     call volume_init(PV)
case ('verpoly')
     PP%yes=.true.
     call string_conf(fid,1,'perturbed_format',2,PP%filetype)
     call string_conf(fid,1,'perturbed_filename',2,PP%fnm)
     call verpoly_init(PP)
case ('none')
end select

close(fid)

!if (W%yes .and. L%yes) then
!   print *, "cann't deal with poly and layered media at same time currently"
!   stop 1
!end if
end subroutine init_media

! -----------------------------  1d structrue ----------------------------
subroutine layer1d_read(L)
type(STRUCT_1D) :: L
integer :: fid,kmax,k
real(DP) :: d2m,v2ms,d2kgm3,d0
character (len=SEIS_STRLEN) :: str,layer_type

if (trim(L%filetype)/='ascii') then
   call error_except('only ascii input accepted for cart1d type')
end if

fid=2001
open(fid,file=trim(L%fnm),status="old")
call string_conf(fid,1,'radial_sampling',2,kmax)
if (kmax<=1) then
   call error_except('sampling point should be larger than 1 for cart1d type')
end if

call string_conf(fid,1,'distance2meter',2,d2m)
call string_conf(fid,1,'velocity2m/s',2,v2ms)
call string_conf(fid,1,'density2kg/m^3',2,d2kgm3)
#ifdef WITHQS
call string_conf(fid,1,'QsINF',2,L%QsINF)
#endif

L%nk=kmax
allocate(L%z(kmax))
allocate(L%d(kmax))
allocate(L%Vp(kmax))
allocate(L%Vs(kmax))
allocate(L%Dp(kmax))
allocate(L%Qs(kmax))

call string_conf(fid,1,'<anchor_data>',1,str)
do k=1,kmax
#ifdef WITHQS
   read(fid,*) L%z(k),L%Dp(k),L%Vp(k),L%Vs(k),L%Qs(k),L%QsF0
#else
   read(fid,*) L%z(k),L%Dp(k),L%Vp(k),L%Vs(k)
#endif
end do
L%Vp=L%Vp*v2ms; L%Vs=L%Vs*v2ms; L%Dp=L%Dp*d2kgm3

! check
if ( any(L%Dp<=0) .or. any(L%Vp<=0) .or. any(L%Vs<0) ) then
   print *, 'media parameter is negtive'
   print *, 'Vp=',L%Vp
   print *, 'Vs=',L%Vs
   print *, 'rho=',L%Dp
   call error_except('cart type read failed')
end if
if ( any(L%Vp**2<=2.0*L%Vs**2) ) then
   print *, 'Vp**2 < 2Vs**2'
   print *, 'Vp=',L%Vp
   print *, 'Vs=',L%Vs
   call error_except('cart type read failed')
end if

call string_conf(fid,1,'radial_meaning',2,layer_type)

select case (trim(layer_type))
case ('axis')
     d0=max(L%z(1),L%z(kmax))
     do k=1,kmax
        L%d(k)=d0-L%z(k)
     end do
case ('thickness')
     L%d(1)=0.0
     do k=2,kmax
        L%d(k)=L%z(k-1)+L%d(k-1)
     end do
case ('depth')
case default
     call error_except("layer_meaning can only take axis, depth or thickness")
end select

if (L%d(1)>L%d(kmax)) then
   L%d =L%d (kmax:1:-1)
   L%Vp=L%Vp(kmax:1:-1)
   L%Vs=L%Vs(kmax:1:-1)
   L%Dp=L%Dp(kmax:1:-1)
end if

L%d=L%d*d2m

close(fid)

#ifdef WITHQS
QsF0=L%QsF0; QsINF=L%QsINF
#endif
end subroutine layer1d_read

! -----------------------------  interface structrue ----------------------------
subroutine interface_read(L)
type(STRUCT_INTERFACE) :: L
integer :: fid,imax,jmax,kmax,i,j,k
real(DP) :: d2m,v2ms,d2kgm3
character (len=SEIS_STRLEN) :: str

if (trim(L%filetype)/='ascii') then
   call error_except('only ascii input accepted for interface type')
end if

fid=2001
open(fid,file=trim(L%fnm),status="old")
call string_conf(fid,1,'number_of_interface',2,kmax)
if (kmax<1) then
   call error_except('at least 1 layer should exist for interface type')
end if

call string_conf(fid,1,'distance2meter',2,d2m)
call string_conf(fid,1,'velocity2m/s',2,v2ms)
call string_conf(fid,1,'density2kg/m^3',2,d2kgm3)
#ifdef WITHQS
call string_conf(fid,1,'QsF0',2,L%QsF0)
call string_conf(fid,1,'QsINF',2,L%QsINF)
#endif
!call string_conf(fid,1,'layer_sealevel',2,layer_sealevel)
call string_conf(fid,1,'horizontal_sampling',2,imax)
call string_conf(fid,1,'horizontal_sampling',3,jmax)
L%ni=imax; L%nj=jmax; L%nk=kmax
allocate(L%x(imax))
allocate(L%y(jmax))
allocate(L%z(imax,jmax,kmax))
allocate(L%h(imax,jmax,kmax-1))
allocate(L%Vp(2,kmax))
allocate(L%Vs(2,kmax))
allocate(L%Dp(2,kmax))
allocate(L%Qs(2,kmax))

call string_conf(fid,1,'<anchor_media>',1,str)
do k=1,kmax
#ifdef WITHQS
   read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k),L%Qs(:,k)
#else
   read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k)
#endif
end do

!L%Vp(1,1)=L%Vp(2,1); L%Vs(1,1)=L%Vs(2,1)
!L%Dp(1,1)=L%Dp(2,1); L%Qs(1,1)=L%Qs(2,1)
L%Vp=L%Vp*v2ms; L%Vs=L%Vs*v2ms; L%Dp=L%Dp*d2kgm3

! check
do k=1,kmax
do i=1,2
   if (k==1 .and. i==1) then
      if (L%Dp(i,k)<0 .or. L%Vp(i,k)<0 .or. L%Vs(i,k)<0 ) then
          call error_except('model parameter should >=0 for (1,1)')
      end if
   elseif ( L%Dp(i,k)<=0 .or. L%Vp(i,k)<=0 .or. L%Vs(i,k)<0 ) then
      print *, 'media parameter is negtive'
      print *, 'Vp=',L%Vp
      print *, 'Vs=',L%Vs
      print *, 'Dp=',L%Dp
      call error_except('interface type read failed')
   elseif ( L%Vp(i,k)**2<=2.0*L%Vs(i,k)**2 ) then
      print *, 'Vp**2 < 2Vs**2'
      print *, 'Vp=',L%Vp
      print *, 'Vs=',L%Vs
      call error_except('interface type read failed')
   end if
end do
end do

call string_conf(fid,1,'<anchor_interface>',1,str)
do j=1,jmax
do i=1,imax     
    read(fid,*) L%x(i),L%y(j),( L%z(i,j,k), k=1,kmax )
end do 
end do

do k=1,kmax-1
   L%h(:,:,k)=L%z(:,:,k)-L%z(:,:,k+1)
end do

do k=1,kmax-1
   if (any(L%h(:,:,k)<0.0)) then
      print *, 'thickness of layer',k,' should be not less than 0'
      print *, minloc(L%h(:,:,k))
      call error_except("thickness btween two interfaces shouldn' be less than 0")
   end if
end do

L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP
L%z=L%z*d2m; L%h=L%h*d2m

close(fid)

#ifdef WITHQS
QsF0=L%QsF0; QsINF=L%QsINF
#endif
end subroutine interface_read

! -----------------------------  layered structrue ----------------------------
subroutine layered_read(L)
type(STRUCT_LAYERED) :: L
integer :: fid,imax,jmax,kmax,i,j,k
real(DP) :: d2m,v2ms,d2kgm3
character (len=SEIS_STRLEN) :: str,layer_type

if (trim(L%filetype)/='ascii') then
   call error_except('only ascii input accepted for layered type')
end if

fid=2001
open(fid,file=trim(L%fnm),status="old")
call string_conf(fid,1,'number_of_layer',2,kmax)
if (kmax<1) then
   call error_except('at least 1 layer should exist for layered type')
end if

call string_conf(fid,1,'distance2meter',2,d2m)
call string_conf(fid,1,'velocity2m/s',2,v2ms)
call string_conf(fid,1,'density2kg/m^3',2,d2kgm3)
#ifdef WITHQS
call string_conf(fid,1,'QsF0',2,L%QsF0)
call string_conf(fid,1,'QsINF',2,L%QsINF)
#endif
!call string_conf(fid,1,'layer_sealevel',2,layer_sealevel)
call string_conf(fid,1,'horizontal_sampling',2,imax)
call string_conf(fid,1,'horizontal_sampling',3,jmax)
L%ni=imax; L%nj=jmax; L%nk=kmax
allocate(L%x(imax))
allocate(L%y(jmax))
allocate(L%d(imax,jmax,kmax))
allocate(L%h(imax,jmax,kmax))
allocate(L%Vp(2,kmax))
allocate(L%Vs(2,kmax))
allocate(L%Dp(2,kmax))
allocate(L%Qs(2,kmax))

call string_conf(fid,1,'<anchor_media>',1,str)
do k=1,kmax
#ifdef WITHQS
   read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k),L%Qs(:,k)
#else
   read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k)
#endif
end do

L%Vp=L%Vp*v2ms; L%Vs=L%Vs*v2ms; L%Dp=L%Dp*d2kgm3

! check
if ( any(L%Dp<=0) .or. any(L%Vp<=0) .or. any(L%Vs<0) ) then
   print *, 'media parameter is negtive'
   print *, 'Vp=',L%Vp
   print *, 'Vs=',L%Vs
   print *, 'Dp=',L%Dp
   call error_except('layered type read failed')
end if
if ( any(L%Vp**2<=2.0*L%Vs**2) ) then
   print *, 'Vp**2 < 2Vs**2'
   print *, 'Vp=',L%Vp
   print *, 'Vs=',L%Vs
   call error_except('layered type read failed')
end if

call string_conf(fid,1,'layer_meaning',2,layer_type)
call string_conf(fid,1,'<anchor_layer>',1,str)
select case (trim(layer_type))
case ('depth')
  do j=1,jmax
  do i=1,imax     
      read(fid,*) L%x(i),L%y(j),( L%d(i,j,k), k=1,kmax )
  end do 
  end do
  L%h(:,:,1)=L%d(:,:,1)
  do k=2,kmax
     L%h(:,:,k)=L%d(:,:,k)-L%d(:,:,k-1)
  end do
case ('thickness')
  do j=1,jmax
  do i=1,imax     
      read(fid,*) L%x(i),L%y(j),( L%h(i,j,k), k=1,kmax )
  end do 
  end do
  L%d(:,:,1)=L%h(:,:,1)
  do k=2,kmax
     L%d(:,:,k)=L%d(:,:,k-1)+L%h(:,:,k)
  end do
case default
     call error_except("layer_meaning can only take depth or thickness")
end select

do k=1,kmax-1
   if (any(L%h(:,:,k)<0.0)) then
      print *, minloc(L%h(:,:,k))
      call error_except('thickness of layer should not be less than 0')
   end if
end do
if (any(L%h(:,:,kmax)<SEIS_ZERO)) then
   print *, minloc(L%h(:,:,kmax))
   call error_except('thickness of lowest layer should larger than 0')
end if

L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP
L%d=L%d*d2m; L%h=L%h*d2m

close(fid)

#ifdef WITHQS
QsF0=L%QsF0; QsINF=L%QsINF
#endif
end subroutine layered_read

! -----------------------------  composite structrue ----------------------------
subroutine composite_read(L)
type(STRUCT_COMPOSITE) :: L
integer :: imax,jmax,kmax

if (trim(L%filetype)/='nc') then
   call error_except('only nc input accepted for composite type')
end if

call nfseis_diminfo(L%fnm,'theta',imax)
call nfseis_diminfo(L%fnm,'phi',jmax)
call nfseis_diminfo(L%fnm,'layer',kmax)

L%ni=imax; L%nj=jmax; L%nk=kmax
allocate(L%x(imax))
allocate(L%y(jmax))
allocate(L%h(imax,jmax,kmax))
allocate(L%Vp(2,imax,jmax,kmax))
allocate(L%Vs(2,imax,jmax,kmax))
allocate(L%Dp(2,imax,jmax,kmax))
allocate(L%Vp_poly_d(kmax))
allocate(L%Vs_poly_d(kmax))
allocate(L%Dp_poly_d(kmax))
#ifdef WITHQS
allocate(L%Qs(2,imax,jmax,kmax))
allocate(L%Qs_poly_d(kmax))
#endif

call nfseis_varget(L%fnm,'theta',L%x,(/1/),(/imax/),(/1/))
call nfseis_varget(L%fnm,'phi',L%y,(/1/),(/jmax/),(/1/))
call nfseis_varget(L%fnm,'thickness',L%h,(/1,1,1/),(/imax,jmax,kmax/),(/1,1,1/) )
call nfseis_varget(L%fnm,'Vp',L%Vp,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
call nfseis_varget(L%fnm,'Vs',L%Vs,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
call nfseis_varget(L%fnm,'rho',L%Dp,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
call nfseis_varget(L%fnm,'Vp_poly_d',L%Vp_poly_d,(/1/),(/kmax/),(/1/) )
call nfseis_varget(L%fnm,'Vs_poly_d',L%Vs_poly_d,(/1/),(/kmax/),(/1/) )
call nfseis_varget(L%fnm,'rho_poly_d',L%Dp_poly_d,(/1/),(/kmax/),(/1/) )
#ifdef WITHQS
! - ys 09/07/2010
call nfseis_varget(L%fnm,'Qs_poly_d',L%Qs_poly_d,(/1/),(/kmax/),(/1/) )
!
call nfseis_varget(L%fnm,'Qs',L%Qs,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
call nfseis_attget(L%fnm,'QsF0',L%QsF0)
call nfseis_attget(L%fnm,'QsINF',L%QsINF)
QsF0=L%QsF0; QsINF=L%QsINF
#endif

L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP

if (any(L%h(:,:,1:kmax-1)<0.0)) then
   call error_except('thickness of layer should not be less than 0')
end if
if (any(L%h(:,:,kmax)<SEIS_ZERO)) then
   print *, minloc(L%h(:,:,kmax))
   call error_except('thickness of lowest layer should larger than 0')
end if
end subroutine composite_read

! -----------------------------  volume structrue ----------------------------
subroutine volume_init(P)
type(STRUCT_VOLUME) :: P
integer :: imax,jmax,kmax

if (trim(P%filetype)/='nc') then
   call error_except('only nc input accepted for volume type')
end if

  call nfseis_diminfo(P%fnm,'theta',imax)
  call nfseis_diminfo(P%fnm,'phi',jmax)
  call nfseis_diminfo(P%fnm,'depth',kmax)
  allocate(P%gx(imax)); P%gx=0.0
  allocate(P%gy(jmax)); P%gy=0.0
  allocate(P%grz(kmax)); P%grz=0.0
  call nfseis_varget(P%fnm,'theta',P%gx,(/1/),(/imax/),(/1/))
  call nfseis_varget(P%fnm,'phi',P%gy,(/1/),(/jmax/),(/1/))
  call nfseis_varget(P%fnm,'depth2sealevel',P%grz,(/1/),(/kmax/),(/1/))
  call nfseis_attget(P%fnm,'sealevel',P%rz0)
  P%gx=P%gx*PI/180.0_DP; P%gy=P%gy*PI/180.0_DP
  P%imax=imax; P%jmax=jmax; P%kmax=kmax
  imax=min(nx*max(2*NGRIDX,1),imax)
  jmax=min(ny*max(2*NGRIDY,1),jmax)
  kmax=min(nz*max(2*NGRIDZ,1),kmax)
  call volume_alloc(P,imax,jmax,kmax)
end subroutine volume_init

subroutine volume_alloc(P,imax,jmax,kmax)
type(STRUCT_VOLUME) :: P
integer,intent(in) :: imax,jmax,kmax

if (associated(P%Vp)) then
if (size(P%Vp,1)<imax .or. size(P%Vp,2)<jmax .or. size(P%Vp,3)<kmax) then
   deallocate(P%x);deallocate(P%y);deallocate(P%rz);
   deallocate(P%Vp); deallocate(P%Vs); deallocate(P%Dp)
end if
end if

P%ni=imax;P%nj=jmax;P%nk=kmax
if (.not. associated(P%Vp)) then
   allocate(P%x(imax)); P%x=0.0
   allocate(P%y(jmax)); P%y=0.0
   allocate(P%rz(kmax)); P%rz=0.0
   allocate(P%Vp(imax,jmax,kmax)); P%Vp=0.0
   allocate(P%Vs(imax,jmax,kmax)); P%Vs=0.0
   allocate(P%Dp(imax,jmax,kmax)); P%Dp=0.0
end if
end subroutine volume_alloc

subroutine volume_read(P)
type(STRUCT_VOLUME) :: P
integer :: i1,i2,j1,j2,k1,k2,imax,jmax,kmax
integer :: indx(1)
indx=maxloc(P%gx,P%gx<=max(x(nx1),P%gx(1)));i1=max(indx(1),1)
indx=minloc(P%gx,P%gx>=min(x(nx2),P%gx(P%imax)));i2=indx(1)
indx=maxloc(P%gy,P%gy<=max(y(ny1),P%gy(1)));j1=max(indx(1),1)
indx=minloc(P%gy,P%gy>=min(y(ny2),P%gy(P%jmax)));j2=indx(1)
k1=1; k2=P%kmax
imax=i2-i1+1; jmax=j2-j1+1; kmax=k2-k1+1
call volume_alloc(P,imax,jmax,kmax)

P%x(1:imax)=P%gx(i1:i2); P%y(1:jmax)=P%gy(j1:j2)
P%rz(1:kmax)=P%grz(k1:k2);
call nfseis_varget( P%fnm,'Vp',P%Vp(1:imax,1:jmax,1:kmax), &
     (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
call nfseis_varget( P%fnm,'Vs',P%Vs(1:imax,1:jmax,1:kmax), &
     (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
call nfseis_varget( P%fnm,'rho',P%Dp(1:imax,1:jmax,1:kmax), &
     (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
end subroutine volume_read

! ---------------------------  verpoly structrue ----------------------------
subroutine verpoly_init(P)
type(STRUCT_VERPOLY) :: P
integer :: imax,jmax,kmax,nmax

if (trim(P%filetype)/='nc') then
   call error_except('only nc input accepted for verpoly type')
end if

   call nfseis_diminfo(P%fnm,'theta',imax)
   call nfseis_diminfo(P%fnm,'phi',jmax)
   call nfseis_diminfo(P%fnm,'depth',kmax)
   call nfseis_diminfo(P%fnm,'polynomial',nmax)
   allocate(P%gx(imax)); P%gx=0.0
   allocate(P%gy(jmax)); P%gy=0.0
   allocate(P%Vp_poly_d(nmax)); P%Vp_poly_d=0.0
   allocate(P%Vs_poly_d(nmax)); P%Vs_poly_d=0.0
   allocate(P%Dp_poly_d(nmax)); P%Dp_poly_d=0.0
   call nfseis_varget(P%fnm,'theta',P%gx,(/1/),(/imax/),(/1/))
   call nfseis_varget(P%fnm,'phi',P%gy,(/1/),(/jmax/),(/1/))
   call nfseis_varget(P%fnm,'Vp_poly_d',P%Vp_poly_d,(/1/),(/nmax/),(/1/))
   call nfseis_varget(P%fnm,'Vs_poly_d',P%Vs_poly_d,(/1/),(/nmax/),(/1/))
   call nfseis_varget(P%fnm,'rho_poly_d',P%Dp_poly_d,(/1/),(/nmax/),(/1/))
   call nfseis_attget(P%fnm,'sealevel',P%rz0)
   P%gx=P%gx*PI/180.0_DP; P%gy=P%gy*PI/180.0_DP
   P%imax=imax; P%jmax=jmax; P%kmax=kmax; P%nmax=nmax
   imax=min(nx*max(2*NGRIDX,1),imax)
   jmax=min(ny*max(2*NGRIDY,1),jmax)
   kmax=min(nz*max(2*NGRIDZ,1),kmax)
   call verpoly_alloc(P,imax,jmax,kmax)
end subroutine verpoly_init

subroutine verpoly_alloc(P,imax,jmax,kmax)
type(STRUCT_VERPOLY) :: P
integer,intent(in) :: imax,jmax,kmax

if (associated(P%rz)) then
if (size(P%rz,1)<imax .or. size(P%rz,2)<jmax .or. size(P%rz,3)<kmax) then
   deallocate(P%x);deallocate(P%y);deallocate(P%rz);
   deallocate(P%Vp_poly_c); deallocate(P%Vs_poly_c); deallocate(P%Dp_poly_c)
end if
end if

P%ni=imax;P%nj=jmax;P%nk=kmax
if (.not. associated(P%rz)) then
   allocate(P%x(imax)); P%x=0.0
   allocate(P%y(jmax)); P%y=0.0
   allocate(P%rz(imax,jmax,kmax)); P%rz=0.0
   allocate(P%Vp_poly_c(imax,jmax,kmax,P%nmax)); P%Vp_poly_c=0.0
   allocate(P%Vs_poly_c(imax,jmax,kmax,P%nmax)); P%Vs_poly_c=0.0
   allocate(P%Dp_poly_c(imax,jmax,kmax,P%nmax)); P%Dp_poly_c=0.0
end if
end subroutine verpoly_alloc

subroutine verpoly_read(P)
type(STRUCT_VERPOLY) :: P
integer :: i1,i2,j1,j2,k1,k2,imax,jmax,kmax,nmax
integer :: indx(1)
indx=maxloc(P%gx,P%gx<=max(x(nx1),P%gx(1)));i1=max(indx(1),1)
indx=minloc(P%gx,P%gx>=min(x(nx2),P%gx(P%imax)));i2=indx(1)
indx=maxloc(P%gy,P%gy<=max(y(ny1),P%gy(1)));j1=max(indx(1),1)
indx=minloc(P%gy,P%gy>=min(y(ny2),P%gy(P%jmax)));j2=indx(1)
k1=1; k2=P%kmax
imax=i2-i1+1; jmax=j2-j1+1; kmax=k2-k1+1; nmax=P%nmax
call verpoly_alloc(P,imax,jmax,kmax)

P%x(1:imax)=P%gx(i1:i2); P%y(1:jmax)=P%gy(j1:j2)
call nfseis_varget( P%fnm,'depth2sealevel',P%rz(1:imax,1:jmax,1:kmax), &
     (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
call nfseis_varget( P%fnm,'Vp_poly_c',P%Vp_poly_c(1:imax,1:jmax,1:kmax,:), &
     (/i1,j1,k1,1/),(/imax,jmax,kmax,nmax/),(/1,1,1,1/) )
call nfseis_varget( P%fnm,'Vs_poly_c',P%Vs_poly_c(1:imax,1:jmax,1:kmax,:), &
     (/i1,j1,k1,1/),(/imax,jmax,kmax,nmax/),(/1,1,1,1/) )
call nfseis_varget( P%fnm,'rho_poly_c',P%Dp_poly_c(1:imax,1:jmax,1:kmax,:), &
     (/i1,j1,k1,1/),(/imax,jmax,kmax,nmax/),(/1,1,1,1/) )
end subroutine verpoly_read

!---------------------------------

subroutine effmedia_eval

real(SP),dimension(-NGRIDX:NGRIDX) :: xvec
real(SP),dimension(-NGRIDY:NGRIDY) :: yvec
real(SP),dimension(-NGRIDZ:NGRIDZ) :: zvec
real(SP) :: dx1,dx2,dy1,dy2,dz1,dz2

real(SP),dimension(2) :: lam, miu
real(SP),dimension(2) :: Vp,Vs,Dp,dVp,dVs,dDp
real(SP) :: rho0,mu0,lam0
logical :: flag_water
integer :: nwater

real(SP) :: x0,y0,z0,d0,x1,y1,z1,z2
integer :: nsamp
integer :: i,j,k,mi,mj,mk,n,m
integer :: i1,i2,j1,j2,k1,k2,k0,n1,n2
integer :: indx(1)

!layered
real(SP) :: h,d
real(SP) :: L1,L2
!verpoly
integer :: wi,wj,wk

!imax=L%ni; jmax=L%nj; kmax=L%nk
nsamp=max(2*NGRIDX,1)*max(2*NGRIDY,1)*max(2*NGRIDZ,1)

do k=nk1,nk2
#ifdef MediaMPI
   if (masternode)  &
#endif
       print *, ' k=',k-nk1+1, ' of',nk
do j=nj1,nj2
do i=ni1,ni2

   rho0=0.0;mu0=0.0;lam0=0.0
   flag_water=.false.; nwater=0

   xvec(0)=x(i); yvec(0)=y(j); zvec(0)=z(k)

   dx1=(x(i)-x(i-1))/2.0/(NGRIDX+1); dx2=(x(i+1)-x(i))/2.0/(NGRIDX+1)
   xvec(-NGRIDX:-1)=(/-NGRIDX:-1/)*dx1+x(i);
   xvec(1:NGRIDX)  =(/ 1:NGRIDX /)*dx2+x(i)

   dy1=(y(j)-y(j-1))/2.0/(NGRIDY+1); dy2=(y(j+1)-y(j))/2.0/(NGRIDY+1)
   yvec(-NGRIDY:-1)=(/-NGRIDY:-1/)*dy1+y(j);
   yvec( 1:NGRIDY )=(/ 1:NGRIDY /)*dy2+y(j)

   dz1=(z(k)-z(k-1))/2.0/(NGRIDZ+1); dz2=(z(k+1)-z(k))/2.0/(NGRIDZ+1)
   zvec(-NGRIDZ:-1)=(/-NGRIDZ:-1/)*dz1+z(k);
   zvec( 1:NGRIDZ )=(/ 1:NGRIDZ /)*dz2+z(k)

do mk=-NGRIDZ,NGRIDZ
   if (NGRIDZ/=0 .and. mk==0) cycle
do mj=-NGRIDY,NGRIDY
   if (NGRIDY/=0 .and. mj==0) cycle
do mi=-NGRIDX,NGRIDX
   if (NGRIDX/=0 .and. mi==0) cycle

   x0=xvec(mi); y0=yvec(mj); z0=zvec(mk);
   d0=max(ztopo-z0,0.0);

! -------------   background model ---------------------------
! cart1d
if (L1D%yes) then
   k1=1; k2=L1D%nk
   if (abs(L1D%d(k1)-d0)<SEIS_ZERO .or. d0<L1D%d(k1)) then
      Vp(:)=L1D%Vp(k1); Vs(:)=L1D%Vs(k1); Dp(:)=L1D%Dp(k1)
   elseif (abs(L1D%d(k2)-d0)<SEIS_ZERO) then
      if (abs(L1D%d(k2-1)-L1D%d(k2))<SEIS_ZERO) then
         Vp(1)=L1D%Vp(k2-1);Vp(2)=L1D%Vp(k2);
         Vs(1)=L1D%Vs(k2-1);Vs(2)=L1D%Vs(k2)
         Dp(1)=L1D%Dp(k2-1);Dp(2)=L1D%Dp(k2)
      else
         Vp(:)=L1D%Vp(k2); Vs(:)=L1D%Vs(k2); Dp(:)=L1D%Dp(k2)
      end if
   elseif (L1D%d(k2)<d0) then
      Vp(:)=L1D%Vp(k2); Vs(:)=L1D%Vs(k2); Dp(:)=L1D%Dp(k2)
   else
       do while (k2-k1>1)
          k0=(k2-k1)/2+k1
          if (abs(L1D%d(k0)-d0)<SEIS_ZERO) then
             k1=k0;k2=k0; exit
          elseif (L1D%d(k0)>d0) then
             k2=k0;
          else
             k1=k0
          end if
       end do

       if (k1==k2) then
          if (abs(L1D%d(k1-1)-L1D%d(k1))<SEIS_ZERO) then
             Vp(1)=L1D%Vp(k1-1); Vp(2)=L1D%Vp(k1);
             Vs(1)=L1D%Vs(k1-1); Vs(2)=L1D%Vs(k1)
             Dp(1)=L1D%Dp(k1-1); Dp(2)=L1D%Dp(k1)
          elseif (abs(L1D%d(k1+1)-L1D%d(k1))<SEIS_ZERO) then
             Vp(1)=L1D%Vp(k1); Vp(2)=L1D%Vp(k1+1);
             Vs(1)=L1D%Vs(k1); Vs(2)=L1D%Vs(k1+1)
             Dp(1)=L1D%Dp(k1); Dp(2)=L1D%Dp(k1+1)
          else
             Vp(:)=L1D%Vp(k1); Vs(:)=L1D%Vs(k1); Dp(:)=L1D%Dp(k1)
          end if
       else
          d=L1D%d(k2); !-d0;
          h=L1D%d(k2)-L1D%d(k1)
          L1=(d-d0)/h; L2=1.0-L1
          Vp(:)=L1D%Vp(k1)*L1+L1D%Vp(k2)*L2
          Vs(:)=L1D%Vs(k1)*L1+L1D%Vs(k2)*L2
          Dp(:)=L1D%Dp(k1)*L1+L1D%Dp(k2)*L2
       end if
   end if
end if
! interface
if (LA%yes) then
   call indx_locate_1d(x0,LA%x,i1,i2,x1)
   call indx_locate_1d(y0,LA%y,j1,j2,y1)

   do n=LA%nk,1,-1

      z1=interp_2d(LA%x(i1:i2),LA%y(j1:j2),LA%z(i1:i2,j1:j2,n),2,2,x1,y1)
      if (n==LA%nk) then
         h=1.0e3
      else
         h =interp_2d(LA%x(i1:i2),LA%y(j1:j2),LA%h(i1:i2,j1:j2,n),2,2,x1,y1)
      end if

      if (h<=SEIS_ZERO .and. n>1) cycle

      if (h<=SEIS_ZERO .and. n==1) then
         if (LA%Vp(1,n)>SEIS_ZERO) then
            n1=1;k1=1;n2=1;k2=1
         end if
         Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
         Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
         Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
         exit
      end if

      if (abs(z1-z0)<=dz2/5.0) then  !  .and. h>SEIS_ZERO
         k2=n; n2=2
         if (n==1) then
            if (LA%Vp(1,n)>SEIS_ZERO) then
               k1=n; n1=1
            else
               k1=n; n1=2
            end if
         end if
         do m=n-1,1,-1
            h=interp_2d(LA%x(i1:i2),LA%y(j1:j2),LA%h(i1:i2,j1:j2,m),2,2,x1,y1)
            if (h>SEIS_ZERO) then
                k1=m+1; n1=1
                exit
            elseif (m==1) then
               if (LA%Vp(1,m)>SEIS_ZERO) then
                  k1=m; n1=1
               else
                  k1=n; n1=2
               end if
            end if
         end do
         Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
         Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
         Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
         exit
      end if

      if (z0<z1) then
         if (n==LA%nk) then
            Vp(:)=LA%Vp(2,n)
            Vs(:)=LA%Vs(2,n)
            Dp(:)=LA%Dp(2,n)
         else
            L2=(z1-z0)/h; L1=1.0-L2
            Vp(:)=LA%Vp(2,n)*L1+LA%Vp(1,n+1)*L2
            Vs(:)=LA%Vs(2,n)*L1+LA%Vs(1,n+1)*L2
            Dp(:)=LA%Dp(2,n)*L1+LA%Dp(1,n+1)*L2
         end if
         exit
      end if

      k1=n; k2=n; n1=2; n2=2

      if (n==1) then
         if (LA%Vp(1,n)>SEIS_ZERO) then
            n1=1;k1=1;n2=1;k2=1
         end if
         Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
         Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
         Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
         exit
      end if

   end do
end if
! layered
if (L3D%yes) then
   d=0.0
   call indx_locate_1d(x0,L3D%x,i1,i2,x1)
   call indx_locate_1d(y0,L3D%y,j1,j2,y1)
   do n=1,L3D%nk
      h=interp_2d(L3D%x(i1:i2),L3D%y(j1:j2),L3D%h(i1:i2,j1:j2,n),2,2,x1,y1)

      if (h<=SEIS_ZERO) cycle

      d=d+h
      if (abs(d-d0)<=dz2/5.0) then
         if (n==L3D%nk) call error_except('please increase the thickness of the lowest layer')
         n1=n;
         do m=n+1,L3D%nk
            h=interp_2d(L3D%x(i1:i2),L3D%y(j1:j2),L3D%h(i1:i2,j1:j2,m),2,2,x1,y1)
            if (h>SEIS_ZERO) then
                n2=m
                exit
            end if
         end do
         Vp(1)=L3D%Vp(2,n1); Vp(2)=L3D%Vp(1,n2)
         Vs(1)=L3D%Vs(2,n1); Vs(2)=L3D%Vs(1,n2)
         Dp(1)=L3D%Dp(2,n1); Dp(2)=L3D%Dp(1,n2)
         exit
      elseif (d>d0) then
         L1=(d-d0)/h; L2=1.0-L1
         Vp(:)=L3D%Vp(1,n)*L1+L3D%Vp(2,n)*L2
         Vs(:)=L3D%Vs(1,n)*L1+L3D%Vs(2,n)*L2
         Dp(:)=L3D%Dp(1,n)*L1+L3D%Dp(2,n)*L2
         exit
      elseif (n==L3D%nk) then
         call error_except('please increase the thickness of the lowest layer')
      end if
   end do
end if
! composite
if (LC%yes) then
   d=0.0
   call indx_locate_1d(x0,LC%x,i1,i2,x1)
   call indx_locate_1d(y0,LC%y,j1,j2,y1)
   do n=1,LC%nk
      h=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%h(i1:i2,j1:j2,n),2,2,x1,y1)

      if (h<=SEIS_ZERO) cycle

      d=d+h
      if (abs(d-d0)<=dz2/5.0) then
         if (n==LC%nk) call error_except('please increase the thickness of the lowest layer')
         n1=n
         do m=n+1,LC%nk
            h=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%h(i1:i2,j1:j2,m),2,2,x1,y1)
            if (h>SEIS_ZERO) then
                n2=m
                exit
            end if
         end do
         Vp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(2,i1:i2,j1:j2,n1),2,2,x1,y1)
         Vs(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(2,i1:i2,j1:j2,n1),2,2,x1,y1)
         Dp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(2,i1:i2,j1:j2,n1),2,2,x1,y1)
         Vp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(1,i1:i2,j1:j2,n2),2,2,x1,y1)
         Vs(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(1,i1:i2,j1:j2,n2),2,2,x1,y1)
         Dp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(1,i1:i2,j1:j2,n2),2,2,x1,y1)
         exit
      elseif (d>d0) then
         Vp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(1,i1:i2,j1:j2,n),2,2,x1,y1)
         Vs(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(1,i1:i2,j1:j2,n),2,2,x1,y1)
         Dp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(1,i1:i2,j1:j2,n),2,2,x1,y1)
         Vp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(2,i1:i2,j1:j2,n),2,2,x1,y1)
         Vs(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(2,i1:i2,j1:j2,n),2,2,x1,y1)
         Dp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(2,i1:i2,j1:j2,n),2,2,x1,y1)
         Vp(:)=Vp(1)+(Vp(2)-Vp(1))/h**LC%Vp_poly_d(n)*(d0-(d-h))**LC%Vp_poly_d(n)
         Vs(:)=Vs(1)+(Vs(2)-Vs(1))/h**LC%Vs_poly_d(n)*(d0-(d-h))**LC%Vs_poly_d(n)
         Dp(:)=Dp(1)+(Dp(2)-Dp(1))/h**LC%Dp_poly_d(n)*(d0-(d-h))**LC%Dp_poly_d(n)
         exit
      elseif (n==LC%nk) then
         call error_except('please increase the thickness of the lowest layer')
      end if
   end do
end if
! volume
if (BV%yes) then
   call indx_locate_1d(x0,BV%x(1:BV%ni),i1,i2,x1)
   call indx_locate_1d(y0,BV%y(1:BV%nj),j1,j2,y1)
   call indx_locate_1d(BV%rz0-z0,BV%rz(1:BV%nk),k1,k2,z1)
   Vp(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
         BV%Vp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
   Vs(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
         BV%Vs(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
   Dp(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
         BV%Dp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
end if
! verpoly
if (BP%yes) then
   indx=minloc(abs(BP%x(1:BP%ni)-x0)); wi=indx(1)
   indx=minloc(abs(BP%y(1:BP%nj)-y0)); wj=indx(1)
   !if (abs(BP%x(wi)-x0)+abs(BP%y(wj)-y0)>SEIS_ZERO) then
   !   print *, 'error'
   !   stop 1
   !else
   indx=minloc(abs(BP%rz(wi,wj,1:BP%nk)-(BP%rz0-z0))); wk=indx(1)
   if (abs(BP%rz(wi,wj,wk)-(BP%rz0-z0))<=dz2/5.0  &
      .or. (BP%rz0-z0)<BP%rz(wi,wj,1)) then
      k1=max(1,wk-1); k2=wk
      z1=BP%rz(wi,wj,k2)-BP%rz(wi,wj,k1); z2=0
      Vp(1)=eval_poly(BP%Vp_poly_c(wi,wj,k1,:),BP%Vp_poly_d,z1)
      Vp(2)=eval_poly(BP%Vp_poly_c(wi,wj,k2,:),BP%Vp_poly_d,z2)
      Vs(1)=eval_poly(BP%Vs_poly_c(wi,wj,k1,:),BP%Vs_poly_d,z1)
      Vs(2)=eval_poly(BP%Vs_poly_c(wi,wj,k2,:),BP%Vs_poly_d,z2)
      Dp(1)=eval_poly(BP%Dp_poly_c(wi,wj,k1,:),BP%Dp_poly_d,z1)
      Dp(2)=eval_poly(BP%Dp_poly_c(wi,wj,k2,:),BP%Dp_poly_d,z2)
   else
      if (BP%rz0-z0<BP%rz(wi,wj,wk)) wk=wk-1
      z1=(BP%rz0-z0)-BP%rz(wi,wj,wk)
      Vp(:)=eval_poly(BP%Vp_poly_c(wi,wj,wk,:),BP%Vp_poly_d,z1)
      Vs(:)=eval_poly(BP%Vs_poly_c(wi,wj,wk,:),BP%Vs_poly_d,z1)
      Dp(:)=eval_poly(BP%Dp_poly_c(wi,wj,wk,:),BP%Dp_poly_d,z1)
   end if
end if

! -------------   perturbed model ---------------------------
! volume
if (PV%yes) then
   call indx_locate_1d(x0,PV%x(1:PV%ni),i1,i2,x1)
   call indx_locate_1d(y0,PV%y(1:PV%nj),j1,j2,y1)
   call indx_locate_1d(PV%rz0-z0,PV%rz(1:PV%nk),k1,k2,z1)
   dVp=interp_3d(PV%x(i1:i2),PV%y(j1:j2),PV%rz(k1:k2), &
         PV%Vp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
   dVs=interp_3d(PV%x(i1:i2),PV%y(j1:j2),PV%rz(k1:k2), &
         PV%Vs(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
   dDp=interp_3d(PV%x(i1:i2),PV%y(j1:j2),PV%rz(k1:k2), &
         PV%Dp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
   Vp=Vp*(1.0+dVp)
   Vs=Vs*(1.0+dVs)
   Dp=Dp*(1.0+dDp)
end if
! verpoly
if (PP%yes) then
   indx=minloc(abs(PP%x(1:PP%ni)-x0)); wi=indx(1)
   indx=minloc(abs(PP%y(1:PP%nj)-y0)); wj=indx(1)
   !if (abs(PP%x(wi)-x0)+abs(PP%y(wj)-y0)>SEIS_ZERO) then
   !   print *, 'error'
   !   stop 1
   !else
   indx=minloc(abs(PP%rz(wi,wj,1:PP%nk)-(PP%rz0-z0))); wk=indx(1)
   if (abs(PP%rz(wi,wj,wk)-(PP%rz0-z0))<=dz2/5.0  &
      .or. (PP%rz0-z0)<PP%rz(wi,wj,1)) then
      k1=max(1,wk-1); k2=wk
      z1=PP%rz(wi,wj,k2)-PP%rz(wi,wj,k1); z2=0
      dVp(1)=eval_poly(PP%Vp_poly_c(wi,wj,k1,:),PP%Vp_poly_d,z1)
      dVp(2)=eval_poly(PP%Vp_poly_c(wi,wj,k2,:),PP%Vp_poly_d,z2)
      dVs(1)=eval_poly(PP%Vs_poly_c(wi,wj,k1,:),PP%Vs_poly_d,z1)
      dVs(2)=eval_poly(PP%Vs_poly_c(wi,wj,k2,:),PP%Vs_poly_d,z2)
      dDp(1)=eval_poly(PP%Dp_poly_c(wi,wj,k1,:),PP%Dp_poly_d,z1)
      dDp(2)=eval_poly(PP%Dp_poly_c(wi,wj,k2,:),PP%Dp_poly_d,z2)
   else
      if (PP%rz0-z0<PP%rz(wi,wj,wk)) wk=wk-1
      z1=(PP%rz0-z0)-PP%rz(wi,wj,wk)
      dVp(:)=eval_poly(PP%Vp_poly_c(wi,wj,wk,:),PP%Vp_poly_d,z1)
      dVs(:)=eval_poly(PP%Vs_poly_c(wi,wj,wk,:),PP%Vs_poly_d,z1)
      dDp(:)=eval_poly(PP%Dp_poly_c(wi,wj,wk,:),PP%Dp_poly_d,z1)
   end if
   Vp=Vp*(1.0+dVp)
   Vs=Vs*(1.0+dVs)
   Dp=Dp*(1.0+dDp)
end if

   ! to lame parameters
   miu=Dp*Vs*Vs
   lam=Vp*Vp*Dp - 2.0*miu

   !accumulate
   rho0=rho0+0.5*(Dp(1)+Dp(2))
   if (miu(1)<=SEIS_ZERO .or. miu(2)<=SEIS_ZERO) then
      flag_water=.true.; nwater=nwater+1
   else
      mu0 = mu0+0.5*(1.0/miu(1)+1.0/miu(2))
   end if
   !mu0 = mu0+0.5*(1.0/miu(1)+1.0/miu(2))
   lam0=lam0+0.5*(1.0/lam(1)+1.0/lam(2))

end do !mi
end do !mj
end do !mk

    rho(i,j,k)=rho0/nsamp
    if (flag_water) then
       if (nwater>=nsamp/2) then
          mu(i,j,k)=0.0
       else
          mu(i,j,k)=(nsamp-nwater)/mu0
       end if
    else
       mu(i,j,k)=nsamp/mu0
    end if
    lambda(i,j,k)=nsamp/lam0

end do !i
end do !j
end do !k
end subroutine effmedia_eval

function eval_poly(c,d,x) result(f)
real(SP),dimension(:),intent(in) :: c,d
real(SP),intent(in) :: x
real(SP) :: f
integer :: nmax,n
nmax=size(d)
!f=pp(1)*x0**3+pp(2)*x0**2+pp(3)*x0+pp(4)
f=0.0
do n=1,nmax
   f=f+c(n) * (x**d(n))
end do
end function eval_poly

subroutine media_extend
  call extend_equal(rho)
  call extend_equal(mu)
  call extend_equal(lambda)
end subroutine media_extend
subroutine extend_equal(w)
   real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2),intent(inout) :: w
   integer i,j,k,n
   ! x1, x2
   do k=nz1,nz2
   do j=ny1,ny2
      do n=1,LenFD
         w(ni1-n,j,k)=w(ni1,j,k)
         w(ni2+n,j,k)=w(ni2,j,k)
      end do
   end do
   end do
   ! y1, y2
   do k=nz1,nz2
   do i=nx1,nx2
      do n=1,LenFD
         w(i,nj1-n,k)=w(i,nj1,k)
         w(i,nj2+n,k)=w(i,nj2,k)
      end do
   end do
   end do
   ! z1, z2
   do j=ny1,ny2
   do i=nx1,nx2
      do n=1,LenFD
         w(i,j,nk1-n)=w(i,j,nk1)
         w(i,j,nk2+n)=w(i,j,nk2)
      end do
   end do
   end do
end subroutine extend_equal

function dist_point2plane(x0,y0,z0,A,B,C) result (L)
real(DP),intent(in) :: x0,y0,z0
real(DP),dimension(SEIS_GEO),intent(in) :: A,B,C
real(DP) L
real(DP),dimension(SEIS_GEO) :: AB,AC,p
real(DP) c1,c2,c3,d

AB=B-A; AC=C-A
call times_product(AB,AC,p)
c1=p(1);c2=p(2);c3=p(3);
d=-dot_product(p,A)
L=abs( (c1*x0+c2*y0+c3*z0+d)/sqrt(dot_product(p,p)) )
end function dist_point2plane

subroutine indx_locate_1d(x0,x,i1,i2,x1)
real(SP),intent(in) :: x0
real(SP),dimension(:),intent(in) :: x
integer,intent(out) :: i1,i2
real(SP),intent(out) :: x1
integer :: i0
!logical,intent(in) :: backward
  !integer :: indx(1)
  !indx=minloc(x,x>=x0)
  !i2=max(indx(1),2)
  !i1=i2-1

i1=1; i2=size(x); x1=x0
if (x0<=x(1)) then
   x1=x(1); i1=1; i2=i1+1
elseif (x0>=x(i2)) then
   x1=x(i2); i1=i2-1
else
  do while (i2-i1>1)
     i0=(i2-i1)/2+i1
     if (x(i0)==x0) then
        i1=i0;i2=i0+1; exit
     elseif (x(i0)>x0) then
        i2=i0;
     else
        i1=i0
     end if
  end do
end if
end subroutine indx_locate_1d

subroutine media_create(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  integer,dimension(SEIS_GEO) ::subs,subc,subt
   filenm=media_fnm_get(n_i,n_j,n_k)
   call nfseis_grid3d_skel(filenm,nx,ny,nz,"media generated by seis3d_media" )
   subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
   call nfseis_grid3d_attput(filenm,        &
        subs,subc,subt,                   &
        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &
        nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &
        ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &
        ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,    &
        (/ngi1,ngi2,ngj1,ngj2,ngk1,ngk2/) )
#ifdef WITHQS
   call nfseis_attput(filenm,'QsF0',QsF0)
   call nfseis_attput(filenm,'QsINF',QsINF)
   call nfseis_grid3d_addvar(filenm,'Qs')
#endif
   call nfseis_grid3d_addvar(filenm,'rho')
   call nfseis_grid3d_addvar(filenm,'lambda')
   call nfseis_grid3d_addvar(filenm,'mu')

   call nfseis_varput(filenm,'rho',rho,subs,subc,subt)
   call nfseis_varput(filenm,'mu',mu,subs,subc,subt)
   call nfseis_varput(filenm,'lambda',lambda,subs,subc,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs,subs,subc,subt)
#endif
end subroutine media_create

#ifdef MediaMPI
subroutine media_exchange
integer :: ierr
integer,dimension(MPI_STATUS_SIZE) :: istatus
! exchange
call MPI_SENDRECV(rho(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),111,        &
                  rho(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),111,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),112, &
                  rho(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),112, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),121,        &
                  rho(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),121,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),122, &
                  rho(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),122, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),131,        &
                  rho(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),131,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),132, &
                  rho(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),132, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(mu (ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),211,        &
                  mu (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),211,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(mu (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),212, &
                  mu (ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),212, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(mu (nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),221,        &
                  mu (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),221,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(mu (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),222, &
                  mu (nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),222, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(mu (nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),231,        &
                  mu (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),231,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(mu (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),232, &
                  mu (nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),232, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(lambda(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),311,        &
                  lambda(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),311,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(lambda(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),312, &
                  lambda(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),312, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(lambda(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),321,        &
                  lambda(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),321,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(lambda(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),322, &
                  lambda(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),322, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(lambda(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),331,        &
                  lambda(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),331,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(lambda(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),332, &
                  lambda(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),332, &
                  SWMPI_COMM,istatus,ierr)
end subroutine media_exchange
#else
subroutine media_exchange
character (len=SEIS_STRLEN) :: filenm
integer,dimension(SEIS_GEO) :: subt
integer,dimension(SEIS_GEO) ::         &
     subs_x1,subc_x1, subs_x2,subc_x2, &
     subs_y1,subc_y1, subs_y2,subc_y2, &
     subs_z1,subc_z1, subs_z2,subc_z2
integer,dimension(LenFD) ::            &
     indx_x1,indx_x2,indx_y1,indx_y2,indx_z1,indx_z2
integer :: n_i,n_j,n_k
integer :: i,j,k

subt=(/ 1,1,1 /)
subs_x1=(/ ni2+1,ny1  ,nz1  /); subc_x1=(/ LenFD , ny    , nz    /)
subs_x2=(/ nx1  ,ny1  ,nz1  /); subc_x2=(/ LenFD , ny    , nz    /)
subs_y1=(/ nx1  ,nj2+1,nz1  /); subc_y1=(/ nx    , LenFD , nz    /)
subs_y2=(/ nx1  ,ny1  ,nz1  /); subc_y2=(/ nx    , LenFD , nz    /)
subs_z1=(/ nx1  ,ny1  ,nk2+1/); subc_z1=(/ nx    , ny    , LenFD /)
subs_z2=(/ nx1  ,ny1  ,nz1  /); subc_z2=(/ nx    , ny    , LenFD /)
indx_x1=(/ (i,i=ni1,ni1+LenFD-1) /)
indx_x2=(/ (i,i=ni2-LenFD+1,ni2) /)
indx_y1=(/ (j,j=nj1,nj1+LenFD-1) /)
indx_y2=(/ (j,j=nj2-LenFD+1,nj2) /)
indx_z1=(/ (k,k=nk1,nk1+LenFD-1) /)
indx_z2=(/ (k,k=nk2-LenFD+1,nk2) /)

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call media_import(n_i,n_j,n_k)

   ! to x1
if (n_i>0) then
   call swmpi_change_fnm(n_i-1,n_j,n_k)
   filenm=media_fnm_get(n_i-1,n_j,n_k)
   call nfseis_varput(filenm,'rho',rho(indx_x1,:,:),subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'mu' ,mu(indx_x1,:,:), subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'lambda',lambda(indx_x1,:,:),subs_x1,subc_x1,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(indx_x1,:,:),subs_x1,subc_x1,subt)
#endif
end if
! to x2
if (n_i<dims(1)-1) then
   call swmpi_change_fnm(n_i+1,n_j,n_k)
   filenm=media_fnm_get(n_i+1,n_j,n_k)
   call nfseis_varput(filenm,'rho',rho(indx_x2,:,:),subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'mu', mu(indx_x2,:,:),subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'lambda',lambda(indx_x2,:,:),subs_x2,subc_x2,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(indx_x2,:,:),subs_x2,subc_x2,subt)
#endif
end if
! to y1
if (n_j>0) then
   call swmpi_change_fnm(n_i,n_j-1,n_k)
   filenm=media_fnm_get(n_i,n_j-1,n_k)
   call nfseis_varput(filenm,'rho',rho(:,indx_y1,:),subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'mu', mu(:,indx_y1,:), subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'lambda',lambda(:,indx_y1,:),subs_y1,subc_y1,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,indx_y1,:),subs_y1,subc_y1,subt)
#endif
end if
! to y2
if (n_j<dims(2)-1) then
   call swmpi_change_fnm(n_i,n_j+1,n_k)
   filenm=media_fnm_get(n_i,n_j+1,n_k)
   call nfseis_varput(filenm,'rho',rho(:,indx_y2,:),subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'mu', mu(:,indx_y2,:), subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'lambda',lambda(:,indx_y2,:),subs_y2,subc_y2,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,indx_y2,:),subs_y2,subc_y2,subt)
#endif
end if
! to k1
if (n_k>0) then
   call swmpi_change_fnm(n_i,n_j,n_k-1)
   filenm=media_fnm_get(n_i,n_j,n_k-1)
   call nfseis_varput(filenm,'rho', rho(:,:,indx_z1),subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'mu', mu(:,:,indx_z1), subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'lambda',lambda(:,:,indx_z1),subs_z1,subc_z1,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,:,indx_z1),subs_z1,subc_z1,subt)
#endif
end if
! to k2
if (n_k<dims(3)-1) then
   call swmpi_change_fnm(n_i,n_j,n_k+1)
   filenm=media_fnm_get(n_i,n_j,n_k+1)
   call nfseis_varput(filenm,'rho', rho(:,:,indx_z2),subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'mu', mu(:,:,indx_z2),  subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'lambda',lambda(:,:,indx_z2),subs_z2,subc_z2,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,:,indx_z2),subs_z2,subc_z2,subt)
#endif
end if
end do
end do
end do
end subroutine media_exchange
#endif

subroutine media_stept_check(n_i,n_j,n_k)
integer,intent(in) :: n_i,n_j,n_k
real(SP) :: Vp,dtLe,dtlocal
integer :: i,j,k

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      !gx(i,j,k)=z(k)*cos(real(y(j),DP))*cos(real(x(i),DP))
      !gy(i,j,k)=z(k)*cos(real(y(j),DP))*sin(real(x(i),DP))
      !gz(i,j,k)=z(k)*sin(real(y(j),DP))
      gx(i,j,k)=z(k)*cos(real(x(i)-PI/2,DP))*cos(real(y(j),DP))
      gy(i,j,k)=z(k)*cos(real(x(i)-PI/2,DP))*sin(real(y(j),DP))
      gz(i,j,k)=z(k)*sin(real(x(i)-PI/2,DP))
   end do
   end do
   end do
   
   do k=nk1,nk2
   do j=nj1,nj2
   do i=ni1,ni2
      Vp=sqrt( (lambda(i,j,k)+2.0*mu(i,j,k))/rho(i,j,k) )
      dtLe=min( dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) )
      dtlocal=1.3/Vp * dtLe 
      if (dtlocal<dtmax) then
         dtmax=dtlocal
         dtnode=(/ n_i,n_j,n_k /)
         dtindx=(/ i,j,k /)
         dtmaxVp=Vp
         dtmaxL=dtLe
      end if
   end do
   end do
   end do
   write(*,"(a,3i5,a,f,a,3i5,a,2f)") "   dtmax in thread", n_i,n_j,n_k, " is", dtmax,  &
       " located on", dtindx," with Vp and dL", dtmaxVp,dtmaxL
end subroutine media_stept_check

subroutine alloc_local
allocate(gx(nx1:nx2,ny1:ny2,nz1:nz2)); gx=0.0_SP
allocate(gy(nx1:nx2,ny1:ny2,nz1:nz2)); gy=0.0_SP
allocate(gz(nx1:nx2,ny1:ny2,nz1:nz2)); gz=0.0_SP
end subroutine alloc_local

subroutine dealloc_media
if (L3D%yes) then
   deallocate(L3D%x)
   deallocate(L3D%y)
   deallocate(L3D%h)
   deallocate(L3D%Vp)
   deallocate(L3D%Vs)
   deallocate(L3D%Dp)
end if
end subroutine dealloc_media

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef MediaMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef MediaMPI
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

end program seis3d_media

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:

