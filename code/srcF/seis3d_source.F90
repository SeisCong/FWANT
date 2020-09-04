program seis3d_source

! This program inits and distributes the seismic source
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:54:24 -0500 (Fri, 16 Jan 2009) $
! $Revision: 68 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod
use string_mod
use math_mod
use nfseis_mod
use para_mod
use grid_mod
use media_mod
use src_mod
use mpi_mod
implicit none

real(SP) :: src_hyper_height
logical,allocatable :: force_flag(:),moment_flag(:)
real(SP),allocatable :: force_axis(:,:), moment_axis(:,:)

character (len=SEIS_STRLEN) :: filenm
integer :: p(1)
real(SP) :: xmin,xmax,ymin,ymax,zmin,zmax
real(SP) :: x0,y0,z0
integer :: i,j,k,n_i,n_j,n_k
integer :: n,m,npt,nfrc,nmom,gi,gj,gk,si,sj,sk
logical :: iflag
integer :: ncid,indxid,axisid,siftid,t0id,stftid,stffid
integer :: fxid,fyid,fzid,mxxid,myyid,mzzid,mxyid,mxzid,myzid

!----------------------------------------------------------------------!

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)
call media_fnm_init(fnm_conf)
call src_fnm_init(fnm_conf)

! load grid
call grid_alloc

! source
call read_src_para(fnm_src_conf)

! find src index
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_coord_import(n_i,n_j,n_k)

   xmin=minval(x); xmax=maxval(x)
   ymin=minval(y); ymax=maxval(y)
   zmin=minval(z); zmax=maxval(z)
   print *, 'xmin=',xmin/pi*180,'xmax=',xmax/pi*180
   print *, 'ymin=',ymin/pi*180,'ymax=',ymax/pi*180
   print *, 'zmin=',zmin,'zmax=',zmax

!force
do npt=1,num_force
   if (force_flag(npt)) cycle

   x0=force_axis(1,npt);y0=force_axis(2,npt);z0=force_axis(3,npt)
   if (x0<xmin .or. x0>xmax .or. y0<ymin .or. y0>ymax) cycle

   p=minloc(abs(x-x0)); i=loct_i(p(1))
   p=minloc(abs(y-y0)); j=loct_j(p(1))

   if( z0 >=zmin .and. z0 <=zmax ) then
      p=minloc(abs(z-z0)); k=loct_k(p(1))
   elseif (z0>src_hyper_height .and. n_k==dims(3)-1) then
      k=nk2; z0=z(k)
   else
      cycle
   end if

   if ( i<ni1 .or. i>ni2 .or. j<nj1 .or. j>nj2 .or. k<nk1 .or. k>nk2 ) cycle

   write(*,"(i10,a,i10,2(i2),a,3(i2))") npt, ' in', n_i,n_j,n_k, ' of ',dims
   force_flag(npt)=.true.
   gi=swmpi_globi(i,n_i); gj=swmpi_globj(j,n_j); gk=swmpi_globk(k,n_k)
   force_indx(:,npt)=(/ gi,gj,gk /)
   force_axis(3,npt)=z0

   ! calculate shift
   force_shift(:,npt)=0.0_SP
   if (x0>x(i)) then
      force_shift(1,npt)=(x0-x(i))/(x(i+1)-x(i))
   else
      force_shift(1,npt)=(x0-x(i))/(x(i)-x(i-1))
   end if
   if (y0>y(j)) then
      force_shift(2,npt)=(y0-y(j))/(y(j+1)-y(j))
   else
      force_shift(2,npt)=(y0-y(j))/(y(j)-y(j-1))
   end if
   if (z0>z(k)) then
      force_shift(3,npt)=(z0-z(k))/(z(k+1)-z(k))
   else
      force_shift(3,npt)=(z0-z(k))/(z(k)-z(k-1))
   end if
end do !n

!moment
do npt=1,num_moment
   if (moment_flag(npt)) cycle

   x0=moment_axis(1,npt);y0=moment_axis(2,npt);z0=moment_axis(3,npt)
   if (x0<xmin .or. x0>xmax .or. y0<ymin .or. y0>ymax) cycle

   p=minloc(abs(x-x0)); i=loct_i(p(1))
   p=minloc(abs(y-y0)); j=loct_j(p(1))

   if( z0 >=zmin .and. z0 <=zmax ) then
      p=minloc(abs(z-z0)); k=loct_k(p(1))
   elseif (z0>src_hyper_height .and. n_k==dims(3)-1) then
      k=nk2; z0=z(k)
   else
      cycle
   end if

   if ( i<ni1 .or. i>ni2 .or. j<nj1 .or. j>nj2 .or. k<nk1 .or. k>nk2 ) cycle

   write(*,"(i10,a,i10,2(i2),a,3(i2))") npt, ' in', n_i,n_j,n_k, ' of ',dims
   moment_flag(npt)=.true.
   gi=swmpi_globi(i,n_i); gj=swmpi_globj(j,n_j); gk=swmpi_globk(k,n_k)
   moment_indx(:,npt)=(/ gi,gj,gk /)
   moment_axis(3,npt)=z0

   ! calculate shift
   moment_shift(:,npt)=0.0_SP
   if (x0>x(i)) then
      moment_shift(1,npt)=(x0-x(i))/(x(i+1)-x(i))
   else
      moment_shift(1,npt)=(x0-x(i))/(x(i)-x(i-1))
   end if
   if (y0>y(j)) then
      moment_shift(2,npt)=(y0-y(j))/(y(j+1)-y(j))
   else
      moment_shift(2,npt)=(y0-y(j))/(y(j)-y(j-1))
   end if
   if (z0>z(k)) then
      moment_shift(3,npt)=(z0-z(k))/(z(k+1)-z(k))
   else
      moment_shift(3,npt)=(z0-z(k))/(z(k)-z(k-1))
   end if
end do !n

end do !n_k
end do !n_j
end do !n_i

! check
iflag=.false.
do n=1,num_force
   if ( .not. force_flag(n)) then
      iflag=.true.
      print *, 'n=',n,'loct=',force_axis(:,n),'indx=',force_indx(:,n)
   end if
end do
do n=1,num_moment
   if ( .not. moment_flag(n)) then
      iflag=.true.
      print *, 'n=',n,'loct=',moment_axis(:,n),'indx=',moment_indx(:,n)
   end if
end do
if (iflag) then
   print *, 'there are some source points out of the computational domain'
   stop 1
end if

! distribute to thread
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   nfrc=0
   do n=1,num_force
      gi=force_indx(1,n);gj=force_indx(2,n);gk=force_indx(3,n)
      if (      gi>=ngx1 .and. gi<=ngx2                          &
          .and. gj>=ngy1 .and. gj<=ngy2                          &
          .and. gk>=ngz1 .and. gk<=ngz2 ) nfrc=nfrc+1
   end do
   nmom=0
   do n=1,num_moment
      gi=moment_indx(1,n);gj=moment_indx(2,n);gk=moment_indx(3,n)
      if (      gi>=ngx1 .and. gi<=ngx2                          &
          .and. gj>=ngy1 .and. gj<=ngy2                          &
          .and. gk>=ngz1 .and. gk<=ngz2 ) nmom=nmom+1
   end do

   filenm=src_fnm_get(n_i,n_j,n_k)
   call srcnode_skel(filenm,nfrc,ntwin_force,nmom,ntwin_moment)
   call nfseis_open(filenm,ncid)
  
   !force
   if (nfrc>0) then
      call nfseis_inq_varid(ncid,'force_indx',indxid)
      call nfseis_inq_varid(ncid,'force_axis',axisid)
      call nfseis_inq_varid(ncid,'force_shift',siftid)
      call nfseis_inq_varid(ncid,'force_start_time',t0id)
      call nfseis_inq_varid(ncid,'force_stf_time',stftid)
      call nfseis_inq_varid(ncid,'force_stf_freq',stffid)
      call nfseis_inq_varid(ncid,'Fx',fxid)
      call nfseis_inq_varid(ncid,'Fy',fyid)
      call nfseis_inq_varid(ncid,'Fz',fzid)
      call nfseis_put(ncid,stftid,frcstf_time,(/1/),(/ntwin_force/),(/1/))
      call nfseis_put(ncid,stffid,frcstf_freq,(/1/),(/ntwin_force/),(/1/))
      m=0
      do n=1,num_force
         gi=force_indx(1,n);gj=force_indx(2,n);gk=force_indx(3,n)
         if (      gi>=ngx1 .and. gi<=ngx2                          &
             .and. gj>=ngy1 .and. gj<=ngy2                          &
             .and. gk>=ngz1 .and. gk<=ngz2 ) then
            m=m+1
            call nfseis_put(ncid,fxid,ForceX(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
            call nfseis_put(ncid,fyid,ForceY(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
            call nfseis_put(ncid,fzid,ForceZ(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
            call nfseis_put(ncid,axisid,force_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
            call nfseis_put(ncid,siftid,force_shift(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
            call nfseis_put(ncid,t0id,force_t0(n),(/m/),(/1/),(/1/))
            call nfseis_put(ncid,indxid,         &
                 (/ swmpi_locli(gi,n_i),         &
                    swmpi_loclj(gj,n_j),         &
                    swmpi_loclk(gk,n_k) /),      &
                 (/1,m/),(/SEIS_GEO,1/),(/1,1/))
         end if
      end do
   end if
   !moment
   if (nmom>0) then
      call nfseis_inq_varid(ncid,'moment_indx',indxid)
      call nfseis_inq_varid(ncid,'moment_axis',axisid)
      call nfseis_inq_varid(ncid,'moment_shift',siftid)
      call nfseis_inq_varid(ncid,'moment_start_time',t0id)
      call nfseis_inq_varid(ncid,'moment_stf_time',stftid)
      call nfseis_inq_varid(ncid,'moment_stf_freq',stffid)
      call nfseis_inq_varid(ncid,'Mxx',mxxid)
      call nfseis_inq_varid(ncid,'Myy',myyid)
      call nfseis_inq_varid(ncid,'Mzz',mzzid)
      call nfseis_inq_varid(ncid,'Mxy',mxyid)
      call nfseis_inq_varid(ncid,'Mxz',mxzid)
      call nfseis_inq_varid(ncid,'Myz',myzid)
      call nfseis_put(ncid,stftid,momstf_time,(/1/),(/ntwin_moment/),(/1/))
      call nfseis_put(ncid,stffid,momstf_freq,(/1/),(/ntwin_moment/),(/1/))
      m=0
      do n=1,num_moment
         gi=moment_indx(1,n);gj=moment_indx(2,n);gk=moment_indx(3,n)
         if (      gi>=ngx1 .and. gi<=ngx2                          &
             .and. gj>=ngy1 .and. gj<=ngy2                          &
             .and. gk>=ngz1 .and. gk<=ngz2 ) then
            m=m+1
            call nfseis_put(ncid,mxxid,MomTxx(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,myyid,MomTyy(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,mzzid,MomTzz(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,mxyid,MomTxy(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,mxzid,MomTxz(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,myzid,MomTyz(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,axisid,moment_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
            call nfseis_put(ncid,siftid,moment_shift(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
            call nfseis_put(ncid,t0id,moment_t0(n),(/m/),(/1/),(/1/))
            call nfseis_put(ncid,indxid,         &
                 (/ swmpi_locli(gi,n_i),         &
                    swmpi_loclj(gj,n_j),         &
                    swmpi_loclk(gk,n_k) /),      &
                 (/1,m/),(/SEIS_GEO,1/),(/1,1/))
         end if
      end do
   end if
   call nfseis_close(ncid)
end do
end do
end do

call src_destroy
call dealloc_local

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

!*************************************************************************
!*                 PART-I  src alloc and dealloc                         *
!*************************************************************************

subroutine alloc_force(npt)
  integer,intent(in) :: npt
  allocate(force_flag(npt)); force_flag=.false.
  allocate(force_axis(SEIS_GEO,npt)); force_axis=0.0_SP
end subroutine alloc_force
subroutine alloc_moment(npt)
  integer,intent(in) :: npt
  allocate(moment_flag(npt)); moment_flag=.false.
  allocate(moment_axis(SEIS_GEO,npt)); moment_axis=0.0_SP
end subroutine alloc_moment

subroutine dealloc_local
  if (allocated(force_flag)) deallocate(force_flag)
  if (allocated(moment_flag)) deallocate(moment_flag)
end subroutine dealloc_local

!*************************************************************************
!*                    PART-II  src init and distrib                      *
!*************************************************************************

subroutine read_src_para(fnm_conf)
character (len=*),intent(in) :: fnm_conf
character (len=SEIS_STRLEN) :: mommech,stf_type,str
real(DP) :: strike,dip,rake
real(SP) :: f0,m0
integer fid,n,m

fid=1001
open(fid,file=trim(fnm_conf),status="old")

call string_conf(fid,1,'src_hyper_height',2,src_hyper_height)

!force
call string_conf(fid,1,"number_of_force_source",2,num_force)
if (num_force>=1) then
   call string_conf(fid,1,"force_stf_window",2,ntwin_force)
   call string_conf(fid,1,"force_stf_type",2,stf_type)
   frcstf_id=stf_name2id(trim(stf_type))
   call src_alloc_force(num_force,ntwin_force)
   call alloc_force(num_force)
   do m=1,ntwin_force
      call string_conf(fid,1,"force_stf_timefactor",m+1,frcstf_time(m))
      call string_conf(fid,1,"force_stf_freqfactor",m+1,frcstf_freq(m))
   end do
   call string_conf(fid,1,"<anchor_force>",1,str)
   do n=1,num_force
   do m=1,ntwin_force
      read(fid,*) force_axis(:,n),force_t0(n), &
          f0,ForceX(m,n),ForceY(m,n),ForceZ(m,n)
      ForceX(m,n)=ForceX(m,n)*f0
      ForceY(m,n)=ForceY(m,n)*f0
      ForceZ(m,n)=ForceZ(m,n)*f0
   end do
      force_axis(1:2,n)=force_axis(1:2,n)*PI/180.0_SP
   end do
end if
!moment
call string_conf(fid,1,"number_of_moment_source",2,num_moment)
if (num_moment>=1) then
   call string_conf(fid,1,"moment_stf_window",2,ntwin_moment)
   call string_conf(fid,1,"moment_stf_type",2,stf_type)
   momstf_id=stf_name2id(trim(stf_type))
   call src_alloc_moment(num_moment,ntwin_moment)
   call alloc_moment(num_moment)
   do m=1,ntwin_moment
      call string_conf(fid,1,"moment_stf_timefactor",m+1,momstf_time(m))
      call string_conf(fid,1,"moment_stf_freqfactor",m+1,momstf_freq(m))
   end do
   call string_conf(fid,1,"moment_mech_input",2,mommech)
   call string_conf(fid,1,"<anchor_moment>",1,str)
   if (trim(mommech)=='moment') then
      do n=1,num_moment
      do m=1,ntwin_moment
         read(fid,*) moment_axis(:,n),moment_t0(n),m0, &
           MomTxx(m,n),MomTyy(m,n),MomTzz(m,n),MomTxy(m,n),MomTxz(m,n),MomTyz(m,n)
         MomTxx(m,n)=m0*MomTxx(m,n);MomTyy(m,n)=m0*MomTyy(m,n);MomTzz(m,n)=m0*MomTzz(m,n)
         MomTxy(m,n)=m0*MomTxy(m,n);MomTxz(m,n)=m0*MomTxz(m,n);MomTyz(m,n)=m0*MomTyz(m,n)
      end do
         moment_axis(1:2,n)=moment_axis(1:2,n)*PI/180.0_SP
      end do
   else
      do n=1,num_moment
      do m=1,ntwin_moment
         read(fid,*) moment_axis(:,n),moment_t0(n),m0,strike,dip,rake
         call angle2moment(strike,dip,rake, &
              MomTxx(m,n),MomTyy(m,n),MomTzz(m,n),MomTxy(m,n),MomTxz(m,n),MomTyz(m,n))
         MomTxx(m,n)=m0*MomTxx(m,n);MomTyy(m,n)=m0*MomTyy(m,n);MomTzz(m,n)=m0*MomTzz(m,n)
         MomTxy(m,n)=m0*MomTxy(m,n);MomTxz(m,n)=m0*MomTxz(m,n);MomTyz(m,n)=m0*MomTyz(m,n)
      end do
         moment_axis(1:2,n)=moment_axis(1:2,n)*PI/180.0_SP
      end do
   end if
end if
close(fid)
end subroutine read_src_para

subroutine angle2moment(strike,dip,rake,Mxx,Myy,Mzz,Mxy,Mxz,Myz)
real(DP),intent(in) :: strike,dip,rake
real(SP),intent(out) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

real(DP) :: strike_pi,dip_pi,rake_pi
real(SP) :: M11,M22,M33,M12,M13,M23

dip_pi=dip/180.0_DP*PI
strike_pi=strike/180.0_DP*PI
rake_pi=rake/180.0_DP*PI
!in Aki and Richard's
M11=-(sin(dip_pi)*cos(rake_pi)*sin(2.0_DP*strike_pi)   &
         +sin(2.0_DP*dip_pi)*sin(rake_pi)*sin(strike_pi)**2)
M22=sin(dip_pi)*cos(rake_pi)*sin(2.0_DP*strike_pi)     &
         -sin(2.0_DP*dip_pi)*sin(rake_pi)*cos(strike_pi)**2
!Mzz=sin(2.0*dip_pi)*sin(rake_pi)
M33=-(M11+M22)
M12=sin(dip_pi)*cos(rake_pi)*cos(2.0_DP*strike_pi)     &
         +0.5_DP*sin(2.0_DP*dip_pi)*sin(rake_pi)*sin(2.0_DP*strike_pi)
M13=-(cos(dip_pi)*cos(rake_pi)*cos(strike_pi)       &
         +cos(2.0_DP*dip_pi)*sin(rake_pi)*sin(strike_pi))
M23=-(cos(dip_pi)*cos(rake_pi)*sin(strike_pi)       &
         -cos(2.0_DP*dip_pi)*sin(rake_pi)*cos(strike_pi))
!to spherical
Mxx= M11; Myy=M22; Mzz=M33
Mxy= -M12; Mxz=M13; Myz=-M23
end subroutine angle2moment

!----------------------------------------------------------------

subroutine srcnode_skel(filenm,nfrc,ntwfrc,nmom,ntwmom)
character (len=*),intent(in) :: filenm
integer,intent(in) :: nfrc,ntwfrc,nmom,ntwmom

integer :: ncid,ierr,oldMode
integer :: geoid,vid,nmomid,ntwmomid,nfrcid,ntwfrcid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'srcnode_init:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in srcnode_skel')
! -- define dim
ierr=nf90_def_dim(ncid,'geo_dimension',SEIS_GEO,geoid)
     call nfseis_except(ierr,'geodim dim in srcnode_skel')
ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_force",nfrc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_moment",nmom)
!force
if (nfrc>0) then
ierr=nf90_def_dim(ncid,'force_point',nfrc,nfrcid)
     call nfseis_except(ierr,'force_point dim in srcnode_skel')
ierr=nf90_def_dim(ncid,'force_time_window' ,ntwfrc,ntwfrcid)
     call nfseis_except(ierr,'force_time_window dim in srcnode_skel')
! -- define variable
ierr=nf90_def_var(ncid,'Fx', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'Fy', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'Fz', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_indx', nf90_int,(/ geoid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_axis', SEISNC_DATATYPE,(/ geoid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_shift', SEISNC_DATATYPE,(/ geoid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_start_time', SEISNC_DATATYPE,(/ nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_stf_time', SEISNC_DATATYPE,(/ ntwfrcid /),vid )
ierr=nf90_def_var(ncid,'force_stf_freq', SEISNC_DATATYPE,(/ ntwfrcid /),vid )
ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_id",frcstf_id)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_type", &
      trim(stf_id2name(frcstf_id)) )
end if
!moment
if (nmom>0) then
ierr=nf90_def_dim(ncid,'moment_point',nmom,nmomid)
     call nfseis_except(ierr,'moment_point dim in srcnode_skel')
ierr=nf90_def_dim(ncid,'moment_time_window' ,ntwmom,ntwmomid)
     call nfseis_except(ierr,'moment_time_window dim in srcnode_skel')
! -- define variable
ierr=nf90_def_var(ncid,'Mxx',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Myy',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Mzz',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Mxy',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Mxz',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Myz',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'moment_indx', nf90_int,(/ geoid, nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_axis', SEISNC_DATATYPE,(/ geoid, nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_shift', SEISNC_DATATYPE,(/ geoid, nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_start_time', SEISNC_DATATYPE,(/ nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_stf_time', SEISNC_DATATYPE,(/ ntwmomid /),vid )
ierr=nf90_def_var(ncid,'moment_stf_freq', SEISNC_DATATYPE,(/ ntwmomid /),vid )
ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_stf_id",momstf_id)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_stf_type", &
      trim(stf_id2name(momstf_id)) )
end if
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in srcnode_skel')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in srcnode_skel')
end subroutine srcnode_skel

end program seis3d_source

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
