program SI_ker2info

! This program is used to calculate kernel distribution.
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2008 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 23:30:18 -0500 (Wed, 14 Jan 2009) $
! $Revision: 53 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

#define VERBOSE

use constants_mod
use string_mod
use para_mod
use io_mod
use media_mod
use nfseis_mod
use src_mod
use mpi_mod
#ifdef KernelInfoMPI
use mpi
#endif

implicit none
integer,parameter :: MAXKSTNM=8,MAXKEVNM=20

character (len=8) :: str_d0,str_d1,str_d2
character (len=10):: str_t0,str_t1,str_t2

character (len=300) :: rcdstr
character (len=SEIS_STRLEN) ::           &
    fnm_main_conf,fnm_ker_conf,          &
    evt_list,sgt_list,pnm_info,          &
    fnm_synx,fnm_syny,fnm_synz,          &
    pnm_ker,fnm_ker
character (len=SEIS_STRLEN) ::           &
    KSTNM,KEVNM
character (len=SEIS_STRLEN) :: filenm
character (len=15) :: idstr

integer,dimension(SEIS_GEO) :: blksiz
integer :: kid

integer :: ncmp,npk,nfq
integer :: picknum(6)
integer :: pickfreq
real(SP),dimension(3,3) :: &
  rot_mat
character (len=2) :: filt_nm
real(SP) :: t1,t2

real(SP),dimension(:,:,:),allocatable :: &
    Kap,Kaq,Kbp,Kbq, &
    Sap,Saq,Sbp,Sbq
integer,dimension(:,:,:),allocatable :: &
    Nap,Naq,Nbp,Nbq, &
    Map,Maq,Mbp,Mbq
real :: kap0,kaq0,kbp0,kbq0

integer :: n_i,n_j,n_k
integer :: i,j,k,m,n,ierr
integer,dimension(SEIS_GEO) :: &
   bsubs,bsubc,bsubt,          &
    subs, subc, subt

integer :: staid,evtid

integer :: ncid
integer,dimension(4) :: &
  Sid,Nid,Mid,stnmid,evnmid

!----------------------------------------------------------------------

#ifdef KernelInfoMPI
call MPI_INIT(ierr)
#endif

call get_conf_name(fnm_conf)

! read kernel conf
fnm_ker_conf='TomoKernel.conf'
call init_kernel(fnm_ker_conf)
fnm_conf=fnm_main_conf

call para_init(fnm_conf)
call swmpi_init(fnm_conf)

#ifdef KernelInfoMPI
  call swmpi_cart_creat
  call swmpi_reinit_para
  call swmpi_datatype
  if (masternode) then
     read(*,"(a)") sgt_list
     read(*,"(a)") pnm_info
     read(*,*) kap0,kaq0,kbp0,kbq0
     write(*,*) sgt_list
     write(*,*) pnm_info
     write(*,*) kap0,kaq0,kbp0,kbq0
  end if
  call MPI_BCAST(sgt_list,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(pnm_info,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kap0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kaq0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kbp0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kbq0,1,MPI_REAL,0,SWMPI_COMM,ierr)
#else
  call swmpi_set_gindx(0,0,0)
  read(*,*) sgt_list
  read(*,*) pnm_info
  read(*,*) kap0,kaq0,kbp0,kbq0
#endif

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

!----------------------------------------------------------------------

#ifndef KernelInfoMPI
  print *, 'input mpi id:'
  read *, n_i,n_j,n_k
  call swmpi_change_fnm(n_i,n_j,n_k)
  call swmpi_set_gindx(n_i,n_j,n_k)
  thisid=(/ n_i,n_j,n_k /)
#else
  n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#endif

call date_and_time(date=str_d0,time=str_t0)
write(idstr,"(a8,3(i2.2))") ' thisid=',n_i,n_j,n_k
#ifdef VERBOSE
  write(*,*) idstr, ' begins from ',str_d0,  &
                  ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
#endif

call io_snap_locate(n_i,n_j,n_k)

!----------------------------------------------------------------------
kernode_if : if (snap_ishere(kid) .and. sgt_out(kid)) then

blksiz=snap_subc(:,kid)  
call alloc_ker_var(blksiz(1),blksiz(2),blksiz(3))

staid=10001; evtid=10002
! loop all stations
open(staid,file=trim(sgt_list),status='old',iostat=ierr)
if (ierr>0) call error_except("sgt_list open err:"//trim(sgt_list))

fnm_ker=get_fnm_snapnode_n('./','kerinfo_',kid,0,n_i,n_j,n_k)
call kernelinfo_skel(trim(pnm_info)//"/"//trim(fnm_ker), &
     blksiz(1),blksiz(2),blksiz(3))

!-----------------------------------------------------------------------------
sta_loop: do

  !read(staid,*,iostat=ierr) KSTNM,evt_list
  read(staid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit sta_loop
  !if (KSTNM(1:1)=='#') cycle
  rcdstr=adjustl(rcdstr)
  if (rcdstr(1:1)=='#') cycle

  n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
  KSTNM=rcdstr(1:n)
  rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
  if (n==0) then
     evt_list=trim(rcdstr)
  else
     evt_list=rcdstr(1:n)
  end if

#ifdef VERBOSE
  call date_and_time(date=str_d1,time=str_t1)
  write(*,*) idstr//trim(KSTNM)//" "//trim(evt_list), &
                    ': begins from ',str_d1,  &
                  ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
#endif

   bsubs=(/ 1,1,1 /)
   bsubc=blksiz
   bsubt=(/ 1,1,1 /)
   subc=bsubc; subt=snap_subt(:,kid)*bsubt
   subs=snap_subs(:,kid)+(bsubs-1)*subt

! loop all event
!-----------------------------------------------------------------------------
open(evtid,file=trim(evt_list),status='old',iostat=ierr)
if (ierr>0) call error_except("evt_list open err:"//trim(evt_list))

evt_loop: do
  read(evtid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  if (len_trim(rcdstr)==300) then
     print *, len_trim(rcdstr)
     print *, rcdstr
     call error_except("beyond maximum character 300 limit")
  end if
  rcdstr=adjustl(rcdstr)
  if (rcdstr(1:1)=='#') cycle

#ifdef VERBOSE
   write(*,*) idstr//trim(rcdstr)
#endif

  n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
  KEVNM=rcdstr(1:n)
  rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  read(rcdstr,*) picknum, rot_mat

  read(evtid,"(a132)") fnm_synx
  read(evtid,"(a132)") fnm_syny
  read(evtid,"(a132)") fnm_synz

  if (all(picknum<=0)) cycle

  call reset_var_zero

! loop component
do ncmp=1,6
   if (picknum(ncmp)<=0) cycle
do npk=1,picknum(ncmp)
   do
   read(evtid,'(a300)',iostat=ierr) rcdstr
   if (ierr<0) call error_except("read event pick error")
   if (len_trim(rcdstr)==300) then
      print *, len_trim(rcdstr)
      print *, rcdstr
      call error_except("beyond maximum character 300 limit")
   end if
   rcdstr=adjustl(rcdstr)
   if (rcdstr(1:1)=='#') cycle
   exit
   end do

   !read(evtid,*) pickfreq(npk),filt_nm,t1,t2,pnm_ker(npk)
   ! pickfreq
#ifdef VERBOSE
   write(*,*) idstr//trim(rcdstr)
#endif
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) pickfreq
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   ! pickfilt
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   filt_nm=rcdstr(1:n)
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   ! t1 t2
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) t1
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) t2
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   ! pnm_ker
   pnm_ker=trim(rcdstr)

   ! kernel nc file
   fnm_ker=get_fnm_snapnode_n('./','kernel_',kid,0,n_i,n_j,n_k)
   fnm_ker=trim(pnm_ker)//"/"//trim(fnm_ker)
   call nfseis_varget(trim(fnm_ker),'phase_Vp',Kap,bsubs,bsubc,bsubt)
   call nfseis_varget(trim(fnm_ker),'amplitude_Vp',Kaq,bsubs,bsubc,bsubt)
   call nfseis_varget(trim(fnm_ker),'phase_Vs',Kbp,bsubs,bsubc,bsubt)
   call nfseis_varget(trim(fnm_ker),'amplitude_Vs',Kbq,bsubs,bsubc,bsubt)

   do k=1,blksiz(3)
   do j=1,blksiz(2)
   do i=1,blksiz(1)
      Sap(i,j,k)=Sap(i,j,k)+abs(Kap(i,j,k))
      Saq(i,j,k)=Saq(i,j,k)+abs(Kaq(i,j,k))
      Sbp(i,j,k)=Sbp(i,j,k)+abs(Kbp(i,j,k))
      Sbq(i,j,k)=Sbq(i,j,k)+abs(Kbq(i,j,k))
      if (abs(Kap(i,j,k))>Kap0) Map(i,j,k)=Map(i,j,k)+1
      if (abs(Kaq(i,j,k))>Kaq0) Maq(i,j,k)=Maq(i,j,k)+1
      if (abs(Kbp(i,j,k))>Kbp0) Mbp(i,j,k)=Mbp(i,j,k)+1
      if (abs(Kbq(i,j,k))>Kbq0) Mbq(i,j,k)=Mbq(i,j,k)+1
   end do
   end do
   end do

end do ! npk
end do ! ncmp

   do k=1,blksiz(3)
   do j=1,blksiz(2)
   do i=1,blksiz(1)
      if (Map(i,j,k)>=1) then
         Nap(i,j,k)=Nap(i,j,k)+1
         ierr=nf90_put_var(ncid,stnmid(1),trim(KSTNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KSTNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,evnmid(1),trim(KEVNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KEVNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,Mid(1),Map(i:i,j,k),                               &
              (/i,j,k,Nap(i,j,k)/),(/1,1,1,1/),(/1,1,1,1/))
      end if
      if (Maq(i,j,k)>=1) then
         Naq(i,j,k)=Naq(i,j,k)+1
         ierr=nf90_put_var(ncid,stnmid(2),trim(KSTNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KSTNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,evnmid(2),trim(KEVNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KEVNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,Mid(2),Maq(i:i,j,k),                               &
              (/i,j,k,Nap(i,j,k)/),(/1,1,1,1/),(/1,1,1,1/))
      end if
      if (Mbp(i,j,k)>=1) then
         Nbp(i,j,k)=Nbp(i,j,k)+1
         ierr=nf90_put_var(ncid,stnmid(3),trim(KSTNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KSTNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,evnmid(3),trim(KEVNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KEVNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,Mid(3),Mbp(i:i,j,k),                               &
              (/i,j,k,Nap(i,j,k)/),(/1,1,1,1/),(/1,1,1,1/))
      end if
      if (Mbq(i,j,k)>=1) then
         Nbq(i,j,k)=Nbq(i,j,k)+1
         ierr=nf90_put_var(ncid,stnmid(4),trim(KSTNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KSTNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,evnmid(4),trim(KEVNM),                          &
              (/1,i,j,k,Nap(i,j,k)/),(/len_trim(KEVNM),1,1,1,1/),(/1,1,1,1,1/))
         ierr=nf90_put_var(ncid,Mid(4),Mbq(i:i,j,k),                               &
              (/i,j,k,Nap(i,j,k)/),(/1,1,1,1/),(/1,1,1,1/))
      end if
   end do
   end do
   end do

end do evt_loop
close(evtid)

end do sta_loop
close(staid)

   ! put
   call nfseis_put(ncid,Nid(1),Nap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,Nid(2),Naq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,Nid(3),Nbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,Nid(4),Nbq,bsubs,bsubc,bsubt)

   call nfseis_put(ncid,Sid(1),Sap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,Sid(2),Saq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,Sid(3),Sbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,Sid(4),Sbq,bsubs,bsubc,bsubt)

   call nfseis_close(ncid)

!-----------------------------------------------------------------------------
end if kernode_if

#ifdef VERBOSE
  call date_and_time(date=str_d2,time=str_t2)
  write(*,*) idstr, ' finish at ',str_d2,  &
                     ',',str_t2(1:2),'/',str_t2(3:4),'/',str_t2(5:10)
#endif

!-----------------------------------------------------------------------------
call dealloc_all

#ifdef KernelInfoMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine alloc_ker_var(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  allocate(Kap(ki,kj,kk)); Kap=0.0
  allocate(Kaq(ki,kj,kk)); Kaq=0.0
  allocate(Kbp(ki,kj,kk)); Kbp=0.0
  allocate(Kbq(ki,kj,kk)); Kbq=0.0
  allocate(Sap(ki,kj,kk)); Sap=0.0
  allocate(Saq(ki,kj,kk)); Saq=0.0
  allocate(Sbp(ki,kj,kk)); Sbp=0.0
  allocate(Sbq(ki,kj,kk)); Sbq=0.0
  allocate(Nap(ki,kj,kk)); Nap=0
  allocate(Naq(ki,kj,kk)); Naq=0
  allocate(Nbp(ki,kj,kk)); Nbp=0
  allocate(Nbq(ki,kj,kk)); Nbq=0
  allocate(Map(ki,kj,kk)); Map=0
  allocate(Maq(ki,kj,kk)); Maq=0
  allocate(Mbp(ki,kj,kk)); Mbp=0
  allocate(Mbq(ki,kj,kk)); Mbq=0
end subroutine alloc_ker_var

subroutine dealloc_all
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
  if (allocated(Sap)) deallocate(Sap)
  if (allocated(Saq)) deallocate(Saq)
  if (allocated(Sbp)) deallocate(Sbp)
  if (allocated(Sbq)) deallocate(Sbq)
  if (allocated(Nap)) deallocate(Nap)
  if (allocated(Naq)) deallocate(Naq)
  if (allocated(Nbp)) deallocate(Nbp)
  if (allocated(Nbq)) deallocate(Nbq)
  if (allocated(Map)) deallocate(Map)
  if (allocated(Maq)) deallocate(Maq)
  if (allocated(Mbp)) deallocate(Mbp)
  if (allocated(Mbq)) deallocate(Mbq)
end subroutine dealloc_all

subroutine reset_var_zero
  !Nap=0; Naq=0; Nbp=0; Nbq=0
  Map=0; Maq=0; Mbp=0; Mbq=0
end subroutine reset_var_zero

!-----------------------------------------------------------------------------
subroutine init_kernel(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer :: fid,n

fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)
call string_conf(fid,1,'snap_id',2,kid)

!call string_conf(fid,1,'SGT_LIST',2,sgt_list)
close(fid)
end subroutine init_kernel

subroutine kernelinfo_skel(filenm,nx,ny,nz)
character (len=*) :: filenm
integer,intent(in) :: nx,ny,nz

integer :: ierr,oldMode
integer :: Iid,Jid,Kid,numid,nstnmid,nevnmid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'kernelinfo_skel:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in kernelinfo_skel')
! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, Iid)
     call nfseis_except(ierr,'I dim in kernelinfo_skel')
ierr=nf90_def_dim(ncid, 'J', ny, Jid)
ierr=nf90_def_dim(ncid, 'K', nz, Kid)
ierr=nf90_def_dim(ncid, 'maxkstnm', MAXKSTNM, nstnmid)
ierr=nf90_def_dim(ncid, 'maxkevnm', MAXKEVNM, nevnmid)
ierr=nf90_def_dim(ncid, 'N', nf90_unlimited, numid)
     call nfseis_except(ierr,'N dim in kernelinfo_skel')
ierr=nf90_put_att(ncid,NF90_GLOBAL,'threshold_Kap',Kap0)
ierr=nf90_put_att(ncid,NF90_GLOBAL,'threshold_Kaq',Kaq0)
ierr=nf90_put_att(ncid,NF90_GLOBAL,'threshold_Kbp',Kbp0)
ierr=nf90_put_att(ncid,NF90_GLOBAL,'threshold_Kbq',Kbq0)
!--
ierr=nf90_def_var(ncid,'Sap',SEISNC_DATATYPE,(/Iid,Jid,Kid/),Sid(1))
ierr=nf90_def_var(ncid,'Saq',SEISNC_DATATYPE,(/Iid,Jid,Kid/),Sid(2))
ierr=nf90_def_var(ncid,'Sbp',SEISNC_DATATYPE,(/Iid,Jid,Kid/),Sid(3))
ierr=nf90_def_var(ncid,'Sbq',SEISNC_DATATYPE,(/Iid,Jid,Kid/),Sid(4))
ierr=nf90_def_var(ncid,'Nap',nf90_int,(/Iid,Jid,Kid/),Nid(1))
ierr=nf90_def_var(ncid,'Naq',nf90_int,(/Iid,Jid,Kid/),Nid(2))
ierr=nf90_def_var(ncid,'Nbp',nf90_int,(/Iid,Jid,Kid/),Nid(3))
ierr=nf90_def_var(ncid,'Nbq',nf90_int,(/Iid,Jid,Kid/),Nid(4))
ierr=nf90_def_var(ncid,'Map',nf90_int,(/Iid,Jid,Kid,numid/),Mid(1))
ierr=nf90_def_var(ncid,'Maq',nf90_int,(/Iid,Jid,Kid,numid/),Mid(2))
ierr=nf90_def_var(ncid,'Mbp',nf90_int,(/Iid,Jid,Kid,numid/),Mid(3))
ierr=nf90_def_var(ncid,'Mbq',nf90_int,(/Iid,Jid,Kid,numid/),Mid(4))
ierr=nf90_def_var(ncid,'KSTNMap',nf90_char,(/nstnmid,Iid,Jid,Kid,numid/),stnmid(1))
ierr=nf90_def_var(ncid,'KSTNMaq',nf90_char,(/nstnmid,Iid,Jid,Kid,numid/),stnmid(2))
ierr=nf90_def_var(ncid,'KSTNMbp',nf90_char,(/nstnmid,Iid,Jid,Kid,numid/),stnmid(3))
ierr=nf90_def_var(ncid,'KSTNMbq',nf90_char,(/nstnmid,Iid,Jid,Kid,numid/),stnmid(4))
ierr=nf90_def_var(ncid,'KEVNMap',nf90_char,(/nevnmid,Iid,Jid,Kid,numid/),evnmid(1))
ierr=nf90_def_var(ncid,'KEVNMaq',nf90_char,(/nevnmid,Iid,Jid,Kid,numid/),evnmid(2))
ierr=nf90_def_var(ncid,'KEVNMbp',nf90_char,(/nevnmid,Iid,Jid,Kid,numid/),evnmid(3))
ierr=nf90_def_var(ncid,'KEVNMbq',nf90_char,(/nevnmid,Iid,Jid,Kid,numid/),evnmid(4))
!--
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in kernelinfo_skel')
end subroutine kernelinfo_skel

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef KernelInfoMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef KernelInfoMPI
  print *, thisid
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

end program SI_ker2info

