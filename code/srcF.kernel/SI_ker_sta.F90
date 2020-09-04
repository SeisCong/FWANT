program SI_ker_sta

! This program calculates finite-frequency sensitivity kernels for multi pairs
! of event and station by swapping event strain tensors when keeping station
! sgt in memory.
! The core of kernel calculation coming from the code of Z Li, C Po, Z Zhang and Y Shen.
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

#define VERBOSE
!#define DEBUG

use constants_mod
use string_mod
use para_mod
use io_mod
use media_mod
use nfseis_mod
use src_mod
use mpi_mod
#ifdef KernelMPI
use mpi
#endif

implicit none

integer,parameter ::   &
  SIG_FILTER_P1  =100, &
  SIG_FILTER_P2  =120, &
  SIG_FILTER_NONE=130

character (len=8) :: str_d0,str_d1,str_d2
character (len=10):: str_t0,str_t1,str_t2

character (len=300) :: rcdstr
character (len=SEIS_STRLEN) ::           &
    fnm_main_conf,fnm_ker_conf,          &
    EVENT_ROOT,EVENT_OUTD,               &
    SGT_ROOT,SGT_STAX,SGT_STAY,SGT_STAZ, &
    evt_list,sgt_list,                   &
    fnm_synx,fnm_syny,fnm_synz,          &
    fnm_ker
character (len=SEIS_STRLEN) ::           &
    KSTANM,KEVNM
character (len=SEIS_STRLEN) :: filenm
character (len=15) :: idstr

character (len=SEIS_STRLEN),dimension(:),allocatable :: pnm_ker

integer :: staid,evtid

logical :: flag_fst
integer,dimension(SEIS_GEO) :: &
  bstart,bcount,bstride,       &
  sub_gs,sub_gc,sub_gt,        &
  subnc_s,subnc_c,subnc_t
integer :: sub_gtinv
integer,dimension(SEIS_GEO) :: blksiz,blkfst,blknum
integer :: kid,knt,nt1,nt2
real(SP) :: kdt,Tmin,Tmax

integer :: num_freq,num_pick
integer :: ncmp,npk,nfq
integer :: picknum(6)
integer,allocatable :: &
  picktidx(:,:),       &
  pickfreq(:),         &
  pickfilt(:)
real(SP),allocatable :: &
  picktwin(:,:)
real(SP),dimension(3,3) :: &
  rot_mat
character (len=2) :: filt_nm
real(SP) :: t1,t2

integer :: NOD,NOD1,NEL,NFACT
real(DP),dimension(:,:),allocatable :: filt_a,filt_b
real(DP),dimension(:),allocatable :: vecx,vecy,vecw

real(SP),dimension(:),pointer :: Vpr,Upr
real(SP),dimension(:),allocatable,target :: T,Vx,Ux,Vy,Uy,Vz,Uz
real(SP),dimension(:),allocatable :: VV,UU

real(SP) :: SGT_M0
real(SP),dimension(:,:,:,:),pointer :: Udapr,Udbpr
real(SP),dimension(:,:,:,:),allocatable :: &
    ExxX,EyyX,EzzX,ExyX,ExzX,EyzX,         &
    ExxY,EyyY,EzzY,ExyY,ExzY,EyzY,         &
    ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ
real(SP),dimension(:,:,:,:),allocatable,target :: &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS
real(SP),dimension(:,:,:,:),allocatable :: &
    Uda,Udb
logical :: flag_udx,flag_udy,flag_udz
logical :: flag_sgtx,flag_sgty,flag_sgtz

real(SP),dimension(:),allocatable :: Ka,Kb
real(SP),dimension(:),allocatable :: E11R,E22R,E33R,E12R,E13R,E23R
real(SP),dimension(:),allocatable :: E11S,E22S,E33S,E12S,E13S,E23S
real(SP),dimension(:,:,:),allocatable :: Kap,Kaq,Kbp,Kbq

integer :: n_i,n_j,n_k
integer :: blki,blkj,blkk
integer :: i,j,k,m,n,mt,ierr,i1,i2,j1,j2,k1,k2
integer :: p(1)
integer,dimension(SEIS_GEO) :: &
   bsubs,bsubc,bsubt,          &
    subs, subc, subt
!#ifdef KernelMPI
!integer,dimension(MPI_STATUS_SIZE) :: istatus
!#endif

integer :: ncidS,ncidF
integer :: TxxSid,TyySid,TzzSid,TxySid,TxzSid,TyzSid
integer :: TxxFid,TyyFid,TzzFid,TxyFid,TxzFid,TyzFid

integer,dimension(:),allocatable :: ncid
integer,dimension(:),allocatable :: kapid,kaqid,kbpid,kbqid

!----------------------------------------------------------------------

#ifdef KernelMPI
call MPI_INIT(ierr)
#endif

call get_conf_name(fnm_conf)

! read kernel conf
fnm_ker_conf='TomoKernel.conf'
call init_kernel(fnm_ker_conf)
fnm_conf=fnm_main_conf

call para_init(fnm_conf)
call swmpi_init(fnm_conf)

#ifdef KernelMPI
call swmpi_cart_creat
call swmpi_reinit_para
call swmpi_datatype
  if (masternode) then
     read(*,"(a)") sgt_list
     write(*,*) sgt_list
  end if
  call MPI_BCAST(sgt_list,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
#else
call swmpi_set_gindx(0,0,0)
  read(*,"(a)") sgt_list
  write(*,*) sgt_list
#endif

call media_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

call init_filter(fnm_ker_conf)

!----------------------------------------------------------------------

#ifndef KernelMPI
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

sub_gs=(bstart-1)*snap_subt(:,kid)+snap_subs(:,kid)
sub_gs(1)=out_i(sub_gs(1))
sub_gs(2)=out_j(sub_gs(2))
sub_gs(3)=out_k(sub_gs(3))
sub_gc=(snap_subc(:,kid)-bstart(:))/bstride(:)+1
do n=1,SEIS_GEO
   if (bcount(n)/=-1 .and. bcount(n)<sub_gc(n)) then
       sub_gc(n)=bcount(n)
   end if
end do
sub_gt=bstride*snap_subt(:,kid)
sub_gtinv=snap_tinv(kid)

call io_snap_locate(n_i,n_j,n_k)

!----------------------------------------------------------------------
kernode_if : if (snap_ishere(kid) .and. sgt_out(kid)) then

do i=1,sub_gc(1)
   i1=sub_gs(1)+(i-1)*sub_gt(1)
   if (i1>=snap_gsubs(1,kid)) exit
end do
do i=sub_gc(1),1,-1
   i2=sub_gs(1)+(i-1)*sub_gt(1)
   if (i2<=snap_gsube(1,kid)) exit
end do
subnc_s(1)=(i1-snap_gsubs(1,kid))/snap_subt(1,kid)+1
subnc_c(1)=(i2-i1)/sub_gt(1)+1

do j=1,sub_gc(2)
   j1=sub_gs(2)+(j-1)*sub_gt(2)
   if (j1>=snap_gsubs(2,kid)) exit
end do
do j=sub_gc(2),1,-1
   j2=sub_gs(2)+(j-1)*sub_gt(2)
   if (j2<=snap_gsube(2,kid)) exit
end do
subnc_s(2)=(j1-snap_gsubs(2,kid))/snap_subt(2,kid)+1
subnc_c(2)=(j2-j1)/sub_gt(2)+1

do k=1,sub_gc(3)
   k1=sub_gs(3)+(k-1)*sub_gt(3)
   if (k1>=snap_gsubs(3,kid)) exit
end do
do k=sub_gc(3),1,-1
   k2=sub_gs(3)+(k-1)*sub_gt(3)
   if (k2<=snap_gsube(3,kid)) exit
end do
subnc_s(3)=(k1-snap_gsubs(3,kid))/snap_subt(3,kid)+1
subnc_c(3)=(k2-k1)/sub_gt(3)+1

subnc_t=bstride

where (blksiz==-1 .or. blksiz>subnc_c)
   blksiz=subnc_c
end where
blknum=(subnc_c+blksiz-1)/blksiz

#ifdef DEBUG
write(100+myid,*) "thisid="
write(100+myid,*) thisid
write(100+myid,*) "sub_gs="
write(100+myid,*) sub_gs
write(100+myid,*) "sub_gc="
write(100+myid,*) sub_gc
write(100+myid,*) "sub_gt="
write(100+myid,*) sub_gt
write(100+myid,*) "sub_gtinv="
write(100+myid,*) sub_gtinv
write(100+myid,*) "subnc_s="
write(100+myid,*) subnc_s
write(100+myid,*) "subnc_c="
write(100+myid,*) subnc_c
write(100+myid,*) "subnc_t="
write(100+myid,*) subnc_t

write(100+myid,*) "blksiz="
write(100+myid,*) blksiz
write(100+myid,*) "blknum="
write(100+myid,*) blknum
#endif

call alloc_ker_var(blksiz(1),blksiz(2),blksiz(3),knt)
call alloc_media_local(blksiz(1),blksiz(2),blksiz(3))
call alloc_pick(num_pick)
call alloc_seismo_var(knt)

kdt=sub_gtinv*stept
T=(/1:knt/)*kdt
Tmin=kdt
Tmax=knt*kdt

staid=10001; evtid=10002
! loop all stations
open(staid,file=trim(sgt_list),status='old',iostat=ierr)
if (ierr>0) call error_except("sgt_list open err:"//trim(sgt_list))

!-----------------------------------------------------------------------------
sta_loop: do

  read(staid,*,iostat=ierr) KSTANM,evt_list
  if (ierr<0) exit sta_loop
  if (KSTANM(1:1)=='#') cycle

#ifdef VERBOSE
  call date_and_time(date=str_d1,time=str_t1)
  write(*,*) idstr//trim(KSTANM)//" "//trim(evt_list), &
                    ': begins from ',str_d1,  &
                  ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
#endif

! loop each block
!-----------------------------------------------------------------------------
do blkk=blkfst(3),blknum(3)
do blkj=blkfst(2),blknum(2)
do blki=blkfst(1),blknum(1)

#ifdef VERBOSE
  call date_and_time(date=str_d1,time=str_t1)
   write(*,"(a7,3(i4,a1,i4.4),a9,3(i2.2))")   &
       ' block:',                             &
       blki,'/',blknum(1),                    &
       blkj,'/',blknum(2),                    &
       blkk,'/',blknum(3),                    &
       ', thisid=',n_i,n_j,n_k
  write(*,*) idstr, ': current time ',str_d1, &
           ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
#endif

   flag_sgtx=.false.; flag_sgty=.false.; flag_sgtz=.false.

   bsubs(1)=(blki-1)*blksiz(1)*subnc_t(1)+subnc_s(1)
   bsubs(2)=(blkj-1)*blksiz(2)*subnc_t(2)+subnc_s(2)
   bsubs(3)=(blkk-1)*blksiz(3)*subnc_t(3)+subnc_s(3)

   bsubc(1)=min(blksiz(1),subnc_c(1)-(bsubs(1)-subnc_s(1))/subnc_t(1))
   bsubc(2)=min(blksiz(2),subnc_c(2)-(bsubs(2)-subnc_s(2))/subnc_t(2))
   bsubc(3)=min(blksiz(3),subnc_c(3)-(bsubs(3)-subnc_s(3))/subnc_t(3))

   bsubt=subnc_t

   subc=bsubc;
   subt=snap_subt(:,kid)*bsubt
   subs=snap_subs(:,kid)+(bsubs-1)*snap_subt(:,kid)

#ifdef DEBUG
write(100+myid,*) "-------------------"
write(100+myid,*) "bli,j,k=",blki,blkj,blkk
write(100+myid,*) "bsubs="
write(100+myid,*) bsubs
write(100+myid,*) "bsubc="
write(100+myid,*) bsubc
write(100+myid,*) "bsubt="
write(100+myid,*) bsubt
write(100+myid,*) "subs="
write(100+myid,*) subs
write(100+myid,*) "subc="
write(100+myid,*) subc
write(100+myid,*) "subt="
write(100+myid,*) subt
#endif

if ( product(bsubc)*knt/=size(ExxS) ) then
   call realloc_ker_var(bsubc(1),bsubc(2),bsubc(3),knt)
   call alloc_media_local(bsubc(1),bsubc(2),bsubc(3))
end if

! load media
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt) 

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

  if (maxval(picknum)>num_pick) then
     num_pick=maxval(picknum)
     call alloc_pick(num_pick)
  end if

  flag_udx=.false.;flag_udy=.false.;flag_udz=.false.
  if (picknum(1)>0 .or. any(picknum(4:6)>0)) flag_udx=.true.
  if (picknum(2)>0 .or. any(picknum(4:6)>0)) flag_udy=.true.
  if (picknum(3)>0 .or. any(picknum(4:6)>0)) flag_udz=.true.

  read(evtid,"(a132)") fnm_synx
  read(evtid,"(a132)") fnm_syny
  read(evtid,"(a132)") fnm_synz
  if (flag_udx) call init_synthetic(fnm_synx,T,Vx,Ux)
  if (flag_udy) call init_synthetic(fnm_syny,T,Vy,Uy)
  if (flag_udz) call init_synthetic(fnm_synz,T,Vz,Uz)

  if (all(picknum<=0)) cycle

  if (abs(T(2)-T(1)-kdt)>SEIS_EQUAL) then
     print *, kdt,T(2),T(1)
     print *, T(2)-T(1),abs(T(2)-T(1))
     !call error_except('dt /= stept*tinv')
  end if

! load station sgtx
if ( flag_udx .and. (.not. flag_sgtx) ) then
   if ((.not. allocated(ExxX)) .or. (product(bsubc)*knt/=size(ExxX))) then
      call alloc_sgtx(bsubc(1),bsubc(2),bsubc(3),knt)
   end if
if (snap_tcnt(kid)<knt) then
   m=0; n=0; mt=snap_tcnt(kid)
   do 
      if (m>=knt) exit
      n=n+1
      filenm=get_fnm_snapnode_n(                                 &
         trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAX), &
         'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
      if (m+mt>knt) mt=knt-m

      ierr=nf90_get_var(ncidF,TxxFid,ExxX(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyyFid,EyyX(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TzzFid,EzzX(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxyFid,ExyX(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxzFid,ExzX(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyzFid,EyzX(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel: read sgtx fail')
      m=m+mt
      call nfseis_close(ncidF)
   end do
else
   filenm=get_fnm_snapnode_n(                                 &
      trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAX), &
      'sgt_',kid,1,n_i,n_j,n_k)
   call nfseis_open(filenm,ncidF)
   call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
   call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
   call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
   call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
   call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
   call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
   ierr=nf90_get_var(ncidF,TxxFid,ExxX,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyyFid,EyyX,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TzzFid,EzzX,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxyFid,ExyX,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxzFid,ExzX,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyzFid,EyzX,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   call nfseis_close(ncidF)
end if
  flag_sgtx=.true.
  call stress2strain(ExxX,EyyX,EzzX,ExyX,ExzX,EyzX,lambda,mu,bsubc)
end if
! load station sgty
if ( flag_udy .and. (.not. flag_sgty) ) then
   if ((.not. allocated(ExxY)) .or. (product(bsubc)*knt/=size(ExxY))) then
      call alloc_sgty(bsubc(1),bsubc(2),bsubc(3),knt)
   end if
if (snap_tcnt(kid)<knt) then
   m=0; n=0; mt=snap_tcnt(kid)
   do 
      if (m>=knt) exit
      n=n+1
      filenm=get_fnm_snapnode_n(                                 &
         trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAY), &
         'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
      if (m+mt>knt) mt=knt-m

      ierr=nf90_get_var(ncidF,TxxFid,ExxY(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyyFid,EyyY(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TzzFid,EzzY(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxyFid,ExyY(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxzFid,ExzY(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyzFid,EyzY(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel: read sgtx fail')
      m=m+mt
      call nfseis_close(ncidF)
   end do
else
   filenm=get_fnm_snapnode_n(                                 &
      trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAY), &
      'sgt_',kid,1,n_i,n_j,n_k)
   call nfseis_open(filenm,ncidF)
   call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
   call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
   call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
   call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
   call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
   call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
   ierr=nf90_get_var(ncidF,TxxFid,ExxY,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyyFid,EyyY,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TzzFid,EzzY,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxyFid,ExyY,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxzFid,ExzY,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyzFid,EyzY,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   call nfseis_close(ncidF)
end if
  flag_sgty=.true.
  call stress2strain(ExxY,EyyY,EzzY,ExyY,ExzY,EyzY,lambda,mu,bsubc)
end if
! load station sgtz
if ( flag_udz .and. (.not. flag_sgtz) ) then
   if ((.not. allocated(ExxZ)) .or. (product(bsubc)*knt/=size(ExxZ))) then
      call alloc_sgtz(bsubc(1),bsubc(2),bsubc(3),knt)
   end if
if (snap_tcnt(kid)<knt) then
   m=0; n=0; mt=snap_tcnt(kid)
   do 
      if (m>=knt) exit
      n=n+1
      filenm=get_fnm_snapnode_n(                                 &
         trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAZ), &
         'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
      if (m+mt>knt) mt=knt-m

      ierr=nf90_get_var(ncidF,TxxFid,ExxZ(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyyFid,EyyZ(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TzzFid,EzzZ(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxyFid,ExyZ(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxzFid,ExzZ(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyzFid,EyzZ(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel: read sgtx fail')
      m=m+mt
      call nfseis_close(ncidF)
   end do
else
   filenm=get_fnm_snapnode_n(                                 &
      trim(SGT_ROOT)//'/'//trim(KSTANM)//'/'//trim(SGT_STAZ), &
      'sgt_',kid,1,n_i,n_j,n_k)
   call nfseis_open(filenm,ncidF)
   call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
   call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
   call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
   call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
   call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
   call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
   ierr=nf90_get_var(ncidF,TxxFid,ExxZ,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyyFid,EyyZ,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TzzFid,EzzZ,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxyFid,ExyZ,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxzFid,ExzZ,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyzFid,EyzZ,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   call nfseis_close(ncidF)
end if
  flag_sgtz=.true.
  call stress2strain(ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ,lambda,mu,bsubc)
end if
! load event sgt
if (snap_tcnt(kid)<knt) then
   m=0; n=0; mt=snap_tcnt(kid)
   do 
      if (m>=knt) exit
      n=n+1
      filenm=get_fnm_snapnode_n(                                    &
         trim(EVENT_ROOT)//'/'//trim(KEVNM)//'/'//trim(EVENT_OUTD), &
         'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidS)
      call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
      call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
      call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
      call nfseis_inq_varid(ncidS, 'Txy', TxySid)
      call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
      call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)
      if (m+mt>knt) mt=knt-m
      ierr=nf90_get_var(ncidS,TxxSid,ExxS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyySid,EyyS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TzzSid,EzzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxySid,ExyS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxzSid,ExzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyzSid,EyzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel: read source sgt fail')
      m=m+mt
      call nfseis_close(ncidS)
   end do
else
   filenm=get_fnm_snapnode_n(                                    &
      trim(EVENT_ROOT)//'/'//trim(KEVNM)//'/'//trim(EVENT_OUTD), &
      'sgt_',kid,1,n_i,n_j,n_k)
   call nfseis_open(filenm,ncidS)
   call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
   call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
   call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
   call nfseis_inq_varid(ncidS, 'Txy', TxySid)
   call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
   call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)
   ierr=nf90_get_var(ncidS,TxxSid,ExxS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyySid,EyyS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TzzSid,EzzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxySid,ExyS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxzSid,ExzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyzSid,EyzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   call nfseis_close(ncidS)
end if
  call stress2strain(ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,lambda,mu,bsubc)

! scatter
call sgt2scatter3(                  &
     ExxX,EyyX,EzzX,ExyX,ExzX,EyzX, &
     ExxY,EyyY,EzzY,ExyY,ExzY,EyzY, &
     ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ, &
     ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
     flag_udx,flag_udy,flag_udz,    &
     lambda,mu,1,knt,kdt,bsubc)

! loop component
do ncmp=1,6

   if (ncmp==4 .and. any(picknum(4:6)>0)) then
      call syn_rotate(Vx,Vy,Vz,rot_mat,1,knt)
      call syn_rotate(Ux,Uy,Uz,rot_mat,1,knt)
      call scatter_rotate(ExxS,EyyS,EzzS,rot_mat,1,knt,bsubc)
      call scatter_rotate(ExyS,EyzS,ExzS,rot_mat,1,knt,bsubc)
   end if

   if (picknum(ncmp)<=0) cycle

   if (ncmp==1) then
      Vpr=>Vx; Upr=>Ux;
      Udapr=>ExxS; Udbpr=>ExyS
   elseif (ncmp==2) then
      Vpr=>Vy; Upr=>Uy;
      Udapr=>EyyS; Udbpr=>EyzS
   elseif (ncmp==3) then
      Vpr=>Vz; Upr=>Uz;
      Udapr=>EzzS; Udbpr=>ExzS
   elseif (ncmp==4) then
      Vpr=>Vx; Upr=>Ux;
      Udapr=>ExxS; Udbpr=>ExyS
   elseif (ncmp==5) then
      Vpr=>Vy; Upr=>Uy;
      Udapr=>EyyS; Udbpr=>EyzS
   elseif (ncmp==6) then
      Vpr=>Vz; Upr=>Uz;
      Udapr=>EzzS; Udbpr=>ExzS
   end if

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
   read(rcdstr,*) pickfreq(npk)
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   if (pickfreq(npk)>num_freq) call error_except('frequency band id overflow')
   ! pickfilt
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   filt_nm=rcdstr(1:n)
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   select case (trim(filt_nm))
   case ('P1')
     pickfilt(npk)=SIG_FILTER_P1
   case ('P2')
     pickfilt(npk)=SIG_FILTER_P2
   case ('N')
     pickfilt(npk)=SIG_FILTER_NONE
   end select
   ! t1 t2
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) t1
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
   read(rcdstr,*) t2
   rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
   if (t1>Tmax) then
      write(*,*) "t1,Tmax=",t1,Tmax
      call warning_print("t1>Tmax, skip this rcd "//trim(rcdstr))
      cycle
   elseif (t1<=Tmin) then
      write(*,*) "t1,Tmin=",t1,Tmin
      call warning_print("t1<Tmin, reset t1 to Tmin"//trim(rcdstr))
      picktidx(1,npk)=1
      picktwin(1,npk)=T(1)
   else
      p=maxloc(T,T<=t1); picktidx(1,npk)=p(1); picktwin(1,npk)=T(p(1));
   end if
   if (t2<Tmin) then
      write(*,*) "t2,Tmin=",t2,Tmin
      call warning_print("t2<Tmin, skip this rcd"//trim(rcdstr))
      cycle
   elseif (t2>=Tmax) then
      write(*,*) "t2,Tmax=",t2,Tmax
      call warning_print("t2>Tmax, reset t2 to Tmax"//trim(rcdstr))
      picktidx(2,npk)=knt
      picktwin(2,npk)=T(knt)
   else
      p=minloc(T,T>=t2); picktidx(2,npk)=p(1); picktwin(2,npk)=T(p(1));
   end if
   ! pnm_ker
   pnm_ker(npk)=trim(rcdstr)

   ! kernel nc file
   if (flag_fst .and. (blki==1 .and. blkj==1 .and. blkk==1)) then
      fnm_ker=get_fnm_snapnode_n('./','kernel_',kid,0,n_i,n_j,n_k)
      call nfseis_grid3d_def(trim(pnm_ker(npk))//"/"//trim(fnm_ker),   &
          snap_subc(1,kid),snap_subc(2,kid), snap_subc(3,kid),ncid(npk), &
          "Finite frequency kernel")
      call nfseis_grid3d_defvar(ncid(npk),'phase_Vp'    ,kapid(npk))
      call nfseis_grid3d_defvar(ncid(npk),'amplitude_Vp',kaqid(npk))
      call nfseis_grid3d_defvar(ncid(npk),'phase_Vs'    ,kbpid(npk))
      call nfseis_grid3d_defvar(ncid(npk),'amplitude_Vs',kbqid(npk))
      call nfseis_grid3d_enddef(ncid(npk))
   else
      fnm_ker=get_fnm_snapnode_n('./','kernel_',kid,0,n_i,n_j,n_k)
      call nfseis_open(trim(pnm_ker(npk))//"/"//trim(fnm_ker),ncid(npk))
      call nfseis_inq_varid(ncid(npk),'phase_Vp'    ,kapid(npk))
      call nfseis_inq_varid(ncid(npk),'amplitude_Vp',kaqid(npk))
      call nfseis_inq_varid(ncid(npk),'phase_Vs'    ,kbpid(npk))
      call nfseis_inq_varid(ncid(npk),'amplitude_Vs',kbqid(npk))
   end if
end do

! cal kernel
do npk=1,picknum(ncmp)
   nfq=pickfreq(npk)

   ! filtering when different frq or filt
   if (npk==1  &
       .or. (nfq/=pickfreq(npk-1)) .or. (pickfilt(npk)/=pickfilt(npk-1)) ) then
   select case (pickfilt(npk))
   case (SIG_FILTER_P1)
     call scatter_filter(Udapr,Uda,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call scatter_filter(Udbpr,Udb,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call syn_filter(Vpr,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
     call syn_filter(Upr,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
   case (SIG_FILTER_P2)
     call scatter_filtfilt(Udapr,Uda,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call scatter_filtfilt(Udbpr,Udb,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
     call syn_filtfilt(Vpr,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
     call syn_filtfilt(Upr,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
   case (SIG_FILTER_NONE)
     Uda=Udapr;Udb=Udbpr;
     VV=Vpr; UU=Upr;
   case default
     print *, "npk,nfq,pickfilt=",npk,nfq,pickfilt(npk)
     call error_except('pickfilt error')
   end select
   end if

#ifdef VERBOSE
#ifdef KernelMPI
if (myid==0) &
#endif
   call export_info(pnm_ker(npk),T,VV,UU, &
        picktwin(:,npk),picktidx(:,npk),  &
        filt_a(:,nfq),filt_b(:,nfq))
#ifdef DEBUG
#ifdef KernelMPI
if (myid==0) &
#endif
   call export_debug(pnm_ker(npk),T,Vpr,Upr)
#endif
#endif

   nt1=picktidx(1,npk);nt2=picktidx(2,npk)
   call scatter2kernel(VV,UU,Uda,Udb,nt1,nt2,kdt,bsubc,Kap,Kaq,Kbp,Kbq)

   ! put
   call nfseis_put(ncid(npk),kapid(npk),Kap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(npk),kaqid(npk),Kaq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(npk),kbpid(npk),Kbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(npk),kbqid(npk),Kbq,bsubs,bsubc,bsubt)
end do ! npk

! close nc file
do npk=1,picknum(ncmp)
   call nfseis_close(ncid(npk))
end do ! npk

end do ! ncmp

end do evt_loop
close(evtid)

end do ! blki
end do ! blkj
end do ! blkk

end do sta_loop
close(staid)

!-----------------------------------------------------------------------------
end if kernode_if

#ifdef VERBOSE
  call date_and_time(date=str_d2,time=str_t2)
  write(*,*) idstr, ' finish at ',str_d2,  &
                     ',',str_t2(1:2),'/',str_t2(3:4),'/',str_t2(5:10)
#endif

!-----------------------------------------------------------------------------
call media_destroy
call dealloc_all

#ifdef KernelMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine alloc_filt_var(nf,NOD1,NEL)
  integer,intent(in) :: nf,NOD1,NEL
  allocate(filt_b(NOD1,nf)); filt_b=0.0_DP
  allocate(filt_a(NOD1,nf)); filt_a=0.0_DP
  allocate(vecx(NEL)); vecx=0.0_DP
  allocate(vecy(NEL)); vecy=0.0_DP
  allocate(vecw(NEL)); vecw=0.0_DP
end subroutine alloc_filt_var

subroutine alloc_pick(npk)
  integer,intent(in) :: npk
  if (allocated(pnm_ker)) deallocate(pnm_ker)
  if (allocated(pickfreq)) deallocate(pickfreq)
  if (allocated(pickfilt)) deallocate(pickfilt)
  if (allocated(picktidx)) deallocate(picktidx)
  if (allocated(picktwin)) deallocate(picktwin)
  if (allocated(ncid ))    deallocate(ncid )
  if (allocated(kapid))    deallocate(kapid)
  if (allocated(kaqid))    deallocate(kaqid)
  if (allocated(kbpid))    deallocate(kbpid)
  if (allocated(kbqid))    deallocate(kbqid)

  allocate(pnm_ker(npk)); 
  allocate(pickfreq(npk)); pickfreq=0
  allocate(pickfilt(npk)); pickfilt=0
  allocate(picktidx(2,npk)); picktidx=0
  allocate(picktwin(2,npk)); picktwin=0.0
  allocate(ncid (npk)); ncid=0
  allocate(kapid(npk)); kapid=0
  allocate(kaqid(npk)); kaqid=0
  allocate(kbpid(npk)); kbpid=0
  allocate(kbqid(npk)); kbqid=0
end subroutine alloc_pick

subroutine alloc_seismo_var(nt)
  integer,intent(in) :: nt
  allocate(T(nt)); T=0.0
  allocate(VV(nt)); VV=0.0
  allocate(UU(nt)); UU=0.0
  allocate(Vx(nt)); Vx=0.0
  allocate(Ux(nt)); Ux=0.0
  allocate(Vy(nt)); Vy=0.0
  allocate(Uy(nt)); Uy=0.0
  allocate(Vz(nt)); Vz=0.0
  allocate(Uz(nt)); Uz=0.0
end subroutine alloc_seismo_var

subroutine alloc_ker_var(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0

  allocate(Uda(ki,kj,kk,kt)); Uda=0.0
  allocate(Udb(ki,kj,kk,kt)); Udb=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0

  allocate(Ka(kt),Kb(kt));Ka=0.0;Kb=0.0
  allocate(E11S(kt)); E11S=0.0
  allocate(E22S(kt)); E22S=0.0
  allocate(E33S(kt)); E33S=0.0
  allocate(E12S(kt)); E12S=0.0
  allocate(E13S(kt)); E13S=0.0
  allocate(E23S(kt)); E23S=0.0
  allocate(E11R(kt)); E11R=0.0
  allocate(E22R(kt)); E22R=0.0
  allocate(E33R(kt)); E33R=0.0
  allocate(E12R(kt)); E12R=0.0
  allocate(E13R(kt)); E13R=0.0
  allocate(E23R(kt)); E23R=0.0
end subroutine alloc_ker_var

subroutine realloc_ker_var(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' realloc_ker_var'
#endif
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
  if (allocated(Uda)) deallocate(Uda)
  if (allocated(Udb)) deallocate(Udb)

  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0

  allocate(Uda(ki,kj,kk,kt)); Uda=0.0
  allocate(Udb(ki,kj,kk,kt)); Udb=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0 
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0 
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0 
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0 
end subroutine realloc_ker_var

subroutine alloc_sgtx(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' alloc_sgtx'
#endif
  if (allocated(ExxX)) deallocate(ExxX)
  if (allocated(EyyX)) deallocate(EyyX)
  if (allocated(EzzX)) deallocate(EzzX)
  if (allocated(ExyX)) deallocate(ExyX)
  if (allocated(ExzX)) deallocate(ExzX)
  if (allocated(EyzX)) deallocate(EyzX)
  allocate(ExxX(ki,kj,kk,kt)); ExxX=0.0
  allocate(EyyX(ki,kj,kk,kt)); EyyX=0.0
  allocate(EzzX(ki,kj,kk,kt)); EzzX=0.0
  allocate(ExyX(ki,kj,kk,kt)); ExyX=0.0
  allocate(ExzX(ki,kj,kk,kt)); ExzX=0.0
  allocate(EyzX(ki,kj,kk,kt)); EyzX=0.0
end subroutine alloc_sgtx
subroutine alloc_sgty(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' alloc_sgty'
#endif
  if (allocated(ExxY)) deallocate(ExxY)
  if (allocated(EyyY)) deallocate(EyyY)
  if (allocated(EzzY)) deallocate(EzzY)
  if (allocated(ExyY)) deallocate(ExyY)
  if (allocated(ExzY)) deallocate(ExzY)
  if (allocated(EyzY)) deallocate(EyzY)
  allocate(ExxY(ki,kj,kk,kt)); ExxY=0.0
  allocate(EyyY(ki,kj,kk,kt)); EyyY=0.0
  allocate(EzzY(ki,kj,kk,kt)); EzzY=0.0
  allocate(ExyY(ki,kj,kk,kt)); ExyY=0.0
  allocate(ExzY(ki,kj,kk,kt)); ExzY=0.0
  allocate(EyzY(ki,kj,kk,kt)); EyzY=0.0
end subroutine alloc_sgty
subroutine alloc_sgtz(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
#ifdef VERBOSE
  write(*,*) idstr,' alloc_sgtz'
#endif
  if (allocated(ExxZ)) deallocate(ExxZ)
  if (allocated(EyyZ)) deallocate(EyyZ)
  if (allocated(EzzZ)) deallocate(EzzZ)
  if (allocated(ExyZ)) deallocate(ExyZ)
  if (allocated(ExzZ)) deallocate(ExzZ)
  if (allocated(EyzZ)) deallocate(EyzZ)
  allocate(ExxZ(ki,kj,kk,kt)); ExxZ=0.0
  allocate(EyyZ(ki,kj,kk,kt)); EyyZ=0.0
  allocate(EzzZ(ki,kj,kk,kt)); EzzZ=0.0
  allocate(ExyZ(ki,kj,kk,kt)); ExyZ=0.0
  allocate(ExzZ(ki,kj,kk,kt)); ExzZ=0.0
  allocate(EyzZ(ki,kj,kk,kt)); EyzZ=0.0
end subroutine alloc_sgtz

subroutine alloc_media_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  if (allocated(lambda)) deallocate(lambda)
  if (allocated(mu)) deallocate(mu)
  allocate(lambda(ki,kj,kk)); lambda=0.0
  allocate(mu(ki,kj,kk)); mu=0.0
end subroutine alloc_media_local

subroutine dealloc_all
  if (allocated(Vx)) deallocate(Vx)
  if (allocated(Ux)) deallocate(Ux)
  if (allocated(T)) deallocate(T)
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxX)) deallocate(ExxX)
  if (allocated(EyyX)) deallocate(EyyX)
  if (allocated(EzzX)) deallocate(EzzX)
  if (allocated(ExyX)) deallocate(ExyX)
  if (allocated(ExzX)) deallocate(ExzX)
  if (allocated(EyzX)) deallocate(EyzX)
  if (allocated(Ka)) deallocate(Ka)
  if (allocated(Kb)) deallocate(Kb)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
end subroutine dealloc_all

!-----------------------------------------------------------------------------
subroutine init_kernel(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer :: fid,n

fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)
call string_conf(fid,1,'snap_id',2,kid)
call string_conf(fid,1,'first_run',2,flag_fst)
do n=1,SEIS_GEO
 call string_conf(fid,1,'sub_start',n+1,bstart(n))
 call string_conf(fid,1,'sub_count',n+1,bcount(n))
 call string_conf(fid,1,'sub_stride',n+1,bstride(n))
end do
do n=1,SEIS_GEO
 call string_conf(fid,1,'block_from',n+1,blkfst(n))
 call string_conf(fid,1,'block_size',n+1,blksiz(n))
end do

call string_conf(fid,1,"SGT_TIMEDIM_SIZE",2,knt)
call string_conf(fid,1,"SGT_M0",2,SGT_M0)
call string_conf(fid,1,'SGT_ROOT',2,SGT_ROOT)
call string_conf(fid,1,'SGT_STAX',2,SGT_STAX)
call string_conf(fid,1,'SGT_STAY',2,SGT_STAY)
call string_conf(fid,1,'SGT_STAZ',2,SGT_STAZ)
!call string_conf(fid,1,'SGT_LIST',2,sgt_list)

call string_conf(fid,1,'EVENT_ROOT',2,EVENT_ROOT)
call string_conf(fid,1,'EVENT_OUTD',2,EVENT_OUTD)

! numbers
call string_conf(fid,1,'maximum_of_pick_per_component',2,num_pick)
call string_conf(fid,1,'number_of_freqband',2,num_freq)

close(fid)
end subroutine init_kernel

subroutine init_filter(fnm_conf)
character (len=*),intent(in) :: fnm_conf

character (len=SEIS_STRLEN) :: fnm_filter,str
integer :: fid,gid,n,m

fid=5002; gid=5003
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'<filter_file>',1,str)
do n=1,num_freq
   read(fid,*) fnm_filter

   if (trim(fnm_filter)/='none') then
      open(gid,file=trim(fnm_filter),status='old')
      read(gid,*) NOD1
      NOD=NOD1-1; NFACT=3*NOD
      NEL=knt+2*NFACT
      if (.not. allocated(filt_a)) then
         call alloc_filt_var(num_freq,NOD1,NEL)
      end if
      do m=1,NOD1
         !read(gid,*) filt_a(m,n),filt_b(m,n)
         read(gid,"(2(f20.16))") filt_a(m,n),filt_b(m,n)
      end do
      close(gid)
   end if
end do
close(fid)
end subroutine init_filter

subroutine init_synthetic(fnm_syn,T,V,U)
character (len=*),intent(in) :: fnm_syn
real(SP),dimension(:),intent(out) :: T,V,U
integer :: fid,n

fid=5001
open(fid,file=trim(fnm_syn),status='old')
  do n=1,knt
     read(fid,*) T(n),V(n),U(n)
  end do
close(fid)
end subroutine init_synthetic

subroutine export_info(pnm,T,V,U,twin,tidx,fcta,fctb)
character (len=*),intent(in) :: pnm
real(SP),dimension(:),intent(in) :: T,V,U
real(SP),dimension(:),intent(in) :: twin
integer,dimension(:),intent(in) :: tidx
real(DP),dimension(:),intent(in) :: fcta,fctb

character (len=SEIS_STRLEN) :: filenm
integer :: n,fid
fid=6001

filenm=trim(pnm)//'/'//'kernel.para.log'
open(fid,file=trim(filenm),status='unknown')
  write(fid,*) twin
  write(fid,*) tidx
  write(fid,"(f20.16)") fcta
  write(fid,"(f20.16)") fctb
close(fid)

filenm=trim(pnm)//'/'//'synthetic.filtered.dat'
open(fid,file=trim(filenm),status='unknown')
do n=1,knt
   write(fid,*) T(n),V(n),U(n)
end do
close(fid)
end subroutine export_info

subroutine export_debug(pnm,T,V,U)
character (len=*),intent(in) :: pnm
real(SP),dimension(:),intent(in) :: T,V,U

character (len=SEIS_STRLEN) :: filenm
integer :: n,fid
fid=6001

filenm=trim(pnm)//'/'//'kernel.debug.log'
open(fid,file=trim(filenm),status='unknown')
  write(fid,*) 'knt=',knt
  write(fid,*) 'NEL=',NEL
  write(fid,*) 'NFACT=',NFACT
  write(fid,*) 'NOD=',NOD
close(fid)

filenm=trim(pnm)//'/'//'synthetic.orig.dat'
open(fid,file=trim(filenm),status='unknown')
do n=1,knt
   write(fid,*) T(n),V(n),U(n)
end do
close(fid)

filenm=trim(pnm)//'/'//'filter.vec.dat'
open(fid,file=trim(filenm),status='unknown')
do n=1,NEL
   write(fid,*) vecx(n),vecy(n)
end do
close(fid)
end subroutine export_debug

subroutine stress2strain(Exx,Eyy,Ezz,Exy,Exz,Eyz,lam,miu,scl)
real(SP),dimension(:,:,:,:),intent(inout) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
real(SP),dimension(:,:,:),intent(in) :: lam,miu
integer,dimension(SEIS_GEO),intent(in) :: scl
real(SP) E1,E2,E3,E0
integer i,j,k,n,n0

n0=size(Exx,4)
do n=1,n0
do k=1,scl(3)
do j=1,scl(2)
do i=1,scl(1)
   E1=(lam(i,j,k)+miu(i,j,k))/(miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   E2=-lam(i,j,k)/(2.0*miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   E3=1.0/miu(i,j,k)
   E0=E2*(Exx(i,j,k,n)+Eyy(i,j,k,n)+Ezz(i,j,k,n))
   Exx(i,j,k,n)=E0-(E2-E1)*Exx(i,j,k,n)
   Eyy(i,j,k,n)=E0-(E2-E1)*Eyy(i,j,k,n)
   Ezz(i,j,k,n)=E0-(E2-E1)*Ezz(i,j,k,n)
   Exy(i,j,k,n)=0.5*E3*Exy(i,j,k,n)
   Exz(i,j,k,n)=0.5*E3*Exz(i,j,k,n)
   Eyz(i,j,k,n)=0.5*E3*Eyz(i,j,k,n)
end do
end do
end do
end do
end subroutine stress2strain

subroutine sgt2scatter(                   &
           ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
           ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
           lam,miu,nt1,nt2,dt,cblk)

integer,intent(in) :: nt1,nt2
real(SP),dimension(:,:,:,:),intent(inout) :: &
     ExxS,EyyS
real(SP),dimension(:,:,:,:),intent(in) ::    &
               EzzS,ExyS,ExzS,EyzS,          &
     ExxR,EyyR,EzzR,ExyR,ExzR,EyzR
real(SP),dimension(:,:,:),intent(in) :: miu,lam
integer,dimension(SEIS_GEO),intent(in) :: cblk
real(SP),intent(in) :: dt

integer :: n,i,j,k
real(kind=SP) :: x1,x2

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)
   E11R(1:nt2)=ExxR(i,j,k,1:nt2);E22R=EyyR(i,j,k,1:nt2);E33R=EzzR(i,j,k,1:nt2)
   E12R(1:nt2)=ExyR(i,j,k,1:nt2);E13R=ExzR(i,j,k,1:nt2);E23R=EyzR(i,j,k,1:nt2)
   E11S(1:nt2)=ExxS(i,j,k,1:nt2);E22S=EyyS(i,j,k,1:nt2);E33S=EzzS(i,j,k,1:nt2)
   E12S(1:nt2)=ExyS(i,j,k,1:nt2);E13S=ExzS(i,j,k,1:nt2);E23S=EyzS(i,j,k,1:nt2)

do n=nt1,nt2
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do

   ExxS(i,j,k,n)=x1*2.0_SP*(lam(i,j,k)+2.0_SP*miu(i,j,k))
   EyyS(i,j,k,n)=2.0_SP*(x2-x1)*2.0_SP*miu(i,j,k)

end do

end do
end do
end do
end subroutine sgt2scatter

subroutine sgt2scatter3(           &
    ExxX,EyyX,EzzX,ExyX,ExzX,EyzX, &
    ExxY,EyyY,EzzY,ExyY,ExzY,EyzY, &
    ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ, &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
    flag_udx,flag_udy,flag_udz,    &
    lam,miu,nt1,nt2,dt,cblk)

logical,intent(in) :: flag_udx,flag_udy,flag_udz
integer,intent(in) :: nt1,nt2
real(SP),dimension(:,:,:,:),intent(inout) :: &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS
real(SP),dimension(:,:,:,:),intent(in) ::    &
    ExxX,EyyX,EzzX,ExyX,ExzX,EyzX,           &
    ExxY,EyyY,EzzY,ExyY,ExzY,EyzY,           &
    ExxZ,EyyZ,EzzZ,ExyZ,ExzZ,EyzZ
real(SP),dimension(:,:,:),intent(in) :: miu,lam
integer,dimension(SEIS_GEO),intent(in) :: cblk
real(SP),intent(in) :: dt

integer :: n,i,j,k
real(kind=SP) :: x1,x2

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)
   E11S(1:nt2)=ExxS(i,j,k,1:nt2);E22S=EyyS(i,j,k,1:nt2);E33S=EzzS(i,j,k,1:nt2)
   E12S(1:nt2)=ExyS(i,j,k,1:nt2);E13S=ExzS(i,j,k,1:nt2);E23S=EyzS(i,j,k,1:nt2)

! Ux
if (flag_udx) then
   E11R(1:nt2)=ExxX(i,j,k,1:nt2);E22R=EyyX(i,j,k,1:nt2);E33R=EzzX(i,j,k,1:nt2)
   E12R(1:nt2)=ExyX(i,j,k,1:nt2);E13R=ExzX(i,j,k,1:nt2);E23R=EyzX(i,j,k,1:nt2)
do n=nt1,nt2
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do
   ExxS(i,j,k,n)=x1*2.0_SP*(lam(i,j,k)+2.0_SP*miu(i,j,k))
   ExyS(i,j,k,n)=2.0_SP*(x2-x1)*2.0_SP*miu(i,j,k)
end do
end if

! Uy
if (flag_udy) then
   E11R(1:nt2)=ExxY(i,j,k,1:nt2);E22R=EyyY(i,j,k,1:nt2);E33R=EzzY(i,j,k,1:nt2)
   E12R(1:nt2)=ExyY(i,j,k,1:nt2);E13R=ExzY(i,j,k,1:nt2);E23R=EyzY(i,j,k,1:nt2)
do n=nt1,nt2
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do
   EyyS(i,j,k,n)=x1*2.0_SP*(lam(i,j,k)+2.0_SP*miu(i,j,k))
   EyzS(i,j,k,n)=2.0_SP*(x2-x1)*2.0_SP*miu(i,j,k)
end do
end if

! Uz
if (flag_udz) then
   E11R(1:nt2)=ExxZ(i,j,k,1:nt2);E22R=EyyZ(i,j,k,1:nt2);E33R=EzzZ(i,j,k,1:nt2)
   E12R(1:nt2)=ExyZ(i,j,k,1:nt2);E13R=ExzZ(i,j,k,1:nt2);E23R=EyzZ(i,j,k,1:nt2)
do n=nt1,nt2
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do
   EzzS(i,j,k,n)=x1*2.0_SP*(lam(i,j,k)+2.0_SP*miu(i,j,k))
   ExzS(i,j,k,n)=2.0_SP*(x2-x1)*2.0_SP*miu(i,j,k)
end do
end if

end do
end do
end do
end subroutine sgt2scatter3

subroutine syn_rotate(VX,Vy,Vz,rot_mat,nt1,nt2)
real(SP),dimension(:),intent(inout) :: Vx,Vy,Vz
real(SP),dimension(:,:),intent(in) :: rot_mat
integer,intent(in) :: nt1,nt2

real(SP) :: v1,v2,v3
integer :: n
do n=nt1,nt2
   v1=Vx(n); v2=Vy(n); v3=Vz(n)
   Vx(n)=v1*rot_mat(1,1)+v2*rot_mat(2,1)+v3*rot_mat(3,1)
   Vy(n)=v1*rot_mat(1,2)+v2*rot_mat(2,2)+v3*rot_mat(3,2)
   Vz(n)=v1*rot_mat(1,3)+v2*rot_mat(2,3)+v3*rot_mat(3,3)
end do
end subroutine syn_rotate

subroutine scatter_rotate(VX,Vy,Vz,rot_mat,nt1,nt2,cblk)
real(SP),dimension(:,:,:,:),intent(inout) :: Vx,Vy,Vz
real(SP),dimension(:,:),intent(in) :: rot_mat
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt1,nt2
real(SP) :: v1,v2,v3
integer :: n,i,j,k
do n=nt1,nt2
do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)
   v1=Vx(i,j,k,n); v2=Vy(i,j,k,n); v3=Vz(i,j,k,n)
   Vx(i,j,k,n)=v1*rot_mat(1,1)+v2*rot_mat(2,1)+v3*rot_mat(3,1)
   Vy(i,j,k,n)=v1*rot_mat(1,2)+v2*rot_mat(2,2)+v3*rot_mat(3,2)
   Vz(i,j,k,n)=v1*rot_mat(1,3)+v2*rot_mat(2,3)+v3*rot_mat(3,3)
end do
end do
end do
end do
end subroutine scatter_rotate

subroutine scatter2kernel(V,U,Sa,Sb, &
           nt1,nt2,dt,cblk,Kap,Kaq,Kbp,Kbq)

real(SP),dimension(:),intent(in) :: V,U
real(SP),dimension(:,:,:,:),intent(in) :: Sa,Sb
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt1,nt2
real(SP),intent(in) :: dt

real(SP),dimension(:,:,:),intent(out) :: Kap,Kaq,Kbp,Kbq

real(SP) :: V0,U0
integer :: n,i,j,k

! normalization
!V0=(0.5*V(nt1)**2+0.5*V(nt2)**2)*kdt; U0=(0.5*U(nt1)**2+0.5*U(nt2)**2)*kdt
!do m=nt1+1,nt2-1
!   V0=V0+V(m)**2.0*kdt; U0=U0+U(m)**2.0*kdt
!end do
!V=V/V0/SGT_M0; U=U/U0/SGT_M0

! normalization
V0=(0.5*V(nt1)**2+0.5*V(nt2)**2)*kdt; U0=(0.5*U(nt1)**2+0.5*U(nt2)**2)*kdt
do n=nt1+1,nt2-1
   V0=V0+V(n)**2.0*dt; U0=U0+U(n)**2.0*dt
end do
!V=V/V0/SGT_M0; U=U/U0/SGT_M0
V0=V0*SGT_M0; U0=U0*SGT_M0

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   Ka(nt1:nt2)=Sa(i,j,k,nt1:nt2);
   Kb(nt1:nt2)=Sb(i,j,k,nt1:nt2);

   ! Forming kernels for tau_p and tau_q.
   Kap(i,j,k)=0.5*dt*(V(nt1)*Ka(nt1)+V(nt2)*Ka(nt2))
   Kaq(i,j,k)=0.5*dt*(U(nt1)*Ka(nt1)+U(nt2)*Ka(nt2))
   Kbp(i,j,k)=0.5*dt*(V(nt1)*Kb(nt1)+V(nt2)*Kb(nt2))
   Kbq(i,j,k)=0.5*dt*(U(nt1)*Kb(nt1)+U(nt2)*Kb(nt2))
   do n=nt1+1,nt2-1
      Kap(i,j,k)=Kap(i,j,k)+V(n)*Ka(n)*dt
      Kaq(i,j,k)=Kaq(i,j,k)+U(n)*Ka(n)*dt
      Kbp(i,j,k)=Kbp(i,j,k)+V(n)*Kb(n)*dt
      Kbq(i,j,k)=Kbq(i,j,k)+U(n)*Kb(n)*dt
   end do
   
   Kap(i,j,k)=Kap(i,j,k)/V0
   Kaq(i,j,k)=Kaq(i,j,k)/U0
   Kbp(i,j,k)=Kbp(i,j,k)/V0
   Kbq(i,j,k)=Kbq(i,j,k)/U0

end do
end do
end do
end subroutine scatter2kernel

subroutine scatter_filtfilt(w0,w,bval,aval,nt,cblk)
real(SP),dimension(:,:,:,:),intent(in) :: w0
real(SP),dimension(:,:,:,:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt
integer :: i,j,k,n

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   vecx(NFACT+1:NFACT+nt)=w0(i,j,k,1:nt)
   do n=1,NFACT
      vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   call filtfilt(NOD,NEL,bval,aval,vecx,vecy,vecw)
   w(i,j,k,1:nt)=vecy(NFACT+1:NFACT+nt)
   
end do
end do
end do
end subroutine scatter_filtfilt

subroutine scatter_filter(w0,w,bval,aval,nt,cblk)
real(SP),dimension(:,:,:,:),intent(in) :: w0
real(SP),dimension(:,:,:,:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,dimension(SEIS_GEO),intent(in) :: cblk
integer,intent(in) :: nt
integer :: i,j,k
integer :: n

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   vecx(NFACT+1:NFACT+nt)=w0(i,j,k,1:nt)
   do n=1,NFACT
      !vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(n)=0.0
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   call filter(NOD,NEL,bval,aval,vecx,vecy)
   w(i,j,k,1:nt)=vecy(NFACT+1:NFACT+nt)
   
end do
end do
end do
end subroutine scatter_filter

subroutine syn_filtfilt(w0,w,bval,aval,nt)
real(SP),dimension(:),intent(in) :: w0
real(SP),dimension(:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,intent(in) :: nt
integer :: n

   vecx(NFACT+1:NFACT+nt)=w0(1:nt)
   do n=1,NFACT
      vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   !do n=1,NFACT
   !   vecx(n)=2*w0(1)-w0(NFACT+1-(n-1))
   !end do
   !do n=NFACT+1,NFACT+nt
   !   vecx(n)=w0(n-NFACT)
   !end do
   !do n=NFACT+nt+1,NEL
   !   vecx(n)=2*w0(nt)-w0(nt-1-(n-(NFACT+nt+1)))
   !end do
   call filtfilt(NOD,NEL,bval,aval,vecx,vecy,vecw)
   w(1:nt)=vecy(NFACT+1:NFACT+nt)

end subroutine syn_filtfilt

subroutine syn_filter(w0,w,bval,aval,nt)
real(SP),dimension(:),intent(in) :: w0
real(SP),dimension(:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,intent(in) :: nt
integer :: n

   vecx(NFACT+1:NFACT+nt)=w0(1:nt)
   !do n=1,NFACT
   !   vecx(n)=2*vecx(1)-vecx(NFACT+1-(n-1))
   !   vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   !end do
   do n=1,NFACT
      !vecx(n)=2*vecx(NFACT+1)-vecx(2*NFACT+1-(n-1))
      vecx(n)=0.0
      vecx(NEL-NFACT+n)=2*vecx(NEL-NFACT)-vecx(NEL-NFACT-n)
   end do
   call filter(NOD,NEL,bval,aval,vecx,vecy)
   w(1:nt)=vecy(NFACT+1:NFACT+nt)

end subroutine syn_filter

subroutine filtfilt(NOD,NEL,b,a,x,y,w)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real(DP),dimension(:) :: x,y,w

call filter(NOD,NEL,b,a,x,y)
w=y(NEL:1:-1)
call filter(NOD,NEL,b,a,w,y)
y=y(NEL:1:-1)
end subroutine filtfilt

subroutine filter(NOD,NEL,b,a,x,y)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real(DP),dimension(:) :: x,y
integer :: n,m

do n=1,NEL
   y(n)=b(1)*x(n)
   do m=2,min(NOD+1,n)
      y(n)=y(n)+b(m)*x(n-m+1)-a(m)*y(n-m+1)
   end do
   y(n)=y(n)/a(1)
end do
end subroutine filter

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef KernelMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef KernelMPI
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

subroutine warning_print(msg)
  character (len=*),intent(in) :: msg
  print *, idstr//" warning :"//trim(msg)
end subroutine warning_print

end program SI_ker_sta

