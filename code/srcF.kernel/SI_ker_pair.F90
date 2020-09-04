program SI_ker_pair

! This program calculates finite-frequency sensitivity kernels for pair of
! event and station.
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
  SIG_FILTER_NONE=130, &
  SIG_FILTER_EMPT=200

character (len=SEIS_STRLEN) ::  &
    fnm_main_conf,fnm_ker_conf, &
    pnm_wave,pnm_sgt,           &
    fnm_syn,fnm_ker
character (len=SEIS_STRLEN) :: filenm
character (len=SEIS_STRLEN),dimension(:,:,:),allocatable :: pnm_ker

integer,dimension(SEIS_GEO) :: blksiz,blknum
integer :: kid,knt,nt1,nt2
real(SP) :: kdt

integer :: num_freq,num_filt,num_twin
integer :: nfq,npt,ntw
logical,dimension(:,:,:),allocatable :: measured
integer,allocatable ::  &
  filt_type(:,:),       &
  twinpix(:,:),         &
  twinidx(:,:,:,:)
real(SP),allocatable :: &
  twinlim(:,:,:,:)

integer :: NOD,NOD1,NEL
real(DP),dimension(:,:),allocatable :: filt_a,filt_b
real(DP),dimension(:),allocatable :: vecx,vecy,vecw

real(SP),dimension(:),allocatable :: V,U,T
real(SP),dimension(:),allocatable :: VV,UU
real(SP) :: sgt_m0

real(SP),dimension(:,:,:,:),allocatable :: &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS ,        &
    ExxR,EyyR,EzzR,ExyR,ExzR,EyzR

real(SP),dimension(:),allocatable :: Ka,Kb
real(SP),dimension(:),allocatable :: E11R,E22R,E33R,E12R,E13R,E23R
real(SP),dimension(:),allocatable :: E11S,E22S,E33S,E12S,E13S,E23S
real(SP),dimension(:,:,:),allocatable :: Kap,Kaq,Kbp,Kbq

integer :: n_i,n_j,n_k
integer :: i,j,k,m,n,mt,ierr
integer,dimension(SEIS_GEO) :: &
   bsubs,bsubc,bsubt,          &
    subs, subc, subt
!#ifdef KernelMPI
!integer,dimension(MPI_STATUS_SIZE) :: istatus
!#endif

integer :: ncidS,ncidF
integer :: TxxSid,TyySid,TzzSid,TxySid,TxzSid,TyzSid
integer :: TxxFid,TyyFid,TzzFid,TxyFid,TxzFid,TyzFid

integer,dimension(:,:,:),allocatable :: ncid
integer,dimension(:,:,:),allocatable :: kapid,kaqid,kbpid,kbqid

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
#else
call swmpi_set_gindx(0,0,0)
#endif

call media_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

!----------------------------------------------------------------------
! seismo on receiver
!#ifdef KernelMPI
!if (seismo_on_this(pnm_obsinfo,id_of_obs,indx_in_obs, &
!                   thisid(1),thisid(2),thisid(3),n)) then
!   filenm=get_fnm_seismo(pnm_obs,thisid(1),thisid(2),thisid(3))
!   i=snap_tinv(kid)/pt_tinv; ! recv type ouput may be denser than snap
!   call retrieve_recvline(filenm,n,varnm_obs,Vz,i,knt,i)
!   do n=0,dims(1)*dims(2)*dims(3)-1
!   if (n/=myid) then
!      call MPI_SEND(Vz,knt,MPI_REAL,n,3,MPI_COMM_WORLD,ierr)
!   end if
!   end do
!else
!   call MPI_RECV(Vz,knt,MPI_REAL,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,istatus,ierr)
!end if
!#endif

call init_synthetic(fnm_syn)
call init_filter(fnm_ker_conf)

kdt=snap_tinv(kid)*stept
if (abs(T(2)-T(1)-kdt)>SEIS_ZERO) then
   print *, kdt,T(2),T(1)
   print *, T(2)-T(1),abs(T(2)-T(1))
   !call error_except('dt /= stept*tinv')
end if

call twinlim2num

#ifdef KernelMPI
if (myid==0) &
#endif
   call export_info

!-----------------------------------------------------------------------------
#ifndef KernelMPI
print *, 'input mpi id:'
read *, n_i,n_j,n_k
call swmpi_change_fnm(n_i,n_j,n_k)
call swmpi_set_gindx(n_i,n_j,n_k)
thisid=(/ n_i,n_j,n_k /)
#else
n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#endif

call io_snap_locate(n_i,n_j,n_k)

kernode_if : if (snap_ishere(kid) .and. sgt_out(kid)) then

where (blksiz==-1)
   blksiz=snap_subc(:,kid)  
end where
blknum=(snap_subc(:,kid)+blksiz-1)/blksiz
call alloc_ker_var(blksiz(1),blksiz(2),blksiz(3),knt)
call alloc_media_local(blksiz(1),blksiz(2),blksiz(3))

! create kernel nc file
call alloc_nc_var(num_freq,num_filt,num_twin)
do nfq=1,num_freq
do npt=1,num_filt
do ntw=1,twinpix(npt,nfq)
   if (.not. measured(ntw,npt,nfq)) cycle
   fnm_ker=get_fnm_snapnode_n('./','kernel_',kid,0,n_i,n_j,n_k)
   call nfseis_grid3d_def(trim(pnm_ker(ntw,npt,nfq))//"/"//trim(fnm_ker),   &
     snap_subc(1,kid),snap_subc(2,kid), snap_subc(3,kid),ncid(ntw,npt,nfq), &
     "Finite frequency kernel")
   call nfseis_grid3d_defvar(ncid(ntw,npt,nfq),'kernel_phase_Vp'    ,kapid(ntw,npt,nfq))
   call nfseis_grid3d_defvar(ncid(ntw,npt,nfq),'kernel_amplitude_Vp',kaqid(ntw,npt,nfq))
   call nfseis_grid3d_defvar(ncid(ntw,npt,nfq),'kernel_phase_Vs'    ,kbpid(ntw,npt,nfq))
   call nfseis_grid3d_defvar(ncid(ntw,npt,nfq),'kernel_amplitude_Vs',kbqid(ntw,npt,nfq))
   call nfseis_grid3d_enddef(ncid(ntw,npt,nfq))
end do
end do
end do

! loop each block
!-----------------------------------------------------------------------------
do k=1,blknum(3)
do j=1,blknum(2)
do i=1,blknum(1)

#ifdef VERBOSE
   write(*,"(a7,3(i4,a1,i4.4),a9,3(i2.2))") &
       ' block:',                           &
       i,'/',blknum(1),                     &
       j,'/',blknum(2),                     &
       k,'/',blknum(3),                     &
       ', thisid=',n_i,n_j,n_k
#endif

   bsubs=(/ (i-1)*blksiz(1)+1,(j-1)*blksiz(2)+1,(k-1)*blksiz(3)+1 /)
   bsubc=(/ min(blksiz(1),snap_subc(1,kid)-bsubs(1)+1),  &
            min(blksiz(2),snap_subc(2,kid)-bsubs(2)+1),  &
            min(blksiz(3),snap_subc(3,kid)-bsubs(3)+1) /)
   bsubt=(/ 1,1,1 /)
   subc=bsubc; subt=snap_subt(:,kid)*bsubt
   subs=snap_subs(:,kid)+(bsubs-1)*subt

if (product(bsubc)*knt/=size(ExxS)) then
   call realloc_ker_var(bsubc(1),bsubc(2),bsubc(3),knt)
   call realloc_media_local(bsubc(1),bsubc(2),bsubc(3))
end if

! load modeling field
if (snap_tcnt(kid)<knt) then
   m=0; n=0; mt=snap_tcnt(kid)
   do 
      if (m>=knt) exit
      n=n+1
      filenm=get_fnm_snapnode_n(pnm_wave,'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidS)
      call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
      call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
      call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
      call nfseis_inq_varid(ncidS, 'Txy', TxySid)
      call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
      call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)

      filenm=get_fnm_snapnode_n(pnm_sgt,'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
      if (m+mt>knt) mt=knt-m

      ierr=nf90_get_var(ncidS,TxxSid,ExxS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyySid,EyyS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TzzSid,EzzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxySid,ExyS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxzSid,ExzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyzSid,EyzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))

      ierr=nf90_get_var(ncidF,TxxFid,ExxR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyyFid,EyyR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TzzFid,EzzR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxyFid,ExyR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxzFid,ExzR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyzFid,EyzR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel_mpi: get equavilent stress fail')
      m=m+mt
      call nfseis_close(ncidS)
      call nfseis_close(ncidF)
   end do
else
   if (i==1 .and. j==1 .and. k==1) then
      filenm=get_fnm_snapnode_n(pnm_wave,'sgt_',kid,1,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidS)
      call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
      call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
      call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
      call nfseis_inq_varid(ncidS, 'Txy', TxySid)
      call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
      call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)
   
      filenm=get_fnm_snapnode_n(pnm_sgt,'sgt_',kid,1,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
   end if
   ierr=nf90_get_var(ncidS,TxxSid,ExxS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyySid,EyyS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TzzSid,EzzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxySid,ExyS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxzSid,ExzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyzSid,EyzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxxFid,ExxR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyyFid,EyyR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TzzFid,EzzR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxyFid,ExyR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxzFid,ExzR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyzFid,EyzR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
end if

! load media
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt) 

! convert to stain
  call stress2strain(ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,lambda,mu,bsubc)
  call stress2strain(ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,lambda,mu,bsubc)

! cal scatter wave
  call sgt2scatter(                   &
       ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
       ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
       lambda,mu,1,knt,kdt,bsubc)

! cal kernel
do nfq=1,num_freq
do npt=1,num_filt

if (.not. measured(1,npt,nfq)) cycle

select case (filt_type(npt,nfq))
case (SIG_FILTER_P1)
  call scatter_filter(ExxS,ExxR,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
  call scatter_filter(EyyS,EyyR,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
  call syn_filter(V,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
  call syn_filter(U,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
case (SIG_FILTER_P2)
  call scatter_filtfilt(ExxS,ExxR,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
  call scatter_filtfilt(EyyS,EyyR,filt_b(:,nfq),filt_a(:,nfq),knt,bsubc)
  call syn_filtfilt(V,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
  call syn_filtfilt(U,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
case (SIG_FILTER_NONE)
  ExxR=ExxS;EyyR=EyyS;
  VV=V; UU=U;
case default
  print *, "npt,nfq,filt_type=",npt,nfq,filt_type(npt,nfq)
  call error_except('filt_type error')
end select

do ntw=1,twinpix(npt,nfq)
   nt1=twinidx(1,ntw,npt,nfq);nt2=twinidx(2,ntw,npt,nfq)
   call scatter2kernel(VV,UU,ExxR,EyyR,nt1,nt2,kdt,bsubc,Kap,Kaq,Kbp,Kbq)
   ! put
   call nfseis_put(ncid(ntw,npt,nfq),kapid(ntw,npt,nfq),Kap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(ntw,npt,nfq),kaqid(ntw,npt,nfq),Kaq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(ntw,npt,nfq),kbpid(ntw,npt,nfq),Kbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid(ntw,npt,nfq),kbqid(ntw,npt,nfq),Kbq,bsubs,bsubc,bsubt)
end do ! ntw

end do ! npt
end do ! nfq

end do
end do
end do ! block

!-----------------------------------------------------------------------------
if (snap_tcnt(kid)>=knt) then
   call nfseis_close(ncidS)
   call nfseis_close(ncidF)
end if

do nfq=1,num_freq
do npt=1,num_filt
do ntw=1,twinpix(npt,nfq)
   call nfseis_close(ncid(ntw,npt,nfq))
end do
end do
end do

end if kernode_if

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

subroutine alloc_number_var(nf,np,nw)
  integer,intent(in) :: nf,np,nw
  !allocate(ptsamp(nf)); ptsamp=1
  !allocate(ptsamp_type(nw,np,nf)); ptsamp_type=SIG_FILTER_EMPT
  allocate(measured(nw,np,nf)); measured=.false.
  allocate(filt_type(np,nf)); filt_type=SIG_FILTER_EMPT
  allocate(twinpix(np,nf)); twinpix=0
  allocate(twinlim(2,nw,np,nf)); twinlim=0.0
  allocate(twinidx(2,nw,np,nf)); twinidx=0
  allocate(pnm_ker(nw,np,nf))
end subroutine alloc_number_var

subroutine alloc_filt_var(nf,NOD1,NEL)
  integer,intent(in) :: nf,NOD1,NEL
  allocate(filt_b(NOD1,nf)); filt_b=0.0_DP
  allocate(filt_a(NOD1,nf)); filt_a=0.0_DP
  allocate(vecx(NEL)); vecx=0.0_DP
  allocate(vecy(NEL)); vecy=0.0_DP
  allocate(vecw(NEL)); vecw=0.0_DP
end subroutine alloc_filt_var

subroutine alloc_nc_var(nf,np,nw)
  integer,intent(in) :: nf,np,nw
  allocate(ncid(nw,np,nf)); ncid=0
  allocate(kapid(nw,np,nf)); kapid=0
  allocate(kaqid(nw,np,nf)); kaqid=0
  allocate(kbpid(nw,np,nf)); kbpid=0
  allocate(kbqid(nw,np,nf)); kbqid=0
end subroutine alloc_nc_var

subroutine alloc_seismo_var(nt)
  integer,intent(in) :: nt
  allocate(V(nt)); V=0.0
  allocate(U(nt)); U=0.0
  allocate(T(nt)); T=0.0
  allocate(VV(nt)); VV=0.0
  allocate(UU(nt)); UU=0.0
end subroutine alloc_seismo_var

subroutine alloc_ker_var(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0
  allocate(ExxR(ki,kj,kk,kt)); ExxR=0.0
  allocate(EyyR(ki,kj,kk,kt)); EyyR=0.0
  allocate(EzzR(ki,kj,kk,kt)); EzzR=0.0
  allocate(ExyR(ki,kj,kk,kt)); ExyR=0.0
  allocate(ExzR(ki,kj,kk,kt)); ExzR=0.0
  allocate(EyzR(ki,kj,kk,kt)); EyzR=0.0

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
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxR)) deallocate(ExxR)
  if (allocated(EyyR)) deallocate(EyyR)
  if (allocated(EzzR)) deallocate(EzzR)
  if (allocated(ExyR)) deallocate(ExyR)
  if (allocated(ExzR)) deallocate(ExzR)
  if (allocated(EyzR)) deallocate(EyzR)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)

  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0
  allocate(ExxR(ki,kj,kk,kt)); ExxR=0.0
  allocate(EyyR(ki,kj,kk,kt)); EyyR=0.0
  allocate(EzzR(ki,kj,kk,kt)); EzzR=0.0
  allocate(ExyR(ki,kj,kk,kt)); ExyR=0.0
  allocate(ExzR(ki,kj,kk,kt)); ExzR=0.0
  allocate(EyzR(ki,kj,kk,kt)); EyzR=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0 
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0 
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0 
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0 
end subroutine realloc_ker_var

subroutine alloc_media_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  allocate(lambda(ki,kj,kk)); lambda=0.0
  allocate(mu(ki,kj,kk)); mu=0.0
end subroutine alloc_media_local
subroutine realloc_media_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  if (allocated(lambda)) deallocate(lambda)
  if (allocated(mu)) deallocate(mu)
  allocate(lambda(ki,kj,kk)); lambda=0.0
  allocate(mu(ki,kj,kk)); mu=0.0
end subroutine realloc_media_local

subroutine dealloc_all
  if (allocated(V)) deallocate(V)
  if (allocated(U)) deallocate(U)
  if (allocated(T)) deallocate(T)
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxR)) deallocate(ExxR)
  if (allocated(EyyR)) deallocate(EyyR)
  if (allocated(EzzR)) deallocate(EzzR)
  if (allocated(ExyR)) deallocate(ExyR)
  if (allocated(ExzR)) deallocate(ExzR)
  if (allocated(EyzR)) deallocate(EyzR)
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
integer :: fid,n,num_meas
integer :: nfq,npt,ntw
real(SP) :: t1,t2
character (len=SEIS_STRLEN) :: str
integer,dimension(:,:),allocatable :: pix

fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)

call string_conf(fid,1,'SGT_ROOT',2,pnm_sgt)
call string_conf(fid,1,"SGT_M0",2,sgt_m0)
call string_conf(fid,1,'WAVE_ROOT',2,pnm_wave)

call string_conf(fid,1,'snap_id',2,kid)
do n=1,SEIS_GEO
 call string_conf(fid,1,'block_size',n+1,blksiz(n))
end do

call string_conf(fid,1,'fnm_syn',2,fnm_syn)

! numbers
call string_conf(fid,1,'number_of_freqband',2,num_freq)
call string_conf(fid,1,'number_of_measured',2,num_meas)
!call string_conf(fid,1,'number_of_filttype',2,num_filt)
!call string_conf(fid,1,'number_of_timewind',2,num_twin)
num_filt=3

! search num_twin
call string_conf(fid,1,'<time_window>',1,str)
num_twin=0
allocate(pix(num_filt,num_freq)); pix=0
do n=1,num_meas
   read(fid,*) nfq,npt
   read(fid,"(a132)") str
   pix(npt,nfq)=pix(npt,nfq)+1
   num_twin=max(num_twin,pix(npt,nfq))
end do
deallocate(pix)

call alloc_number_var(num_freq,num_filt,num_twin)

call string_conf(fid,1,'<time_window>',1,str)
do n=1,num_meas
   !read(fid,*) nfq,npt,t1,t2,str
   read(fid,*) nfq,npt,t1,t2
   read(fid,"(a132)") str
   ntw=twinpix(npt,nfq)+1
   twinpix(npt,nfq)=ntw
   measured(ntw,npt,nfq)=.true.
   pnm_ker(ntw,npt,nfq)=str
   twinlim(:,ntw,npt,nfq)=(/ t1,t2 /)
end do
filt_type(1,:)=SIG_FILTER_P1
filt_type(2,:)=SIG_FILTER_P2
filt_type(3,:)=SIG_FILTER_NONE

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
      NOD=NOD1-1; NEL=knt+2*NOD1+NOD
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

subroutine init_synthetic(fnm_syn)
character (len=*),intent(in) :: fnm_syn
integer :: fid,n

fid=5001
open(fid,file=trim(fnm_syn),status='old')
  read(fid,*) knt
  call alloc_seismo_var(knt)
  do n=1,knt
     read(fid,*) T(n),V(n),U(n)
  end do
close(fid)
end subroutine init_synthetic

subroutine twinlim2num
integer :: p(1)
integer :: nfq,npt,ntw

do nfq=1,num_freq
do npt=1,num_filt
do ntw=1,twinpix(npt,nfq)
   if (.not. measured(ntw,npt,nfq)) cycle
   p=maxloc(T,T<=twinlim(1,ntw,npt,nfq));
   twinidx(1,ntw,npt,nfq)=p(1)
   twinlim(1,ntw,npt,nfq)=T(p(1))

   p=minloc(T,T>=twinlim(2,ntw,npt,nfq));
   twinidx(2,ntw,npt,nfq)=p(1)
   twinlim(2,ntw,npt,nfq)=T(p(1))
end do
end do
end do
end subroutine twinlim2num

subroutine export_info
character (len=SEIS_STRLEN) :: filenm
integer :: nfq,npt,ntw,n,fid
fid=6001
do nfq=1,num_freq
do npt=1,num_filt
   if (.not. measured(1,npt,nfq)) cycle

   do ntw=1,twinpix(npt,nfq)
      filenm=trim(pnm_ker(ntw,npt,nfq))//'/kernel.para.log'
      open(fid,file=trim(filenm),status='unknown')
      write(fid,*) twinidx(:,ntw,npt,nfq),twinlim(:,ntw,npt,nfq)
      if (filt_type(npt,nfq)/=SIG_FILTER_NONE) then
         write(fid,"(f20.16)") filt_a(:,nfq)
         write(fid,"(f20.16)") filt_b(:,nfq)
      end if
      close(fid)
   end do

   if (filt_type(npt,nfq)==SIG_FILTER_NONE) cycle

   call syn_filtfilt(V,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
   call syn_filtfilt(U,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
   filenm=trim(pnm_ker(1,npt,nfq))//'/synthetic.filtfilt.dat'
   open(fid,file=trim(filenm),status='unknown')
   do n=1,knt
      write(fid,*) T(n),VV(n),UU(n)
   end do
   close(fid)

   call syn_filter(V,VV,filt_b(:,nfq),filt_a(:,nfq),knt)
   call syn_filter(U,UU,filt_b(:,nfq),filt_a(:,nfq),knt)
   filenm=trim(pnm_ker(1,npt,nfq))//'/synthetic.filter.dat'
   open(fid,file=trim(filenm),status='unknown')
   do n=1,knt
      write(fid,*) T(n),VV(n),UU(n)
   end do
   close(fid)

end do
end do
end subroutine export_info

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
!V=V/V0/sgt_m0; U=U/U0/sgt_m0

! normalization
V0=(0.5*V(nt1)**2+0.5*V(nt2)**2)*kdt; U0=(0.5*U(nt1)**2+0.5*U(nt2)**2)*kdt
do n=nt1+1,nt2-1
   V0=V0+V(n)**2.0*dt; U0=U0+U(n)**2.0*dt
end do
!V=V/V0/sgt_m0; U=U/U0/sgt_m0
V0=V0*sgt_m0; U0=U0*sgt_m0

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
integer :: i,j,k

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   vecx(NOD+1:NOD+nt)=w0(i,j,k,1:nt)
   call filtfilt(NOD,NEL,bval,aval,vecx,vecy,vecw)
   w(i,j,k,1:nt)=vecy(NOD+1:NOD+nt)
   
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

do k=1,cblk(3)
do j=1,cblk(2)
do i=1,cblk(1)

   vecx(NOD+1:NOD+nt)=w0(i,j,k,1:nt)
   call filter(NOD,NEL,bval,aval,vecx,vecy)
   w(i,j,k,1:nt)=vecy(NOD+1:NOD+nt)
   
end do
end do
end do
end subroutine scatter_filter

subroutine syn_filtfilt(w0,w,bval,aval,nt)
real(SP),dimension(:),intent(in) :: w0
real(SP),dimension(:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,intent(in) :: nt

   vecx(NOD+1:NOD+nt)=w0(1:nt)
   call filtfilt(NOD,NEL,bval,aval,vecx,vecy,vecw)
   w(1:nt)=vecy(NOD+1:NOD+nt)

end subroutine syn_filtfilt

subroutine syn_filter(w0,w,bval,aval,nt)
real(SP),dimension(:),intent(in) :: w0
real(SP),dimension(:),intent(out) :: w
real(DP),dimension(:),intent(in) :: bval,aval
integer,intent(in) :: nt

   vecx(NOD+1:NOD+nt)=w0(1:nt)
   call filter(NOD,NEL,bval,aval,vecx,vecy)
   w(1:nt)=vecy(NOD+1:NOD+nt)

end subroutine syn_filter

subroutine filtfilt(NOD,NEL,b,a,x,y,w)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real(DP),dimension(:) :: x,y,w

x(1:NOD)=0.0; x(NEL-2*(NOD+1)+1:NEL)=0.0
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

y(1:NOD)=0.0
do n=NOD+1,NEL
   y(n)=b(1)*x(n)
   do m=2,NOD+1
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

end program SI_ker_pair

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
