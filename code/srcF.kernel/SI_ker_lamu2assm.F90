program SI_ker2assm

! This program is used to assemble kernel value on inversion block.
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
#define PhaseDelay
!#define AmplitudeReduction

use constants_mod
use string_mod
use para_mod
use io_mod
use media_mod
use nfseis_mod
use src_mod
use mpi_mod
#ifdef KernelAssmMPI
use mpi
#endif

implicit none

type SNAPBD
     integer :: n_i,n_j,n_k
     integer,dimension(SEIS_GEO) :: indxs,indxe,indxc
     integer,dimension(SEIS_GEO) :: subs,sube,subt,subc
end type SNAPBD

integer :: nthd
type(SNAPBD),allocatable,dimension(:) :: info_list

character (len=8) :: str_d0,str_d1,str_d2
character (len=10):: str_t0,str_t1,str_t2

character (len=300) :: rcdstr
character (len=SEIS_STRLEN) ::   &
    fnm_main_conf,fnm_ker_conf,  &
    evt_list,sgt_list,pnm_spool, &
    fnm_synx,fnm_syny,fnm_synz,  &
    pnm_ker,fnm_ker,fnm_invblk
character (len=SEIS_STRLEN) :: &
    KSTNM,KEVNM
character (len=SEIS_STRLEN) :: filenm,filenm_part
character (len=15) :: idstr

integer,dimension(SEIS_GEO) :: blksiz
integer :: kid

integer :: ncmp,npk
integer :: picknum(6)
integer :: pickfreq
real(SP),dimension(3,3) :: &
  rot_mat
character (len=2) :: filt_nm
real(SP) :: t1,t2

real(SP),dimension(:,:,:),allocatable :: &
    Kap,Kaq,Kbp,Kbq

integer :: cur_bdap,cur_bdaq,cur_bdbp,cur_bdbq
integer :: max_bdap,max_bdaq,max_bdbp,max_bdbq
integer :: num_blk,max_pix
real :: kap0,kaq0,kbp0,kbq0
real :: Cap0,Caq0,Cbp0,Cbq0
real(SP),dimension(:),allocatable :: &
    Bap,Baq,Bbp,Bbq,                 &
    Sap,Saq,Sbp,Sbq,                 &
    kerval
real(SP),dimension(:,:),allocatable :: volm_blk
integer,dimension(:),allocatable :: npix_blk,keridx
integer,dimension(:,:,:),allocatable :: pixs_blk
real(SP) :: V

integer :: n_i,n_j,n_k
integer :: i,j,k,m,n,i1,i2,j1,j2,k1,k2,ierr
integer,dimension(SEIS_GEO) :: bsubs,bsubc,bsubt

integer :: staid,evtid

!----------------------------------------------------------------------

#ifdef KernelAssmMPI
call MPI_INIT(ierr)
#endif

call get_conf_name(fnm_conf)

! read kernel conf
fnm_ker_conf='TomoKernel.conf'
call init_kernel(fnm_ker_conf)
fnm_conf=fnm_main_conf

call para_init(fnm_conf)
call swmpi_init(fnm_conf)

#ifdef KernelAssmMPI
  call swmpi_cart_creat
  call swmpi_reinit_para
  call swmpi_datatype
  if (masternode) then
     read(*,"(a)") sgt_list
     read(*,"(a)") fnm_invblk
     read(*,"(a)") pnm_spool
     read(*,*) Cap0,Caq0,Cbp0,Cbq0
     read(*,*) kap0,kaq0,kbp0,kbq0

     write(*,*) sgt_list
     write(*,*) fnm_invblk
     write(*,*) pnm_spool
     write(*,*) Cap0,Caq0,Cbp0,Cbq0
     write(*,*) kap0,kaq0,kbp0,kbq0
  end if
  call MPI_BCAST(sgt_list,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(fnm_invblk,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(pnm_spool,SEIS_STRLEN,MPI_CHARACTER,0,SWMPI_COMM,ierr)
  call MPI_BCAST(Cap0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(Caq0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(Cbp0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(Cbq0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kap0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kaq0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kbp0,1,MPI_REAL,0,SWMPI_COMM,ierr)
  call MPI_BCAST(kbq0,1,MPI_REAL,0,SWMPI_COMM,ierr)
#else
  call swmpi_set_gindx(0,0,0)
  read(*,"(a)") sgt_list
  read(*,"(a)") fnm_invblk
  read(*,"(a)") pnm_spool
  read(*,*) Cap0,Caq0,Cbp0,Cbq0
  read(*,*) kap0,kaq0,kbp0,kbq0
  write(*,*) sgt_list
  write(*,*) fnm_invblk
  write(*,*) pnm_spool
  write(*,*) Cap0,Caq0,Cbp0,Cbq0
  write(*,*) kap0,kaq0,kbp0,kbq0
#endif

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

blksiz=snap_subc(:,kid)
!bsubs=(/ 1,1,1 /); bsubc=blksiz; bsubt=(/ 1,1,1 /)

!-- read inversion cell information
call init_inv_block(fnm_invblk)
call alloc_local(blksiz(1),blksiz(2),blksiz(3))

!call io_snap_locate(n_i,n_j,n_k)
call locate_snap(kid,nthd)

#ifdef PhaseDelay
write(*,*) "assemble phase-delay kernel ..."
#endif
#ifdef AmplitudeReduction
write(*,*) "assemble amplitude-reduction kernel ..."
#endif

call date_and_time(date=str_d0,time=str_t0)
write(idstr,"(a8,3(i2.2))") ' thisid=',0,0,0
#ifdef VERBOSE
  write(*,*) ' begins from ',str_d0,  &
                  ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
#endif

staid=10001; evtid=10002
! loop all stations
!----------------------------------------------------------------------
open(staid,file=trim(sgt_list),status='old',iostat=ierr)
if (ierr>0) call error_except("sgt_list open err:"//trim(sgt_list))

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
  write(*,*) trim(KSTNM)//" "//trim(evt_list), &
                    ': begins from ',str_d1,  &
                  ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
#endif

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
   write(*,*) trim(rcdstr)
#endif

  n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
  KEVNM=rcdstr(1:n)
  rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  read(rcdstr,*) picknum, rot_mat
  read(evtid,"(a132)") fnm_synx
  read(evtid,"(a132)") fnm_syny
  read(evtid,"(a132)") fnm_synz

  if (all(picknum<=0)) cycle

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
   write(*,*) trim(rcdstr)
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

   ! read kernel
   do n=1,nthd
      n_i=info_list(n)%n_i;n_j=info_list(n)%n_j;n_k=info_list(n)%n_k
      call swmpi_change_fnm(n_i,n_j,n_k)
      call swmpi_set_gindx(n_i,n_j,n_k)
   
      i1=info_list(n)%indxs(1);j1=info_list(n)%indxs(2);k1=info_list(n)%indxs(3);
      i2=info_list(n)%indxe(1);j2=info_list(n)%indxe(2);k2=info_list(n)%indxe(3);
      bsubc=info_list(n)%indxc; bsubs=(/1,1,1/); bsubt=(/1,1,1/)

      fnm_ker=get_fnm_snapnode_n('./','kernel_lamu_',kid,0,n_i,n_j,n_k)
      fnm_ker=trim(pnm_ker)//"/"//trim(fnm_ker)
#ifdef PhaseDelay
      call nfseis_varget(trim(fnm_ker),'phase_la',Kap(i1:i2,j1:j2,k1:k2),bsubs,bsubc,bsubt)
      call nfseis_varget(trim(fnm_ker),'phase_mu',Kbp(i1:i2,j1:j2,k1:k2),bsubs,bsubc,bsubt)
#endif
#ifdef AmplitudeReduction
      call nfseis_varget(trim(fnm_ker),'amplitude_la',Kaq(i1:i2,j1:j2,k1:k2),bsubs,bsubc,bsubt)
      call nfseis_varget(trim(fnm_ker),'amplitude_mu',Kbq(i1:i2,j1:j2,k1:k2),bsubs,bsubc,bsubt)
#endif
   end do

#ifdef PhaseDelay
   Bap=0.0; Bbp=0.0;
#endif
#ifdef AmplitudeReduction
   Baq=0.0; Bbq=0.0
#endif
   do n=1,num_blk
   do m=1,npix_blk(n)
      i=pixs_blk(1,m,n);j=pixs_blk(2,m,n);k=pixs_blk(3,m,n)
      V=volm_blk(m,n)
#ifdef PhaseDelay
      Bap(n)=Bap(n)+Kap(i,j,k)*V
      Bbp(n)=Bbp(n)+Kbp(i,j,k)*V
#endif
#ifdef AmplitudeReduction
      Baq(n)=Baq(n)+Kaq(i,j,k)*V
      Bbq(n)=Bbq(n)+Kbq(i,j,k)*V
#endif
   end do
#ifdef PhaseDelay
      Sap(n)=Sap(n)+abs(Bap(n))
      Sbp(n)=Sbp(n)+abs(Bbp(n))
#endif
#ifdef AmplitudeReduction
      Saq(n)=Saq(n)+abs(Baq(n))
      Sbq(n)=Sbq(n)+abs(Bbq(n))
#endif
   end do
   
   filenm_part=splash2underline(pnm_ker)

#ifdef PhaseDelay
   Bap=Bap*Cap0
   filenm=trim(pnm_spool)//'/'//trim(filenm_part)//'.Kap.unformat'
   call export_G_compact(Bap,Kap0,filenm,cur_bdap,max_bdap)

   Bbp=Bbp*Cbp0
   filenm=trim(pnm_spool)//'/'//trim(filenm_part)//'.Kbp.unformat'
   call export_G_compact(Bbp,Kbp0,filenm,cur_bdbp,max_bdbp)
#endif

#ifdef AmplitudeReduction
   Baq=Baq*Caq0
   filenm=trim(pnm_spool)//'/'//trim(filenm_part)//'.Kaq.unformat'
   call export_G_compact(Baq,Kaq0,filenm,cur_bdaq,max_bdaq)

   Bbq=Bbq*Cbq0
   filenm=trim(pnm_spool)//'/'//trim(filenm_part)//'.Kbq.unformat'
   call export_G_compact(Bbq,Kbq0,filenm,cur_bdbq,max_bdbq)
#endif

#ifdef VERBOSE
   !write(*,*) "max_band=",max_bdap,max_bdaq,max_bdbp,max_bdbq
   write(*,*) "G_band=",cur_bdap,cur_bdaq,cur_bdbp,cur_bdbq, &
              maxval( (/max_bdap,max_bdaq,max_bdbp,max_bdbq/) )
#endif

end do ! npk
end do ! ncmp

end do evt_loop
close(evtid)

end do sta_loop
close(staid)

call export_inv_info

#ifdef PhaseDelay
Sap=Sap*Cap0
filenm=trim(pnm_spool)//'/'//'inv_G_sum.Kap.unformat'
call export_G_sum(Sap,filenm)
Sbp=Sbp*Cbp0
filenm=trim(pnm_spool)//'/'//'inv_G_sum.Kbp.unformat'
call export_G_sum(Sbp,filenm)
#endif
#ifdef AmplitudeReduction
Saq=Saq*Cap0
filenm=trim(pnm_spool)//'/'//'inv_G_sum.Kaq.unformat'
call export_G_sum(Saq,filenm)
Sbq=Sbq*Cbp0
filenm=trim(pnm_spool)//'/'//'inv_G_sum.Kbq.unformat'
call export_G_sum(Sbq,filenm)
#endif

!-----------------------------------------------------------------------------

#ifdef VERBOSE
  call date_and_time(date=str_d2,time=str_t2)
  write(*,*) ' finish at ',str_d2,  &
                     ',',str_t2(1:2),'/',str_t2(3:4),'/',str_t2(5:10)
#endif

!-----------------------------------------------------------------------------
call dealloc_all

#ifdef KernelAssmMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine alloc_inv_block(nblk,npix)
  integer,intent(in) :: nblk,npix
  allocate(volm_blk(npix,nblk)); volm_blk=0.0
  allocate(npix_blk(nblk)); npix_blk=0
  allocate(pixs_blk(SEIS_GEO,npix,nblk)); pixs_blk=0
  allocate(keridx(nblk)); keridx=0
  allocate(kerval(nblk)); kerval=0.0
  allocate(Sap(nblk)); Sap=0.0
  allocate(Saq(nblk)); Saq=0.0
  allocate(Sbp(nblk)); Sbp=0.0
  allocate(Sbq(nblk)); Sbq=0.0
  allocate(Bap(nblk)); Bap=0.0
  allocate(Baq(nblk)); Baq=0.0
  allocate(Bbp(nblk)); Bbp=0.0
  allocate(Bbq(nblk)); Bbq=0.0
end subroutine alloc_inv_block

subroutine alloc_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  allocate(info_list(dims(1)*dims(2)*dims(3)))
  allocate(Kap(ki,kj,kk)); Kap=0.0
  allocate(Kaq(ki,kj,kk)); Kaq=0.0
  allocate(Kbp(ki,kj,kk)); Kbp=0.0
  allocate(Kbq(ki,kj,kk)); Kbq=0.0
end subroutine alloc_local

subroutine dealloc_all
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
  if (allocated(Bap)) deallocate(Bap)
  if (allocated(Baq)) deallocate(Baq)
  if (allocated(Bbp)) deallocate(Bbp)
  if (allocated(Bbq)) deallocate(Bbq)
  if (allocated(Sap)) deallocate(Sap)
  if (allocated(Saq)) deallocate(Saq)
  if (allocated(Sbp)) deallocate(Sbp)
  if (allocated(Sbq)) deallocate(Sbq)
  if (allocated(info_list)) deallocate(info_list)
  if (allocated(volm_blk)) deallocate(volm_blk)
  if (allocated(npix_blk)) deallocate(npix_blk)
  if (allocated(pixs_blk)) deallocate(pixs_blk)
end subroutine dealloc_all

!-----------------------------------------------------------------------------
subroutine init_kernel(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer :: fid

fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)
call string_conf(fid,1,'snap_id',2,kid)

!call string_conf(fid,1,'SGT_LIST',2,sgt_list)
close(fid)
end subroutine init_kernel

subroutine init_inv_block(fnm_invblk)
character (len=*),intent(in) :: fnm_invblk

call nfseis_diminfo(fnm_invblk,'block',num_blk)
call nfseis_diminfo(fnm_invblk,'max_pix',max_pix)

call alloc_inv_block(num_blk,max_pix)

call nfseis_varget(fnm_invblk,'volume',volm_blk, &
     (/1,1/),(/max_pix,num_blk/),(/1,1/))
call nfseis_varget(fnm_invblk,'num_of_pix',npix_blk, &
     (/1/),(/num_blk/),(/1/))
call nfseis_varget(fnm_invblk,'indx_of_pix',pixs_blk, &
     (/1,1,1/),(/SEIS_GEO,max_pix,num_blk/),(/1,1,1/))
! initialize
max_bdap=0
max_bdaq=0
max_bdbp=0
max_bdbq=0
end subroutine init_inv_block

subroutine locate_snap(id,nthd)
integer :: id,nthd
integer :: n_i1,n_i2,n_j1,n_j2,n_k1,n_k2
integer :: n_i,n_j,n_k
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

function splash2underline(pnm) result(fnm)
character (len=*),intent(in) :: pnm
character (len=SEIS_STRLEN) :: fnm
integer :: n

fnm=pnm
do n=1,len_trim(fnm)
   if (fnm(n:n)/="." .and. fnm(n:n)/='/') exit
   fnm(n:n)=" "
end do
fnm=adjustl(fnm)
do n=1,len_trim(fnm)
   if (fnm(n:n)=='/') fnm(n:n)="_"
end do
end function splash2underline

subroutine export_G_compact(kerblk,ker0,filenm,cur_bdwid,max_bdwid)
character (len=*),intent(in) :: filenm
real,dimension(:),intent(in) :: kerblk
real,intent(in) :: ker0
integer,intent(out) :: cur_bdwid
integer,intent(inout) :: max_bdwid

integer :: fid,m,n,ncoef

fid=1001
ncoef=0
do n=1,num_blk
   if (abs(kerblk(n))>ker0) then
      ncoef=ncoef+1
      keridx(ncoef)=n
      kerval(ncoef)=kerblk(n)
   end if
end do
cur_bdwid=ncoef
max_bdwid=max(max_bdwid,ncoef)

open(fid,file=trim(filenm),status='unknown',form='unformatted')
write(fid) num_blk,ker0
write(fid) ncoef
write(fid) (keridx(m),kerval(m),m=1,ncoef)
close(fid)

end subroutine export_G_compact

subroutine export_inv_info
integer :: fid
character (len=SEIS_STRLEN) :: filenm

fid=1002
filenm=trim(pnm_spool)//'/'//'inv_G_band.conf'
open(fid,file=trim(filenm),status='unknown')
write(fid,"(a,i)") "max_band_width_Kap = ",max_bdap
write(fid,"(a,i)") "max_band_width_Kaq = ",max_bdaq
write(fid,"(a,i)") "max_band_width_Kbp = ",max_bdbp
write(fid,"(a,i)") "max_band_width_Kbq = ",max_bdbq
close(fid)
end subroutine export_inv_info

subroutine export_G_sum(kerval,filenm)
character (len=*),intent(in) :: filenm
real,dimension(:),intent(in) :: kerval

integer :: fid,m

fid=1001
open(fid,file=trim(filenm),status='unknown',form='unformatted')
write(fid) num_blk
write(fid) (kerval(m),m=1,num_blk)
close(fid)
end subroutine export_G_sum

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef KernelAssmMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef KernelAssmMPI
  print *, thisid
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

end program SI_ker2assm

