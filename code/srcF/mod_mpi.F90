module mpi_mod

! This module is used for mpi parallel
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-18 14:50:53 -0500 (Sun, 18 Jan 2009) $
! $Revision: 74 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"
!#define WithoutFreeSurface

use mpi

use constants_mod
use para_mod
use string_mod

implicit none

private
public :: swmpi_init,                          &
          swmpi_change_fnm,                    &
          swmpi_cart_creat,                    &
          swmpi_reinit_para,                   &
          swmpi_set_gindx,                     &
          swmpi_datatype,                      &
          swmpi_rename_fnm,                    &
          swmpi_globi,swmpi_globj,swmpi_globk, &
          swmpi_locli,swmpi_loclj,swmpi_loclk, &
          point_in_thisnode,                   &
          swmpi_time_init,                     &
          swmpi_time_write,                    &
          swmpi_time_end,                      &
          set_mpi_subfix,                      &
          swmpi_except

integer,public :: ndims
integer,dimension(SEIS_GEO),public :: dims
integer,public :: SWMPI_COMM !communicator include topo information
integer,public :: myid
integer,dimension(SEIS_GEO),public :: thisid  !the thread coords in CART topo
integer,dimension(SEIS_GEO,2),public :: neigid

integer,public :: DTypeXS,DTypeXL, DTypeYS,DTypeYL, DTypeZS,DTypeZL

logical,public :: freenode,masternode
logical,dimension(SEIS_GEO,2),public :: absnode
character (len=SEIS_STRLEN),public :: fnm_swmpi

logical,dimension(SEIS_GEO) :: periods
logical :: reorder

integer,dimension(MPI_STATUS_SIZE) :: istatus
integer psize !the processor's size in MPI_COMM_WORLD
integer ierr
character (len=MPI_MAX_PROCESSOR_NAME) :: hostname !the host name
integer hostlen !the hostname's char length

character (len=8) :: str_d0,str_d1,str_d2
character (len=10):: str_t0,str_t1,str_t2
real (kind=8) :: wtime0,wtime1,wtime2
integer :: wtimen0

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine swmpi_init(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer fid,n
fid=1001
open(fid,file=trim(fnm_conf),status="old")
  !call string_conf(fid,1,'ndims',2,ndims)
  ndims=3
  do n=1,ndims
     call string_conf(fid,1,'dims',n+1,dims(n))
     !call string_conf(fid,1,'periods',n+1,periods(n))
  end do
  !call string_conf(fid,1,'reorder',2,reorder)
  periods=.false.
  reorder=.true.
close(fid)
end subroutine swmpi_init
#ifdef MPI1DTOPO
subroutine swmpi_cart_creat
   integer m,n
   SWMPI_COMM=MPI_COMM_WORLD
   call MPI_COMM_RANK(SWMPI_COMM, &
                      myid,       &
                      ierr)
   thisid(3)=myid/(dims(1)*dims(2))
   thisid(2)=myid/dims(1)
   thisid(1)=mod(myid,dims(1))

   if (thisid(1)==0) then
      neigid(1,1)=MPI_PROC_NULL
   else
      neigid(1,1)=myid-1
   end if
   if (thisid(1)==dims(1)-1) then
      neigid(1,2)=MPI_PROC_NULL
   else
      neigid(1,2)=myid+1
   end if

   if (thisid(2)==0) then
      neigid(2,1)=MPI_PROC_NULL
   else
      neigid(2,1)=myid-dims(1)
   end if
   if (thisid(2)==dims(2)-1) then
      neigid(2,2)=MPI_PROC_NULL
   else
      neigid(2,2)=myid+dims(1)
   end if

   if (thisid(3)==0) then
      neigid(3,1)=MPI_PROC_NULL
   else
      neigid(3,1)=myid-dims(1)*dims(2)
   end if
   if (thisid(3)==dims(3)-1) then
      neigid(3,2)=MPI_PROC_NULL
   else
      neigid(3,2)=myid+dims(1)*dims(2)
   end if

   call MPI_COMM_SIZE(SWMPI_COMM,      &
                      psize,           &
                      ierr)
   call MPI_GET_PROCESSOR_NAME(        &
                      hostname,        &
                      hostlen,         &
                      ierr)
   fnm_swmpi=set_mpi_prefix(thisid(1), &
                            thisid(2), &
                            thisid(3))
   ! swmpi010101_coordx.nc
   !set flags
   freenode=.false.
#ifndef WithoutFreeSurface
   if (neigid(3,2) == MPI_PROC_NULL) freenode=.true.
#endif

   absnode=.false.
   do n=1,SEIS_GEO
   do m=1,2
      if (neigid(n,m)==MPI_PROC_NULL) absnode(n,m)=.true.
   end do
   end do
   !absnode(SEIS_GEO,2)=.false.  ! z2 is inner or free

   masternode=.false.
   if (myid==0) masternode=.true.

end subroutine swmpi_cart_creat
#else
subroutine swmpi_cart_creat
   integer m,n
   call MPI_CART_CREATE(MPI_COMM_WORLD,         &
                        ndims,                  &
                        dims,                   &
                        periods,                &
                        reorder,                &
                        SWMPI_COMM,             &
                        ierr)
   call MPI_COMM_RANK(SWMPI_COMM,               &
                      myid,                     &
                      ierr)
   call MPI_CART_COORDS(SWMPI_COMM,             &
                        myid,                   &
                        ndims,                  &
                        thisid,                 &
                        ierr)
   call MPI_CART_SHIFT(SWMPI_COMM,              &
                       0,1,                     &
                       neigid(1,1),neigid(1,2), &
                       ierr)
   call MPI_CART_SHIFT(SWMPI_COMM,              &
                       1,1,                     &
                       neigid(2,1),neigid(2,2), &
                       ierr)
   call MPI_CART_SHIFT(SWMPI_COMM,              &
                       2,1,                     &
                       neigid(3,1),neigid(3,2), &
                       ierr)
   call MPI_COMM_SIZE(SWMPI_COMM,               &
                      psize,                    &
                      ierr)
   call MPI_GET_PROCESSOR_NAME(                 &
                      hostname,                 &
                      hostlen,                  &
                      ierr)
   fnm_swmpi=set_mpi_prefix(thisid(1),          &
                            thisid(2),          &
                            thisid(3))
   ! swmpi010101_coordx.nc
   !set flags
   freenode=.false.
#ifndef WithoutFreeSurface
   if (neigid(3,2) == MPI_PROC_NULL) freenode=.true.
#endif

   absnode=.false.
   do n=1,SEIS_GEO
   do m=1,2
      if (neigid(n,m)==MPI_PROC_NULL) absnode(n,m)=.true.
   end do
   end do
   !absnode(SEIS_GEO,2)=.false.  ! z2 is inner or free

   masternode=.false.
   if (myid==0) masternode=.true.

end subroutine swmpi_cart_creat
#endif
subroutine swmpi_reinit_para
   ! init ngi1 etc
   ngi1=swmpi_globi(ni1,thisid(1))
   ngi2=swmpi_globi(ni2,thisid(1))
   ngj1=swmpi_globj(nj1,thisid(2))
   ngj2=swmpi_globj(nj2,thisid(2))
   ngk1=swmpi_globk(nk1,thisid(3))
   ngk2=swmpi_globk(nk2,thisid(3))
   ngx1=swmpi_globi(nx1,thisid(1))
   ngx2=swmpi_globi(nx2,thisid(1))
   ngy1=swmpi_globj(ny1,thisid(2))
   ngy2=swmpi_globj(ny2,thisid(2))
   ngz1=swmpi_globk(nz1,thisid(3))
   ngz2=swmpi_globk(nz2,thisid(3))
   point_in_this=(/ ngi1,ngi2,ngj1,ngj2,ngk1,ngk2 /)
   if (thisid(1)==0)        point_in_this(1)=ngx1
   if (thisid(1)==dims(1)-1) point_in_this(2)=ngx2
   if (thisid(2)==0)        point_in_this(3)=ngy1
   if (thisid(2)==dims(2)-1) point_in_this(4)=ngy2
   if (thisid(3)==0)        point_in_this(5)=ngz1
   if (thisid(3)==dims(3)-1) point_in_this(6)=ngz2

   NTPI=ni*dims(1); NTPJ=nj*dims(2); NTPK=nk*dims(3)
   NTPX=NTPI+(nx-ni); NTPY=NTPJ+(ny-nj); NTPZ=NTPK+(nz-nk)
   npi1=thisid(1)*ni+1; npi2=(thisid(1)+1)*ni
   npj1=thisid(2)*nj+1; npj2=(thisid(2)+1)*nj
   npk1=thisid(3)*nk+1; npk2=(thisid(3)+1)*nk
end subroutine swmpi_reinit_para

subroutine swmpi_set_gindx(n_i,n_j,n_k)
   integer,intent(in) :: n_i,n_j,n_k
   ! init ngi1 etc
   ngi1=swmpi_globi(ni1,n_i)
   ngi2=swmpi_globi(ni2,n_i)
   ngj1=swmpi_globj(nj1,n_j)
   ngj2=swmpi_globj(nj2,n_j)
   ngk1=swmpi_globk(nk1,n_k)
   ngk2=swmpi_globk(nk2,n_k)
   ngx1=swmpi_globi(nx1,n_i)
   ngx2=swmpi_globi(nx2,n_i)
   ngy1=swmpi_globj(ny1,n_j)
   ngy2=swmpi_globj(ny2,n_j)
   ngz1=swmpi_globk(nz1,n_k)
   ngz2=swmpi_globk(nz2,n_k)
   point_in_this=(/ ngi1,ngi2,ngj1,ngj2,ngk1,ngk2 /)
   if (n_i==0)         point_in_this(1)=ngx1
   if (n_i==dims(1)-1) point_in_this(2)=ngx2
   if (n_j==0)         point_in_this(3)=ngy1
   if (n_j==dims(2)-1) point_in_this(4)=ngy2
   if (n_k==0)         point_in_this(5)=ngz1
   if (n_k==dims(3)-1) point_in_this(6)=ngz2
   NTPI=ni*dims(1); NTPJ=nj*dims(2); NTPK=nk*dims(3)
   NTPX=NTPI+(nx-ni); NTPY=NTPJ+(ny-nj); NTPZ=NTPK+(nz-nk)
   npi1=n_i*ni+1; npi2=(n_i+1)*ni
   npj1=n_j*nj+1; npj2=(n_j+1)*nj
   npk1=n_k*nk+1; npk2=(n_k+1)*nk
end subroutine swmpi_set_gindx

!*************************************************************************
!* nx,ny,nz, et al, should be inited before this subroutine              *
!*************************************************************************
subroutine swmpi_datatype
   !---- mpi type definition  ----
   call MPI_TYPE_VECTOR(ny*nz,            &
                        LenFDS,           &
                        nx,               &
                        SEISMPI_DATATYPE, &
                        DTypeXS,          &
                        ierr)
   call MPI_TYPE_VECTOR(ny*nz,            &
                        LenFDL,           &
                        nx,               &
                        SEISMPI_DATATYPE, &
                        DTypeXL,          &
                        ierr)
   call MPI_TYPE_VECTOR(nz,               &
                        nx*LenFDS,        &
                        nx*ny,            &
                        SEISMPI_DATATYPE, &
                        DTypeYS,          &
                        ierr)
   call MPI_TYPE_VECTOR(nz,               &
                        nx*LenFDL,        &
                        nx*ny,            &
                        SEISMPI_DATATYPE, &
                        DTypeYL,          &
                        ierr)
   call MPI_TYPE_VECTOR(LenFDS,           &
                        nx*ny,            &
                        nx*ny,            &
                        SEISMPI_DATATYPE, &
                        DTypeZS,          &
                        ierr)
   call MPI_TYPE_VECTOR(LenFDL,           &
                        nx*ny,            &
                        nx*ny,            &
                        SEISMPI_DATATYPE, &
                        DTypeZL,          &
                        ierr)
   call MPI_TYPE_COMMIT(DTypeXS,ierr)
   call MPI_TYPE_COMMIT(DTypeYS,ierr)
   call MPI_TYPE_COMMIT(DTypeZS,ierr)
   call MPI_TYPE_COMMIT(DTypeXL,ierr)
   call MPI_TYPE_COMMIT(DTypeYL,ierr)
   call MPI_TYPE_COMMIT(DTypeZL,ierr)
end subroutine swmpi_datatype

subroutine swmpi_time_init(filenm,ntime)
   character (len=*),intent(in) :: filenm
   integer,intent(in) :: ntime
   integer fid
   fid=1009
   if (myid/=0) return
   wtime0=MPI_WTIME()
   call date_and_time(date=str_d0,time=str_t0)

if (ntime==0) then
   open(fid,file=trim(filenm),status="unknown")
   write(fid,*) "# seis3d_wave run time log"
else
   open(fid,file=trim(filenm),status="old",position="append")
   write(fid,*)
   write(fid,*) "# seis3d_wave restart from ntime=",ntime
end if
     write(fid,*)
     write(fid,*) 'the program begins from ',str_d0,  &
                  ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
     write(fid,*)
     write(fid,*)  'each step uses time'
     !write(fid,'(4a10)') 'step','hour','minute','second'
   close(fid)
   wtime1=wtime0
   wtimen0=ntime
end subroutine swmpi_time_init
subroutine swmpi_time_write(ntime,filenm)
   character (len=*),intent(in) :: filenm
   integer,intent(in) :: ntime
   integer :: fid
   !integer h,m
   real (kind=8) :: s,s0,s1
   if (myid/=0) return
   fid=1009
   wtime2=MPI_WTIME()
   s=wtime2-wtime1
   !h=int(s/3600)
   !m=int((s-h*3600)/60)
   !s=s-h*3600-m*60
   s0=(wtime2-wtime0)/3600.0
   s1=s0/(ntime-wtimen0)*(nt-wtimen0)
   open(fid,file=trim(filenm),status="old",position="append")
     write(fid,'(i6,a10,g12.5,a4,g12.5,a17,g12.5,a11)') &
             ntime,'step uses',s,'s,  ',                &
             s0,' hours passed and', s1-s0, ' hours left'
   !  write(fid,'(i6,a10,i5,a1,i4,a1,f7.3,a1)')   &
   !          ntime,'time is',h,'h',m,'m',s,'s'
   close(fid)

   wtime1=wtime2
end subroutine swmpi_time_write
subroutine swmpi_time_end(filenm)
   character (len=*),intent(in) :: filenm
   integer fid
   integer d,h,m
   real (kind=8) :: s
   if (myid/=0) return
   fid=1009
   s=MPI_WTIME()-wtime0
   d =int(s/3600.0/24.0)
   h =int((s-d*3600*24)/3600)
   m =int(s-d*3600*24-h*3600)/60
   call date_and_time(date=str_d1,time=str_t1)
   open(fid,file=trim(filenm),status="old",position="append")
     write(fid,*) '------------------------------------'
     write(fid,*) 'the program'
     write(fid,*) '  begins from ',str_d0,  &
                     ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
     write(fid,*) '  finish at ',str_d1,  &
                     ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
     write(fid, * )  'the total run time is'
     write(fid,'(3(i10,a10))') d,'day',h,'hour',m,'minute'
   close(fid)
end subroutine swmpi_time_end
function set_mpi_prefix(i,j,k) result(filenm)
    integer,intent(in) :: i,j,k
    character (len=12) :: filenm
    character (len=2) :: str1,str2,str3
    write(str1,"(i2.2)") i
    write(str2,"(i2.2)") j
    write(str3,"(i2.2)") k
    filenm  ='swmpi'//str1//str2//str3//'_'
end function set_mpi_prefix
function set_mpi_subfix(i,j,k) result(filenm)
    integer,intent(in) :: i,j,k
    character (len=9) :: filenm
    character (len=2) :: str1,str2,str3
    write(str1,"(i2.2)") i
    write(str2,"(i2.2)") j
    write(str3,"(i2.2)") k
    filenm  ='mpi'//str1//str2//str3
end function set_mpi_subfix
subroutine swmpi_change_fnm(i,j,k)
  integer,intent(in) :: i,j,k
  character (len=SEIS_STRLEN) :: str1,str2,str3
  write(str1,"(i2.2)") i
  write(str2,"(i2.2)") j
  write(str3,"(i2.2)") k
  fnm_swmpi ='swmpi'//trim(adjustl(str1))  &
                    //trim(adjustl(str2))  &
                    //trim(adjustl(str3))  &
                    //'_'
end subroutine swmpi_change_fnm

function swmpi_rename_fnm(pnm_input,fnm_in) result(fnm_out)
  character (len=*),intent(in),optional :: fnm_in,pnm_input
  character (len=SEIS_STRLEN) :: fnm_out
  if (present(pnm_input) .and. present(fnm_in)) then
  fnm_out=trim(pnm_input)//'/'//trim(fnm_swmpi)//trim(fnm_in)
  elseif (present(fnm_in)) then
  fnm_out=trim(fnm_swmpi)//trim(fnm_in)
  else
    print *, 'Error: no fnm_in passed into swmpi_rename_fnm'
    stop 1
  end if
  !fnm_out=trim(fnm_swmpi)//trim(fnm_in)
end function swmpi_rename_fnm

function swmpi_globi(i,n_i) result(gi)
integer,intent(in) :: i,n_i
integer :: gi
gi=(i-ni1+1)+n_i*ni+(ni1-1) ! term of (ni1-1) is for clarity
end function swmpi_globi
function swmpi_globj(j,n_j) result(gj)
integer,intent(in) :: j,n_j
integer :: gj
gj=(j-nj1+1)+n_j*nj+(nj1-1)
end function swmpi_globj
function swmpi_globk(k,n_k) result(gk)
integer,intent(in) :: k,n_k
integer :: gk
gk=(k-nk1+1)+n_k*nk+(nk1-1)
end function swmpi_globk

function swmpi_locli(gi,n_i) result(i)
integer,intent(in) :: gi,n_i
integer :: i
!i=mod(gi-(ni1-nx1),ni)+(ni1-nx1)
i=gi-n_i*ni
end function swmpi_locli
function swmpi_loclj(gj,n_j) result(j)
integer,intent(in) :: gj,n_j
integer :: j
!j=mod(gj-(nj1-ny1),nj)+(nj1-ny1)
j=gj-n_j*nj
end function swmpi_loclj
function swmpi_loclk(gk,n_k) result(k)
integer,intent(in) :: gk,n_k
integer :: k
!k=mod(gk-(nk1-nz1),nk)+(nk1-nz1)
k=gk-n_k*nk
end function swmpi_loclk

function point_in_thisnode(i,j,k) result(iflag)
integer,intent(in) :: i,j,k
logical iflag
integer glob_i1,glob_i2,glob_j1,glob_j2,glob_k1,glob_k2
glob_i1=swmpi_globi(ni1,thisid(1)); glob_i2=swmpi_globi(ni2,thisid(1))
glob_j1=swmpi_globj(nj1,thisid(2)); glob_j2=swmpi_globj(nj2,thisid(2))
glob_k1=swmpi_globk(nk1,thisid(3)); glob_k2=swmpi_globk(nk2,thisid(3))
iflag=.false.
if (      i>=glob_i1 .and. i<=glob_i2  &
    .and. j>=glob_j1 .and. j<=glob_j2  &
    .and. k>=glob_k1 .and. k<=glob_k2 ) then
   iflag=.true.
end if
end function point_in_thisnode

subroutine swmpi_except(msg)
    character (len=*),intent(in) :: msg
    integer :: ierr
    print *, trim(msg)
    call MPI_ABORT(SWMPI_COMM,1,ierr)
end subroutine swmpi_except

end module mpi_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
