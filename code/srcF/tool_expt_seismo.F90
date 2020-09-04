program tool_expt_seismo

! This program retrieve seismogram data from output node nc file
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 23:30:18 -0500 (Wed, 14 Jan 2009) $
! $Revision: 53 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

use constants_mod
use string_mod, only : string_conf
use math_mod
use para_mod
use mpi_mod
use nfseis_mod
use grid_mod
use io_mod

implicit none

integer itype,id,indx,nt2
integer i,j,k,n1,n2,n3,n_i,n_j,n_k,npt

character (len=SEIS_STRLEN) :: fnmstr,pnm_expt,fnm_expt

real,dimension(:),allocatable :: &
     Vx,Vy,Vz,T

!---------------------------------------------

pnm_expt='./export'

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)

write(*,*) "Input type: recv(=1), line(=2) or snap(=3)? and maximum time count"
read(*,*) itype, nt2

if (itype==1) then
   write(*,*) "Input point index"
   read(*,*) indx
   id=0
   fnm_expt= trim(pnm_expt)//"/"//'recv_'       &
           //"pt"//trim(io_out_pattern(indx,5)) &
           //".vel"
elseif (itype==2) then
   write(*,*) "Input line id and point index"
   read(*,*) id,indx
   fnm_expt= trim(pnm_expt)//"/"                 &
           //'line_'//trim(io_out_pattern(id,3)) &
           //"_pt"//trim(io_out_pattern(indx,5)) &
           //".vel"
elseif (itype==3) then
   write(*,*) "Input snap id and point index"
   read(*,*) id, i,j,k
   fnm_expt= trim(pnm_expt)//"/"                 &
           //'snap_'//trim(io_out_pattern(id,3)) &
           //"_n1"//trim(io_out_pattern(i,5)) &
           //"_n2"//trim(io_out_pattern(j,5)) &
           //"_n3"//trim(io_out_pattern(k,5)) &
           //".vel"
else
   write(*,*) "Wrong type", itype
   stop 1
end if

allocate( Vx(nt)); Vx=0.0
allocate( Vy(nt)); Vy=0.0
allocate( Vz(nt)); Vz=0.0
allocate( T(nt)); T=0.0

if (itype==1 .or. itype==2) then

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
if (seismo_on_this(pnm_grid,id,indx,n_i,n_j,n_k,npt)) then
   fnmstr=get_fnm_seismo(pnm_out,n_i,n_j,n_k)
   call retrieve_recvline(fnmstr,npt,'Vx',Vx,1,nt2,1,T)
   call retrieve_recvline(fnmstr,npt,'Vy',Vy,1,nt2,1)
   call retrieve_recvline(fnmstr,npt,'Vz',Vz,1,nt2,1)
   exit
endif
end do
end do
end do

else

   call locate_snap_indx(id,i,j,k,n_i,n_j,n_k,n1,n2,n3)
   fnmstr=trim(pnm_out)                         &
       //'/'//trim(io_enum('snap_',id))         &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_time(fnmstr,T,nt2)
   call retrieve_snap_seis(fnmstr,n1,n2,n3,'Vx',Vx,nt2)
   call retrieve_snap_seis(fnmstr,n1,n2,n3,'Vy',Vy,nt2)
   call retrieve_snap_seis(fnmstr,n1,n2,n3,'Vz',Vz,nt2)

end if

call export_ascii(fnm_expt,T,Vx,Vy,Vz,nt2)
!call export_hfin_bin(fnm_expt,tdelta,t,Vx,Vy,Vz,m)

deallocate( Vx)
deallocate( Vy)
deallocate( Vz)
deallocate( T)

print *, 'export seismo finished'

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine export_hfin_bin( filenm,stept,t,Vx,Vy,Vz,nt)
character (len=*) :: filenm
real stept
real,dimension(:) :: t,Vx,Vy,Vz
integer nt,fid,i
fid=1001
open(fid,file=trim(filenm),status='unknown',form='unformatted')
	write (fid) nt
	write (fid) stept
	do i=1,nt
	write (fid) t(i),Vx(i),Vy(i),Vz(i)
	end do
close(fid)
end subroutine export_hfin_bin
subroutine export_ascii(filenm,t,Vx,Vy,Vz,nt)
character (len=*) :: filenm
real,dimension(:) :: t,Vx,Vy,Vz
integer nt,fid,i
fid=1001
open(fid,file=trim(filenm),status='unknown')
	do i=1,nt
	write (fid,*) t(i),Vx(i),Vy(i),Vz(i)
	end do
close(fid)
end subroutine export_ascii

subroutine locate_snap_indx(id,i,j,k,n_i,n_j,n_k,n1,n2,n3)
integer id,i,j,k,n_i,n_j,n_k,n1,n2,n3
integer gi,gj,gk

gi=(i-1)*snap_subt(1,id)+snap_subs(1,id)
gj=(j-1)*snap_subt(2,id)+snap_subs(2,id)
gk=(k-1)*snap_subt(3,id)+snap_subs(3,id)

n_i=(gi-ni1)/ni
n_j=(gj-nj1)/nj
n_k=(gk-nk1)/nk

call swmpi_change_fnm(n_i,n_j,n_k)
call swmpi_set_gindx(n_i,n_j,n_k)

n1=(gi-max(ngi1,snap_subs(1,id)))/snap_subt(1,id)+1
n2=(gj-max(ngj1,snap_subs(2,id)))/snap_subt(2,id)+1
n3=(gk-max(ngk1,snap_subs(3,id)))/snap_subt(3,id)+1
end subroutine locate_snap_indx

!subroutine retrieve_snap(id,i,j,k,nt2)
!logical iflag
!tdelta=snap_tinv(id)*stept
!n=0; m=0
!do
!   if (m>=nt2) exit
!   n=n+1
!   fnm_nc=get_fnm_stressnode_n(pnm_out,'snap_',id,n,n_i,n_j,n_k)
!   inquire(file=trim(fnm_nc),exist=iflag)
!   if (.not. iflag .and. n==1) then
!      print *, trim(fnm_nc), ' does not exist'
!      stop 1
!   end if
!   if (.not. iflag) exit
!
!   call nfseis_diminfo(fnm_nc,'time',mt)
!   if (mt+m>nt2) mt=nt2-m
!   call nfseis_varget(fnm_nc,'Vx',Vx(m+1:m+mt),  &
!        (/n1,n2,n3,1/),(/1,1,1,mt/),(/1,1,1,1/))
!   call nfseis_varget(fnm_nc,'Vy',Vy(m+1:m+mt),  &
!        (/n1,n2,n3,1/),(/1,1,1,mt/),(/1,1,1,1/))
!   call nfseis_varget(fnm_nc,'Vz',Vz(m+1:m+mt),  &
!        (/n1,n2,n3,1/),(/1,1,1,mt/),(/1,1,1,1/))
!   call nfseis_varget(fnm_nc,'time',t(m+1:m+mt), &
!        (/1/),(/mt/),(/1/))
!   m=m+mt
!end do
!end subroutine retrieve_snap

end program tool_expt_seismo

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
