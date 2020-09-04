program inv_ckb_syn

implicit none
integer,parameter :: DP = KIND(1.0D0)
integer,parameter :: SP = KIND(1.0)

character (len=*),parameter ::  &
  hline="----------------------------------------------------"
character (len=300) :: rcdstr,synstr
character (len=25) :: rhsstr

character (len=132) :: fnm_G,fnm_list,fnm_syn,fnm_ckb

real(DP) :: ker_thres,data_thres
real(DP) :: rhs,tmp
real(SP) :: rhs_min,rhs_max
real(SP),dimension(3) :: ker0

integer :: num_cmp, num_blk,num_mval

integer :: i,j,k,m,n,ncoef,nit,nael
integer :: listid,fid,gid,ierr

integer,dimension(:),allocatable :: indx0,indx
real(SP),dimension(:),allocatable :: coef_sp
real(DP),dimension(:),allocatable :: coef_dp
real(DP),dimension(:),allocatable :: mval

write(6,*)'  # of inverted structure components'
read(5,*) num_cmp
write(*,*) num_cmp

print *,'  file name of checkerboard='
read(*,'(a)') fnm_ckb
write(*,*) " ",trim(fnm_ckb)

print *,'  file list of Gd from record='
read(*,'(a)') fnm_list
write(*,*) " ",trim(fnm_list)

print *,'  file list of Gd predicted='
read(*,'(a)') fnm_syn
write(*,*) " ",trim(fnm_syn)

write(*,*)'  kernel and measurement threshold'
read(*,*) ker_thres, data_thres
write(*,*) " ",ker_thres, data_thres

!--------------------------- init ckb---------------------------
listid=10001; gid=10002; fid=1003

open(fid,file=trim(fnm_ckb),status='old',iostat=ierr)
if (ierr>0) call error_except("fnm_ckb open err:"//trim(fnm_ckb))
read(fid,*) num_mval

call alloc_var(num_mval)

do n=1,num_mval
   read(fid,*) mval(n)
end do
close(fid)

!----------------------- read in G matrix ------------------------
rhs_min=1.e7; rhs_max=-1.e7

open(listid,file=trim(fnm_list),status='old',iostat=ierr)
if (ierr>0) call error_except("fnm_list open err:"//trim(fnm_list))

open(fid,file=trim(fnm_syn),status='unknown',iostat=ierr)
if (ierr>0) call error_except("fnm_syn open err:"//trim(fnm_syn))

list_loop: do
  read(listid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  rcdstr=adjustl(rcdstr); if (rcdstr(1:1)=='#') cycle
  synstr=""

  ! kernel
  nael=0
  do k=1,num_cmp
     n=index(rcdstr(1:len_trim(rcdstr))," ")
     fnm_G=rcdstr(1:n)
     synstr=trim(synstr)//" "//rcdstr(1:n)
     rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

     open(gid,file=trim(fnm_G),status='old',iostat=ierr,form='unformatted')
     if (ierr>0) call error_except("open kernel err:"//trim(fnm_G))
     read(gid) num_blk,ker0(k)
     read(gid) ncoef
     read(gid) (indx0(nael+i),coef_sp(nael+i),i=1,ncoef)
     close(gid)

     do i=1,ncoef
        indx0(nael+i)=indx0(nael+i)+num_blk*(k-1)
     end do
     nael=nael+ncoef
  end do
  ! measure, weight and station id
  read(rcdstr,*) rhs
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  ! check measurement
  if (abs(rhs)>data_thres) then
     write(*,*) "omit large rhs=",rhs
     cycle
     !call error_except("measurement is too large")
  end if

  ! check kernel value
  m=0
  do i=1,nael
     j=indx0(i)
     tmp=real(coef_sp(i),DP)
     if (dabs(tmp) > ker_thres)then
        m=m+1
        indx(m)=j
        coef_dp(m)=tmp
     end if
  end do

  ! synthetic data
  tmp=0.0_DP
  do i=1,m
     j=indx(i)
     tmp=tmp+coef_dp(i)*mval(j)
  end do
  rhs=real(tmp,SP)
  write(rhsstr,*) rhs

  synstr=trim(synstr)//" "//trim(rhsstr)//" "//trim(rcdstr)
  write(fid,"(a)") trim(synstr)
  
  rhs_min=min(rhs_min,rhs)
  rhs_max=max(rhs_max,rhs)

end do list_loop

write(*,*) "    rhs_min,rhs_max"
write(*,*) rhs_min,rhs_max
close(listid)
close(fid)

write(*,*) hline

!----------------------------------------------------------
contains
!----------------------------------------------------------

subroutine alloc_var(nx)
  integer :: nx

  allocate(indx0(nx)); indx0=0
  allocate(indx(nx)); indx=0

  allocate(coef_sp(nx));coef_sp=0.0
  allocate(coef_dp(nx));coef_dp=0.0_DP
  allocate(mval(nx));mval=0.0_DP
end subroutine alloc_var

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
  print *, trim(msg)
  stop 1
end subroutine error_except

end program inv_ckb_syn

