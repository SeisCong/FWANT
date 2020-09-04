program solver_damp_smooth

use LSQR_mod
implicit none

!_________________________________________________________________________
!solver_WS.f:  
!
!  WS is for an ad hoc weighting<weig_sta> imposed on STATION CORRECTION
!
!   output files:
!    "try.sum"-- summary of fitting results including
!     VARIANCE REDUCTION, MISFIT, MODEL VARIANCE and MODEL NORM
!    "try.xyz"-- the solved solution of m which is organized in arry
!     (1,mvar)
!
!   input files and parameters:
!    fnm_list-- kernel file list and measurement,weight,and station indx
!
!    the parameter <mode (imultis)> controls several inversion options that
!    mode=1--> Simple damped LS inversion of Gm=d according to the original
!              pixel parameterization
!        =2--> MultiScaled inversion; 3D linear interpolatory wavelet
!              check solver_LSQR_pixel.f for 3D Haar or cubic splines
!              in subroutine <face_lift> (bior1.1 or bior3.1)
!              check solver_LSQR_node.f for linear interpolatory wavelet
!        =3--> Convulutional quelling that behaves similar as MINIMUN CURVATURE
!              but is more flexible by convoling with  GAUSSIAN function.
!              The parameter <sigma0> controls the one sigma # of grids of the
!              imposed a priori smoothness preference, bigger <sigma0> enforces
!              smoother solution!
!
!    Be sure to compare with <solver_SVD.f> for other insights
!
!   08/13 Y. Shen: change mdata and station location files
!                   change the defination of variance (fit in the program).
!   08/12/2008 W. Zhang: port to fortran90 format and only keep simple mode

character (len=*),parameter ::  &
  hline="----------------------------------------------------"
character (len=300) :: rcdstr

character (len=132) :: fnm_G,fnm_list,fnm_smooth
character (len=10) :: mode0(4),mode

real(DP) :: ker_thres,data_thres
real(DP) :: lamb,eta,weig_sta,weig_evt,weig_loc,tmp
real(SP) :: rhs,weig_dat
real(SP) :: rhs_min,rhs_max,rhs_min_weig,rhs_max_weig
real(SP),dimension(3) :: ker0

integer :: imultis
integer :: num_cmp,nsta,nevt,nloc
integer :: nrow_smoth,nbd_smoth,nel_smoth
integer ::                                     &
  num_blk,num_dat,num_mval,num_xval,num_siz1d, &
  num_sta,num_evt,num_loc,num_nel,num_row,num_bd
integer,allocatable :: idsta(:),idevt(:),idloc(:)
real(DP),allocatable :: wgloc(:)

integer :: i,j,k,m,n,ncoef,nit,nael
integer :: listid,gid,ierr

integer,dimension(:),allocatable :: indx0,indx
real(SP),dimension(:),allocatable :: coef_sp
real(DP),dimension(:),allocatable :: coef_dp
real(DP),dimension(:),allocatable :: mval,errx

mode0=(/'Simple    ','Multiscale','Quelling  ','HybridVQLM'/)

! ------------------------ input parameters ------------------------
write(*,*) hline
write(*,*)' mode:'
!write(*,*)'     0 x= Resynthesizing to different level using'
!write(*,*)'          previous results!'
write(*,*)'     1 = Simple damping;'
!write(*,*)'     2 x= MultiScaled transformation:'
!write(*,*)'     3 x= Convolution quelling.'
!write(*,*)'     4 x= vertically quelling and laterally multiscale.'
read(*,*) imultis
mode=mode0(imultis)

if (imultis/=1) call error_except("only Simple damping valid")
 
print *,'  file list of G matrix, data and weighting='
read(*,'(a)') fnm_list
write(*,*) " ",trim(fnm_list)

print *,'  file name of smoothing operator='
read(*,'(a)') fnm_smooth
write(*,*) " ",trim(fnm_smooth)

write(6,*)'  # of inverted structure components'
read(5,*) num_cmp
write(*,*) num_cmp

write(*,*)'  # of station terms per measurement and weighting'
read(*,*) nsta,weig_sta
write(*,*) nsta,weig_sta

write(*,*)'  # of event terms per measurement and weighting'
read(*,*) nevt,weig_evt
write(*,*) nevt,weig_evt

write(*,*)'  # of location terms per measurement and weighting'
read(*,*) nloc,weig_loc
write(*,*) nloc,weig_loc

write(*,*)' dampping and smoothing'
read(*,*) lamb,eta
write(*,*) lamb,eta

write(6,*)'  kernel and measurement threshold'
read(5,*) ker_thres, data_thres
write(*,*) " ",ker_thres, data_thres

!--------------------------- initial ---------------------------
call init_smoothing(fnm_smooth,num_cmp,nrow_smoth,nbd_smoth,nel_smoth)

call init_data_damp(fnm_list,num_cmp,nsta,nevt,nloc,    &
     num_dat,num_sta,num_evt,num_loc,num_mval,num_xval, &
     num_siz1d,num_row,num_bd,num_nel)

num_row=num_row+nrow_smoth
num_nel=num_nel+nel_smoth
num_bd=max(num_bd,nbd_smoth)
num_siz1d=max(num_siz1d,num_row)

write(*,*) '  # of cmp, sta, evt, loc per data'
write(*,*) num_cmp,nsta,nevt,nloc
write(*,*) '  # of dat, mval, xval, size'
write(*,*) num_dat,num_mval,num_xval,num_siz1d
write(*,*) '  # of row, band, nel'
write(*,*) num_row,num_bd,num_nel

call alloc_var(nsta,nevt,nloc,num_xval,num_bd)
call lsqr_init(num_row,num_xval,num_siz1d,num_nel)

listid=10001; gid=10002
rhs_min=1.e7; rhs_max=-1.e7
rhs_min_weig=1.e7; rhs_max_weig=-1.e7
!----------------------- read in G matrix ------------------------
open(listid,file=trim(fnm_list),status='old',iostat=ierr)
if (ierr>0) call error_except("fnm_list open err:"//trim(fnm_list))

list_loop: do
  read(listid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  rcdstr=adjustl(rcdstr); if (rcdstr(1:1)=='#') cycle

  ! kernel
  nael=0
  do k=1,num_cmp
     n=index(rcdstr(1:len_trim(rcdstr))," ")
     fnm_G=rcdstr(1:n)
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
  read(rcdstr,*) rhs, weig_dat
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  ! check measurement
  if (abs(rhs)>data_thres) then
     write(*,*) "omit large rhs=",rhs
     cycle
     !call error_except("measurement is too large")
  end if

  do k=1,nsta
      read(rcdstr,*) idsta(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  end do

  do k=1,nevt
      read(rcdstr,*) idevt(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  end do

  do k=1,nloc
      read(rcdstr,*) idloc(k), wgloc(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  end do

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
  
  ! station term
  do k=1,nsta
     m=m+1
     indx(m)=num_mval+abs(idsta(k))
     coef_dp(m)=sign(weig_sta,real(idsta(k),DP))
  end do
  
  ! event term
  do k=1,nevt
     m=m+1
     indx(m)=num_mval+num_sta+abs(idevt(k))
     coef_dp(m)=sign(weig_evt,real(idevt(k),DP))
  end do

  ! location term
  do k=1,nloc
     m=m+1
     indx(m)=num_mval+num_sta+num_evt+idloc(k)
     coef_dp(m)=wgloc(k)*weig_loc
  end do

  rhs_min=min(rhs_min,rhs)
  rhs_max=max(rhs_max,rhs)

  ! covariance matrix
  do i=1,m
     coef_dp(i)=coef_dp(i)*weig_dat
  end do
  rhs=rhs*weig_dat

  rhs_min_weig=min(rhs_min_weig,rhs)
  rhs_max_weig=max(rhs_max_weig,rhs)

  call lldrow(coef_dp,indx,m,real(rhs,DP))
end do list_loop

write(*,*) "    rhs_min,rhs_max"
write(*,*) rhs_min,rhs_max
write(*,*) "    rhs_min_weig,rhs_max_weig"
write(*,*) rhs_min_weig,rhs_max_weig
close(listid)

call lsqr_check("  ---- data, station and event ----")

!--------------------------- Dampping ------------------------------
do i=1,num_xval
!do i=1,num_mval
   !write(*,*) i
   coef_dp(1)=lamb
   indx(1)=i
   call lldrow(coef_dp,indx,1,0.0_DP)
end do

call lsqr_check("  ---- data, station and event + dampping ----")

!--------------------------- smoothing ------------------------------
open(gid,file=trim(fnm_smooth),status='old',iostat=ierr)
if (ierr>0) call error_except("open smoothing file err:"//trim(fnm_smooth))
read(gid,*) num_blk,nrow_smoth,nbd_smoth,nel_smoth
do n=1,nrow_smoth
   read(gid,*) ncoef,(indx0(i),i=1,ncoef),(coef_sp(i),i=1,ncoef)
   do k=1,num_cmp
      do i=1,ncoef
         indx(i)=indx0(i)+(k-1)*num_blk
         coef_dp(i)=coef_sp(i)*eta
      end do
      call lldrow(coef_dp,indx,ncoef,0.0_DP)
   end do
end do
close(gid)

call lsqr_check("  ---- data, station and event + dampping + smoothing ----")

!--------------------------- solving -------------------------
write(*,*) "    begin solving"
call slsqr(nit)
call lgtsol(mval,errx,num_xval)
write(*,*) "    finish solving"
!do i=num_mval+1,num_mval+num_sta
!   mval(i)=mval(i)*weig_sta
!   errx(i)=errx(i)*weig_sta
!end do
!do i=num_mval+num_sta+1,num_mval+num_sta+num_evt
!   mval(i)=mval(i)*weig_evt
!   errx(i)=errx(i)*weig_evt
!end do

!--------------------------- output ----------------------------
call export_result

call result_evaluate(fnm_list,fnm_smooth)

call lsqr_destroy

write(*,*) hline

!----------------------------------------------------------
contains
!----------------------------------------------------------

subroutine alloc_var(nst,nev,nlo,nx,nbd)
  integer :: nst,nev,nlo,nx,nbd
  if (nst>0) then
     allocate(idsta(nst)); idsta=0.0
  end if
  if (nev>0) then
     allocate(idevt(nev)); idevt=0.0
  end if
  if (nlo>0) then
     allocate(idloc(nlo)); idloc=0.0
     allocate(wgloc(nlo)); wgloc=0.0_DP
  end if

  allocate(indx0(nbd)); indx0=0
  allocate(indx(nbd)); indx=0

  allocate(coef_sp(nbd));coef_sp=0.0
  allocate(coef_dp(nbd));coef_dp=0.0_DP
  allocate(mval(nx));mval=0.0_DP
  allocate(errx(nx));errx=0.0_DP
end subroutine alloc_var

subroutine init_smoothing(filenm,ncmp,max_row,max_bd,max_nel)
character (len=*),intent(in) :: filenm
integer,intent(in) :: ncmp
integer,intent(out) :: max_row,max_bd,max_nel
integer :: fid,nblk,ierr

fid=5001;
max_row=0; max_bd=0; max_nel=0
open(fid,file=trim(filenm),status='old',iostat=ierr)
if (ierr>0) call error_except("init_smoothing open err:"//trim(filenm))
read(fid,*) nblk,max_row,max_bd,max_nel
close(fid)
max_bd =max_bd *ncmp
max_row=max_row*ncmp
max_nel=max_nel*ncmp
end subroutine init_smoothing

subroutine init_data_damp(filenm,ncmp,nst,nev,nlo,  &
    max_dat,max_st,max_ev,max_lo,max_mval,max_xval, &
    max_siz,max_row,max_bd,max_nel)
character (len=*),intent(in) :: filenm
integer,intent(in) :: ncmp,nst,nev,nlo
integer,intent(out) ::                     &
  max_dat,max_st,max_ev,max_lo,max_mval,max_xval, &
  max_siz,max_row,max_bd,max_nel
character (len=132) :: fnm_k
integer :: idst(nst),idev(nev),idlo(nlo)
real    :: weig(nlo)
integer :: fid,gid,k,n,nblk,nael,ierr
real(SP) :: ker0

fid=5001; gid=5002
max_st=0
max_ev=0
max_lo=0
max_dat=0
max_mval=0 ! model parameterization
max_xval=0 ! m plus station and event terms
max_siz=0  ! max size of vector
max_row=0
max_bd=0
max_nel=0

! number related with kernels
open(fid,file=trim(filenm),status='old',iostat=ierr)
if (ierr>0) call error_except("init_G_scale open err:"//trim(filenm))
do
  read(fid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  rcdstr=adjustl(rcdstr); if (rcdstr(1:1)=='#') cycle

  max_dat=max_dat+1
  nael=0
  do k=1,ncmp
     n=index(rcdstr(1:len_trim(rcdstr))," ")
     fnm_k=rcdstr(1:n)
     !write(*,*) trim(fnm_k)
     rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
     open(gid,file=trim(fnm_k),status='old',iostat=ierr,form='unformatted')
     if (ierr>0) call error_except("init_G_scale kfile err:"//trim(filenm))
     read(gid) nblk, ker0
     read(gid) ncoef
     close(gid)
     nael=nael+ncoef
  end do

  read(rcdstr,*) rhs, weig_dat
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

  do k=1,nst
      read(rcdstr,*) idst(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      max_st=max(max_st,abs(idst(k)))
  end do

  do k=1,nev
      read(rcdstr,*) idev(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      max_ev=max(max_ev,abs(idev(k)))
  end do

  do k=1,nlo
      read(rcdstr,*) idlo(k), weig(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      max_lo=max(max_lo,abs(idlo(k)))
  end do

  max_bd=max(max_bd,nael)
  if(max_nel+nael<0) then
    print *,'serious error: "max_nel" overflow ',nael,max_nel 
    stop
  endif
  max_nel=max_nel+nael
end do
close(fid)

max_mval=nblk*ncmp
max_xval=max_mval+max_st+max_ev+max_lo
max_nel=max_nel+(nst+nev+nlo)*max_dat +max_xval
max_row=max_dat+max_xval
max_bd=max_bd+nst+nev+nlo
max_siz=max(max_xval,max_row)
end subroutine init_data_damp

!----------------------------------------------------------
subroutine result_evaluate(fnm_list,fnm_smo)
character (len=*),intent(in) :: fnm_list,fnm_smo
real(DP) :: fit,fit1,dvar,fit_weig,dvar_weig,fit_kai2
real(DP) :: mean,vmin,vmax,cov,norm,smocov
real(DP) :: vmin_sta,vmax_sta,vmin_evt,vmax_evt,vmin_loc,vmax_loc
real(DP) :: tmp,tmp1,tmp_weig
real(SP) :: rhs,rhs_weig
integer :: i,j,k,n,m,m1,nblk,nael,ierr,fid,gid,ndat
real(SP) :: ker0

!------------- data variance reduction ------------------
fit=0.0_DP; fit1=0.0_DP; dvar=0.0_DP
fit_weig=0.0_DP; dvar_weig=0.0_DP; fit_kai2=0.0_DP
ndat=0
vmin_loc=1.d7; vmax_loc=-1.d7

fid=5001; gid=5002
open(fid,file=trim(fnm_list),status='old',iostat=ierr)
open(76,file='try.prd')
do
  read(fid,'(a300)',iostat=ierr) rcdstr
  if (ierr<0) exit
  rcdstr=adjustl(rcdstr); if (rcdstr(1:1)=='#') cycle

  ndat=ndat+1

  ! kernel
  nael=0
  do k=1,num_cmp
     n=index(trim(rcdstr(1:len_trim(rcdstr)))," ")
     fnm_G=rcdstr(1:n)
     rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)

     open(gid,file=trim(fnm_G),status='old',iostat=ierr,form='unformatted')
     if (ierr>0) call error_except("open kernel err:"//trim(fnm_G))
     read(gid) nblk,ker0
     read(gid) ncoef
     read(gid) (indx0(nael+i),coef_sp(nael+i),i=1,ncoef)
     close(gid)

     do i=1,ncoef
        indx0(nael+i)=indx0(nael+i)+nblk*(k-1)
     end do
     nael=nael+ncoef
  end do

  read(rcdstr,*) rhs, weig_dat
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  ! check measurement
  if (abs(rhs)>data_thres) then
     !write(*,*) "omit large rhs=",rhs
     cycle
  end if

  do k=1,nsta
      read(rcdstr,*) idsta(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  end do

  do k=1,nevt
      read(rcdstr,*) idevt(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  end do

  do k=1,nloc
      read(rcdstr,*) idloc(k), wgloc(k)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
      n=index(rcdstr(1:len_trim(rcdstr))," ");rcdstr(1:)=rcdstr(n+1:); rcdstr=adjustl(rcdstr)
  end do

  rhs_weig=rhs*weig_dat

  ! data variance
  dvar=dvar+dble(rhs)*dble(rhs)
  dvar_weig=dvar_weig+dble(rhs_weig)*dble(rhs_weig)

  ! synthetic data
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
  m1=m
  ! station term
  do k=1,nsta
     m=m+1
     j=num_mval+abs(idsta(k))
     indx(m)=j
     coef_dp(m)=sign(weig_sta,real(idsta(k),DP))
  end do
  ! event term
  do k=1,nevt
     m=m+1
     j=num_mval+num_sta+abs(idevt(k))
     indx(m)=j
     coef_dp(m)=sign(weig_evt,real(idevt(k),DP))
  end do
  ! location term
  do k=1,nloc
     m=m+1
     indx(m)=num_mval+num_sta+num_evt+idloc(k)
     coef_dp(m)=wgloc(k)*weig_loc
  end do
  do k=1,nloc
     j=num_mval+num_sta+num_evt+idloc(k)
     vmin_loc=dmin1(mval(j)*wgloc(k)*weig_loc,vmin_loc)
     vmax_loc=dmax1(mval(j)*wgloc(k)*weig_loc,vmax_loc)
  end do

  tmp=0.0_DP
  do i=1,m
     j=indx(i)
     tmp=tmp+coef_dp(i)*mval(j)
  end do
  tmp1=0.0_DP
  do i=1,m1
     j=indx(i)
     tmp1=tmp1+coef_dp(i)*mval(j)
  end do
  tmp_weig=tmp*weig_dat

  fit=fit+(tmp-dble(rhs))*(tmp-dble(rhs))
  fit1=fit1+(tmp1-dble(rhs))*(tmp1-dble(rhs))
  fit_weig=fit_weig+(tmp_weig-dble(rhs_weig))*(tmp_weig-dble(rhs_weig))
  write(76,'(e12.5,3(2x,e12.5))') tmp1,dble(rhs),tmp_weig, dble(rhs_weig)
end do
close(fid)
close(76)

! kai^2 in Montelli, et al. GJI(2004),158:637
fit_kai2 = fit_weig/ndat
!  fit=(1.-dsqrt(fit/dvar))*100.
!  Y. Shen, 8/13/2004
!  I define variance as deltT*deltT, thus variance reduction is
!  (1 - {sum[deltT(final)]**2/sum[deltT(original)]**2})*100.  
fit = (1. - (fit/dvar))*100.0
fit1= (1. - (fit1/dvar))*100.0
fit_weig = (1. - (fit_weig/dvar_weig))*100.0

!------------------------- smooth --------------------------
smocov=0.0_DP
open(fid,file=trim(fnm_smo),status='old',iostat=ierr)
read(fid,*) num_blk,nrow_smoth,nbd_smoth,nel_smoth
do n=1,nrow_smoth
   read(fid,*) ncoef,(indx0(i),i=1,ncoef),(coef_sp(i),i=1,ncoef)
   do k=1,num_cmp
      tmp=0.0_DP
      do i=1,ncoef
         j=indx0(i)+(k-1)*num_blk
         tmp=tmp+coef_sp(i)*mval(j)
      end do
      smocov=smocov+tmp*tmp
   end do
end do
smocov = dsqrt(smocov/dble(num_mval))
close(fid)

!------------------------- statistics --------------------------
!  Y.Shen, 08/16/04
!  much of the model space is unsampled or poorly samples.  It is
!  more meaningful to evaluate model norm only at the grids that are
!  sampled as measured by Gt*G (see Hung, Shen, Chiao, 2004).
!
!  This should be done separately (BUILDG/gtgdiag.f) and the result to be read in
!       open(81,file='/home/yang/FFT/africa/BUILDG/diag.out')
!       do k=1,nvar
!       read(81,*) ndum,gtgdiag(k)
!       end do
!       close(81)
!  calculate model norm etc.

mean=0.d0; cov=0.d0; norm=0.d0
vmin=1.d7; vmax=-1.d7;
vmin_sta=1.d7; vmax_sta=-1.d7
vmin_evt=1.d7; vmax_evt=-1.d7

do j=1,num_mval
   cov=cov+errx(j)*errx(j)
   norm=norm+mval(j)*mval(j)
   mean=mean+mval(j)
   vmin=dmin1(mval(j),vmin)
   vmax=dmax1(mval(j),vmax)
end do
do j=num_mval+1,num_mval+num_sta
   vmin_sta=dmin1(mval(j)*weig_sta,vmin_sta)
   vmax_sta=dmax1(mval(j)*weig_sta,vmax_sta)
end do
do j=num_mval+num_sta+1,num_mval+num_sta+num_evt
   vmin_evt=dmin1(mval(j)*weig_evt,vmin_evt)
   vmax_evt=dmax1(mval(j)*weig_evt,vmax_evt)
end do
!norm=dsqrt(norm)/dble(nvar)
!cov=dsqrt(cov)/dble(nvar)
!mean=mean/dble(nvar)
norm = dsqrt(norm/dble(num_mval))
cov  = dsqrt(cov/dble(num_mval))
mean =mean/dble(num_mval)
!write(*,*) 'number of points with gtgdiag > ',gtgthres, nk

!------------------------- output --------------------------
! summary of fitting results
open(77,file='try.sum')
write(77,*) "mode = ", mode
!write(77,*)' Remember that model variance can only be compared'
!write(77,*)'  under the same mode !!'
write(77,*) "damping = ", lamb
write(77,*) "smoothing = ", eta
write(77,*) "percent of data_var reduction = ", fit
write(77,*) "percent of data_var reduction of structure = ", fit1
write(77,*) "percent of data_var reduction w weight = ", fit_weig
write(77,*) "kai^2 = ", fit_kai2*ndat
write(77,*) "number of data = ", ndat
write(77,*) "kai^2/N = ", fit_kai2
write(77,*) "error_variance (model_var?) = ", cov
write(77,*) "smooth_norm = ", smocov
write(77,*) "model_norm (model_var?) = ", norm
write(77,*) "model_mean = ", mean
write(77,*) "model_min = ", vmin
write(77,*) "model_max = ", vmax
write(77,*) "station_term_min = ", vmin_sta
write(77,*) "station_term_max = ", vmax_sta
write(77,*) "event_term_min = ", vmin_evt
write(77,*) "event_term_max = ", vmax_evt
write(77,*) "location_term_min = ", vmin_loc
write(77,*) "location_term_max = ", vmax_loc
close(77)

write(*,*) "damping, smoothing, data_var_reduc, data_var_reduc_of_structure, model_norm"
write(*,*) lamb, eta, fit, fit1, norm
write(*,*) ""
end subroutine result_evaluate

!----------------------------------------------------------
subroutine export_result
integer :: i,j

! output result
open(78,file='try.xyz')
  write(78,'(1e12.3)')(mval(i),i=1,num_mval)
close(78)

open(78,file='try.err')
  write(78,'(1e12.3)')(errx(i),i=1,num_mval)
close(78)

! station term
open(79,file='try.sta')
do i=1,num_sta
   j=num_mval+i
   write(79,'(i5,f10.5,f10.5)') i,mval(j),mval(j)*weig_sta
end do
close(79)

! event term
open(80,file='try.evt')
do i=1,num_evt
   j=num_mval+num_sta+i
   write(80,'(i5,f10.5,f10.5)') i,mval(j),mval(j)*weig_evt
end do
close(80)

! location term
open(81,file='try.loc')
do i=1,num_loc
   j=num_mval+num_sta+num_evt+i
   write(81,'(i5,f10.5)') i,mval(j)
end do
close(81)

end subroutine export_result

end program solver_damp_smooth

