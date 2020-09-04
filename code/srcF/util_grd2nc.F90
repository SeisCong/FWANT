program grd2nc

! This program converts the file in PLOT3D format to netcdf file
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 13:02:04 -0500 (Fri, 16 Jan 2009) $
! $Revision: 72 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#define USE_GETARG

use netcdf
implicit none
integer,parameter :: SEIS_STRLEN=132
real,parameter :: SEIS_ZERO=1e-15
character (len=SEIS_STRLEN) :: fnm_plot3d,fnm_nc,fnm_grd
integer n,imax,jmax,kmax,i0,j0,k0
#ifdef USE_GETARG
integer numarg,iargc
#endif

type GRDSTCT
     integer ni,nj,nk
     integer i0,j0,k0
     real*8,dimension(:,:,:),pointer :: x,y,z
end type GRDSTCT
type (GRDSTCT),dimension(:),allocatable :: grd_list

i0=1;j0=1;k0=1;

#ifdef USE_GETARG
   numarg = iargc( )
   if (numarg<1) then
      print *, "useage: ./util_grd2nc fnm_grd i0 j0 k0"
      stop 1
   end if
   call getarg(1,fnm_plot3d)
   if (numarg>=4) then
      call getarg(2,fnm_nc); read(i0,*) fnm_nc
      call getarg(3,fnm_nc); read(j0,*) fnm_nc
      call getarg(4,fnm_nc); read(k0,*) fnm_nc
   end if
#else
   print *, 'please input the PLOT3D file name with/without subfix'
   read *, fnm_plot3d
#endif

n=len_trim(fnm_plot3d)

if (n>4 .and. fnm_plot3d(n-2:n)=='grd') then
   fnm_nc  =fnm_plot3d(1:n-4)//'.nc'
   fnm_grd =fnm_plot3d(1:n-4)//'.grd'
else
   fnm_nc  =trim(fnm_plot3d)//'.nc'
   fnm_grd =trim(fnm_plot3d)//'.grd'
end if

call read_grd2nc(trim(fnm_grd),imax,jmax,kmax)

call export_grd2nc(trim(fnm_nc),imax,jmax,kmax,i0,j0,k0)

call destroy_grdlist

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine read_grd2nc(fnm_grd,imax,jmax,kmax)
  character (len=*),intent(in) :: fnm_grd
  integer imax,jmax,kmax
  integer imin,jmin,kmin

  integer i,j,k,m,n,nmax,ni,nj,nk,mi,mj,mk,fid
  real*8 x0,y0,z0

  fid=1003
  open(fid,file=trim(fnm_grd),status='old')
  ! read subdomain counts
  read(fid,*) nmax
  allocate(grd_list(nmax))

  ! read subdomain dimensions
  do n=1,nmax
     read(fid,*)  ni,nj,nk
     grd_list(n)%ni=ni; grd_list(n)%nj=nj; grd_list(n)%nk=nk
     grd_list(n)%i0=1; grd_list(n)%j0=1; grd_list(n)%k0=1
     allocate(grd_list(n)%x(ni,nj,nk))
     allocate(grd_list(n)%y(ni,nj,nk))
     allocate(grd_list(n)%z(ni,nj,nk))
  end do

  ! read coord data
  do n=1,nmax
     ni=grd_list(n)%ni; nj=grd_list(n)%nj; nk=grd_list(n)%nk
     read(fid,*) ( ((grd_list(n)%x(i,j,k), i=1,ni ), j=1,nj),k=1,nk ),  &
                 ( ((grd_list(n)%y(i,j,k), i=1,ni ), j=1,nj),k=1,nk ),  &
                 ( ((grd_list(n)%z(i,j,k), i=1,ni ), j=1,nj),k=1,nk )
  end do
  close(fid)

  ! calculate indx
  do n=1,nmax
     ni=grd_list(n)%ni; nj=grd_list(n)%nj; nk=grd_list(n)%nk
  do m=n+1,nmax
     mi=grd_list(m)%ni; mj=grd_list(m)%nj; mk=grd_list(m)%nk

     !x1y1z1 -> x2 or y2 or z2
     x0=grd_list(m)%x(1,1,1); y0=grd_list(m)%y(1,1,1); z0=grd_list(m)%z(1,1,1)
       !x2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(ni,:,:),grd_list(n)%y(ni,:,:),grd_list(n)%z(ni,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+ni-1
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !y2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,nj,:),grd_list(n)%y(:,nj,:),grd_list(n)%z(:,nj,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0+nj-1
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !z2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,nk),grd_list(n)%y(:,:,nk),grd_list(n)%z(:,:,nk), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0+nk-1
     end if

     !x2y1z1 -> x1 or y2 or z2
     x0=grd_list(m)%x(mi,1,1); y0=grd_list(m)%y(mi,1,1); z0=grd_list(m)%z(mi,1,1)
       !x1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(1,:,:),grd_list(n)%y(1,:,:),grd_list(n)%z(1,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0-mi+1
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !y2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,nj,:),grd_list(n)%y(:,nj,:),grd_list(n)%z(:,nj,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0+nj-1
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !z2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,nk),grd_list(n)%y(:,:,nk),grd_list(n)%z(:,:,nk), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0+nk-1
     end if

     !x1y2z1 -> x2 or y1 or z2
     x0=grd_list(m)%x(1,mj,1); y0=grd_list(m)%y(1,mj,1); z0=grd_list(m)%z(1,mj,1)
       !x2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(ni,:,:),grd_list(n)%y(ni,:,:),grd_list(n)%z(ni,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+ni-1
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !y1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,1,:),grd_list(n)%y(:,1,:),grd_list(n)%z(:,1,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0-mj+1
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !z2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,nk),grd_list(n)%y(:,:,nk),grd_list(n)%z(:,:,nk), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0+nk-1
     end if

     !x2y2z1 -> x1 or y1 or z2
     x0=grd_list(m)%x(mi,mj,1); y0=grd_list(m)%y(mi,mj,1); z0=grd_list(m)%z(mi,mj,1)
       !x1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(1,:,:),grd_list(n)%y(1,:,:),grd_list(n)%z(1,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0-mi+1
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !y1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,1,:),grd_list(n)%y(:,1,:),grd_list(n)%z(:,1,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0-mj+1
         grd_list(m)%k0=grd_list(n)%k0+k-1
     end if
       !z2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,nk),grd_list(n)%y(:,:,nk),grd_list(n)%z(:,:,nk), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0+nk-1
     end if

     ! z2 other

     !x1y1z2 -> x2 or y2 or z1
     x0=grd_list(m)%x(1,1,mk); y0=grd_list(m)%y(1,1,mk); z0=grd_list(m)%z(1,1,mk)
       !x2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(ni,:,:),grd_list(n)%y(ni,:,:),grd_list(n)%z(ni,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+ni-1
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !y2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,nj,:),grd_list(n)%y(:,nj,:),grd_list(n)%z(:,nj,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0+nj-1
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !z1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,1),grd_list(n)%y(:,:,1),grd_list(n)%z(:,:,1), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0-mk+1
     end if

     !x2y1z2 -> x1 or y2 or z1
     x0=grd_list(m)%x(mi,1,mk); y0=grd_list(m)%y(mi,1,mk); z0=grd_list(m)%z(mi,1,mk)
       !x1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(1,:,:),grd_list(n)%y(1,:,:),grd_list(n)%z(1,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0-mi+1
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !y2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,nj,:),grd_list(n)%y(:,nj,:),grd_list(n)%z(:,nj,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0+nj-1
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !z1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,1),grd_list(n)%y(:,:,1),grd_list(n)%z(:,:,1), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0+j-1
         grd_list(m)%k0=grd_list(n)%k0-mk+1
     end if

     !x1y2z2 -> x2 or y1 or z1
     x0=grd_list(m)%x(1,mj,mk); y0=grd_list(m)%y(1,mj,mk); z0=grd_list(m)%z(1,mj,mk)
       !x2
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(ni,:,:),grd_list(n)%y(ni,:,:),grd_list(n)%z(ni,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+ni-1
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !y1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,1,:),grd_list(n)%y(:,1,:),grd_list(n)%z(:,1,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0-mj+1
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !z1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,1),grd_list(n)%y(:,:,1),grd_list(n)%z(:,:,1), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-1
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0-mk+1
     end if

     !x2y2z2 -> x1 or y1 or z1
     x0=grd_list(m)%x(mi,mj,mk); y0=grd_list(m)%y(mi,mj,mk); z0=grd_list(m)%z(mi,mj,mk)
       !x1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(1,:,:),grd_list(n)%y(1,:,:),grd_list(n)%z(1,:,:), &
          j,k) ) then
         grd_list(m)%i0=grd_list(n)%i0-mi+1
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !y1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,1,:),grd_list(n)%y(:,1,:),grd_list(n)%z(:,1,:), &
          i,k) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0-mj+1
         grd_list(m)%k0=grd_list(n)%k0+k-mk
     end if
       !z1
     if ( is_in_plane(x0,y0,z0, &
          grd_list(n)%x(:,:,1),grd_list(n)%y(:,:,1),grd_list(n)%z(:,:,1), &
          i,j) ) then
         grd_list(m)%i0=grd_list(n)%i0+i-mi
         grd_list(m)%j0=grd_list(n)%j0+j-mj
         grd_list(m)%k0=grd_list(n)%k0-mk+1
     end if

  end do !m
  end do !n

  ! find i[jk]max[min]
  imin=1;jmin=1;kmin=1
  imax=0;jmax=0;kmax=0
  do n=1,nmax
     if (imin>grd_list(n)%i0) imin=grd_list(n)%i0
     if (jmin>grd_list(n)%j0) jmin=grd_list(n)%j0
     if (kmin>grd_list(n)%k0) kmin=grd_list(n)%k0

     if (imax<grd_list(n)%i0+grd_list(n)%ni-1) imax=grd_list(n)%i0+grd_list(n)%ni-1
     if (jmax<grd_list(n)%j0+grd_list(n)%nj-1) jmax=grd_list(n)%j0+grd_list(n)%nj-1
     if (kmax<grd_list(n)%k0+grd_list(n)%nk-1) kmax=grd_list(n)%k0+grd_list(n)%nk-1
  end do
  do n=1,nmax
     grd_list(n)%i0=grd_list(n)%i0-imin+1
     grd_list(n)%j0=grd_list(n)%j0-jmin+1
     grd_list(n)%k0=grd_list(n)%k0-kmin+1
  end do
  imax=imax-imin+1
  jmax=jmax-jmin+1
  kmax=kmax-kmin+1
end subroutine read_grd2nc

function is_in_plane(x0,y0,z0,x,y,z,i,j) result( isIn )
real*8 :: x0,y0,z0
real*8,dimension(:,:) :: x,y,z
logical isIn
integer i,j

integer ni,nj

isIn=.false.
ni=size(x,1); nj=size(x,2)
do j=1,nj
do i=1,ni
   if ( (x(i,j)-x0)**2+(y(i,j)-y0)**2+(z(i,j)-z0)**2 <= SEIS_ZERO ) then
      isIn=.true.
      return
   end if
end do
end do
end function is_in_plane

subroutine destroy_grdlist
  integer n,nmax
  nmax=ubound(grd_list,1)
  do n=1,nmax
     deallocate(grd_list(n)%x)
     deallocate(grd_list(n)%y)
     deallocate(grd_list(n)%z)
  end do
  deallocate(grd_list)
end subroutine destroy_grdlist

subroutine export_grd2nc(fnm_nc,imax,jmax,kmax,i0,j0,k0)
  character (len=*),intent(in) :: fnm_nc
  integer imax,jmax,kmax,i0,j0,k0

  integer n,nmax,ni,nj,nk
  integer ncid,ierr,oldMode,xdimid,ydimid,zdimid,xid,yid,zid

  nmax=ubound(grd_list,1)

  ! creat netcdf file
  ierr=nf90_create(path = trim(fnm_nc),cmode= nf90_clobber, ncid = ncid )
       call nfseis_handle_err(ierr)
  !--
  ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
  ierr=nf90_def_dim(ncid,"I",imax,xdimid)
  ierr=nf90_def_dim(ncid,"J",jmax,ydimid)
  ierr=nf90_def_dim(ncid,"K",kmax,zdimid)
  ! define variable
  ierr=nf90_def_var(ncid,'x',nf90_float,(/ xdimid,ydimid,zdimid /), xid )
       call nfseis_handle_err(ierr)
  ierr=nf90_def_var(ncid,'y',nf90_float,(/ xdimid,ydimid,zdimid /), yid )
       call nfseis_handle_err(ierr)
  ierr=nf90_def_var(ncid,'z',nf90_float,(/ xdimid,ydimid,zdimid /), zid )
       call nfseis_handle_err(ierr)
  ! define att
  ierr=nf90_put_att(ncid,NF90_GLOBAL,'indx',(/ i0,j0,k0 /) )
       call nfseis_handle_err(ierr)
  ! end define
  ierr=nf90_enddef(ncid)
  ! put data
  do n=1,nmax
     ni=grd_list(n)%ni; nj=grd_list(n)%nj; nk=grd_list(n)%nk
     i0=grd_list(n)%i0; j0=grd_list(n)%j0; k0=grd_list(n)%k0
     ierr=nf90_put_var(ncid,xid,grd_list(n)%x(1:ni,1:nj,1:nk),              &
          (/ i0,j0,k0 /), (/ ni,nj,nk /) )
          call nfseis_handle_err(ierr)
     ierr=nf90_put_var(ncid,yid,grd_list(n)%y(1:ni,1:nj,1:nk),              &
          (/ i0,j0,k0 /), (/ ni,nj,nk /) )
          call nfseis_handle_err(ierr)
     ierr=nf90_put_var(ncid,zid,grd_list(n)%z(1:ni,1:nj,1:nk),              &
          (/ i0,j0,k0 /), (/ ni,nj,nk /) )
          call nfseis_handle_err(ierr)
  end do
  !---- close files -----
  ierr=nf90_close(ncid)
end subroutine export_grd2nc

subroutine nfseis_handle_err(ierr)
    integer,intent(in) :: ierr
    if (ierr /= nf90_noerr) then
       print *, trim(nf90_strerror(ierr))
       stop "Stopped"
    end if
end subroutine nfseis_handle_err

end program grd2nc

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
