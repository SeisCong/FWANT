module media_mod

! This module contains the variables for 3D medium
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-16 12:39:55 -0500 (Fri, 16 Jan 2009) $
! $Revision: 60 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

use constants_mod
use string_mod
use math_mod
use para_mod
use mpi_mod
use nfseis_mod

implicit none
private
public ::         &
  media_fnm_init, &
  media_fnm_get,  &
  media_destroy,  &
  media_alloc,    &
  media_import

!-----------------------------------------------------------------------------

real(SP),dimension(:,:,:),allocatable,public :: rho,mu,lambda
real(SP),dimension(:,:,:),allocatable,public :: Qs
real(SP),public :: QsF0,QsINF
character (len=SEIS_STRLEN),public ::       &
     fnm_media_conf, pnm_media
integer :: ierr

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

!*************************************************************************
!*                            alloc and dealloc                          *
!*************************************************************************
subroutine media_alloc
  allocate( mu(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  mu=0.0_SP
  allocate( lambda(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  lambda=0.0_SP
  allocate( rho(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  rho=0.0_SP
#ifdef WITHQS
  allocate( Qs(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Qs=0.0_SP
#endif
end subroutine media_alloc
subroutine media_destroy
  if (allocated(mu)) deallocate( mu )
  if (allocated(lambda)) deallocate( lambda )
  if (allocated(rho)) deallocate( rho )
  if (allocated(Qs)) deallocate( Qs )
end subroutine media_destroy

!*************************************************************************
!*                             media io                                  *
!*************************************************************************
subroutine media_fnm_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
    call string_conf(fid,1,'MEDIA_CONF',2,fnm_media_conf)
    call string_conf(fid,1,'MEDIA_ROOT',2,pnm_media)
  close(fid)
end subroutine media_fnm_init

function media_fnm_get(n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_media)//'/'//'media'//'_'//set_mpi_subfix(n_i,n_j,n_k)//'.nc'
end function media_fnm_get

subroutine media_import(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  integer,dimension(SEIS_GEO) :: subs,subc,subt
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt)
  call nfseis_varget( filenm, 'rho', rho, subs,subc,subt)
#ifdef WITHQS
  call nfseis_varget( filenm, 'Qs', Qs, subs,subc,subt)
  call nfseis_attget( filenm, 'QsF0', QsF0)
  call nfseis_attget( filenm, 'QsINF', QsINF)
#endif
end subroutine media_import
  
end module media_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
