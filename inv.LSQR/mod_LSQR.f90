module LSQR_mod
!----------------------------------------------------------------
!                    LSQR package
!----------------------------------------------------------------
! call lsqr_init
! call lldrow (coef,idx,ncoef,rhs)
! call slsqr(nit)
! call lgtsol(dd,errx,nnn)

! common storage for coefficient matrix
!     row order with coefficients, col index, and number
!     in row
! definition of variables
!     mrow,ncol = dimensions of system of equations (mrow by ncol)
!     nel = total number of active elements in A
!     ra = coefficient matrix values
!     ja = column indices
!     na = array of number of active elements per row

implicit none
private

public ::       &
  lsqr_init,    &
  lsqr_destroy, &
  lldrow,       &
  slsqr,        &
  lgtsol,       &
  lsqr_check,   &
  error_except

integer,parameter,public :: DP = KIND(1.0D0)
integer,parameter,public :: SP = KIND(1.0)
real(DP),parameter,public :: CONS_EQ=1.0e-20_DP
integer,parameter,public :: INT_KIND=8; ! integer length

!integer :: ncol,max_row !,max_nel
integer(kind=INT_KIND)::max_nel,nel,ncol,max_row,mrow
!integer :: mrow !,nel
real(DP) :: res2

integer(kind=INT_KIND),allocatable :: ja(:),na(:)
real(DP),dimension(:),allocatable ::  &
  ra,x,b,u,v,w, q,sig

!----------------------------------------------------------------
contains
!----------------------------------------------------------------

! ... initialize for conjugate gradient least squares
subroutine lsqr_init(m,n,siz,num)
  integer,intent(in) :: m,n,siz !,num
  integer(kind=INT_KIND),intent(in):: num
  integer :: ierr

  mrow=0; ncol=n; nel=0
  max_row=m; max_nel=num;

  allocate(na(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: na,no enough memory space")
  na=0
  allocate(b(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: b,no enough memory space")
  b=0.0_DP
  allocate(u(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: u,no enough memory space")
  u=0.0_DP
  allocate(q(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: q,no enough memory space")
  q=0.0_DP

  allocate(v(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: v,no enough memory space")
  v=0.0_DP
  allocate(w(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: w,no enough memory space")
  w=0.0_DP
  allocate(x(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: x,no enough memory space")
  x=0.0_DP
  allocate(sig(siz),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: sig,no enough memory space")
  sig=0.0_DP

  allocate(ja(max_nel),stat=ierr); 
  if (ierr/=0) call error_except("lsqr_init: ja,no enough memory space")
  ja=0
  allocate(ra(max_nel),stat=ierr);
  if (ierr/=0) call error_except("lsqr_init: ra,no enough memory space")
  ra=0.0_DP
end subroutine lsqr_init

subroutine lsqr_destroy
  integer(kind=INT_KIND) :: ierr,i

  if (allocated(na)) deallocate(na);
  if (allocated(b)) deallocate(b);
  if (allocated(u)) deallocate(u);
  if (allocated(q)) deallocate(q);

  if (allocated(v)) deallocate(v);
  if (allocated(w)) deallocate(w);
  if (allocated(x)) deallocate(x);
  if (allocated(sig)) deallocate(sig);

  if (allocated(ja)) deallocate(ja)
  if (allocated(ra)) deallocate(ra)
end subroutine lsqr_destroy

! ... add one equation of data into row of matrix
!     for conjugate gradient method (slsqr) of solution
subroutine lldrow(coef,idx,ncoef,rhs)
integer,intent(in) :: idx(:)
integer,intent(in) ::ncoef
real(DP),intent(in) :: coef(:),rhs
integer(kind=INT_KIND) :: i

if (ncoef<=0) return

if (nel+ncoef>max_nel) then
   write(*,*) mrow+1, nel, ncoef, max_nel
   call error_except("exceed max_nel")
end if
mrow=mrow+1
b(mrow) = rhs
na(mrow)= ncoef

ra(nel+1:nel+ncoef)=coef(1:ncoef)
ja(nel+1:nel+ncoef)=idx(1:ncoef)
nel=nel+ncoef
!debug
! write(*,*) 'nel in lldrow, maxnel: ',nel,max_nel
end subroutine lldrow

subroutine lsqr_check(msg)
character (len=*),intent(in) :: msg
write(*,*) msg
write(*,*) "    max_row,mrow=",max_row,mrow
write(*,*) "    max_nel,nel=",max_nel,nel
end subroutine lsqr_check

subroutine slsqr(itct)
! ... conjugate gradient routine for non-square least squares
!     formal parameters:
!       atol = ATOL from Paige and Saunders
!       btol = BTOL           "
!       conlim = CONLIM       "
!       itct = iteration count (returned)

!      parameter(atol=1.d-6,btol=1.d-6,conlim=0.d0)
! Y.Shen, 08/23/04
!parameter(atol=1.d-7,btol=1.d-7,conlim=0.d0)

!real(DP),parameter :: &
!  atol=1.0e-20_DP,     &
!  btol=1.0e-20_DP,     &
!  conlim=0.0_DP
real(DP),parameter :: &
  atol=1.0e-7_DP,     &
  btol=1.0e-7_DP,     &
  conlim=0.0_DP
integer(kind=INT_KIND),parameter :: maxit=30000
integer(kind=INT_KIND),intent(out) :: itct

integer(kind=INT_KIND) :: i
real(DP) :: alpha,beta,phibar,rhobar,phi,rho,c,s,theta
real(DP) :: temp,anorm,test1,test2,rnorm

!     maxit should be set approximatly 4*ncol

! ... initial housekeeping
      itct = 0
      anorm = dotp(ra,ra,nel)
      anorm = dsqrt(anorm)

! ... set initial values for start of iterations
      write(*,*) ncol,mrow,nel
      do i=1,mrow
         u(i) = b(i)
      enddo
      call lnrliz(beta,u,mrow)

! ... compute A transpose u and put results into v
      call laty(v,u)
      call lnrliz(alpha,v,ncol)

! ... move v into w
      do i=1,ncol
         w(i) = v(i)
      enddo

! ... zero out x and sig vectors
      do i=1,ncol
         x(i) = 0.0_DP
         sig(i) = 0.0_DP
      enddo

! ... assign phibar and rhobar
      phibar = beta
      rhobar = alpha

!****************************************
!           main iteration loop         *
!****************************************
loop_iteration: do

! ... compute Av, subtract alpha*u, and assign to u
      call lax(q,v)
      !call lupdat(u,q,'-',alpha,u,mrow)
      call lupdat(u,q,'-',alpha,mrow)
      call lnrliz(beta,u,mrow)

! ... compute A transpose u, subtract beta*v, and assign to v
      call laty(q,u)
      !call lupdat(v,q,'-',beta,v,ncol)
      call lupdat(v,q,'-',beta,ncol)
      call lnrliz(alpha,v,ncol)

! ... do steps (a) through (g) in Paige and Saunders
      rho = dsqrt(rhobar**2+beta**2)

      if(rho.eq.0.d0) return
      !if(abs(rho) <= CONS_EQ) return

      c = rhobar/rho
      s = beta/rho
      theta = s*alpha
      rhobar = -c*alpha
      phi = c*phibar
      phibar = s*phibar

! ... update solution vector and w vector
      do i=1,ncol
         temp = w(i)/rho
         x(i) = x(i)+phi*temp
         w(i) = v(i)-theta*temp
         sig(i) = sig(i)+temp**2
      enddo

! ... loop back if convergence not attained
! ... test for convergence
!
!       calculate r (tmp1)
!
!      call lax(tmp1,x)
!      call lupdat(tmp1,b,'-',1.0,tmp1,mrow)
!      rnorm=dsqrt(dotp(tmp1,tmp1,ncol))
      rnorm=phibar

!     criteria 1 in P&S p.54
      test1=rnorm/anorm/dsqrt(dotp(x,x,ncol))

! ... A(transpose)r (tmp2)
!
!      call laty(tmp2,tmp1)
!
!     criteria 2 in P&S p.54
!
!
!      test2=dsqrt(dotp(tmp2,tmp2,ncol))/anorm/rnorm
      test2=alpha*dabs(c)/anorm

      itct = itct+1
!      write(*,*)itct,test1,atol,test2
      if (mod(itct,6000).eq.0)then
         write(*,'(19hiter., test1,test2:,4i11,3d13.4)') &
                  itct,mrow,ncol,nel,test1,test2,phibar
      end if
      if ((test1.le.atol).or.(test2.le.btol).or.(itct.eq.maxit)) then
         write (*,'(19hiter., test1,test2:,4i11,3d13.4)') &
                 itct,mrow,ncol,nel,test1,test2,phibar
         res2 = phibar
         return
      endif
end do loop_iteration

end subroutine slsqr

subroutine lgtsol(solvec,serr,nvar)
! ... extract solution vector from conjugate gradient least squares
! ... formal parameters:
!       solvec = solution vector
!       serr = standard error estimates
  real(DP) :: solvec(:),serr(:)
  integer :: nvar

  real(DP) :: temp
  integer(kind=INT_KIND) :: i

  temp = res2/dble(max0(mrow-nvar,1))
  do i=1,nvar
     solvec(i) = x(i)
     serr(i) = dsqrt(temp*sig(i))
  enddo
end subroutine lgtsol

! ... normalize the vector 'vect' to unit magnitude and
!     return the scale factor in 'v0'
subroutine lnrliz(v0,vect,n)
integer(kind=INT_KIND),intent(in) :: n
real(DP) :: v0,vect(:)

integer(kind=INT_KIND) :: i

v0 = dotp(vect,vect,n)
v0 = dsqrt(v0)

if (abs(v0) <= CONS_EQ) return

do i=1,n
   vect(i) = vect(i)/v0
enddo
end subroutine lnrliz

!subroutine lupdat(a,b,flag,scal,c,n)
!! ... update vector A by adding scaled version of vector B
!!     character flag indicates addition or subtraction
!  real(DP),intent(out) :: a(:)
!  real(DP),intent(in) :: b(:),c(:),scal
!  integer,intent(in) :: n
!  character*2 flag
!
!  integer :: i
!  if (flag(1:1).eq.'+') then
!     do i=1,n
!        a(i) = b(i)+c(i)*scal
!     enddo
!  else
!     do i=1,n
!        a(i) = b(i)-c(i)*scal
!     enddo
!  end if
!end subroutine lupdat
subroutine lupdat(a,b,flag,scal,n)
! ... update vector A by scaled then adding vector B
!     character flag indicates addition or subtraction
  real(DP),intent(inout) :: a(:)
  real(DP),intent(in) :: b(:),scal
  integer(kind=INT_KIND),intent(in) :: n
  character (len=*) :: flag

  integer(kind=INT_KIND) :: i
  if (flag(1:1).eq.'+') then
     do i=1,n
        a(i) = b(i)+a(i)*scal
     enddo
  else
     do i=1,n
        a(i) = b(i)-a(i)*scal
     enddo
  end if
end subroutine lupdat

function dotp(a,b,n) result(cov2)
  integer(kind=INT_KIND),intent(in) :: n
  real(DP),intent(in) :: a(:),b(:)
  real(DP) :: cov2

  integer(kind=INT_KIND) :: i
  cov2=0.0_DP
  do i=1,n
     cov2=cov2+a(i)*b(i)
  enddo
end function dotp

subroutine lax(yy,xx)
! ... Form matrix product AX and put result in vector Y
!     X is input in the dummy array "xx" which is assumed
!c       to  be of length "ncol";  the output Y is returned in
!       the dummy array "yy" of length "mrow"
  real(DP) :: xx(:),yy(:),y0
  integer(kind=INT_KIND) :: i,j,l,l1,l2
  l2 = 0
  do i=1,mrow
     y0 = 0.d0
     l1 = l2+1; l2 = l2+na(i)
     do l=l1,l2
        j = ja(l)
        y0 = y0+ra(l)*xx(j)
     enddo
     yy(i) = y0
  enddo
end subroutine lax

subroutine laty(xx,yy)
! ... form matrix product A(transpose)Y and put result in vector X
!     The input vector Y of length "mrow" is passed in the dummy
!       array "yy" and the results are returned in the dummy array "xx"
!       which is of length "ncol"
  real(DP) :: xx(:),yy(:),yi
  integer(kind=INT_KIND) :: i,j,l,l1,l2
  l2 = 0
  do i=1,ncol
     xx(i) = 0.0_DP
  enddo
  do i=1,mrow
     yi = yy(i)
     l1 = l2+1; l2 = l2+na(i)
     do l=l1,l2
        j = ja(l)
        xx(j) = xx(j)+ra(l)*yi
     enddo
  enddo
end subroutine laty

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
  print *, trim(msg)
  stop 1
end subroutine error_except

end module LSQR_mod

