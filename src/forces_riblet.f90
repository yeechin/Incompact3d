!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!=======================================================================
! This program computes the drag coefficients of a riblet channel. 
!=======================================================================

module forces_riblet
  USE decomp_2d
  implicit none

  integer :: iforces_riblet
!  real(mytype),save,allocatable,dimension(:,:,:) :: ux01, uy01, ux11, uy11, ppi1
!  real(mytype),allocatable,dimension(:) :: xld,xrd,yld,yud
!  integer,allocatable,dimension(:) :: icvlf,icvrt,jcvlw,jcvup
!  integer,allocatable,dimension(:) :: icvlf_lx,icvrt_lx,icvlf_ly,icvrt_ly
!  integer,allocatable,dimension(:) :: jcvlw_lx,jcvup_lx,jcvlw_ly,jcvup_ly

!  character(len=*), parameter :: io_restart_forces = "restart-forces-io", &
!       resfile = "restart-forces"
  
contains

  subroutine init_forces_riblet

    USE decomp_2d
    USE decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    USE param
    USE variables
    implicit none

    integer :: iv,stp1,stp2,h

!    call alloc_x(ux01)
!    call alloc_x(uy01)
!    call alloc_x(ux11)
!    call alloc_x(uy11)
!    call alloc_x(ppi1)

 !   ux01 = zero
!    uy01 = zero
!    ux11 = zero
!    uy11 = zero

!    allocate(icvlf(nvol),icvrt(nvol),jcvlw(nvol),jcvup(nvol))
!    allocate(icvlf_lx(nvol),icvrt_lx(nvol),icvlf_ly(nvol),icvrt_ly(nvol))
!    allocate(jcvlw_lx(nvol),jcvup_lx(nvol),jcvlw_ly(nvol),jcvup_ly(nvol))

    !     Definition of the Control Volume
    !*****************************************************************
    !! xld,xrd,yld,yud: limits of control volume (!!don't use cex and cey anymore!!)

!    do iv=1,nvol
       ! ok for istret=0 (!!to do for istret=1!!)
!       icvlf(iv) = nint(xld(iv)/dx)+1
!       icvrt(iv) = nint(xrd(iv)/dx)+1
!       if (istret.eq.0) then 
!         jcvlw(iv) = nint(yld(iv)/dy)+1
!         jcvup(iv) = nint(yud(iv)/dy)+1
!       else
!         stp1=0
!         stp2=0
!         do h = 1, ny-1  
!           if ((-yp(h+1)-yp(h)+two*yld(iv)).lt.(yld(iv)-yp(h)).and.(stp1.eq.0)) then 
!             jcvlw(iv) = h+1
!             stp1=1
!           endif
!           if ((-yp(h+1)-yp(h)+two*yud(iv)).lt.(yud(iv)-yp(h)).and.(stp2.eq.0)) then
!             jcvup(iv) = h
!             stp2=1 
!           endif
!         enddo
!       endif
!       icvlf_lx(iv) = icvlf(iv)
!       icvrt_lx(iv) = icvrt(iv)
!       jcvlw_lx(iv) = max(jcvlw(iv)+1-xstart(2),1)
!       jcvup_lx(iv) = min(jcvup(iv)+1-xstart(2),xsize(2))
!       icvlf_ly(iv) = max(icvlf(iv)+1-ystart(1),1)
!       icvrt_ly(iv) = min(icvrt(iv)+1-ystart(1),ysize(1))
!       jcvlw_ly(iv) = jcvlw(iv)
!       jcvup_ly(iv) = jcvup(iv)
!    enddo

    if (nrank==0) then
       write(*,*) '========================Forces============================='
       write(*,*) '                        '
       write(*,*) '  Calculate drag coefficient of a riblet surface '
       write(*,*) '                       '

       write(*,*) '==========================================================='
    endif

!    call decomp_2d_init_io(io_restart_forces)
!    call decomp_2d_register_variable(io_restart_forces, "ux01", 1, 0, 0, mytype)
!    call decomp_2d_register_variable(io_restart_forces, "uy01", 1, 0, 0, mytype)
!    call decomp_2d_register_variable(io_restart_forces, "ux11", 1, 0, 0, mytype)
!    call decomp_2d_register_variable(io_restart_forces, "uy11", 1, 0, 0, mytype)
    
  end subroutine init_forces_riblet


  subroutine force_riblet(ux1,ep1)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm_param

    use var, only : ta1, tb1, tc1, td1, di1
    use var, only : ux2, uy2, ta2, tb2, tc2, td2, di2

    implicit none
    character(len=30) :: filename, filename2
!    integer :: nzmsize
!    integer                                             :: i, iv, j, k, kk, code, jj
!    integer                                             :: nvect1,nvect2,nvect3

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ep1

!    real(mytype), dimension(ysize(1),ysize(2),ysize(3)) :: ppi2

!    real(mytype), dimension(nz) :: yLift,xDrag
!    real(mytype) :: yLift_mean,xDrag_mean

!    real(mytype), dimension(ny-1) :: del_y

!    real(mytype), dimension(nz) :: tunstxl, tunstyl
!    real(mytype), dimension(nz) :: tconvxl,tconvyl
!    real(mytype), dimension(nz) :: tpresxl,tpresyl
!    real(mytype), dimension(nz) :: tdiffxl,tdiffyl

!    real(mytype), dimension(nz) :: tunstx, tunsty
!    real(mytype), dimension(nz) :: tconvx,tconvy
!    real(mytype), dimension(nz) :: tpresx,tpresy
!    real(mytype), dimension(nz) :: tdiffx,tdiffy

!    real(mytype) :: uxmid,uymid,prmid
!    real(mytype) :: dudxmid,dudymid,dvdxmid,dvdymid
!    real(mytype) :: fac,tsumx,tsumy
!    real(mytype) :: fcvx,fcvy,fprx,fpry,fdix,fdiy
!    real(mytype) :: xmom,ymom
	
    integer       ::i,j,k,kk
	integer       ::np_riblet
	integer :: code
	real(mytype), dimension(ny*nz/n_riblet) :: umean_riblet !streamwise and riblet-wise average of u
	real(mytype), dimension(ny*nz/n_riblet) :: umean_riblet_recv !receive buffer
	real(mytype), dimension((ny+1)/2,nz/n_riblet/2+1) :: umean_riblet_half
	real(mytype)    ::dudyw_normal_int, dudzw_normal_int
	real(mytype)    ::xDrag_coeff
	
	
	np_riblet=nz/n_riblet
	umean_riblet=zero
	umean_riblet_recv=zero
	do j=1,xsize(2)
		do k=1,xsize(3)
			do i=1,xsize(1)
				kk=xstart(3)-1+k !global index in z-direction
				kk=mod(kk,np_riblet)
				if(kk==0)kk=np_riblet
				umean_riblet(xstart(2)-1+j+(kk-1)*ny)=umean_riblet(xstart(2)-1+j+(kk-1)*ny)+ux1(i,j,k)/nx/n_riblet
			end do
		end do
	end do
	call MPI_ALLREDUCE(umean_riblet,umean_riblet_recv,ny*nz/n_riblet,real_type,MPI_SUM,MPI_COMM_WORLD,code)
	
	umean_riblet_half=zero
	do j=1,(ny+1)/2
		do k=1,np_riblet/2+1
			if(k==1) then
				umean_riblet_half(j,k)=(umean_riblet_recv(j+(k-1)*ny)+umean_riblet_recv(ny+1-j+(k-1)*ny))/two
			else
				umean_riblet_half(j,k)=(umean_riblet_recv(j+(k-1)*ny)+umean_riblet_recv(j+(np_riblet+1-k)*ny) &
			           +umean_riblet_recv(ny+1-j+(k-1)*ny)+umean_riblet_recv(ny+1-j+(np_riblet+1-k)*ny))/four
		    end if
		end do
	end do
  call dudy_wall(umean_riblet_half,dudyw_normal_int)
  call dudz_wall(umean_riblet_half,dudzw_normal_int)
  xDrag_coeff = two*xnu*(dudyw_normal_int+dudzw_normal_int)/(two/three)**2
  
  if(nrank.eq.0) then
!	  write(*,*) t,xDrag_coeff
	  write(138,*)t,xDrag_coeff
	  call flush(138)
  end if
  if(mod(itime,icheckpoint).eq.0)then
	  if(nrank.eq.0) then
		  write(filename,"('forces.dat',I7.7)") itime
		  call system("cp forces.dat " //filename)
	  end if
  end if
  return
  end subroutine force_riblet

  subroutine dudy_wall(u,dudyw_normal_int)
	  USE param
	  USE complex_geometry
	  USE decomp_2d
	  USE variables
	  USE ibm, only : polint
	  
	  implicit none
	  
	  real(mytype),dimension((ny+1)/2,nz/n_riblet/2+1)  ::u
	  integer                                             :: i, j, k
	  integer                                             :: jy
	  
	  ! The size of dudyw is specific for riblet flow
	  real(mytype),dimension(nz/n_riblet/2+1)         ::dudyw
	  real(mytype)                                    ::dudyw_normal_int
	  
	  real(mytype),dimension(2)                           ::xpol,ypol
	  real(mytype)                                        ::dypol, deltay
	  
	  real(mytype),dimension(10)                          :: xa, ya
	  integer                                             :: ia, na, np_half_riblet, jpif, nypif
!	  real(mytype)               :: zeromach
	  
!      zeromach=one
!      do while ((one + zeromach / two) .gt. one)
!         zeromach = zeromach/two
!      end do
!      zeromach = ten*zeromach
	  
	  np_half_riblet=nz/n_riblet/2+1
	  
	  do k=1,np_half_riblet
		  ia=0
		  nypif=npif
		  ia=ia+1
		  xa(ia)=yf(1,1,k)
		  ya(ia)=zero
		  jy=1
		  do while(yp(jy).lt.yf(1,1,k))
			  jy=jy+1
		  end do
		  if(nyfpif(1,1,k).lt.npif) nypif=nyfpif(1,1,k)
		  do jpif=1,nypif
			  ia=ia+1
			  if(izap.eq.1) then !zapping
				  xa(ia)=yp(jy+jpif)
				  ya(ia)=u(jy+jpif,k)
			  else !no zapping
				  xa(ia)=yp(jy+jpif-1)
				  ya(ia)=u(jy+jpif-1,k)
			  end if
		  end do
		  na = ia
		  deltay=yp(2)-yp(1)
		  xpol(1)=yf(1,1,k)-deltay/ten
		  xpol(2)=yf(1,1,k)+deltay/ten
		  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
		  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
	      dudyw(k)=(ypol(2)-ypol(1))/deltay*five
		  
!		  call normal_riblet((k-1)*dz,normaly,normalz)
!		  if(normaly.lt.zeromach)then
!			  dudyw(k)=zero
!		  else
!			  dudyw(k)=dudyw(k)*normaly
!		      dudyw(k)=dudyw(k)*sqrt(1+(normalz/normaly)**2)
		      if(k==1.or.k==np_half_riblet) dudyw(k)=dudyw(k)*zpfive
	  end do
	  
	  dudyw_normal_int=sum(dudyw)/real(np_half_riblet-1)
	  return
  end subroutine dudy_wall


  ! The subroutine computes dudz on the upper and lower riblet surfaces in a channel riblet flow  
    subroutine dudz_wall(u,dudzw_normal_int)
  	  USE param
  	  USE complex_geometry
  	  USE decomp_2d
  	  USE variables
	  USE ibm, only : polint
	  
  	  implicit none
	  
  	  real(mytype),dimension((ny+1)/2,nz/n_riblet/2+1)  ::u
  	  integer                                             :: i, j, k
  	  integer                                             :: kz
	  
  	  ! The size of dudyw is specific for riblet flow
  	  real(mytype),dimension((ny+1)/2) ::dudzw
	  
  	  real(mytype),dimension(2)                           ::xpol,ypol
  	  real(mytype)                                        ::dypol,deltaz
	  
  	  real(mytype),dimension(10)                          :: xa, ya
  	  integer                                             :: ia, na, kpif,nzpif
!	  real(mytype)                                        :: normaly,normalz
	  real(mytype)                                        :: W_riblet, dudzw_normal_int
	  
!	  real(mytype)               :: zeromach
	  
!      zeromach=one
!      do while ((one + zeromach / two) .gt. one)
!         zeromach = zeromach/two
!      end do
!      zeromach = ten*zeromach
       W_riblet = zlz/n_riblet ! Width of riblet
	  
	  dudzw=zero

	  do j=2,(ny+1)/2
		  if(nobjz(1,j).ne.0) then
			  ia=0
			  nzpif=npif
			  ia=ia+1
			  xa(ia)=zf(1,1,j)
			  ya(ia)=zero
			  
			  kz=(zf(1,1,j)+dz)/dz+1
			  if(nzfpif(1,1,j).lt.npif) nzpif=nzfpif(1,1,j)
			  do kpif=1,nzpif
				  ia=ia+1
				  if(izap.eq.1) then
					  xa(ia)=(kz-1)*dz+kpif*dz
					  ya(ia)=u(j,kz+kpif)
				  else
					  xa(ia)=(kz-1)*dz+(kpif-1)*dz
					  ya(ia)=u(j,kz+kpif-1)
				  end if
			  end do
			  na = ia
			  deltaz=dz
			  xpol(1)=zf(1,1,j)-deltaz/ten
			  xpol(2)=zf(1,1,j)+deltaz/ten
			  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
			  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
			  dudzw(j)=(ypol(2)-ypol(1))/deltaz*five
			  
!			  call normal_riblet(zf(1,1,j),normaly,normalz)
			  
!			  if(normalz.lt.zeromach) then
!				  dudzw(j)=zero
!			  else
!				  dudzw(j)=dudzw(j)*normalz	  
		  end if
	  end do
	  dudzw(1)=zero
	  
	  dudzw_normal_int=zero
	  do j=2,(ny+1)/2-1
		  dudzw_normal_int=dudzw_normal_int+dudzw(j)*(yp(j+1)-yp(j-1))*zpfive
	  end do
	  dudzw_normal_int=dudzw_normal_int/W_riblet*two
  	  return
    end subroutine dudz_wall 



  ! get the normal direction on the riblet surface, only for the first half riblet
 ! subroutine normal_riblet(z,normaly,normalz)
!	  USE param
!	  USE complex_geometry
!	  USE decomp_2d
!	  USE variables
	  
!	  implicit none
!	  real(mytype)          ::z,normaly,normalz
!	  real(mytype)          ::W_riblet,zm
!	  real(mytype)          ::dydz
!	  real(mytype)          ::temp
	  
!  	if(A_riblet.eq.six.and.B_riblet.eq.six) then ! Perfect balanced polynomials
!  		zstar = zpfive
!  		a = eight
!  		b = eight
!  	else
!  		zstar = (three-zpfive*B_riblet)/(A_riblet-B_riblet)
!  		a = (two*(A_riblet+B_riblet)*zstar+A_riblet-B_riblet)/three/zstar
!  		b = ((A_riblet+B_riblet)*zstar-b_riblet)/onepfive/(zstar-zpfive)
!  	endif

!	W_riblet = zlz/n_riblet ! Width of riblet
!	zm = z/W_riblet !Non-dimensionalized zm
	
!	if(zm.le.zstar) then !Belongs to the first polynomial
!		dydz = three*a*zm**2 - two*A_riblet*zm
!	else
!		dydz = three*b*(zm-zpfive)**2 + two*B_riblet*(zm-zpfive)
!	endif
!	dydz = two*gamma_riblet*dydz
	
!	normaly=1
!	normalz=-dydz
!	
!	temp = sqrt(normaly**2+normalz**2)
!	normaly=normaly/temp
!	normalz=normalz/temp ! unit normal vector pointing to fluid
!  end subroutine normal_riblet
!  
!
!
!
  ! The following subroutines are for the whole riblet shape, however, currently, we only need to 
  ! calculate derivatives on the averaged half riblet
! The subroutine computes dudy on the upper and lower riblet surfaces in a channel riblet flow
!  subroutine dudy_wall(u,dudyw)
!	  USE param
!	  USE complex_geometry
!	  USE decomp_2d
!	  USE variables
!	  
!	  implicit none
!	  
!	  real(mytype),dimension(ysize(1),ysize(2),ysize(3))  ::u
!	  integer                                             :: i, j, k
!	  integer                                             :: jy
!	  
!	  ! The size of dudyw is specific for riblet flow
!	  real(mytype),dimension(2,ysize(1),ysize(3))         ::dudyw
!	  
!	  real(mytype),dimension(2)                           ::xpol,ypol
!	  real(mytype)                                        ::dypol
!	  
!	  real(mytype),dimension(10)                          :: xa, ya
!	  integer                                             :: ia, na
!	  integer                                             :: remp ! 0 for lower riblet, 1 for upper riblet
!	  
!	  do k=1,ysize(3)
!		  do i=1,ysize(1)
!			  if(nobjy(i,k).ne.0) then
!				  ia=0
!				  do j=1,nobjy(i,k) ! for channel riblet, nobjy=2
!					  !The first boundary
!					  nypif = npif
!					  ia=ia+1
!					  xa(ia)=yi(j,i,k)
!					  ya(ia)=0.
!					  if(yi(j,i,k).gt.0.)then ! Object immersed
!						  jy=1
!						  do while(yp(jy).lt.yi(j,i,k))
!							  jy=jy+1
!						  end do
!						  jy=jy-1
!						  jpoli=jy+1
!						  if(nyipif(j,i,k).lt.npif) nypif=nyipif(j,i,k)
!						  do jpif=1,nypif
!							  ia=ia+1
!							  if(izap.eq.1)then !zapping
!								  xa(ia)=yp(jy-jpif)
!								  ya(ia)=u(i,jy-jpif,k)
!						      else ! no zapping
!								  xa(ia)=yp(jy-jpif+1)
!								  ya(ia)=u(i,jy-jpif+1,k)
!							  end if
!						  end do
!					  else
!						  jpoli=1
!					  end if
!					  
!					  !The second boundary
!					  nypif=npif
!					  ia=ia+1
!					  xa(ia)=yf(j,i,k)
!					  ya(ia)=0
!					  if(yf(j,i,k).lt.yly) then !Object immersed
!						  jy=1
!						  do while(yp(jy).lt.yf(j,i,k))
!							  jy=jy+1
!						  end do
!						  jpolf=jy-1
!						  if(nyfpif(j,i,k).lt.npif) nypif=nyfpif(j,i,k)
!						  do jpif=1,nypif
!							  ia=ia+1
!							  if(izap.eq.1) then!zapping
!								  xa(ia)=yp(jy+jpif)
!								  ya(ia)=u(i,jy+jpif,k)
!							  else ! no zapping
!									  xa(ia)=yp(jy+jpif-1)
!									  ya(ia)=u(i,jy+jpif-1,k)
!							  end if
!						  end do
!					  else
!						  jpolf=ny
!					  end if
!					  na=ia
!					  ! specific for riblet flow, the lower boundary and the upper boundary
!					  if(j==1) then
!						  xpol(1)=yf(j,i,k)-dy/ten
!						  xpol(2)=yf(j,i,k)+dy/ten
!						  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
!						  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
!					      dudyw(j,i,k)=(ypol(2)-ypol(1))/dy*five
!					  else if(j==2) then
!						  xpol(1)=yi(j,i,k)-dy/ten
!						  xpol(2)=yi(j,i,k)+dy/ten
!						  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
!						  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
!					      dudyw(j,i,k)=(ypol(2)-ypol(1))/dy*five
!					  end if
!					  ia=0
!				  end do
!			  end if
!		  end do
!	  end do
!	  
!	  return
!  end subroutine dudy_wall 

! The subroutine computes dudz on the upper and lower riblet surfaces in a channel riblet flow  
!  subroutine dudz_wall(u,dudzw)
!	  USE param
!	  USE complex_geometry
!	  USE decomp_2d
!	  USE variables
!	  
!	  implicit none
!	  
!	  real(mytype),dimension(ysize(1),ysize(2),ysize(3))  ::u
!	  integer                                             :: i, j, k
!	  integer                                             :: kz
!	  
!	  ! The size of dudyw is specific for riblet flow
!	  real(mytype),dimension(nobjmax*2,zsize(1),ysize(2)) ::dudzw
!	  
!	  real(mytype),dimension(2)                           ::xpol,ypol
!	  real(mytype)                                        ::dypol
!	  
!	  real(mytype),dimension(10)                          :: xa, ya
!	  integer                                             :: ia, na
!	  integer                                             :: remp ! 0 for lower riblet, 1 for upper riblet
!	  
!	  do j=1,zsize(2)
!		  do i=1,zsize(1)
!			  if(nobjz(i,j).ne.0) then
!				  remp=0
!				  if(j>(ny+1)/2) remp=1
!				  ia=0
!				  do k=1,nobjz(i,j) ! for channel riblet, nobjy=n_riblet+1 because of riblet tips
!					  !are set at the spanwise boundary
!					  
!					  !The first boundary
!					  nzpif = npif
!					  ia=ia+1
!					  xa(ia)=zi(k,i,j)
!					  ya(ia)=0.
!					  if(zi(j,i,k).gt.0.)then ! Object immersed
!						  kz=zi(k,i,j)/dz+1
!						  kpoli=kz+1
!						  if(nzipif(k,i,j).lt.npif) nzpif=nzipif(k,i,j)
!						  do kpif=1,nzpif
!							  ia=ia+1
!							  if(izap.eq.1)then !zapping
!								  xa(ia)=(kz-1)*dz-kpif*dz
!								  ya(ia)=u(i,j,kz-kpif)
!					      else ! no zapping
!							  xa(ia)=(kz-1)*dz-(kpif-1)*dz
!							  ya(ia)=u(i,j,kz-kpif+1)
!						  end if
!					  end do
!				  else
!					  kpoli=1
!				  end if
!				  
!				  !The second boundary
!				  nzpif=npif
!				  ia=ia+1
!				  xa(ia)=zf(k,i,j)
!				  ya(ia)=0
!				  if(zf(k,i,j).lt.zlz) then !Object immersed
!					  kz=(zf(k,i,j)+dz)/dz+1
!					  kpolf=kz-1
!					  if(nzfpif(k,i,j).lt.npif) nzpif=nzfpif(k,i,j)
!					  do kpif=1,nzpif
!						  ia=ia+1
!						  if(izap.eq.1) then!zapping
!							  xa(ia)=(kz-1)*dz+kpif*dz
!							  ya(ia)=u(i,j,kz+kpif)
!						  else ! no zapping
!							  xa(ia)=(kz-1)*dz+(kpif-1)*dz
!							  ya(ia)=u(i,j,kz+kpif-1)
!						  end if
!					  end do
!				  else
!					  kpolf=nz
!				  end if
!				  na=ia
!				  ! specific for riblet flow
!				  
!				  if(k==1) then
!					  xpol(1)=zf(k,i,j)-dz/ten
!					  xpol(2)=zf(k,i,j)+dz/ten
!					  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
!					  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
!				      dudzw(k,i,j)=(ypol(2)-ypol(1))/dz*five
!				  else if(k==nobjz(i,j)) then
!					  xpol(1)=zi(k,i,j)-dz/ten
!					  xpol(2)=zi(k,i,j)+dz/ten
!					  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
!					  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
!				      dudzw(2*(nobjz(i,j)-1),i,j)=(ypol(2)-ypol(1))/dz*five
!				  else
!					  xpol(1)=zi(k,i,j)-dz/ten
!					  xpol(2)=zi(k,i,j)+dz/ten
!					  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
!					  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
!				      dudzw(2*(k-1),i,j)=(ypol(2)-ypol(1))/dz*five
!					  
!					  xpol(1)=zf(k,i,j)-dz/ten
!					  xpol(2)=zf(k,i,j)+dz/ten
!					  call polint(xa,ya,na,xpol(1),ypol(1),dypol)
!					  call polint(xa,ya,na,xpol(2),ypol(2),dypol)
!				      dudzw(2*k-1,i,j)=(ypol(2)-ypol(1))/dz*five 
!				  end if
!				  ia=0
!			  end do
!		  end if
!	  end do
!  end do
!  
!  return
!end subroutine dudz_wall 
!
!! get the normal direction on the riblet surface
!subroutine normal_riblet(z,normaly,normalz,remp)
!  USE param
!  USE complex_geometry
!  USE decomp_2d
!  USE variables
!  
!  implicit none
!  real(mytype)          ::z,normaly,normalz
!  real(mytype)          ::W_riblet,zm
!  real(mytype)          ::dydz
!  real(mytype)          ::temp
!  integer              :: remp !0 for the lower wall, 1 for the upper wall
!  integer              :: iflag ! 0 for \hat{z}*<=0.5; 1 for \hat{z}>0.5
!  
!	if(A_riblet.eq.six.and.B_riblet.eq.six) then ! Perfect balanced polynomials
!		zstar = zpfive
!		a = eight
!		b = eight
!	else
!		zstar = (three-zpfive*B_riblet)/(A_riblet-B_riblet)
!		a = (two*(A_riblet+B_riblet)*zstar+A_riblet-B_riblet)/three/zstar
!		b = ((A_riblet+B_riblet)*zstar-b_riblet)/onepfive/(zstar-zpfive)
!	endif
!
!iflag = 0
!
!W_riblet = zlz/n_riblet ! Width of riblet
!zm = mod(zm, W_riblet)/W_riblet !Non-dimensionalized zm
!if(zm>zpfive) then
!	zm = one - zm
!	iflag = 1
!end if
!if(zm.le.zstar) then !Belongs to the first polynomial
!	dydz = three*a*zm**2 - two*A_riblet*zm
!else
!	dydz = three*b*(zm-zpfive)**2 + two*B_riblet*(zm-zpfive)
!endif
!dydz = two*gamma_riblet*dydz
!
!normaly=1
!normalz=-dydz
!
!if(iflag==1) normalz=-normalz! right half of the riblet
!if(remp==1) normaly=-normaly !upper riblet surface
!
!temp = sqrt(normaly**2+normalz**2)
!normaly=normaly/temp
!normalz=normalz/temp ! unit normal vector pointing to fluid
!end subroutine normal_riblet
!
end module forces_riblet
