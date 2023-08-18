!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module channel

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_channel, boundary_conditions_channel, postprocess_channel, &
            visu_channel, visu_channel_init, momentum_forcing_channel, &
            geomcomplex_channel, critR

contains
  !############################################################################
  subroutine init_channel (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    real(mytype) :: ftent
    integer :: k,j,i,fh,ierror,ii,is,it,code, jj

    if (idir_stream /= 1 .and. idir_stream /= 3) then
       if (nrank == 0) then
          write(*,*) '!! ERROR in imposing sorce term for momentum !!'
          write(*,*) '!! idir_stream ', idir_stream
          write(*,*) '!! idir_stream has to be:'
          write(*,*) '!! - 1 for streamwise direction in X'
          write(*,*) '!! - 3 for streamwise direction in Z'
          write(*,*) '!! Y is not supported and other values do not make sense'
          write(*,*) '!! Calculation will be now stop'
        endif
        call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
    endif

    if (iscalar==1) then
       if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) then
          write(*,*) 'Imposing linear temperature profile'
       end if
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
             if (istret/=0) y=yp(j+xstart(2)-1)
             do i=1,xsize(1)
                phi1(i,j,k,:) = one - y/yly
             enddo
          enddo
       enddo

       phi1(:,:,:,:) = zero !change as much as you want
       if ((nclyS1 == 2).and.(xstart(2) == 1)) then
         !! Generate a hot patch on bottom boundary
         phi1(:,1,:,:) = one
       endif
       if ((nclySn == 2).and.(xend(2) == ny)) then
         phi1(:,xsize(2),:,:) = zero
       endif
    endif
!
    ux1=zero
    uy1=zero
    uz1=zero
    byx1=zero;byy1=zero;byz1=zero
    ! if to decide type of initialization to apply 
    if (iin == 0) then ! laminar flow
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
             if (istret/=0) y=yp(j+xstart(2)-1)-yly*half
             um=exp(-zptwo*y*y)
             do i=1,xsize(1)
                if (idir_stream == 1) then
                   ux1(i,j,k)=one-y*y
                   uy1(i,j,k)=zero
                   uz1(i,j,k)=sin(real(i-1,mytype)*dx)+cos(real(k-1,mytype)*dz)
                else
                   uz1(i,j,k)=one-y*y
                   uy1(i,j,k)=zero
                   ux1(i,j,k)=zero
                endif
             enddo
          enddo
       enddo     
    elseif (iin <= 2) then ! Traditional init to turbulent flows using random numbers + lam profile
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
       !modulation of the random noise + initial velocity profile
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
             if (istret/=0) y=yp(j+xstart(2)-1)-yly*half
             um=exp(-zptwo*y*y)
             do i=1,xsize(1)
                if (idir_stream == 1) then
                   ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+one-y*y
                   uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
                   uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
                else
                   uz1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+one-y*y
                   uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
                   ux1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)
                endif
             enddo
          enddo
       enddo
    elseif (iin == 4) then ! SEM
       call sem_init_channel(ux1, uy1, uz1)
    endif
   
    !INIT FOR G AND U=MEAN FLOW + NOISE 
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
             uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
             uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
          enddo
       enddo
    enddo

    return
  end subroutine init_channel
  !############################################################################
  !############################################################################
  subroutine boundary_conditions_channel (ux,uy,uz,phi)

    use param
    use var, only : di2
    use variables
    use decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    if (.not. cpg ) then ! if not constant pressure gradient
       if (idir_stream == 1) then
          call channel_cfr(ux,two/three)
       else
          call channel_cfr(uz,two/three)
       endif
    end if

    if (iscalar /= 0) then
       if (iimplicit <= 0) then
          if ((nclyS1 == 2).and.(xstart(2) == 1)) then
             !! Generate a hot patch on bottom boundary
             phi(:,1,:,:) = one
          endif
          if ((nclySn == 2).and.(xend(2) == ny)) THEN
             phi(:,xsize(2),:,:) = zero
          endif
       else
          !
          ! Implicit boundary conditions are usually given in input file
          ! It is possible to modify g_sc here
          ! It is not possible to modify alpha_sc and beta_sc here
          !
          ! Bottom temperature if alpha_sc(:,1)=1 and beta_sc(:,1)=0 (default)
          !if (nclyS1.eq.2) g_sc(:,1) = one
          ! Top temperature if alpha_sc(:,2)=1 and beta_sc(:,2)=0 (default)
          !if (nclySn.eq.2) g_sc(:,2) = zero
       endif
    endif

  end subroutine boundary_conditions_channel
  !############################################################################
  !!
  !!  SUBROUTINE: channel_cfr
  !!      AUTHOR: Kay Schäfer
  !! DESCRIPTION: Inforces constant flow rate without need of data transposition
  !!
  !############################################################################
  subroutine channel_cfr (ux, constant)

    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux
    real(mytype), intent(in) :: constant

    integer :: code, i, j, k, jloc
    real(mytype) :: can, ub, uball, coeff

    ub = zero
    uball = zero
    coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

    do k = 1, xsize(3)
       do jloc = 1, xsize(2)
          j = jloc + xstart(2) - 1
          do i = 1, xsize(1)
            ub = ub + ux(i,jloc,k) / ppy(j)
          enddo
       enddo
    enddo

    ub = ub * coeff

    call MPI_ALLREDUCE(ub,uball,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)

    can = - (constant - uball)

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) 'UT', uball, can

    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          ux(i,j,k) = ux(i,j,k) - can
        enddo
      enddo
    enddo

  end subroutine channel_cfr
  !############################################################################
  !############################################################################
  subroutine postprocess_channel(ux1,uy1,uz1,pp3,phi1,ep1)

    use var, ONLY : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_channel
  subroutine visu_channel_init(visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_channel_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_channel
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs channel-specific visualization
  !!
  !############################################################################
  subroutine visu_channel(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, only : Rx1,Ry1,Rz1
    use var, ONLY : nzmsize
    use visu, only : write_field
    
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    integer, intent(in) :: num

    integer :: i,j,k

    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

    !Q=-0.5*(ta1**2+te1**2+di1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
    call write_field(di1, ".", "critq", num, flush = .true.) ! Reusing temporary array, force flush
    Rx1 = zero; Ry1 = zero; Rz1 = zero
    do k=1,xsize(3)
      do j=1,xsize(2)
        do i=1,xsize(1)
          call critR(ta1(i,j,k),td1(i,j,k),tg1(i,j,k),     &
                     tb1(i,j,k),te1(i,j,k),th1(i,j,k),     &
                     tc1(i,j,k),tf1(i,j,k),ti1(i,j,k),     &
                     Rx1(i,j,k),Ry1(i,j,k),Rz1(i,j,k))
         enddo
       enddo
     enddo
     call write_field(Rx1, ".", "Rx", num)
     call write_field(Ry1, ".", "Ry", num)
     call write_field(Rz1, ".", "Rz", num)

  end subroutine visu_channel

  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: critR
  !!      AUTHOR: Yiqian Wang
  !! DESCRIPTION: Calculate the Liutex vector
  !!
  !############################################################################
  subroutine critR(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,Rx,Ry,Rz)
	  
	  implicit none
	  
	  real(mytype),intent(in) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	  real(mytype),intent(out) :: Rx,Ry,Rz
	  
	  real(mytype) :: a(3,3)
	  
	  real(mytype) :: aa, bb, cc
	  real(mytype) :: delta
	  real(mytype) :: tt(3,3)
	  
	  complex(mytype) :: eig1c, eig2c
	  real(mytype) :: eig3r
	  
	  real(mytype) :: qq, rr
	  real(mytype) :: aaaa, bbbb
	  
	  real(mytype) :: vr(3)
	  real(mytype) :: temp,Rmag
	  
	  real(mytype) :: delta1, delta2, delta3
	  
      a(1,1) = dudx
      a(1,2) = dudy
      a(1,3) = dudz
      a(2,1) = dvdx
      a(2,2) = dvdy
      a(2,3) = dvdz
      a(3,1) = dwdx
      a(3,2) = dwdy
      a(3,3) = dwdz
	  
      !-----------------------------------------------------------------------
      ! Cubic Formula
      ! Reference: Numerical Recipes in FORTRAN 77, Second Edition
      ! 5.6 Quadratic and Cubic Equations
      ! Page 179
      !-----------------------------------------------------------------------

      ! cubic equation
      ! x**3 + aa * x**2 + bb * x + cc = 0

      ! coefficients of characteristic equation 
      aa = -(a(1,1)+a(2,2)+a(3,3))

      tt = matmul(a,a)

      bb = -0.5*(tt(1,1)+tt(2,2)+tt(3,3)-(a(1,1)+a(2,2)+a(3,3))**2)

      cc = -(a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))                            &
             -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))                           &
             +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1)))

      ! discriminant of characteristic equation
      delta = 18*aa*bb*cc-4*aa**3*cc+aa**2*bb**2-4*bb**3-27*cc**2

      qq = (aa**2-3*bb)/9.0
      rr = (2*aa**3-9*aa*bb+27*cc)/54.0

      ! delta = rr**2 - qq**3
      ! alleviate round error
      delta = -delta/108

      if(delta > 0.0) then ! one real root and two complex conjugate roots

        aaaa = -sign(1.0_mytype, rr)*(abs(rr)+sqrt(delta))**(1.0_mytype/3.0_mytype)

        if(aaaa == 0.0) then
          bbbb = 0.0
        else
          bbbb = qq/aaaa
        end if

        eig1c = cmplx(-0.5*(aaaa+bbbb)-aa/3.0, 0.5*sqrt(3.0)*(aaaa-bbbb))
        eig2c = cmplx(real(eig1c), -aimag(eig1c))
        eig3r = aaaa+bbbb-aa/3.0

        ! real right eigenvector

        delta1 = (a(1,1)-eig3r)*(a(2,2)-eig3r) - a(2,1)*a(1,2)
        delta2 = (a(2,2)-eig3r)*(a(3,3)-eig3r) - a(2,3)*a(3,2)
        delta3 = (a(1,1)-eig3r)*(a(3,3)-eig3r) - a(1,3)*a(3,1)

        if(delta1 == 0.0 .and. delta2 == 0.0 .and. delta3 == 0.0) then
          write(*,*) 'ERROR: delta1 = delta2 = delta3 = 0.0'
          write(*,*) a(1,1)-eig3r,  a(1,2),       a(1,3)
          write(*,*) a(2,1),        a(2,2)-eig3r, a(2,3)
          write(*,*) a(3,1),        a(3,2),       a(3,3)-eig3r
          stop
        end if

        if(abs(delta1) >= abs(delta2) .and.                                  &
           abs(delta1) >= abs(delta3)) then

          vr(1) = (-(a(2,2)-eig3r)*a(1,3) +         a(1,2)*a(2,3))/delta1
          vr(2) = (         a(2,1)*a(1,3) - (a(1,1)-eig3r)*a(2,3))/delta1
          vr(3) = 1.0

        else if(abs(delta2) >= abs(delta1) .and.                             &
                abs(delta2) >= abs(delta3)) then

          vr(1) = 1.0
          vr(2) = (-(a(3,3)-eig3r)*a(2,1) +         a(2,3)*a(3,1))/delta2
          vr(3) = (         a(3,2)*a(2,1) - (a(2,2)-eig3r)*a(3,1))/delta2

        else if(abs(delta3) >= abs(delta1) .and.                             &
                abs(delta3) >= abs(delta2)) then

           vr(1) = (-(a(3,3)-eig3r)*a(1,2) +         a(1,3)*a(3,2))/delta3
           vr(2) = 1.0
           vr(3) = (         a(3,1)*a(1,2) - (a(1,1)-eig3r)*a(3,2))/delta3

        else

          write(*,*) 'ERROR: '
          write(*,*) delta1, delta2, delta3
          stop

        end if

        temp = sqrt(vr(1)**2+vr(2)**2+vr(3)**2)

        vr(1) = vr(1)/temp
        vr(2) = vr(2)/temp
        vr(3) = vr(3)/temp
		
		! Explicit formula by Yiqian Wang Journal of Hydrodynamics 2019
		temp = (dwdy-dvdz)*vr(1)+(dudz-dwdx)*vr(2)+(dvdx-dudy)*vr(3)
		vr = sign(1.0_mytype,temp)*vr
		temp = abs(temp)
!     When Liutex is basically the same as the vorticity in the direction of 
!     the real eigenvector, the term inside sqrt() could possibly be a small
!     negative. 
		Rmag = temp - sqrt(abs(temp**2-4.0_mytype*(aimag(eig1c))**2))
		Rx = Rmag * vr(1)
		Ry = Rmag * vr(2)
		Rz = Rmag * vr(3)
		
!		if(temp**2-4*(aimag(eig1c)**2)<0.0) then
!		write(*,*)Rx,Rmag,vr(1),vr(2),vr(3),temp**2-4*(aimag(eig1c)**2)
!		endif
		
	else !three real eigenvalues
		Rx = zero
		Ry = zero
		Rz = zero
	endif
	  

	  
  end subroutine critR
  !############################################################################
  !############################################################################
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Applies rotation for t < spinup_time.
  !!
  !############################################################################
  subroutine momentum_forcing_channel(dux1, duy1, duz1, ux1, uy1, uz1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    if (cpg) then
        !! fcpg: add constant pressure gradient in streamwise direction
        if (idir_stream == 1) then
           dux1(:,:,:,1) = dux1(:,:,:,1) + fcpg !* (re/re_cent)**2
        else
           duz1(:,:,:,1) = duz1(:,:,:,1) + fcpg !* (re/re_cent)**2
        endif
    endif

    ! To update to take into account possible flow in z dir
    if (itime < spinup_time .and. iin <= 2) then
       if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
          write(*,*) 'Rotating turbulent channel at speed ',wrotation
       dux1(:,:,:,1) = dux1(:,:,:,1) - wrotation*uy1(:,:,:)
       duy1(:,:,:,1) = duy1(:,:,:,1) + wrotation*ux1(:,:,:)
    endif

  end subroutine momentum_forcing_channel
  !############################################################################
  !############################################################################
  subroutine geomcomplex_channel(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,yp,remp)

    use decomp_2d, only : mytype
    use param, only : zero, one, two, ten
    use ibm

    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: remp
    integer                    :: j
    real(mytype)               :: ym
    real(mytype)               :: zeromach
    real(mytype)               :: h

    epsi(:,:,:) = zero
    h = (yly - two) / two

    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = ten*zeromach

    do j=nyi,nyf
       ym=yp(j)
       if ((ym.le.h).or.(ym.ge.(h+two))) then
          epsi(:,j,:)=remp
       endif
    enddo

    return
  end subroutine geomcomplex_channel
  !############################################################################
  !############################################################################
  subroutine sem_init_channel(ux1, uy1, uz1)

    implicit none

    ! Arguments
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    ! Local variables
    integer :: i, j, k, ii, jj
    real(mytype) :: x, y, z, ftent, um
    integer ( kind = 4 ), parameter :: nsemini = 1000 ! For the moment we fix it but after this can go in the input file
    real(mytype), dimension(3,nsemini) :: eddy, posvor
    real(mytype)     :: volsemini, rrand, ddx, ddy, ddz, lsem, upr, vpr, wpr
    real(mytype), dimension(3) :: dim_min, dim_max

     dim_min(1) = zero
     dim_min(2) = zero
     dim_min(3) = zero
     dim_max(1) = xlx
     dim_max(2) = yly
     dim_max(3) = zlz
     volsemini = xlx * yly * zlz
     ! 3 int to get different random numbers
     do jj = 1, nsemini
        ! Vortex Position
        do ii = 1, 3
           call random_number(rrand)
           posvor(ii,jj) = dim_min(ii)+(dim_max(ii)-dim_min(ii))*rrand
        enddo
        ! Eddy intensity
        do ii = 1, 3
           call random_number(rrand)
           if (rrand <= zpfive) then
              eddy(ii,jj) = -one
           else
              eddy(ii,jj) = +one
           endif 
        enddo
     enddo
     ! Loops to apply the fluctuations 
     do k = 1, xsize(3)
        z = real((k+xstart(3)-1-1),mytype)*dz
        do j = 1, xsize(2)
           if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
           if (istret/=0) y=yp(j+xstart(2)-1)
           do i = 1, xsize(1)
              x = real(i-1,mytype)*dx
              lsem = 0.15_mytype ! For the moment we keep it constant
              upr = zero
              vpr = zero
              wpr = zero
              do jj = 1, nsemini
                 ddx = abs(x-posvor(1,jj))
                 ddy = abs(y-posvor(2,jj))
                 ddz = abs(z-posvor(3,jj))
                 if (ddx < lsem .and. ddy < lsem .and. ddz < lsem) then
                    ! coefficients for the intensity of the fluctuation
                    ftent = (one-ddx/lsem)*(one-ddy/lsem)*(one-ddz/lsem)
                    ftent = ftent / (sqrt(two/three*lsem))**3
                    upr = upr + eddy(1,jj) * ftent
                    vpr = vpr + eddy(2,jj) * ftent
                    wpr = wpr + eddy(3,jj) * ftent
                 endif
              enddo
              upr = upr * sqrt(volsemini/nsemini)
              vpr = vpr * sqrt(volsemini/nsemini)
              wpr = wpr * sqrt(volsemini/nsemini)
              ! 
              um  = one-(y-yly*half)**2 ! we can use a better arroximation 
              if (idir_stream == 1) then
                 ux1(i,j,k)=upr*sqrt(two/three*init_noise*um) + um
                 uy1(i,j,k)=vpr*sqrt(two/three*init_noise*um)
                 uz1(i,j,k)=wpr*sqrt(two/three*init_noise*um)
              else
                 uz1(i,j,k)=upr*sqrt(two/three*init_noise*um) + um
                 uy1(i,j,k)=vpr*sqrt(two/three*init_noise*um)
                 ux1(i,j,k)=wpr*sqrt(two/three*init_noise*um)
              endif
           enddo
        enddo
     enddo

  end subroutine sem_init_channel
  !############################################################################
end module channel
