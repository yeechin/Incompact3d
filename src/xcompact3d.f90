!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use var
  use case

  use transeq, only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier, only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
       calc_divu_constraint, solve_poisson, cor_vel
  use tools, only : restart, simu_stats, apply_spatial_filter, read_inflow
  use turbine, only : compute_turbines
  use ibm_param
  use ibm, only : body
  use genepsi, only : genepsi3d

  implicit none


  call init_xcompact3d()

  do itime=ifirst,ilast
     !t=itime*dt
     t=t0 + (itime0 + itime + 1 - ifirst)*dt
     call simu_stats(2)

     if (iturbine.ne.0) call compute_turbines(ux1, uy1, uz1)

     if (iin.eq.3.and.mod(itime,ntimesteps)==1) then
        call read_inflow(ux_inflow,uy_inflow,uz_inflow,itime/ntimesteps)
     endif

     if ((itype.eq.itype_abl.or.iturbine.ne.0).and.(ifilter.ne.0).and.(ilesmod.ne.0)) then
        call filter(C_filter)
        call apply_spatial_filter(ux1,uy1,uz1,phi1)
     endif

     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

        if (imove.eq.1) then ! update epsi for moving objects
          if ((iibm.eq.2).or.(iibm.eq.3)) then
             call genepsi3d(ep1)
          else if (iibm.eq.1) then
             call body(ux1,uy1,uz1,ep1)
          endif
        endif
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)
#ifdef DEBG
        call check_transients()
#endif
        
        if (ilmn) then
           !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
           call velocity_to_momentum(rho1,ux1,uy1,uz1)
        endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        if (ilmn) then
           call momentum_to_velocity(rho1,ux1,uy1,uz1)
           !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
           !! Note - all other solvers work on velocity always
        endif
        
        call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

     enddo !! End sub timesteps

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

     call simu_stats(3)
	 if (itime.ge.initstat) then ! Statistics of Liutex vector
     if (itype.eq.itype_channel.or.itype.eq.itype_channel_riblet) then
       call update_liutex_channel(ux1,uy1,uz1)
     endif
     endif

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1,Rx1,Ry1,Rz1,vorx1,vory1,vorz1)

  enddo !! End time loop

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_init
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case
  use sandbox, only : init_sandbox
  use forces
  use forces_riblet

  use var

  use navier, only : calc_divu_constraint
  use tools, only : test_speed_min_max, test_scalar_min_max, &
       restart, &
       simu_stats, compute_cfldiff, &
       init_inflow_outflow

  use param, only : ilesmod, jles,itype
  use param, only : irestart

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : nstat, nvisu, nprobe, ilist

  use les, only: init_explicit_les
  use turbine, only: init_turbines

  use visu, only : visu_init, visu_ready

  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body

  use probes, only : init_probes

  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
    
  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'Xcompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'Xcompact3d is run with the provided file -->', trim(InputFN)
  endif

#ifdef ADIOS2
  if (nrank .eq. 0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_io_init()
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  if (itype.eq.itype_cyl.and.iforces.eq.1) then
     call init_forces()
     if (irestart==1) then
        call restart_forces(0)
     endif
  endif
  
  if (itype.eq.itype_channel_riblet.and.iforces_riblet.eq.1) then
     call init_forces_riblet()
  endif

  !####################################################################
  ! initialise visu
  if (ivisu.ne.0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if
  ! compute diffusion number of simulation
  call compute_cfldiff()
  !####################################################################
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     if (itype == itype_sandbox) then
        call init_sandbox(ux1,uy1,uz1,ep1,phi1,1)
     end if
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif

  if ((ioutflow.eq.1).or.(iin.eq.3)) then
     call init_inflow_outflow()
  end if

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

  call init_probes()

  if (iturbine.ne.0) call init_turbines(ux1, uy1, uz1)

  if (itype==2) then
     if(nrank.eq.0)then
        open(42,file='time_evol.dat',form='formatted')
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        open(38,file='forces.dat',form='formatted')
     endif
  endif
  if(itype==itype_channel_riblet) then
	  if(nrank.eq.0) then
		  open(138,file='forces.dat',form='formatted')
	  end if
  end if

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise

  use tools, only : simu_stats
  use param, only : itype, jles, ilesmod
  use probes, only : finalize_probes
  use visu, only : visu_finalise
  use les, only: finalise_explicit_les

  implicit none

  integer :: ierr
  
  if (itype==2) then
     if(nrank.eq.0)then
        close(42)
     endif
  endif
  if (itype==5) then
     if(nrank.eq.0)then
        close(38)
     endif
  endif
  
  call simu_stats(4)
  call finalize_probes()
  call visu_finalise()
  if (ilesmod.ne.0) then
     if (jles.gt.0) call finalise_explicit_les()
  endif
  call decomp_2d_io_finalise()
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d

subroutine check_transients()

  use decomp_2d, only : mytype
  use mpi
  use var
  
  implicit none

  real(mytype) :: dep, dep1
  integer :: code
   
  dep=maxval(abs(dux1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX dux1 ', dep1
 
  dep=maxval(abs(duy1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duy1 ', dep1
 
  dep=maxval(abs(duz1))
  call MPI_ALLREDUCE(dep,dep1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
  if (nrank == 0) write(*,*)'## MAX duz1 ', dep1
  
end subroutine check_transients

subroutine update_liutex_channel(ux1, uy1, uz1)

  use decomp_2d
  use variables
  use channel, only : critR
  use param
  use var, only : ux2, uy2, uz2, ux3, uy3, uz3
  USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
  USE var, only : Rx1,Ry1,Rz1
  USE var, only: vorx1, vory1, vorz1
  use var, ONLY : nzmsize
!  use visu, only : write_field

  use ibm_param, only : ubcx,ubcy,ubcz


  implicit none

  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

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

Rx1 = zero; Ry1 = zero; Rz1 = zero
  do k=1,xsize(3)
     do j=1,xsize(2)
        do i=1,xsize(1)
		  call critR(ta1(i,j,k),td1(i,j,k),tg1(i,j,k),          &
		             tb1(i,j,k),te1(i,j,k),th1(i,j,k),          &
					 tc1(i,j,k),tf1(i,j,k),ti1(i,j,k),          &
					 Rx1(i,j,k),Ry1(i,j,k),Rz1(i,j,k))
        enddo
     enddo
  enddo

  vorx1 = tf1 - th1 ! modified by Yiqian Wang
  vory1 = tg1 - tc1
  vorz1 = tb1 - td1

! output Liutex field to check correctness of the implementation
! call write_field(Rx1, ".", "Rx", 123456)
! call write_field(Ry1, ".", "Ry", 123456)
! call write_field(Rz1, ".", "Rz", 123456)
end subroutine update_liutex_channel
