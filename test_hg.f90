Program test_hg
!$$$  program documentation block
!         .           .            .
!  program name: test_enkf
!    programmer: da,cheng        org: umd      date: 2016-Mar-23
!
!  purpose:
!
!  revision history:
!    2016-Mar-23     da    - creator
!
!  file dependencies:
!
!  attributes: 
!    language: fortran 90
!    machine : 
!
!
!$$$ end documentation block
  Use mod_type,           Only : rdef
  Use mod_lorenz63_fwd,   Only : lorenz63_rk4, t_lorenz63
  Use mod_rnorm,          Only : rnorm, set_random_seed
  Use mod_lorenz63_3dvar, Only : lorenz63_3dvar
  Use mod_lorenz63_letkf, Only : lorenz63_letkf
  Use mod_lorenz63_enkf,  Only : lorenz63_enkf
  Use mod_inflation,      Only : inflate_rtps, inflate_rtpp
  Implicit none

  Integer,parameter :: nx = 3
  Integer,parameter :: filter_LETKF  = 1
  Integer,parameter :: filter_POENKF = 2

  Integer                :: filter_opt
  Logical                :: loutm = .false.  ! output ensemble members
  Integer                :: nn 
  Integer                :: initopt 
  Integer                :: ncycle, nspin
  Type(t_lorenz63)       :: xt
  Real(rdef),allocatable :: xb(:,:), xa(:,:)
  Real(rdef),allocatable :: dx(:,:)
  Real(rdef)             :: Pb(nx,nx), Pa(nx,nx)
  Real(rdef)             :: xbm(nx), xam(nx), xam_enkf(nx), xam_3dvar(nx)
  Real(rdef),allocatable :: yo(:,:)
  Real(rdef)             :: erro(nx)
  Real(rdef)             :: mseb(nx,nx), msea(nx,nx)
  Real(rdef)             :: infl_param
  Real(rdef)             :: rdef1d(nx)
  Real(rdef)             :: errb(3,3)
  Integer                :: i, j, k


!----------------------------------------------------------------------------------
! 0. configuration for lorenz63 model
!----------------------------------------------------------------------------------
  xt%dt   = 0.01d0   ! length for one time step
  xt%nt   = 8       ! generate nt-step forecast
  xt%para = (/ 10.0d0, 28.0d0, 8.0d0/3.0d0 /)  

  nspin   = 4000     ! spin up the model with nspin steps, Then the final state 
                    ! is used as the initial condition
  ncycle  = 4000     ! total n-step forecast cycles
  Write(6,*) "input obs error:"
  Read(*,*) erro(1)
  erro(:) = erro(1)

!----------------------------------------------------------------------------------
! 1. spin up the lorenz 63 model
!----------------------------------------------------------------------------------
  xt%x0   = (/  0.1d0,  0.1d0, 0.1d0 /)
  Write(6,*) "spin up lorenz63 model: dt =", xt%dt, ", nspin =", nspin
  Call lorenz63_rk4( xt%x0(:), xt%xn(:), xt%para, nspin, xt%dt )
  xt%x0(:) = xt%xn(:)
  Write(6,*) "initial condition: x=", xt%x0(:)

!----------------------------------------------------------------------------------
! 2. generate simulated observations
!----------------------------------------------------------------------------------
  Write(6,*) "generate simulated obs: err(yo) = ", erro(:)
  Allocate( yo(nx,ncycle) )
  rdef1d = xt%x0(:)
  !Call set_random_seed( 0 )
  Write(6,*) "Write simulated obs (yo) into fort.10010"
  Do i = 1, ncycle
     Call lorenz63_rk4( rdef1d, rdef1d, xt%para, xt%nt, xt%dt )
     Do j = 1, nx
        yo(j,i) = rdef1d(j) + erro(j)*rnorm()
     Enddo
     Write(10010,*) ( yo(j,i), j = 1, nx )
  Enddo

!----------------------------------------------------------------------------------
! 3. generate initial perurbation
!----------------------------------------------------------------------------------
  Write(6,*) "input ensemble number nn"
  Read(*,*) nn
  Write(6,*) "nmember=", nn
  Allocate( xb(nx,nn), xa(nx,nn) )
  Allocate( dx(nx,nn) )
  ! generate initial ensembles based on historical trajectories
  Do k = 1, nn
     xb(:,k) = (/0.1d0, 0.1d0, 0.1d0/) 
     !Call lorenz63_rk4( xb(:,k), xb(:,k), xt%para, nspin+10*k*xt%nt, xt%dt )
     Call lorenz63_rk4( xb(:,k), xb(:,k), xt%para, 2*nspin+5*k*xt%nt, xt%dt )
  Enddo

!----------------------------------------------------------------------------------
! 3. load B
!----------------------------------------------------------------------------------
  Write(6,*) "load background error covariance matrix B(3x3) from fort.1040"
  Do i = 1, nx
     Read(1040,*) ( errb(i,j), j = 1, nx )
     Write(6,*)   ( errb(i,j), j = 1, nx )
  Enddo
  Write(6,*) "rawB="
  !errb = 1000000000*errb
  !errb(1,2:3) = 0.d0
  !errb(2,1) = 0.d0; errb(2,3) = 0.d0;
  !errb(3,1:2) = 0.d0
  Do i = 1, nx
     Write(6,*)   ( errb(i,j), j = 1, nx )
  Enddo

  pause "continue"

!----------------------------------------------------------------------------------
! 4. run EnKF at every cycle
!----------------------------------------------------------------------------------
  Write(6,*) "input inflation_param="
  Read(*,*) infl_param
  Write(6,*) "infl_param=", infl_param

  !Write(6,*) "Select filter "
  !Write(6,*) "LETKF=",filter_LETKF
  !Write(6,*) "POENKF=",filter_POENKF
  !Read(*,*) filter_opt
  !If ( filter_opt == filter_LETKF ) then
  !   Write(6,*) "using LETKF"
  !Elseif ( filter_opt == filter_POENKF ) then
  !   Write(6,*) "using PO-ENKF"
  !Endif

  !Write(6,*) "run EnSRF at every cycle"
  Write(6,*) "Write truth (xt)     into fort.10020"
  Write(6,*) "Write background(xb) into fort.10030"
  Write(6,*) "Write analysis (xa)  into fort.10040"

  Do i = 1, ncycle

     Write(6,*) "----------------------------------"
     Write(6,*) "i=", i

     ! true trajectory
     Call lorenz63_rk4( xt%x0(:), xt%xn(:), xt%para, xt%nt, xt%dt )

     ! ensemble forecast
     Do k = 1, nn
        Call lorenz63_rk4( xb(:,k), xb(:,k), xt%para, xt%nt, xt%dt )
        if ( loutm ) Then
           Write(2000+k,*) (xb(j,k),j=1,nx)
        Endif
     Enddo

     ! calculate ensemble mean & Pb
     Do j = 1, nx
        xbm(j) = SUM( xb(j,:) )/nn
        dx(j,:) = xb(j,:) - xbm(j)
     Enddo
     Pb = MATMUL( dx, TRANSPOSE(dx) )/(nn-1)

     ! Calculate local bred vector dimension
 
     ! assimilating obs
     If ( .true. ) Then  ! LETKF
        Write(6,*) "using LETKF"
        Call lorenz63_letkf( nn, &
                          xb, &
                          (/.true., .true., .true./), &
                          yo(:,i), &
                          erro(:), &
                          xa, &
                          xam )
        !Call lorenz63_enkf( nn, &
        !                  xb, &
        !                  (/.true., .true., .true./), &
        !                  yo(:,i), &
        !                  erro(:), &
        !                  xa, &
        !                  xam )

        Call inflate_rtpp( nx, nn, xb, xa, infl_param )

        xam_enkf = xam

        ! remove mean
        DO j = 1, nx
           do k = 1, nn
              xa(j,k) = xa(j,k) - xam(j)
           enddo
        Enddo

        ! 3d-var
        Call lorenz63_3dvar( xam_enkf, &
                             errb, &
                             (/.true., .true., .true./), &
                             yo(:,i), &
                             erro(:), &
                             xam_3dvar )

        xam(:) = 0.5*xam_enkf(:) + 0.5*xam_3dvar(:)
        ! add back
        DO j = 1, nx
           do k = 1, nn
              xa(j,k) = xa(j,k) + xam(j)
           enddo
        Enddo
 
     Endif

     If ( loutm ) Then
        Do k = 1, nn
           Write(4000+k,*) (xa(j,k),j=1,nx)
        Enddo
     Endif

     ! calculate ensemble Pa
     Do j = 1, nx
        dx(j,:) = xa(j,:) - xam(j)
     Enddo
     Pa = MATMUL( dx, TRANSPOSE(dx) )/(nn-1)


     Write(10020,*) ( xt%xn(j), j = 1, nx )
     Write(10030,*) ( xbm(j), j = 1, nx )
     Write(10040,*) ( xam(j), j = 1, nx )
     Write(10050,*) (Pb(j,j), j = 1, nx)
     Write(10060,*) (Pa(j,j), j = 1, nx)

     Write(6,*) "Tr(Pa)/Tr(Pb)=", ( Pa(1,1)+Pa(2,2)+Pa(3,3) )/( Pb(1,1)+Pb(2,2)+Pb(3,3) )

    ! Relaxation to Prior Spread
    !Call inflate_rtps( nx, nn, xb, xa, infl_param )
    !Call inflate_rtpp( nx, nn, xb, xa, infl_param )

    ! Xa members after relaxation
     if ( loutm ) Then
        Do k = 1, nn
           Write(6000+k,*) (xa(j,k),j=1,nx)
        Enddo
     Endif

    ! mean xa & Pa after relaxation
     !Do j = 1, nx
     !   xam(j) = SUM( xa(j,:) )/nn
     !   dx(j,:) = xa(j,:) - xam(j)
     !Enddo
     !Pa = MATMUL( dx, TRANSPOSE(dx) )/(nn-1)
     !Write(10065,*) ( Pa(j,j), j = 1, nx ) 

     xt%x0(:) = xt%xn(:)
     xb(:,:)   = xa(:,:)


  Enddo

  Deallocate( yo ) 
  Deallocate( xb, xa )

Endprogram

