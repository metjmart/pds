!=================================================================
!=================================================================
!=================================================================
!      =====                                           =====
!      =====              Module netNudge              =====
!      =====                                           =====
!=================================================================
!=================================================================
!=================================================================

!=======================================================================
! This module reads a list of WRF input files, then creates WRF
! grid/analysis nudging files in NetCDF format.
!
! Revision history:
!    Date          Programmer                  Description of change
!    ----          ----------                  ---------------------
!    July 2012     Craig  Mattocks, UM-RSMAS   Wrote original code
!=======================================================================
! To build executable:
!
!    setenv NETCDF /usr/local/netcdf4
!    setenv NETCDF_INC "-I${NETCDF}/include"
!    setenv NETCDF_LIB "-L${NETCDF}/lib -lnetcdf -lnetcdff"
!
!    Compiler flags (Intel Fortran compiler)
!       setenv F90 ifort
!       setenv FFLAGS "-O0 -free -g -check all -warn all -traceback -fp-stack-check -fpe0"
!       setenv FFLAGS "-O3 -free -xHost -ip -fp-model fast=2 -no-heap-arrays -no-prec-div -no-prec-sqrt -ftz -align all -m64"
!
!    Compiler flags (Portland Group Fortran compiler)
!       setenv F90 pgf90
!       setenv FFLAGS "-O3 -Mfree -fastsse -Mvect=noaltcode -Msmartalloc -Mprefetch=distance:8 -Mfprelaxed"
!
!    ${F90} -o netNudge ${FFLAGS} ${NETCDF_INC} ${NETCDF_LIB} netNudge.f90
!    testNetNudge wrfinput_d01_t1 wrfinput_d01_t2 etc.
!
! Check output:
!     ncdump -v U_NDG_NEW -f fortran wrffdda_d01
!=======================================================================
module netNudge

   !--------------------------------------------
   ! Force explicit declaration of all variables
   !--------------------------------------------
   implicit none

   !--------------------------------------------------------------
   ! Retain values of class instance variables between invocations
   !--------------------------------------------------------------
   save

   !-----------
   ! Visibility
   !-----------
   private                            ! Everything is private by default
   public :: RUN_UNIT_TESTS
   public :: deBug, xBug, yBug, zBug, tBug

   !----------------
   ! Nudging control
   !----------------
   logical, parameter :: rampDown = .true.   ! Ramp nudging down to
                                             ! previous values at end of
                                             ! simulation?
   integer, parameter :: rampTime = 1        ! Nudging time slot to ramp
                                             ! down to at end of simulation

   !----------------------------------
   ! Global class variables for NetCDF
   !----------------------------------
   character (len = *), parameter :: title       = "OUTPUT FROM REAL_EM V3.4 PREPROCESSOR"
   character (len = *), parameter :: institution = "University of Miami - RSMAS"
   character (len = *), parameter :: contact     = "dnolan@rsmas.miami.edu"
   character (len = *), parameter :: source      = "netNudge module"
   character (len = *), parameter :: history     = "Generated for WRF ideal em_quarter_ss simulation"
   character (len = *), parameter :: references  = "http://www.mmm.ucar.edu/wrf/OnLineTutorial/Compile/arw_compile3_ideal.htm"
   character (len = *), parameter :: model       = "WRF-ARW"
   character (len = *), parameter :: version     = "3.4"
   character (len = *), parameter :: Conventions = "CF-1.5"
   character (len = *), parameter :: memoryOrder = "XYZ"
   integer            , parameter :: fieldType   = 104

   integer            , save :: diff_opt, km_opt, damp_opt
   integer            , save :: mp_physics, cu_physics
   integer            , save :: ra_lw_physics, ra_sw_physics
   integer            , save :: sf_sfclay_physics, sf_surface_physics, bl_pbl_physics
   integer            , save :: surface_input_source, sst_update
   integer            , save :: grid_fdda ,  gfdda_interval_m,  gfdda_end_h
   integer            , save :: grid_sfdda, sgfdda_interval_m, sgfdda_end_h
   integer            , save :: hypsometric_opt
   integer            , save :: grid_id, parent_id
   integer            , save :: i_parent_start, j_parent_start
   integer            , save :: parent_grid_ratio
   integer            , save :: num_land_cat, isWater, isLake
   integer            , save :: isIce, isUrban, iSoilWater
   real               , save :: dx, dy, dt
   real               , save :: dampcoef, khdif, kvdif
   real               , save :: cen_lat, cen_lon
   real               , save :: truelat1, truelat2
   real               , save :: moad_cen_lat, stand_lon
   real               , save :: pole_lat, pole_lon
   real               , save :: map_proj
   character (len = 4), save :: mminlu

   !-----------
   ! Debug info
   !-----------
   logical, parameter :: deBug = .true.   ! Print debug info?
   integer, parameter :: xBug = 10
   integer, parameter :: yBug = 10
   integer, parameter :: zBug =  1
   integer, parameter :: tBug =  1

   !----------------
   ! Class functions
   !----------------
   contains

      !=================================================================
      ! This subroutine reads in a list of files provided by a user in
      ! the command line.
      !
      ! Input:
      !    none
      !
      ! Output:
      !    nf              Number of files
      !    files           List of file names
      !    nudge_int_min   user-defined nudging time interval (minutes)
      !=================================================================
      subroutine userInput (nf, files, nudge_int_min)
         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         implicit none

         !---------------
         ! File variables
         !---------------
         !Note : command arguments must be passed as "file1 file2 file3 nudge_interval_minutes" for example
         integer, intent(out) :: nf, nudge_int_min
         character (len = *), dimension(:), allocatable, intent(out) :: files
         character*8 :: nudge_int_min_char
         integer :: n, nArgs

         !----------------------------------------
         ! Determine how many files user has input
         !----------------------------------------
         nArgs = command_argument_count()
         nf = nArgs - 1         !Note, currently number of files is (number of args - 1)...
                                !if we add more arguments to the command line, this will have to change.

         !------------------------------------------------------------
         ! Read list of files and define nudge time interval (minutes)
         !------------------------------------------------------------
         if (nArgs > 0) then
            if (allocated(files)) deallocate(files)
            allocate(files(nf))
            do n = 1, nf
               call get_command_argument(n, files(n))
            end do
            call get_command_argument(nArgs, nudge_int_min_char)
            read(nudge_int_min_char,'(i10)') nudge_int_min
         else
            write (*,'(a)') "Error: no WRF input files provided!"
            stop
         end if
      end subroutine userInput

      !=================================================================
      ! Input:
      !    n                     Ordinal number of file to read in list
      !                             of WRF input files
      !    file                  Name of WRF input file
      !
      ! Output:
      !    west_east             x dimension
      !    south_north           y dimension
      !    bottom_top            z dimension
      !    Time                  t dimension
      !    DateStrLen            Number of characters in a date string
      !    start_date            Start date of 1st WRF input file
      !    simulation_start_date Simulation start date of 1st WRF input
      !                             file
      !    Times                 Array of WRF date character strings
      !    U                     4D array containing x-wind component
      !                             (m/s)
      !    V                     4D array containing y-wind component
      !                             (m/s)
      !    T                     4D array containing perturbation
      !                             potential temperature (theta-t0) (K)
      !    Q                     4D array containing water vapor mixing ratio
      !                             (QVAPOR) (kg H20/kg air)
      !    PH                    4D array containing perturbation
      !                             geopotential (m2/s2)
      !    MU                    3D array containing perturbation dry
      !                             air mass in column (Pa)
      !    other                 Values of global class NetCDF variables
      !                             at top of this module modified when
      !                             1st WRF input file is read (n = 1)
      !=================================================================
      subroutine readWRFinput (n, file,                                 &
                               west_east, south_north, bottom_top, Time,&
                               DateStrLen,                              &
                               start_date, simulation_start_date,       &
                               Times, U, V, T, Q, PH, MU)
         !--------------------------------
         ! Import Fortran 90 NetCDF module
         !--------------------------------
         use netcdf

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         implicit none

         !---------------
         ! NetCDF file ID
         !---------------
         integer :: ncID

         !---------------------
         ! Number of dimensions
         !---------------------
         integer :: nDims

         !--------------------
         ! Number of variables
         !--------------------
         integer :: nVars

         !----------------------------
         ! Number of global attributes
         !----------------------------
         integer :: nGlobalAtts

         !----------------------------------
         ! Number of unlimited dimension IDs
         !----------------------------------
         integer :: nUnlimitedDimIDs

         !---------------
         ! File variables
         !---------------
         integer            , intent(in) :: n
         character (len = *), intent(in) :: file

         !------------------
         ! Dimension lengths
         !------------------
         integer, intent(out) :: west_east
         integer, intent(out) :: south_north
         integer, intent(out) :: bottom_top
         integer, intent(out) :: Time
         integer, intent(out) :: DateStrLen

         !---------------
         ! Date variables
         !---------------
         character (len = *), intent(out) :: start_date
         character (len = *), intent(out) :: simulation_start_date

         !---------------
         ! WRF input data
         !---------------
         character (len = *), allocatable, dimension(:)      , intent(out) :: Times
         real               , allocatable, dimension(:,:,:,:), intent(out) ::  U
         real               , allocatable, dimension(:,:,:,:), intent(out) ::  V
         real               , allocatable, dimension(:,:,:,:), intent(out) ::  T
         real               , allocatable, dimension(:,:,:,:), intent(out) ::  Q
         real               , allocatable, dimension(:,:,:,:), intent(out) :: PH
         real               , allocatable, dimension(:,:  ,:), intent(out) :: MU

         !--------------
         ! Dimension IDs
         !--------------
!        character(len = nf90_max_name) :: TimeDimName
         integer ::        Time_dim
         integer ::  DateStrLen_dim
         integer ::   west_east_dim
         integer :: south_north_dim
         integer ::  bottom_top_dim

         !-------------
         ! Variable IDs
         !-------------
         integer :: Times_id
         integer ::     U_id
         integer ::     V_id
         integer ::     T_id
         integer ::     Q_id
         integer ::    PH_id
         integer ::    MU_id

         !----------------------------------------------
         ! Rank (number of dimensions) for each variable
         !----------------------------------------------
         integer, parameter :: Times_rank = 2
         integer, parameter ::     U_rank = 4
         integer, parameter ::     V_rank = 4
         integer, parameter ::     T_rank = 4
         integer, parameter ::     Q_rank = 4
         integer, parameter ::    PH_rank = 4
         integer, parameter ::    MU_rank = 4

         !----------------
         ! Variable shapes
         !----------------
!        integer, dimension (Times_rank) :: Times_dims
!        integer, dimension (    U_rank) ::     U_dims
!        integer, dimension (    V_rank) ::     V_dims
!        integer, dimension (    T_rank) ::     T_dims
!        integer, dimension (    Q_rank) ::     Q_dims
!        integer, dimension (   PH_rank) ::    PH_dims
!        integer, dimension (   MU_rank) ::    MU_dims

         !--------------------
         ! Open WRF input file
         !--------------------
         if (deBug) write (*,*) "readWRFinput: Opening WRF input file: ", file
         call check(nf90_open(trim(file), nf90_nowrite, ncID))

         !-----------------------------------------------
         ! Read general characteristics of WRF input file
         !-----------------------------------------------
         call check( nf90_inquire(ncID, nDims, nVars, nGlobalAtts, nUnlimitedDimIDs) )
         if (deBug) then
            write (*,*) "=========================="
            write (*,*) "WRF netCDF input file info"
            write (*,*) "=========================="
            write (*,*) "ncID              = ", ncID
            write (*,*) "Dimensions        = ", nDims
            write (*,*) "Variables         = ", nVars
            write (*,*) "Global attributes = ", nGlobalAtts
            write (*,*) "Unlimited Dim IDs = ", nUnlimitedDimIDs
         end if

         !-----------------------------------
         ! Read dimensions of WRF data arrays
         !-----------------------------------
         call check( nf90_inq_dimid(ncID, "Time", Time_dim) )
         call check( nf90_inquire_dimension(ncID, Time_dim, len = Time) )
         if (deBug) write (*, '("readWRFinput: Time dimension     = ", i3)') Time

         call check( nf90_inq_dimid(ncID, "DateStrLen", DateStrLen_dim) )
         call check( nf90_inquire_dimension(ncID, DateStrLen_dim, len = DateStrLen) )
         if (deBug) write (*, '("readWRFinput: Date string length = ", i3)') DateStrLen

         call check( nf90_inq_dimid(ncID, "west_east", west_east_dim) )
         call check( nf90_inquire_dimension(ncID, west_east_dim, len = west_east) )
         if (deBug) write (*, '("readWRFinput: x    dimension     = ", i3)') west_east

         call check( nf90_inq_dimid(ncID, "south_north", south_north_dim) )
         call check( nf90_inquire_dimension(ncID, south_north_dim, len = south_north) )
         if (deBug) write (*, '("readWRFinput: y    dimension     = ", i3)') south_north

         call check( nf90_inq_dimid(ncID, "bottom_top", bottom_top_dim) )
         call check( nf90_inquire_dimension(ncID, bottom_top_dim, len = bottom_top) )
         if (deBug) write (*, '("readWRFinput: z    dimension     = ", i3)') bottom_top

         !------------------------------------
         ! Allocate memory for WRF data arrays
         !------------------------------------
         if (.not. allocated(Times)) allocate(Times(Time))
         if (.not. allocated( U   )) allocate( U(west_east, south_north, bottom_top, Time))
         if (.not. allocated( V   )) allocate( V(west_east, south_north, bottom_top, Time))
         if (.not. allocated( T   )) allocate( T(west_east, south_north, bottom_top, Time))
         if (.not. allocated( Q   )) allocate( Q(west_east, south_north, bottom_top, Time))
         if (.not. allocated(PH   )) allocate(PH(west_east, south_north, bottom_top, Time))
         if (.not. allocated(MU   )) allocate(MU(west_east, south_north            , Time))
         if (deBug) write (*, '("readWRFinput: Memory allocated for WRF input data arrays")')

         !--------------
         ! Read WRF data
         !--------------
         if (deBug) write (*, '("readWRFinput: Reading WRF input data...")')
         if (deBug) write (*, '("readWRFinput: Sample WRF input array values:")')

         call check( nf90_inq_varid(ncID, "Times" , Times_id) )
         call check( nf90_get_var  (ncID, Times_id, Times   ) )
         if (deBug) write (*, '("readWRFinput:    Times = ", a   )') Times(tBug)

         call check( nf90_inq_varid(ncID, "U" , U_id) )
         call check( nf90_get_var  (ncID, U_id, U   ) )
         if (deBug) write (*, '("readWRFinput:    U     = ", f8.3)')  U(xBug, yBug, zBug, tBug)

         call check( nf90_inq_varid(ncID, "V" , V_id) )
         call check( nf90_get_var  (ncID, V_id, V   ) )
         if (deBug) write (*, '("readWRFinput:    V     = ", f8.3)')  V(xBug, yBug, zBug, tBug)

         call check( nf90_inq_varid(ncID, "T" , T_id) )
         call check( nf90_get_var  (ncID, T_id, T   ) )
         if (deBug) write (*, '("readWRFinput:    T     = ", f8.3)')  T(xBug, yBug, zBug, tBug)

         call check( nf90_inq_varid(ncID, "QVAPOR", Q_id) )
         call check( nf90_get_var  (ncID, Q_id, Q   ) )
         if (deBug) write (*, '("readWRFinput:    Q     = ", f8.3)')  Q(xBug, yBug, zBug, tBug)

         call check( nf90_inq_varid(ncID, "PH" , PH_id) )
         call check( nf90_get_var  (ncID, PH_id, PH   ) )
         if (deBug) write (*, '("readWRFinput:    PH    = ", f8.3)') PH(xBug, yBug, zBug, tBug)

         call check( nf90_inq_varid(ncID, "MU" , MU_id) )
         call check( nf90_get_var  (ncID, MU_id, MU   ) )
         if (deBug) write (*, '("readWRFinput:    MU    = ", f8.3)') MU(xBug, yBug      , tBug)

         !---------------------------
         ! Read WRF global attributes
         !---------------------------
         if (n == 1) then
            if (deBug) write (*, '("readWRFinput: Reading WRF global attributes...")')

            call check( nf90_get_att(ncID, nf90_global, "START_DATE", start_date) )
            call check( nf90_get_att(ncID, nf90_global, "SIMULATION_START_DATE", simulation_start_date) )

            call check( nf90_get_att(ncID, nf90_global, "DIFF_OPT", diff_opt) )
            call check( nf90_get_att(ncID, nf90_global, "KM_OPT"  , km_opt  ) )
            call check( nf90_get_att(ncID, nf90_global, "DAMP_OPT", damp_opt) )

            call check( nf90_get_att(ncID, nf90_global, "MP_PHYSICS"   , mp_physics) )
            call check( nf90_get_att(ncID, nf90_global, "CU_PHYSICS"   , cu_physics) )

            call check( nf90_get_att(ncID, nf90_global, "RA_LW_PHYSICS", ra_lw_physics) )
            call check( nf90_get_att(ncID, nf90_global, "RA_SW_PHYSICS", ra_sw_physics) )

            call check( nf90_get_att(ncID, nf90_global, "SF_SFCLAY_PHYSICS" , sf_sfclay_physics ) )
            call check( nf90_get_att(ncID, nf90_global, "SF_SURFACE_PHYSICS", sf_surface_physics) )
            call check( nf90_get_att(ncID, nf90_global, "BL_PBL_PHYSICS"    , bl_pbl_physics    ) )

            call check( nf90_get_att(ncID, nf90_global, "SURFACE_INPUT_SOURCE", surface_input_source) )
            call check( nf90_get_att(ncID, nf90_global, "SST_UPDATE"          , sst_update          ) )

            call check( nf90_get_att(ncID, nf90_global, "GRID_FDDA"        , grid_fdda        ) )
            call check( nf90_get_att(ncID, nf90_global, "GFDDA_INTERVAL_M" , gfdda_interval_m ) )
            call check( nf90_get_att(ncID, nf90_global, "GFDDA_END_H"      , gfdda_end_h      ) )

            call check( nf90_get_att(ncID, nf90_global, "GRID_SFDDA"       , grid_sfdda       ) )
            call check( nf90_get_att(ncID, nf90_global, "SGFDDA_INTERVAL_M", sgfdda_interval_m) )
            call check( nf90_get_att(ncID, nf90_global, "SGFDDA_END_H"     , sgfdda_end_h     ) )

            call check( nf90_get_att(ncID, nf90_global, "HYPSOMETRIC_OPT", hypsometric_opt) )

            call check( nf90_get_att(ncID, nf90_global, "GRID_ID"          , grid_id          ) )
            call check( nf90_get_att(ncID, nf90_global, "PARENT_ID"        , parent_id        ) )
            call check( nf90_get_att(ncID, nf90_global, "I_PARENT_START"   , i_parent_start   ) )
            call check( nf90_get_att(ncID, nf90_global, "J_PARENT_START"   , j_parent_start   ) )
            call check( nf90_get_att(ncID, nf90_global, "PARENT_GRID_RATIO", parent_grid_ratio) )

            call check( nf90_get_att(ncID, nf90_global, "NUM_LAND_CAT", num_land_cat) )
            call check( nf90_get_att(ncID, nf90_global, "ISWATER"     , isWater     ) )
            call check( nf90_get_att(ncID, nf90_global, "ISLAKE"      , isLake      ) )
            call check( nf90_get_att(ncID, nf90_global, "ISICE"       , isIce       ) )
            call check( nf90_get_att(ncID, nf90_global, "ISURBAN"     , isUrban     ) )
            call check( nf90_get_att(ncID, nf90_global, "ISOILWATER"  , iSoilWater  ) )

            call check( nf90_get_att(ncID, nf90_global, "DX", dx) )
            call check( nf90_get_att(ncID, nf90_global, "DY", dy) )
            call check( nf90_get_att(ncID, nf90_global, "DT", dt) )

            call check( nf90_get_att(ncID, nf90_global, "DAMPCOEF", dampcoef) )
            call check( nf90_get_att(ncID, nf90_global, "KHDIF"   , khdif   ) )
            call check( nf90_get_att(ncID, nf90_global, "KVDIF"   , kvdif   ) )

            call check( nf90_get_att(ncID, nf90_global, "CEN_LAT", cen_lat) )
            call check( nf90_get_att(ncID, nf90_global, "CEN_LON", cen_lon) )

            call check( nf90_get_att(ncID, nf90_global, "TRUELAT1", truelat1) )
            call check( nf90_get_att(ncID, nf90_global, "TRUELAT2", truelat2) )

            call check( nf90_get_att(ncID, nf90_global, "MOAD_CEN_LAT", moad_cen_lat) )
            call check( nf90_get_att(ncID, nf90_global, "STAND_LON"   , stand_lon   ) )

            call check( nf90_get_att(ncID, nf90_global, "POLE_LAT", pole_lat) )
            call check( nf90_get_att(ncID, nf90_global, "POLE_LON", pole_lon) )
            call check( nf90_get_att(ncID, nf90_global, "MAP_PROJ", map_proj) )

            call check( nf90_get_att(ncID, nf90_global, "MMINLU"  , mminlu) )
         end if

         if (deBug) write (*, '("readWRFinput: Finished readWRFinput")')
      end subroutine readWRFinput

      !=================================================================
      ! This subroutine writes a WRF nudging file in NetCDF format.
      !
      ! Input:
      !    west_east             x dimension
      !    south_north           y dimension
      !    bottom_top            z dimension
      !    Time                  t dimension
      !    DateStrLen            Number of characters in a date string
      !    start_date            Start date of 1st WRF input file
      !    simulation_start_date Simulation start date of 1st WRF input
      !                             file
      !    Times                 Array of WRF date character strings
      !    *_NDG_OLD             4D arrays of nudging values for OLD
      !                             time slot
      !    *_NDG_NEW             4D arrays of nudging values for NEW
      !                             time slot
      !    other                 Global class NetCDF variables at top of
      !                             this module
      !
      ! Output:
      !    WRF analysis/grid nudging file named wrffdda_d<nn>
      !
      ! Check output:
      !    ncdump -v U_NDG_NEW -f fortran wrffdda_d01
      !=================================================================
      subroutine writeNudge (west_east, south_north, bottom_top, Time,  &
                             DateStrLen,                                &
                             start_date, simulation_start_date,         &
                             Times,                                     &
                              U_NDG_OLD,  U_NDG_NEW,                    &
                              V_NDG_OLD,  V_NDG_NEW,                    &
                              T_NDG_OLD,  T_NDG_NEW,                    &
                              Q_NDG_OLD,  Q_NDG_NEW,                    &
                             PH_NDG_OLD, PH_NDG_NEW,                    &
                             MU_NDG_OLD, MU_NDG_NEW)
         !--------------------------------
         ! Import Fortran 90 NetCDF module
         !--------------------------------
         use netcdf

         !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         implicit none

         !---------------
         ! NetCDF file ID
         !---------------
         integer :: ncID

         !---------------
         ! File variables
         !---------------
         character (len =  *), parameter  :: file_prefix = "wrffdda_d"
!        character (len =  *), parameter  :: file_suffix = ".nc"
         character (len = 11)             :: file_name

         !------------------
         ! Old pre-fill flag
         !------------------
         integer :: oldFillMode   ! <--- Craig Mattocks suggested adding this line for faster netcdf generation (no variable prefilling).  Feb 25, 2014

         !------------------
         ! Dimension lengths
         !------------------
         integer, intent(in) :: west_east
         integer, intent(in) :: south_north
         integer, intent(in) :: bottom_top
         integer, intent(in) :: Time
         integer, parameter  :: one_stag = 1
         integer, intent(in) :: DateStrLen

         !---------------
         ! Date variables
         !---------------
         character (len = *), intent(in) :: start_date
         character (len = *), intent(in) :: simulation_start_date
         character (len = DateStrLen)    :: creation_date

         !-----------------
         ! WRF nudging data
         !-----------------
         character (len = *), dimension(Time), intent(in) :: Times

         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  U_NDG_OLD
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  V_NDG_OLD
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  T_NDG_OLD
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  Q_NDG_OLD
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) :: PH_NDG_OLD
         real, dimension(west_east, south_north, one_stag  , Time), intent(in) :: MU_NDG_OLD

         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  U_NDG_NEW
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  V_NDG_NEW
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  T_NDG_NEW
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) ::  Q_NDG_NEW
         real, dimension(west_east, south_north, bottom_top, Time), intent(in) :: PH_NDG_NEW
         real, dimension(west_east, south_north, one_stag  , Time), intent(in) :: MU_NDG_NEW

         !--------------
         ! Dimension IDs
         !--------------
         integer ::        Time_dim
         integer ::   west_east_dim
         integer :: south_north_dim
         integer ::  bottom_top_dim
         integer ::    one_stag_dim
         integer ::  DateStrLen_dim

         !-------------
         ! Variable IDs
         !-------------
         integer ::      Times_id

         integer ::  U_NDG_OLD_id
         integer ::  V_NDG_OLD_id
         integer ::  T_NDG_OLD_id
         integer ::  Q_NDG_OLD_id
         integer :: PH_NDG_OLD_id
         integer :: MU_NDG_OLD_id

         integer ::  U_NDG_NEW_id
         integer ::  V_NDG_NEW_id
         integer ::  T_NDG_NEW_id
         integer ::  Q_NDG_NEW_id
         integer :: PH_NDG_NEW_id
         integer :: MU_NDG_NEW_id

         !----------------------------------------------
         ! Rank (number of dimensions) for each variable
         !----------------------------------------------
         integer, parameter ::      Times_rank = 2

         integer, parameter ::  U_NDG_OLD_rank = 4
         integer, parameter ::  V_NDG_OLD_rank = 4
         integer, parameter ::  T_NDG_OLD_rank = 4
         integer, parameter ::  Q_NDG_OLD_rank = 4
         integer, parameter :: PH_NDG_OLD_rank = 4
         integer, parameter :: MU_NDG_OLD_rank = 4

         integer, parameter ::  U_NDG_NEW_rank = 4
         integer, parameter ::  V_NDG_NEW_rank = 4
         integer, parameter ::  T_NDG_NEW_rank = 4
         integer, parameter ::  Q_NDG_NEW_rank = 4
         integer, parameter :: PH_NDG_NEW_rank = 4
         integer, parameter :: MU_NDG_NEW_rank = 4

         !----------------
         ! Variable shapes
         !----------------
         integer, dimension (     Times_rank) ::      Times_dims

         integer, dimension ( U_NDG_OLD_rank) ::  U_NDG_OLD_dims
         integer, dimension ( V_NDG_OLD_rank) ::  V_NDG_OLD_dims
         integer, dimension ( T_NDG_OLD_rank) ::  T_NDG_OLD_dims
         integer, dimension ( Q_NDG_OLD_rank) ::  Q_NDG_OLD_dims
         integer, dimension (PH_NDG_OLD_rank) :: PH_NDG_OLD_dims
         integer, dimension (MU_NDG_OLD_rank) :: MU_NDG_OLD_dims

         integer, dimension ( U_NDG_NEW_rank) ::  U_NDG_NEW_dims
         integer, dimension ( V_NDG_NEW_rank) ::  V_NDG_NEW_dims
         integer, dimension ( T_NDG_NEW_rank) ::  T_NDG_NEW_dims
         integer, dimension ( Q_NDG_NEW_rank) ::  Q_NDG_NEW_dims
         integer, dimension (PH_NDG_NEW_rank) :: PH_NDG_NEW_dims
         integer, dimension (MU_NDG_NEW_rank) :: MU_NDG_NEW_dims

         !------------------
         ! Attribute vectors
         !------------------
         integer :: intval (1)
         real    :: realval(1)

         !-----------------------------
         ! Size of data chunks to write
         !-----------------------------
         integer, dimension(1), parameter :: start1d = (/ 1 /)
         integer, dimension(2), parameter :: start2d = (/ 1, 1 /)
!        integer, dimension(3), parameter :: start3d = (/ 1, 1, 1 /)
         integer, dimension(4), parameter :: start4d = (/ 1, 1, 1, 1 /)

         integer, dimension(1), parameter :: count1d = (/ 1 /)
         integer, dimension(2)            :: count2d
!        integer, dimension(3)            :: count3d
         integer, dimension(4)            :: count4d

         count2d(1) = DateStrLen
         count2d(2) = Time

         count4d(1) = west_east
         count4d(2) = south_north
         count4d(3) = bottom_top
         count4d(4) = Time

         !-------------------------------------------------
         ! Create a NetCDF output file -- enter define mode
         !-------------------------------------------------
         write (file_name,'(a, i2.2)') file_prefix, grid_id
         if (deBug) print *, "writeNudge: creating NetCDF output file '", file_name, "'"
         call check( nf90_create(file_name, nf90_clobber, ncID) )
    !    call check( nf90_create(file_name, nf90_netcdf4, ncID) )

         !---------------------------------------------------------
         ! Optimize write efficiency by deactivating fill (default)
         !---------------------------------------------------------
         call check( nf90_set_fill(ncID, nf90_noFill, oldFillMode) )   ! <- Craig suggested adding this as part of the 'no-prefill' method of speeding up the code

         !------------------------------------------
         ! Generate a file creation date (timestamp)
         !------------------------------------------
         creation_date = now()
         if (deBug) print *, "writeNudge: file creation date is ", creation_date

         !------------------------------------------------
         ! Define dimensions.
         ! The record (time) dimension is defined to have
         ! unlimited length - it can grow as needed.
         !------------------------------------------------
         if (deBug) print *, "writeNudge: defining dimensions"

         call check( nf90_def_dim(ncID, "Time"       , nf90_unlimited , Time_dim       ) )
!        call check( nf90_def_dim(ncID, "Time"       , Time           , Time_dim       ) )
         call check( nf90_def_dim(ncID, "DateStrLen" , DateStrLen     , DateStrLen_dim ) )
         call check( nf90_def_dim(ncID, "west_east"  , west_east      , west_east_dim  ) )
         call check( nf90_def_dim(ncID, "south_north", south_north    , south_north_dim) )
         call check( nf90_def_dim(ncID, "bottom_top" , bottom_top     , bottom_top_dim ) )
         call check( nf90_def_dim(ncID, "one_stag"   , one_stag       , one_stag_dim   ) )

         !-----------------
         ! Define variables
         !-----------------
         if (deBug) print *, "writeNudge: defining variables"

         Times_dims(1) = DateStrLen_dim
         Times_dims(2) = Time_dim
         call check( nf90_def_var(ncID, "Times", nf90_char, Times_dims, Times_id) )

         U_NDG_OLD_dims(1) = west_east_dim
         U_NDG_OLD_dims(2) = south_north_dim
         U_NDG_OLD_dims(3) = bottom_top_dim
         U_NDG_OLD_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "U_NDG_OLD", nf90_float, U_NDG_OLD_dims, U_NDG_OLD_id) )

         V_NDG_OLD_dims(1) = west_east_dim
         V_NDG_OLD_dims(2) = south_north_dim
         V_NDG_OLD_dims(3) = bottom_top_dim
         V_NDG_OLD_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "V_NDG_OLD", nf90_float, V_NDG_OLD_dims, V_NDG_OLD_id) )

         T_NDG_OLD_dims(1) = west_east_dim
         T_NDG_OLD_dims(2) = south_north_dim
         T_NDG_OLD_dims(3) = bottom_top_dim
         T_NDG_OLD_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "T_NDG_OLD", nf90_float, T_NDG_OLD_dims, T_NDG_OLD_id) )

         Q_NDG_OLD_dims(1) = west_east_dim
         Q_NDG_OLD_dims(2) = south_north_dim
         Q_NDG_OLD_dims(3) = bottom_top_dim
         Q_NDG_OLD_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "Q_NDG_OLD", nf90_float, Q_NDG_OLD_dims, Q_NDG_OLD_id) )

         PH_NDG_OLD_dims(1) = west_east_dim
         PH_NDG_OLD_dims(2) = south_north_dim
         PH_NDG_OLD_dims(3) = bottom_top_dim
         PH_NDG_OLD_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "PH_NDG_OLD", nf90_float, PH_NDG_OLD_dims, PH_NDG_OLD_id) )

         MU_NDG_OLD_dims(1) = west_east_dim
         MU_NDG_OLD_dims(2) = south_north_dim
         MU_NDG_OLD_dims(3) = one_stag_dim
         MU_NDG_OLD_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "MU_NDG_OLD", nf90_float, MU_NDG_OLD_dims, MU_NDG_OLD_id) )

         U_NDG_NEW_dims(1) = west_east_dim
         U_NDG_NEW_dims(2) = south_north_dim
         U_NDG_NEW_dims(3) = bottom_top_dim
         U_NDG_NEW_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "U_NDG_NEW", nf90_float, U_NDG_NEW_dims, U_NDG_NEW_id) )

         V_NDG_NEW_dims(1) = west_east_dim
         V_NDG_NEW_dims(2) = south_north_dim
         V_NDG_NEW_dims(3) = bottom_top_dim
         V_NDG_NEW_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "V_NDG_NEW", nf90_float,  V_NDG_NEW_dims, V_NDG_NEW_id) )

         T_NDG_NEW_dims(1) = west_east_dim
         T_NDG_NEW_dims(2) = south_north_dim
         T_NDG_NEW_dims(3) = bottom_top_dim
         T_NDG_NEW_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "T_NDG_NEW", nf90_float, T_NDG_NEW_dims, T_NDG_NEW_id) )

         Q_NDG_NEW_dims(1) = west_east_dim
         Q_NDG_NEW_dims(2) = south_north_dim
         Q_NDG_NEW_dims(3) = bottom_top_dim
         Q_NDG_NEW_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "Q_NDG_NEW", nf90_float, Q_NDG_NEW_dims, Q_NDG_NEW_id) )

         PH_NDG_NEW_dims(1) = west_east_dim
         PH_NDG_NEW_dims(2) = south_north_dim
         PH_NDG_NEW_dims(3) = bottom_top_dim
         PH_NDG_NEW_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "PH_NDG_NEW", nf90_float, PH_NDG_NEW_dims, PH_NDG_NEW_id) )

         MU_NDG_NEW_dims(1) = west_east_dim
         MU_NDG_NEW_dims(2) = south_north_dim
         MU_NDG_NEW_dims(3) = one_stag_dim
         MU_NDG_NEW_dims(4) = Time_dim
         call check( nf90_def_var(ncID, "MU_NDG_NEW", nf90_float, MU_NDG_NEW_dims, MU_NDG_NEW_id) )

         !-------------------------
         ! Assign global attributes
         !-------------------------
         if (deBug) print *, "writeNudge: assigning global attributes"

         call check( nf90_put_att(ncID, nf90_global, "TITLE"                , title                ) )
         call check( nf90_put_att(ncID, nf90_global, "INSTITUTION"          , institution          ) )
         call check( nf90_put_att(ncID, nf90_global, "CONTACT"              , contact              ) )
         call check( nf90_put_att(ncID, nf90_global, "SOURCE"               , source               ) )
         call check( nf90_put_att(ncID, nf90_global, "HISTORY"              , history              ) )
         call check( nf90_put_att(ncID, nf90_global, "REFERENCES"           , references           ) )
         call check( nf90_put_att(ncID, nf90_global, "MODEL"                , model                ) )
         call check( nf90_put_att(ncID, nf90_global, "VERSION"              , version              ) )
         call check( nf90_put_att(ncID, nf90_global, "CONVENTIONS"          , Conventions          ) )
         call check( nf90_put_att(ncID, nf90_global, "CREATION_DATE"        , creation_date        ) )
         call check( nf90_put_att(ncID, nf90_global, "START_DATE"           , start_date           ) )
         call check( nf90_put_att(ncID, nf90_global, "SIMULATION_START_DATE", simulation_start_date) )

         intval(1) = west_east + 1
         call check( nf90_put_att(ncID, nf90_global, "WEST-EAST_GRID_DIMENSION", intval) )

         intval(1) = south_north + 1
         call check( nf90_put_att(ncID, nf90_global, "SOUTH-NORTH_GRID_DIMENSION", intval) )

         intval(1) = bottom_top + 1
         call check( nf90_put_att(ncID, nf90_global, "BOTTOM-TOP_GRID_DIMENSION", intval) )

         realval(1) = dx
         call check( nf90_put_att(ncID, nf90_global, "DX", realval) )

         realval(1) = dy
         call check( nf90_put_att(ncID, nf90_global, "DY", realval) )

         call check( nf90_put_att(ncID, nf90_global, "GRIDTYPE", 'C') )

         intval(1) = diff_opt
         call check( nf90_put_att(ncID, nf90_global, "DIFF_OPT", intval) )

         intval(1) = km_opt
         call check( nf90_put_att(ncID, nf90_global, "KM_OPT", intval) )

         intval(1) = damp_opt
         call check( nf90_put_att(ncID, nf90_global, "DAMP_OPT", intval) )

         realval(1) = dampcoef
         call check( nf90_put_att(ncID, nf90_global, "DAMPCOEF", realval) )

         realval(1) = khdif
         call check( nf90_put_att(ncID, nf90_global, "KHDIF", realval) )

         realval(1) = kvdif
         call check( nf90_put_att(ncID, nf90_global, "KVDIF", realval) )

         intval(1) = mp_physics
         call check( nf90_put_att(ncID, nf90_global, "MP_PHYSICS", intval) )

         intval(1) = ra_lw_physics
         call check( nf90_put_att(ncID, nf90_global, "RA_LW_PHYSICS", intval) )

         intval(1) = ra_sw_physics
         call check( nf90_put_att(ncID, nf90_global, "RA_SW_PHYSICS", intval) )

         intval(1) = sf_sfclay_physics
         call check( nf90_put_att(ncID, nf90_global, "SF_SFCLAY_PHYSICS", intval) )

         intval(1) = sf_surface_physics
         call check( nf90_put_att(ncID, nf90_global, "SF_SURFACE_PHYSICS", intval) )

         intval(1) = bl_pbl_physics
         call check( nf90_put_att(ncID, nf90_global, "BL_PBL_PHYSICS", intval) )

         intval(1) = cu_physics
         call check( nf90_put_att(ncID, nf90_global, "CU_PHYSICS", intval) )

         intval(1) = surface_input_source
         call check( nf90_put_att(ncID, nf90_global, "SURFACE_INPUT_SOURCE", intval) )

         intval(1) = sst_update
         call check( nf90_put_att(ncID, nf90_global, "SST_UPDATE", intval) )

         intval(1) = grid_fdda
         call check( nf90_put_att(ncID, nf90_global, "GRID_FDDA", intval) )

         intval(1) = gfdda_interval_m
         call check( nf90_put_att(ncID, nf90_global, "GFDDA_INTERVAL_M", intval) )

         intval(1) = gfdda_end_h
         call check( nf90_put_att(ncID, nf90_global, "GFDDA_END_H", intval) )

         intval(1) = grid_sfdda
         call check( nf90_put_att(ncID, nf90_global, "GRID_SFDDA", intval) )

         intval(1) = sgfdda_interval_m
         call check( nf90_put_att(ncID, nf90_global, "SGFDDA_INTERVAL_M", intval) )

         intval(1) = sgfdda_end_h
         call check( nf90_put_att(ncID, nf90_global, "SGFDDA_END_H", intval) )

         intval(1) = hypsometric_opt
         call check( nf90_put_att(ncID, nf90_global, "HYPSOMETRIC_OPT", intval) )

         intval(1) = 1
         call check( nf90_put_att(ncID, nf90_global, "WEST-EAST_PATCH_START_UNSTAG", intval) )

         intval(1) = west_east
         call check( nf90_put_att(ncID, nf90_global, "WEST-EAST_PATCH_END_UNSTAG", intval) )

         intval(1) = 1
         call check( nf90_put_att(ncID, nf90_global, "WEST-EAST_PATCH_START_STAG", intval) )

         intval(1) = west_east + 1
         call check( nf90_put_att(ncID, nf90_global, "WEST-EAST_PATCH_END_STAG", intval) )

         intval(1) = 1
         call check( nf90_put_att(ncID, nf90_global, "SOUTH-NORTH_PATCH_START_UNSTAG", intval) )

         intval(1) = south_north
         call check( nf90_put_att(ncID, nf90_global, "SOUTH-NORTH_PATCH_END_UNSTAG", intval) )

         intval(1) = 1
         call check( nf90_put_att(ncID, nf90_global, "SOUTH-NORTH_PATCH_START_STAG", intval) )

         intval(1) = south_north + 1
         call check( nf90_put_att(ncID, nf90_global, "SOUTH-NORTH_PATCH_END_STAG", intval) )

         intval(1) = 1
         call check( nf90_put_att(ncID, nf90_global, "BOTTOM-TOP_PATCH_START_UNSTAG", intval) )

         intval(1) = bottom_top
         call check( nf90_put_att(ncID, nf90_global, "BOTTOM-TOP_PATCH_END_UNSTAG", intval) )

         intval(1) = 1
         call check( nf90_put_att(ncID, nf90_global, "BOTTOM-TOP_PATCH_START_STAG", intval) )

         intval(1) = bottom_top + 1
         call check( nf90_put_att(ncID, nf90_global, "BOTTOM-TOP_PATCH_END_STAG", intval) )

         intval(1) = grid_id
         call check( nf90_put_att(ncID, nf90_global, "GRID_ID", intval) )

         intval(1) = parent_id
         call check( nf90_put_att(ncID, nf90_global, "PARENT_ID", intval) )

         intval(1) = i_parent_start
         call check( nf90_put_att(ncID, nf90_global, "I_PARENT_START", intval) )

         intval(1) = j_parent_start
         call check( nf90_put_att(ncID, nf90_global, "J_PARENT_START", intval) )

         intval(1) = parent_grid_ratio
         call check( nf90_put_att(ncID, nf90_global, "PARENT_GRID_RATIO", intval) )

         realval(1) = dt
         call check( nf90_put_att(ncID, nf90_global, "DT", realval) )

         realval(1) = cen_lat
         call check( nf90_put_att(ncID, nf90_global, "CEN_LAT", realval) )

         realval(1) = cen_lon
         call check( nf90_put_att(ncID, nf90_global, "CEN_LON", realval) )

         realval(1) = truelat1
         call check( nf90_put_att(ncID, nf90_global, "TRUELAT1", realval) )

         realval(1) = truelat2
         call check( nf90_put_att(ncID, nf90_global, "TRUELAT2", realval) )

         realval(1) = moad_cen_lat
         call check( nf90_put_att(ncID, nf90_global, "MOAD_CEN_LAT", realval) )

         realval(1) = stand_lon
         call check( nf90_put_att(ncID, nf90_global, "STAND_LON", realval) )

         realval(1) = pole_lat
         call check( nf90_put_att(ncID, nf90_global, "POLE_LAT", realval) )

         realval(1) = pole_lon
         call check( nf90_put_att(ncID, nf90_global, "POLE_LON", realval) )

         intval(1) = map_proj
         call check( nf90_put_att(ncID, nf90_global, "MAP_PROJ", intval) )

         call check( nf90_put_att(ncID, nf90_global, "MMINLU", mminlu) )

         intval(1) = num_land_cat
         call check( nf90_put_att(ncID, nf90_global, "NUM_LAND_CAT", intval) )

         intval(1) = isWater
         call check( nf90_put_att(ncID, nf90_global, "ISWATER", intval) )

         intval(1) = isLake
         call check( nf90_put_att(ncID, nf90_global, "ISLAKE", intval) )

         intval(1) = isIce
         call check( nf90_put_att(ncID, nf90_global, "ISICE", intval) )

         intval(1) = isUrban
         call check( nf90_put_att(ncID, nf90_global, "ISURBAN", intval) )

         intval(1) = iSoilWater
         call check( nf90_put_att(ncID, nf90_global, "ISOILWATER", intval) )

         !-------------------------------
         ! Assign per-variable attributes
         !-------------------------------
         intval(1) = fieldType
         call check( nf90_put_att(ncID, U_NDG_OLD_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, U_NDG_OLD_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, U_NDG_OLD_id, "description", "OLD X WIND FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, U_NDG_OLD_id, "units"      , "m s-1") )
         call check( nf90_put_att(ncID, U_NDG_OLD_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, U_NDG_OLD_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, V_NDG_OLD_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, V_NDG_OLD_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, V_NDG_OLD_id, "description", "OLD Y WIND FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, V_NDG_OLD_id, "units"      , "m s-1") )
         call check( nf90_put_att(ncID, V_NDG_OLD_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, V_NDG_OLD_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, T_NDG_OLD_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, T_NDG_OLD_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, T_NDG_OLD_id, "description", "OLD PERT POT TEMP FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, T_NDG_OLD_id, "units"      , "K") )
         call check( nf90_put_att(ncID, T_NDG_OLD_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, T_NDG_OLD_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, Q_NDG_OLD_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, Q_NDG_OLD_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, Q_NDG_OLD_id, "description", "OLD WATER VAPOR MIX RATIO FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, Q_NDG_OLD_id, "units"      , "kg/kg") )
         call check( nf90_put_att(ncID, Q_NDG_OLD_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, Q_NDG_OLD_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, PH_NDG_OLD_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, PH_NDG_OLD_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, PH_NDG_OLD_id, "description", "OLD PERT GEOPOTENTIAL FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, PH_NDG_OLD_id, "units"      , "kg/kg") )
         call check( nf90_put_att(ncID, PH_NDG_OLD_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, PH_NDG_OLD_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, MU_NDG_OLD_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, MU_NDG_OLD_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, MU_NDG_OLD_id, "description", "OLD PERT COLUMN DRY MASS FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, MU_NDG_OLD_id, "units"      , "Pa") )
         call check( nf90_put_att(ncID, MU_NDG_OLD_id, "stagger"    , "Z") )
         call check( nf90_put_att(ncID, MU_NDG_OLD_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, U_NDG_NEW_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, U_NDG_NEW_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, U_NDG_NEW_id, "description", "NEW X WIND FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, U_NDG_NEW_id, "units"      , "m s-1") )
         call check( nf90_put_att(ncID, U_NDG_NEW_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, U_NDG_NEW_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, V_NDG_NEW_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, V_NDG_NEW_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, V_NDG_NEW_id, "description", "NEW Y WIND FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, V_NDG_NEW_id, "units"      , "m s-1") )
         call check( nf90_put_att(ncID, V_NDG_NEW_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, V_NDG_NEW_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, T_NDG_NEW_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, T_NDG_NEW_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, T_NDG_NEW_id, "description", "NEW PERT POT TEMP FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, T_NDG_NEW_id, "units"      , "K") )
         call check( nf90_put_att(ncID, T_NDG_NEW_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, T_NDG_NEW_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, Q_NDG_NEW_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, Q_NDG_NEW_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, Q_NDG_NEW_id, "description", "NEW WATER VAPOR MIX RATIO FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, Q_NDG_NEW_id, "units"      , "kg/kg") )
         call check( nf90_put_att(ncID, Q_NDG_NEW_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, Q_NDG_NEW_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, PH_NDG_NEW_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, PH_NDG_NEW_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, PH_NDG_NEW_id, "description", "NEW PERT GEOPOTENTIAL FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, PH_NDG_NEW_id, "units"      , "kg/kg") )
         call check( nf90_put_att(ncID, PH_NDG_NEW_id, "stagger"    , char(0)) )
         call check( nf90_put_att(ncID, PH_NDG_NEW_id, "coordinates", "XLONG XLAT") )

!        intval(1) = fieldType
         call check( nf90_put_att(ncID, MU_NDG_NEW_id, "FieldType"  , intval) )
         call check( nf90_put_att(ncID, MU_NDG_NEW_id, "MemoryOrder", memoryOrder) )
         call check( nf90_put_att(ncID, MU_NDG_NEW_id, "description", "NEW PERT COLUMN DRY MASS FOR FDDA GRID NUDGING") )
         call check( nf90_put_att(ncID, MU_NDG_NEW_id, "units"      , "Pa") )
         call check( nf90_put_att(ncID, MU_NDG_NEW_id, "stagger"    , "Z") )
         call check( nf90_put_att(ncID, MU_NDG_NEW_id, "coordinates", "XLONG XLAT") )

         !------------------
         ! Leave define mode
         !------------------
         if (deBug) print *, "writeNudge: leaving NetCDF define mode"
         call check( nf90_enddef(ncID) )

         !---------------------------------
         ! Write data to NetCDF output file
         !---------------------------------
         if (deBug) print *, "writeNudge: writing Times"
         call check( nf90_put_var(ncID, Times_id, Times, start=start2d, count=count2d) )

         if (deBug) print *, "writeNudge: writing OLD data for WRF FDDA grid nudging"
!        call check( nf90_put_var(ncID,  U_NDG_OLD_id,  U_NDG_OLD, start=start4d, count=count4d) )
         call check( nf90_put_var(ncID,  U_NDG_OLD_id,  U_NDG_OLD) )
         call check( nf90_put_var(ncID,  V_NDG_OLD_id,  V_NDG_OLD) )
         call check( nf90_put_var(ncID,  T_NDG_OLD_id,  T_NDG_OLD) )
         call check( nf90_put_var(ncID,  Q_NDG_OLD_id,  Q_NDG_OLD) )
         call check( nf90_put_var(ncID, PH_NDG_OLD_id, PH_NDG_OLD) )
         call check( nf90_put_var(ncID, MU_NDG_OLD_id, MU_NDG_OLD) )

         if (deBug) print *, "writeNudge: writing NEW data for WRF FDDA grid nudging"
         call check( nf90_put_var(ncID,  U_NDG_NEW_id,  U_NDG_NEW) )
         call check( nf90_put_var(ncID,  V_NDG_NEW_id,  V_NDG_NEW) )
         call check( nf90_put_var(ncID,  T_NDG_NEW_id,  T_NDG_NEW) )
         call check( nf90_put_var(ncID,  Q_NDG_NEW_id,  Q_NDG_NEW) )
         call check( nf90_put_var(ncID, PH_NDG_NEW_id, PH_NDG_NEW) )
         call check( nf90_put_var(ncID, MU_NDG_NEW_id, MU_NDG_NEW) )

         !-------------------------
         ! Close NetCDF output file
         !-------------------------
			call check( nf90_close(ncID) )
         if (deBug) print *, "writeNudge: all data written, NetCDF file closed"
      end subroutine writeNudge

      !=================================================================
      ! This subroutine reads WRF input files, then creates a WRF grid/
      ! analysis nudging file in NetCDF format.
      !
      ! Input:
      !    none
      !
      ! Output:
      !    WRF nudging file named wrffdda_d<nn> created
      !
      ! Notes:
      !    Ramping of values at end of nudging interval is controlled by
      !       global nudging control parameters 'rampDown' and 'rampTime'
      !       at the top of this module
      !=================================================================
      subroutine nudge()
        !--------------------------------------------
         ! Force explicit declaration of all variables
         !--------------------------------------------
         implicit none

         !---------------
         ! File variables
         !---------------
         character (len = 16), dimension(:), allocatable :: files
         integer :: nf, n, nm1, nudge_int_min

         !------------------
         ! Dimension lengths
         !------------------
         integer :: west_east
         integer :: south_north
         integer :: bottom_top
         integer :: Time
         integer, parameter :: one_stag = 1

         !------------------------
         ! Date and time variables
         !------------------------
         integer :: DateStrLen
         integer, parameter :: dateStringLength =  19
         character (len = dateStringLength) start_date
         character (len = dateStringLength) simulation_start_date
         character (len = dateStringLength) date, ss_date
         character (len = dateStringLength), allocatable, dimension(:) :: allTimes
         real(8) :: timeWindow
         integer :: nRamp

         !---------------
         ! WRF input data
         !---------------
         character (len = dateStringLength), allocatable, dimension(:) :: Times
         real, allocatable, dimension(:,:,:,:) ::  U
         real, allocatable, dimension(:,:,:,:) ::  V
         real, allocatable, dimension(:,:,:,:) ::  T
         real, allocatable, dimension(:,:,:,:) ::  Q
         real, allocatable, dimension(:,:,:,:) :: PH
         real, allocatable, dimension(:,:,  :) :: MU

         !-----------------
         ! WRF nudging data
         !-----------------
         real, allocatable, dimension(:,:,:,:) ::  U_NDG_OLD
         real, allocatable, dimension(:,:,:,:) ::  V_NDG_OLD
         real, allocatable, dimension(:,:,:,:) ::  T_NDG_OLD
         real, allocatable, dimension(:,:,:,:) ::  Q_NDG_OLD
         real, allocatable, dimension(:,:,:,:) :: PH_NDG_OLD
         real, allocatable, dimension(:,:,:,:) :: MU_NDG_OLD

         real, allocatable, dimension(:,:,:,:) ::  U_NDG_NEW
         real, allocatable, dimension(:,:,:,:) ::  V_NDG_NEW
         real, allocatable, dimension(:,:,:,:) ::  T_NDG_NEW
         real, allocatable, dimension(:,:,:,:) ::  Q_NDG_NEW
         real, allocatable, dimension(:,:,:,:) :: PH_NDG_NEW
         real, allocatable, dimension(:,:,:,:) :: MU_NDG_NEW

         !----------------------------------------------
         ! Read list of WRF input files provided by user
         !----------------------------------------------
         call userInput (nf, files, nudge_int_min)

         if (deBug) then
            write (*, '(i3, " WRF input files:")') nf
            do n = 1, nf
               write (*, '(a)') files(n)
            end do
         end if

         !----------------------------------
         ! Allocate memory for nudging times
         !----------------------------------
         if (.not. allocated(allTimes)) allocate(allTimes(nf))

         !--------------------------
         ! Loop over WRF input files
         !--------------------------
         do n = 1, nf
            !-------------------------
            ! Read each WRF input file
            !-------------------------
            call readWRFinput (n, trim(files(n)),                       &
                               west_east, south_north, bottom_top, Time,&
                               DateStrLen,                              &
                               date, ss_date,                           &
                               Times,                                   &
                               U, V, T, Q, PH, MU)

            !-----------------------------------
            ! Allocate memory for nudging arrays
            !-----------------------------------
            if (n == 1) then
               if (.not. allocated( U_NDG_OLD)) allocate( U_NDG_OLD(west_east, south_north, bottom_top, nf))
               if (.not. allocated( V_NDG_OLD)) allocate( V_NDG_OLD(west_east, south_north, bottom_top, nf))
               if (.not. allocated( T_NDG_OLD)) allocate( T_NDG_OLD(west_east, south_north, bottom_top, nf))
               if (.not. allocated( Q_NDG_OLD)) allocate( Q_NDG_OLD(west_east, south_north, bottom_top, nf))
               if (.not. allocated(PH_NDG_OLD)) allocate(PH_NDG_OLD(west_east, south_north, bottom_top, nf))
               if (.not. allocated(MU_NDG_OLD)) allocate(MU_NDG_OLD(west_east, south_north, one_stag  , nf))

               if (.not. allocated( U_NDG_NEW)) allocate( U_NDG_NEW(west_east, south_north, bottom_top, nf))
               if (.not. allocated( V_NDG_NEW)) allocate( V_NDG_NEW(west_east, south_north, bottom_top, nf))
               if (.not. allocated( T_NDG_NEW)) allocate( T_NDG_NEW(west_east, south_north, bottom_top, nf))
               if (.not. allocated( Q_NDG_NEW)) allocate( Q_NDG_NEW(west_east, south_north, bottom_top, nf))
               if (.not. allocated(PH_NDG_NEW)) allocate(PH_NDG_NEW(west_east, south_north, bottom_top, nf))
               if (.not. allocated(MU_NDG_NEW)) allocate(MU_NDG_NEW(west_east, south_north, one_stag  , nf))
            end if

            !--------------------------------------------
            ! Increment date + time for each nudging time
            !--------------------------------------------
            if (n == 1) then
               start_date = date
               simulation_start_date = ss_date
               gfdda_interval_m = nudge_int_min
               timeWindow = gfdda_interval_m * 60.d0 ! Nudging time window (seconds)
               allTimes(1) = Times(1)
            else
               allTimes(n) = upDate(allTimes(1), (n - 1)*timeWindow)
            end if
            if (deBug) write(*, '("nudge: allTimes(",i2,") = ", a19)') n, allTimes(n)

            !--------------------------------------------------------
            ! Append WRF input data to nudging arrays.
            ! Note: assumes only 1 time level in each WRF input file.
            !--------------------------------------------------------
            nm1 = max(n-1, 1)

             U_NDG_OLD(:,:,:,n  ) =  U(:,:,:,1)
             V_NDG_OLD(:,:,:,n  ) =  V(:,:,:,1)
             T_NDG_OLD(:,:,:,n  ) =  T(:,:,:,1)
             Q_NDG_OLD(:,:,:,n  ) =  Q(:,:,:,1)
            PH_NDG_OLD(:,:,:,n  ) = PH(:,:,:,1)
            MU_NDG_OLD(:,:,1,n  ) = MU(:,:,  1)

             U_NDG_NEW(:,:,:,nm1) =  U(:,:,:,1)
             V_NDG_NEW(:,:,:,nm1) =  V(:,:,:,1)
             T_NDG_NEW(:,:,:,nm1) =  T(:,:,:,1)
             Q_NDG_NEW(:,:,:,nm1) =  Q(:,:,:,1)
            PH_NDG_NEW(:,:,:,nm1) = PH(:,:,:,1)
            MU_NDG_NEW(:,:,1,nm1) = MU(:,:,  1)

            if (n == nf) then
               if (rampDown) then
                  !------------------------------------
                  ! Ramp nudging down to initial values
                  ! at end of simulation
                  !------------------------------------
                  nRamp = min( max(rampTime, 1), nf )
                   U_NDG_NEW(:,:,:,n) =  U_NDG_OLD(:,:,:,nRamp)
                   V_NDG_NEW(:,:,:,n) =  V_NDG_OLD(:,:,:,nRamp)
                   T_NDG_NEW(:,:,:,n) =  T_NDG_OLD(:,:,:,nRamp)
                   Q_NDG_NEW(:,:,:,n) =  Q_NDG_OLD(:,:,:,nRamp)
                  PH_NDG_NEW(:,:,:,n) = PH_NDG_OLD(:,:,:,nRamp)
                  MU_NDG_NEW(:,:,1,n) = MU_NDG_OLD(:,:,1,nRamp)
               end if
            end if
         end do   ! n = 1, nf

         !---------------------------------------
         ! Adjust size and content of Times array
         !---------------------------------------
         if (allocated(Times)) deallocate(Times)
         allocate(Times(nf))
         Times(:) = allTimes(1:nf)

         !---------------------------------------------
         ! Write data to WRF analysis/grid nudging file
         !---------------------------------------------
         call writeNudge (west_east, south_north, bottom_top, nf,       &
                          DateStrLen,                                   &
                          start_date, simulation_start_date,            &
                          Times,                                        &
                           U_NDG_OLD,  U_NDG_NEW,                       &
                           V_NDG_OLD,  V_NDG_NEW,                       &
                           T_NDG_OLD,  T_NDG_NEW,                       &
                           Q_NDG_OLD,  Q_NDG_NEW,                       &
                          PH_NDG_OLD, PH_NDG_NEW,                       &
                          MU_NDG_OLD, MU_NDG_NEW)

         !------------
         ! Free memory
         !------------
         if (allocated(Times))      deallocate(Times)

         if (allocated( U_NDG_OLD)) deallocate( U_NDG_OLD)
         if (allocated( V_NDG_OLD)) deallocate( V_NDG_OLD)
         if (allocated( T_NDG_OLD)) deallocate( T_NDG_OLD)
         if (allocated( Q_NDG_OLD)) deallocate( Q_NDG_OLD)
         if (allocated(PH_NDG_OLD)) deallocate(PH_NDG_OLD)
         if (allocated(MU_NDG_OLD)) deallocate(MU_NDG_OLD)

         if (allocated( U_NDG_NEW)) deallocate( U_NDG_NEW)
         if (allocated( V_NDG_NEW)) deallocate( V_NDG_NEW)
         if (allocated( T_NDG_NEW)) deallocate( T_NDG_NEW)
         if (allocated( Q_NDG_NEW)) deallocate( Q_NDG_NEW)
         if (allocated(PH_NDG_NEW)) deallocate(PH_NDG_NEW)
         if (allocated(MU_NDG_NEW)) deallocate(MU_NDG_NEW)

         if (allocated( U)) deallocate( U)
         if (allocated( V)) deallocate( V)
         if (allocated( T)) deallocate( T)
         if (allocated( Q)) deallocate( Q)
         if (allocated(PH)) deallocate(PH)
         if (allocated(MU)) deallocate(MU)

         if (allocated(files)) deallocate(files)
      end subroutine nudge

      !=============================================================
      ! Check for NetCDF errors.
      !
      ! Input:
      !     anyErrors   NetCDF error code
      !
      ! Output:
      !     If an error is encountered, prints NetCDF error code and
      !     explanation, then terminates program; otherwise does
      !     nothing.
      !=============================================================
      subroutine check (anyErrors)
         use netcdf
         implicit none
         integer, intent(in) :: anyErrors

         if (anyErrors /= nf90_noerr) then
            write (*,*) " Error number ", anyErrors
            write (*,*) trim(nf90_strerror(anyErrors))
            stop 2
         end if
      end subroutine check

      !============================================================
      ! This function creates a WRF-standard date character string.
      !
      ! Input:
      !    yyyy      year    (4 digits)
      !    mm        month   (2 digits)
      !    dd        day     (2 digits)
      !    hh        hour    (2 digits)
      !    mn        minutes (2 digits)
      !    ss        seconds (2 digits)
      !
      ! Output:
      !    wrfDate   Date string (yyyy-mm-dd_hh:mn:ss)
      !============================================================
      elemental character(len = 19) function wrfDate(yyyy, mm, dd, hh, mn, ss)
         integer, intent(in) :: yyyy
         integer, intent(in) :: mm
         integer, intent(in) :: dd
         integer, intent(in) :: hh
         integer, intent(in) :: mn
         integer, intent(in) :: ss
         write (unit = wrfDate,														&
         		 fmt = '(i4.4, 2("-",i2.2), "_", i2.2, 2(":",i2.2))')		&
               yyyy, mm, dd, hh, mn, ss
      end function wrfDate

      !===========================================================
      ! This function creates a WRF-standard date character string
      ! that represents the current time.
      !
      ! Input:
      !    none
      !
      ! Output:
      !    now       date string (yyyy-mm-dd_hh:mn:ss)
      !===========================================================
      character(len = 19) function now()
         character (len =  8) :: d
         character (len = 10) :: t
         call date_and_time(date = d, time = t)
         now = d(1:4) // '-' // d(5:6) // '-' // d(7:8) // '_' //       &
               t(1:2) // ':' // t(3:4) // ':' // t(5:6)
      end function now

      !=================================================================
      ! This function increments a WRF date by a given amount of time.
      !
      ! Input:
      !    befDate   Date before incrementing time (yyyy-mm-dd_hh:mn:ss)
      !    dt        Time increment (seconds)
      !
      ! Output:
      !    upDate    Date after  incrementing time (yyyy-mm-dd_hh:mn:ss)
      !=================================================================
      character(len = 19) function upDate(befDate, dt)
         character (len = 19), intent(in) :: befDate
         real(8)             , intent(in) :: dt
         real(8) :: aftTime
         integer :: yRef, yyyy, mm, dd, hh, mn, ss

         call extractDate(befDate, yyyy, mm, dd, hh, mn, ss)
         yRef = yyyy
         aftTime = d2s(yyyy, mm, dd, hh, mn, ss, yRef) + dt
         call s2d(aftTime, yRef, yyyy, mm, dd, hh, mn, ss)
         upDate = wrfDate(yyyy, mm, dd, hh, mn, ss)
      end function upDate

      !==============================================================
      ! This subroutine extracts year, month, day, hour, minutes, and
      ! seconds from a date character string.
      !
      ! Input:
      !     date     Date string in WRF format (yyyy-mm-dd_hh:mn:ss)
      !
      ! Output:
      !     yyyy     year    (4 digits)
      !     mm       month   (2 digits)
      !     dd       day     (2 digits)
      !     hh       hour    (2 digits)
      !     mn       minutes (2 digits)
      !     ss       seconds (2 digits)
      !==============================================================
      subroutine extractDate(date, yyyy, mm, dd, hh, mn, ss)
         character(19), intent(in ) :: date
         integer      , intent(out) :: yyyy
         integer      , intent(out) :: mm
         integer      , intent(out) :: dd
         integer      , intent(out) :: hh
         integer      , intent(out) :: mn
         integer      , intent(out) :: ss
         read (date, '(i4.4, 5(1x,i2.2))') yyyy, mm, dd, hh, mn, ss
      end subroutine extractDate

      !===================================================
      ! This subroutine converts time, given in cumulative
      ! seconds since year 'yRef', to a Gregorian date.
      ! NOTE: This is the bug-fixed version.
      !
      ! Input:
      !    cs        Cumulative seconds
      !    yRef      yRef    (4 digits) = reference year
      !
      ! Output:
      !    Gregorian date:
      !    yyyy      year    (4 digits)
      !    mm        month   (2 digits)
      !    dd        day     (2 digits)
      !    hh        hour    (2 digits)
      !    mn        minutes (2 digits)
      !    ss        seconds (2 digits)
      !===================================================
      subroutine s2d(cs, yRef, yyyy, mm, dd, hh, mn, ss)
         real*8  , intent(in ) :: cs
         integer , intent(in ) :: yRef
         integer , intent(out) :: yyyy
         integer , intent(out) :: mm
         integer , intent(out) :: dd
         integer , intent(out) :: hh
         integer , intent(out) :: mn
         integer , intent(out) :: ss

         integer   :: jdy
         integer   :: monthDays(12)
         integer   :: month
         integer*8 :: i
         integer*8 :: js
         real*8    :: ds
         real*8    :: sec

         monthDays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

         yyyy = yRef
         sec  = cs
         js   = sec
         ds   = sec - js
         i    = js
         do while(i > 0)
            if (isLeapYear(yyyy)) then
               i = i - 366 * 86400
            else
               i = i - 365 * 86400
            end if
            if (i > 0) then
               yyyy = yyyy + 1
               js = i               ! Seconds remaining after years counted
            end if
         end do

         jdy = int(js / 86400.)
         js  = js - jdy * 86400.    ! Seconds remaining after days counted
         jdy = jdy + 1              ! Increment day to ordinal date

         hh = int(js / 3600.)
         js = js - hh * 3600.       ! Seconds remaining after hours counted

         mn = int(js / 60.)
         js = js - mn * 60.         ! Seconds remaining after minutes counted
         ss = js

         sec = ds                   ! Fraction of second remaining (if any)

         if (isLeapYear(yyyy)) then ! Adjust February for leap year
            monthDays(2) = 29
         else
            monthDays(2) = 28
         end if

         mm = 1                     ! Compute month and day
         dd = jdy
         do month = 1, 12
            if (jdy > monthDays(month)) then
               jdy = jdy - monthDays(month)
               mm = mm + 1
               dd = jdy
            else
               exit
            end if
         end do
         if (mm > 12) then
            yyyy = yyyy + 1
            mm   = 1
         end if
      end subroutine s2d

      !====================================================
      ! Convert time, given in standard Gregorian notation,
      ! into Julian seconds.
      !
      ! Input:
      !    yyyy     year    (4 digits)
      !    mm       month   (2 digits)
      !    dd       day     (2 digits)
      !    hh       hour    (2 digits)
      !    mn       minutes (2 digits)
      !    ss       seconds (2 digits)
      !    yRef     yRef    (4 digits) - reference year
      !
      ! Output:
      !    d2s      cumulative seconds since  year 'yRef'
      !====================================================
      elemental real*8 function d2s(yyyy,mm,dd,hh,mn,ss, yRef)
         integer, intent(in) :: yyyy
         integer, intent(in) :: mm
         integer, intent(in) :: dd
         integer, intent(in) :: hh
         integer, intent(in) :: mn
         integer, intent(in) :: ss
         integer, intent(in) :: yRef

         integer :: monthDays(12)
         integer :: jyr, jdy
         integer :: i

         monthDays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

         jyr = 0
         if (yyyy > yRef) then
            do i = yRef, yyyy-1
               if (isLeapYear(i)) then
                  jyr = jyr + 366
               else
                  jyr = jyr + 365
               end if
            end do
         end if

         jdy = dd - 1

         if (isLeapYear(yyyy)) then
            monthDays(2) = 29
         else
            monthDays(2) = 28
         end if

         do i = 1, mm-1
            jdy = jdy + monthDays(i)
         end do

         d2s = jyr * 86400.d0                                           &
             + jdy * 86400.d0                                           &
             + hh  *  3600.d0                                           &
             + mn  *    60.d0                                           &
             + ss
      end function d2s

      !======================================================
      ! is this a leap year?
      !
      ! Input:
      !    year     year  (4 digits)
      !
      ! on output:
      !    .true.   if year is a leap year; .false. otherwise
      !======================================================
      elemental logical function isLeapYear(year)
         integer, intent(in) :: year
         if ( (mod(year,   4) == 0) .and.                               &
              (mod(year, 100) /= 0) .or.                                &
              (mod(year, 400) == 0) ) then
            isLeapYear = .true.
         else
            isLeapYear = .false.
         end if
      end function isLeapYear

      !===========================
      ! Run a suite of unit tests.
      !
      ! Input:
      !     none
      !
      ! Output:
      !     Unit tests are run.
      !===========================
      subroutine RUN_UNIT_TESTS()
         !--------
         ! Test #1
         !--------
         call nudge()
      end subroutine RUN_UNIT_TESTS

end module netNudge

!=============================
! Test program
! testNetNudge wrfinput_d01.nc
!=============================
program testNetNudge
   use netNudge
   implicit none
   call RUN_UNIT_TESTS()
!  call netNudge_mp_run_unit_tests()
end program testNetNudge
