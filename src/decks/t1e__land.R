! Trappist1e Landplanet deck, Hannah Woodward 2024
! Based on P1_THAI01e_Ben1.R M40 Desert World for THAI M. Way 2020/05
! - 4x5x40 layers modelE version (model top 0.1mb); 13 layers in the ocean
! - atmospheric composition: GHG.CO2_400ppm.txt
! - land: flat topography with roughness length 1cm
! - ocean: none
! - initial conditions:
!    - atmosphere: 300K isothermal temp, zero water vapor, 1000mb -> 0.1mb
! - radiation: planet (SOCRATES)
! Guide on creating a new planet:
! https://docs.google.com/document/d/1zrXQyEjXLRldWdyiZn2JxxQXz6dU8bayqOltf6eIJgo/edit

Preprocessor Options
!! ----- Update PLANET_PARAMS value to match planet config in
!! ----- $HOME/$MODELDIR/model/shared/PlanetParams_mod.F90
#define PLANET_PARAMS Trappist1e

!! --- Other options
#define USE_ENT
#define NEW_IO
#define PBL_USES_GCM_TENDENCIES        ! Not used in Mars model
#define USE_PLANET_RAD
#define GISS_RAD_OFF
!#define CACHED_SUBDD                  ! for sub-daily diagnostics (model output)
End Preprocessor Options

Object modules: (in order of decreasing priority)
       ! resolution-specific source codes
Atm72x46                               ! horizontal resolution is 72x46 -> 4x5deg
AtmL40p STRAT_DUM                      ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_M
FFT72                                  ! Fast Fourier Transform
IO_DRV                                 ! new i/o
!SUBDD                                 ! for sub-daily diagnostics (model output)

     ! GISS dynamics w/o gravity wave drag
ATMDYN MOMEN2ND                        ! atmospheric dynamics
QUS_DRV QUS3D                          ! advection of Q/tracers

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                                 ! model geometry
DIAG_ZONAL GCDIAGb                     ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                          ! diagn/post-processing output
MODEL_COM                              ! calendar, timing variables
MODELE_DRV                             ! ModelE cap
MODELE                                 ! initialization and main loop
ATM_COM                                ! main atmospheric variables
ATM_DRV                                ! driver for atmosphere-grid components
ATMDYN_COM                             ! atmospheric dynamics
ATM_UTILS                              ! utilities for some atmospheric quantities
QUS_COM QUSDEF                         ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM         ! clouds modules
SURFACE SURFACE_LANDICE FLUXES         ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM        ! land surface and soils + snow model
VEG_DRV                                ! vegetation
! VEG_COM VEGETATION                   ! old vegetation
ENT_DRV ENT_COM   ! + Ent              ! new vegetation
PBL_COM PBL_DRV PBL                    ! atmospheric pbl
ATURB_E1                               ! turbulence in whole atmosphere
LAKES_COM LAKES                        ! lake modules
SEAICE SEAICE_DRV                      ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV        ! land ice modules
ICEDYN_DRV ICEDYN                      ! ice dynamics modules
Zenith	   			                   ! shared zenith calculations
RAD_COM RAD_DRV RADIATION              ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo    ! radiation and albedo
DIAG_COM DIAG DEFACC                   ! diagnostics DEFACC_OLD is before David's fix to diagnostics
OCN_DRV                                ! driver for ocean-grid components
OCEAN OCNML                            ! using qflux ocean, but in reality not

! planet radiation source files
planet_rad planet_alb lw_control sw_control

SparseCommunicator_mod              ! sparse gather/scatter module

Components:
shared MPI_Support solvers giss_LSM
dd2d
socrates
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
!! ----- ATMOSPHERE + SURFACE INITIAL CONDITIONS (CAN LEAVE AS IS)
AIC=AIC_M40_1bar_300K_q0.nc	           ! 300K isothermal temp, zero water vapor, 1000mb
! comment this out to avoid grav wave drag.
!ZVAR=ZVAR4X5.nc                       ! topographic variation for gwdrag

!! ----- CUSTOM TOPOGRAPHY
TOPO=Z72X46N_flat_land.nc              ! flat surface
VEG=V72X46.1.cor2_no_crops03_alb35.nc                  ! provides surface albedo
SOILIC=soilic_wetness0.nc              ! wetness(1-6)=0, temp(1-6)=27C/300K, snow=0
SOIL=allsand.nc		                   ! 100% sand with sl(slope)=0
TOP_INDEX=planet/desert_world/stdev_72x46_desertworld.nc
ROUGHL=z0m_72x46_desertworld_1cm.nc    ! 1cm (0.01m) roughness length

!! ----- RADIATION (CAN LEAVE AS IS)
RADN1=sgpgxg.table8                    ! rad.tables and history files
RADN3=miescatpar.abcdv2
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
ISCCP=ISCCP.tautables
MSU_wts=MSU.RSS.weights.data
REG=REG4X5                             ! special regions-diag

!! ----- GHG GAS CONCENTRATIONS
GHG=GHG.CO2_400ppm.txt

Label and Namelist:
t1e__land (1880 atm.,the current modelE version)

DTFIX=300

!! ----- CHECK ATMOS COMPOSITION SCALING FACTORS
&&PARAMETERS
O2X=0.
NO2X=0.
N2OX=0.
CFC11X=0.
CFC12X=0.
N2CX=0.
XGHGX=0.
YGHGX=0.
SO2X=0.
O3X=0.
CH4X=1.
CO2X=1.

! turn off both horizontal (latitudinal) and vertical distributions of GHGs
l_uniform_ghg=1

!! ----- STELLAR CONFIG
solar_spec_dir='/home/app/socrates/stellar_spectra'
spectral_dir='/home/app/socrates/spectral_files'
!! --- find available rad files using `ls /home/app/socrates/spectral_files`
spectral_file_lw='sp_lw_dsa_ar10bar/sp_lw_15_dsa_ar10bar'
spectral_file_sw='sp_sw_dsa_ar10bar/sp_sw_43_dsa_ar10bar'
!! --- find available host using `ls /home/app/socrates/stellar_spectra`
solar_spec='trappist1'

!! ----- PLANET CONFIG
planet_s0=900.0
planetName='Trappist1e'
eccentricity=0.
obliquity=0.
siderealorbitalPeriod=527040.0d0  ! unit seconds; 6.1 days
siderealRotationPeriod=527040.0d0 ! unit seconds; 6.1 days
quantizeYearLength='False'

! calculate initial surface pressure for hydrostatic consistency with topography
initial_psurf_from_topo=1

MAXCTOP=1.
wsn_max=2.
minGroundTemperature=-250.0d0
maxGroundTemperature=300.0d0
vegCO2X_off=1                             ! ocn is prognostic
variable_lk=1                             ! Variable lakes on
mincolmass=1000.
maxcolmass=20000000000.
land_CO2_bc_flag=0

! linear damping timescales (sec) for winds in the top 3 layers.
rtau=320000.,270000.,220000.,170000.,120000.,70000.

! parameters set for coupled ocean runs:
KOCEAN=0            ! disable ocn

! Add SUBDD output
!SUBDD='olrcs lwcrf_toa incsw_toa srnf_toa'
!SUBDD1='tgrnd swds lwus lwds t u v w swhr lwhr'
!SUBDD2='swup_toa_clrsky'
!SUBDD3='p_surf z p_3d'
!NSUBDD=12
!nday_subdd=48
!write_daily_files=1
minCalendarDaysPerYear=1

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1         ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

xCDpbl=1.
!cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! Tuning parameters as of 2016/08/17 (D.S. Amundsen)
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=0.50   ! below 850mb and MC regions; then tune this to get rad.balance
radiusl_multiplier=0.97

! Cloud inhomogeneity correction
KCLDEP=1    ! use a constant value for CLDEPS
EPSCON=0.12 ! use CLDEPS=0.12

H2OstratX=1.

H2ObyCH4=0.         ! deactivates strat.H2O generated by CH4
KSIALB=0            ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=1850
crops_yr=-1         ! if -1, crops in VEG-file is used
volc_yr=-1
od_cdncx=0.         ! don't include 1st indirect effect
cc_cdncx=0.         ! don't include 2nd indirect effect (used 0.0036)
aer_rad_forc=0
cloud_rad_forc=1

! parameters that control the Shapiro filter
DT_XUfilter=225.    ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=225.    ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.      ! Shapiro filter on V in N-S direction
DT_YUfilter=0.      ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
DTsrc=1756.8
DT=225.
NIsurf=1            ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0      ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1    ! calls radiation twice; use =0 to save cpu time
KCOPY=1             ! saving acc + rsf
KRSF=1200           ! save one checkpoint restart file every century
isccp_diags=1       ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13            ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13            ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48             ! to get daily energy history use nda4=24*3600/DTsrc
Nssw=2              ! until diurnal diags are fixed, Nssw has to be even
Ndisk=480
nrad=2
ndigits_year=5
&&END_PARAMETERS

 &INPUTZ
 YEARI=00001,MONTHI=1,DATEI=1,HOURI=0!, IYEAR1=0001 ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=00501,MONTHE=1,DATEE=1,HOURE=0, KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=00002,MONTHE=1,DATEE=1,HOURE=0,IWRITE=1,JWRITE=1,
/
