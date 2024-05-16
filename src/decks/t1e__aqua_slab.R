! Trappist1e Aquaplanet slab deck, Hannah Woodward 2024
! Based on P1_THAI01e_Hab1.R M40 aqua planet for THAI. Thomas Fauchez 2018/10
! - 4x5x40 layers modelE version (model top 0.1mb); 13 layers in the ocean
! - atmospheric composition: GHG.CO2_400ppm.txt
! - ocean: q-flux
! - radiation: planet (SOCRATES)
! - initial conditions:
!    - atmosphere: 300K isothermal temp, zero water vapor, 1000mb -> 0.1mb
!    - ocean: 300K isothermal temp, 100m depth, 13 layers
! Additional modules required (included in podman image):
! - ATMDYN_HIGHLOWTEMP.f
! - SEAICE_DRVLOWTEMP.f
! Guide on creating a new planet:
! https://docs.google.com/document/d/1zrXQyEjXLRldWdyiZn2JxxQXz6dU8bayqOltf6eIJgo/edit

Preprocessor Options
!! ----- Update PLANET_PARAMS value to match planet config in
!! ----- $HOME/$MODELDIR/model/shared/PlanetParams_mod.F90
#define PLANET_PARAMS Trappist1e

!! --- Other options
!#define TRACERS_ON                    ! include tracers code
#define USE_ENT
#define NEW_IO
#define PBL_USES_GCM_TENDENCIES
#define USE_PLANET_RAD
#define GISS_RAD_OFF
!#define DEBUG_RADIATION_KEEP
!#define CACHED_SUBDD                  ! for sub-daily diagnostics (model output)
#define MORE_CLOUD_DIAGS_3D
End Preprocessor Options

Object modules: (in order of decreasing priority)
       ! resolution-specific source codes
Atm72x46                               ! horizontal resolution is 72x46 -> 4x5deg
AtmL40p STRAT_DUM                      ! vertical resolution is 40 layers -> 0.1mb
DIAG_RES_M
!ORES_5x4_L13                          ! ocean horiz res 4x5deg, 13 vert layers
FFT72                                  ! Fast Fourier Transform
IO_DRV                                 ! new i/o
!SUBDD                                 ! for sub-daily diagnostics (model output)

     ! GISS dynamics w/o gravity wave drag
ATMDYN_HIGHLOWTEMP MOMEN2ND            ! atmospheric dynamics
!QUS_DRV QUS3D                         ! advection of Q/tracers
QUS_DRV TQUS_DRV

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
SEAICE SEAICE_DRVLOWTEMP               ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV        ! land ice modules
ICEDYN_DRV ICEDYN                      ! ice dynamics modules
Zenith                                 ! shared zenith calculations
RAD_COM RAD_DRV RADIATION              ! radiation modules
RAD_UTILS ALBEDO READ_AERO ocalbedo    ! radiation and albedo
DIAG_COM DIAG DEFACC                   ! diagnostics
OCN_DRV                                ! driver for ocean-grid components
OCEAN OCNML                            ! qflux ocean modules

! planet radiation source files
planet_rad planet_alb lw_control sw_control

Components:
shared MPI_Support solvers giss_LSM
dd2d
socrates
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
!! ----- ATMOSPHERE, OCEAN, ICE, GLACIAL INITIAL CONDITIONS (CAN LEAVE AS IS)
AIC=AIC_M40_1bar_300K_q0.nc            ! 300K isothermal temp, zero water vapor, 1000mb pressure
GIC=GIC.E046D3M20A.1DEC1955.aqua.nc    ! initial conditions (ground) aquaplanet (ssi=0)
GLMELT=GLMELT_4X5.OCN.nc	           ! glacial melt distribution

! presc. climatological sea ice
!ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc

OHT=zero_OHT_4x5_100m.nc
OCNML=zero_OCNML_4x5_100m.nc

CDN=CD4X500S.ext.nc                    ! surf.drag coefficient
! comment this out to avoid grav wave drag.
!ZVAR=ZVAR4X5.nc                       ! topographic variation for gwdrag

!! ----- CUSTOM TOPOGRAPHY
TOPO=Z72X46N_gas.1_nocasp_aqua2.nc
OSTRAITS=OSTRAITS_72x46_aqua.nml
RVR=RD_modelE_M_btub004C.nc
NAMERVR=RD_modelE_M.names_btub0.txt
VEG=V72X46.1.cor2_no_cropsYANG.ext.nc  ! albedo=0.2 no veg
SOIL=S4X50093YANG.ext.nc               ! 50/50 clay/sand mix
TOP_INDEX=top_index_72x46_a.ij.ext.nc  ! only used if #define DO_TOPMODEL_RUNOFF
SOILCARB_global=soilcarb_top30cm_4x5.nc

!soil_textures=soil_textures_top30cm.nc

!! ----- RADIATION (CAN LEAVE AS IS)
RADN1=sgpgxg.table8                    ! rad.tables and history files
RADN3=miescatpar.abcdv2
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
ISCCP=ISCCP.tautables
MSU_wts=MSU.RSS.weights.data
REG=REG4X5                             ! special regions-diag

!! ----- GHG GAS CONCENTRATIONS
GHG=GHG.CO2_400ppm.txt

!! ----- UPDATE THE LABEL
Label and Namelist:
t1e__aqua_slab (1880 atm.,the current modelE version)


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
vegCO2X_off=1
variable_lk=1                          ! No dynamic lakes
mincolmass=1000.
maxcolmass=20000000000.

! linear damping timescales (sec) for winds in the top 6 layers
rtau=320000.,270000.,220000.,170000.,120000.,70000.

Kvflxo=0            ! usually set to 1 only during a prescr.ocn run by editing "I"
ocn_cycl=1          ! =0 if ocean varies from year to year

! parameters set for coupled ocean runs:
KOCEAN=1            ! ocn is prognostic

! Add SUBDD output
!SUBDD=' '          ! no sub-daily frequency diags
!SUBDD='olrcs lwcrf_toa incsw_toa swup_toa srnf_toa tsavg swds'
!SUBDD1='lwus lwds FOOPN t u v w swhr lwhr q rh'
!SUBDD2='totcld wtrcld icecld  icf wcf icmmr wcmmr icsiz wcsiz'
!SUBDD3='lwp iwp swup_toa_clrsky tgrnd p_surf z p_3d'
!NSUBDD=12          ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
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

!PTLISO=15.         ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.

! Tuning parameters as of 2016/10/05 (D.S. Amundsen)
U00a=.50            ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=0.5            ! below 850mb and MC regions; then tune this to get rad.balance
wmu_multiplier = 2.0

H2OstratX=1.

H2ObyCH4=0.         ! de-activates strat.H2O generated by CH4
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
KRSF=1200           ! save one checkpoint restart file every 100 years/orbits
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
!nmonav=1
!keep_params=1
!init_topog_related=0
!init_flake=0
&&END_PARAMETERS

 &INPUTZ
 YEARI=00001,MONTHI=1,DATEI=1,HOURI=0!, IYEAR1=0001 ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=00501,MONTHE=1,DATEE=1,HOURE=0, KDIAG=13*0,
 ISTART=2,IRANDI=0, YEARE=00002,MONTHE=1,DATEE=1,HOURE=0,IWRITE=1,JWRITE=1,
/
