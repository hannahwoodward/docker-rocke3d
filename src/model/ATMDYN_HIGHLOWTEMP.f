#include "rundeck_opts.h"

      module ATMDYN
      use geom, only : imh,fim,byim
      implicit none
      private

      public init_ATMDYN, DYNAM
     &     ,FILTER
     &     ,COMPUTE_DYNAM_AIJ_DIAGNOSTICS
     &     ,AFLUX, COMPUTE_MASS_FLUX_DIAGS
     &     ,SDRAG,shap1,shap1d
     &     ,fcuva,fcuvb

C**** Variables used in DIAG5 calculations
!@var FCUVA,FCUVB fourier coefficients for velocities
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: FCUVA,FCUVB

!@var aflux_topo_adjustments whether to adjust uphill air mass fluxes
!@+   around steep orography
      logical :: aflux_topo_adjustments=.true.

!@var pfilter_using_slp whether to zonally filter surface pressure using
!@+   diagnosed SLP.  If false, the filter adjusts surface pressure
!@+   to smooth the zonal PGF
      logical :: pfilter_using_slp=.true.

      contains

      SUBROUTINE init_ATMDYN
      USE DOMAIN_DECOMP_ATM, only: grid
      use domain_decomp_1d, only : am_i_root
      use resolution, only : lm,mfixs,mdrya
      use model_com, only : dtsrc
      use constant, only : planet_name
      use dynamics
      use Dictionary_mod
      implicit none
      call get_param( "DT", DT )
C**** NIdyn=dtsrc/dt(dyn) has to be a multiple of 2
      NIdyn = 2*nint(.5*dtsrc/dt)
      if (is_set_param("DT") .and. nint(DTsrc/dt).ne.NIdyn) then
        if (AM_I_ROOT())
     *       write(6,*) 'DT=',DT,' has to be changed to',DTsrc/NIdyn
        if(trim(planet_name).eq.'Earth') then
          ! We do not allow automatic rounding of timesteps for Earth,
          ! since such a situation is usually due to a rundeck typo
          call stop_model('INPUT: DT inappropriately set',255)
        endif
      end if
      DT = DTsrc/NIdyn
      call set_param( "DT", DT, 'o' )         ! copy DT into DB

c      call setDtParam('dt', dt, dtSrcUsed)

      ! Defaults for bounds on surface pressure (i.e. column mass)
      ! minimum: a multiple of the mass in the constant-pressure domain
      !          (leaving enough for the variable-mass layers)
      ! maximum: a multiple of global mean column mass
      ! the following ratios were taken from the older version of ADVECM
      mincolmass = mfixs*(350d0/150d0)
      maxcolmass = mdrya*(1160d0/984d0)
      if(is_set_param('mincolmass'))
     &     call get_param('mincolmass',mincolmass)
      if(is_set_param('maxcolmass'))
     &     call get_param('maxcolmass',maxcolmass)

      if(trim(planet_name).ne.'Earth') then
        aflux_topo_adjustments = .false.
        ! SLP diagnosis currently uses Earth-specific lapse rates
        pfilter_using_slp = .false.
      endif

      call sync_param( "NFILTR", NFILTR ) !!
      call sync_param( "DT_XVfilter", DT_XVfilter )
      call sync_param( "DT_XUfilter", DT_XUfilter )
      call sync_param( "DT_YVfilter", DT_YVfilter )
      call sync_param( "DT_YUfilter", DT_YUfilter )
      call sync_param( "MFILTR", MFILTR )
      call sync_param( "do_polefix", do_polefix )
C**** Determine if FLTRUV is called.
      QUVfilter = .false.
      if (DT_XUfilter>0. .or. DT_XVfilter>0. .or.
     *    DT_YUfilter>0. .or. DT_YVfilter>0.)  QUVfilter = .true.
      if (QUVfilter) then
         if (DT_XUfilter > 0. .and. DT_XUfilter < DT) then
             DT_XUfilter = DT
             WRITE(6,*) "DT_XUfilter too small; reset to :",DT_XUfilter
             call set_param( "DT_XUfilter", DT_XUfilter, 'o' )
         end if
         if (DT_XVfilter > 0. .and. DT_XVfilter < DT) then
             DT_XVfilter = DT
             WRITE(6,*) "DT_XVfilter too small; reset to :",DT_XVfilter
             call set_param( "DT_XVfilter", DT_XVfilter, 'o' )
         end if
         if (DT_YUfilter > 0. .and. DT_YUfilter < DT) then
             DT_YUfilter = DT
             WRITE(6,*) "DT_YUfilter too small; reset to :",DT_YUfilter
             call set_param( "DT_YUfilter", DT_YUfilter, 'o' )
         end if
         if (DT_YVfilter > 0. .and. DT_YVfilter < DT) then
             DT_YVfilter = DT
             WRITE(6,*) "DT_YVfilter too small; reset to :",DT_YVfilter
             call set_param( "DT_YVfilter", DT_YVfilter, 'o' )
         end if
      end if
c Warn if polar fixes requested for a model not having a half polar box
c     if(do_polefix.eq.1 .and. jm.ne.46) then
c        do_polefix = 0
c        write(6,*) 'Polar fixes are currently applicable only to'//
c    &           'models having a half polar box; no fixes applied'
c     endif

      CALL AVRX
      ALLOCATE( FCUVA(0:IMH, grid%j_strt_halo:grid%j_stop_halo, LM, 2),
     &          FCUVB(0:IMH, grid%j_strt_halo:grid%j_stop_halo, LM, 2))
      end SUBROUTINE init_ATMDYN

c      subroutine setDtParam(tName, tParam, dtSrc)
c      USE Dictionary_mod
c      use BaseTime_mod
c      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
c      character(len=*) :: tName
c      type (BaseTime), intent(in) :: dtSrc
c      real(8), intent(inout) :: tParam
c      real(8) :: tOld
c
c      tOld = tParam
c      call get_param(tName, tParam)
c      tParam = dtSrc%convertToReal()/nint(dtSrc%convertToReal()/tParam)
c      call set_param( tName, tParam, 'o' )
c
c      if (abs(tParam-tOld) .gt. 1.0e-15) then
c        if (AM_I_ROOT()) then
c          write(6,*) trim(tName),' has changed from ', tOld,' to ',
c     *      tParam
c        end if
c      end if
c      end subroutine setDtParam


      SUBROUTINE DYNAM
!@sum  DYNAM Integrate dynamic terms
!@vers 2013/10/31
!@auth Original development team
      Use CONSTANT,   Only: by3,byGRAV,RGAS,SHA,kg2mb,UNDEF_VAL
      Use RESOLUTION, Only: IM,JM,LM,LS1, MFIXs
      USE MODEL_COM, only : DTsrc
      Use ATM_COM,    Only: MA,U,V,T,Q,QCL,QCI,MASUM, MUs,MVs,MWs, GZ, P
      Use GEOM,       Only: AXYP
      USE SOMTQ_COM, only : tmom,mz
      Use DYNAMICS,   Only: MU,MV,MW, pu,pv,sd,dut,dvt
     &    ,cos_limit,nidyn,dt,mrch,nstep,quvfilter,USE_UNR_DRAG
      Use DIAG_COM,   Only: MODD5K,NDA5K,NDAA
      USE DOMAIN_DECOMP_ATM, only: grid
      Use DOMAIN_DECOMP_1D,  Only: GetDomainBounds, GLOBALSUM, SOUTH,
     *                             HALO_UPDATE, HALO_UPDATE_COLUMN
      USE MOMENTS, only : advecv
      IMPLICIT NONE

      Real*8,Dimension (LM, IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &   MEVEN,MODD1,MODD3
      Real*8,Dimension (IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &   MSUMODD, PA,PB,PC, FPEU,FPEV, AM1,AM2
      Real*8,Dimension (IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &   MMA,TZ,        !  even leap frog arrays
     &   UT,VT,TT,TZT,  !  odd leap frog arrays
     &   UX,VX,         !  initial forward step arrays
     &   PIJL, UNRDRAG_x,UNRDRAG_y

      REAL*8 DTFS,DTLF, DAMSUM
      INTEGER I,J,L,IP1,IM1   !@var I,J,L,IP1,IM1  loop variables
      INTEGER NS, NSOLD,MODDA    !? ,NIdynO

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      MSUMODD(:,[grid%J_STRT_HALO,grid%J_STOP_HALO])=UNDEF_VAL

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

!?    NIdynO=MOD(NIdyn,2)   ! NIdyn odd is currently not an option
      DTFS=DT*2./3.
      DTLF=2.*DT
      NS = NIDYN  ;  NSOLD = NIDYN
      MUs(:,:,:) = 0  ;  MVs(:,:,:) = 0  ;  MWs(:,:,:) = 0

C**** Leap-frog re-initialization: IF (NS.LT.NIdyn)
  300 CONTINUE
      Do J=J_0,J_1  ;  Do I=1,IM
         MASUM(I,J) = Sum (MA(:,I,J))  ;  EndDo  ;  EndDo
      Call HALO_UPDATE_COLUMN (GRID, MA,    From=SOUTH)
      Call HALO_UPDATE        (GRID, MASUM, From=SOUTH)
      UX(:,:,:) = U(:,:,:)  ;  UT(:,:,:) = U(:,:,:)
      VX(:,:,:) = V(:,:,:)  ;  VT(:,:,:) = V(:,:,:)
      TZ(:,:,:) = TMOM(MZ,:,:,:)

!**** Initial forward step:  MODD3 = MA + .667*DT*F(U,V,MA)
      MRCH=0
#          ifdef NUDGE_ON
           Call NUDGE_PREP
           Call NUDGE (UX,VX,DTFS)
#          endif
      Call AFLUX  (NS,   U,V,MA,MASUM,    MA   ,MASUM)
      Call ADVECM (DTFS,         MA,      MODD3,MSUMODD)
      Call GWDRAG (DTFS, U,V,       UX,VX,MODD3, T,TZ, .True.)
      Call VDIFF  (DTFS, U,V,       UX,VX,MODD3, T)
      Call ADVECV (DTFS, U,V,MA, MA,UX,VX,MODD3)
      Call PGF    (DTFS, U,V,MA,    UX,VX,MODD3, T,TZ)
       CALL CALC_PIJL (LM,P,PIJL)
       PU(:,:,:) = MU(:,:,:)*kg2mb
       PV(:,:,:) = MV(:,:,:)*kg2mb
       SD(:,:,:) = MW(:,:,:)*kg2mb
       PB(:,:)   = (MSUMODD(:,:) - MFIXs)*kg2mb
      call isotropuv(ux,vx,COS_LIMIT)

!**** Initial backward step:  MODD1 = MA + DT*F(UX,VX,MODD3)
      MRCH=-1
#          ifdef NUDGE_ON
           Call NUDGE (UT,VT,DT)
#          endif
      Call AFLUX  (NS, UX,VX,MODD3,MSUMODD,  MA,   MASUM)
      Call ADVECM (DT,              MA,      MODD1,MSUMODD)
      Call GWDRAG (DT, UX,VX,          UT,VT,MODD1, T,TZ, .False.)
      Call VDIFF  (DT, UX,VX,          UT,VT,MODD1, T)
      Call ADVECV (DT, UX,VX,MODD3, MA,UT,VT,MODD1)
      Call PGF    (DT, UX,VX,MODD3,    UT,VT,MODD1, T,TZ)
       CALL CALC_PIJL (LS1-1,PB,PIJL)
       PU(:,:,:) = MU(:,:,:)*kg2mb
       PV(:,:,:) = MV(:,:,:)*kg2mb
       SD(:,:,:) = MW(:,:,:)*kg2mb
       PA(:,:)   = (MSUMODD(:,:) - MFIXs)*kg2mb
      call isotropuv(ut,vt,COS_LIMIT)
      GO TO 360

!**** Odd leap frog step:  MODD3 = MODD1 + 2*DT*F(U,V,MA)
  340 MRCH=-2
#          ifdef NUDGE_ON
           Call NUDGE (UT,VT,DTLF)
#          endif
      Call AFLUX  (NS,   U,V,MA,MASUM,       MODD1,MSUMODD)
      Call ADVECM (DTLF,         MODD1,      MODD3,MSUMODD)
      Call GWDRAG (DTLF, U,V,          UT,VT,MODD3, T,TZ, .False.)
      Call VDIFF  (DTLF, U,V,          UT,VT,MODD3, T)
      Call ADVECV (DTLF, U,V,MA, MODD1,UT,VT,MODD3)
      Call PGF    (DTLF, U,V,MA,       UT,VT,MODD3, T,TZ)
       CALL CALC_PIJL (LS1-1,P,PIJL)
       PU(:,:,:) = MU(:,:,:)*kg2mb
       PV(:,:,:) = MV(:,:,:)*kg2mb
       SD(:,:,:) = MW(:,:,:)*kg2mb
       PB(:,:)   = (MSUMODD(:,:) - MFIXs)*kg2mb
      call isotropuv(ut,vt,COS_LIMIT)
      PA(:,:) = PB(:,:)     ! LOAD PB TO PA
      MODD1(:,:,:) = MODD3(:,:,:)
      NS = NS - 1

!**** Even leap frog step:  MA = MEVEN + 2*DT*F(UT,VT,MODD1)
  360      MODD5K = Mod (NSTEP+4-NS + NDA5K*NIDYN, NDA5K*NIDYN+2)
      MRCH=2
      MEVEN(:,:,:) = MA(:,:,:)
       PC(:,:) = P(:,:)      ! LOAD P TO PC
#          ifdef NUDGE_ON
           Call NUDGE (U,V,DTLF)
#          endif
      Call AFLUX  (NS,   UT,VT,MODD1,MSUMODD,   MEVEN,MASUM)
      Call ADVECM (DTLF,              MEVEN,    MA,MASUM)
      Call GWDRAG (DTLF, UT,VT,             U,V,MA, T,TZ, .False.)
      Call VDIFF  (DTLF, UT,VT,             U,V,MA, T)
      Call ADVECV (DTLF, UT,VT,MODD1, MEVEN,U,V,MA)
       CALL CALC_PIJL (LS1-1,PA,PIJL)
       PU(:,:,:) = MU(:,:,:)*kg2mb
       PV(:,:,:) = MV(:,:,:)*kg2mb
       SD(:,:,:) = MW(:,:,:)*kg2mb
        P(:,:)   = (MASUM(:,:) - MFIXs)*kg2mb
            MODDA = Mod (NSTEP+4-NS + NDAA*NIDYN, NDAA*NIDYN+2)  ! strat
         IF(MODDA.LT.MRCH) CALL DIAGA0   ! strat
C**** ACCUMULATE MASS FLUXES FOR TRACERS and Q
      MUs(:,:,:) = MUs(:,:,:) + PU(:,:,:)
      MVs(:,:,:) = MVs(:,:,:) + PV(:,:,:)
      MWs(:,:,1:LM-1) = MWs(:,:,1:LM-1) + SD(:,:,:)
C**** ADVECT Q AND T
       TT(:,:,:) =  T(:,:,:)
      TZT(:,:,:) = TZ(:,:,:)
      Do L=1,LM
         MMA(:,:,L) = MEVEN(L,:,:)*AXYP(:,:)  ;  EndDo
      Call AADVT (DTLF, MMA,T,TMOM, .False., FPEU,FPEV)
!     save z-moment of temperature in contiguous memory for later
      TZ(:,:,:) = TMOM(MZ,:,:,:)
       PC(:,:)   = .5*( P(:,:)  + PC(:,:))
       TT(:,:,:) = .5*( T(:,:,:)+ TT(:,:,:))
      TZT(:,:,:) = .5*(TZ(:,:,:)+TZT(:,:,:))

      CALL CALC_PIJL(LS1-1,PC,PIJL)
c      CALL CALC_PIJL(LS1-1,PA,PIJL) ! true leapfrog
      Call PGF    (DTLF, UT,VT,MODD1,       U,V,MA, TT,TZT)

      Call COMPUTE_MASS_FLUX_DIAGS (GZ, MU,MV, DT)

      CALL CALC_AMPK(LS1-1)
      call isotropuv(u,v,COS_LIMIT)
      if (USE_UNR_DRAG==0) CALL SDRAG (DTLF)
         If (Mod(NSTEP+4-NS+NDAA*NIDYN,NDAA*NIDYN+2) < MRCH)  Then
           CALL DIAGA
           CALL DIAGB
           CALL EPFLUX (U,V,T,P)
         ENDIF
C**** Restart after 8 steps due to divergence of solutions
      NS = NS - 1
      If (NSOLD-NS < 8 .and. NS > 1)  GoTo 340
      NSOLD=NS
      If (NS > 1)  GoTo 300

      if (USE_UNR_DRAG==1) then
         Call UNRDRAG (P,U,V,T,TZ,UNRDRAG_x,UNRDRAG_y)
         U(:,:,:) = U(:,:,:) + UNRDRAG_x(:,:,:) * DTsrc
         V(:,:,:) = V(:,:,:) + UNRDRAG_y(:,:,:) * DTsrc  ;  EndIf

!**** Convert summed mass fluxes from (kg/s) to (kg)
      MUs(:,:,:) = MUs(:,:,:) * DTLF
      MVs(:,:,:) = MVs(:,:,:) * DTLF
      MWs(:,:,1:LM-1) = MWs(:,:,1:LM-1) * DTLF

c apply east-west filter to U and V once per physics timestep
      Call FLTRUV
c apply north-south filter to U and V once per physics timestep
      call conserv_amb_ext(u,am1) ! calculate ang. mom. before filter
      call fltry2(u,1d0) ! 2nd arg could be set using DT_YUfilter
      call fltry2(v,1d0) ! 2nd arg could be set using DT_YVfilter
      call conserv_amb_ext(u,am2) ! calculate ang. mom. after filter
      am2(:,j_0stg:j_1stg) = am1(:,j_0stg:j_1stg)-am2(:,j_0stg:j_1stg)
      if(have_south_pole) am2(:,1) = 0.
      call globalsum(grid,am2,damsum,all=.true.)
      call add_am_as_solidbody_rotation(u,damsum) ! maintain global ang. mom.

      RETURN
      END SUBROUTINE DYNAM


      Subroutine COMPUTE_MASS_FLUX_DIAGS (GZ, MU,MV, DT)
      use RESOLUTION, only: IM, LM
      USE DOMAIN_DECOMP_ATM, only: grid
      use DOMAIN_DECOMP_1D, only: halo_update, SOUTH, getDomainBounds
      use DIAG_COM, only: AIJ => AIJ_loc, IJ_FGZU, IJ_FGZV

      Real*8,Intent(In) :: GZ(:,GRID%J_STRT_HALO:,:),
     *                     MU(:,GRID%J_STRT_HALO:,:),
     *                     MV(:,GRID%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      integer :: J_0S, J_1S
      integer :: J_0STG, J_1STG
      integer :: I, IP1, L, J

      call getDomainBounds(grid, J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG,
     &                           J_STRT_SKP =J_0S,  J_STOP_SKP =J_1S)

!     Call HALO_UPDATE (GRID, GZ, From=SOUTH)   haloed in PGF
      DO J=J_0S,J_1S ! eastward transports
      DO L=1,LM
         I=IM
         DO IP1=1,IM
            AIJ(I,J,IJ_FGZU)=AIJ(I,J,IJ_FGZU)+
     &           (GZ(I,J,L)+GZ(IP1,J,L))*MU(I,J,L)*DT ! use DT=DTLF/2
            I=IP1
         END DO
      END DO
      END DO
      DO J=J_0STG,J_1STG ! northward transports
      DO L=1,LM
         DO I=1,IM
            AIJ(I,J,IJ_FGZV)=AIJ(I,J,IJ_FGZV)+
     &           (GZ(I,J-1,L)+GZ(I,J,L))*MV(I,J,L)*DT ! use DT=DTLF/2
         END DO
      END DO
      END DO

      end subroutine compute_mass_flux_diags


      Subroutine COMPUTE_DYNAM_AIJ_DIAGNOSTICS (MUs, MVs, DT)
      use CONSTANT,      only: BY3
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      use DIAG_COM, only: AIJ => AIJ_loc,
     &     IJ_FGZU, IJ_FGZV, IJ_FMV, IJ_FMU
      use RESOLUTION, only: IM,JM,LM

      real*8, intent(in) :: MUs(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: MVs(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      integer :: I, IP1, J, L
      integer :: J_0STG, J_1STG, J_0S, J_1S
      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
      real*8 :: dtlf

      dtlf = 2.*dt

      call getDomainBounds(grid, J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG,
     &                           J_STRT_SKP =J_0S,  J_STOP_SKP =J_1S,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

      do j=J_0STG,J_1STG
      do L=1,LM
         AIJ(:,J,IJ_FMV)  = AIJ(:,J,IJ_FMV )+MVs(:,J,L)*DTLF
      enddo
      enddo

      do j=J_0S,J_1S
      do l=1,lm
         AIJ(:,J,IJ_FMU) = AIJ(:,J,IJ_FMU)+MUs(:,J,L)*DTLF
      enddo
      enddo

      if (haveLatitude(grid, J=1)) then
         do l=1,lm
            AIJ(:,1,IJ_FMU)  = AIJ(:, 1,IJ_FMU )+MUs(:, 1,L)*DTLF*BY3
         enddo
      endif
      if(haveLatitude(grid, J=JM)) then
         do l=1,lm
            AIJ(:,JM,IJ_FMU) = AIJ(:,JM,IJ_FMU )+MUs(:,JM,L)*DTLF*BY3
         enddo
      endif
      end subroutine COMPUTE_DYNAM_AIJ_DIAGNOSTICS


      Subroutine AFLUX (NS, U,V,MA,MASUM, ME,MESUM)
!@sum  AFLUX Calculates horizontal/vertical air mass fluxes
!@+    Input:  NS = decrementing leap-frog time step counter
!@+            MA,MASUM (kg/m^2), U,V (m/s) = step center values
!@+            ME,MESUM (kg/m^2) = mass at beginning of time step
!@+    Output: MU,MV (kg/s) = horizontal mass fluxes
!@+            MW (kg/s) = MW(L-1) + (CONV-CONVs + MM-MMs)*MFRAC/MFRACs
!@+            CONV (kg/s) =  horizontal mass convergence
!@+            SPA (.5 m/s) = filtered (U+U) defined on eastern cell edge
!!@auth Original development team
      Use RESOLUTION, Only: IM,JM,LM,LS1, MFIX,MFIXs,MFRAC
      Use ATM_COM,    Only: ZATMO
      Use GEOM,       Only: IMAXJ, DXYP,DYP,DXV, POLWT
      Use DYNAMICS,   Only: DT, MU,MV,MW,CONV, SPA,DO_POLEFIX
      Use DOMAIN_DECOMP_ATM, Only: GRID
      Use DOMAIN_DECOMP_1D,  Only: GetDomainBounds, SOUTH,NORTH,
     &                             HALO_UPDATE
      IMPLICIT NONE
      Integer,Intent(In) :: NS
      Real*8,Intent(InOut):: U(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *                       V(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *                   MA(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *                   MASUM(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      Real*8,Intent(In)::ME(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *                   MESUM(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)

!**** Local variables
      Real*8,Parameter :: TWOby3 = 2/3d0
      Integer :: I,J,L,Ip1,Im1, J1,J1XP,J1H,J1V, JN,JNXP,JNH,JNV
      Logical :: QSP,QNP
      Real*8  :: DUMMYS(IM),MUS,MVS,PBS,MVSA(LM), zNSxDT,
     *           DUMMYN(IM),MUN,MVN,PBN,MVNA(LM),
     *           USV0(IM,2,LM),VSV0(IM,2,LM), M,CONVs,MVARs

       zNSxDT = 1 / (NS*DT)
!****                             +---------+
!**** COMPUTATION OF MASS FLUXES  |  MA,T  MU    PRIMARY GRID ROW (J)
!**** ARAKAWA'S SCHEME B          +---MV---U,V   VELOCITY GRID ROW (J)
!****
      Call GetDomainBounds (GRID,
     &   J_STRT      = J1  , J_STOP      = JN  , !  primary row limits
     &   J_STRT_SKP  = J1XP, J_STOP_SKP  = JNXP, !  primary row limits excluding poles
     &   J_STRT_HALO = J1H , J_STOP_HALO = JNH , !  halo primary row limits
     &   J_STRT_STGR = J1V , J_STOP_STGR = JNV , !  velocity row limits
     &   HAVE_SOUTH_POLE = QSP, HAVE_NORTH_POLE = QNP)

      CALL HALO_UPDATE(grid, U,FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid, V,FROM=NORTH+SOUTH)

!**** Use interpolated velocity values at poles
      If (J1V==2)  Then
         USV0(:,1,:) = U(:,2,:)
         VSV0(:,1,:) = V(:,2,:)
         U(:,2,:) = POLWT*U(:,2,:) + (1-POLWT)*U(:,3,:)
         V(:,2,:) = POLWT*V(:,2,:) + (1-POLWT)*V(:,3,:)  ;  EndIf
      If (JNV==JM)  Then
         USV0(:,2,:) = U(:,JM,:)
         VSV0(:,2,:) = V(:,JM,:)
         U(:,JM,:) = POLWT*U(:,JM,:) + (1-POLWT)*U(:,JM-1,:)
         V(:,JM,:) = POLWT*V(:,JM,:) + (1-POLWT)*V(:,JM-1,:)  ;  EndIf

C**** Compute MU (kg/s) = eastward mass flux at non-polar points
      Do 2166 L=1,LM
      Do 2154 J=J1XP,JNXP
      DO 2154 I=1,IM
 2154 SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
      CALL AVRX (SPA(1,J1H,L))
      I=IM
      Do 2166 J=J1XP,JNXP
      DO 2165 IP1=1,IM
      MU(I,J,L) = .25*DYP(J)*SPA(I,J,L)*(MA(L,I,J)+MA(L,Ip1,J))
 2165 I=IP1
 2166 CONTINUE
C**** Compute MV (kg/s) = northward mass flux
      Do 2172 L=1,LM
      IM1=IM
      DO 2172 J=J1V,JNV
      DO 2170 I=1,IM
      MV(I,J,L) = .25*DXV(J)*(V(I,J,L)+V(Im1,J,L))*
     *            (MA(L,I,J)+MA(L,I,J-1))
 2170 IM1=I
 2172 CONTINUE

!**** Restore uninterpolated values of U,V at poles
      If (J1V==2)  Then
         U(:,2,:) = USV0(:,1,:)
         V(:,2,:) = VSV0(:,1,:)  ;  EndIf
      If (JNV==JM)  Then
         U(:,JM,:) = USV0(:,2,:)
         V(:,JM,:) = VSV0(:,2,:)  ;  EndIf

C**** Compute MU*3 at poles
      If (QSP)  Then
         Do L=1,LM
            MUS = Sum( U(:,2,L))*byIM*.25*DYP(2)*MA(L,1,1)
            MVS = Sum(MV(:,2,L))*byIM
            MVSA(L) = MVS
            DUMMYS(1) = 0
            Do I=2,IM
               DUMMYS(I) = DUMMYS(I-1) + (MV(I,2,L)-MVS)  ;  EndDo
            PBS = Sum(DUMMYS(:))*byIM
            SPA(:,1,L) = 4*(PBS-DUMMYS(:)+MUS) / (DYP(2)*MA(L,1,1))
             MU(:,1,L) = 3*(PBS-DUMMYS(:)+MUS)  ;  EndDo  ;  EndIf
      If (QNP)  Then
         Do L=1,LM
            MUN = Sum( U(:,JM,L))*byIM*.25*DYP(JM-1)*MA(L,1,JM)
            MVN = Sum(MV(:,JM,L))*byIM
            MVNA(L) = MVN
            DUMMYN(1) = 0
            Do I=2,IM
               DUMMYN(I) = DUMMYN(I-1) + (MV(I,JM,L)-MVN)  ;  EndDo
            PBN = Sum(DUMMYN(:))*byIM
            SPA(:,JM,L) = 4*(DUMMYN(:)-PBN+MUN) / (DYP(JM-1)*MA(L,1,JM))
             MU(:,JM,L) = 3*(DUMMYN(:)-PBN+MUN)  ;  EndDo  ;  EndIf

      if(do_polefix.eq.1) then
c To maintain consistency with subroutine ADVECV,
c adjust pu at the pole if no corner fluxes are used there
c in ADVECV.
         If (QSP)  MU(:,1 ,:) = MU(:,1 ,:)*TWOby3
         If (QNP)  MU(:,JM,:) = MU(:,JM,:)*TWOby3  ;  EndIf

      if(aflux_topo_adjustments) then

!**** Modify eastward uphill air mass fluxes around steep topography
      Do 310 J=J1XP,JNXP
      I = IM
      Do 310 Ip1=1,IM
      If (ZATMO(I,J) == ZATMO(Ip1,J))  GoTo 310
      If (ZATMO(I,J) <  ZATMO(Ip1,J))
     *   Then  ;  M = MASUM(I,J)
                  Do L=1,LS1-1
                     If (MU(I,J,L) <= 0)  GoTo 310
                     M = M - MA(L,I,J)
                     If (M <= MASUM(Ip1,J))  GoTo 310
                     MU(I,J,L+1) = MU(I,J,L+1) + MU(I,J,L)
                     MU(I,J,L) = 0  ;  EndDo
         Else  ;  M = MASUM(Ip1,J)
                  Do L=1,LS1-1
                     If (MU(I,J,L) >= 0)  GoTo 310
                     M = M - MA(L,Ip1,J)
                     If (M <= MASUM(I,J))  GoTo 310
                     MU(I,J,L+1) = MU(I,J,L+1) + MU(I,J,L)
                     MU(I,J,L) = 0  ;  EndDo  ;  EndIf
  310 I = Ip1

!**** Modify northward uphill air mass fluxes around steep topography
!**** Exceptional J loop boundaries, not:  Do J=J1V,JNV
      Do 320 J=Max(J1XP,3),JNXP
      Do 320 I=1,IM
      If (ZATMO(I,J-1) == ZATMO(I,J))  GoTo 320
      If (ZATMO(I,J-1) <  ZATMO(I,J))
     *   Then  ;  M = MASUM(I,J-1)
                  Do L=1,LS1-1
                     If (MV(I,J,L) <= 0)  GoTo 320
                     M = M - MA(L,I,J-1)
                     If (M <= MASUM(I,J))  GoTo 320
                     MV(I,J,L+1) = MV(I,J,L+1) + MV(I,J,L)
                     MV(I,J,L) = 0  ;  EndDo
         Else  ;  M = MASUM(I,J)
                  Do L=1,LS1-1
                     If (MV(I,J,L) >= 0)  GoTo 320
                     M = M - MA(L,I,J)
                     If (M <= MASUM(I,J-1))  GoTo 320
                     MV(I,J,L+1) = MV(I,J,L+1) + MV(I,J,L)
                     MV(I,J,L) = 0  ;  EndDo  ;  EndIf
  320 Continue

      endif ! aflux_topo_adjustments

!**** Compute CONV (kg/s) = horizontal mass convergence
      Call HALO_UPDATE (GRID,MU,From=SOUTH+NORTH) ! full halos needed later
      Call HALO_UPDATE (GRID,MV,From=SOUTH+NORTH)
      DO 2400 L=1,LM
      Do 1510 J=J1XP,JNXP
      IM1=IM
      DO 1510 I=1,IM
      CONV(I,J,L) = MU(IM1,J,L)-MU(I,J,L) + MV(I,J,L)-MV(I,J+1,L)
 1510 IM1=I
      If (QSP)  CONV(1,1,L) = -MVSA(L)
      If (QNP)  CONV(1,JM,L) = MVNA(L)
 2400 CONTINUE

C**** Compute MW (kg/s) = downward vertical mass flux
      DO 2435 J=J1,JN
      Do 2435 I=1,IMAXJ(J)
         CONVs = Sum(CONV(I,J,:))
         MVARs = MESUM(I,J) - MFIXs
         MW(I,J,LM-1) = CONV(I,J,LM) - CONVs*MFRAC(LM) +
     +      (ME(LM,I,J) - MFIX(LM) - MVARs*MFRAC(LM))*DXYP(J)*zNSxDT
         Do 2430 L=LM-2,1,-1
 2430       MW(I,J,L) = CONV(I,J,L+1) - CONVs*MFRAC(L+1) + MW(I,J,L+1) +
     +       (ME(L+1,I,J) - MFIX(L+1) - MVARs*MFRAC(L+1))*DXYP(J)*zNSxDT
 2435 CONTINUE
      Do L=1,LM-1
         If (QSP)  MW(2:IM,1 ,L) = MW(1,1 ,L)
         If (QNP)  MW(2:IM,JM,L) = MW(1,JM,L)  ;  EndDo

      RETURN
      EndSubroutine AFLUX


      SUBROUTINE ADVECM (DT1, MOLD,MNEW,MSUM)
!@sum  ADVECM Calculates updated column pressures using mass fluxes
!@auth Original development team
      Use RESOLUTION, Only: IM,JM,LM, MTOP
      Use ATM_COM,    Only: ZATMO,U,V,T,Q
      Use GEOM,       Only: IMAXJ,byDXYP
      Use DYNAMICS,   Only: MRCH,MW,CONV
      use dynamics, only : mincolmass,maxcolmass
      USE DOMAIN_DECOMP_ATM, only: grid
      Use DOMAIN_DECOMP_1D,  Only: GetDomainBounds, HALO_UPDATE_COLUMN,
     *                             HALO_UPDATE, GLOBALMAX, NORTH,SOUTH
      IMPLICIT NONE
      Real*8,Intent(In) :: DT1,
     *                   MOLD(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      Real*8,Intent(Out) :: MSUM(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *                   MNEW(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)

!**** Local variables
      Integer :: I,J,L,Im1, J1,JN
      INTEGER :: n_exception, n_exception_all
      Logical :: QSP,QNP
      Real*8  :: MNEWs

      Call GetDomainBounds (GRID,
     &   J_STRT = J1, J_STOP = JN, !  primary row limits
     &   HAVE_SOUTH_POLE = QSP, HAVE_NORTH_POLE = QNP)

!**** Compute MNEW and MSUM, the new atmospheric mass distribution
      ! 1st pass count warning/termination events
      ! This avoids the need for 2 halo fills during normal
      ! execution.
      n_exception = 0
      Do J=J1,JN
      Do I=1,IMAXJ(J)
         MNEW(1,I,J) = MOLD(1,I,J) +
     +      DT1*(CONV(I,J,1) + MW(I,J,1))*byDXYP(J)
         Do L=2,LM-1
            MNEW(L,I,J) = MOLD(L,I,J) +
     +         DT1*(CONV(I,J,L) + MW(I,J,L) - MW(I,J,L-1))*byDXYP(J)
            EndDo
         MNEW(LM,I,J) = MOLD(LM,I,J) +
     +         DT1*(CONV(I,J,LM) - MW(I,J,LM-1))*byDXYP(J)
         MSUM(I,J) = Sum (MNEW(:,I,J))
         if(MSUM(I,J)+MTOP > MAXCOLMASS .or.
     &      MSUM(I,J)+MTOP < MINCOLMASS) then
           n_exception = 1
           if(MSUM(I,J)+MTOP > MAXCOLMASS*(1200d0/1160d0) .or.
     &        MSUM(I,J)+MTOP < MINCOLMASS*(250d0/350d0)) then
             n_exception = 2  ! stop execution if too far out of range
           endif
         endif
      EndDo
      EndDo

      Call GLOBALMAX (grid, n_exception, n_exception_all)
      IF (n_exception_all > 0) Then
        ! need halos
        CALL HALO_UPDATE(grid, U, FROM=NORTH)
        CALL HALO_UPDATE(grid, V, FROM=NORTH)
        ! 2nd pass report problems
        Do J=J1,JN  ;  Do I =1,IMAXJ(J)
          if(MSUM(I,J)+MTOP > MAXCOLMASS .or.
     &       MSUM(I,J)+MTOP < MINCOLMASS) then
            Im1 = I-1  ;  If (Im1==0) Im1 = IM
            Write (6,990) I,J,MRCH,ZATMO(I,J),DT1,
     *           (L,U(IM1,J,L),U(I,J,L),U(IM1,J+1,L),U(I,J+1,L),
     *              V(IM1,J,L),V(I,J,L),V(IM1,J+1,L),V(I,J+1,L),
     *              MNEW(L,I,J),MOLD(L,I,J),T(I,J,L),Q(I,J,L)*1000,
     *            L=LM,1,-1)
            Write (6,*) "Pressure diagnostic error"
          endif
        EndDo  ;  EndDo
      endif

      if(n_exception_all==2)
     &     call stop_model('ADVECM: Mass diagnostic error',11)

      If (QSP)  Then
         Do I=2,IM
            MNEW(:,I,1) = MNEW(:,1,1)
            MSUM(I,1)   = MSUM(1,1)  ;  EndDo  ;  EndIf
      If (QNP)  Then
         Do I=2,IM
            MNEW(:,I,JM) = MNEW(:,1,JM)
            MSUM(I,JM)   = MSUM(1,JM)  ;  EndDo  ;  EndIf

      Call HALO_UPDATE_COLUMN (GRID, MNEW, From=SOUTH)
      Call HALO_UPDATE        (GRID, MSUM, From=SOUTH)
      Call MAtoP (MNEW)
      Return
  990 Format (/'0PRESSURE DIAGNOSTIC  I,J,MRCH,ZATMO,DT=',3I4,2F10.2/
     *  '  L     U(I-1,J)     U(I,J)   U(I-1,J+1)    U(I,J+1)',
     *      '    V(I-1,J)     V(I,J)   V(I-1,J+1)    V(I,J+1)',
     *      '       MNEW        MOLD      T(I,J)      Q*1000' /
     *  (I3,11F12.3,F12.6))
      EndSubroutine ADVECM


      Subroutine PGF (DT1, U,V,MAM, UT,VT,MAFTER, S0,SZ)
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
!**** Input: DT1 = time step (s)
!****        U,V = mean horizontal velocity during time step (m/s)
!****        MAM = mean mass distribution during time step (kg/m^2)
!****      S0,SZ = potential temperature and vertical gradient (K)
!****     MAFTER = mass distribution at end of time step (kg/m^2)
!**** Output: UT,VT = velocity updated by pressure gradient force (m/s)

!**** R (J/kg*C) = gas constant = 287 for dry air
!**** K          = exponent of exner function = R/SHA
!**** M (kg/m^2) = vertical coordinate = air mass above the level
!**** DM(kg/m^2) = layer mass difference = MAM
!**** P (Pa)     = pressure = M*GRAV
!**** DP(Pa)     = layer pressure difference = PD - PU
!**** A (m^3/kg) = specific volume = R*T / P
!**** S (K)      = potential temperature = S0 - SZ*2*(M-M0)/(MD-MU) =
!****            = S0 - SZ*2*(P-P0)/(PD-PU) = S0 - SZ*2*(P-P0)/DP =
!****            = S0+SZ*2*P0/DP - SZ*2*P/DP
!**** T (K)      = temperature = S * P(mb)^K = S*.01^K * P^K =
!****            = [(S0+SZ*2*P0/DP)*.01^K - P*(SZ*2/DP)*.01^K]*P^K =
!****            = (X - P*Y)*P^K = X*P^K - Y*P^(K+1)

!**** Integral of A*dM from MU to MD (from top to bottom of layer)
!**** Int[A*dM] = Int[R*T*dP/P*G] = R*Int{[X*P^(K-1) - Y*P^K]*dP}/G =
!**** = R*{X*P^K/K - Y*P^(K+1)/(K+1)}/G from PU to PD =
!**** = R*{X*(PD^K-PU^K)/K - Y*[PD^(K+1)-PU^(K+1)]/(K+1)}/G

!**** Compute DGZ thickness everwhere in a layer from layer bottom
!**** G*dZ = - A*dP = - (R*T/P)*dP = - R*[X*P^(K-1) - Y*P^K]*dP
!**** DGZ = - Int{R*[X*P^(K-1) - Y*P^K]*dP} from PD to P =
!****     = - R*{X*(P^K-PD^K)/K - Y*[P^(K+1)-PD^(K+1)]/(K+1)}
!****     = R*{X*(PD^K-P^K)/K - Y*[PD^(K+1)-P^(K+1)]/(K+1)}

!**** DGZup = R*{X*(PD^K-PU^K)/K - Y*[PD^(K+1)-PU^(K+1)]/(K+1)}
!**** Int[A*dM] = DGZup/G

!**** Compute mass weighted average value of DGZ in a layer
!**** DGZave = Int{DGZ*dP}/DP from PD to PU =
!**** = R*Int({X*(PD^K-P^K)/K - Y*[PD^(K+1)-P^(K+1)]/(K+1)}*dP)/DP =
!**** = R*{X*[P*PD^K - P^(K+1)/(K+1)]/K -
!****    - Y*[P*PD^(K+1) - P^(K+2)/(K+2)]/(K+1)}/DP =
!**** = R*(X*{(PD-PU)*PD^K - [PD^(K+1)-PU^(K+1)]/(K+1)}/K -
!****    - Y*{(PD-PU)*PD^(K+1) - [PD^(K+2)-PU^(K+2)]/(K+2)}/(K+1))/DP =
!**** = R*(X*{DP*PD^K - [PD^(K+1)-PU^(K+1)]/(K+1)}/K -
!****    - Y*{DP*PD^(K+1) - [PD^(K+2)-PU^(K+2)]/(K+2)}/(K+1))/DP =

!**** GZave(L) = GZATMO + Sum[DGZup(1:L-1)] + DGZave(L)

!**** PGFU (kg*m/s^2) =
!**** = {Mean[Int(A*dM)] * dP/dX + Mean(DM) * dGZave/dX} * DX*DY =
!**** = {Mean[Int(A*dM)] * dP + Mean(DM) * dGZave} * DY
!**** DUT (kg*m/s) = dTIME * PolarFiltered(PGFU)

      Use CONSTANT,   Only: GRAV,RGAS,KAPA,byGRAV,
     *                      zK=>byKAPA,zKp1=>byKAPAp1,zKp2=>byKAPAp2
      Use RESOLUTION, Only: IM,JM,LM, MTOP
      Use ATM_COM,    Only: ZATMO
      Use DIAG_COM,   Only: MODD5K
      Use GEOM,       Only: IMAXJ,DXV,DYV,DXYS,DXYN, ACOR,ACOR2
      Use ATM_COM,    Only: GZ,PHI
      Use DYNAMICS,   Only: DUT,DVT,AdM=>SPA,DO_POLEFIX,MRCH
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      Use DOMAIN_DECOMP_1D,  Only: SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      IMPLICIT NONE
      Real*8 :: DT1
      Real*8,Dimension(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *   U,V, UT,VT, S0,SZ
      Real*8,Dimension(LM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *   MAM,MAFTER
!**** Local variables
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *   P,PGFU
      Real*8 :: DT4,DGZU(LM),DGZA(LM),
     *          M,PU,PKU,PKPU,PKPPU,DP,zDP,X,Y,PD,PKD,PKPD,PKPPD,GZD,
     *          FACTOR,FLUX, VMASS
     *  ,HUNDREDTHeKAPA
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO):: FD,RFDUX
      INTEGER I,J,L,IM1,IP1,IPOLE  !@var I,J,IP1,IM1,L,IPOLE loop variab.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      DT4=DT1/4.
      HUNDREDTHeKAPA = .01d0**KAPA

!**** Loop over grid columns
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
!**** Integrate pressures from the top down
         M   = MTOP
         PU  = M*GRAV
         PKU = PU**KAPA  ;  PKPU = PKU*PU  ;  PKPPU = PKPU*PU
         Do L=LM,1,-1
            DP  = MAM(L,I,J)*GRAV
            zDP = 1 / DP
            Y   = SZ(I,J,L)*2*zDP*HUNDREDTHeKAPA
            X   = S0(I,J,L)*HUNDREDTHeKAPA + Y*(PU+.5*DP)
            PD  = PU + DP
            PKD = PD**KAPA  ;  PKPD = PKD*PD  ;  PKPPD = PKPD*PD
!           AdM = RGAS*(X*(PKD-PKU)*zK - Y*(PKPD-PKPU)*zKp1)/GRAV
            DGZU(L) = RGAS*(X*(PKD-PKU)*zK - Y*(PKPD-PKPU)*zKp1)
            DGZA(L) = RGAS*(X*(DP*PKD - (PKPD-PKPU)*zKp1)*zK -
     -                      Y*(DP*PKPD - (PKPPD-PKPPU)*zKp2)*zKp1)*zDP
            AdM(I,J,L) = DGZU(L)*byGRAV
              P(I,J,L) = GRAV*(M + .5*MAM(L,I,J))
            M   = M + MAM(L,I,J)
            PU  = PD
            PKU = PKD  ;  PKPU = PKPD  ;  PKPPU=PKPPD  ;  EndDo
!**** Integrate altitude from the bottom up
         GZD = ZATMO(I,J)
         Do L=1,LM
            GZ(I,J,L) = GZD + DGZA(L)
            GZD = GZD + DGZU(L)  ;  EndDo  ;  EndDo  ;  EndDo

C**** SET POLAR VALUES FROM THOSE AT I=1
      IF (haveLatitude(grid, J=1)) THEN
        DO L=1,LM
            AdM(2:IM,1,L) = AdM(1,1,L)
              P(2:IM,1,L) =   P(1,1,L)
             GZ(2:IM,1,L) =  GZ(1,1,L)  ;  EndDo  ;  EndIf
      IF (haveLatitude(grid, J=JM)) THEN
        DO L=1,LM
            AdM(2:IM,JM,L) = AdM(1,JM,L)
              P(2:IM,JM,L) =   P(1,JM,L)
             GZ(2:IM,JM,L) =  GZ(1,JM,L)  ;  EndDo  ;  EndIf

C****
C**** PRESSURE GRADIENT FORCE
C****
C**** NORTH-SOUTH DERIVATIVE AFFECTS THE V-COMPONENT OF MOMENTUM
C
      CALL HALO_UPDATE(grid, P,   FROM=SOUTH)
      Call HALO_UPDATE (GRID, AdM, From=SOUTH)
      Call HALO_UPDATE (GRID, GZ,  From=SOUTH)
      PHI(:,:,:) = GZ(:,:,:)
      DO 3236 L=1,LM
      DO 3236 J=J_0STG,J_1STG
      FACTOR = DT4*DXV(J)
      IM1=IM
      DO 3234 I=1,IM
      FLUX = ((AdM(I,J,L)+AdM(I,J-1,L))*( P(I,J,L)- P(I,J-1,L)) +
     +        (MAM(L,I,J)+MAM(L,I,J-1))*(GZ(I,J,L)-GZ(I,J-1,L))) *
     *       FACTOR
      DVT(I,J,L)  =DVT(I,J,L)  -FLUX
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUX
 3234 IM1=I
 3236 CONTINUE
C
C**** SMOOTHED EAST-WEST DERIVATIVE AFFECTS THE U-COMPONENT
C
      DO L=1,LM
         If (J_0STG == 2)   PGFU(:,1 ,L) = 0
         If (J_1STG == JM)  PGFU(:,JM,L) = 0
        I=IM

        DO J=Max(2,J_0STG-1),J_1STG
          DO IP1=1,IM
            PGFU(I,J,L) =
     =         (AdM(Ip1,J,L)+AdM(I,J,L))*( P(Ip1,J,L)- P(I,J,L)) +
     +         (MAM(L,Ip1,J)+MAM(L,I,J))*(GZ(Ip1,J,L)-GZ(I,J,L))
            I=IP1
          END DO
        END DO

        Call AVRX (PGFU(1,J_0H,L), JRANGE=(/Max(2,J_0H),Min(JM-1,J_1)/))

        DO J=J_0STG,J_1STG
            FACTOR = -DT4*DYV(J)
            DUT(:,J,L) = DUT(:,J,L) + FACTOR*(PGFU(:,J,L)+PGFU(:,J-1,L))
        END DO
      END DO

c correct for erroneous dxyv at the poles
      if(do_polefix.eq.1) then
         do ipole=1,2
            if(haveLatitude(grid,J=2) .and. ipole.eq.1) then
               j = 2
            else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
               j = JM
            else
               cycle
            endif
            dut(:,j,:) = dut(:,j,:)*acor
            dvt(:,j,:) = dvt(:,j,:)*acor2
         enddo
      endif
C
C**** CALL DIAGNOSTICS
      IF(MRCH.GT.0) THEN
         IF(MODD5K.LT.MRCH) CALL DIAG5D (6,MRCH,DUT,DVT)
         CALL DIAGCD (grid,3,U,V,DUT,DVT,DT1)
      ENDIF
C****
!**** Undo scaling performed at beginning of ADVECV
C****
      Do 3525 J=J_0STG,J_1STG
      I=IM
      DO 3525 IP1=1,IM
      Do L=1,LM
         VMASS = .5*((MAFTER(L,I,J-1)+MAFTER(L,Ip1,J-1))*DXYN(J-1) +
     +               (MAFTER(L,I,J  )+MAFTER(L,Ip1,J  ))*DXYS(J))
         UT(I,J,L) = UT(I,J,L) + DUT(I,J,L) / VMASS
         VT(I,J,L) = VT(I,J,L) + DVT(I,J,L) / VMASS  ;  EndDo
 3525 I = IP1

      RETURN
      END SUBROUTINE PGF

      SUBROUTINE AVRX(X,jrange)
!@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
!@auth Original development team
      USE RESOLUTION, only : im,jm
      USE GEOM, only : dlon,dxp,dyp,bydyp
      !USE DYNAMICS, only : xAVRX
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds
      USE MOMENTS, only : moment_enq_order
      USE constant, only : byrt2
      IMPLICIT NONE
      REAL*8, INTENT(INOUT), optional ::
     &     X(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      Integer, Intent(In), optional :: jrange(2)
      REAL*8, ALLOCATABLE, SAVE  :: DRAT(:)
      REAL*8, SAVE ::  BYSN(IMH)
      REAL*8, DIMENSION(0:IMH) :: AN,BN
CCC   INTEGER, SAVE :: NMIN(grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER, ALLOCATABLE, SAVE :: NMIN(:)
      INTEGER J,N
      LOGICAL, SAVE :: init = .false.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, J0, J1
      REAL*8, SAVE :: xAVRX
      INTEGER order

      if ( present(X) ) goto 1000

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO = J_0H, J_STOP_HALO = J_1H,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
C
      call moment_enq_order(order)
      if ( order==4 ) then
        xAVRX = byrt2
      else if ( order==2 ) then
        xAVRX = 1.d0
      else
        call stop_model("unsupported scheme order in AVRX0",255)
      endif
C
      IF (.NOT. init) THEN
        init = .true.
C       CALL FFT0(IM)
        j0 = MAX(1,J_0H)
        j1 = MIN(JM,J_1H)
        ALLOCATE(DRAT(j0:j1), NMIN(j0:j1))
        DO N=1,IMH
          BYSN(N)=xAVRX/SIN(.5*DLON*N)
        END DO
        DO J=j0,j1
          DRAT(J) = DXP(J)*BYDYP(3)
          DO N=IMH,1,-1
            IF(BYSN(N)*DRAT(J) .GT.1.) THEN
              NMIN(J) = N+1
              EXIT
            ENDIF
          END DO
        END DO
      END IF
      RETURN
C****
!!!      ENTRY AVRX (X)
 1000 continue
C****

      If (Present(jrange)) Then
        j0 = jrange(1)
        j1 = jrange(2)
      Else
        call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
        j0=J_0S
        j1=J_1S
      End If

      DO J=j0,j1
        IF (DRAT(J).GT.1) CYCLE
        CALL FFT (X(1,J),AN,BN)
        DO N=NMIN(J),IMH-1
          AN(N)=BYSN(N)*DRAT(J) * AN(N)
          BN(N)=BYSN(N)*DRAT(J) * BN(N)
        END DO
        AN(IMH) = BYSN(IMH)*DRAT(J) * AN(IMH)
        CALL FFTI(AN,BN,X(1,J))
      END DO

      RETURN
      END SUBROUTINE AVRX

      SUBROUTINE FILTER
!@sum  FILTER Performs 8-th order shapiro filter in zonal direction
!@auth Original development team
!@calls SHAP1D
C****
C**** MFILTR=1  SMOOTH P USING SEA LEVEL PRESSURE FILTER
C****        2  SMOOTH T USING TROPOSPHERIC STRATIFICATION OF TEMPER
C****        3  SMOOTH P AND T
C****
      Use CONSTANT,   Only: byGRAV,RGAS,SHA,KAPA
      USE RESOLUTION, only : ls1,ptop,psf,pmtop
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only : itime
      USE ATM_COM, only : t,p,q,qcl,qci,zatmo
      USE GEOM, only : areag,dxyp,byim
      USE SOMTQ_COM, only : tmom,qmom
      USE ATM_COM, only : pk
      USE DYNAMICS, only : COS_LIMIT,mfiltr,sig
#ifdef TRACERS_ON
      USE TRACER_COM, only: NTM,trm,trmom
      use OldTracer_mod, only: trname,ITIME_TR0
#endif
      USE FLUXES, only : atmsrf
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds, GLOBALSUM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: X,Y
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *        POLD, PRAT
      REAL*8 PSUMO,PSUMN,PDIF,AKAP,PS,ZS
      REAL*8, EXTERNAL :: SLP
      INTEGER I,J,L,N  !@var I,J,L  loop variables
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: KEJ,PEJ
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      REAL*8 initialTotalEnergy, finalTotalEnergy
      real*8 getTotalEnergy ! external for now
      real*8, dimension(im) :: rhosrf,pgfx

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      IF (MOD(MFILTR,2).NE.1) GO TO 200
C**** Initialise total energy (J/m^2)
      initialTotalEnergy = getTotalEnergy()

      ! Save old pressure
      do j=j_0s,j_1s
        pold(:,j)=p(:,j)
      enddo

      if(pfilter_using_slp) then
C****
C**** SEA LEVEL PRESSURE FILTER ON P
C****
      DO J=J_0S,J_1S
        DO I=1,IM
          PS=P(I,J)+PTOP
          ZS=ZATMO(I,J)*BYGRAV
          X(I,J)=SLP(PS,atmsrf%TSAVG(I,J),ZS)
          Y(I,J)=X(I,J)/PS
        END DO
      END DO
      CALL SHAP1D (8,X)
      call isotropslp(x,COS_LIMIT)
      DO J=J_0S,J_1S
        PSUMO=0.
        PSUMN=0.
        DO I=1,IM
          PSUMO=PSUMO+P(I,J)
          P(I,J)=X(I,J)/Y(I,J)-PTOP
C**** reduce large variations (mainly due to topography)
          P(I,J)=MIN(MAX(P(I,J),0.99d0*POLD(I,J)),1.01d0*POLD(I,J))
          PSUMN=PSUMN+P(I,J)
        END DO
        PDIF=(PSUMN-PSUMO)*BYIM
        DO I=1,IM
          P(I,J)=P(I,J)-PDIF
        END DO
      END DO

      else

        ! Smooth the approximate form of zonal PGF instead.
        ! This code is from ATMDYN2.f.
        do j=j_0s,j_1s
          do i=1,im
            rhosrf(i) = ((p(i,j)+ptop)**(1.-kapa))/(rgas*t(i,j,1))
          enddo
          do i=1,im-1
            pgfx(i) = (p(i+1,j)-p(i,j))+
     &           .5*(rhosrf(i+1)+rhosrf(i))*(zatmo(i+1,j)-zatmo(i,j))
          enddo
          i=im
             pgfx(i) = (p(1,j)-p(i,j))+
     &           .5*(rhosrf(1)+rhosrf(i))*(zatmo(1,j)-zatmo(i,j))
          pgfx = pgfx - sum(pgfx)*byim
          x(1,j) = 0.
          do i=2,im
            x(i,j) = x(i-1,j) + pgfx(i-1)
          enddo
          y(:,j) = x(:,j)
        enddo
        call shap1d (8,x)
        call isotropslp(x,COS_LIMIT)
        do j=j_0s,j_1s
          p(:,j) = p(:,j) + (x(:,j)-y(:,j))
        enddo

      endif ! slp-based filter or not

C**** Scale mixing ratios (incl moments) to conserve mass/heat
      DO J=J_0S,J_1S
        DO I=1,IM
          PRAT(I,J)=POLD(I,J)/P(I,J)
        END DO
      END DO
      DO L=1,LS1-1
      DO J=J_0S,J_1S
      DO I=1,IM
c adjust pot. temp. to maintain unchanged absolute temp.
        T(I,J,L)= T(I,J,L)*
     &       ((POLD(I,J)*SIG(L)+PTOP)/(P(I,J)*SIG(L)+PTOP))**KAPA
        Q(I,J,L)= Q(I,J,L)*PRAT(I,J)
!       WM(I,J,L)=WM(I,J,L)*PRAT(I,J)
        QCL(I,J,L)=QCL(I,J,L)*PRAT(I,J)
        QCI(I,J,L)=QCI(I,J,L)*PRAT(I,J)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*PRAT(I,J)
      END DO
      END DO
      END DO
#ifdef TRACERS_ON
C**** In general, only an air tracer is affected by the filter
C**** This fix conserves tracer concentration, BUT NOT MASS!
C****   It is not wanted for most tracers.
C**** Thus, this code slows model and should be removed if not wanted
C**** Instead of if(trname...) could use n=n_air and avoid the IF-test
C**** But if n_air=0 this will cause problems...
      do n=1,ntm
      if (trname(n).ne.'Air' .and. trname(n).ne.'CO2n') cycle
!     if (itime.lt.itime_tr0(n)) cycle   !probably not needed
      DO L=1,LS1-1
        DO J=J_0S,J_1S
          DO I=1,IM
             trm(I,J,L,n)=  trm(I,J,L,n)/PRAT(I,J)
             trmom(:,I,J,L,n)=trmom(:,I,J,L,n)/PRAT(I,J)
      end do; end do; end do
      end do
#endif
      CALL CALC_AMPK(LS1-1)

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)

  200 IF (MFILTR.LT.2) RETURN
C****
C**** TEMPERATURE STRATIFICATION FILTER ON T
C****
      AKAP=KAPA-.205d0    ! what is this number?
      DO L=1,LM
        IF(L.LT.LS1) THEN
          DO J=J_0S,J_1S
            Y(:,J)=(SIG(L)*P(:,J)+PTOP)**AKAP
            X(:,J)=T(:,J,L)*Y(:,J)
          END DO
          CALL SHAP1D (8,X)
          DO J=J_0S,J_1S
            T(:,J,L)=X(:,J)/Y(:,J)
          END DO
        ELSE
          DO J=J_0S,J_1S
            X(:,J)=T(:,J,L)
          END DO
          CALL SHAP1D (8,X)
          DO J=J_0S,J_1S
            T(:,J,L)=X(:,J)
          END DO
        END IF
      END DO
C
      RETURN
      END SUBROUTINE FILTER

      subroutine fltry2(q3d,strength)
!@sum  fltry2 noise reduction filter for a velocity-type field
!@sum  at secondary latitudes
      use resolution, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      use domain_decomp_1d, only : halo_update,north,south
      implicit none
      integer, parameter :: nshap=8
      real*8, parameter :: by4ton=1./(4.**nshap)
      real*8 :: dt
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: q3d
      real*8 :: strength
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: yn
      real*8 yvby4ton
      integer i,j,l,jj,n,j_0stg,j_1stg,j_1f
      real*8, dimension(im) :: yjm1,yj
      logical :: have_south_pole,have_north_pole
c      real*8 :: by4ton
c      by4ton=1./(4.**nshap)

      call getDomainBounds(grid, j_strt_stgr=j_0stg, j_stop_stgr=j_1stg,
     &         have_south_pole=have_south_pole,
     &         have_north_pole=have_north_pole)


      yvby4ton = min(strength,1d0)*by4ton*((-1)**(nshap))

      if(have_north_pole) then
        j_1f=jm-1
      else
        j_1f=j_1stg
      endif

      do l=1,lm
        do j=j_0stg,j_1stg
          yn(:,j,l)=q3d(:,j,l)
        enddo
      enddo
      do n=1,nshap
        call halo_update(grid, yn, from=north+south)
        do l=1,lm
          if(have_south_pole) then ! pole-crossing conditions
            yjm1(1:im/2)    = -yn(im/2+1:im,2,l)
            yjm1(im/2+1:im) = -yn(1:im/2,2,l)
          else
            yjm1(:)   = yn(:,j_0stg-1,l)
          endif
          do j=j_0stg,j_1f
            do i=1,im
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)+yn(i,j+1,l)
              yjm1(i) = yj(i)
            enddo
          enddo
          if(have_north_pole) then ! pole-crossing conditions
            j=jm
            do i=1,im/2
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)-yn(i+im/2,j,l)
            enddo
            do i=im/2+1,im
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)-yj(i-im/2)
            enddo
          endif
        enddo                 ! l
      enddo                   ! nshap
      do l=1,lm
        do j=j_0stg,j_1stg
          q3d(:,j,l) = q3d(:,j,l) -yn(:,j,l)*yvby4ton
        enddo
      enddo
      return
      end subroutine fltry2

      Subroutine FLTRUV
!@sum  FLTRUV Filters 2 gridpoint noise from the velocity fields
!@auth Original development team
      USE CONSTANT, only : sha
      USE RESOLUTION, only : im,jm,lm
      Use ATM_COM,    Only: MA,U,V,T, PK
      USE GEOM, only : dxyn,dxys
      USE DYNAMICS, only : dt,mrch,ang_uv, COS_LIMIT,do_polefix
     &  ,DT_XUfilter,DT_XVfilter,DT_YVfilter,DT_YUfilter
C**********************************************************************
C**** FILTERING IS DONE IN X-DIRECTION WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES.
C**********************************************************************
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *     DUT,DVT,USAVE,VSAVE
      Real*8 :: X(IM),YV(Max(2*JM,IM)), ANGM,MMUV(IM),MMUVs
      REAL*8 XUby4toN,XVby4toN,YVby4toN,YUby4toN
      REAL*8 :: DT1=0.
      INTEGER I,J,K,L,N,IP1  !@var I,J,L,N  loop variables
      REAL*8 YV2,YVJ,YVJM1,X1,XI,XIM1
      INTEGER, PARAMETER :: NSHAP=8  ! NSHAP MUST BE EVEN
      REAL*8, PARAMETER :: BY16=1./16., by4toN=1./(4.**NSHAP)
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      INTEGER :: II
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
C****
      USAVE=U ; VSAVE=V
      if (DT_XUfilter.gt.0.) then
        XUby4toN = (DT/DT_XUfilter)*by4toN
      else
        XUby4toN = 0.
      end if
      if (DT_XVfilter.gt.0.) then
        XVby4toN = (DT/DT_XVfilter)*by4toN
      else
        XVby4toN = 0.
      end if
C****
C**** Filtering in east-west direction
C****
      DO 350 L=1,LM
C**** Filter U component of velocity
      DO 240 J=J_0STG,J_1STG
      DO 210 I=1,IM
  210 X(I) = U(I,J,L)
      DO 230 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 220 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  220 XIM1 = XI
  230 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 240 I=1,IM
  240 U(I,J,L) = U(I,J,L) - X(I)*XUby4toN
C**** Filter V component of velocity
      DO 340 J=J_0STG,J_1STG
      DO 310 I=1,IM
  310 X(I) = V(I,J,L)
      DO 330 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 320 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  320 XIM1 = XI
  330 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 340 I=1,IM
  340 V(I,J,L) = V(I,J,L) - X(I)*XVby4toN
  350 CONTINUE

c This routine is now called by subroutine DYNAM
c      if(do_polefix.eq.1) then
c         call isotropuv(u,v,COS_LIMIT)
c      endif

C**** Conserve angular momentum along latitudes
      DO L=1,LM
        DO J=J_0STG,J_1STG
          ANGM=0.
          MMUVs = 0
          I=IM
          DO IP1=1,IM
             MMUV(I) = .5*((MA(L,Ip1,J-1)+MA(L,I,J-1))*DXYN(J-1) +
     +                     (MA(L,Ip1,J  )+MA(L,I,J  ))*DXYS(J))
             MMUVs = MMUVs + MMUV(I)
             ANGM = ANGM - MMUV(I)*(U(I,J,L)-USAVE(I,J,L))
            I=IP1
          END DO
          If (ANG_UV == 1)  U(:,J,L) = U(:,J,L) + ANGM/MMUVs
          DUT(:,J,L) = (U(:,J,L)-USAVE(:,J,L))*MMUV(:)
          DVT(:,J,L) = (V(:,J,L)-VSAVE(:,J,L))*MMUV(:)
        END DO
      END DO

      USAVE(:,:,:) = .5*(USAVE(:,:,:)+U(:,:,:))
      VSAVE(:,:,:) = .5*(VSAVE(:,:,:)+V(:,:,:))
      Call DIAGCD (GRID,5,USAVE,VSAVE,DUT,DVT,DT1)

      RETURN
      END SUBROUTINE FLTRUV

      subroutine isotropslp(slp,coscut)
      use RESOLUTION, only : im,jm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds
      use GEOM, only : cosp,dxp
      use DYNAMICS, only : dt
      implicit none
      real*8, parameter :: k=1d3
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: slp
      real*8 :: coscut,fac
      integer :: ipole,j,jcut,jinc,jp
      integer :: j_0s, j_1s, hemi

      call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      Do j = j_0s, j_1s
        If (far_from_pole(j, cosp, coscut)) Cycle
        fac = k*dt/(dxp(j)*dxp(j))
        Call shap1(slp(1,j),im,fac)
      enddo

      return
      end subroutine isotropslp

      subroutine isotropuv(u,v,coscut)
!@sum  isotropuv isotropizes the velocity field in the near-polar row(s)
!@auth M. Kelley
      USE RESOLUTION, only : im,jm,lm
      use DYNAMICS, only : dt
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds
      USE GEOM, only : cosv,dxv,cosi=>cosiv,sini=>siniv
      implicit none
      real*8, parameter :: klo=1d3,khi=1d7
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *  U, V
      real*8 :: coscut,fac,k
      real*8, dimension(im) :: ua,va
      real*8, dimension(0:imh) :: an,bn
      integer :: i,j,l,hemi,jp,jinc,jcut,ipole
      Integer :: J_0STG, J_1STG

      call getDomainBounds(grid, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)

      do l=1,lm
        do j = J_0STG, J_1STG
          hemi = Hemisphere(j)
          If (far_from_pole(j, cosv, coscut)) Cycle

c compute xy velocities
          do i=1,im
            ua(i) = cosi(i)*u(i,j,l)-hemi*sini(i)*v(i,j,l)
            va(i) = cosi(i)*v(i,j,l)+hemi*sini(i)*u(i,j,l)
          enddo
c filter the xy velocities

          k = maxval(abs(u(:,j,l)))*2.*dt/dxv(j)
          if(k.lt.0.5) then
            k = klo
          else if(k.gt.1d0) then
            k = khi
          else
            k = klo + 2d0*(k-0.5)*(khi-klo)
          endif
          fac = k*dt/(dxv(j)*dxv(j))
          call shap1(ua,im,fac)
          call shap1(va,im,fac)
          if(at_pole(j)) then   ! really strong filtering right at the pole
            call fft(ua,an,bn)
            an(2:imh) = 0.
            bn(2:imh) = 0.
            call ffti(an,bn,ua)
            call fft(va,an,bn)
            an(2:imh) = 0.
            bn(2:imh) = 0.
            call ffti(an,bn,va)
          endif
c convert xy velocities back to polar coordinates
          do i=1,im
            u(i,j,l) = cosi(i)*ua(i)+hemi*sini(i)*va(i)
            v(i,j,l) = cosi(i)*va(i)-hemi*sini(i)*ua(i)
          enddo
        enddo                   ! j
      enddo                     ! l

      return
      end subroutine isotropuv

      Integer Function Hemisphere(j)
      Use GEOM, only: FJEQ
      Integer :: j

      If (J < FJEQ) Then
        hemisphere = -1
      Else
        hemisphere = +1
      End If
      End Function Hemisphere

      ! Detect whether at the pole on staggered grid
      Logical Function at_pole(j)
      Use resolution, only : jm
      Integer :: j
      If (j == jm .or. j == 2) Then
        at_pole = .true.
      else
        at_pole = .false.
      end if
      End Function at_pole

      Logical Function far_from_pole(j, cosj, coscut)
      Use resolution, only: JM
        Integer :: j
        Real*8 :: cosj(JM), coscut

        far_from_pole = (cosj(j) >= coscut)

      End Function far_from_pole

      subroutine shap1(x,im,fac)
      implicit none
      integer :: im
      real*8, dimension(im) :: x
      real*8 :: fac,facby4,x1,xim1,xi
      integer :: i,n,nn
      n = int(fac) + 1
      facby4 = fac*.25d0/n
      do nn=1,n
      x1 = x(1)
      xim1 = x(im)
      do i=1,im-1
         xi = x(i)
         x(i) = x(i) + facby4*(xim1-xi-xi+x(i+1))
         xim1 = xi
      enddo
      i = im
      x(i) = x(i) + facby4*(xim1-x(i)-x(i)+x1)
      enddo
      return
      end subroutine shap1

      SUBROUTINE SHAP1D (NORDER,X)
!@sum  SHAP1D Smoothes in zonal direction use n-th order shapiro filter
!@auth Original development team
      USE RESOLUTION, only :im,jm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds
      IMPLICIT NONE
!@var NORDER order of shapiro filter (must be even)
      INTEGER, INTENT(IN) :: NORDER
      REAL*8, INTENT(INOUT),
     *        DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: X
      REAL*8, DIMENSION(IM)::XS
      REAL*8 by4toN,XS1,XSIM1,XSI
      INTEGER I,J,N   !@var I,J,N  loop variables
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      by4toN=1./4.**NORDER
      DO J=J_0S,J_1S
        XS(:) = X(:,J)
        DO N=1,NORDER
          XS1=XS(1)
          XSIM1=XS(IM)
          DO I=1,IM-1
            XSI=XS(I)
            XS(I)=XSIM1-XSI-XSI+XS(I+1)
            XSIM1=XSI
          END DO
          XS(IM)=XSIM1-XS(IM)-XS(IM)+XS1
        END DO
        X(:,J)=X(:,J)-XS(:)*by4toN
      END DO
      RETURN
      END SUBROUTINE SHAP1D

      SUBROUTINE SDRAG(DT1)
!@sum  SDRAG puts a drag on the winds in the top layers of atmosphere
!@auth Original Development Team
      Use CONSTANT,   Only: GRAV,RGAS,SHA
      USE RESOLUTION, only : ls1
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only : itime
      Use ATM_COM,    Only: MA,U,V,T, PEDN,PK
      USE GEOM, only : cosv,imaxj,kmaxj,idij,idjj,rapj,dxyv,dxyn,dxys
     *     ,rapvs,rapvn
      USE DIAG_COM, only : ajl=>ajl_loc,jl_dudtsdrg
      USE DYNAMICS, only : x_sdrag,csdragl,lsdrag
     *     ,lpsdrag,ang_sdrag,Wc_Jdrag,wmax,vsdragl
      use dynamics, only : l1_rtau,rtau,linear_sdrag
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      IMPLICIT NONE

!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var L(P)SDRAG lowest level at which SDRAG_lin is applied (near poles)
C**** SDRAG_const is applied above PTOP (150 mb) and below the SDRAG_lin
C**** regime (but not above P_CSDRAG)
      Real*8 :: X,MAUV,WL,TL,RHO,CDN,MMUV(LM),DU
!@var DUT,DVT change in momentum (kg*m/s)
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *        DUT,DVT
      INTEGER I,J,L,IP1,K,Lmax
      logical cd_lin
!@var ang_mom is the sum of angular momentun at layers LS1 to LM
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)    ::
     *        ANG_MOM, SUM_MMUV
!@var wmaxp =.75*wmax,the imposed limit for stratospheric winds (m/s)
      real*8 wmaxp,wmaxj,xjud
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      ANG_MOM(:,:) = 0  ;  DUT(:,:,:) = 0  ;  DVT(:,:,:) = 0

      if(linear_sdrag) then

        do l=l1_rtau,lm
        x = dt1/rtau(l)
        do j=j_0stg, j_1stg
        i=im
        do ip1=1,im
           MAUV = (MA(L,Ip1,J-1)+MA(L,I,J-1))*RAPVN(J-1) +
     +            (MA(L,Ip1,J  )+MA(L,I,J  ))*RAPVS(J)
c**** adjust diags for possible difference between dt1 and dtsrc
c        call inc_ajl(i,j,l,jl_dudtsdrg,-u(i,j,l)*x) ! for a-grid only
          ajl(j,l,jl_dudtsdrg) = ajl(j,l,jl_dudtsdrg) -u(i,j,l)*x
           DUT(I,J,L) = - X*MAUV*DXYV(J)*U(I,J,L)
           DVT(I,J,L) = - X*MAUV*DXYV(J)*V(I,J,L)
          u(i,j,l)=u(i,j,l)*(1.-x)
          v(i,j,l)=v(i,j,l)*(1.-x)
          ang_mom(i,j) = ang_mom(i,j) - dut(i,j,l)
          i=ip1
        enddo
        enddo
        enddo

      else ! default method

      wmaxp = wmax*3.d0/4.d0
      DO L=LS1,LM
      DO J=J_0STG, J_1STG
      cd_lin=.false.
      IF( L.ge.LSDRAG .or.
     *   (L.ge.LPSDRAG.and.COSV(J).LE..15) ) cd_lin=.true.
      wmaxj=wmax
      if(COSV(J).LE..15) wmaxj=wmaxp
      I=IM
      DO IP1=1,IM
        TL=T(I,J,L)*PK(L,I,J)   ! not quite correct - should be on UV grid
C**** check T to make sure it stayed within physical bounds
        if (TL.lt.50..or.TL.gt.423.) then
          write(99,'(a,i8,3i4,a,3f10.2)')
     *    ' SDRAG:',itime,i,j,l,'  T,U,V=',TL,U(I,J,L),V(I,J,L)
          call stop_model('Stopped in ATMDYN::SDRAG',11)
        end if
        RHO=PEDN(L+1,I,J)/(RGAS*TL)   ! not quite correct - should be on UV grid
        WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
        xjud=1.
        if(Wc_JDRAG.gt.0.) xjud=(Wc_JDRAG/(Wc_JDRAG+min(WL,wmaxj)))**2
C**** WL is restricted to Wmax by adjusting X, if necessary;
C**** the following is equivalent to first reducing (U,V), if necessary,
C**** then finding the drag and applying it to the reduced winds
                    CDN=CSDRAGl(l)*xjud
        IF (cd_lin) CDN=(X_SDRAG(1)+X_SDRAG(2)*min(WL,wmaxj))*xjud
         MAUV = (MA(L,Ip1,J-1)+MA(L,I,J-1))*RAPVN(J-1) +
     +          (MA(L,Ip1,J  )+MA(L,I,J  ))*RAPVS(J)
         X = DT1*RHO*CDN*Min(WL,WMAXJ)*GRAV*VSDRAGL(L) / MAUV
        if (wl.gt.wmaxj) X = 1. - (1.-X)*wmaxj/wl
C**** adjust diags for possible difference between DT1 and DTSRC
c        call inc_ajl(i,j,l,JL_DUDTSDRG,-U(I,J,L)*X) ! for a-grid only
        ajl(j,l,jl_dudtsdrg) = ajl(j,l,jl_dudtsdrg) -u(i,j,l)*x
         DUT(I,J,L) = - X*MAUV*DXYV(J)*U(I,J,L)
         DVT(I,J,L) = - X*MAUV*DXYV(J)*V(I,J,L)
         ANG_MOM(I,J) = ANG_MOM(I,J) -DUT(I,J,L)
        U(I,J,L)=U(I,J,L)*(1.-X)
        V(I,J,L)=V(I,J,L)*(1.-X)
        I=IP1
      END DO
      END DO
      END DO

      endif ! linear sdrag or not

C*
C***  Add the lost angular momentum uniformly back in if ang_sdrag>0
C***  only below 150mb if ang_sdrag=1, into whole column if ang_sdrag>1
C*
      if (ang_sdrag.gt.0) then
        lmax=ls1-1
        if (ang_sdrag.gt.1) lmax=lm
        SUM_MMUV(:,:) = 0
        do j = J_0STG,J_1STG
        I=IM
        do ip1 = 1,im
          do l = 1,lmax
            MMUV(L) = .5*((MA(L,Ip1,J-1)+MA(L,I,J-1))*DXYN(J-1) +
     +                    (MA(L,Ip1,J  )+MA(L,I,J  ))*DXYS(J))
            SUM_MMUV(I,J) = SUM_MMUV(I,J) + MMUV(L)
          end do
C*
          do l = 1,lmax
            DU = ANG_MOM(I,J) / SUM_MMUV(I,J)
            DUT(I,J,L) = DUT(I,J,L) + DU*MMUV(L)
c            call inc_ajl(i,j,l,JL_DUDTSDRG,du) ! for a-grid only
            ajl(j,l,jl_dudtsdrg) = ajl(j,l,jl_dudtsdrg) +du
            U(I,J,L)=U(I,J,L) + du
          end do
          I=IP1
        end do
        end do
      end if

C**** conservation diagnostic
C**** (technically we should use U,V from before but this is ok)
      CALL DIAGCD (grid,4,U,V,DUT,DVT,DT1)

      RETURN
      END SUBROUTINE SDRAG

      end module ATMDYN

      subroutine add_am_as_solidbody_rotation(u,dam)
      use constant, only : radius,mb2kg
      use resolution, only : pstrat
      use resolution, only : im,jm,lm
      use atm_com, only : p
      use geom, only : cosv,dxyn,dxys,fim
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      use domain_decomp_1d, only : globalsum
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: u
      real*8 :: dam
      integer :: j,l
      real*8 :: u0,xintsum
      real*8, dimension(grid%j_strt_halo:grid%j_stop_halo) ::
     &     psumj,xintj
      integer :: j_0stg, j_1stg, j_0, j_1
      logical :: have_south_pole, have_north_pole

      call getDomainBounds(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_stgr=j_0stg, j_stop_stgr=j_1stg,
     &               have_south_pole=have_south_pole,
     &               have_north_pole=have_north_pole)

      do j=j_0stg-1,j_1
        psumj(j) = sum(p(:,j))+fim*pstrat
      enddo
      do j=j_0stg,j_1stg
        xintj(j) = cosv(j)*cosv(j)*
     &       (psumj(j-1)*dxyn(j-1)+psumj(j)*dxys(j))
      enddo
      if(have_south_pole) xintj(1)=0.
      call globalsum(grid,xintj,xintsum,all=.true.)
      u0 = dam/(radius*mb2kg*xintsum)
      do l=1,lm
      do j=j_0stg,j_1stg
        u(:,j,l) = u(:,j,l) + u0*cosv(j)
      enddo
      enddo
      return
      end subroutine add_am_as_solidbody_rotation

      SUBROUTINE conserv_AMB_ext(U,AM)
      USE CONSTANT, only : omega,radius,mb2kg
      USE RESOLUTION, only : ls1,psfmpt,pstrat
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : p
      USE DYNAMICS, only : dsig
      USE GEOM, only : cosv,dxyn,dxys,dxyv,byaxyp
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, SOUTH, HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : CHECKSUM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: U
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
      INTEGER :: I,IP1,J,L
      REAL*8 :: PSJ,PSIJ,UE,UEDMS,FACJ

      INTEGER :: J_0S, J_1S, J_0STG, J_1STG, J_0, J_1, I_0, I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1,
     *               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)


C****
C**** ANGULAR MOMENTUM ON B GRID
C****
      CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      DO J=J_0STG,J_1STG
      PSJ=(2.*PSFMPT*DXYV(J))
      UE=RADIUS*OMEGA*COSV(J)
      UEDMS=2.*UE*PSTRAT*DXYV(J)
      FACJ=.5*COSV(J)*RADIUS*mb2kg
      I=IM
      DO IP1=1,IM
        PSIJ=(P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J))*DXYS(J)
        AM(I,J)=0.
        DO L=1,LS1-1
          AM(I,J)=AM(I,J)+U(I,J,L)*DSIG(L)
        END DO
        AM(I,J)=AM(I,J)*PSIJ
        DO L=LS1,LM
          AM(I,J)=AM(I,J)+U(I,J,L)*PSJ*DSIG(L)
        END DO
        AM(I,J)=(UEDMS+UE*PSIJ+AM(I,J))*FACJ
        I=IP1
      END DO
      END DO

      RETURN
C****
      END SUBROUTINE conserv_AMB_ext

      SUBROUTINE conserv_AM(AM)
!@sum  conserv_AM calculates A-grid column-sum atmospheric angular momentum,
!@sum  per unit area
!@auth Gary Russell/Gavin Schmidt
      USE RESOLUTION, only : im
      USE ATM_COM, only : u
      USE GEOM, only : byaxyp
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
      INTEGER :: I,J
      INTEGER :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1)

C****
C**** ANGULAR MOMENTUM ON B GRID
C****
      call conserv_AMB_ext(U,AM)

c move to A grid
      call regrid_btoa_ext(am)

c scale by area
      DO J=J_0,J_1
        DO I=I_0,I_1
          am(I,J)=am(I,J)*BYAXYP(I,J)
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE conserv_AM

      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates A-grid column-sum atmospheric kinetic energy,
!@sum  (J/m2)
!@auth Gary Russell/Gavin Schmidt
      USE CONSTANT, only : mb2kg
      USE RESOLUTION, only : ls1,psfmpt
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : p,u,v
      USE DYNAMICS, only : dsig
      USE GEOM, only : dxyn,dxys,dxyv,byaxyp
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      USE DOMAIN_DECOMP_1D, only : CHECKSUM, HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : SOUTH
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RKE
      INTEGER :: I,IP1,J,L
      INTEGER :: J_0STG,J_1STG, J_0, J_1, I_0, I_1
      REAL*8 :: PSJ,PSIJ

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1,
     *               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)

C****
C**** KINETIC ENERGY ON B GRID
C****

      CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      DO J=J_0STG,J_1STG
      PSJ=(2.*PSFMPT*DXYV(J))
      I=IM
      DO IP1=1,IM
        PSIJ=(P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J))*DXYS(J)
        RKE(I,J)=0.
        DO L=1,LS1-1
          RKE(I,J)=RKE(I,J)+
     &         (U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*DSIG(L)
        END DO
        RKE(I,J)=RKE(I,J)*PSIJ
        DO L=LS1,LM
          RKE(I,J)=RKE(I,J)+
     &         (U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*DSIG(L)*PSJ
        END DO
        RKE(I,J)=0.25*RKE(I,J)*mb2kg
        I=IP1
      END DO
      END DO

c move to A grid
      call regrid_btoa_ext(rke)

c scale by area
      DO J=J_0,J_1
        DO I=I_0,I_1
          rke(I,J)=rke(I,J)*BYAXYP(I,J)
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE conserv_KE

      SUBROUTINE calc_kea_3d(kea)
!@sum  calc_kea_3d calculates square of wind speed on the A grid
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v
c      USE GEOM, only : ravps,ravpn
      USE DOMAIN_DECOMP_ATM, only: grid
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: KEA
      INTEGER :: I,J,L
      DO L=1,LM
      DO J=GRID%J_STRT_STGR,GRID%J_STOP_STGR
      DO I=1,IM
        KEA(I,J,L)=.5*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
      ENDDO
      ENDDO
      ENDDO
      call regrid_btoa_3d(kea)
      RETURN

      END SUBROUTINE calc_kea_3d

      subroutine recalc_agrid_uv
!@sum recalc_agrid_uv Computes u_a,v_a from u and v
!@var u x-component at secondary grids (B_grid)
!@var v y-component at secondary grids (B_grid)
!@var u_a x-component at primary grids (A_grid)
!@var v_a y-component at primary grids (A_grid)
!@auth Ye Cheng
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v,ua=>ualij,va=>valij
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      Use DOMAIN_DECOMP_1D,  Only: NORTH
      USE DOMAIN_DECOMP_1D, only : halo_update
      USE GEOM, only : imaxj,idij,idjj,kmaxj,rapj,cosiv,siniv
      implicit none
      real*8, dimension(im) :: ra
      integer, dimension(im) :: idj
      real*8 :: HEMI,u_t,v_t,rak,ck,sk,uk,vk
      integer :: i,j,l,k,idik,idjk,kmax

      integer :: J_0S,J_1S
      logical :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      call getDomainBounds(grid, J_STRT_SKP=J_0S,   J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE    )
!     polar boxes

C**** Update halos of U and V
      CALL HALO_UPDATE(grid,u, from=NORTH)
      CALL HALO_UPDATE(grid,v, from=NORTH)


      if (HAVE_SOUTH_POLE) then
        J=1
        KMAX=KMAXJ(J)
        HEMI=-1.
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
      end if              !south pole
!
      if (HAVE_NORTH_POLE) then
        J=JM
        KMAX=KMAXJ(J)
        HEMI=1.
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJJ(K,J)
              RAK=RAPJ(K,J)
              ck=cosiv(k)
              sk=siniv(k)
              uk=u(idik,idjk,L)
              vk=v(idik,idjk,L)
              u_t=u_t+rak*(uk*ck-hemi*vk*sk)
              v_t=v_t+rak*(vk*ck+hemi*uk*sk)
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
      end if                !north pole

!     non polar boxes
C**** Update halos of u and v. (Needed bcs. IDJJ(3:4,J_1S)=J_1S+1)
C     ---> done by calling routine...

c      CALL HALO_UPDATE(grid, u, FROM=NORTH)
c      CALL HALO_UPDATE(grid, v, FROM=NORTH)
      DO J=J_0S,J_1S
        KMAX=KMAXJ(J)
        DO K=1,KMAX
          IDJ(K)=IDJJ(K,J)
          RA(K)=RAPJ(K,J)
        END DO
        DO I=1,IMAXJ(J)
          DO L=1,LM
            u_t=0.d0; v_t=0.d0
            DO K=1,KMAX
              IDIK=IDIJ(K,I,J)
              IDJK=IDJ(K)
              RAK=RA(K)
              u_t=u_t+u(IDIK,IDJK,L)*RAK
              v_t=v_t+v(IDIK,IDJK,L)*RAK
            END DO
            ua(l,i,j)=u_t
            va(l,i,j)=v_t
          END DO
        END DO
      END DO
C****
      return
      end subroutine recalc_agrid_uv

      subroutine regrid_atov_1d(u_a,v_a,uv1d)
      USE RESOLUTION, only : im,jm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : halo_update,SOUTH
      USE GEOM, only : rapvs,rapvn,cosiv,siniv
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo)  ::
     &          u_a,v_a
      real*8, dimension(2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)),
     &        intent(out) :: uv1d
      real*8 :: hemi
      integer :: i,ip1,j,n
      real*8, dimension(im) :: usouth,vsouth,unorth,vnorth
      integer :: j_0stg, j_1stg
      j_0stg = grid%j_strt_stgr
      j_1stg = grid%j_stop_stgr
      CALL HALO_UPDATE(grid,U_A,from=SOUTH)
      CALL HALO_UPDATE(grid,V_A,from=SOUTH)
      j=j_0stg-1
      if(j.eq.1) then
        hemi = -1.
        usouth(:)=2.*(u_a(1,j)*cosiv(:)+v_a(1,j)*siniv(:)*hemi)
        vsouth(:)=2.*(v_a(1,j)*cosiv(:)-u_a(1,j)*siniv(:)*hemi)
      else
        i=im
        do ip1=1,im
          usouth(i)=(u_a(i,j)+u_a(ip1,j))
          vsouth(i)=(v_a(i,j)+v_a(ip1,j))
          i=ip1
        enddo
      endif
      n = 0
      do j=j_0stg,j_1stg
        if(j.lt.jm) then
          i=im
          do ip1=1,im
            unorth(i)=(u_a(i,j)+u_a(ip1,j))
            vnorth(i)=(v_a(i,j)+v_a(ip1,j))
            i=ip1
          enddo
        else
          hemi = +1.
          unorth(:)=2.*(u_a(1,j)*cosiv(:)+v_a(1,j)*siniv(:)*hemi)
          vnorth(:)=2.*(v_a(1,j)*cosiv(:)-u_a(1,j)*siniv(:)*hemi)
        endif
        do i=1,im
          n = n + 1
          uv1d(n) = rapvn(j-1)*usouth(i)+rapvs(j)*unorth(i)
          n = n + 1
          uv1d(n) = rapvn(j-1)*vsouth(i)+rapvs(j)*vnorth(i)
          usouth(i) = unorth(i)
          vsouth(i) = vnorth(i)
        enddo
      enddo
      return
      end subroutine regrid_atov_1d

      subroutine get_nuv(nuv)
      use resolution, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: nuv
      nuv = 2*im*(1+grid%j_stop_stgr-grid%j_strt_stgr)
      return
      end subroutine get_nuv

      subroutine get_vpkey_of_n(n,vpkey)
      implicit none
      integer :: n,vpkey
      vpkey = 1+(n-1)/2
      return
      end subroutine get_vpkey_of_n

      subroutine get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
      use resolution, only : im
      use geom, only : rapvn,rapvs
      implicit none
      integer :: n
      integer, dimension(4) :: ilist,jlist
      real*8, dimension(4) :: wts
      integer :: nnbr
      integer :: iv,jv,ivp1
      call get_ivjv_of_n(n,iv,jv)
      nnbr = 4
      ivp1 = iv+1 - im*(iv/im)
      ilist(1:4) = (/ iv, ivp1, iv, ivp1 /)
      jlist(1:4) = (/ jv-1, jv-1, jv, jv /)
      wts(1:4) = (/ rapvn(jv-1), rapvn(jv-1), rapvs(jv), rapvs(jv) /)
      return
      end subroutine get_regrid_info_for_n

      subroutine get_uv_of_n(n,uv)
      use resolution, only : im,lm
      use atm_com, only : u,v
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: iv,jv
      call get_ivjv_of_n(n,iv,jv)
      if(mod(n,2).eq.1) then
        uv(1:lm) = u(iv,jv,1:lm)
      else
        uv(1:lm) = v(iv,jv,1:lm)
      endif
      return
      end subroutine get_uv_of_n

      subroutine store_uv_of_n(n,uv)
      use resolution, only : im,lm
      use atm_com, only : u,v
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: iv,jv
      call get_ivjv_of_n(n,iv,jv)
      if(mod(n,2).eq.1) then
        u(iv,jv,1:lm) = uv(1:lm)
      else
        v(iv,jv,1:lm) = uv(1:lm)
      endif
      return
      end subroutine store_uv_of_n

      subroutine get_ivjv_of_n(n,iv,jv)
      use resolution, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: n
      integer :: iv,jv
      integer :: nv,njm1
      nv = 1+(n-1)/2
      njm1 = (nv-1)/im
      jv = grid%j_strt_stgr + njm1
      iv = nv - njm1*im
      return
      end subroutine get_ivjv_of_n

      subroutine replicate_uv_to_agrid(ur,vr,k,ursp,vrsp,urnp,vrnp)
      use resolution, only : im,jm,lm
      USE ATM_COM, only : u,v
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : hasSouthPole, hasNorthPole
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     UR,VR
      real*8, dimension(im,lm) :: ursp,vrsp,urnp,vrnp
      integer :: i,j,l
      integer :: J_0S,J_1S
      if(k.ne.4)
     &     call stop_model('incorrect k in replicate_uv_to_agrid',255)
      J_0S = GRID%J_STRT_SKP
      J_1S = GRID%J_STOP_SKP
      do j=j_0s,j_1s
      do i=2,im
      do l=1,lm
        ur(1,l,i,j) = u(i-1,j  ,l)
        vr(1,l,i,j) = v(i-1,j  ,l)
        ur(2,l,i,j) = u(i  ,j  ,l)
        vr(2,l,i,j) = v(i  ,j  ,l)
        ur(3,l,i,j) = u(i-1,j+1,l)
        vr(3,l,i,j) = v(i-1,j+1,l)
        ur(4,l,i,j) = u(i  ,j+1,l)
        vr(4,l,i,j) = v(i  ,j+1,l)
      enddo ! l
      enddo ! i
      i = 1
      do l=1,lm
        ur(1,l,i,j) = u(im ,j  ,l)
        vr(1,l,i,j) = v(im ,j  ,l)
        ur(2,l,i,j) = u(i  ,j  ,l)
        vr(2,l,i,j) = v(i  ,j  ,l)
        ur(3,l,i,j) = u(im ,j+1,l)
        vr(3,l,i,j) = v(im ,j+1,l)
        ur(4,l,i,j) = u(i  ,j+1,l)
        vr(4,l,i,j) = v(i  ,j+1,l)
      enddo ! l
      enddo ! j
      if(hasSouthPole(grid)) then
        ursp(:,:) = u(:,2,:)
        vrsp(:,:) = v(:,2,:)
      endif
      if(hasNorthPole(grid)) then
        urnp(:,:) = u(:,jm,:)
        vrnp(:,:) = v(:,jm,:)
      endif
      return
      end subroutine replicate_uv_to_agrid

      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
      use resolution, only : im,jm,lm
      USE ATM_COM, only : u,v
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE_BLOCK,SOUTH,
     &     hasSouthPole, hasNorthPole
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     DU,DV
      real*8, dimension(im,lm) :: dusp,dvsp,dunp,dvnp
      integer :: i,j,l
      integer :: J_0STG,J_1STG
      if(k.ne.4) call stop_model(
     &     'incorrect k in avg_replicated_duv_to_vgrid',255)
      J_0STG = GRID%J_STRT_STGR
      J_1STG = GRID%J_STOP_STGR
      CALL HALO_UPDATE_BLOCK(GRID,DU,FROM=SOUTH)
      CALL HALO_UPDATE_BLOCK(GRID,DV,FROM=SOUTH)
c
c copy circumpolar data into the appropriate spots in du,dv
c
      if(hasSouthPole(grid)) then
        j=1
        do i=2,im
        do l=1,lm
          du(3,l,i,j) = dusp(i-1,l)
          du(4,l,i,j) = dusp(i  ,l)
          dv(3,l,i,j) = dvsp(i-1,l)
          dv(4,l,i,j) = dvsp(i  ,l)
        enddo
        enddo
        i=1
        do l=1,lm
          du(3,l,i,j) = dusp(im ,l)
          du(4,l,i,j) = dusp(i  ,l)
          dv(3,l,i,j) = dvsp(im ,l)
          dv(4,l,i,j) = dvsp(i  ,l)
        enddo
#ifndef ALT_CLDMIX_UV
c compensate for the factor of 2 in ravj(1).  change ravj(1) later.
        du(3:4,:,:,j) = du(3:4,:,:,j)*.5
        dv(3:4,:,:,j) = dv(3:4,:,:,j)*.5
#endif
      endif
      if(hasNorthPole(grid)) then
        j=jm
        do i=2,im
        do l=1,lm
          du(1,l,i,j) = dunp(i-1,l)
          du(2,l,i,j) = dunp(i  ,l)
          dv(1,l,i,j) = dvnp(i-1,l)
          dv(2,l,i,j) = dvnp(i  ,l)
        enddo
        enddo
        i=1
        do l=1,lm
          du(1,l,i,j) = dunp(im ,l)
          du(2,l,i,j) = dunp(i  ,l)
          dv(1,l,i,j) = dvnp(im ,l)
          dv(2,l,i,j) = dvnp(i  ,l)
        enddo
#ifndef ALT_CLDMIX_UV
c compensate for the factor of 2 in ravj(jm).  change ravj(jm) later.
        du(1:2,:,:,j) = du(1:2,:,:,j)*.5
        dv(1:2,:,:,j) = dv(1:2,:,:,j)*.5
#endif
      endif
c
c now do the averaging
c
      do j=j_0stg,j_1stg
      do i=1,im-1
      do l=1,lm
        u(i,j,l)=u(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (du(4,l,i,j-1)+du(3,l,i+1,j-1)+du(2,l,i,j)+du(1,l,i+1,j))
        v(i,j,l)=v(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (dv(4,l,i,j-1)+dv(3,l,i+1,j-1)+dv(2,l,i,j)+dv(1,l,i+1,j))
      enddo ! l
      enddo ! i
      i = im
      do l=1,lm
        u(i,j,l)=u(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (du(4,l,i,j-1)+du(3,l,1,j-1)+du(2,l,i,j)+du(1,l,1,j))
        v(i,j,l)=v(i,j,l)+
#ifdef ALT_CLDMIX_UV
     &       .25*
#endif
     &       (dv(4,l,i,j-1)+dv(3,l,1,j-1)+dv(2,l,i,j)+dv(1,l,1,j))
      enddo ! l
      enddo ! j
      return
      end subroutine avg_replicated_duv_to_vgrid

      SUBROUTINE regrid_btoa_3d(x)
      USE RESOLUTION, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH,
     &     hasSouthPole, hasNorthPole
      USE GEOM, only : byim
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: X
      INTEGER :: I,IM1,J,L
      REAL*8 :: XIM1J,XIJ
      call halo_update(grid,x,from=north)
      DO L=1,LM
      if(hasSouthPole(grid)) then
        x(:,1,l) = sum(x(:,2,l))*byim
      endif
      DO J=GRID%J_STRT_SKP,GRID%J_STOP_SKP
      IM1=IM
      XIM1J = x(im1,j,l)
      DO I=1,IM
        XIJ = x(i,j,l)
        X(I,J,L)=.25*(XIM1J+XIJ+
     &       x(im1,j+1,l)+x(i,j+1,l))
        XIM1J = XIJ
        IM1=I
      ENDDO
      ENDDO
      if(hasNorthPole(grid)) then
        x(:,jm,l) = sum(x(:,jm,l))*byim
      endif
      ENDDO
      RETURN
      END SUBROUTINE regrid_btoa_3d

      subroutine regrid_btoa_ext(x)
c regrids scalar x_bgrid*dxyv -> x_agrid*dxyp
      USE RESOLUTION, only : im,jm
      USE GEOM, only : rapvs,rapvn,dxyp,dxyv
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, HALO_UPDATE, NORTH,
     &     hasSouthPole, hasNorthPole
      USE GEOM, only : byim
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: X
      INTEGER :: I,IM1,J
      INTEGER :: J_0S,J_1S
      REAL*8 :: XIM1J,XIJ
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      call halo_update(grid,x,from=north)
      if(hasSouthPole(grid)) then
        X(:,1) = SUM(X(:,2))*BYIM*(DXYP(1)/DXYV(2))
      endif
      DO J=J_0S,J_1S
      IM1=IM
      XIM1J = X(IM1,J)
      DO I=1,IM
        XIJ = X(I,J)
c        X(I,J) = .25*(XIM1J+X(I,J)+X(IM1,J+1)+X(I,J+1))
        X(I,J) = (
     &       (XIM1J+X(I,J))*RAPVS(J)
     &      +(X(IM1,J+1)+X(I,J+1))*RAPVN(J) )
        XIM1J = XIJ
        IM1 = I
      ENDDO
      ENDDO
      if(hasNorthPole(grid)) then
        X(:,JM) = SUM(X(:,JM))*BYIM*(DXYP(JM)/DXYV(JM))
      endif
      return
      end subroutine regrid_btoa_ext


      Subroutine DIAGCD (GRID,M,UX,VX,DUT,DVT,DT1)
!@sum  DIAGCD Keeps track of the conservation properties of angular
!@+    momentum and kinetic energy inside dynamics routines
!@auth Gary Russell
      Use CONSTANT,   Only: RADIUS,OMEGA
      use resolution, only : im,jm,lm
      USE MODEL_COM, only : mdiag,mdyn
      USE GEOM, only : cosv, ravpn,ravps,bydxyp,fim,byim
      USE DIAG_COM, only : consrv=>consrv_loc
      Use DYNAMICS,   Only: CONV
      USE DOMAIN_DECOMP_1D, only : CHECKSUM, HALO_UPDATE, DIST_GRID
      Use DOMAIN_DECOMP_1D,  Only: GetDomainBounds, SOUTH
      USE GETTIME_MOD
      IMPLICIT NONE
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCD IS BEING CALLED
C**** M=1  AFTER ADVECTION IN DYNAMICS
C****   2  AFTER CORIOLIS FORCE IN DYNAMICS
C****   3  AFTER PRESSURE GRADIENT FORCE IN DYNAMICS
C****   4  AFTER STRATOS DRAG IN DYNAMICS
C****   5  AFTER FLTRUV IN DYNAMICS
C****   6  AFTER GRAVITY WAVE DRAG IN DYNAMICS
C****
      TYPE (DIST_GRID), INTENT(IN) :: grid
!@var M index denoting from where DIAGCD is called
      INTEGER, INTENT(IN) :: M
!@var DT1 current time step
      REAL*8, INTENT(IN) :: DT1
!@var UX,VX current velocities
      REAL*8, INTENT(IN),
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        UX,VX
!@var DUT,DVT current momentum changes (kg*m/s)
      REAL*8, INTENT(IN),
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        DUT,DVT
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PI
     &     ,DAMB,DKEB
      INTEGER :: I,J,L,N,IP1
      LOGICAL dopit
      REAL*8 :: DUTI,DUTIL,RKEI,RKEIL,begin
      INTEGER, DIMENSION(6) ::
     *     NAMOFM=(/2,3,4,5,6,7/), NKEOFM=(/14,15,16,17,18,19/)

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(BEGIN)

      call getDomainBounds(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

!****
!**** Mass change by advection
!****
      IF (M.eq.1) THEN
        dopit=.true.
        If (HAVE_SOUTH_POLE)  PI(1)  = FIM*Sum(CONV(1,1 ,:))
        If (HAVE_NORTH_POLE)  PI(JM) = FIM*Sum(CONV(1,JM,:))
        DO J=J_0S,J_1S
           PI(J) = Sum(CONV(:,J,:))
        END DO
      ELSE
        PI=0.
        dopit=.false.
      END IF
C****
C**** CHANGE OF ANGULAR MOMENTUM AND KINETIC ENERGY BY VARIOUS
C**** PROCESSES IN DYNAMICS
C****
C****

      CALL HALO_UPDATE(grid, PI, FROM=SOUTH)

      DO J=J_0STG,J_1STG
        DUTIL=0.
        RKEIL=0.
        DO L=1,LM
          DUTI=0.
          RKEI=0.
          DO I=1,IM
            DUTI=DUTI+DUT(I,J,L)
            RKEI=RKEI+(UX(I,J,L)*DUT(I,J,L)+VX(I,J,L)*DVT(I,J,L))
          END DO
          DUTIL=DUTIL+DUTI
          RKEIL=RKEIL+RKEI
        END DO
        if (dopit) DUTIL=DUTIL+2.*DT1*RADIUS*OMEGA*COSV(J)*
     *       (PI(J-1)*RAVPN(J-1)+PI(J)*RAVPS(J))
         DAMB(J) = DUTIL*COSV(J)*RADIUS
         DKEB(J) = RKEIL
      END DO
C****

c
c regrid to primary latitudes
c
      call regrid_to_primary_1d(damb)
      N=NAMOFM(M)
      DO J=J_0,J_1
        CONSRV(J,N)=CONSRV(J,N)+DAMB(J)*BYDXYP(J)*BYIM
      ENDDO
      call regrid_to_primary_1d(dkeb)
      N=NKEOFM(M)
      DO J=J_0,J_1
        CONSRV(J,N)=CONSRV(J,N)+DKEB(J)*BYDXYP(J)*BYIM
      ENDDO
      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAGCD


      subroutine regrid_to_primary_1d(x)
      USE RESOLUTION, only : jm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, NORTH,
     &     hasSouthPole, hasNorthPole

      implicit none
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: X
      integer :: j
      CALL HALO_UPDATE(grid, X, FROM=NORTH)
      if(hasSouthPole(grid)) X(1)=0.
      DO J=GRID%J_STRT,GRID%J_STOP_SKP
        X(J)=.5*(X(J)+X(J+1))
      ENDDO
      if(hasNorthPole(grid)) X(JM)=.5*X(JM)
      return
      end subroutine regrid_to_primary_1d

      SUBROUTINE DIAG5D (M5,NDT,DUT,DVT)
      use resolution, only : im,jm,lm
      USE MODEL_COM, only : MDIAG,MDYN
      USE DYNAMICS, only : dsig
      USE GC_COM, only : speca,nspher,klayer,jeq
      USE DIAG_COM, only : imh,fim
      USE ATMDYN, only : FCUVA,FCUVB
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      USE DOMAIN_DECOMP_1D, only : GLOBALSUM, WRITE_PARALLEL
      USE GETTIME_MOD
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        DUT,DVT  !  (kg*m/s)

      INTEGER :: M5,NDT

      REAL*8, DIMENSION(IMH+1) :: X
      REAL*8, DIMENSION(0:IMH) :: FA,FB
      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE
      REAL*8, DIMENSION
     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part
      REAL*8 BEGIN

      INTEGER :: J,J45N,KUV,KSPHER,L,MKE,N,NM
      INTEGER :: J_0STG,J_1STG

      call getDomainBounds(GRID,J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG)

      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      MKE=M5

      GO TO (810,810,810,100,100,  100,810),M5
C****  810 WRITE (6,910) M5
  810 CALL WRITE_PARALLEL(M5, UNIT=6, format=
     & "('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)")
C****  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)
      call stop_model('INCORRECT VALUE OF M5 WHEN CALLING DIAG5D',255)
C****
C**** KINETIC ENERGY
C****
C**** TRANSFER RATES FOR KINETIC ENERGY IN THE DYNAMICS
  100 CALL GETTIME(BEGIN)
      KE(:,:)=0.
      KE_part(:,:,:)=0.

      DO L=1,LM
        DO J=J_0STG,J_1STG
          KSPHER=KLAYER(L)
          IF (J > JEQ) KSPHER= KSPHER+1
          DO KUV=1,2 ! loop over u,v
            IF(KUV.EQ.1) CALL FFT(DUT(1,J,L),FA,FB)
            IF(KUV.EQ.2) CALL FFT(DVT(1,J,L),FA,FB)
            DO N=1,NM
               X(N) = .5*FIM *
     &          (FA(N-1)*FCUVA(N-1,J,L,KUV)+FB(N-1)*FCUVB(N-1,J,L,KUV))
            ENDDO
            X(1)=X(1)+X(1)
            X(NM)=X(NM)+X(NM)
            IF (J.NE.JEQ) KE_part(:,J,KSPHER)=KE_part(:,J,KSPHER)+
     &                                        X(:)*DSIG(L)
            IF (J.EQ.J45N) THEN     ! 45 N
               KE_part(:,J,KSPHER+2)=KE_part(:,J,KSPHER+2)+X(:)*DSIG(L)
            ELSE IF (J.EQ.JEQ) THEN ! EQUATOR
              DO N=1,NM
                KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                X(N)*DSIG(L)
                KE_part(N,J,KSPHER  )=KE_part(N,J,KSPHER  )+
     &                                .5D0*X(N)*DSIG(L)       ! CONTRIB TO SH
                KE_part(N,J,KSPHER+1)=KE_part(N,J,KSPHER+1)+
     &                                .5D0*X(N)*DSIG(L)       ! CONTRIB TO NH
              ENDDO
              IF (KUV.EQ.2) KSPHER=KSPHER+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      CALL GLOBALSUM(grid, KE_part(1:NM,:,1:NSPHER), KE(1:NM,1:NSPHER),
     &   ALL=.TRUE.)

      DO 180 KSPHER=1,NSPHER
      DO 180 N=1,NM
  180 SPECA(N,MKE,KSPHER)=SPECA(N,MKE,KSPHER)+KE(N,KSPHER)/NDT
C****
      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAG5D

      SUBROUTINE DIAG5F(UX,VX)
C**** FOURIER COEFFICIENTS FOR CURRENT WIND FIELD
C****
      use resolution, only : im,jm,lm
      USE MODEL_COM, only : IDACC,MDIAG,MDYN
      USE DIAG_COM, only : ia_d5f,imh
      USE ATMDYN, only : FCUVA,FCUVB
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE GETTIME_MOD
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        UX,VX
      INTEGER :: J,L
      INTEGER :: J_0STG, J_1STG
      REAL*8 BEGIN

      call getDomainBounds(GRID, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)
      CALL GETTIME(BEGIN)
      IDACC(ia_d5f)=IDACC(ia_d5f)+1
      DO L=1,LM
         DO J=J_0STG,J_1STG
            CALL FFT(UX(1,J,L),FCUVA(0,J,L,1),FCUVB(0,J,L,1))
            CALL FFT(VX(1,J,L),FCUVA(0,J,L,2),FCUVB(0,J,L,2))
         ENDDO
      ENDDO
      CALL TIMEOUT(BEGIN,MDIAG,MDYN)

      RETURN
      END SUBROUTINE DIAG5F

      SUBROUTINE QDYNAM
!@sum  QDYNAM is the driver to integrate dynamic terms by the method
!@+          of pre-computing Courant limits using mean fluxes
!@+    It replaces CALL AADVT (MMA,Q,QMOM, SD,PU,PV, DTLF,.TRUE.,
!@vers 2013/03/27
!@auth J. Lerner
      use resolution, only : im,jm,lm
      USE ATM_COM, only : q
      USE SOMTQ_COM, only : qmom
      USE DIAG_COM, only: byim
      USE GC_COM, only: agc=>agc_loc
      USE GCDIAG, only : jl_totntlh,jl_zmfntlh,jl_totvtlh,jl_zmfvtlh
      Use ATM_COM, Only: PS,MB,MMA,MWs
      USE TRACER_ADV, only:
     *    AADVQ,AADVQ0,sbf,sbm,sfbm,scf,scm,sfcm,ncyc
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      USE DOMAIN_DECOMP_1D, only : halo_update, south, north
      IMPLICIT NONE
      Real*8 :: byNCYC,byMMA
      INTEGER I,J,L   !@var I,J,L loop variables

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)


      CALL CALC_AMP(PS,MB)
      CALL HALO_UPDATE(grid, MB, FROM=SOUTH+NORTH) ! for convenience later
      CALL AADVQ0 (1._8)  ! uses the fluxes MUs,MVs,MWs from DYNAM
C****
C**** convert from concentration to mass units
C****
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
        Q(I,J,L)=Q(I,J,L)*MB(I,J,L)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*MB(I,J,L)
      enddo; enddo; enddo
C**** ADVECT
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
      CALL AADVQ (Q,QMOM, .TRUE. ,'q       ')
        byncyc = 1./ncyc
        AGC(:,:,jl_totntlh) = AGC(:,:,jl_totntlh) + sbf(:,:)
        AGC(:,:,jl_zmfntlh) = AGC(:,:,jl_zmfntlh)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        AGC(:,:,jl_totvtlh) = AGC(:,:,jl_totvtlh) + scf(:,:)
        AGC(:,:,jl_zmfvtlh)  = AGC(:,:,jl_zmfvtlh)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc
C****
C**** convert from mass to concentration units (using updated MMA)
C****
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
        byMMA = 1 / MMA(I,J,L)
        Q(I,J,L) = Q(I,J,L)*byMMA
        QMOM(:,I,J,L) = QMOM(:,I,J,L)*byMMA
      enddo; enddo; enddo

#ifndef TRACERS_ON
c Unscale the vertical mass flux accumulation for use by column physics.
c Switch the sign convention back to "positive downward".
      MWs(:,:,:) = -MWs(:,:,:)*NCYC
#else
c TRDYNAM will do the unscaling
#endif

      RETURN
      END SUBROUTINE QDYNAM
c      end module ATMDYN_QDYNAM

#ifdef TRACERS_ON
      SUBROUTINE TrDYNAM
!@sum  TrDYNAM is the driver to integrate tracer dynamic terms
!@auth J. Lerner
      use resolution, only : im,jm,lm
      USE MODEL_COM, only: itime
      USE TRACER_COM, only: trm, trmom, NTM
      use OldTracer_mod, only: itime_tr0, trname, t_qlimit
      USE TRACER_ADV
#ifndef SKIP_TRACER_DIAGS
      USE TRDIAG_COM, only: TAJLN=>TAJLN_loc, TAIJN=>TAIJN_LOC,
     *     jlnt_nt_tot,jlnt_nt_mm,jlnt_vt_tot,jlnt_vt_mm,
     *     tij_uflx,tij_vflx
#endif
      USE ATM_COM, only : MWs
      IMPLICIT NONE
      REAL*8 byncyc
      INTEGER N

C**** uses the fluxes MUs,MVs,MWs from DYNAM and QDYNAM
      DO N=1,NTM
        IF (itime.LT.itime_tr0(N)) cycle
        sfbm = 0.; sbm = 0.; sbf = 0.
        sfcm = 0.; scm = 0.; scf = 0.
        safv = 0.; sbfv = 0.

        CALL AADVQ (TRM(:,:,:,n),TrMOM(:,:,:,:,n),t_qlimit(n),trname(n))

C**** Flux diagnostics
#ifndef SKIP_TRACER_DIAGS
        byncyc = 1./ncyc
        TAJLN(:,:,jlnt_nt_tot,n) = TAJLN(:,:,jlnt_nt_tot,n) + sbf(:,:)
        TAJLN(:,:,jlnt_nt_mm, n) = TAJLN(:,:,jlnt_nt_mm, n)
     &    + sbm(:,:)*sfbm(:,:)*byim*byncyc
        TAJLN(:,:,jlnt_vt_tot,n) = TAJLN(:,:,jlnt_vt_tot,n) + scf(:,:)
        TAJLN(:,:,jlnt_vt_mm, n) = TAJLN(:,:,jlnt_vt_mm, n)
     &    + scm(:,:)*sfcm(:,:)*byim*byncyc

#ifdef TRACERS_WATER
C**** vertically integrated atmospheric fluxes
        TAIJN(:,:,tij_uflx,n) = TAIJN(:,:,tij_uflx,n) + safv(:,:)
        TAIJN(:,:,tij_vflx,n) = TAIJN(:,:,tij_vflx,n) + sbfv(:,:)
#endif
#endif

      ENDDO

c Unscale the vertical mass flux accumulation for use by column physics.
c Switch the sign convention back to "positive downward".
      MWs(:,:,:) = -MWs(:,:,:)*NCYC

      RETURN
      END SUBROUTINE TrDYNAM
#endif

      module UNRDRAG_COM
      !@sum  UNRDRAG_COM model variables for (alternative) gravity wave drag
      !@auth Tiehan Zhou / Marvin A. Geller
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR
      USE RESOLUTION, only: IM, JM
      implicit none
      save
      !@var r8: kind parameter of real*8
            integer, parameter :: r8 = selected_real_kind(12)
      !@var Z4var: Subgrid-scale orographic variances multiplied by 4 at uv grids
      !@+          of 2.5 x 2 resolution. It is real*4 rather than real*8.
      !@+          Its unit is m^2.
            real :: Z4var(IM, 2:JM)
      !@var Eke_by2: It is a tunable parameter from Eq.(3.1b) in McFarlane (JAS, 1987).
      !@+            Its unit is m^(-1).
            real(r8) :: Eke_by2 = 5.5E-6_r8
      !
      !@  Following parameters/variables are related to calculating nonorogrphic drag.
      !@+ The parameterization is described in Alexander and Dunkerton (JAS, 1999).
      !
      !   A. Parameters related to gravity wave source (assumed Gaussian shape in C):
      !
      !@var flag = 1 for B1 ( peak flux at c0 = 0 )
      !@+   flag = 0 for B2 ( peak flux at ci = 0 )
            integer, parameter :: flag = 0
      !@var Bt: sum of |momentum flux| for all +/-c (kg/m/s^2)
            real(r8) :: Bt(JM,INT_DAYS_PER_YEAR)
      !@var N_Kh: number of horizontal wavenumbers
            integer, parameter :: N_Kh = 1
      !@var Bm: amplitude for the spectrum (m^2/s^2) ~ u'w'
            real(r8), parameter :: Bm(N_Kh) = (/0.01_r8/)
      !@var Cw: half-width for the spectrum in C (m/s)
            real(r8), parameter :: Cw(N_Kh) = (/10.0_r8/)
      !@var [C_inf, C_sup]: the range of phase velocities for the Gaussian exp(-(C/Cw)**2)
            real(r8), parameter :: C_inf(N_Kh) = -3.0_r8 * Cw
            real(r8), parameter :: C_sup(N_Kh) =  3.0_r8 * Cw
      !@var N_C: number of C samples in source spectrum
            integer, parameter :: N_C = 100
      !@var dc: spectral resolution (m/s)
            real(r8), parameter :: dc(N_Kh) = (C_sup - C_inf)/(N_C - 1)
      !@var C: horizontal phase speed grid
            real(r8) :: C(N_C, N_Kh)
      !@var IZ0: vertical grid index of GW source (function of latitude)
            integer::IZ0(2:JM)
      !@var Wavelenth: wavelength of included gravity wave
      !@+              Unit of Wavelenth is km.
            real(r8), parameter :: Wavelenth(N_Kh) = (/100.0_r8/)
      !
      !    B. Other parameters and variables:
      !
      !@var N_Az: number of azimuths which is even.
            integer, parameter :: N_Az = 4
      !@var Kh: horizontal wave number grid
            real(r8) :: Kh(N_Kh)
      !@var Ah1, Ah2: used to compute components of velocity in azimuthal directions
            real(r8) :: Ah1(N_Az), Ah2(N_Az)
      !@param aLn2: ln(2.0)
            real(r8), parameter :: aLn2 = 0.69314718055994529_r8
      !@var L_min:
            integer :: L_min
      end module UNRDRAG_COM
      subroutine UNRDRAG (PB,U,V,T,SZ,UNRDRAG_x,UNRDRAG_y)
      !@sum  UNRDRAG is the driver for (alternative) gravity wave drag
      !@auth Tiehan Zhou / Marvin A. Geller
      USE UNRDRAG_COM
      USE CONSTANT, only : grav, bygrav, kapa, rgas
      USE GEOM, only: RAPVS, RAPVN
      USE RESOLUTION, only : ls1,psfmpt,ptop
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only: modelEclock
      USE DYNAMICS, only : sig,dsig,sige
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, Only : getDomainBounds
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      Use DOMAIN_DECOMP_1D,  Only: SOUTH

      implicit none
      real(r8), dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                   U, V, T, SZ, UNRDRAG_x, UNRDRAG_y
      real(r8), dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PB
      intent(inout) :: PB, T, SZ
      intent(in) :: U, V
      intent(out) :: UNRDRAG_x, UNRDRAG_y
      real(r8), parameter :: byrgas = 1.0_r8/rgas
      real(r8), parameter :: dkapa = 1.0_r8 - kapa
      real(r8), parameter :: g_sq = grav * grav
      real(r8), parameter :: byPSFMPT = 1.0_r8/PSFMPT
      real(r8), dimension(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *                   T2, SZ2
      real(r8), dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: PB2
      real(r8), dimension(LM) :: dp, P_mid
      real(r8), dimension(LM+1) :: P_edge
      real(r8), dimension(LM) :: uc, vc, rho, bvf_sq, drag_x, drag_y
      real(r8), dimension(2:LM) :: ue, ve, rhoe, bvfe
      real(r8), dimension(2:LM) :: hb, U_Comp
      real(r8), dimension(LM) :: GW, GWF_X, GWF_Y
      real(r8) :: bvf_sqe
      !@var Eps: intermittency factor
      real(r8) :: Eps
      integer :: I, J, L, IP1
      integer :: Ikh, IC, IAZ
      integer :: IC0, MC
      real :: h_4sq
      real(r8) :: byPB2, by_dp_sum
      real(r8) :: Bsum
      real(r8) :: Ugw_S, Vgw_S, U_Comp_S
      real(r8) :: SGN, x
      real(r8) :: B(N_C,N_Az,N_kh)
      real(r8) :: C_actual(N_C)
      !
      !Extract domain decomposition info
      !
      integer :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,             !&
     &         J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,   !&
     &         J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,     !&
     &         J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,     !&
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,            !&
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      if (HAVE_SOUTH_POLE) then
      do I = 2, IM
      PB(I,1) = PB(1,1)
      end do
      do L = 1, LM
      do I = 2, IM
      T(I,1,L) = T(1,1,L)
      SZ(I,1,L) = SZ(1,1,L)
      end do
      end do
      end if

      if (HAVE_NORTH_POLE) then
      do I = 2, IM
      PB(I,JM) = PB(1,JM)
      end do
      do L = 1, LM
      do I = 2, IM
      T(I,JM,L) = T(1,JM,L)
      SZ(I,JM,L) = SZ(1,JM,L)
      end do
      end do
      end if

      call HALO_UPDATE(GRID, T,  from=SOUTH)
      call HALO_UPDATE(GRID, SZ, from=SOUTH)
      call HALO_UPDATE(GRID, PB, from=SOUTH)

      do L= LS1, LM
      dp(L) = PSFMPT * DSIG(L)
      P_mid(L) = SIG(L) * PSFMPT + PTOP
      P_edge(L) = SIGE(L) * PSFMPT + PTOP
      end do
      P_edge(LM+1) = SIGE(LM+1) * PSFMPT + PTOP

      do J = J_0S, J_1
         I = IM
         do IP1 = 1, IM
            PB2(I,J) = (PB(I,J-1) + PB(IP1,J-1))*RAPVN(J-1)     !&
     &                  + (PB(I,J) + PB(IP1,J))*RAPVS(J)
            do L= 1, LM
            T2(L,I,J) = 0.25_r8 *
     *          (T(I,J-1,L) + T(IP1,J-1,L) + T(I,J,L) + T(IP1,J,L))
            SZ2(L,I,J) = 0.25_r8 *
     *       (SZ(I,J-1,L) + SZ(IP1,J-1,L) + SZ(I,J,L) + SZ(IP1,J,L))
            end do
         I = IP1
         end do
      end do

      Latitude:  do J = J_0STG, J_1STG
      Longitude: do I = 1, IM
      byPB2 = 1.0_r8/PB2(I,J)
      h_4sq = Z4var(I,J)
      !
      !Following loop calculates dp(1:LS1-1), P_mid(1:LS1-1), P_edge(1:LS1-1)
      !
      do L= 1, LS1-1
      dp(L) = PB2(I,J) * DSIG(L)
      P_mid(L) = SIG(L) * PB2(I,J) + PTOP
      P_edge(L) = SIGE(L) * PB2(I,J) + PTOP
      end do
      !
      !Following loop calculates rho(:) ,bvf_sq(:), uc(:), vc(:) at the middle levels.
      !
      do L= 1, LM
      rho(L) = P_mid(L)**dkapa / T2(L,I,J) * byrgas
      bvf_sq(L) = 2.0_r8*SZ2(L,I,J) / (dp(L)*T2(L,I,J)) * g_sq*rho(L)
      uc(L) = U(I,J,L)
      vc(L) = V(I,J,L)
      end do
      !
      !Following loop calculates rhoe(:) ,bvf_sqe(:), ue(:), ve(:) at the edge levels.
      !
      do L= 2, LM
      by_dp_sum = 1.0_r8 / (dp(L) + dp(L-1))
      rhoe(L) = (rho(L) * dp(L-1) + rho(L-1) * dp(L)) * by_dp_sum
      bvf_sqe = (bvf_sq(L) * dp(L-1) + bvf_sq(L-1) * dp(L)) * by_dp_sum
      if (bvf_sqe <= 0.0_r8) then
          bvfe(L) = -1.0
      else
          bvfe(L) = sqrt (bvf_sqe)
      end if
      ue(L) = (uc(L) * dp(L-1) + uc(L-1) * dp(L)) * by_dp_sum
      ve(L) = (vc(L) * dp(L-1) + vc(L-1) * dp(L)) * by_dp_sum
      hb(L) = P_edge(L) / rhoe(L) * bygrav
      end do

      if (h_4sq <= 0.0) then
          !For efficiency, orogrphic variances over oceans were set to a negative value.
          drag_x(:) = 0.0_r8
          drag_y(:) = 0.0_r8
      else
          call orographic_drag (ue, ve, rhoe, bvfe, h_4sq, Eke_by2
     *                           , drag_x, drag_y)
      end if
          UNRDRAG_x(I,J,1:LS1-1) = drag_x(1:LS1-1) * byPB2
          UNRDRAG_y(I,J,1:LS1-1) = drag_y(1:LS1-1) * byPB2
          UNRDRAG_x(I,J,LS1:LM) = drag_x(LS1:LM) * byPSFMPT
          UNRDRAG_y(I,J,LS1:LM) = drag_y(LS1:LM) * byPSFMPT

      !...Calculating Eps
      Eps = 0.0_r8
      do Ikh = 1, N_Kh
         Bsum = 0.0_r8
         do IC = 1, N_C
            Bsum = Bsum + Bm(Ikh) * exp(-(C(IC,Ikh)/Cw(Ikh))**2 * aLn2)
         end do
         Eps = Eps + Bsum * dc(Ikh) * ( rhoe(IZ0(J)) * 100.0_r8 )
                                                    !!!100.0_r8 arises from the units of P and rho.
      end do
      Eps = Bt(J,modelEclock%getDayOfYear()) / Eps
      !...Calculating source spectra (function of azimuth, horizontal wave number)
      do IAZ = 1, N_Az
         Ugw_S = ue(IZ0(J))
         Vgw_S = ve(IZ0(J))
         U_Comp_S = Ugw_S * Ah1(IAZ) + Vgw_S * Ah2(IAZ)
         do Ikh = 1, N_Kh
            do IC = 1, N_C
               x = C(IC,Ikh)*(1-flag) + (C(IC,Ikh) - U_Comp_S)*flag
               if ( x == 0.0_r8 ) then
                  SGN = 0.0_r8
               else
                  SGN = sign(1.0_r8, x)
               end if
               B(IC,IAZ,Ikh) = SGN * Bm(Ikh) *
     *            exp( -(C(IC,Ikh)/Cw(Ikh))**2 * aLn2 ) * rhoe(IZ0(J))
            end do
         end do
      end do

      GWF_X(:) = 0.0_r8
      GWF_Y(:) = 0.0_r8
      !...Loop over azimuth
      do IAZ = 1, N_Az
         U_Comp(:) = ue(:) * Ah1(IAZ) + ve(:) * Ah2(IAZ)
         !...Loop over horizontal wave number
         do Ikh = 1, N_kh
            IC = N_C
            IC0 = 1
            MC = 0
            do while ( B(IC,IAZ,Ikh) > 0.0_r8 )
               MC = MC + 1
               IC0 = IC
               IC = IC - 1
               if ( IC == 0) exit
            end do
         !...if flag == 0 then parameterize spectrum with c_actual = c - ucomp(iz0)
         C_actual(:) = C(:,Ikh) * real(flag, r8)  !&
     &          + ( C(:,Ikh) + U_Comp(IZ0(J)) ) * real(1 - flag, r8)
         call nonorographic_drag (C_actual(IC0), dc(Ikh), B(IC,IAZ,Ikh)
     &           , Eps, Kh(Ikh), hb, rhoe, U_Comp, bvfe, MC, IZ0(J), GW)
         GWF_X(:) = GWF_X(:) + GW(:) * Ah1(IAZ)
         GWF_Y(:) = GWF_Y(:) + GW(:) * Ah2(IAZ)
         end do      !horizontal wavenumber grid
      end do      !azimuth grid
      do L = L_min, LM
      if (L < LS1) then
         UNRDRAG_x(I,J,L) = UNRDRAG_x(I,J,L) + GWF_X(L) * byPB2
         UNRDRAG_y(I,J,L) = UNRDRAG_y(I,J,L) + GWF_Y(L) * byPB2
      else
         UNRDRAG_x(I,J,L) = UNRDRAG_x(I,J,L) + GWF_X(L) * byPSFMPT
         UNRDRAG_y(I,J,L) = UNRDRAG_y(I,J,L) + GWF_Y(L) * byPSFMPT
      end if
      end do

      end do Longitude
      end do Latitude
      end subroutine UNRDRAG

      subroutine init_UNRDRAG
      !@sum  init_UNRDRAG initializes parameters for (alternative) gravity wave drag
      !@auth Tiehan Zhou / Marvin A. Geller
      USE RESOLUTION, only: JM, LM, PLbot
      USE CONSTANT, only : pi, twopi
      USE GEOM, only: LAT_DG
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR
      USE UNRDRAG_COM, only: Z4var, Bt
      USE UNRDRAG_COM, only: r8, N_C, C_inf, dc, C, IZ0, N_Kh, Wavelenth
      USE UNRDRAG_COM, only: Kh, Ah1, Ah2, N_Az, aLn2, L_min
      USE FILEMANAGER, only: openunit, closeunit
      implicit none
      integer :: iu_Z4var, I, IAZ, J, IT
      real(r8) :: x, Phi
      real(r8) :: Bt_Smax, Bt_Nmax
      real(r8) :: Bt_Tmax
      character(Len=80) :: Title
      call openunit("Z4var", iu_Z4var, .true., .true.)
      read(iu_Z4var) Title, Z4var
      call closeunit(iu_Z4var)

      Bt_Smax = 6.0_r8 * 0.001_r8
      Bt_Nmax = 0.5_r8 * 0.001_r8
      do IT = 1, INT_DAYS_PER_YEAR
         x = cos(twopi * real(IT-16, r8) / real(INT_DAYS_PER_YEAR, r8))
         do J = 1, JM
            if ( LAT_DG(J,2) <= 1.0E-8 .and. x <= 0.0_r8 ) then
               Bt(J,IT) = -Bt_Smax *
     *          exp(-((LAT_DG(J,2) + 60.0_r8)/15.0_r8)**2 * aLn2 ) * x
            elseif ( LAT_DG(J,2) > 1.0E-8 .and. x >= 0.0_r8 ) then
               Bt(J,IT) =  Bt_Nmax *
     *          exp(-((LAT_DG(J,2) - 60.0_r8)/15.0_r8)**2 * aLn2 ) * x
            else
               Bt(J,IT) = 0.0_r8
            end if
         end do
      end do

      Bt_Tmax = 0.5_r8 * 0.001_r8
      do IT = 1, INT_DAYS_PER_YEAR
         x = cos(twopi * real(IT-16, r8) / real(INT_DAYS_PER_YEAR, r8))
         Phi = -10.0_r8 * x
         do J = 1, JM
            Bt(J,IT) = Bt(J,IT) + Bt_Tmax *
     *          exp(-( (LAT_DG(J,2) - Phi)/5.0_r8 )**2 * aLn2 ) *
     *          0.25_r8 * ( 3.0_r8 - x )
         end do
      end do

      Bt = Bt + 1.0_r8 * 0.001_r8

      do I = 1, N_C
         C(I, :) = C_inf(:) + real(I - 1, r8) * dc(:)
      end do
      do I = 1, N_Kh
      Kh(I) = twopi / (Wavelenth(I) * 1000.0_r8)
                            !!!Factor 1000.0 arises from the unit of Wavelenth.
      end do
      do IAZ = 1, N_Az
      x = twopi / real(N_Az, r8) * real(IAZ - 1, r8)
      Ah1(IAZ) = cos(x)
      Ah2(IAZ) = sin(x)
      end do
      I = 1
      do while ( PLbot(I) >= 100.0_r8 )
         I = I + 1
         if ( I == LM + 2 ) exit
      end do
         IZ0(:) = I - 1
      L_min = minval(IZ0)
      end subroutine init_UNRDRAG

      subroutine orographic_drag (u,v,rho, bvf,h_4sq,coef,drag_x,drag_y)
      !@sum   orographic_drag
      !@auth Tiehan Zhou / Marvin A. Geller
      USE CONSTANT, only : grav
      USE RESOLUTION, only: LM
      USE DYNAMICS, only : bydsig
      USE UNRDRAG_COM, only : r8
      implicit none
      real(r8), intent(in) :: coef
      real(r8), dimension(2:LM), intent(in) :: u, v, rho, bvf
      real(r8), dimension(LM), intent(out) :: drag_x, drag_y
      real, intent(in) :: h_4sq
      !@param  byFc_sq: inverse Froude number squared
      real(r8), parameter :: byFc_sq = 0.5_r8
      real(r8), parameter :: Fc_sq = 1.0_r8/byFc_sq
      real(r8) :: flux(2:LM+1)
      real(r8) :: wind(2:LM)
      real(r8) :: he_sq, u0_sq, u0, flux_temp
      real(r8) :: Fr_sq, const, Fr_sq_min
      real(r8) :: drag, proj_x, proj_y
      integer :: L
      do L = 1, LM
      flux(L+1) = 0.0_r8
      drag_x(L) = 0.0_r8
      drag_y(L) = 0.0_r8
      end do
      Fr_sq_min = Fc_sq
      u0_sq = u(2) * u(2) + v(2) * v(2)
      u0 = sqrt(u0_sq)
      if (u0 < 1.0E-8_r8) then
          return
      else
          wind(2) = u0
          do L = 3, LM
          wind(L) = (u(L) * u(2) + v(L) * v(2)) / u0
          end do
      end if

      if (bvf(2) <= 0.0_r8) then
          return
      else
          he_sq = min(real(h_4sq,r8), byFc_sq * u0_sq / (bvf(2)*bvf(2)))
      end if

      const = he_sq * rho(2) * bvf(2) * wind(2)
      flux(2) = -coef * const
      flux_temp = flux(2) * byFc_sq

      do L = 3, LM
      if (wind(L) < 1.0E-8_r8) then
          flux(L) = 0.0_r8
          exit
      elseif (bvf(L) <= 0.0_r8) then
          flux(L) = flux(L-1)
      else
          Fr_sq = rho(L) * wind(L)**3 / ( const * bvf(L) )
          if (Fr_sq >= Fr_sq_min) then
          flux(L) = flux(L-1)
          else
          flux(L) = flux_temp * Fr_sq
          Fr_sq_min = Fr_sq
          end if
      end if
      end do

      proj_x = u(2) / u0
      proj_y = v(2) / u0

      do L = 2, LM
      drag = -grav * (flux(L+1) - flux(L)) * byDSIG(L)
      drag_x(L) = drag * proj_x
      drag_y(L) = drag * proj_y
      end do
      end subroutine orographic_drag

      subroutine nonorographic_drag (c,dc,b,eps,kh,hb,rho
     *                                ,u,bf,nc,iz0, gwfrc)
      !@sum   nonorographic_drag
      !@auth Tiehan Zhou / Marvin A. Geller
      USE CONSTANT, only : grav, by3
      USE RESOLUTION, only: LM
      USE DYNAMICS, only : bydsig
      USE UNRDRAG_COM, only : r8
      !===============================================
      !...AD parameterizaion with arbitrary tabulated
      !...momentum flux (phase speed) spectrum: rho*u'*v' = b(c)*dc
      !...computes force gwfrc in a single azimuthal direction
      !
      !...Input:
      !...c(1:nc) phase speed (m/s)
      !...b(1:nc) mom flux density (kg/m**2/s)
      !...kh=horizontal wavenumber (1/m)
      !...eps=intermittency factor
      !...rho(2:LM)=density profile (kg/m**3)
      !...u(2:LM)=wind profile (m/s)
      !...bf(2:LM)=bouyancy frequency (1/s)
      !...iz0=source grid level
      !
      !...gwfrc(1:LM)=output force
      !

      implicit none

      !...arguements
      integer, intent(in) :: nc, iz0
      real(r8), intent(in), dimension(nc) :: c, b                  !global variables
      real(r8), intent(in), dimension(2:LM) :: hb, rho, u, bf   !global variables
      real(r8), intent(out) :: gwfrc(LM)                        !global variables
      real(r8), intent(in) :: kh, dc                            !global variables
      real(r8), intent(in) :: eps                               !intermittency factor
      !...The following are global variables
      real(r8) :: crfl(iz0:LM+1)      !phase speed of turning point at z
      integer :: n_rfl(iz0:LM+1)
      real(r8) :: dc3(iz0:LM+1), total_flux(iz0:LM+1)
      real(r8) :: k2, alp2, cm
      !!!...The following are local variables
      integer :: i, j, number, n1, n2, n3, n_end
      real(r8)::unsat              !unsat > 0 if unsaturated at z
      real(r8)::const, dc1, dc2, s1, s2, ratio, dc_end
      !------------------------------------------------------------------
      gwfrc(:) = 0.0_r8  !necessary when the subroutine 'unsaturated' returns quickly.

      if ( nc == 0 ) return

      k2 = kh * kh
      total_flux(:) = 0.0_r8

      !find the turning point, which must decrease with altitude
      cm = 2000.0
      do i = iz0, LM
      alp2 = 1.0_r8 / (4.0_r8 * hb(i) * hb(i) )
      if (bf(i) < 0.0_r8) then
         crfl(i) = u(i)
      else
         crfl(i) = u(i) + bf(i) / sqrt( k2 + alp2 )
      end if
      crfl(i) = min ( crfl(i), cm )
      cm = crfl(i)
      end do
      n_rfl(:) = int( ( crfl(:) - c(1) ) / dc ) + 1
      dc3(:) = crfl(:) - ( c(1) + dc * real ( n_rfl(:) - 1 ) )
      i = iz0
      if ( crfl(i) < c(nc) .and. c(1) < crfl(i) ) then
         n_end = n_rfl(i)
         dc_end = dc3(i)
      else if ( crfl(i) >= c(nc) ) then
         n_end = nc - 1
         dc_end = dc
      else
         return
      end if

      !...find unsaturated waves in the spectrum at z(i)
      if (bf(i) > 0.0_r8) then
         const = kh * rho(i) / bf(i) / 2.0_r8
         number = 0
         SOURCE: do j = 1, n_end
            unsat = const * ( c(j) - u(i) )**3 - b(j)
            if ( number == 0 .and. unsat > 0.0_r8 ) then
               number = 1
               if ( j == 1 ) then
                  n1 = j
                  dc1 = 0.0_r8
               else
                  n1 = j -1
                  s2 = const * ( c(n1) - u(i) )**3 - b(n1)
                  s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
                  ratio = s2 / ( s2 - s1 )
                  dc1 = ratio * dc
               end if
            else if ( number == 1 .and. unsat <= 0.0_r8 ) then
               number = 2
               n2 = j - 1
               s2 = const * ( c(n2) - u(i) )**3 - b(n2)
               s1 = const * ( c(n2+1) - u(i) )**3 - b(n2+1)
               ratio = s2 / ( s2 - s1 )
               dc2 = ratio * dc
            else if( ( number == 2 ) .and. ( unsat > 0.0_r8 ) ) then
               number = 3
               n3 = j - 1
               exit SOURCE
            end if
         end do SOURCE
      else
         number = 1
         n1 = 1
         dc1 = 0.0_r8
      end if

      if ( number == 0 ) return

      BAND_1_2: if ( number == 1 ) then
         n2 = n_end
         dc2 = dc_end
         call unsaturated ( n1, dc1, n2, dc2, i )
      else if ( number >= 2 ) then BAND_1_2
         call unsaturated ( n1, dc1, n2, dc2, i )
         BAND_2:if ( number == 3 ) then
            n1 = n3
            s2 = const * ( c(n1) - u(i) )**3 - b(n1)
            s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
            ratio = s2 / ( s2 - s1 )
            dc1 = ratio * dc
            n2 = n_end
            dc2 = dc_end
            call unsaturated ( n1, dc1, n2, dc2, i )
         end if BAND_2
      end if BAND_1_2

      total_flux(LM) = total_flux(LM-2) * by3
      total_flux(LM-1) = total_flux(LM-2) * 2.0_r8 * by3

      do i = iz0, LM
         gwfrc(i) = -grav * ( total_flux(i+1) - total_flux(i) ) *
     *                                          byDSIG(i) * eps
      end do

      contains

      recursive subroutine unsaturated ( nz1, dcz1, nz2, dcz2, level)
      implicit none
      integer, intent(in)::nz1, nz2, level
      real(r8), intent(in)::dcz1, dcz2
      integer::number, n1, n2, n3, n_end, n_start
      integer::i, j
      real(r8)::const, unsat, dc1, dc2, s1, s2, ratio, flux_rfl, dc_end

      i = level
      do j = nz1, nz2 - 1
         total_flux(i) = total_flux(i) + b(j) * dc
      end do
      total_flux(i) = total_flux(i) - b(nz1) * dcz1 + b(nz2) * dcz2

      if (level == LM) return

      i = level + 1
      if ( c(nz1) + dcz1 < crfl(i) .and. crfl(i) < c(nz2) + dcz2 ) then
         n_end = n_rfl(i)
         dc_end = dc3(i)
      else if ( crfl(i) >= c(nz2) + dcz2 ) then
         n_end = nz2
         dc_end = dcz2
      else
         n_end = nz1
         dc_end = dcz1
      end if

      flux_rfl = 0.0_r8
      do j = n_end, nz2 - 1
      flux_rfl =  flux_rfl + b(j) * dc
      end do
      flux_rfl = flux_rfl - b(n_end) * dc_end + b(nz2) * dcz2
      total_flux(:level) = total_flux(:level) - flux_rfl

      if ( crfl(i) <= c(nz1) + dcz1 ) return

      if (bf(i) > 0.0_r8) then
         const = kh * rho(i) / bf(i) / 2.0_r8
         number = 0
         if ( nz1 == 1 .and. dcz1 == 0.0_r8 ) then
            n_start = nz1
         else
            n_start = nz1 + 1
         end if
         find_bands: do j = n_start, n_end
            unsat = const * ( c(j) - u(i) )**3 - b(j)
            upward: if ( number == 0 .and. unsat > 0.0_r8 ) then
               number = 1
               if ( j == 1 ) then
                  n1 = j
                  dc1 = 0.0_r8
               else
                  n1 = j - 1
                  s2 = const * ( c(n1) - u(i) )**3 - b(n1)
                  s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
                  ratio = s2 / ( s2 - s1 )
                  dc1 = ratio * dc
                  if ( n1 == nz1 .and. s1 * s2 <= 0.0_r8 ) then
                     dc1 = max ( dc1, dcz1 )
                  else if ( n1 == nz1 .and. s1 * s2 > 0.0_r8 ) then
                     dc1 = dcz1
                  end if
               end if
            else if ( number == 1 .and. unsat <= 0.0_r8 ) then upward
               number = 2
               n2 = j - 1
               s2 = const * ( c(n2) -u(i) )**3 - b(n2)
               s1 = const * ( c(n2+1) - u(i) )**3 - b(n2+1)
               ratio = s2 / ( s2 - s1 )
               dc2 = ratio * dc
               if ( n2 == n_end ) dc2 = min ( dc2, dc_end )
            else if ( number == 2 .and. unsat > 0.0_r8 ) then upward
               number = 3
               n3 = j - 1
               exit find_bands
            end if upward
         end do find_bands
      else
         number = 1
         n1 = nz1
         dc1 = dcz1
      end if

      if ( number == 0 ) return

      outgoing_1_2: if ( number == 1 ) then
         n2 = n_end
         dc2 = dc_end
         call unsaturated ( n1, dc1, n2, dc2, i )
      else if ( number >= 2 ) then outgoing_1_2
         call unsaturated ( n1, dc1, n2, dc2, i )
         outgoing_2: if ( number == 3 ) then
            n1 = n3
            s2 = const * ( c(n1) - u(i) )**3 - b(n1)
            s1 = const * ( c(n1+1) - u(i) )**3 - b(n1+1)
            ratio = s2 / ( s2 - s1 )
            dc1 = ratio * dc
            n2 = n_end
            dc2 = dc_end
            call unsaturated ( n1, dc1, n2, dc2, i )
         end if outgoing_2
      end if outgoing_1_2

      end subroutine unsaturated
      end subroutine nonorographic_drag

