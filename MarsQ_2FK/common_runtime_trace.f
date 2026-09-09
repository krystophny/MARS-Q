C=======================================================================
C DEFAULT-OFF ORDERED TRAPPED FULL-BOUNCE GEOMETRY PACKET.
C
C KBTIME and KPHI produce one lower-to-upper zero-orbit-width leg.  The
C return leg retraces those source-native positions in reverse while physical
C time keeps increasing.  This writer emits that exact discrete reflection.
C It does not infer a finite-orbit toroidal position or split the KG/KH cosine
C into leg-resolved complex amplitudes.
C=======================================================================
      SUBROUTINE WRITEKJPTRAPPEDORBIT(JS,JS_MAT,KGRID,RLAM)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ToolBox
      IMPLICIT NONE

      INTEGER JS,JS_MAT,KGRID,J,IOUT,NHALF,FID,JHALF,LEG,
     &        ENDPOINT,ENDPOINTOWNER,PERIODICCLOSURE,KSTATUS
      REAL*8 RLAM,DPSIS,RADIAL,TAUEND,HCHIFACTOR,VPABS,
     &       CHI,PHI,TAU,TAUFRACTION,BOUNCEANGLE,BOUNCEPHASE,
     &       VPSTATE,DTAUMEASURE,ORIENTSTATE,ORIENTVPAR
      REAL*8 HCHIHALF(NCHIT+2),VPABSHALF(NCHIT+2)
      LOGICAL OTRACE,OEXIST,OFAILED
      CHARACTER*128 PATH

      IF (KGRID.NE.1.AND.KGRID.NE.2) RETURN
      CALL KELLTRACESELECT(JS,KGRID,OTRACE)
      IF (.NOT.OTRACE) RETURN
C     The reflected geometry is exact only for the zero-FOW KJPCOEFF lane.
      IF (IFOWT.NE.0) RETURN
      IF (.NOT.ALLOCATED(RCHIK).OR..NOT.ALLOCATED(RPHIK).OR.
     &    .NOT.ALLOCATED(RTK).OR..NOT.ALLOCATED(RHK).OR.
     &    .NOT.ALLOCATED(RJBK)) RETURN

      IF (KGRID.EQ.1) THEN
         DPSIS = DPSIDS(JS)
         RADIAL = CS(JS)
      ELSE
         DPSIS = DPSIDSM(JS)
         RADIAL = CSM(JS)
      ENDIF
      IF (DPSIS.EQ.0.D0.OR.RLAM.LE.0.D0) RETURN

      NHALF = NCHI2+2
      IF (NHALF.LT.2) RETURN
      TAUEND = RTK(NHALF)
      IF (TAUEND.LE.0.D0) RETURN
      IF (ABS(OMEGAB*TAUEND-PI).GT.
     &    1.D-12*MAX(1.D0,ABS(OMEGAB*TAUEND),PI)) RETURN
      DO J=2,NHALF
         IF (RTK(J).LT.RTK(J-1)) RETURN
      ENDDO

      DO J=1,NHALF
         IF (RJBK(J).NE.0.D0) THEN
            HCHIHALF(J) = DPSIS/RJBK(J)
         ELSE
            HCHIHALF(J) = 0.D0
         ENDIF
         VPABSHALF(J) = 0.D0
         IF (J.GT.1.AND.J.LT.NHALF) THEN
            VPABS = 1.D0-RLAM/RHK(J)
            IF (VPABS.GT.0.D0) VPABSHALF(J) = SQRT(VPABS)
         ENDIF
      ENDDO

      WRITE(PATH,'("ELL_TRACE_JS",I4.4,"_G",I1,"_ORBIT.OUT")')
     &      JS,KGRID
C$OMP CRITICAL(ELL_TRACE_WRITE)
      INQUIRE(FILE=PATH,EXIST=OEXIST)
      FID=ASSIGNFREEFILEUNIT()
      OPEN(FID,FILE=PATH,STATUS='UNKNOWN',POSITION='APPEND',
     &     ACTION='WRITE')
      IF (.NOT.OEXIST) THEN
         WRITE(FID,'(A)') '# schema: iter-tc24-mars-trapped-full-orbit-v1'
         WRITE(FID,'(A)') '# source: KJPCOEFF trapped KPARTICLE=0 IFOWT=0'
         WRITE(FID,'(A,I8)') '# js = ',JS
         WRITE(FID,'(A,I8)') '# js_mat = ',JS_MAT
         WRITE(FID,'(A,I8)') '# kgrid = ',KGRID
         WRITE(FID,'(A,ES24.16)') '# rho_pol = ',RADIAL
         WRITE(FID,'(A)') '# orbit_span = full_bounce_reflected_native_zero_fow'
         WRITE(FID,'(A)') '# start_end_point = lower_turning_point'
         WRITE(FID,'(A)') '# intermediate_turning_point = upper_turning_point'
         WRITE(FID,'(A)') '# endpoint_policy = own_initial_lower_and_upper_once; final_lower_is_closure'
         WRITE(FID,'(A)') '# time_orientation = increasing_native_RTK_on_both_legs'
         WRITE(FID,'(A)') '# position_return = exact_reverse_of_native_RCHIK_RPHIK_samples'
         WRITE(FID,'(A)') '# delta_tau_measure = positive_native_RTK_cell_measure; first_row_zero'
         WRITE(FID,'(A)') '# complex_amplitude = unavailable_not_inferred_from_folded_KG_KH_cosine'
         WRITE(FID,'(A)') '# columns: js js_mat kgrid kparticle sample_index half_sample_index leg '
     &      //'lambda ell chi phi tau tau_fraction bounce_angle bounce_phase rho_pol '
     &      //'b0_over_b b_norm jb dpsids hchi vpar_state delta_tau_measure '
     &      //'orientation_state orientation_vpar endpoint_flag endpoint_owner '
     &      //'periodic_closure'
      ENDIF

      OFAILED = .FALSE.
      DO IOUT=1,2*NHALF-1
         CALL KTRAPPEDFULLPOINT(NHALF,IOUT,DBLE(KNTVELL),TAUEND,
     &      OMEGAB,RCHIK,RPHIK,RTK,HCHIHALF,VPABSHALF,JHALF,LEG,
     &      CHI,PHI,TAU,TAUFRACTION,BOUNCEANGLE,BOUNCEPHASE,
     &      HCHIFACTOR,VPSTATE,DTAUMEASURE,ORIENTSTATE,ORIENTVPAR,
     &      ENDPOINT,ENDPOINTOWNER,PERIODICCLOSURE,KSTATUS)
         IF (KSTATUS.NE.0) THEN
            OFAILED = .TRUE.
            EXIT
         ENDIF
         WRITE(FID,1000) JS,JS_MAT,KGRID,0,IOUT,JHALF,LEG,
     &      RLAM,DBLE(KNTVELL),CHI,PHI,TAU,TAUFRACTION,BOUNCEANGLE,
     &      BOUNCEPHASE,RADIAL,RHK(JHALF),B0K/RHK(JHALF),RJBK(JHALF),
     &      DPSIS,HCHIFACTOR,VPSTATE,DTAUMEASURE,ORIENTSTATE,
     &      ORIENTVPAR,ENDPOINT,ENDPOINTOWNER,PERIODICCLOSURE
      ENDDO
      CLOSE(FID)
C$OMP END CRITICAL(ELL_TRACE_WRITE)
      IF (OFAILED) RETURN

      RETURN
 1000 FORMAT(7I8,18(1X,E24.16),3(1X,I3))
      END SUBROUTINE WRITEKJPTRAPPEDORBIT

C=======================================================================
C DEFAULT-OFF COMMON ORBIT RESPONSE PACKET.
C
C This writer exposes the native trapped KJPCOEFF integrands at the exact
C quadrature points used by KG and KH.  It is diagnostic-only: no production
C array is modified and no torque/sign/phase is fitted.  The selected
C surface and integer ell are the same entries in ELL_TRACE.REQUEST used by
C the existing KG/KH/action traces.  The request must therefore use a scalar
C KNTVELL, not the assembled 999 selector.
C
C Interior rows carry the regular and endpoint-subtracted G/H integrands; the
C two endpoint rows carry the analytic singular add-backs.  Summing the rows
C and applying the recorded G/H normalisations reproduces the native factor
C formulas.  RCHIK/RPHIK/RTK are MARS's ordered trapped half-bounce chart.
C `rho_pol` is CS (full mesh) or CSM (half mesh); a physical s_tor map is not
C invented here and must be joined from the accepted radial map.  Schema v2
C additionally exports the chart factor hchi=B.grad(chi)/B and the physical
C one-sided parallel orientation.  This is diagnostic metadata only.
C=======================================================================
      SUBROUTINE WRITEKJPCOMMONTRACE(JS,JS_MAT,KGRID,RLAM)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ToolBox
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER JS,JS_MAT,KGRID,J,K,L,FID
      REAL*8 RLAM,OMEGAE,DPSIS,RADIAL,TAUEND,GNORM,HNORM,
     &       ARG,SQRTARG,PGPHASE,PHPHASE,GP1,GP2,GP3,
     &       HP1,HP6,HP2,HP3,HP4,HP5,HP7,CTMPL,CTMPU,
     &       PHASE0,DIFFERCHI,DIFPI,HCHIFACTOR,VPSTATE,
     &       ORIENTSTATE,ORIENTVPAR
      COMPLEX*16 PG,PH,FLG,FUG,FLX1,FUX1,FLX2,FUX2,
     &       FLQ1,FUQ1,FLQ2,FUQ2,FLQ3,FUQ3,FLDP,FUDP,
     &       GPARA,GPERP,GDPHI,HX1,HX2,HQ1,HQ2,HQ3,HDP,ZERO,
     &       GENDP,HENDX1,HENDX2,HENDQ1,HENDQ2,HENDQ3,HENDDP,
     &       CALPHA_LOCAL
      LOGICAL OTRACE,OEXIST
      CHARACTER*128 PATH

      IF (KGRID.NE.1.AND.KGRID.NE.2) RETURN
      CALL KELLTRACESELECT(JS,KGRID,OTRACE)
      IF (.NOT.OTRACE) RETURN
      IF (.NOT.ALLOCATED(RCHIK).OR..NOT.ALLOCATED(RPHIK).OR.
     &    .NOT.ALLOCATED(RTK).OR..NOT.ALLOCATED(RHK).OR.
     &    .NOT.ALLOCATED(RJBK)) RETURN

      IF (KGRID.EQ.1) THEN
         OMEGAE = ROT(JS)
         DPSIS  = DPSIDS(JS)
         RADIAL = CS(JS)
      ELSE
         OMEGAE = ROTM(JS)
         DPSIS  = DPSIDSM(JS)
         RADIAL = CSM(JS)
      ENDIF
      IF (DPSIS.EQ.0.D0) RETURN
      IF (RLAM.LE.0.D0) RETURN

      TAUEND = RTK(NCHI2+2)
      IF (TAUEND.LE.0.D0) RETURN
      GNORM = RCHIHK/4.D0/PI
      HNORM = RCHIHK/4.D0*OMEGAB/PI
      PHASE0 = 4.D0*SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK
      IF (IPERTURB.NE.0) THEN
C        CALPHA is zero in the executed perturbation mode, as in KH.
         CALPHA_LOCAL = 0.D0
      ELSEIF (V2XKEY.EQ.1 .OR. V2XKEY.EQ.3) THEN
         CALPHA_LOCAL = 0.D0
      ELSE
         CALPHA_LOCAL = CI/(OMEGA-RNTOR*OMEGAE)
      ENDIF
      ZERO = DCMPLX(0.D0,0.D0)

      WRITE(PATH,'("ELL_TRACE_JS",I4.4,"_G",I1,"_COMMON.OUT")')
     &      JS,KGRID
C$OMP CRITICAL(ELL_TRACE_WRITE)
      INQUIRE(FILE=PATH,EXIST=OEXIST)
      FID=ASSIGNFREEFILEUNIT()
      OPEN(FID,FILE=PATH,STATUS='UNKNOWN',POSITION='APPEND',
     &     ACTION='WRITE')
      IF (.NOT.OEXIST) THEN
         WRITE(FID,'(A)') '# schema: iter-tc24-mars-common-orbit-trace-v2'
         WRITE(FID,'(A)') '# source: KJPCOEFF trapped KPARTICLE=0'
         WRITE(FID,'(A,I8)') '# js = ',JS
         WRITE(FID,'(A,I8)') '# js_mat = ',JS_MAT
         WRITE(FID,'(A,I8)') '# kgrid = ',KGRID
         WRITE(FID,'(A,ES24.16)') '# rho_pol = ',RADIAL
         WRITE(FID,'(A)') '# radial_coordinate = CS/CSM=rho_pol; join s_tor externally'
         WRITE(FID,'(A)') '# position_coordinates = Boozer(chi,phi,tau)'
         WRITE(FID,'(A)') '# orbit_class = trapped'
         WRITE(FID,'(A)') '# orbit_span = half_bounce'
         WRITE(FID,'(A)') '# start_point = lower_turning_point'
         WRITE(FID,'(A)') '# end_point = upper_turning_point'
         WRITE(FID,'(A)') '# endpoint_bounce_angle = pi'
         WRITE(FID,'(A)') '# time_orientation = increasing_native_RTK'
         WRITE(FID,'(A)') '# phase_gauge = chi=RCHIK(1), phi=RPHIK(1)=0, tau=RTK(1)=0'
         WRITE(FID,'(A)') '# tau = native MARS normalized bounce-time coordinate'
         WRITE(FID,'(A)') '# endpoint_flag: -1 lower add-back, 0 interior, +1 upper add-back'
         WRITE(FID,'(A)') '# state_velocity_convention = vpar_state = sign(hchi)*(v_parallel/v)'
         WRITE(FID,'(A)') '# chart_factor = hchi = B_dot_grad_chi/B = dpsids/jb'
         WRITE(FID,'(A)') '# orientation_convention = sign(v_parallel) = sign(vpar_state*hchi)'
         WRITE(FID,'(A)') '# orientation_zero = 0 when vpar_state or hchi is zero'
         WRITE(FID,'(A)') '# columns: js js_mat kgrid kparticle m_index sample_index lambda m ell '
     &      //'chi phi tau tau_fraction bounce_angle rho_pol b0_over_b b_norm jb '
     &      //'gphase_re gphase_im gpara_re gpara_im gperp_re gperp_im '
     &      //'gdphi_re gdphi_im hphase_re hphase_im hx1_re hx1_im hx2_re hx2_im '
     &      //'hq1_re hq1_im hq2_re hq2_im hq3_re hq3_im hdp_re hdp_im '
     &      //'g_normalization h_normalization dpsids hchi vpar_state '
     &      //'orientation_state orientation_vpar endpoint_flag'
      ENDIF

      DO L=1,MLMAX
         IF (NINT(RLM(L)).NE.KNTVELL) CYCLE
         DO K=1,MSMAX
C           KG endpoint coefficients (singular terms and add-backs).
            FLG = RJBK(1)/2.D0*SQRT(RLAM)/SQRT(HPL)*
     &            EXP(-CI*RM(K,2)*CHIL)
            FUG = RJBK(NCHI2+2)/2.D0*SQRT(RLAM)/SQRT(-HPU)*
     &            EXP(CI*(-RNTOR*RPHIK(NCHI2+2)-RM(K,2)*CHIU))*
     &            COS(RLM(L)*PI)

C           KH endpoint coefficients, copied term-for-term from KH.
            FLX1 = (RX1BK(1)-RX1RK(1)*CALPHA_LOCAL)*
     &             SQRT(RLAM)*EXP(CI*RM(K,2)*CHIL)*2.D0/SQRT(HPL)
            FLX2 = RX2K(1)*SQRT(RLAM)*EXP(CI*RM(K,2)*CHIL)*2.D0/SQRT(HPL)
            FLQ1 = RQ1K(1)*RLAM*SQRT(RLAM)*EXP(CI*RM(K,2)*CHIL)*
     &             2.D0/SQRT(HPL)
            FLQ2 = RQ2K(1)*RLAM*SQRT(RLAM)*EXP(CI*RM(K,2)*CHIL)*
     &             2.D0/SQRT(HPL)
            FLQ3 = RQ3K*RLAM*SQRT(RLAM)*EXP(CI*RM(K,2)*CHIL)*
     &             2.D0/SQRT(HPL)
            FLDP = RJBK(1)/DPSIS*SQRT(RLAM)*EXP(CI*RM(K,2)*CHIL)*
     &             2.D0/SQRT(HPL)
            FUX1 = (RX1BK(NCHI2+2)-RX1RK(NCHI2+2)*CALPHA_LOCAL)*
     &             SQRT(RLAM)*EXP(CI*(RM(K,2)*CHIU+
     &             RNTOR*RPHIK(NCHI2+2)))*2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
            FUX2 = RX2K(NCHI2+2)*SQRT(RLAM)*
     &             EXP(CI*(RM(K,2)*CHIU+RNTOR*RPHIK(NCHI2+2))) *
     &             2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
            FUQ1 = RQ1K(NCHI2+2)*RLAM*SQRT(RLAM)*
     &             EXP(CI*(RM(K,2)*CHIU+RNTOR*RPHIK(NCHI2+2))) *
     &             2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
            FUQ2 = RQ2K(NCHI2+2)*RLAM*SQRT(RLAM)*
     &             EXP(CI*(RM(K,2)*CHIU+RNTOR*RPHIK(NCHI2+2))) *
     &             2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
            FUQ3 = RQ3K*RLAM*SQRT(RLAM)*
     &             EXP(CI*(RM(K,2)*CHIU+RNTOR*RPHIK(NCHI2+2))) *
     &             2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
            FUDP = RJBK(NCHI2+2)/DPSIS*SQRT(RLAM)*
     &             EXP(CI*(RM(K,2)*CHIU+RNTOR*RPHIK(NCHI2+2))) *
     &             2.D0*COS(RLM(L)*PI)/SQRT(-HPU)

            GENDP = FLG*4.D0*SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK
            HENDX1 = FLX1*PHASE0
            HENDX2 = FLX2*PHASE0
            HENDQ1 = FLQ1*PHASE0
            HENDQ2 = FLQ2*PHASE0
            HENDQ3 = FLQ3*PHASE0
            HENDDP = FLDP*PHASE0

C           Lower endpoint add-back row.  Emit only the lower coefficient;
C           the upper coefficient is emitted in the separate upper row below.
C           This makes the row sum equal the one native (FLG+FUG) add-back.
            IF (RJBK(1).NE.0.D0) THEN
               HCHIFACTOR = DPSIS/RJBK(1)
            ELSE
               HCHIFACTOR = 0.D0
            ENDIF
            VPSTATE = 0.D0
            ORIENTSTATE = 0.D0
            ORIENTVPAR = 0.D0
            CALL WRITEKJPCOMMONROW(FID,JS,JS_MAT,KGRID,0,K,1,RLAM,
     &       RM(K,2),RLM(L),RCHIK(1),RPHIK(1),RTK(1),0.D0,
     &       OMEGAB*RTK(1),RADIAL,RHK(1),B0K/RHK(1),RJBK(1),ZERO,ZERO,
     &       GENDP,GENDP,ZERO,HENDX1,HENDX2,HENDQ1,HENDQ2,HENDQ3,HENDDP,
     &       GNORM,HNORM,DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,ORIENTVPAR,-1)

            DO J=2,NCHI2+1
               ARG=1.D0-RLAM/RHK(J)
               IF (ARG.LE.0.D0) CYCLE
               SQRTARG=SQRT(ARG)
               PGPHASE=-RNTOR*RPHIK(J)-RM(K,2)*RCHIK(J)
               PG=EXP(CI*PGPHASE)*COS(RLM(L)*OMEGAB*RTK(J))
               GP1=RJBK(J)*SQRTARG
               GP2=0.5D0*RJBK(J)*RLAM/RHK(J)/SQRTARG
               GP3=0.5D0*RJBK(J)/SQRTARG
               CTMPL=1.D0/SQRT(DIFPI(RCHIK(J)-CHIL))
               CTMPU=1.D0/SQRT(DIFPI(CHIU-RCHIK(J)))
               GPARA=GP1*PG
               GPERP=GP2*PG-CTMPL*FLG-CTMPU*FUG
               GDPHI=GP3*PG-CTMPL*FLG-CTMPU*FUG

               PHPHASE=RNTOR*RPHIK(J)+RM(K,2)*RCHIK(J)
               PH=EXP(CI*PHPHASE)*2.D0*COS(RLM(L)*OMEGAB*RTK(J))
               HP1=((1.D0-RLAM/RHK(J))*RX1PK(J)+
     &              (2.D0-RLAM/RHK(J))*RX1BK(J))/SQRTARG
               HP6=(2.D0-RLAM/RHK(J))*RX1RK(J)/SQRTARG
               HP2=(2.D0-RLAM/RHK(J))*RX2K(J)/SQRTARG
               HP3=RLAM*RQ1K(J)/SQRTARG
               HP4=RLAM*RQ2K(J)/SQRTARG
               HP5=RLAM*RQ3K/SQRTARG
               HP7=RJBK(J)/DPSIS/SQRTARG
               HX1=PH*(HP1-HP6*CALPHA_LOCAL)-CTMPL*FLX1-CTMPU*FUX1
               HX2=PH*HP2-CTMPL*FLX2-CTMPU*FUX2
               HQ1=PH*HP3-CTMPL*FLQ1-CTMPU*FUQ1
               HQ2=PH*HP4-CTMPL*FLQ2-CTMPU*FUQ2
               HQ3=PH*HP5-CTMPL*FLQ3-CTMPU*FUQ3
               HDP=PH*HP7-CTMPL*FLDP-CTMPU*FUDP
               IF (RJBK(J).NE.0.D0) THEN
                  HCHIFACTOR = DPSIS/RJBK(J)
               ELSE
                  HCHIFACTOR = 0.D0
               ENDIF
               VPSTATE = SQRTARG
               IF (VPSTATE.EQ.0.D0 .OR. HCHIFACTOR.EQ.0.D0) THEN
                  ORIENTSTATE = 0.D0
                  ORIENTVPAR = 0.D0
               ELSE
                  ORIENTSTATE = SIGN(1.D0,VPSTATE)
                  ORIENTVPAR = SIGN(1.D0,VPSTATE*HCHIFACTOR)
               ENDIF
               CALL WRITEKJPCOMMONROW(FID,JS,JS_MAT,KGRID,0,K,J,RLAM,
     &          RM(K,2),RLM(L),RCHIK(J),RPHIK(J),RTK(J),
     &          RTK(J)/TAUEND,OMEGAB*RTK(J),RADIAL,RHK(J),B0K/RHK(J),RJBK(J),
     &          PG,GPARA,GPERP,GDPHI,PH,HX1,HX2,HQ1,HQ2,HQ3,HDP,
     &          GNORM,HNORM,DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,ORIENTVPAR,0)
            ENDDO

C           Upper endpoint add-back row.  Emit only the upper coefficient so
C           the two endpoint rows are not a duplicated complete add-back.
            GENDP = FUG*4.D0*SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK
            HENDX1 = FUX1*PHASE0
            HENDX2 = FUX2*PHASE0
            HENDQ1 = FUQ1*PHASE0
            HENDQ2 = FUQ2*PHASE0
            HENDQ3 = FUQ3*PHASE0
            HENDDP = FUDP*PHASE0
            IF (RJBK(NCHI2+2).NE.0.D0) THEN
               HCHIFACTOR = DPSIS/RJBK(NCHI2+2)
            ELSE
               HCHIFACTOR = 0.D0
            ENDIF
            VPSTATE = 0.D0
            ORIENTSTATE = 0.D0
            ORIENTVPAR = 0.D0
            CALL WRITEKJPCOMMONROW(FID,JS,JS_MAT,KGRID,0,K,NCHI2+2,
     &       RLAM,RM(K,2),RLM(L),RCHIK(NCHI2+2),
     &       RPHIK(NCHI2+2),RTK(NCHI2+2),1.D0,
     &       OMEGAB*RTK(NCHI2+2),RADIAL,RHK(NCHI2+2),B0K/RHK(NCHI2+2),
     &       RJBK(NCHI2+2),ZERO,ZERO,GENDP,GENDP,ZERO,HENDX1,HENDX2,
     &       HENDQ1,HENDQ2,HENDQ3,HENDDP,GNORM,HNORM,DPSIS,HCHIFACTOR,
     &       VPSTATE,ORIENTSTATE,ORIENTVPAR,1)
         ENDDO
      ENDDO
      CLOSE(FID)
C$OMP END CRITICAL(ELL_TRACE_WRITE)

      RETURN
      END SUBROUTINE WRITEKJPCOMMONTRACE

C=======================================================================
C One row writer keeps the schema and complex-pair ordering in one place.
C Six integer fields precede twelve scalar fields and eleven complex pairs;
C the two normalization scalars and five orientation fields make forty-one
C real fields in total.
C=======================================================================
      SUBROUTINE WRITEKJPCOMMONROW(FID,JS,JS_MAT,KGRID,KPARTICLE,
     & MIDX,JIDX,RLAM,MVAL,ELL,CHI,PHI,TAU,TAUFRACTION,BOUNCE_ANGLE,
     & RADIAL,RHVAL,BVAL,JBVAL,PG,GPARA,GPERP,GDPHI,PH,HX1,HX2,HQ1,HQ2,
     & HQ3,HDP,
     & GNORM,HNORM,DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,ORIENTVPAR,
     & ENDPOINT)

      IMPLICIT NONE
      INTEGER FID,JS,JS_MAT,KGRID,KPARTICLE,MIDX,JIDX,ENDPOINT
      REAL*8 RLAM,MVAL,ELL,CHI,PHI,TAU,TAUFRACTION,BOUNCE_ANGLE,
     &       RADIAL,RHVAL,BVAL,JBVAL,GNORM,HNORM,DPSIS,HCHIFACTOR,
     &       VPSTATE,ORIENTSTATE,ORIENTVPAR
      COMPLEX*16 PG,GPARA,GPERP,GDPHI,PH,HX1,HX2,HQ1,HQ2,HQ3,HDP
      WRITE(FID,1000) JS,JS_MAT,KGRID,KPARTICLE,MIDX,JIDX,
     & RLAM,MVAL,ELL,CHI,PHI,TAU,TAUFRACTION,BOUNCE_ANGLE,RADIAL,RHVAL,BVAL,
     & JBVAL,REAL(PG),AIMAG(PG),REAL(GPARA),AIMAG(GPARA),
     & REAL(GPERP),AIMAG(GPERP),REAL(GDPHI),AIMAG(GDPHI),REAL(PH),
     & AIMAG(PH),REAL(HX1),AIMAG(HX1),REAL(HX2),AIMAG(HX2),REAL(HQ1),
     & AIMAG(HQ1),REAL(HQ2),AIMAG(HQ2),REAL(HQ3),AIMAG(HQ3),REAL(HDP),
     & AIMAG(HDP),GNORM,HNORM,DPSIS,HCHIFACTOR,VPSTATE,
     & ORIENTSTATE,ORIENTVPAR,ENDPOINT
 1000 FORMAT(6I8,41(1X,E24.16),1X,I3)
      END SUBROUTINE WRITEKJPCOMMONROW

C=======================================================================
C DEFAULT-OFF LEG-RESOLVED KG/KH SOURCE PACKET.
C
C KG and KH currently fold the two trapped velocity legs into cos(theta).
C This packet writes the two pre-fold complex legs at the same native
C quadrature cells, including the analytic endpoint add-backs.  The row
C carries the source surface/orbit identifiers and the native K/M/ell
C indices; no correspondence, phase, radial, or field quantity is inferred.
C It is intentionally upstream of KJPFILL and the DWK cache.  The existing
C cache schema has already summed lambda/cell/leg contributions, so this
C producer packet is the exact first source boundary for a future cache
C extension; it does not pretend to provide a downstream consumer.
C=======================================================================
      SUBROUTINE WRITEKJPLEGTRACE(JS,JS_MAT,KGRID,RLAM)

      USE RCOMDM
      USE DIMENSIM
      USE GLOBALM
      USE KINETICM
      USE ToolBox
      IMPLICIT NONE
      INCLUDE 'compam.inc'

      INTEGER JS,JS_MAT,KGRID,J,K,L,M,FID,LEG,ENDPOINT
      INTEGER ORBITID,CELLID
      REAL*8 RLAM,OMEGAE,DPSIS,RADIAL,GNORM,HNORM,PHASE0,
     &       ARG,SQRTARG,THETA,PGPHASE,PHPHASE,CTMPL,CTMPU,
     &       GP1,GP2,GP3,HP1,HP6,HP2,HP3,HP4,HP5,HP7,
     &       HCHIFACTOR,VPSTATE,ORIENTSTATE,ORIENTVPAR
      COMPLEX*16 FLG,FUG,FLX1,FLX2,FLQ1,FLQ2,FLQ3,FLDP,
     &       FUX1,FUX2,FUQ1,FUQ2,FUQ3,FUDP,GLEG,HLEG,
     &       GPARA,GPERP,GDPHI,HX1,HX2,HQ1,HQ2,HQ3,HDP,ZERO,CALPHA
      LOGICAL OTRACE,OEXIST
      CHARACTER*128 PATH
      REAL*8 DIFFERCHI,DIFPI

      IF (KGRID.NE.1.AND.KGRID.NE.2) RETURN
      CALL KELLTRACESELECT(JS,KGRID,OTRACE)
      IF (.NOT.OTRACE) RETURN
      IF (IFOWT.NE.0) RETURN
      IF (.NOT.ALLOCATED(RCHIK).OR..NOT.ALLOCATED(RPHIK).OR.
     &    .NOT.ALLOCATED(RTK).OR..NOT.ALLOCATED(RHK).OR.
     &    .NOT.ALLOCATED(RJBK)) RETURN
      IF (RLAM.LE.0.D0) RETURN

      IF (KGRID.EQ.1) THEN
         OMEGAE=ROT(JS)
         DPSIS=DPSIDS(JS)
         RADIAL=CS(JS)
      ELSE
         OMEGAE=ROTM(JS)
         DPSIS=DPSIDSM(JS)
         RADIAL=CSM(JS)
      ENDIF
      IF (DPSIS.EQ.0.D0) RETURN
      IF (RTK(NCHI2+2).LE.0.D0) RETURN

      GNORM=RCHIHK/4.D0/PI
      HNORM=RCHIHK/4.D0*OMEGAB/PI
      PHASE0=4.D0*SQRT(DIFFERCHI(CHIU,CHIL))/RCHIHK
      ZERO=DCMPLX(0.D0,0.D0)
      CALPHA=CI/(OMEGA-RNTOR*OMEGAE)
      IF (IPERTURB.NE.0) CALPHA=ZERO
      IF (V2XKEY.EQ.1 .OR. V2XKEY.EQ.3) CALPHA=ZERO
      ORBITID=2*JS+KGRID

      WRITE(PATH,'("ELL_TRACE_JS",I4.4,"_G",I1,"_LEGS.OUT")')
     &      JS,KGRID
C$OMP CRITICAL(ELL_TRACE_WRITE)
      INQUIRE(FILE=PATH,EXIST=OEXIST)
      FID=ASSIGNFREEFILEUNIT()
      OPEN(FID,FILE=PATH,STATUS='UNKNOWN',POSITION='APPEND',
     &     ACTION='WRITE')
      IF (.NOT.OEXIST) THEN
         WRITE(FID,'(A)') '# schema: iter-tc24-mars-leg-resolved-kg-kh-v1'
         WRITE(FID,'(A)') '# source: KG/KH trapped KPARTICLE=0 before cosine fold'
         WRITE(FID,'(A)') '# orbit_id = 2*js+kgrid; lambda remains source-bound'
         WRITE(FID,'(A)') '# cell_id = native RCHIK/RTK sample index'
         WRITE(FID,'(A)') '# leg = -1, +1; sum of both legs restores native cosine terms'
         WRITE(FID,'(A)') '# endpoint = -1 lower add-back, 0 interior, +1 upper add-back'
         WRITE(FID,'(A)') '# downstream_cache = not yet carried; KJPFILL is first boundary'
         WRITE(FID,'(A)') '# fields: js js_mat kgrid kparticle orbit_id cell_id k m sample_index leg endpoint '
     &      //'lambda k_value m_value ell_value chi phi tau tau_fraction theta rho_pol '
     &      //'rho_h b0_over_b b_norm jb g_norm h_norm dpsids hchi vpar_state '
     &      //'orientation_state orientation_vpar gpara_re gpara_im gperp_re '
     &      //'gperp_im gdphi_re gdphi_im hx1_re hx1_im hx2_re hx2_im hq1_re '
     &      //'hq1_im hq2_re hq2_im hq3_re hq3_im hdp_re hdp_im'
      ENDIF

      DO L=1,MLMAX
         IF (NINT(RLM(L)).NE.KNTVELL) CYCLE
         DO K=1,MSMAX
         DO M=1,MSMAX
            FLG=RJBK(1)/2.D0*SQRT(RLAM)/SQRT(HPL)*
     &          EXP(-CI*RM(K,2)*CHIL)
            FUG=RJBK(NCHI2+2)/2.D0*SQRT(RLAM)/SQRT(-HPU)*
     &          EXP(CI*(-RNTOR*RPHIK(NCHI2+2)-RM(K,2)*CHIU))*
     &          COS(RLM(L)*PI)
            JLOOP: DO J=1,NCHI2+2
               IF (J.EQ.1) THEN
                  ENDPOINT=-1
               ELSEIF (J.EQ.NCHI2+2) THEN
                  ENDPOINT=1
               ELSE
                  ENDPOINT=0
               ENDIF
               CELLID=J
               IF (RJBK(J).NE.0.D0) THEN
                  HCHIFACTOR=DPSIS/RJBK(J)
               ELSE
                  HCHIFACTOR=0.D0
               ENDIF
               VPSTATE=0.D0
               ORIENTSTATE=0.D0
               ORIENTVPAR=0.D0

               IF (ENDPOINT.NE.0) THEN
                  THETA=RLM(L)*OMEGAB*RTK(J)
                  LEGLOOP: DO LEG=-1,1,2
                  IF (ENDPOINT.LT.0) THEN
                     GPARA=ZERO
                     GPERP=FLG*4.D0*SQRT(DIFFERCHI(CHIU,CHIL))/
     &                  RCHIHK/2.D0
                     GDPHI=GPERP
                     HX1=(RX1BK(1)-RX1RK(1)*CALPHA)*SQRT(RLAM)*
     &                  EXP(CI*RM(M,2)*CHIL)*2.D0/SQRT(HPL)*PHASE0/2.D0
                     HX2=RX2K(1)*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &                  2.D0/SQRT(HPL)*PHASE0/2.D0
                     HQ1=RQ1K(1)*RLAM*SQRT(RLAM)*
     &                  EXP(CI*RM(M,2)*CHIL)*2.D0/SQRT(HPL)*PHASE0/2.D0
                     HQ2=RQ2K(1)*RLAM*SQRT(RLAM)*
     &                  EXP(CI*RM(M,2)*CHIL)*2.D0/SQRT(HPL)*PHASE0/2.D0
                     HQ3=RQ3K*RLAM*SQRT(RLAM)*
     &                  EXP(CI*RM(M,2)*CHIL)*2.D0/SQRT(HPL)*PHASE0/2.D0
                     HDP=RJBK(1)/DPSIS*SQRT(RLAM)*
     &                  EXP(CI*RM(M,2)*CHIL)*2.D0/SQRT(HPL)*PHASE0/2.D0
                  ELSE
                     GPARA=ZERO
                     GPERP=FUG*4.D0*SQRT(DIFFERCHI(CHIU,CHIL))/
     &                  RCHIHK/2.D0
                     GDPHI=GPERP
                     HX1=(RX1BK(NCHI2+2)-RX1RK(NCHI2+2)*CALPHA)*SQRT(RLAM)*
     &                  EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &                  2.D0*COS(RLM(L)*PI)/SQRT(-HPU)*PHASE0/2.D0
                     HX2=RX2K(NCHI2+2)*SQRT(RLAM)*
     &                  EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &                  2.D0*COS(RLM(L)*PI)/SQRT(-HPU)*PHASE0/2.D0
                     HQ1=RQ1K(NCHI2+2)*RLAM*SQRT(RLAM)*
     &                  EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &                  2.D0*COS(RLM(L)*PI)/SQRT(-HPU)*PHASE0/2.D0
                     HQ2=RQ2K(NCHI2+2)*RLAM*SQRT(RLAM)*
     &                  EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &                  2.D0*COS(RLM(L)*PI)/SQRT(-HPU)*PHASE0/2.D0
                     HQ3=RQ3K*RLAM*SQRT(RLAM)*
     &                  EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &                  2.D0*COS(RLM(L)*PI)/SQRT(-HPU)*PHASE0/2.D0
                     HDP=RJBK(NCHI2+2)/DPSIS*SQRT(RLAM)*
     &                  EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &                  2.D0*COS(RLM(L)*PI)/SQRT(-HPU)*PHASE0/2.D0
                  ENDIF
                  CALL WRITEKJPLEGROW(FID,JS,JS_MAT,KGRID,ORBITID,CELLID,
     &               K,M,RM(K,2),RM(M,2),RLM(L),RLAM,RCHIK(J),RPHIK(J),RTK(J),
     &               RTK(J)/RTK(NCHI2+2),THETA,J,LEG,ENDPOINT,GPARA,GPERP,
     &               GDPHI,HX1,HX2,
     &               HQ1,HQ2,HQ3,HDP,GNORM,HNORM,RADIAL,RHK(J),B0K/RHK(J),
     &               RJBK(J),DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,
     &               ORIENTVPAR)
                  ENDDO LEGLOOP
                  CYCLE JLOOP
               ENDIF

               ARG=1.D0-RLAM/RHK(J)
               IF (ARG.LE.0.D0) CYCLE JLOOP
               SQRTARG=SQRT(ARG)
               THETA=RLM(L)*OMEGAB*RTK(J)
               PGPHASE=-RNTOR*RPHIK(J)-RM(K,2)*RCHIK(J)
               PHPHASE=RNTOR*RPHIK(J)+RM(M,2)*RCHIK(J)
               GP1=RJBK(J)*SQRTARG
               GP2=0.5D0*RJBK(J)*RLAM/RHK(J)/SQRTARG
               GP3=0.5D0*RJBK(J)/SQRTARG
               CTMPL=1.D0/SQRT(DIFPI(RCHIK(J)-CHIL))
               CTMPU=1.D0/SQRT(DIFPI(CHIU-RCHIK(J)))
               HP1=((1.D0-RLAM/RHK(J))*RX1PK(J)+
     &              (2.D0-RLAM/RHK(J))*RX1BK(J))/SQRTARG
               HP6=(2.D0-RLAM/RHK(J))*RX1RK(J)/SQRTARG
               HP2=(2.D0-RLAM/RHK(J))*RX2K(J)/SQRTARG
               HP3=RLAM*RQ1K(J)/SQRTARG
               HP4=RLAM*RQ2K(J)/SQRTARG
               HP5=RLAM*RQ3K/SQRTARG
               HP7=RJBK(J)/DPSIS/SQRTARG
               FLX1=(RX1BK(1)-RX1RK(1)*CALPHA)*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &              2.D0/SQRT(HPL)
               FLX2=RX2K(1)*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &              2.D0/SQRT(HPL)
               FLQ1=RQ1K(1)*RLAM*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &              2.D0/SQRT(HPL)
               FLQ2=RQ2K(1)*RLAM*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &              2.D0/SQRT(HPL)
               FLQ3=RQ3K*RLAM*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &              2.D0/SQRT(HPL)
               FLDP=RJBK(1)/DPSIS*SQRT(RLAM)*EXP(CI*RM(M,2)*CHIL)*
     &              2.D0/SQRT(HPL)
               FUX1=RX1BK(NCHI2+2)*SQRT(RLAM)*
     &              EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &              2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
               FUX2=RX2K(NCHI2+2)*SQRT(RLAM)*
     &              EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &              2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
               FUQ1=RQ1K(NCHI2+2)*RLAM*SQRT(RLAM)*
     &              EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &              2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
               FUQ2=RQ2K(NCHI2+2)*RLAM*SQRT(RLAM)*
     &              EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &              2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
               FUQ3=RQ3K*RLAM*SQRT(RLAM)*
     &              EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &              2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
               FUDP=RJBK(NCHI2+2)/DPSIS*SQRT(RLAM)*
     &              EXP(CI*(RM(M,2)*CHIU+RNTOR*RPHIK(NCHI2+2)))*
     &              2.D0*COS(RLM(L)*PI)/SQRT(-HPU)
               DO LEG=-1,1,2
                  GLEG=EXP(CI*(PGPHASE+LEG*THETA))/2.D0
                  HLEG=EXP(CI*(PHPHASE+LEG*THETA))
                  GPARA=GP1*GLEG
                  GPERP=GP2*GLEG-
     &               (CTMPL*FLG+CTMPU*FUG)/2.D0
                  GDPHI=GP3*GLEG-
     &               (CTMPL*FLG+CTMPU*FUG)/2.D0
                  HX1=HLEG*(HP1-HP6*CALPHA)-
     &               (CTMPL*FLX1+CTMPU*FUX1)/2.D0
                  HX2=HLEG*HP2-(CTMPL*FLX2+CTMPU*FUX2)/2.D0
                  HQ1=HLEG*HP3-(CTMPL*FLQ1+CTMPU*FUQ1)/2.D0
                  HQ2=HLEG*HP4-(CTMPL*FLQ2+CTMPU*FUQ2)/2.D0
                  HQ3=HLEG*HP5-(CTMPL*FLQ3+CTMPU*FUQ3)/2.D0
                  HDP=HLEG*HP7-(CTMPL*FLDP+CTMPU*FUDP)/2.D0
                  VPSTATE=SQRTARG
                  IF (HCHIFACTOR.EQ.0.D0) THEN
                     ORIENTSTATE=0.D0
                     ORIENTVPAR=0.D0
                  ELSE
                     ORIENTSTATE=1.D0
                     ORIENTVPAR=SIGN(1.D0,HCHIFACTOR)
                  ENDIF
                  CALL WRITEKJPLEGROW(FID,JS,JS_MAT,KGRID,ORBITID,CELLID,
     &               K,M,RM(K,2),RM(M,2),RLM(L),RLAM,RCHIK(J),RPHIK(J),RTK(J),
     &               RTK(J)/RTK(NCHI2+2),THETA,J,LEG,ENDPOINT,GPARA,GPERP,
     &               GDPHI,HX1,HX2,HQ1,HQ2,HQ3,HDP,GNORM,HNORM,RADIAL,
     &               RHK(J),B0K/RHK(J),RJBK(J),DPSIS,HCHIFACTOR,VPSTATE,
     &               ORIENTSTATE,ORIENTVPAR)
               ENDDO
            ENDDO JLOOP
         ENDDO
         ENDDO
      ENDDO
      CLOSE(FID)
C$OMP END CRITICAL(ELL_TRACE_WRITE)
      RETURN
      END SUBROUTINE WRITEKJPLEGTRACE

      SUBROUTINE WRITEKJPLEGROW(FID,JS,JS_MAT,KGRID,ORBITID,CELLID,K,M,KVAL,
     & MVAL,ELL,RLAM,CHI,PHI,TAU,TAUFRACTION,THETA,J,LEG,ENDPOINT,GPARA,GPERP,
     & GDPHI,HX1,HX2,HQ1,HQ2,HQ3,HDP,GNORM,HNORM,RADIAL,RHVAL,BVAL,JBVAL,
     & DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,ORIENTVPAR)
      IMPLICIT NONE
      INTEGER FID,JS,JS_MAT,KGRID,ORBITID,CELLID,K,M,J,LEG,ENDPOINT
      REAL*8 KVAL,MVAL,ELL,RLAM,CHI,PHI,TAU,TAUFRACTION,THETA,GNORM,HNORM,
     &       RADIAL,RHVAL,BVAL,JBVAL,DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,
     &       ORIENTVPAR
      COMPLEX*16 GPARA,GPERP,GDPHI,HX1,HX2,HQ1,HQ2,HQ3,HDP
      WRITE(FID,1000) JS,JS_MAT,KGRID,0,ORBITID,CELLID,K,M,J,LEG,ENDPOINT,
     & RLAM,KVAL,MVAL,ELL,CHI,PHI,TAU,TAUFRACTION,THETA,RADIAL,RHVAL,BVAL,JBVAL,
     & GNORM,HNORM,DPSIS,HCHIFACTOR,VPSTATE,ORIENTSTATE,ORIENTVPAR,
     & REAL(GPARA),AIMAG(GPARA),REAL(GPERP),AIMAG(GPERP),REAL(GDPHI),
     & AIMAG(GDPHI),REAL(HX1),AIMAG(HX1),REAL(HX2),AIMAG(HX2),REAL(HQ1),
     & AIMAG(HQ1),REAL(HQ2),AIMAG(HQ2),REAL(HQ3),AIMAG(HQ3),REAL(HDP),
     & AIMAG(HDP)
 1000 FORMAT(11I8,38(1X,E24.16))
      END SUBROUTINE WRITEKJPLEGROW
