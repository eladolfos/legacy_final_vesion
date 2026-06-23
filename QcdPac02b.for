      FUNCTION ALAMBD(AMU)
C                                                   -=-=- alambd
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Qcdpar
C Description: Many callable functions to return QCD parameters
C ListOfFiles: alambd nfl alamf amass amhatf nfltot amumin ch 
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Qcdpar.f,v 1.1 97/12/21 20:34:52 wkt Exp $ 
C #Log:	Qcdpar.f,v $
c Revision 1.1  97/12/21  20:34:52  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      ALAMBD = ALAMF(NFL(AMU))

      RETURN
      END
C
C***************************************************************
C
      FUNCTION ALAMF(N)
C                                                   -=-=- alamf
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                 Returns the value of LambdaCWZ in the energy range
C                 with N "light" quarks.
 
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / IOUNIT / NIN, NOUT, NWRT
      LOGICAL SET
C
      IF (.NOT.SET) CALL LAMCWZ
      IF ((N.LT.0) .OR. (N.GT.9)) THEN
         WRITE (NOUT, *) ' N IS OUT OF RANGE IN ALAMF'
         ALAMF=0.
      ELSE
         ALAMF = ALAM(MAX(N, NF-NHQ))
      ENDIF
      RETURN
      END
C
C***************************************************************
C
      FUNCTION ALEPI (AMU, NEF)
C                                                   -=-=- alepi
C                   Returns ALPHA/PI using the Effective Lamda appropriate for
C                             NEF flavors without regard to the value of AMU.
C                   Appropriate for Renormalization Schemes with fixed NEF.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
C
      IF(.NOT.SET) CALL LAMCWZ
 
ccpy

	print*,' NORDER, NEF, AMU =',NORDER, NEF, AMU

 
      ALM = ALAM(NEF)
      ALEPI = ALPQCD (NORDER, NEF, AMU/ALM, IRT)
 
      IF     (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALEPI', 'MU', AMU,
     >             ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2, NWRT, 'ALPI > 3; Be aware!', 'ALEPI', ALEPI,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEF, ALM, AMU
      ENDIF
 
      RETURN
      END
C
C                      **************************************
C
      SUBROUTINE ALFSET (QS, ALFS)
C                                                   -=-=- alfset
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Setalf
C Description: Routines to set the value of alpha_s and lambda_qcd
C ListOfFiles: alfset setlam lamcwz setl1
C ====================================================================
 
C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Setalf.f,v 1.1 97/12/21 20:34:55 wkt Exp $ 
C #Log:	Setalf.f,v $
c Revision 1.1  97/12/21  20:34:55  wkt
c Initial revision
c 

C                                 Given the value of Alpha(strong), ALFS, at
C                                 the scale QS, this routine determines the
C                                 effective # of flavors and Effective Lamda
C                                 at QS appropriate for the current value of
C                                 NORDER; and then call SET1 to setup the
C                                 whole package for subsequent use.
C                                 Calculates Alpha_s at loop order
C                                 specified by current value of NORDER.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      EXTERNAL RTALF
      COMMON / RTALFC / ALFST, JORD, NEFF
 
      DATA ALAM, BLAM, ERR / 0.01, 10.0, 0.02 /
 
      QST   = QS
      ALFST = ALFS
      CALL PARQCD (2, 'ORDR', ORDR, IR1)
      JORD  = ORDR
 
      NEFF = NFL(QS)
 
      EFLLN  = QZBRNT (RTALF, ALAM, BLAM, ERR, IR2)
      EFFLAM = QS / EXP (EFLLN)
 
        print *, "Call ALFSET", NEFF, EFFLAM
      CALL SETL1 (NEFF, EFFLAM)
 
      END
 
C
C**************************************************************
C
      FUNCTION ALPHEM()
C                                                   -=-=- alphem
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns the value of electromagnetic interaction coupling.
      COMMON /COMALP/A
      ALPHEM=A
      RETURN
C                      ******************
      END
C
C                                                          =-=-= Setalf
      FUNCTION ALPI (AMU)
C                                                   -=-=- alpi
C  These comments are enclosed in the lead subprogram to survive forsplit.

C ====================================================================
C GroupName: Alphas
C Description: Various callable functions for alpha_s and alpha_em
C ListOfFiles: alpi alpior alepi alpqcd g alphem
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Alphas.f,v 1.3 98/08/16 10:33:20 wkt Exp $ 
C #Log:	Alphas.f,v $
c Revision 1.3  98/08/16  10:33:20  wkt
c Warning in AlpQcd suppressed (cf. comments); numbers corrected.
c 
c Revision 1.2  98/08/11  21:39:31  wkt
c cross-line string argument corrected.
c 
c Revision 1.1  97/12/21  20:34:17  wkt
c Initial revision
c 

C               Returns effective g**2/(4pi**2) = alpha/pi.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
C                      Use the following as subroutine argument of type
C                      set by IMPLICIT statement:
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /

CCPY March 2014
      Common
     > / PdfSwh / Iset, IpdMod, Iptn0, NuIni
CCPY
	real*8 CT14Alphas
	External CT14Alphas
CJI add lhapdf
      integer initlhapdf, iset_save
      common/lhapdf/initlhapdf, iset_save

      if(initlhapdf == 1) then
          pi = 4.0d0*atan(1.0d0)
          call lhapdf_alphasq(1, iset_save, amu, alpi)
          alpi = alpi/pi
          return
      elseIF (Iset>=11900 .and. Iset<=11952) then !ZL
CCPY March 2014
        PI=4.0D0*ATAN(1.0D0)
        ALPI=CT14Alphas(AMU)/PI 
        RETURN 
      ENDIF

      IF(.NOT.SET) CALL LAMCWZ
 
      NEFF = NFL(AMU)
      ALM  = ALAM(NEFF)
      ALPI = ALPQCD (NORDER, NEFF, AMU/ALM, IRT)
CZLdebug        print *, "here", ALPI*PI, NORDER, AMU, ALM
 
      IF (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALPI', 'MU', AMU,
     >              ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEFF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2, NWRT, 'ALPI > 3; Be aware!', 'ALPI', ALPI,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEFF, ALM, AMU
      ENDIF
 
      RETURN
      END
C
C************
C
      FUNCTION ALPIOR (AMU, NL)
C                                                   -=-=- alpior
C               Returns effective g**2/(4pi**2) = alpha/pi.
C               Use formula with NL loops for beta function.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
 
C                      Use the following as subroutine argument of type
C                      set by IMPLICIT statement:
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
 
      IF (.NOT.SET) CALL LAMCWZ
 
      NEFF = NFL(AMU)
      ALM  = ALAM(NEFF)
      ALPIOR = ALPQCD (NL, NEFF, AMU/ALM, IRT)
 
      IF (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALPIOR', 'MU', AMU,
     >              ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEFF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2,NWRT,'ALPIOR > 3; Be aware!','ALPIOR',ALPIOR,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEFF, ALM, AMU
      ENDIF
      END
C
C*****************************************************************
C
      FUNCTION ALPQCD (IRDR, NF, RML, IRT)
C                                                   -=-=- alpqcd
 
C                                 Returns the QCD alpha/pi for RML = MU / LAMDA
C                                 using the standard perturbative formula for
C                                 NF flavors and to IRDR th order in 1/LOG(RML)
 
C                                 Return Code:  IRT
C                                                0:   O.K.
C                                                1:   Mu < Lamda; returns 99.
C                                                2:   Alpha > 10 ; be careful!
C                                                3:   IRDR out of range
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
      PARAMETER (CG = 3.0, TR = 0.5, CF = 4.0/3.0)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA IW1, Iw2 / 0, 0 /
 
      IRT = 0
 
      IF (IRDR .LT. 1 .OR. IRDR .GT. 2) THEN
        WRITE(NOUT, *)
     >  'Order parameter out of range in ALPQCD; IRDR = ', IRDR
        IRT = 3
        STOP
      ENDIF
 
      B0 = (11.* CG  - 2.* NF) / 3.
      B1 = (34.* CG**2 - 10.* CG * NF - 6.* CF * NF) / 3.
      RM2 = RML ** 2

C           AlpQcd is used mainly as a mathematical function, for
C            inversion, as well as evaluation of the physical Alpi.
C           Warning should be deferred to the calling program.
 
      IF (RM2 .LE. 1.) THEN
         IRT = 1
C         CALL QWARN(IW1, NWRT,
C     >    'RM2 (=MU/LAMDA) < 1. not allowed in ALPQCD; Alpha->99.',
C     >    'RM2', RM2, D1, BIG, 1)
         ALPQCD = 99
         RETURN
      ENDIF
 
      ALN = LOG (RM2)
      AL = 4./ B0 / ALN
 
      IF (IRDR .GE. 2) AL = AL * (1.- B1 * LOG(ALN) / ALN / B0**2)
 
      IF (AL .GE. 3.) THEN
         IRT = 2
C         CALL QWARN(IW2, NWRT, 'ALPQCD > 3. in ALPQCD', 'ALPQCD', AL,
C     >              D0, D1, 1)
      ENDIF
 
      ALPQCD = AL
 
      RETURN
C                       *********************
      END
C
      FUNCTION AMASS(I)
C                                                   -=-=- amass
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                     Returns mass of parton I.
C
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMQMS/ VALMAS(9)
      LOGICAL SET
      
      AMASS = 0D0
      II = IABS(I)
      IF (II.GE.1 .and. II.LE.6)  THEN
          AMASS = VALMAS(II)
      ENDIF
      
      RETURN
      END
C
C***********************************************************
C
      FUNCTION AMHATF(I)
C                                                   -=-=- amhatf
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                Returns the boundary in mass scale between the regions
C                with I & (I-1) effective "light" quarks.
C
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / IOUNIT / NIN, NOUT, NWRT
      LOGICAL SET
C
      IF (.NOT.SET) CALL LAMCWZ
      IF ((I.LE.0).OR.(I.GT.9)) THEN
         WRITE (NOUT,*) 'I IS OUT OF RANGE IN AMHATF'
         AMHATF = 0
      ELSE
         AMHATF = AMHAT(I)
      ENDIF
      RETURN
      END
C
C********************************************************************
C
      FUNCTION AMUMIN()
C                                                   -=-=- amumin
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns the minimum mu allowed by perturbative theory.
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NFL, NORDER, SET
      LOGICAL SET
      IF (.NOT.SET) CALL LAMCWZ
      AMUMIN = AMN
      RETURN
      END
C
C************************************************************************
C
      FUNCTION ANOM(Q1,Q2,GARRAY,N)
C                                                   -=-=- anom
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns integral from Q1 to Q2 of
C                 (dmu/mu) * gam(g(mu))
C              where N terms in anomalous dim. gam(g) are used:
C                 gam(g) = sum(i=1 to N) of gamma(i)*((g**2/(4*pi**2))**i).
C                 GAMMA(i) = sum(j=1 to i) of gammnf(i,j)*(nlfl**(j-1))
C                   except GAMMA(1)=gammnf(1,1)+gammnf(1,2)*nlfl
C                 gammnf(i,j) = GARRAY(1+j+i*(i-1)/2)
C                   except gammnf(1,1) = GARRAY(1) & gammnf(1,2)= GARRAY(2)
C                       where nlfl is the number of "light" flavor.
       DIMENSION GAMMA(5),GARRAY(10)
       COMMON /QCDPAR/ AL, NF, NORDER, SET
       COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
       COMMON /IOUNIT/ NIN, NOUT, NWRT
       LOGICAL SET
C
      IF(.NOT.SET) CALL LAMCWZ
      ANOM=0.
      IF (N.LE.0.) RETURN
      IF ((Q1.LE.AMN).OR.(Q2.LE.AMN)) THEN
         WRITE (NOUT, *) 'Q1 OR/AND Q2 IS TOO SMALL IN ANOM'
         RETURN
      ENDIF
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
C
  90  T1= LOG(Q2)
      B1=FLOAT(33-2*(NF-NMMIN+1))/12.
      GAMMA(1)=GARRAY(1)+GARRAY(2)*(NF-NMMIN+1)
      DO 200 J=NMMIN+1,NMMAX+1
         T2=T1
         B1=B1+1./6.
         GAMMA(1)=GAMMA(1)-GARRAY(2)
         IF (J.EQ.(NMMAX+1)) THEN
            T1= LOG(Q1)
         ELSE
            T1= LOG(AMHAT(NF-J+1))
         ENDIF
         ANOM=ANOM+0.5*GAMMA(1)* LOG((T2- LOG(ALAM(NF+1-J)))/
     >        (T1- LOG(ALAM(NF+1-J))))/B1
 200  CONTINUE
      RETURN
      END
C
C********************
C
      FUNCTION CH(I)
C                                                   -=-=- ch
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                     Returns charge of parton I.
C                        I=0 is gluon.
C                        I>0 is quark.
C                        I<0 is antiquark.
C                     See  BLOCK DATA for code.
C
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMQCH/ VALCH(9)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      LOGICAL SET
      CH=0.
      IF (IABS(I).GT.NF) THEN
          WRITE (NOUT, *) 'I IS OUT OF RANGE IN FUNCTION CH'
          print *, I
      ELSEIF (I.NE.0)  THEN
          CH = VALCH(IABS(I))
          IF (I.LT.0)  CH = -CH
      ENDIF
      RETURN
      END
C
C**********************************************************************
C
C                                                          =-=-= Parqcd
      SUBROUTINE CNVL1 (IRDR, JRDR, NF, VLAM)
C                                                   -=-=- cnvl1
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Setaux
C Description: Auxilary functions for lambda and alpha conversions
C ListOfFiles: cnvl1 zcnvlm trnlam zbrlam rtalf
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Setaux.f,v 1.1 97/12/21 20:34:57 wkt Exp $ 
C #Log:	Setaux.f,v $
c Revision 1.1  97/12/21  20:34:57  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL ZCNVLM
C  Auxiliary routine for SETLAM.
C                  Given Lamda (NF) = VLAM at order IRDR, this subroutine
C                  finds the corresponding Lamda at order JRDR so that the
C                  resulting Alpha(Mu, NRDR) remains approximately the same
C                  within the range of Neff = NF.
 
      COMMON / LAMCNV / AMU, ULAM, NFL, IRD, JRD
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA ALM, BLM, ERR, AMUMIN / 0.001, 2.0, 0.02, 1.5 /
 
      IRD = IRDR
      JRD = JRDR
      ULAM = VLAM
 
      CALL PARQCD(2, 'NFL', ANF, IRT)
      NTL = NFLTOT()
      IF (NF .GT. NTL) THEN
         WRITE (NOUT, *) ' NF .GT. NTOTAL in CNVL1; set NF = NTOTAL'
         WRITE (NOUT, *) ' NF, NTOTAL = ', NF, NTL
         NF = NTL
      ENDIF
C                                                First match at NFth threshold
      NFL = NF
      AMU = AMHATF(NF)
      AMU = MAX (AMU, AMUMIN)
      VLM1 = QZBRNT (ZCNVLM, ALM, BLM, ERR, IR1)
C                                          Match again at the next threshold
      IF (NF .LT. NTL) THEN
        AMU = AMHATF(NF+1)
        AMU = MAX (AMU, AMUMIN)
        VLM2 = QZBRNT(ZCNVLM, ALM, BLM, ERR, IR2)
      ELSE
        VLM2 = VLM1
      ENDIF
C                              Take the average and return new value of VLAM
      VLAM = (VLM1 + VLM2) / 2
 
      RETURN
C                        ****************************
      END
 
      SUBROUTINE EVOLUF(FQ2,FQ1,Q2,Q1,N)
C                                                   -=-=- evoluf
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                        ----------------------------
C?????????????????????????TO CHECK ????????????????????????????
C                        ----------------------------
C            Returns moments of parton distribution function at Q2 from Q1.
C            by integrating
C                D(F(Q,N,I))/D(LN(Q**2))=-GAMA(I,J)*F(Q,N,J)
C            where N is rank of the moments.
C            GAMA is aN NF X NF  matrix.
C            F(I=1) corresponds to gluon.
C            F(I=2 to NF+1) corresponds to NF quarks.
C
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      LOGICAL SET
      DIMENSION FQ1(11),FQ2(11),U(11,11),UINV(11,11),
     >          GAMAD(11,11),TEMP(11,11)
      DATA U,UINV,GAMAD/363*0./
C
      IF(.NOT.SET) CALL LAMCWZ
      IF((Q1.LE.AMN).OR.(Q2.LE.AMN)) GOTO 300
      IF(N.LT.2) GOTO 600
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
C
      T2=2.* LOG(Q1)
      B1=FLOAT(33-2*(NF-NMMAX-1))/4.
      SUM=0.
      DO 10 I=2,N
         SUM=SUM+1./FLOAT(I)
 10   CONTINUE
      GAGG=.75-9./FLOAT((N+1)*(N+2))-9./FLOAT(N*(N-1))
     >        +FLOAT(NF-NMMAX-1)/2.+SUM*9.
      GAFG=3./FLOAT(N)-6./FLOAT((N+1)*(N+2))
      GAGF=1./FLOAT(N+1)+2./FLOAT(N*(N-1))
      GAFF=1.-2./FLOAT(N*(N+1))+4.*SUM
      FQ2(1)=FQ1(1)
      DO 20 I=2,NF+1
         DO 30 J=3,NF+1
            UINV(J,I)=1./FLOAT(NF)
  30     CONTINUE
         FQ2(I)=FQ1(I)
         UINV(I,I)=UINV(I,I)-1.
         U(I,I)=-1.
         U(2,I)=1.
         U(I,1)=1.
         U(I,2)=1.
 20   CONTINUE
C
      DO 100 K=NMMAX+1,NMMIN+1,-1
         T1=T2
         B1=B1-.5
         GAGG=GAGG+.5
         IF(K.EQ.(NMMIN+1)) THEN
            T2=2.* LOG(Q2)
         ELSE
            T2=2.* LOG(AMHAT(NF-K+2))
         ENDIF
         A=SQRT((GAGG-GAFF)**2+4.*NF*GAGF*GAFG)
         U(1,1)=.5*(GAGG-GAFF+A)
         U(1,2)=.5*(GAGG-GAFF-A)
         DO 40 L=2,NF+1
            UINV(1,L)=-U(1,2)/A/FLOAT(NF)
            UINV(2,L)=U(1,1)/A/FLOAT(NF)
  40     CONTINUE
         UINV(1,1)=1./A
         UINV(1,2)=-1./A
         D=2* LOG(ALAM(NF+1-K))
         C=(T2-D)/(T1-D)
         DO 50 M=3,NF+1
            GAMAD(M,M)=C**(-GAFF/B1)
  50     CONTINUE
         GAMAD(1,1)=C**(-(U(1,1)+GAFF)/B1)
         GAMAD(2,2)=C**(-(U(1,2)+GAFF)/B1)
         CALL MTMUL(NF+1,NF+1,NF+1,U,GAMAD,TEMP)
         CALL MTMUL(NF+1,NF+1,NF+1,TEMP,UINV,TEMP)
         !CALL MTMUL(NF+1,NF+1,1,TEMP,FQ2,FQ2)
 100  CONTINUE
      RETURN
C
300   WRITE (NOUT,400)
400   FORMAT('Q1 OR/AND Q2 IS TOO SMALL IN EVOLUF')
600   DO 500 I=1,11
500      FQ2(I)=0.
      RETURN
      END
C
C                             *************************
      FUNCTION G (AMU)
C                                                   -=-=- g
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C               Returns effective coupling.
      PARAMETER (PI = 3.1415927)
      G=2.*PI*SQRT(ALPI(AMU))
      RETURN
C                      ******************
      END
C
      SUBROUTINE LAMCWZ
C                                                   -=-=- lamcwz
C                       Set /CWZPRM/ from /QCDPAR/
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
        print *, "Call LAMCWZ", NF,AL
      CALL SETL1 (NF, AL)
      END
C
C***********************
C
      SUBROUTINE MTMUL(L,M,N,A,B,C)
C                                                   -=-=- mtmul
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Qcdmsc
C Description: Miscellaneous functions and subroutines for the pac.
C ListOfFiles: mtmul qwarn qzbrnt sud anom evoluf
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Qcdmsc.f,v 1.1 97/12/21 20:34:49 wkt Exp $ 
C #Log:	Qcdmsc.f,v $
c Revision 1.1  97/12/21  20:34:49  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  TO DO MATRIX MULTIPLICATION
      DIMENSION A(11,11),B(11,11),C(11,11)
      DO 20 I=1,L
      DO 20 K=1,N
      C(I,K)=0.
      DO 20 J=1,M
20       C(I,K)=C(I,K)+A(I,J)*B(J,K)
      RETURN
      END
C
C***********************************************************
C
      FUNCTION NAMQCD(NNAME)
C                                                   -=-=- namqcd
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                   Find integer code corresponding to NAME. If no
C                   match,then return zero.
C
      CHARACTER NNAME*(*), NAME*8
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
C              ==== MACHINE DEPENDENCE
C              ==== ASSUMES CHARACTER CODES FOR '0' TO '9' ARE CONSECUTIVE
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
C
C                Use local variable to avoid problems with short
C                passed argument:
      NAME = NNAME
      NAMQCD=0
      IF ( (NAME .EQ. 'ALAM') .OR. (NAME .EQ. 'LAMB') .OR.
     1        (NAME .EQ. 'LAM') .OR. (NAME .EQ. 'LAMBDA') )
     2             NAMQCD=1
      IF ( (NAME .EQ. 'NFL') .OR. (NAME(1:3) .EQ. '#FL') .OR.
     1        (NAME .EQ. '# FL') )
     2             NAMQCD=2
      DO 10 I=1, 9
         IF (NAME .EQ. 'M'//CHAR(I+IASC0))
     1             NAMQCD=I+2
10       CONTINUE
      DO 20 I= 0, NF
         IF (NAME .EQ. 'LAM'//CHAR(I+IASC0))
     1             NAMQCD=I+13
20       CONTINUE
      IF (NAME(:3).EQ.'ORD' .OR. NAME(:3).EQ.'NRD') NAMQCD = 24
      RETURN
      END
C
C***************************************************************
C
C                                                          =-=-= Qcdmsc
      FUNCTION NFL(AMU)
C                                                   -=-=- nfl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C           NFL returns the number of 'light' flavors.
C
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
CJI LHAPDF
      integer initlhapdf, iset_save
      common/lhapdf/ initlhapdf, iset_save

      if(initlhapdf == 1) then
          call getthreshold(4, amass4)
          call getthreshold(5, amass5)
          if(amu > amass5) then
              nfl = 5
              return
          elseif(amu > amass4) then
              nfl = 4
              return
          else
              nfl = 3
              return
          endif
      endif
C
      IF (.NOT. SET) CALL LAMCWZ
      NFL = NF - NHQ
      IF ((NFL .EQ. NF) .OR. (AMU .LE. AMN)) GOTO 20
      DO 10 I = NF - NHQ + 1, NF
         IF (AMU .GE. AMHAT(I)) THEN
            NFL = I
         ELSE
            GOTO 20
         ENDIF
10       CONTINUE
20    RETURN
      END
C
C***************************************************************
C
      FUNCTION NFLTOT()
C                                                   -=-=- nfltot
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C          Returns the total number of flavors.
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
      NFLTOT=NF
      RETURN
      END
C
C***********************************************************
C
      SUBROUTINE PARQCD(IACT,NAME,VALUE,IRET)
C                                                   -=-=- parqcd
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: Parqcd
C Description: Input-output routines for setting and reading QCD parameters
C ListOfFiles: parqcd qcdin qcdout qcdset qcdget namqcd
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Parqcd.f,v 1.1 97/12/21 20:34:45 wkt Exp $ 
C #Log:	Parqcd.f,v $
c Revision 1.1  97/12/21  20:34:45  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C             Actions: 0     type list of variables on unit VALUE.
C                      1     set variable with name NAME to VALUE, if
C                             it exists, else set IRET to 0.
C                      2     find value of variable. If it does not exist,
C                             set IRET to 0.
C                      3     request values of all parameters from terminal.
C                      4     type list of all values on unit VALUE
C
C             IRET =   0     variable not found.
C                      1     successful search
C                      2     variable found, but bad value.
C                      3     bad value for IACT.
C                      4     no variable search (i.e., IACT is 0,3,or 4).
C
C             NAME is assumed upper-case.
C             if necessary, VALUE is converted to integer by NINT(VALUE)
C
      INTEGER IACT,IRET
      CHARACTER*(*) NAME
C
        print *, "PARQCD", IACT,NAME,VALUE
      IRET=1
      IF (IACT.EQ.0) THEN
         WRITE (NINT(VALUE), *)  'LAM(BDA), NFL, ORD(ER), Mi, ',
     >               '(i in 1 to 9), LAMi (i in 1 to NFL)'
         IRET=4
      ELSEIF (IACT.EQ.1) THEN
         CALL QCDSET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.2) THEN
         CALL QCDGET (NAME,VALUE,IRET)
      ELSEIF (IACT.EQ.3) THEN
         CALL QCDIN
         IRET=4
      ELSEIF (IACT.EQ.4) THEN
         CALL QCDOUT(NINT(VALUE))
         IRET=4
      ELSE
         IRET=3
      ENDIF
 
      RETURN
      END
C
C*******************************************************************
C
      SUBROUTINE QCDGET(NAME,VALUE,IRET)
C                                                   -=-=- qcdget
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Sets VALUE to the value of variable named NAME.
C            IRET=  0          variable not found.
C                   1          success.
C
C            NAME is assumed to be an upper-case character variable.
C
      CHARACTER*(*) NAME
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / COMQMS / VALQMS(9)
      LOGICAL SET
      PARAMETER (PI=3.1415927, EULER=0.57721566)
C
      ICODE = NAMQCD(NAME)
      IRET = 1
      IF (ICODE .EQ. 1) THEN
         VALUE = AL
      ELSEIF (ICODE .EQ. 2) THEN
         VALUE = NF
      ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 12))  THEN
         VALUE = VALQMS(ICODE - 2)
      ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
         VALUE = ALAM(ICODE - 13)
      ELSEIF (ICODE .EQ. 24) THEN
         VALUE = NORDER
      ELSE
         IRET=0
      ENDIF
      END
C
C *****
C
      SUBROUTINE QCDIN
C                                                   -=-=- qcdin
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                 Inputs QCD parameters.
C
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      DIMENSION VALMAS(9)
C              ==== MACHINE DEPENDENCE
C              ==== ASSUMES CHARACTER CODES FOR '0' TO '9' ARE CONSECUTIVE
C              NOTE THAT ICHAR('0') IS NON-STANDARD FORTRAN (ACCORDING TO
C              MICROSOFT).
      CHARACTER ONECH*(1)
      ONECH = '0'
      IASC0 = ICHAR(ONECH)
C
      CALL QCDGET ('LAM',ALAM,IRET1)
      CALL QCDGET ('NFL',ANF,IRET2)
      CALL QCDGET ('ORDER',ORDER,IRET3)
      NF = NINT(ANF)
      NORDER = NINT(ORDER)
 1    WRITE (NOUT, *) 'LambdaMSBAR, # Flavors, loop order ?'
      READ (NIN,*, IOSTAT = IRET) ALAM, NF, NORDER
      ORDER = NORDER
      ANF = NF
      IF (IRET .LT. 0) GOTO 22
      IF (IRET .EQ. 0) THEN
         CALL QCDSET ('LAM',ALAM,IRET1)
         CALL QCDSET ('NFL',ANF,IRET2)
         CALL QCDSET ('ORDER',ORDER,IRET3)
         ENDIF
      IF ((IRET.NE.0) .OR. (IRET1.NE.1) .OR. (IRET2.NE.1)
     >     .OR. (IRET3.NE.1)) THEN
         WRITE (NOUT, *) 'Bad value(s), try again.'
         GOTO 1
         ENDIF
      DO 20, I = 1, NF
         CALL QCDGET('M'//CHAR(I+IASC0),VALMAS(I),IRET1)
 10      WRITE (NOUT, '(1X,A,I2,A)') 'Mass of Quark', I, '?'
         READ (NIN,*, IOSTAT=IRET) VALMAS(I)
         IF (IRET .LT. 0) GOTO 22
         IF (IRET .EQ. 0)
     >      CALL QCDSET('M'//CHAR(I+IASC0),VALMAS(I),IRET1)
         IF ((IRET .NE. 0) .OR. (IRET1 .NE. 1)) THEN
            WRITE (NOUT, *) 'Bad value, try again.'
            GOTO 10
            ENDIF
 20      CONTINUE
      RETURN
C
 22   WRITE (NOUT, *) 'END OF FILE ON INPUT'
      WRITE (NOUT, *)
      RETURN
      END
C
C *****
C
      SUBROUTINE QCDOUT(NOUT)
C                                                   -=-=- qcdout
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Prints out the values of parameters to unit NOUT.
C
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /COMQMS/ VALQMS(9)
      LOGICAL SET
C
      IF (.NOT. SET) CALL LAMCWZ
      WRITE (NOUT,110) AL, NF, NORDER
 110   FORMAT(
     1 ' Lambda (MSBAR) =',G13.5,', NFL (total # of Flavors) =',I3,
     2 ', Order (loops) =', I2)
      WRITE (NOUT,120) (I,VALQMS(I),I=1,NF)
 120   FORMAT (3(' M', I1, '=', G13.5, :, ','))
      IF (NHQ .GT. 0)
     1   WRITE (NOUT,130) (I, ALAMF(I), I = NF-NHQ, NF)
 130   FORMAT (: ' ! Effective lambda given number of light quarks:'/
     >    (2(' ! ', I1, ' quarks => lambda = ', G13.5 : '; ')) )
      RETURN
      END
C
C***************************************************************
C
      SUBROUTINE QCDSET (NAME,VALUE,IRET)
C                                                   -=-=- qcdset
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Assign the variable whose name is specified by NAME the value VALUE
C             IRET=  0      variable not found.
C                    1      success.
C                    2      variable found, but bad value.
C             NAME is assumed upper-case and VALUE is real.
C             If necessary, VALUE is converted to integer by NINT(VALUE).
C
      CHARACTER*(*) NAME
      COMMON / COMQMS / VALMAS(9)
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
      PARAMETER (PI=3.1415927, EULER=0.57721566)
C
      IVALUE = NINT(VALUE)
      ICODE  = NAMQCD(NAME)
                print *, "Call QCD SET", ICODE, VALUE
      IF (ICODE .EQ. 0) THEN
         IRET=0
      ELSE
         IRET = 1
         SET = .FALSE.
         IF (ICODE .EQ. 1) THEN
            IF (VALUE.LE.0) GOTO 12
            AL=VALUE
         ELSEIF (ICODE .EQ. 2) THEN
            IF ( (IVALUE .LT. 0) .OR. (IVALUE .GT. 9)) GOTO 12
            NF = IVALUE
         ELSEIF ((ICODE .GE. 3) .AND. (ICODE .LE. 11))  THEN
            IF (VALUE .LT. 0) GOTO 12
C                                  When the mass of a quark is changed,
C     we reset all values of lambda by holding alpha(min(oldmass,newmass))
C           the same. This is only a prescription. It preserves the value                           
C                                  of alpha below this threshold,                   
            Scle = Min (Value , VALMAS(ICODE - 2))
            AlfScle = Alpi(Scle) * Pi
            VALMAS(ICODE - 2) = VALUE
            Call AlfSet (Scle, AlfScle)
         ELSEIF ((ICODE .GE. 13) .AND. (ICODE .LE. 13+NF))  THEN
            IF (VALUE .LE. 0) GOTO 12
            CALL SETL1 (ICODE-13, VALUE)
         ELSEIF (ICODE .EQ. 24)  THEN
            IF ((IVALUE .LT. 1) .OR. (IVALUE .GT. 2)) GOTO 12
            NORDER = IVALUE
         ENDIF
         IF (.NOT. SET) CALL LAMCWZ
      ENDIF
      RETURN
C
C              Illegal value
 12   IRET=2
      RETURN
      END
C
C************************************************
C
      SUBROUTINE QWARN (IWRN, NWRT1, MSG, NMVAR, VARIAB,
C                                                   -=-=- qwarn
     >                  VMIN, VMAX, IACT)
 
C     Subroutine to handle warning messages.  Writes the (warning) message
C     and prints out the name and value of an offending variable to SYS$OUT
C     the first time, and to output file unit # NWRT1 in subsequent times.
C
C     The switch IACT decides whether the limits (VMIN, VMAX) are active or
C     not.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
 
      CHARACTER*(*) MSG, NMVAR

      Save Iw
 
      IW = IWRN
      VR = VARIAB
 
      If (Iw .LT. 100) Then
         WRITE (NWRT1,'(I5, 3X,A/ 1X,A,'' = '',1PD16.7)') IW, MSG,
     >                  NMVAR, VR
      Else
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' !! Error # > 100 !! ; better check file unit #', NWRT1
      EndIf   

      IF  (IW .EQ. 0) THEN
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' Complete set of warning messages on file unit #', NWRT1
         IF (IACT .EQ. 1) THEN
         WRITE (NOUT,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         WRITE (NWRT1,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         ENDIF
      ENDIF
 
      IWRN = IW + 1
 
      RETURN
C                         *************************
      END

      FUNCTION QZBRNT(FUNC, X1, X2, TOLIN, IRT)
C                                                   -=-=- qzbrnt
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)

      external func
 
      TOL = ABS(TOLIN)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        WRITE (NOUT, *) 'Root must be bracketed for QZBRNT.'
        IRT = 1
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          QZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE (NOUT, *) 'QZBRNT exceeding maximum iterations.'
      IRT = 2
      QZBRNT=B
      RETURN
C**************************************************
      END
C
      FUNCTION RTALF (EFLLN)
C                                                   -=-=- rtalf
C  Auxiliary function for ALFSET, which solves equation RTALF=0.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589)
      COMMON / RTALFC / ALFST, JORD, NEFF
 
      EFMULM = EXP (EFLLN)
      TEM1 = PI / ALFST
      TEM2 = 1. / ALPQCD (JORD, NEFF, EFMULM, I)
 
      RTALF = TEM1 - TEM2
C************************************************************ 
      END

C                                                          =-=-= Qcdpar
      SUBROUTINE SETL1  (NEF, VLAM)
C                                                   -=-=- setl1
C     Given LAMDA = VLAM for NEF flavors:
C                    (i) fills the array  ALAM (0:NF) with effective LAMDAs;
C                    (ii) fills the array AMHAT (NF) with threshold masses;
C                    (iii) count the # of "heavy quarks" (QMS > EFFLAM);
C                    (iv) fix the parameter AMN defined as MAX (ALAM),
C                         times safety factor;
C                    (v) set AL in / QCDPAR / equal to ALAM (NF);
C                    (vi) let SET = .TRUE.
C       Uses formula with NORDER (1 or 2) -- see /QCDPAR/
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      LOGICAL SET
 
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / COMQMS / QMS(9)
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      IF (NEF .LT. 0 .OR. NEF .GT. NF) THEN
        WRITE(NOUT,*)'NEF out of range in SETL1, NEF, NF =',NEF,NF
        STOP
      ENDIF
C             Mass Thresholds are given by the Quark masses in the CWZ scheme
      AMHAT(0) = 0.
      DO 5 N = 1, NF
         AMHAT(N) = QMS(N)
    5    CONTINUE
      ALAM(NEF) = VLAM
      DO 10 N = NEF, 1, -1
         CALL TRNLAM(NORDER, N, -1, IR1)
   10    CONTINUE
      DO 20 N = NEF, NF-1
         CALL TRNLAM(NORDER, N, 1, IR1)
   20    CONTINUE
C=========================                Find first light quark:
      DO 30, N = NF, 1, -1
         IF ((ALAM(N) .GE. 0.7 * AMHAT(N))
     >       .OR. (ALAM(N-1) .GE. 0.7 * AMHAT(N)))THEN
            NHQ = NF - N
            GOTO 40
            ENDIF
   30    CONTINUE
      NHQ = NF
   40 CONTINUE
      DO 50, N = NF-NHQ, 1, -1
         AMHAT(N) = 0
         ALAM(N-1) = ALAM(N)
   50    CONTINUE
C========================               Find minimum mu
      AMN = ALAM(NF)
      DO 60, N = 0, NF-1
         IF (ALAM(N) .GT. AMN)  AMN = ALAM(N)
   60    CONTINUE
      AMN = AMN * 1.0001
      AL = ALAM(NF)
      SET = .TRUE.
      RETURN
C**************************************************************
      END

C                                                          =-=-= Setaux
      SUBROUTINE SETLAM (NEF, WLAM, IRDR)
C                                                   -=-=- setlam
C     The values of LAMBDA=WLAM with NEF effective flavors is given. The
C     coupling is assumed to be given by the IRDR formula.
C     First lambda is converted to a value that gives approximately the
C     same value of alpha_s when the formula with the current value of
C     NORDER is used.  Then SETL1 is called to update the rest of the
C     internal arrays.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET

 
      IF ((NEF .LT. 0) .OR. (NEF .GT. NF)) THEN
         WRITE(NOUT,*)'NEF out of range in SETLAM, NEF, NF=', NEF,NF
         STOP
      ENDIF
C                         Adjust Lamda value if ORDER parameters donot match.
C                         NORDER is not changed.
      VLAM = WLAM
      IF (IRDR .NE. NORDER) CALL CNVL1 (IRDR, NORDER, NEF, VLAM)
        print *, "Call SETLAM", NEF, VLAM
      CALL SETL1 (NEF, VLAM)
      END
C
C************************************************************
C
      Subroutine SetQCD
C                                                   -=-=- setqcd
C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: Setqcd
C Description: Set up the qcdpac of programs, initiate common blocks
C ListOfFiles: setqcd
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Setqcd.f,v 1.1 97/12/21 20:35:00 wkt Exp $ 
C #Log:	Setqcd.f,v $
c Revision 1.1  97/12/21  20:35:00  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      External DatQCD

      Dummy = 0.

      Return
C                        ****************************
      END

      BLOCK DATA DATQCD
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /COMQCH/ VALQCH(9)
      COMMON /COMQMS/ VALQMS(9)
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMALP/ ALPHA
      LOGICAL SET
C
      DATA AL, NF, NORDER, SET / .226, 5, 2, .FALSE. /
      DATA VALQCH/ 0.66666667, -0.33333333,
     >  -0.33333333, 0.66666667,
     >  -0.33333333, 0.66666667,
     >  3*0./
      DATA VALQMS/  2*0.001, 0.2, 1.3, 4.5, 174., 3*0./
      DATA ALPHA/  7.29927E-3 /
 
C                       ******************************
      END

C                                                          =-=-= Alphas
      FUNCTION SUD(Q1,Q2,GARRAY,N)
C                                                   -=-=- sud
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C              Returns integral from Q1 to Q2 of
C                 (dmu/mu) * gam(g(mu)) * ln(Q2/mu)
C              where N terms in anomalous dim. gam(g) are used:
C                 gam(g) = sum(i=1 to N) of GAMMA(i)*((g**2/(4*pi**2))**i).
C                 GAMMA(i) = sum(j=1 t0 i) of gammnf(i,j)*(nlfl**(j-1))
C                 gammnf(i,j) = GARRAY(j+i*(i-1)/2)
C                        where nlfl is the number of "light" flavors.
C
       DIMENSION GAMMA(5),GARRAY(10)
       COMMON / QCDPAR / AL, NF, NORDER, SET
       COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
       COMMON / IOUNIT / NIN, NOUT, NWRT
       LOGICAL SET
       GAMMNF(K,L)=GARRAY(L+K*(K-1)/2)
C
      IF(.NOT.SET) CALL LAMCWZ
      SUD=0.
      IF(N.LE.0.) RETURN
      IF((Q1.LE.AMN).OR.(Q2.LE.AMN)) THEN
         WRITE(NOUT, *) 'Q1 OR/AND Q2 IS TOO SMALL IN SUD'
         RETURN
      ENDIF
C
      NMMIN=NF-NFL(Q2)
      NMMAX=NF-NFL(Q1)
      T10=2.* LOG(Q2)
      B1=FLOAT(33-2*(NF-NMMIN+1))/12.
      B2B1SQ=FLOAT(153-19*(NF-NMMIN+1))/(24.*B1*B1)
      GAMMA(2)=GAMMNF(2,1)+GAMMNF(2,2)*(NF-NMMIN+1)
      DO 200 J=NMMIN+1,NMMAX+1
         T20=T10
         B2B1SQ=(B2B1SQ*B1*B1+19./24.)/(B1+1./6.)**2
         B1=B1+1./6.
         IF (J.EQ.(NMMAX+1)) THEN
            T10=2.* LOG(Q1)
         ELSE
            T10=2.* LOG(AMHAT(NF-J+1))
         ENDIF
         ALLAM=2.* LOG(ALAM(NF+1-J))
         TQ2=2.* LOG(Q2)-ALLAM
         T1=T10-ALLAM
         T2=T20-ALLAM
         ALNT2= LOG(T2)
         ALNT1= LOG(T1)
         SUD=SUD+GARRAY(1)*0.25/B1*(TQ2*(ALNT2-ALNT1
     1             +B2B1SQ*((ALNT2+1.)/T2-(ALNT1+1.)/T1))
     2             +T1-T2+B2B1SQ*(ALNT2**2-ALNT1**2)/2.)
C
         IF (N.GE.2) THEN
            GAMMA(2)=GAMMA(2)-GAMMNF(2,2)
            SUD=SUD+0.25*GAMMA(2)/(B1*B1)*(ALNT1-ALNT2
     >                                +TQ2*(1./T1-1./T2))
         ENDIF
200   CONTINUE
      RETURN
      END
C
C**********************
C
      SUBROUTINE TRNLAM (IRDR, NF, IACT, IRT)
C                                                   -=-=- trnlam
 
C     This routine transforms LAMDA (N) to LAMDA (N+IACT) where IACT = 1/-1
C     The transformation is obtained by requiring the coupling constant to
C                be continuous at the scale Mu = Mass of the (N+1)th quark.
 
C                                         IRT is an return code.
C                                            (0 for OK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      EXTERNAL ZBRLAM
 
      DATA ALM0, BLM0, RERR / 0.01, 10.0, 0.0001 /
      DATA IR1, SML / 0, 1.E-5 /
 
      IRT = 0
 
      N = NF
      JRDR = IRDR
      JACT = IACT
      VLAM = ALAM(N)
 
      IF (JACT .GT. 0) THEN
         N1 = N + 1
         THMS = AMHAT(N1)
         ALM = LOG (THMS/VLAM)
         BLM = BLM0
      ELSE
         N1 = N -1
         THMS = AMHAT(N)
         ALM = ALM0
         THMS = MAX (THMS, SML)
         BLM = LOG (THMS/VLAM)
      ENDIF
C                          Fix up for light quark:
      IF (VLAM .GE. 0.7 * THMS) THEN
         IF (JACT . EQ. 1) THEN
            AMHAT(N1) = 0
         ELSE
            AMHAT(N) = 0
         ENDIF
         IRT = 4
         ALAM(N1) = VLAM
         RETURN
      ENDIF
 
C             QZBRNT is the root-finding function to solve ALPHA(N) = ALPHA(N1)
C             Since 1/Alpha is roughly linear in Log(Mu/Lamda), we use the
C             former in ZBRLAM and the latter as the function variable.
      IF (ALM .GE. BLM) THEN
         WRITE (NOUT, *) 'TRNLAM has ALM >= BLM: ', ALM, BLM
         WRITE (NOUT, *) 'I do not know how to continue'
         STOP
         ENDIF
      VMULM = THMS/VLAM
      ERR = RERR * LOG (VMULM)
      WLLN = QZBRNT (ZBRLAM, ALM, BLM, ERR, IR1)
      ALAM(N1) = THMS / EXP (WLLN)
 
      IF (IR1 .NE. 0) THEN
         WRITE (NOUT, *) 'QZBRNT failed to find VLAM in TRNLAM; ',
     >        'NF, VLAM =', NF, VLAM
         WRITE (NOUT, *) 'I do not know how to continue'
        STOP
      ENDIF
      RETURN
      END
C                             *************************
 
      FUNCTION ZBRLAM (WLLN)
C                                                   -=-=- zbrlam
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      WMULM = EXP (WLLN)
      TEM1 = 1./ ALPQCD(JRDR, N1, WMULM, I)
      TEM2 = 1./ ALPQCD(JRDR, N,  VMULM, I)
 
      ZBRLAM = TEM1 - TEM2
 
      END
 
 
C************************************************************
 
      FUNCTION ZCNVLM (VLAM)
C                                                   -=-=- zcnvlm
C  Auxiliary function for CNVL1, which solves ZCNVLAM=0.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / LAMCNV / AMU, ULAM, NFL, IRD, JRD
 
      ZCNVLM= ALPQCD (IRD,NFL,AMU/ULAM,I) - ALPQCD (JRD,NFL,AMU/VLAM,I)
 
      END
C
C**************************************************************
C
