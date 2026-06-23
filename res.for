CCPY Jan 2021: Changing the input card format (Legacy.in) and subroutine SetC1_C4
CCPY Oct 2019: Changing the value of NRES1 could change the result of calculation
CCPY May 2007: Add A3 for gg-->H process.  
C March 7, 2005: Add heavy quark amss effect in the CSS piece for 
c (b+bbar --> HB) process. The keyword is iHQMass.
C Error: March 1, 2005: Correct the error of C(1) for H+ and HB;
c    SET I_CXFF=1 FOR NLO CALCULATION 
C Feb 10, 2005, B(2) for bb-->HB process is added, and eliminate 
C "omega" in HB process which is not needed if the input
C B quark mass is the MS-bar running mass at M_B scale. 
C Error: April 14, 2004 (correct ALEPASY for A0, etc)
C Feb 2004, Add (INONPERT.EQ.5) 
C Jan 2004, AN ERROR IN B2 IS CORRECTED 
C Oct 2002: I adjust the values of 
C       ACC_SUDNPT=1.D-4
C       ACC_RERR=1.D-4
C AND CHANGE THE LINE
C       IF(DABS(RES_1(I)/RESUM_1).LT.ACC_SUDNPT) GOTO 210
C   
CCPY April 2002: correct the ERROR in B2, inside WSETUP:
C -----------------------------------------------------------------------
C  %%%%%%%%%%%%%%%%%           VERSION 5.2             %%%%%%%%%%%%%%%%
C -----------------------------------------------------------------------
CCPY  Oct. 14, 1995 -- Add virtual gluon "resonance"
CSM   Feb.  8, 1996 -- Add gg initial state
CCPY  May 1996 -- a major bug was found for (P N) process
C                 inside WCXFCXF(B_SUD,CFCF)
CsB   Mar. 27, 1997
C     This file and the resumed numerical results were tested against C.-P.'s
C     version. This file contains all information what is in C.-P.'s
C     and some additional comments. The pion initial state is also contained
C     here but not in C.-P.'s file.
C     The numerical results are identical
C     for ppBar, pp and pN processes
C     for W+/-, Z0, A0, AA, AG and H0 bosons
C     for Q = 80 GeV, y = 0 and 2.5 and QT = 2.4 and 10 GeV's
C     (for any of the above combinations).
C     Remains to be checked:
C     - the pion initial state process was not tested against another
C       code or the literature,
C     - W+/- production in pN collision is zero,
C     - agreement of GG and GL options with C.-P.
C
CsB   July 17, 1997
C     CFlag_C3 = 5 option is added.
C
CsB   Nov. 17, 1997 q Q -> Z Z process is added.
C

C -----------------------------------------------------------------------
      SUBROUTINE WSETUP(NF_EFF)
C -----------------------------------------------------------------------
C FOR W+, W-, Z0, A0 PRODUCTION
C DEFINE CONSTANT VARIABLES.
C UNIT:   GeV AND pb
C THIS ONLY WORKS FOR  C1/C2=2*EXP(-EULER)  
      IMPLICIT NONE
      INCLUDE 'common.for'

      INTEGER NF_EFF
      REAL*8 SCALING,H1,H2
      integer i

CJI Sept 2013: changed this BPIECE to be included in the legacy.in file
C     under the variable named N_SUD_B to match the desired form for the
C     header file used when running ResBos
C     For regular usage, bpiece=T
C      LOGICAL BPIECE
C      DATA BPIECE/.TRUE./
C      DATA BPIECE/.FALSE./


      DATA ZEROJ0 /0.0, 2.4048256, 5.5200781, 8.6537279, 11.7915344,
     >14.9309177, 18.0710640, 21.2116366, 24.3524715, 27.4934791,
     >30.6346065, 33.7758202, 36.9170984, 40.0584258, 43.1997917,
     >46.3411884, 49.4826099, 52.6240518, 55.7655108, 58.9069839,
     >62.0484692, 65.1899648, 68.3314693,
     >69.9022245, 73.6145006,
     >77.7560256, 80.8975559, 84.0390908, 87.1806298, 90.3221726/
!     >93.4637188, 96.6052680, 99.7468199, 102.888374, 106.029931,
!     >109.171490, 112.313050, 115.454613, 118.596177, 121.737742,
!     >124.879309, 128.020877, 131.162446, 134.304017, 137.445588,
!     >140.587160, 143.728734, 146.870308, 150.011883, 153.153458,
!     >156.295034/

      REAL*8 H1_C2
      real*8 BETA2
      real*8 A4c
CJI Jan 2021: Include canonical values
      real*8 A1c,A2c,A3c,B1c,B2c,B3c,C1_IN,C2_IN,C3_IN
      common /canonical/ a1c,a2c,a3c,b1c,b2c,b3c,C1_IN,C2_IN,C3_IN

CJI September 2016, add in the choice of scheme
      CHARACTER*3 SCHEME
      COMMON / RESUMSCHEME / SCHEME
      
C THE ACCURACY IN DOING B-SPACE INTEGRATIONS IN    SUDNPT
CCPY      ACC_SUDNPT=1.D-3
C      ACC_SUDNPT=1.D-4
      ACC_SUDNPT=1.D-5
      ACC_AERR=1.D-16
CCPY      ACC_RERR=1.D-3
CBY      
      ACC_RERR=1.D-4
      ACC_RES1=1.D-3
C SCALING IS USED TO SCALE AS LOG(ENERGY OF THE COLLIDER)
C WE ASSUME COMPARED WITH TEVATRON
      SCALING=(1.0/3.0)*DLOG(ECM)/DLOG(1800.0D0)

C FOR NONPERT
      YNPT_CUT=50.0*SCALING
CCPY March 2019: Add more cycles of ZEROJ0
C FOR RESUM
      NRES1=22
C      NRES1=25
C FOR PDF
C FOR TEVATRON  1.8 TEV
C      PREV_UNDER=1.D8
C FOR SSC  40 TEV
      PREV_UNDER=1.D0
C FOR ZJSUDNPT
      GOAROUND=.FALSE.
C      GOAROUND=.TRUE.
C FOR HMRSE, HMRSB PDF
C FIT3=TRUE FOR 3RD ORDER POLY. FIT. OTHERWISE, LINEAR FIT
      FIT3=.FALSE.
C
      DEBUG=.FALSE.
C      DEBUG=.TRUE.

C TO TEST THE EFFECTS DUE TO LEPASY=1 IN C-FUNCTIONS
      TESTASY=.FALSE.
CCPY      TESTASY=.TRUE.

C TO SPEED UP
C THIS HAS BEEN CHECKED FOR W AND Z AT TEVATRON.
C      ACC_RES1=0.05
c      YNPT_CUT=250.0*SCALING
c      NRES1=INT(7*SCALING)

      P00=1.0D0
      P02=-1.0D0*9.0D0/2.0D0/(8.0D0**2)
      P04=1.0D0*9.0D0*25.0D0*49.0D0/24.0D0/(8.0D0**4)
      Q01=-1.0D0/8.0D0
      Q03=1.0D0*9.0D0*25.0D0/6.0D0/(8.0D0**3)
      EULER=0.577215664901D0
      ZETA3=1.2020569031595942854D0

      NC=3
      CA=NC
      CF=(NC**2-1.0D0)/2.0D0/NC
      NF=NF_EFF
      TR=NF/2.0D0
CPN
      Call SetC1_C4    ! Set up the constants C1..C4

CSM -- moved BETA here -- no good reason
      BETA1=(11.0D0*CA-4.0D0*TR)/12.0D0
CZL -- add beta2 for scale dependent part of A3
CJI Switch the beta2 to be consistent with alpha_s/pi expansion as 
C   beta1 above is
      BETA2=(17*NC**2-5*NC*NF-3*CF*NF)/24.0d0
CSM -- no reason to redefine B0
C      B0=C1/C2
C
CsB >>>
C      A1=CF
C      A2=0.5D0*CF*((67.0D0/18.0D0-PI2/6.0D0)*CA-10.0D0/9.0D0*TR)
C      B1=   -3.0D0/2.0D0*CF
C      IF(BPIECE) THEN
C      	B2=0.5D0*CF**2*(PI2-3.0D0/4.0D0-12.0D0*ZETA3)+
C     >   0.5D0*CF*CA*(11.0D0*PI2/9.0D0-193.0D0/12.0D0+6.0D0*ZETA3)+
C     >   0.5D0*CF*TR*(-4.0D0*PI2/9.0D0+17.0D0/3.0D0)
C      ELSE
C        B2=0.0
C      ENDIF
      A1c = CF
      A1 = CF
CsB A2_here = A2_PRD
CJI Modify to match my notes
      A2c = 0.5D0*CF*((67.0D0/18.0D0-PI2/6.0D0)*CA-10.0D0/9.0D0*TR)
      A2 = A2c-2.d0*A1c*beta1*DLog(B0/C1)
CZL---------------------------------------------------
CZL copied from GG case in the following code
      IF(N_SUD_A.ge.3) THEN
        A3c= CF**2*NF*0.5D0*(ZETA3-55.0/48.0)-CF*NF**2/108.0+
     >      CA**2*CF*(11.0*ZETA3/24.0+11.0*pi**4/720.0-67.0*pi**2/216.0
     >                +245.0/96.0)+
     >      CA*CF*NF*(-7.0*ZETA3/12.0+5.0*pi**2/108.0-209.0/432.0)
!ZL only valid for certain combination of C1...C4
CJI Modify to match my notes
        A3=A3c
     &     -2*A2c*2*BETA1*log(B0/C1)
     &     +(2*BETA1)**2*A1*log(B0/C1)**2
     &     -A1*2.d0*BETA2*log(B0/C1)
CJI Add in A4 from 1808.08981
        if(n_sud_a.ge.4) then
        A4c = (20702 - 5171.9*nf + 195.5722*nf*nf +
     &         3.272344*nf*nf*nf)/4**4
        A4 = A4c
        endif
        ENDIF
CZL---------------------------------------------------
CsB B1_Here = B1_PRD
      B1c = CF*(-3.0D0/2.0D0 )
      B1 = B1c - A1c*2.d0*DLog(C2*B0/C1)
      IF(N_SUD_B.ge.2) THEN
CsB B2_here = 2*B2_PRD
        If( (TYPE_V.EQ.'W+' .or. TYPE_V.EQ.'W-'
     > .or. TYPE_V.EQ.'Z0' .or. TYPE_V.EQ.'A0')
     > .or. (Type_V.Eq.'W+_RA_UDB' .or. Type_V.Eq.'W-_RA_DUB'
     >     .or. Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
CJI 2013: Added in ZU and ZD to match c++ code
     >     .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') )then
CCPY Jan 2004, AN ERROR IN B2 IS CORRECTED FROM 
C    >    -((67.d0/18.d0-Pi/6.d0)*NC - 5.d0/9.d0*NF)*DLog(C2*b0/C1) +
C TO
C    >    -((67.d0/18.d0-Pi2/6.d0)*NC - 5.d0/9.d0*NF)*DLog(C2*b0/C1) +

         B2c=0.5D0*CF**2*(PI2-3.0D0/4.0D0-12.0D0*ZETA3)+
     >   0.5D0*CF*CA*(11.0D0*PI2/9.0D0-193.0D0/12.0D0+6.0D0*ZETA3)+
     >   0.5D0*CF*TR*(-4.0D0*PI2/9.0D0+17.0D0/3.0D0)
CsB 2*beta1_here = beta0_PRD
CCPY April 2002: correct the ERROR in B2, inside WSETUP:
          B2c=0.5D0*B2c
CJI Modify to match my notes
         B2= B2c+beta1*(2*A1*dlog(b0/C1)**2-2*A1*dlog(C2)**2+
     >   2*B1c*dlog(C2)) - A2c*2*log(b0*C2/C1)
!     >   2.d0*CF*(
!     >    -((67.d0/18.d0-Pi2/6.d0)*NC - 5.d0/9.d0*NF)*DLog(C2*B0/C1) +
!     >     2.d0*beta1*( DLog(B0/C1)**2-DLog(C2)**2-3.d0/2.d0*Dlog(C2) )
!     >           )

CJI September 2016, add in the choice of scheme
          if(SCHEME.eq.'CFG') then
              B2 = B2 - 0.5*CF*beta1*(PI2-8d0)
          endif
        ElseIf(TYPE_V.EQ.'H+'.OR.TYPE_V.EQ.'HB') then
CCPY Feb 2005: added for canonical choice only
          B2=0.5D0*CF**2*(PI2-3.0D0/4.0D0-12.0D0*ZETA3)+
     >   0.5D0*CF*CA*(11.0D0*PI2/9.0D0-193.0D0/12.0D0+6.0D0*ZETA3)+
     >   0.5D0*CF*TR*(-4.0D0*PI2/9.0D0+17.0D0/3.0D0) 
          B2=0.5D0*B2
C          print*,' B2 (for DY) =',B2
          B2=B2+BETA1*CF*3.0
C          print*,' B2 (for HB) =',B2      
        ElseIf(TYPE_V.EQ.'AA') then
CCPY Nov 2005: added for canonical choice only
C process independent part
          B2=0.5D0*CF**2*(PI2-3.0D0/4.0D0-12.0D0*ZETA3)+
     >   0.5D0*CF*CA*(-11.0D0*PI2/9.0D0-17.0D0/12.0D0+6.0D0*ZETA3)+
     >   0.5D0*CF*TR*(4.0D0*PI2/9.0D0+1.0D0/3.0D0) 
          B2=0.5D0*B2 + BETA1*CF*PI2/6.0D0
c process dependent part          
CCPY PUT IN AN "EFFECTIVE" C(1) TERM FOR 'AA' PRODUCTION
          H1_C2=0.6D0
          B2=B2 + 2.0d0*BETA1*H1_C2
 
         Else
          B2=0.0
        End If
      ELSE
        B2=0.0
      ENDIF
CsB <<<
      IF(DEBUG) THEN
       PRINT*,'PI,PI2,CF,CA,NC,TR,NF,BETA1,EULER,ZETA3'
       PRINT*,PI,PI2,CF,CA,NC,TR,NF,BETA1,EULER,ZETA3
       PRINT*,'A1,A2,B1,B2,C1,C2,C3,C4,B0'
       PRINT*,A1,A2,B1,B2,C1,C2,C3,C4,B0
       PRINT*,'P00,P02,P04,Q01,Q03'
       PRINT*,P00,P02,P04,Q01,Q03
      ENDIF

      B3c = (7358.86-721.516*nf_eff
     >   +20.5951*nf_eff**2)/4**3
      if(SCHEME.eq.'CFG') then
          H1 = 0.5*CF*(PI2-8)
          H2 = CF*CA*(59*ZETA3/18-1535/192d0+215*PI2/216d0-PI2**2/240d0)
     >       + CF**2/4d0*(-15*ZETA3+511/16d0-67*PI2/12d0+17*PI2**2/45d0)
     >       + CF*NF/864d0*(192*ZETA3+1143-152*PI2)
          B3c = B3c - (beta2*H1+2*beta1*(H2-0.5*H1**2))
      endif

CJI B3 scale dependence from my notes
      B3=B3c-2*A3c*log(b0*C2/C1)
     >+2*beta2*(A1c*log(b0/C1)**2+log(C2)*(B1c-A1c*log(C2)))
     >-4.0/3.0*beta1**2*(2*A1c*log(b0/C1)**3
     >        +log(C2)**2*(2*A1c*log(C2)-3*B1c))
     >+4*beta1*(A2c*log(b0/C1)**2+log(C2)*(B2c-A2c*log(C2)))

      RETURN
      END

C ---------------------------------------------------------------------------
      SUBROUTINE HSETUP(NF_EFF)
C ---------------------------------------------------------------------------
C FOR HIGGS PRODUCTION
C DEFINE CONSTANT VARIABLES.
C UNIT:   GeV AND pb
C THIS ONLY WORKS FOR  C1=2*EXP(-EULER)   AND   C2=1
      IMPLICIT NONE
      INCLUDE 'common.for'
      integer i
      REAL*8 HZEROJ0
      DIMENSION HZEROJ0(0:29)

      INTEGER NF_EFF

      DATA HZEROJ0 /0.0, 2.4048256, 5.5200781, 8.6537279, 11.7915344,
     >14.9309177, 18.0710640, 21.2116366, 24.3524715, 27.4934791,
     >30.6346065, 33.7758202, 36.9170984, 40.0584258, 43.1997917,
     >46.3411884, 49.4826099, 52.6240518, 55.7655108, 58.9069839,
     >62.0484692, 65.1899648, 68.3314693,
     >69.9022245, 73.6145006,
     >77.7560256, 80.8975559, 84.0390908, 87.1806298, 90.3221726/

      REAL*8 H1_C2

      real*8 BETA2
      
      DO 11 I=0,29
       ZEROJ0(I)=HZEROJ0(I)
11    CONTINUE
C THE ACCURACY IN DOING B-SPACE INTEGRATIONS IN    SUDNPT
CCPY      ACC_SUDNPT=1.D-3
      ACC_SUDNPT=1.D-4
C      ACC_SUDNPT=1.D-5
      ACC_AERR=1.D-16
CCPY      ACC_RERR=1.D-3
      ACC_RERR=1.D-4
C      ACC_RERR=1.D-5
      ACC_RES1=1.D-3
C FOR NONPERT
      YNPT_CUT=100.0
CCPY March 2019: Add more cycles of ZEROJ0
C FOR RESUM
      NRES1=22
C      NRES1=25
C FOR PDF
      PREV_UNDER=1.D0
C FOR ZJSUDNPT
      GOAROUND=.FALSE.
C      GOAROUND=.TRUE.
C FOR HMRSE, HMRSB PDF
      FIT3=.FALSE.
C
      DEBUG=.FALSE.
C      DEBUG=.TRUE.

C TO SPEED UP
C      ACC_RES1=0.05
C      YNPT_CUT=25.0*SCALING
C      NRES1=INT(7*SCALING)

      P00=1.0D0
      P02=-1.0D0*9.0D0/2.0D0/(8.0D0**2)
      P04=1.0D0*9.0D0*25.0D0*49.0D0/24.0D0/(8.0D0**4)
      Q01=-1.0D0/8.0D0
      Q03=1.0D0*9.0D0*25.0D0/6.0D0/(8.0D0**3)
      EULER=0.577215664901D0
      ZETA3=1.2020569032D0

      NC=3
      CA=NC
      CF=(NC**2-1.0D0)/2.0D0/NC
      NF=NF_EFF
      TR=NF/2.0D0

CPN                                            Calculate the constants C1..C4
      Call SetC1_C4

CsB___A1 and B1
      A1=CA
      BETA1=(11.0D0*CA-4.0D0*TR)/12.0D0
      BETA2=(17*NC**2-5*NC*NF-3*CF*NF)/6.d0
      B1=-2.0*BETA1

CsB___A2 is included for G G -> A A and G G -> H0
      If (N_SUD_A.Eq.1) then
        A2 = 0.d0
CCHEN!  
        A3 = 0.d0
      Else If (N_SUD_A.GE.2) then
        A2=0.5D0*CF*((67.0D0/18.0D0-PI2/6.0D0)*CA-10.0D0/9.0D0*TR)
!ZL only valid for certain combination of C1...C4
        A2=A2*CA/CF
     &      -A1*2*BETA1*log(B0/C1)
CCHEN!  
        A3 = 0.d0
        
        IF(N_SUD_A.EQ.3) THEN
CCHEN! MARCH 2007 ADD IN A3 FOR GG->H, NOT DONE FOR B3 YET!
          A3= CF**2*NF*0.5D0*(ZETA3-55.0/48.0)-CF*NF**2/108.0+
     >      CA**2*CF*(11.0*ZETA3/24.0+11.0*pi**4/720.0-67.0*pi**2/216.0
     >                +245.0/96.0)+
     >      CA*CF*NF*(-7.0*ZETA3/12.0+5.0*pi**2/108.0-209.0/432.0)
!ZL only valid for certain combination of C1...C4
          A3=A3*CA/CF
     &        -2*A2*2*BETA1*log(B0/C1)
     &        -(2*BETA1)**2*A1*log(B0/C1)**2
     &        -A1/2.d0*BETA2*log(B0/C1)
        ENDIF
      End if       

CJI Sept 2013: Seperate the A and B piece to match the desired header format
      IF(N_SUD_B.Eq.1) then
        B2=0.0
      ELSE IF (N_SUD_B.Eq.2) then 
CCPY FEB 2004, INCLUDING B2
CCPY ONLY SETUP FOR CANONICAL C's (IFLAG_C3=1)
        IF(TYPE_V.EQ.'H0') THEN
          B2 = CA*CA*(23.0/24.0+11.0*pi**2/18.0-3.0/2.0*zeta3) +
     >         CF*Tr - 2*CA*Tr*(1.0/12.0+1.0*pi**2/9.0) -
     >         11.0*CF*CA/8.0 
!ZL only valid for certain combination of C1...C4
          B2=B2+2*BETA1*B1*log(C2)

CCPY SEPT 2009        ELSEIF(TYPE_V.EQ.'AG') THEN
        ELSEIF(TYPE_V.EQ.'AG'.OR.TYPE_V.EQ.'ZG') THEN
CCPY Nov 2005: added for canonical choice only
C process independent part
          B2 = CA*CA*(8.0/3.0+3.0*zeta3) -
     >         CF*Tr - 4.0*CA*Tr/3.0
          B2 = -0.5*B2 + BETA1*CA*PI2/6.0d0
c process dependent part          
CCPY PUT IN AN "EFFECTIVE" C(1) TERM FOR 'AG' PRODUCTION
          H1_C2=6.64D0
          B2=B2 + 2.0d0*BETA1*H1_C2
          
        ELSE
          B2=0.D0
        ENDIF         
      End If


      IF(DEBUG) THEN
       PRINT*,'PI,PI2,CF,CA,NC,TR,NF,BETA1,EULER,ZETA3'
       PRINT*,PI,PI2,CF,CA,NC,TR,NF,BETA1,EULER,ZETA3
       PRINT*,'A1,A2,B1,B2,C1,C2,C3,C4,B0'
       PRINT*,A1,A2,B1,B2,C1,C2,C3,C4,B0
       PRINT*,'P00,P02,P04,Q01,Q03'
       PRINT*,P00,P02,P04,Q01,Q03
      ENDIF

      RETURN
      END ! HSETUP

C ---------------------------------------------------------------------------
      SUBROUTINE VJSETUP(NF_EFF)
C ---------------------------------------------------------------------------
C FOR HIGGS + 1 Jet PRODUCTION
C DEFINE CONSTANT VARIABLES.
C UNIT:   GeV AND pb
C THIS ONLY WORKS FOR  C1=2*EXP(-EULER)   AND   C2=1
      IMPLICIT NONE
      INCLUDE 'common.for'
      integer i
      REAL*8 HZEROJ0
      DIMENSION HZEROJ0(0:29)

      Integer HasJet
      REAL*8 D1s, R,t,ptj,u
      COMMON /Jet/ u,D1s, R,t,ptj, HasJet

      INTEGER NF_EFF

      DATA HZEROJ0 /0.0, 2.4048256, 5.5200781, 8.6537279, 11.7915344,
     >14.9309177, 18.0710640, 21.2116366, 24.3524715, 27.4934791,
     >30.6346065, 33.7758202, 36.9170984, 40.0584258, 43.1997917,
     >46.3411884, 49.4826099, 52.6240518, 55.7655108, 58.9069839,
     >62.0484692, 65.1899648, 68.3314693,
     >69.9022245, 73.6145006,
     >77.7560256, 80.8975559, 84.0390908, 87.1806298, 90.3221726/

      REAL*8 H1_C2

      real*8 BETA2

      integer iflavor
      common/flavor/iflavor

      DO 11 I=0,29
       ZEROJ0(I)=HZEROJ0(I)
11    CONTINUE
C THE ACCURACY IN DOING B-SPACE INTEGRATIONS IN    SUDNPT
CCPY      ACC_SUDNPT=1.D-3
      ACC_SUDNPT=1.D-4
C      ACC_SUDNPT=1.D-5
      ACC_AERR=1.D-16
CCPY      ACC_RERR=1.D-3
      ACC_RERR=1.D-4
C      ACC_RERR=1.D-5
      ACC_RES1=1.D-3
C FOR NONPERT
      YNPT_CUT=100.0
CCPY March 2019: Add more cycles of ZEROJ0
C FOR RESUM
      NRES1=22
C      NRES1=25
C FOR PDF
      PREV_UNDER=1.D0
C FOR ZJSUDNPT
      GOAROUND=.FALSE.
C      GOAROUND=.TRUE.
C FOR HMRSE, HMRSB PDF
      FIT3=.FALSE.
C
      DEBUG=.FALSE.
C      DEBUG=.TRUE.

C TO SPEED UP
C      ACC_RES1=0.05
C      YNPT_CUT=25.0*SCALING
C      NRES1=INT(7*SCALING)

      P00=1.0D0
      P02=-1.0D0*9.0D0/2.0D0/(8.0D0**2)
      P04=1.0D0*9.0D0*25.0D0*49.0D0/24.0D0/(8.0D0**4)
      Q01=-1.0D0/8.0D0
      Q03=1.0D0*9.0D0*25.0D0/6.0D0/(8.0D0**3)
      EULER=0.577215664901D0
      ZETA3=1.2020569032D0

      NC=3
      CA=NC
      CF=(NC**2-1.0D0)/2.0D0/NC
      NF=NF_EFF
      TR=NF/2.0D0

CPN                                            Calculate the constants C1..C4
      Call SetC1_C4
      BETA1=(11.0D0*CA-4.0D0*TR)/12.0D0
      BETA2=(17*NC**2-5*NC*NF-3*CF*NF)/6.d0

      If(Type_V.Eq.'HJ') THEN
      if(iflavor.eq.0) then
CsB___A1 and B1
      A1=CA
      B1=-2.0*BETA1
      B2=0

CsB___A2 is included for G G -> A A and G G -> H0
      If (N_SUD_A.Eq.1) then
        A2 = 0.d0
CCHEN!  
        A3 = 0.d0
      Else If (N_SUD_A.GE.2) then
        A2=0.5D0*CF*((67.0D0/18.0D0-PI2/6.0D0)*CA-10.0D0/9.0D0*TR)
!ZL only valid for certain combination of C1...C4
        A2=A2*CA/CF
     &      -A1*2*BETA1*log(B0/C1)
CCHEN!  
        A3 = 0.d0
        
        IF(N_SUD_A.EQ.3) THEN
CCHEN! MARCH 2007 ADD IN A3 FOR GG->H, NOT DONE FOR B3 YET!
          A3= CF**2*NF*0.5D0*(ZETA3-55.0/48.0)-CF*NF**2/108.0+
     >      CA**2*CF*(11.0*ZETA3/24.0+11.0*pi**4/720.0-67.0*pi**2/216.0
     >                +245.0/96.0)+
     >      CA*CF*NF*(-7.0*ZETA3/12.0+5.0*pi**2/108.0-209.0/432.0)
!ZL only valid for certain combination of C1...C4
          A3=A3*CA/CF
     &        -2*A2*2*BETA1*log(B0/C1)
     &        -(2*BETA1)**2*A1*log(B0/C1)**2
     &        -A1/2.d0*BETA2*log(B0/C1)
        ENDIF
      End if       

      D1s = CA
      else if(iflavor.eq.1) then
      A1 = CF/2d0+CA/2d0
      B1 = -BETA1-3d0/4d0*CF-0.5*CA*Log(u/t)+0.5*CF*Log(u/t)
      D1s = CF
      endif 
      A2 = 0.5D0*CF*((67.0D0/18.0D0-PI2/6.0D0)*CA-10.0D0/9.0D0*TR -
     - 4.d0*beta1*DLog(B0/C1))
      ELSE
          A2 =0
      ENDIF
CZL---------------------------------------------------
CZL copied from GG case in the following code
      IF(N_SUD_A.GE.3) THEN
        A3= CF**2*NF*0.5D0*(ZETA3-55.0/48.0)-CF*NF**2/108.0+
     >      CA**2*CF*(11.0*ZETA3/24.0+11.0*pi**4/720.0-67.0*pi**2/216.0
     >                +245.0/96.0)+
     >      CA*CF*NF*(-7.0*ZETA3/12.0+5.0*pi**2/108.0-209.0/432.0)
!ZL only valid for certain combination of C1...C4
        A3=A3
     &     -2*A2*2*BETA1*log(B0/C1)
     &     -(2*BETA1)**2*A1*log(B0/C1)**2
     &     -A1/2.d0*BETA2*log(B0/C1)
        ENDIF
CZL---------------------------------------------------
CsB B1_Here = B1_PRD
      B1 = CF*(-3.0D0/2.0D0 - 2.d0*DLog(C2*B0/C1))

      D1s = CA/2.0

      IF(DEBUG) THEN
       PRINT*,'PI,PI2,CF,CA,NC,TR,NF,BETA1,EULER,ZETA3'
       PRINT*,PI,PI2,CF,CA,NC,TR,NF,BETA1,EULER,ZETA3
       PRINT*,'A1,A2,B1,B2,C1,C2,C3,C4,B0'
       PRINT*,A1,A2,B1,B2,C1,C2,C3,C4,B0
       PRINT*,'P00,P02,P04,Q01,Q03'
       PRINT*,P00,P02,P04,Q01,Q03
      ENDIF

      RETURN
      END ! VJSETUP

C--------------------------------------------------------------------------
      subroutine SetC1_C4
C--------------------------------------------------------------------------
CPN   Jan. 1999
      implicit NONE
      include 'common.for'

CJI September 2016, add in the choice of scheme
      CHARACTER*3 SCHEME
      COMMON / RESUMSCHEME / SCHEME
      real*8 A1c,A2c,A3c,B1c,B2c,B3c,C1_IN,C2_IN,C3_IN
      common /canonical/ a1c,a2c,a3c,b1c,b2c,b3c,C1_IN,C2_IN,C3_IN

      double precision tag

      B0=2.0D0*DEXP(-EULER) ! the numeric value of this is: B0 = 1.12291897

CCPY SCHEME = 'CSS' or 'CFG' 
CFix the ratio of C2/C1, there is no C4 dependence in CSS framework for DY production
       IF(SCHEME.EQ.'CSS') C4=1.D0

CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
      IF(IFLAG_C3.EQ.1) THEN
C THIS IS THE CANONICAL CHOICE FOR W-BOSON
CZL
       C2=1.0D0
       C1=C2*B0
CCPY JUNE 2013       C3=C2*B0
       C3=B0
       C4=C2
      ELSEIF(IFLAG_C3.EQ.2) THEN
C The renormalization scale is set to Q/2 both for the resummed and the fixed
C order pieces.
       C2 = 0.5D0
       C1=C2*B0
       C3 = 0.5*B0
       C4 = C2
CsB The reason we don't want to set C4 different from 1 is that we don't want to
C   spoil th ecancellation between the Asymptotic and the Perturbative pieces in
C   the small QT region.
      ELSEIF(IFLAG_C3.EQ.3) THEN
C ARBITRARY CHOICE FOR W-BOSON
       C2=0.25D0
       C1=C2*B0
       C3=B0
       C4=C2
      ELSEIF(IFLAG_C3.EQ.4) THEN
C ARBITRARY CHOICE FOR W-BOSON
CCPY JUNE 2013       C2=1.5D0
C       C1=C2*B0
C       C3=B0
C       C4=C2
C
       C2=2.0D0
       C1=C2*B0
       C3=2*B0
       C4=C2
      ELSEIF(IFLAG_C3.EQ.5) THEN
C Checking the case: ratio of C2/C1=2.
       C2=1.d0
       C1=0.5d0
       C3=B0
       C4=C2
      ELSEIF(IFLAG_C3.EQ.6) THEN
       C2=1.d0
       C1=C2*B0
       C3=2*B0
       C4=1d0
      ELSEIF(IFLAG_C3. EQ. 99) THEN
CCPY This assignment of IFLAG_C3=99 should be recognized by legacy.in
        C2=C2_IN
        C1=C1_IN*B0
        C3=C3_IN*B0
        C4=1.D0
      ELSE
       WRITE(NWRT,*) ' NO SUCH IFLAG_C3'
       CALL QUIT
      ENDIF


      Return
      End !SetC1_4

C ---------------------------------------------------------------------------
      SUBROUTINE SET_ZFF2
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      integer i,j

      Real*8 Sw_Eff,Cw_Eff,Tw_Eff, T3u,T3d,Qu,Qd

      QU= 2.d0/3.d0
      QD=-1.d0/3.d0

      IF(TYPE_V.EQ.'A0') THEN
C COUPLINGS OF A0-F-F.  (g*Sin_w)*(\gamma_\mu)*(electric charge)
        ZFF2(1)=(2.0/3.0)**2
        ZFF2(2)=(-1.0/3.0)**2
        ZFF2(3)=(-1.0/3.0)**2
        ZFF2(4)=(2.0/3.0)**2
        ZFF2(5)=(-1.0/3.0)**2
        ZFF2(6)=(2.0/3.0)**2
CCPY Error: April 14, 2004 (correct ALEPASY for A0, etc)
        IF(LEPASY.EQ.1) THEN
          DO I=1,6
           ZFF2(I)=0.D0
          ENDDO
        ENDIF
      ELSEIF(TYPE_V.EQ.'AA') THEN
        ZFF2(1)=(2.0/3.0)**4
        ZFF2(2)=(-1.0/3.0)**4
        ZFF2(3)=(-1.0/3.0)**4
        ZFF2(4)=(2.0/3.0)**4
        ZFF2(5)=(-1.0/3.0)**4
        ZFF2(6)=(2.0/3.0)**4
CCPY Error: April 14, 2004 (correct ALEPASY for A0, etc)
        IF(LEPASY.EQ.1) THEN
          DO I=1,6
           ZFF2(I)=0.D0
          ENDDO
        ENDIF

      ELSEIF(TYPE_V.EQ.'ZZ') THEN
CsB_____See: Ohnemus, Owens PRD43 (91) 3627 (3)
C       We write Qq^4 = g+^4/2 + g-^4/2 with g+ = g- = Qq for qQ->AA,
C       and replace g+ and g- for ZZ.
        Sw_Eff = Sqrt(SW2_EFF)
        Cw_Eff = Sqrt(1.d0 - SW2_EFF)
        Tw_Eff = Sw_Eff/Cw_Eff
        T3u =  1.d0/2.d0
        T3d = -1.d0/2.d0
        ZFF2(1) = (T3u/Sw_Eff/Cw_Eff-Qu*Tw_Eff)**4/2.d0 +
     +            (-Qu*Tw_Eff)**4/2.d0
        ZFF2(2) = (T3d/Sw_Eff/Cw_Eff-Qd*Tw_Eff)**4/2.d0 +
     +            (-Qd*Tw_Eff)**4/2.d0
        ZFF2(3)=ZFF2(2)
        ZFF2(4)=ZFF2(1)
        ZFF2(5)=ZFF2(2)
        ZFF2(6)=ZFF2(1)
      ELSEIF(TYPE_V.EQ.'WW_UUB'.or.TYPE_V.EQ.'WW_DDB') THEN
        DO I=1,6
          ZFF2(I)=1d0/2d0
        ENDDO
CCPY The case of LEPASY=1 needs to be checked.
      ELSEIF(TYPE_V.EQ.'Z0') THEN
C COUPLINGS OF Z0-F-F.  -(g/2/Cos_w)*(\gamma_\mu)*(V^I-A^I*\gamma_5)
cpn Couplings depend on lepasy 
ccpy 
       DO I=1,6
        IF(I.EQ.1 .OR. I.EQ.4 .OR. I.EQ.6) THEN
         ZFF_A(I)=1.0D0/2.0D0
         ZFF_V(I)=1.0D0/2.0D0-4.0D0/3.0D0*SW2_EFF
        ELSE
         ZFF_A(I)=-1.0D0/2.0D0
         ZFF_V(I)=-1.0D0/2.0D0+2.0D0/3.0D0*SW2_EFF
        ENDIF
        IF(LEPASY.EQ.1) THEN
          ZFF2(I)=2.0D0*ZFF_A(I)*ZFF_V(I)
        ELSE
          ZFF2(I)=ZFF_A(I)**2+ZFF_V(I)**2
        ENDIF
       ENDDO
ccpy      
CJI Sept 2013: Added in ZU and ZD to match c++ code 
      ELSEIF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     >    .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') then
        DO I=1,6
          ZFF2(I)=1.0D0
        ENDDO
C
      ELSEIF(TYPE_V.EQ.'GL') THEN
C COUPLINGS OF GL-F-F.  (g_strong)*(\gamma_\mu)*(color matrix)
        DO I=1,6
          ZFF2(I)=1.0D0
        ENDDO
CCPY Error: April 14, 2004 (correct ALEPASY for A0, etc)
        IF(LEPASY.EQ.1) THEN
          DO I=1,6
           ZFF2(I)=0.D0
          ENDDO
        ENDIF
      ENDIF

      RETURN
      END
      


C ---------------------------------------------------------------------------
      SUBROUTINE STANDARD
C ---------------------------------------------------------------------------
C SET UP PARAMETERS FOR THE ELECTROWEAK SECTOR
C USING ON-SHELL SCHEME
C FLAVOR LABEL:
C UP=1, DN=2, ST=3, CH=4 BT=5 TP=6, GL=0
C UB=-1, DB=-2, SB=-3, CB=-4 BB=-5 TB=-6

      IMPLICIT NONE
      INCLUDE 'common.for'
      integer i,j

      REAL*8 ALPI_MT,EW_DELKA
      REAL*8 KM_L,KM_A,KM_RHO,KM_ETA
      REAL*8 CW2,EW_A,EW_DELA,EW_DELR,
     >EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW,
     >QU,QD,XLU,XLD,XRU,XRD

      Real*8 GWeak,GWEAK2, Temp,
     >       Sw_Eff,Cw_Eff,Tw_Eff, T3u,T3d

      Real*8 GAMH,GAMH2,XMQX, A_W, A02

      Real*8 COUP_XMT,COUP_XMC,COUP_MB,CR,CR_CB,UR_TC,XLTC,TWO

      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC

CCAO! 
      REAL*8 GAMW_NLO,GAMZ_NLO,ALFA_EM_MZ
      REAL*8 AEM_MZ,AEM
      REAL*8 ZCOUPL_OLD,ZCOUPL_MZ
      REAL*8 ALPHAS_MW,ALPHAS_MZ
      REAL*8 QE,VUB,XLE,XRE
      REAL*8 GAMW_MZ,GAMZ_MZ
      REAL*8 SW2_EFF_OLD,SW2_EFF_MZ
      REAL*8 BR_Z_EE,BR_W_EN

      Logical First
      Data First /.True./

CCAO!      Data PI/3.14159 26535 89793 24D0/
      PI=4.0D0*ATAN(1.0D0)      
      PI2=PI**2
      R2=DSQRT(2.0D0)

C UNIT  CONVERT ENERGY FROM  GeV  TO  TeV, THEN UNIT=1.0D-3
c WARNING!!! You must keep UNIT=1.0D0 in this version of the code!!!!!!
      UNIT=1.0D0

C_____CONVERT FROM  GeV  TO pb
      HBARC2=UNIT**2*0.38937966D0*1.0D9

C_____INPUT DATA
      GMU=1.1663787D-5/UNIT**2
      EW_ALFA=1.0D0/137.035999074D0

C_____W, Z photon, top and H masses are read from input file
      MW2=MW**2
      MZ2=MZ**2

C_____EW_DELA=0.0595 +/- 0.0009
CCAO!      EW_DELA=0.0595+40.0D0/9.0D0*EW_ALFA/PI*DLOG(xMZ/91.187D0)
C          EW_DELA=0.02767+0.03150-0.00007=0.0591
      EW_DELA=0.0591+40.0D0/9.0D0*EW_ALFA/PI*DLOG(MZ/91.187D0/UNIT)

C_____EW_ALFA AT  MZ
      EW_AMZ=EW_ALFA/(1.0D0-EW_DELA)
      EW_A=PI*EW_ALFA/R2/GMU*UNIT**2
      ALPI_MT=ALPI(MT)
      
      IF(FIRST) THEN
        PRINT*,' MT, PI*ALPI(XMTOP)=',MT, ALPI_MT*PI
      ENDIF

C_____Rho and Delta rho
      EW_DELRHO=3.0D0*GMU*MT**2/8.0D0/R2/PI2
      EW_DELRHO=EW_DELRHO*(1.D0-2.0/3.0*(PI2/3.0+1.0)*ALPI_MT)
      EW_RHO=1.0D0/(1.0D0-EW_DELRHO)

C_____ON-SHELL SCHEME
      SW2=1.0D0-MW2/MZ2
      SW=DSQRT(SW2)
      CW2=1.0D0-SW2
      CW=DSQRT(CW2)

CCPY      EW_DELKA=CW2/SW2*EW_DELRHO
C      EW_DELR=EW_DELA-EW_DELKA
C      print*,' old ew_delr,EW_DELA,EW_DELKA =',ew_delr,EW_DELA,EW_DELKA
      
CPN___Calculate delta_r using mW
      A02 = pi*ew_amz/r2/Gmu
      ew_delr = 1 - A02/mW2/sw2 

      EW_DELKA =  EW_DELA - EW_DELR
      IF(FIRST) THEN
      print*,' new ew_delr,EW_DELKA  =',ew_delr,EW_DELKA 
      ENDIF
      
CCAO!      SW2_EFF=SW2*(1.0d0 + EW_DELKA)
      SW2_EFF_OLD=SW2*(1.0d0 + EW_DELKA)

CCPY THIS IS FOR MZ=91.187 and  MW=80.385
CJI Change SW2_EFF_MZ to match that of resbos code
CJI      SW2_EFF_MZ=0.23143
      SW2_EFF_MZ=0.231676
      SW2_EFF=SW2_EFF_MZ

      Sw_Eff = Sqrt(SW2_EFF)
      Cw_Eff = Sqrt(1.d0 - SW2_EFF)
      Tw_Eff = Sw_Eff/Cw_Eff

      IF(FIRST) THEN
      WRITE(*,*) ' MW, MZ =',MW,MZ
      WRITE(*,*) ' SW2, SW2_EFF_OLD, SW2_EFF_MZ, SW2_EFF =',
     > SW2,SW2_EFF_OLD,SW2_EFF_MZ,SW2_EFF
      ENDIF
      
C_____alpha = e^2/(4 Pi)
C      AEM=1.D0/128.0
      AEM = Sqrt(2.d0)*MW2*GMU*SW2/Pi

      AEM_MZ=1.0/128.937
      ALFA_EM_MZ=EW_AMZ
      IF(FIRST) THEN
      WRITE(*,*)'1/AEM , 1/EW_AMZ, 1/AEM_MZ =',1/AEM,1/EW_AMZ,1/AEM_MZ
      ENDIF
C_____Weak couplings: g_w and Fermi constant
      GWEAK2=AEM*4.d0*PI/SW2
      GWEAK = Sqrt(GWEAK2)
      GFERMI=PI*AEM/SQRT(2.d0)/SW2/MW2

CC_____L/R Quark couplings
      QU= 2.d0/3.d0
      QD=-1.d0/3.d0
      XLU=1.d0-2.d0*SW2_eff*QU
      XLD=-1.d0-2.d0*SW2_eff*QD
      XRU=-2.d0*QU*SW2_eff
      XRD=-2.d0*QD*SW2_eff

      QE=-1.0
      XLE=-1.D0-2.*SW2_EFF*QE
      XRE=-2.D0*QE*SW2_EFF

C WCOUPL IS (g**2)/8
      WCOUPL=R2*MW2*GMU/2.0D0
C WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS

C ZCOUPL IS (g/Cos_w)**2/4
      ZCOUPL=R2*MZ2*GMU*(1.0+EW_DELRHO)

C ACOUPL IS (e**2)/2
      ACOUPL=R2*MW2*GMU*2.0D0*SW2
C WCOUPL IS (g**2)/8
      WCOUPL=R2*MW2*GMU/2.0D0
C WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS
C ZCOUPL IS (g/Cos_w)**2/4
      ZCOUPL_OLD=R2*MZ2*GMU*(1.0+EW_DELRHO)
      ZCOUPL_MZ=PI*ALFA_EM_MZ/SW2/CW2
      ZCOUPL=ZCOUPL_MZ
      
      IF(FIRST) THEN
      WRITE(*,*) 'ZCOUPL_OLD, ZCOUPL_MZ, ZCOUPL =',
     > ZCOUPL_OLD, ZCOUPL_MZ, ZCOUPL
      ENDIF
      
C HCOUPL IS FOR HIGGS
      HCOUPL=R2*GMU*PI/576.0D0

C HPCOUPL is 1/2 CR^2/4 for q H+ Q
      If (Type_V.Eq.'H+') then
      IF (I_MODEL.EQ.1) THEN  ! TOPCOLOR MODEL
        TANB = 3.0
        VEV = 246.0
CsB     Moved top mass factor into Subroutine SResum, since Subroutine
C       Standard does not know about the Q_V dependence (it is called only
C       once in the beginning).
        COUP_XMT=1.d0
        CR = TANB*(SQRT(2.0)*COUP_XMT/VEV)
        UR_TC=0.2
        CR_CB = CR*UR_TC
!        PREFAC = PI/12/SS * CR_CB**2
        HpCOUPL = 1./8. * CR_CB**2
        IF(FIRST) THEN
          Write(*,*) ' '
          WRITE(*,*) ' Q_V,VEV,TANB,COUP_XMT,CR,UR_TC,CR_CB '
          WRITE(*,*) Q_V,VEV,TANB,COUP_XMT,CR,UR_TC,CR_CB
          FIRST=.FALSE.
        ENDIF
      ELSE IF (I_MODEL.EQ.2) THEN ! FOR 2HDMIII:
        XLTC = 1.5D0*SQRT(2.0)
        VEV=246.D0
C COUPLINGS OF T-BBAR-H^+
CsB     Moved top mass factor into Subroutine SResum, since Subroutine
C       Standard does not know about the Q_V dependence (it is called only
C       once in the beginning).
        COUP_XMT=1.d0
        COUP_XMC=1.d0
C COUPLINGS OF C-BBAR-H^+
        CR=SQRT(COUP_XMT*COUP_XMC)/VEV
        CR_CB=XLTC*CR
!        PREFAC = PI/12/SS * CR_CB**2
        HpCOUPL = 1./8. * CR_CB**2
        IF(FIRST) THEN
          WRITE(11,*) ' Q_V,VEV,COUP_XMT,COUP_XMC,CR,XLTC,CR_CB '
          WRITE(11,*) Q_V,VEV,COUP_XMT,COUP_XMC,CR,XLTC,CR_CB
          FIRST=.FALSE.
        ENDIF
      ELSE
        WRITE(*,*) ' NO SUCH MODEL '
        Print*, 'Type_V, I_MODEL = ', Type_V,I_MODEL
        CALL QUIT
      ENDIF
      End If

C HBCOUPL is 1/2 CR^2/4 for b A0 B
      If (Type_V.Eq.'HB') then
        TANB = 50.0
        VEV = 246.220569
CsB     Moved top mass factor into Subroutine SResum, since Subroutine
C       Standard does not know about the Q_V dependence (it is called only
C       once in the beginning).
        COUP_MB=1.d0
        CR = TANB*(COUP_MB/VEV)
!        PREFAC = PI/12/SS * CR**2
        HbCOUPL = 1./8. * CR**2
CsB     ADD A FACTOR OF 2 TO ACCOUNT FOR LEFT- AND RIGHT-HANDED CHIRALITIES
        TWO=2.0
        HbCOUPL = TWO * HbCOUPL
        IF(FIRST) THEN
          Write(*,*) ' '
          WRITE(11,*) ' Q_V,VEV,TANB,COUP_MB,CR,TWO,HbCOUPL '
          WRITE(11,*)  Q_V,VEV,TANB,COUP_MB,CR,TWO,HbCOUPL
          FIRST=.FALSE.
        ENDIF
      End If

C_____Vector boson widths
      GAMW=GFERMI*MW**3/6./PI/SQRT(2.)*(9.)
      GAMW2=GAMW**2
C NLO QCD CORRECTION
CCPY      ALPHAS_MW=ALPI(MW)*PI      

      ALPHAS_MW=ALPI(MW)*PI
      IF(FIRST) THEN
      PRINT*,' MW, ALPHAS_MW =',MW, ALPHAS_MW 
      ENDIF
      
      GAMW_NLO=GAMW*(1.0+2.0*ALPHAS_MW/PI/3.0)
      GAMW_MZ=2.093*(MW/80.385)**3
      
      IF(FIRST) THEN
      WRITE(*,*) ' ALPHAS_MW, GAMW, GAMW_NLO, GAMW_MZ =',
     >            ALPHAS_MW, GAMW, GAMW_NLO, GAMW_MZ  
      ENDIF
      
      TEMP=MZ**3/DSQRT(2.d0)/12.d0/PI*GFERMI
      GAMZ=TEMP*(3.d0+3.d0*(XLE**2+XRE**2)+
     +     6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))
      GAMZ2=GAMZ**2
      
C NLO QCD CORRECTION
CCPY      ALPHAS_MZ=ALPI(XMZ)*PI

      ALPHAS_MZ=ALPI(MZ)*PI

      IF(FIRST) THEN
      PRINT*,' MZ, ALPHAS_MZ =',MZ, ALPHAS_MZ 
      ENDIF
      
      GAMZ_NLO=TEMP*(3.d0+3.d0*(XLE**2+XRE**2)+
     >     (6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))*
     >     (1.0+ALPHAS_MZ/PI))
      GAMZ_MZ=2.4960*(MZ/91.1875)**3
     
      IF(FIRST) THEN
      WRITE(*,*) ' ALPHAS_MZ, GAMZ, GAMZ_NLO, GAMZ_MZ =',
     >            ALPHAS_MZ, GAMZ, GAMZ_NLO, GAMZ_MZ        
      ENDIF
      
CCPY 
      BR_W_EN=1.0/9.0/(1.0+2.0*ALPHAS_MW/PI/3.0)
      BR_Z_EE=TEMP*(XLE**2+XRE**2)/GAMZ_NLO

C_____Higgs width
C      Higgs mass is read from the input file
      IF(FIRST) THEN
      Print*, ' Higgs mass =', MH
      ENDIF
      
      GAMH=0.D0
      DO 57 I=1,2
        XMQX = XMBOT
        IF (I.EQ.2) XMQX = XMBOT
        IF(MH.GT.(2.*XMQX))THEN
          GAMH=GAMH+gmu*3.*MH*XMQX**2/4./SQRT(2.)/PI*
     1 (1.-4.*XMQX**2/MH**2)**1.5
        END IF
57    CONTINUE
      IF(MH.GT.(2.*MW))THEN
        A_W=4.*MW**2/MH**2
        GAMH=GAMH+gmu*MH**3/32./SQRT(2.)/PI*
     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      END IF
      IF(MH.GT.(2.*MZ))THEN
        A_W=4.*MZ**2/MH**2
        GAMH=GAMH+gmu*MH**3/64./SQRT(2.)/PI*
     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      END IF
      GAMH2=GAMH**2
      IF(FIRST) THEN
      Print*, ' Gamma(Higgs) = ', GAMH
      ENDIF
      
C (KM MATRIX ELEMENT)**2
C FLAVOR LABEL:
C UP=1, DN=2, ST=3, CH=4 BT=5 TP=6, GL=0
C UB=-1, DB=-2, SB=-3, CB=-4 BB=-5 TB=-6
C TOP QUARK CHANNEL IS SWITCHED OFF
C THIS ONLY WORKS IN    QCD   CALCULATIONS
CCAO!      KM_L=0.22D0
CCAO!      KM_A=1.0D0
CCAO!      KM_RHO=-0.4D0
CCAO!      KM_ETA=0.22D0
      KM_L=0.224D0
      KM_A=0.83
      KM_ETA=0.35D0
      VUB=35.7D-4      
C      KM_RHO=0.155D0
      KM_RHO=DSQRT((VUB/KM_A/KM_L**3)**2-KM_ETA**2)
      IF(FIRST) THEN
      WRITE(*,*) ' CKM PARAMETERS:    KM_L, KM_A, KM_ETA, KM_RHO = ',
     >         KM_L, KM_A, KM_ETA, KM_RHO     
      ENDIF
      
CCAO!      WRITE(71,*) ' CKM PARAMETERS:   KM_L, KM_A, KM_ETA, KM_RHO = ',
CCAO!     >         KM_L, KM_A, KM_ETA, KM_RHO     

CCPY  CKM PARAMETERS:    KM_L,   KM_A,       KM_ETA,  KM_RHO =   
CCPY                     0.224  0.829999983  0.35   0.154760903

       DO I=1,6
       DO J=1,6
        IF(I.EQ.J) THEN
         VKM(I,J)=1.0D0
        ELSE
         VKM(I,J)=0.0D0
        ENDIF
       ENDDO
      ENDDO


CCPY APRIL 2004
C      VKM(1,2)=(1.0D0-KM_L**2/2.0D0)**2
C      VKM(1,3)=(KM_L)**2
C      VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
C      VKM(4,2)=VKM(1,3)
C      VKM(4,3)=VKM(1,2)
C      VKM(4,5)=(KM_A*KM_L**2)**2
C      VKM(2,1)=VKM(1,2)
C      VKM(3,1)=VKM(1,3)
C      VKM(5,1)=VKM(1,5)
C      VKM(2,4)=VKM(4,2)
C      VKM(3,4)=VKM(4,3)
C      VKM(5,4)=VKM(4,5)

      VKM(1,2)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0)**2
      VKM(1,3)=(KM_L)**2
      VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
      VKM(4,2)=(-KM_L+0.5*KM_A**2*KM_L**5*
     >          DSQRT((1.0-2.0*KM_RHO)**2+KM_ETA**2))**2
      VKM(4,3)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0*
     >         (1.0+4.0*KM_A**2))**2
      VKM(4,5)=(KM_A*KM_L**2)**2
      VKM(2,1)=VKM(1,2)
      VKM(3,1)=VKM(1,3)
      VKM(5,1)=VKM(1,5)
      VKM(2,4)=VKM(4,2)
      VKM(3,4)=VKM(4,3)
      VKM(5,4)=VKM(4,5)

      LEPASY=0
      CALL SET_ZFF2
      
      IF(FIRST) THEN
       PRINT*,'GMU,SW2,CW2,SW2_EFF,EW_A,EW_ALFA,EW_DELA,EW_DELR,'
       Print*,' EW_DELKA,EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW '
       PRINT*,GMU,SW2,CW2,SW2_EFF,EW_A,EW_ALFA,EW_DELA,EW_DELR
       Print*,EW_DELKA,EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW
       PRINT*,'MW,MZ,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL'
       PRINT*,MW,MZ,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL
       PRINT*,'ZFF2'
       PRINT*,ZFF2
       PRINT*, 'VKM(1,2),VKM(1,3),VKM(1,5),SUM_1'
       PRINT*,VKM(1,2),VKM(1,3),VKM(1,5), 
     >              VKM(1,2)+VKM(1,3)+VKM(1,5)
       PRINT*,'VKM(4,2),VKM(4,3),VKM(4,5),SUM_4'
       PRINT*,VKM(4,2),VKM(4,3),VKM(4,5), 
     >              VKM(4,2)+VKM(4,3)+VKM(4,5)
       PRINT*,'1/ALFA_EM_MZ  (ALPHA_{EM} AT MZ) =',1.0/ALFA_EM_MZ
       PRINT*,' SW2_EFF =',SW2_EFF
       PRINT*, 'GAMW_NLO, GAMZ_NLO =',GAMW_NLO, GAMZ_NLO 

       PRINT*,'GFermi,Gmu'
       PRINT*, GFermi,Gmu
       PRINT*, 'BR_Z_EE =', BR_Z_EE
       PRINT*, 'BR_W_EN =', BR_W_EN
       FIRST=.FALSE.
      ENDIF

      RETURN
      END

CC ---------------------------------------------------------------------------
      SUBROUTINE STANDARD_042004
C ---------------------------------------------------------------------------
C SET UP PARAMETERS FOR THE ELECTROWEAK SECTOR
C USING ON-SHELL SCHEME
C FLAVOR LABEL:
C UP=1, DN=2, ST=3, CH=4 BT=5 TP=6, GL=0
C UB=-1, DB=-2, SB=-3, CB=-4 BB=-5 TB=-6

      IMPLICIT NONE
      INCLUDE 'common.for'
      integer i,j

      REAL*8 ALPI_MT,EW_DELKA
      REAL*8 KM_L,KM_A,KM_RHO,KM_ETA
      REAL*8 CW2,EW_A,EW_DELA,EW_DELR,
     >EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW,
     >QU,QD,XLU,XLD,XRU,XRD

      Real*8 GWeak,GWEAK2, Temp,
     >       Sw_Eff,Cw_Eff,Tw_Eff, T3u,T3d

      Real*8 GAMH,GAMH2,XMQX, A_W, A02

      Real*8 COUP_XMT,COUP_XMC,COUP_MB,CR,CR_CB,UR_TC,XLTC,TWO

      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC

      Logical First
      Data First /.True./

CsB___Parameters
      PI = 3.14159 26535 89793 24D0
      PI2=PI**2
      R2=DSQRT(2.0D0)

C UNIT  CONVERT ENERGY FROM  GeV  TO  TeV, THEN UNIT=1.0D-3
c WARNING!!! You must keep UNIT=1.0D0 in this version of the code!!!!!!
      UNIT=1.0D0

C_____CONVERT FROM  GeV  TO pb
      HBARC2=UNIT**2*0.38937966D0*1.0D9

C_____INPUT DATA
      GMU=1.166389D-5/UNIT**2
      EW_ALFA=1.0D0/137.0359895D0

C_____W, Z photon, top and H masses are read from input file
      MW2=MW**2
      MZ2=MZ**2

C_____EW_DELA=0.0595 +- 0.0009
      EW_DELA=0.0595+40.0D0/9.0D0*EW_ALFA/PI*DLOG(MZ/91.187D0/UNIT)

C_____EW_ALFA AT  MZ
      EW_AMZ=EW_ALFA/(1.0D0-EW_DELA)
      EW_A=PI*EW_ALFA/R2/GMU*UNIT**2
      ALPI_MT=0.12D0/PI

C_____Rho and Delta rho
      EW_DELRHO=3.0D0*GMU*MT**2/8.0D0/R2/PI2
      EW_DELRHO=EW_DELRHO*(1.D0-2.0/3.0*(PI2/3.0+1.0)*ALPI_MT)
      EW_RHO=1.0D0/(1.0D0-EW_DELRHO)

C_____ON-SHELL SCHEME
      SW2=1.0D0-MW2/MZ2
      SW=DSQRT(SW2)
      CW2=1.0D0-SW2
      CW=DSQRT(CW2)
      EW_DELKA=CW2/SW2*EW_DELRHO
      SW2_EFF=SW2*(1.0+EW_DELKA)


CPN___Calculate delta_r using mW
      A02 = pi*ew_alfa/r2/Gmu
      ew_delr = 1 - A02/mW2/sw2 


C_____L/R Quark couplings
      QU= 2.d0/3.d0
      QD=-1.d0/3.d0
      XLU=1.d0-2.d0*SW2*QU
      XLD=-1.d0-2.d0*SW2*QD
      XRU=-2.d0*QU*SW2
      XRD=-2.d0*QD*SW2

C WCOUPL IS (g**2)/8
      WCOUPL=R2*MW2*GMU/2.0D0
C WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS

C ZCOUPL IS (g/Cos_w)**2/4
      ZCOUPL=R2*MZ2*GMU*(1.0+EW_DELRHO)

C HCOUPL IS FOR HIGGS
      HCOUPL=R2*GMU*PI/576.0D0

C HPCOUPL is 1/2 CR^2/4 for q H+ Q
      If (Type_V.Eq.'H+') then
      IF (I_MODEL.EQ.1) THEN  ! TOPCOLOR MODEL
        TANB = 3.0
        VEV = 246.0
CsB     Moved top mass factor into Subroutine SResum, since Subroutine
C       Standard does not know about the Q_V dependence (it is called only
C       once in the beginning).
        COUP_XMT=1.d0
        CR = TANB*(SQRT(2.0)*COUP_XMT/VEV)
        UR_TC=0.2
        CR_CB = CR*UR_TC
!        PREFAC = PI/12/SS * CR_CB**2
        HpCOUPL = 1./8. * CR_CB**2
        IF(FIRST) THEN
          Write(*,*) ' '
          WRITE(*,*) ' Q_V,VEV,TANB,COUP_XMT,CR,UR_TC,CR_CB '
          WRITE(*,*) Q_V,VEV,TANB,COUP_XMT,CR,UR_TC,CR_CB
          FIRST=.FALSE.
        ENDIF
      ELSE IF (I_MODEL.EQ.2) THEN ! FOR 2HDMIII:
        XLTC = 1.5D0*SQRT(2.0)
        VEV=246.D0
C COUPLINGS OF T-BBAR-H^+
CsB     Moved top mass factor into Subroutine SResum, since Subroutine
C       Standard does not know about the Q_V dependence (it is called only
C       once in the beginning).
        COUP_XMT=1.d0
        COUP_XMC=1.d0
C COUPLINGS OF C-BBAR-H^+
        CR=SQRT(COUP_XMT*COUP_XMC)/VEV
        CR_CB=XLTC*CR
!        PREFAC = PI/12/SS * CR_CB**2
        HpCOUPL = 1./8. * CR_CB**2
        IF(FIRST) THEN
          WRITE(11,*) ' Q_V,VEV,COUP_XMT,COUP_XMC,CR,XLTC,CR_CB '
          WRITE(11,*) Q_V,VEV,COUP_XMT,COUP_XMC,CR,XLTC,CR_CB
          FIRST=.FALSE.
        ENDIF
      ELSE
        WRITE(*,*) ' NO SUCH MODEL '
        Print*, 'Type_V, I_MODEL = ', Type_V,I_MODEL
        CALL QUIT
      ENDIF
      End If

C HBCOUPL is 1/2 CR^2/4 for b A0 B
      If (Type_V.Eq.'HB') then
        TANB = 50.0
        VEV = 246.0
CsB     Moved top mass factor into Subroutine SResum, since Subroutine
C       Standard does not know about the Q_V dependence (it is called only
C       once in the beginning).
        COUP_MB=1.d0
        CR = TANB*(COUP_MB/VEV)
!        PREFAC = PI/12/SS * CR**2
        HbCOUPL = 1./8. * CR**2
CsB     ADD A FACTOR OF 2 TO ACCOUNT FOR LEFT- AND RIGHT-HANDED CHIRALITIES
        TWO=2.0
        HbCOUPL = TWO * HbCOUPL
        IF(FIRST) THEN
          Write(*,*) ' '
          WRITE(11,*) ' Q_V,VEV,TANB,COUP_MB,CR,TWO,HbCOUPL '
          WRITE(11,*)  Q_V,VEV,TANB,COUP_MB,CR,TWO,HbCOUPL
          FIRST=.FALSE.
        ENDIF
      End If

C_____Vector boson widths
      GAMW=gmu*MW**3/6.d0/PI/DSQRT(2.d0)*(9.d0)
      GAMW2=GAMW**2
      TEMP=MZ**3/DSQRT(2.d0)/12.d0/PI*gmu
      GAMZ=TEMP*(3.d0+3.d0*(1.d0-4.d0*SW2+8.d0*SW2**2)+
     +     6.d0*(XLU**2+XRU**2)+9.*(XLD**2+XRD**2))
      GAMZ2=GAMZ**2

C_____Higgs width
C      Higgs mass is read from the input file
      Print*, ' Higgs mass =', MH
      GAMH=0.D0
      DO 57 I=1,2
        XMQX = XMBOT
        IF (I.EQ.2) XMQX = XMBOT
        IF(MH.GT.(2.*XMQX))THEN
          GAMH=GAMH+gmu*3.*MH*XMQX**2/4./SQRT(2.)/PI*
     1 (1.-4.*XMQX**2/MH**2)**1.5
        END IF
57    CONTINUE
      IF(MH.GT.(2.*MW))THEN
        A_W=4.*MW**2/MH**2
        GAMH=GAMH+gmu*MH**3/32./SQRT(2.)/PI*
     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      END IF
      IF(MH.GT.(2.*MZ))THEN
        A_W=4.*MZ**2/MH**2
        GAMH=GAMH+gmu*MH**3/64./SQRT(2.)/PI*
     1 SQRT(1.-A_W)*(4.-4.*A_W+3.*A_W**2)
      END IF
      GAMH2=GAMH**2
      Print*, ' Gamma(Higgs) = ', GAMH

C (KM MATRIX ELEMENT)**2
C TOP QUARK CHANNEL IS SWITCHED OFF
C THIS ONLY WORKS IN    QCD   CALCULATIONS
      KM_L=0.22D0
      KM_A=1.0D0
      KM_RHO=-0.4D0
      KM_ETA=0.22D0

      DO I=1,6
       DO J=1,6
        IF(I.EQ.J) THEN
         VKM(I,J)=1.0D0
        ELSE
         VKM(I,J)=0.0D0
        ENDIF
       ENDDO
      ENDDO


CCPY APRIL 2004
C      VKM(1,2)=(1.0D0-KM_L**2/2.0D0)**2
C      VKM(1,3)=(KM_L)**2
C      VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
C      VKM(4,2)=VKM(1,3)
C      VKM(4,3)=VKM(1,2)
C      VKM(4,5)=(KM_A*KM_L**2)**2
C      VKM(2,1)=VKM(1,2)
C      VKM(3,1)=VKM(1,3)
C      VKM(5,1)=VKM(1,5)
C      VKM(2,4)=VKM(4,2)
C      VKM(3,4)=VKM(4,3)
C      VKM(5,4)=VKM(4,5)

      VKM(1,2)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0)**2
      VKM(1,3)=(KM_L)**2
      VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
      VKM(4,2)=(-KM_L+0.5*KM_A**2*KM_L**5*
     >          DSQRT((1.0-2.0*KM_RHO)**2+KM_ETA**2))**2
      VKM(4,3)=(1.0D0-KM_L**2/2.0D0-KM_L**4/8.0D0*
     >         (1.0+4.0*KM_A**2))**2
      VKM(4,5)=(KM_A*KM_L**2)**2
      VKM(2,1)=VKM(1,2)
      VKM(3,1)=VKM(1,3)
      VKM(5,1)=VKM(1,5)
      VKM(2,4)=VKM(4,2)
      VKM(3,4)=VKM(4,3)
      VKM(5,4)=VKM(4,5)

      IF(TYPE_V.EQ.'A0') THEN
C COUPLINGS OF A0-F-F.  (g*Sin_w)*(\gamma_\mu)*(electric charge)
        ZFF2(1)=(2.0/3.0)**2
        ZFF2(2)=(-1.0/3.0)**2
        ZFF2(3)=(-1.0/3.0)**2
        ZFF2(4)=(2.0/3.0)**2
        ZFF2(5)=(-1.0/3.0)**2
        ZFF2(6)=(2.0/3.0)**2
      ELSEIF(TYPE_V.EQ.'AA') THEN
        ZFF2(1)=(2.0/3.0)**4
        ZFF2(2)=(-1.0/3.0)**4
        ZFF2(3)=(-1.0/3.0)**4
        ZFF2(4)=(2.0/3.0)**4
        ZFF2(5)=(-1.0/3.0)**4
        ZFF2(6)=(2.0/3.0)**4
      ELSEIF(TYPE_V.EQ.'ZZ') THEN
CsB_____See: Ohnemus, Owens PRD43 (91) 3627 (3)
C       We write Qq^4 = g+^4/2 + g-^4/2 with g+ = g- = Qq for qQ->AA,
C       and replace g+ and g- for ZZ.
        Sw_Eff = Sqrt(SW2_EFF)
        Cw_Eff = Sqrt(1.d0 - SW2_EFF)
        Tw_Eff = Sw_Eff/Cw_Eff
        T3u =  1.d0/2.d0
        T3d = -1.d0/2.d0
        ZFF2(1) = (T3u/Sw_Eff/Cw_Eff-Qu*Tw_Eff)**4/2.d0 +
     +            (-Qu*Tw_Eff)**4/2.d0
        ZFF2(2) = (T3d/Sw_Eff/Cw_Eff-Qd*Tw_Eff)**4/2.d0 +
     +            (-Qd*Tw_Eff)**4/2.d0
        ZFF2(3)=ZFF2(2)
        ZFF2(4)=ZFF2(1)
        ZFF2(5)=ZFF2(2)
        ZFF2(6)=ZFF2(1)
      ELSEIF(TYPE_V.EQ.'WW_UUB'.or.TYPE_V.EQ.'WW_DDB') THEN
        DO I=1,6
          ZFF2(I)=1d0/2d0
        ENDDO

      ELSEIF(TYPE_V.EQ.'Z0') THEN
C COUPLINGS OF Z0-F-F.  -(g/2/Cos_w)*(\gamma_\mu)*(V^I-A^I*\gamma_5)
cpn Couplings depend on lepasy 
ccpy 
       DO I=1,6
        IF(I.EQ.1 .OR. I.EQ.4 .OR. I.EQ.6) THEN
         ZFF_A(I)=1.0D0/2.0D0
         ZFF_V(I)=1.0D0/2.0D0-4.0D0/3.0D0*SW2_EFF
        ELSE
         ZFF_A(I)=-1.0D0/2.0D0
         ZFF_V(I)=-1.0D0/2.0D0+2.0D0/3.0D0*SW2_EFF
        ENDIF
        IF(LEPASY.EQ.1) THEN
          ZFF2(I)=2.0D0*ZFF_A(I)*ZFF_V(I)
        ELSE
          ZFF2(I)=ZFF_A(I)**2+ZFF_V(I)**2
        ENDIF
       ENDDO
ccpy      
CJI Sept 2013: Added in ZU and ZD to match c++ code 
      ELSEIF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     >     .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') then
        DO I=1,6
          ZFF2(I)=1.0D0
        ENDDO
C
      ELSEIF(TYPE_V.EQ.'GL') THEN
C COUPLINGS OF GL-F-F.  (g_strong)*(\gamma_\mu)*(color matrix)
        DO I=1,6
          ZFF2(I)=1.0D0
        ENDDO
      ENDIF

      IF(DEBUG) THEN
       PRINT*,'GMU,SW2,CW2,SW2_EFF,EW_A,EW_ALFA,EW_DELA,EW_DELR,'
       Print*,' EW_DELKA,EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW '
       PRINT*,GMU,SW2,CW2,SW2_EFF,EW_A,EW_ALFA,EW_DELA,EW_DELR
       Print*,EW_DELKA,EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW
       PRINT*,'MW,MZ,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL'
       PRINT*,MW,MZ,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL
       PRINT*,'ZFF2'
       PRINT*,ZFF2
       PRINT*, 'VKM(1,2),VKM(1,3),VKM(1,5),SUM_1'
       PRINT*,VKM(1,2),VKM(1,3),VKM(1,5), 
     >              VKM(1,2)+VKM(1,3)+VKM(1,5)
       PRINT*,'VKM(4,2),VKM(4,3),VKM(4,5),SUM_4'
       PRINT*,VKM(4,2),VKM(4,3),VKM(4,5), 
     >              VKM(4,2)+VKM(4,3)+VKM(4,5)
      ENDIF

      RETURN
      END

C **************************************************************************

      SUBROUTINE STANDARD_OLD
C SET UP PARAMETERS FOR THE ELECTROWEAK SECTOR
C USE (EW_ALFA,GMU,MZ) SCHEME
C FLAVOR LABEL:
C UP=1, DN=2, ST=3, CH=4 BT=5 TP=6, GL=0
C UB=-1, DB=-2, SB=-3, CB=-4 BB=-5 TB=-6

      IMPLICIT NONE
      INCLUDE 'common.for'
      integer i,j

      REAL*8 KM_L,KM_A,KM_RHO,KM_ETA
      REAL*8 CW2,EW_A,EW_DELA,EW_DELR,
     >EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW

      PI=4.0D0*ATAN(1.0D0)
      PI2=PI**2

C UNIT  CONVERT ENERGY FROM  GeV  TO  TeV, THEN UNIT=1.0D-3
c WARNING!!! You must keep UNIT=1.0D0 in this version of the code!!!!!!
      UNIT=1.0D0
C CONVERT FROM  GeV  TO pb
      HBARC2=UNIT**2*0.38937966D0*1.0D9

C INPUT DATA
      GMU=1.166389D-5/UNIT**2
      EW_ALFA=1.0D0/137.0359895D0

      MZ2=MZ**2
      R2=DSQRT(2.0D0)
      EW_DELA=0.0602+40.0D0/9.0D0*EW_ALFA/PI*DLOG(MZ/92.0D0/UNIT)

C EW_ALFA AT  MZ
      EW_AMZ=EW_ALFA/(1.0D0-EW_DELA)
      EW_A=PI*EW_ALFA/R2/GMU*UNIT**2
      EW_DELRHO=3.0D0*GMU*MT**2/8.0D0/R2/PI2
      EW_RHO=1.0D0/(1.0D0-EW_DELRHO)
C WE IGNORE THE CONTRIBUTION FROM    EW_DELRHO  IN  EW_DELR
      EW_DELR=EW_DELA
      MW2=MZ2*0.5D0*(1.0D0+DSQRT(1.0D0-4.0D0*EW_A/MZ2/(1.0D0-EW_DELR)))
      IF(MW.LT.75.0D0 .OR. MW.GT.85.0D0) THEN
       MW=DSQRT(MW2)
      ENDIF
      SW2=1.0D0-MW2/MZ2
      SW=DSQRT(SW2)
      CW2=1.0D0-SW2
      CW=DSQRT(CW2)
      SW2_EFF=0.5D0*(1.0D0-DSQRT(1.0D0-4.0D0*EW_A/MZ2
     >/EW_RHO/(1.0D0-EW_DELA)))

C WCOUPL IS (g**2)/8
      WCOUPL=R2*MW2*GMU/2.0D0
C ZCOUPL IS (g/Cos_w)**2/4
      ZCOUPL=R2*MZ2*GMU
C HCOUPL IS FOR HIGGS
      HCOUPL=R2*GMU*PI/576.0D0

C (KM MATRIX ELEMENT)**2
C TOP QUARK CHANNEL IS SWITCHED OFF
C THIS ONLY WORKS IN    QCD   CALCULATIONS
      KM_L=0.22D0
      KM_A=1.0D0
      KM_RHO=-0.4D0
      KM_ETA=0.22D0

      DO I=1,6
       DO J=1,6
        IF(I.EQ.J) THEN
         VKM(I,J)=1.0D0
        ELSE
         VKM(I,J)=0.0D0
        ENDIF
       ENDDO
      ENDDO

      VKM(1,2)=(1.0D0-KM_L**2/2.0D0)**2
      VKM(1,3)=(KM_L)**2
      VKM(1,5)=(KM_A*KM_L**3)**2*(KM_RHO**2+KM_ETA**2)
      VKM(4,2)=VKM(1,3)
      VKM(4,3)=VKM(1,2)
      VKM(4,5)=(KM_A*KM_L**2)**2
      VKM(2,1)=VKM(1,2)
      VKM(3,1)=VKM(1,3)
      VKM(5,1)=VKM(1,5)
      VKM(2,4)=VKM(4,2)
      VKM(3,4)=VKM(4,3)
      VKM(5,4)=VKM(4,5)

      IF(TYPE_V.EQ.'A0') THEN
C COUPLINGS OF A0-F-F.  (g*Sin_w)*(\gamma_\mu)*(electric charge)
        ZFF2(1)=(2.0/3.0)**2
        ZFF2(2)=(-1.0/3.0)**2
        ZFF2(3)=(-1.0/3.0)**2
        ZFF2(4)=(2.0/3.0)**2
        ZFF2(5)=(-1.0/3.0)**2
        ZFF2(6)=(2.0/3.0)**2
      ELSEIF(TYPE_V.EQ.'AA') THEN
        ZFF2(1)=(2.0/3.0)**4
        ZFF2(2)=(-1.0/3.0)**4
        ZFF2(3)=(-1.0/3.0)**4
        ZFF2(4)=(2.0/3.0)**4
        ZFF2(5)=(-1.0/3.0)**4
        ZFF2(6)=(2.0/3.0)**4
      ELSEIF(TYPE_V.EQ.'Z0') THEN
C COUPLINGS OF Z0-F-F.  -(g/2/Cos_w)*(\gamma_\mu)*(V^I-A^I*\gamma_5)
       DO I=1,6
        IF(I.EQ.1 .OR. I.EQ.4 .OR. I.EQ.6) THEN
         ZFF_A(I)=1.0D0/2.0D0
         ZFF_V(I)=1.0D0/2.0D0-4.0D0/3.0D0*SW2_EFF
        ELSE
         ZFF_A(I)=-1.0D0/2.0D0
         ZFF_V(I)=-1.0D0/2.0D0+2.0D0/3.0D0*SW2_EFF
        ENDIF
        IF(LEPASY.EQ.1) THEN
          ZFF2(I)=2.0D0*ZFF_A(I)*ZFF_V(I)
        ELSE
          ZFF2(I)=ZFF_A(I)**2+ZFF_V(I)**2
        ENDIF
       ENDDO
      ELSEIF(TYPE_V.EQ.'GL') THEN
C COUPLINGS OF GL-F-F.  (g_strong)*(\gamma_\mu)*(color matrix)
        DO I=1,6
          ZFF2(I)=1.0D0
        ENDDO
      ENDIF

      IF(DEBUG) THEN
       PRINT*,' GMU,SW2,CW2,SW2_EFF,EW_A,EW_ALFA,EW_DELA,EW_DELR,'
       Print*,' EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW'
       PRINT*, GMU,SW2,CW2,SW2_EFF,EW_A,EW_ALFA,EW_DELA,EW_DELR
       Print*, EW_RHO,EW_DELRHO,R2,EW_AMZ,SW,CW
       PRINT*,'MW,MZ,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL'
       PRINT*,MW,MZ,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL
       PRINT*,'ZFF2'
       PRINT*,ZFF2
      ENDIF

      RETURN
      END

C **************************************************************************

      SUBROUTINE QCD (nOrder)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON/SETUP1/ DEBUG,NIN,NOUT,NWRT
      LOGICAL DEBUG !,ASK

      IF(DEBUG) THEN
       PRINT*,' '
       CALL QCDOUT(NOUT)
      ENDIF

C
C     ALPI (AMU)        : Alpha / Pi             at scale  Amu
C     NFL  (AMU)        : # of effective flavors at scale  Amu
C     ALAMBD (AMU)      : Effective Lamda        at scale  Amu
C     AMUMIN () : Minimum value of Mu allowed by perturbative theory
C
C     ALEPI (AMU, NEF): Alpha / Pi at scale Amu using effective lamda of
C                 NEF # of light quarks (rather than NFL(AMU))
C                 (appropriate for renormalization sch. with fixed NEF)
C     ALPIOR (AMU, NORD): Alpha / Pi, using formula with NORD loops
C                                                at scale  Amu
C     ALAMF (I) : Effective Lamda for I "light quark" flavors
C     AMHATF(I) : Mass Threshold for  I "light quark" flavors
C     AMASS (I) : Mass of                quark flavor I
C     CH    (i)       : Charge of              quark flavor I
C     NFLTOT () : Total # of quark flavors

!      IF(DEBUG) THEN

       AMU=91.18d0

C10    PRINT*,' '
C      PRINT*,' AMU ='
C      READ*,AMU

      PRINT*,'AMU,ALPI(AMU),NFL(AMU),ALAMBD(AMU),AMUMIN(),NFLTOT()'
      PRINT*,AMU,ALPI(AMU),NFL(AMU),ALAMBD(AMU),AMUMIN(),NFLTOT()

c      PRINT*,'ALPI(AMU),ALEPI(AMU,NFL(AMU)),ALPIOR(AMU,1),ALPIOR(AMU,2)'
c      PRINT*,ALPI(AMU),ALEPI(AMU,NFL(AMU)),ALPIOR(AMU,1),ALPIOR(AMU,2)



c      PRINT*,'ALEPI(AMU,5),ALAMF(5),AMHATF(5),AMASS(5),CH(5)'
c      PRINT*,ALEPI(AMU,5),ALAMF(5),AMHATF(5),AMASS(5),CH(5)

c      PRINT*,'ALEPI(AMU,4),ALAMF(4),AMHATF(4),AMASS(4),CH(4)'
c      PRINT*,ALEPI(AMU,4),ALAMF(4),AMHATF(4),AMASS(4),CH(4)

C      IF(ASK('MORE')) GOTO 10

!      ENDIF

      RETURN
      END


C THIS IS RESUM.FOR
CCPY JULY 1993
C ADD DRELL-YAN (PHOTON) PROCESS
      subroutine sresum(vmas,vpt,rapin,gees,q0in,fresum,pma,pert,asy,
     >                  ier_res)
c      real*8 function fresum(vmas,vpt,rapin,gees,q0in)
      IMPLICIT NONE
      INCLUDE 'common.for'
      logical DUMB

      INTEGER NF_EFF
      REAL*8 RES_1,RESUM_1,RESUM_2,RESUM,B_E
      DIMENSION RES_1(1:22)
      LOGICAL FIT3_TEMP,SPEED

      EXTERNAL ZJSUDNPT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ2NT
      integer ier_res
      INTEGER IRUN,NRUN,i
      real*8 vmas,vpt,rapin,gees(3),q0in,fresum,pma,pert,asy
      Real*8 Coup_M, Omega
      Double Precision RUN_MASS_TOT
      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC

      Integer KinCorr
      Common / CorrectKinematics / KinCorr

      real*8 pyalem
      external pyalem

      real*8 Z_beta0, Z_beta1

CJI Jan 2015: Add in hj calculation
      Integer HasJet
      REAL*8 D1s, R, t, ptj, u
      REAL*8 H0Approx, H2, s_alpi
      External H0Approx
      Common /Jet/ u,D1s, R, t, ptj, HasJet
      REAL*8 KK,dsigmadt,Kgg,tau,GHJ2

CJI September 2016, add in the choice of scheme
      CHARACTER*3 SCHEME
      COMMON / RESUMSCHEME / SCHEME

      ier_res=0

CCPY Error: April 14, 2004 (correct ALEPASY for A0, etc)
      CALL SET_ZFF2
      
      
C S_ INDICATES THE PARAMETERS USED IN CALCULATING SUDAKOV FORM FACTORS

      Q_V=vmas
      S_Q=Q_V

      S_LAM=ALAMBD(S_Q)
      NF_EFF=NFL(S_Q)
c this setup depends on NF_EFF and thus S_Q, and so must be done repetively
      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR.
     .   TYPE_V.EQ.'ZG' .OR. TYPE_V.EQ.'GG') THEN ! g g initial state
       CALL HSETUP(NF_EFF)
      ELSE IF(HasJet.Eq.1) THEN
       CALL VJSETUP(NF_EFF)
      ELSE                                        ! q Q initial state
       CALL WSETUP(NF_EFF)
      ENDIF
c
      S_BETA1=(33.0-2.0*NF_EFF)/12.0

c the "*UNIT", intended to simply unit conversion and to make a useful
c   check doesn't work anymore because of the program restructuring.
c   UNIT has been set to unity in the module WHATIS
c      ECM=ECM*UNIT
      ECM=ECM
      ECM2=ECM**2
c
C
CCPY
      IF (Type_V.Eq.'W+_RA_UDB' .or. Type_V.Eq.'W-_RA_DUB' 
     > .or. Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
CJI September 2013: added new names for Z and W to match c++ version
     > .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD'
     > .or. Type_V.Eq.'WU' .or. Type_V.Eq.'WD'
     > .or. Type_V.Eq.'ZJU' .or. Type_V.Eq.'ZJD')THEN
        NO_SIGMA0 = 1
      ENDIF
        
      IF (NO_SIGMA0.EQ.1) THEN
         CONST=1.0D0/2.0D0
      ELSE
        IF(TYPE_V.EQ.'W+'.OR. TYPE_V.EQ.'HP') THEN
C THIS IS THE \sigma_0 IN MY NOTES
          CONST=PI*WCOUPL*2.0/3.0

        ELSEIF(TYPE_V.EQ.'W-'.OR. TYPE_V.EQ.'HM') THEN
C THIS IS THE \sigma_0 IN MY NOTES
          CONST=PI*WCOUPL*2.0/3.0

        ELSEIF(TYPE_V.EQ.'Z0' .OR. TYPE_V.EQ.'HZ') THEN
C THIS IS THE \sigma_0 IN MY NOTES
          CONST=PI*ZCOUPL/3.0

        ELSEIF(TYPE_V.EQ.'H0') THEN
C THIS IS THE (\sigma_0)*(Q_V)^2 IN MY NOTES
!ZL introduce the compensate scale dependence from ALPI
          Z_beta0=4*S_BETA1
          Z_beta1=34*3**2/3.d0-2*NF_EFF*(4/3.d0+5*3/3.d0)
          CONST=HCOUPL*(ALPI(C2*S_Q)*Q_V)**2
     &     *(1+Z_beta0/2d0*ALPI(S_Q)*log(C2**2)+1d0/16d0*ALPI(S_Q)**2
     &       *(Z_beta0**2*log(C2**2)**2-2*Z_beta1*log(C2**2)))

        Else If (TYPE_V.EQ.'H+') then
          CONST=Pi*HpCOUPL*2.0/3.0 ! = Pi/3 CR^2/4
CsB       Moved top mass factor here from Subroutine Standard, since Q_V is not
C         defined there.
          IF (I_MODEL.EQ.1) THEN  ! TOPCOLOR MODEL
            IF(I_RUNMASS.EQ.1) THEN
              COUP_M=RUN_MASS_TOT(Q_V,6)
            ELSE
              COUP_M=MT
            ENDIF
          ELSE IF (I_MODEL.EQ.2) THEN ! FOR 2HDMIII:
            IF(I_RUNMASS.EQ.1) THEN
              COUP_M=Sqrt(RUN_MASS_TOT(Q_V,6)*RUN_MASS_TOT(Q_V,4))
            ELSE
              COUP_M=Sqrt(MT*qmas_c)
            ENDIF
          End If
          CF = 4.d0/3.d0
          CONST=CONST*COUP_M**2

        Else If (TYPE_V.EQ.'HB') then
          CONST=Pi*HbCOUPL*2.0/3.0 ! = Pi/3 CR^2/4
CsB       Moved top mass factor here from Subroutine Standard, since Q_V is not
C         defined there.
          IF (I_RUNMASS.EQ.1) THEN
            COUP_M=RUN_MASS_TOT(Q_V,5)
          ELSE
            COUP_M=XMBOT
          End If
          CF = 4.d0/3.d0

C          If (LTOpt.GT.-1.AND.I_RUNMASS.EQ.1) then
CCPY Feb 10, 2005, eliminate "omega" in HB process which 
C is not needed if the input
C B quark mass is the MS-bar running mass at M_B scale. 
C In C(1) function, we do not include the 3.d0*DLog((Q_V/XMBOT)**2)
C term of Omega, where 
C  Omega = 3.d0*DLog((Q_V/XMBOT)**2) + 4.d0.
C The 2nd term, 4, in Omega is due to the conversion from pole mass to MS-bar
C running mass. Thus, it is not correct to include the following line
c COUP_M = COUP_M*( 1.d0 + CF*Alpi(Q_V)/4.d0*Omega )
C when a MS-bar running mass is used.
C          End If

          CONST=CONST*COUP_M**2
!          Print*, ' CONST,HbCOUPL = ', CONST,HbCOUPL
C          Print*, ' Q_V,COUP_M = ', Q_V,COUP_M

        ELSEIF(TYPE_V.EQ.'A0') THEN
C THIS IS THE \sigma_0 IN MY NOTES
C ASSUME EW_ALFA IS NOT RUNNING
C          EW_ALFA=1.0D0/137.0359895D0
C          CONST=EW_ALFA
C          CONST=(CONST**2)*4.0*PI/9.0/Q_V/Q_V
          CONST=PI**2*pyalem(q_v**2)*4.0/3.0

CCPY For fitting or the plotter, we also include the factor for the 
C decay of gamma^* into a lepton pair
	  CONST=CONST*pyalem(q_v**2)/3.0/pi/q_v/q_v;
          CONST=PI*ACOUPL*2.0/3.0 ! yfu test 2024.10.13 (override)
        ELSEIF(TYPE_V.EQ.'AA') THEN
C From notes 'How to do qqB -> gamma gamma...' pp1:
C Const = sigma_0 = 4/3 Pi^2 alpha
C With alpha = e^2/(4 Pi)
          CONST = 4.d0/3.d0*PI**2*pyalem(q_v**2)

        ELSEIF(TYPE_V.EQ.'ZZ') THEN
CsB________From Ohnemus, Owens PRD43 (91) 3627 (3): the overall factor is the
C          same for qQ->AA and qQ->ZZ, only g+ and g- differ.
CZL          CONST = 4.d0/3.d0*PI**2*pyalem(q_v**2)
          CONST = 4.d0/3.d0*PI**2*ACOUPL/(2d0*PI)

        ELSEIF(TYPE_V.EQ.'WW_UUB'.or.TYPE_V.EQ.'WW_DDB') THEN
          CONST = 4.d0/3.d0*PI**2*ACOUPL/(2d0*PI)
CZL          CONST = 4.d0/3.d0*PI**2*pyalem(q_v**2)

        ELSEIF(TYPE_V.EQ.'GL') THEN
C THIS IS THE \sigma_0 IN MY NOTES
           CONST=8.0/9.0*PI**3*ALPI(C2*S_Q)
!ZL introduce the compensate scale dependence from ALPI
          Z_beta0=4*S_BETA1
          Z_beta1=34*3**2/3.d0-2*NF_EFF*(4/3.d0+5*3/3.d0)
          CONST=CONST
     &     *(1+Z_beta0/2d0*ALPI(S_Q)*log(C2**2)+1d0/16d0*ALPI(S_Q)**2
     &       *(Z_beta0**2*log(C2**2)**2-2*Z_beta1*log(C2**2)))

        ELSEIF(TYPE_V.EQ.'AG') THEN
C THIS IS THE \sigma_0 IN MY NOTES
           CONST=2*PI/8.0/8.0*PI**2*ALPI(C2*S_Q)**2
!ZL introduce the compensate scale dependence from ALPI
          Z_beta0=4*S_BETA1
          Z_beta1=34*3**2/3.d0-2*NF_EFF*(4/3.d0+5*3/3.d0)
!          CONST=CONST
!     &     *(1+Z_beta0/2d0*ALPI(S_Q)*log(C2**2)+1d0/16d0*ALPI(S_Q)**2
!     &       *(Z_beta0**2*log(C2**2)**2-2*Z_beta1*log(C2**2)))

        ELSEIF(TYPE_V.EQ.'ZG') THEN
CsB The overall constant is the same as for GG -> AA
           CONST=2*PI/8.0/8.0*PI**2*ALPI(C2*S_Q)**2
!ZL introduce the compensate scale dependence from ALPI
          Z_beta0=4*S_BETA1
          Z_beta1=34*3**2/3.d0-2*NF_EFF*(4/3.d0+5*3/3.d0)
          CONST=CONST
     &     *(1+Z_beta0/2d0*ALPI(S_Q)*log(C2**2)+1d0/16d0*ALPI(S_Q)**2
     &       *(Z_beta0**2*log(C2**2)**2-2*Z_beta1*log(C2**2)))

        ELSEIF(TYPE_V.EQ.'GG') THEN
CSM
           CONST=PI**3/16.0*ALPI(C2*S_Q)
!ZL introduce the compensate scale dependence from ALPI
          Z_beta0=4*S_BETA1
          Z_beta1=34*3**2/3.d0-2*NF_EFF*(4/3.d0+5*3/3.d0)
          CONST=CONST
     &     *(1+Z_beta0/2d0*ALPI(S_Q)*log(C2**2)+1d0/16d0*ALPI(S_Q)**2
     &       *(Z_beta0**2*log(C2**2)**2-2*Z_beta1*log(C2**2)))
        ELSEif(TYPE_V.eq.'HJ') THEN
          KK = 1d0! + ALPI(S_Q)*11d0/4d0
          tau = 4d0 * mt**2/mh**2
          GHJ2 = asin(1d0/sqrt(tau))**2
          GHJ2 = 1d0+(1d0-tau)*GHJ2
          GHJ2 = 4d0*dsqrt(2d0)*(ALPI(125d0)/4d0)**2*GMU*4.0/9.0!tau**2*GHJ2**2
        ELSE
          PRINT*,' WRONG TYPE_V'
          CALL QUIT
        ENDIF
      ENDIF

      DUMB=.false.
      if(DUMB) then
      WRITE(NOUT,*) ' ECM,IBEAM,LEPASY'
      WRITE(NOUT,*) ECM,IBEAM,LEPASY
      WRITE(NOUT,*) ' ISET,INONPERT,IFLAG_C3,C1,C2,C3,C4'
      WRITE(NOUT,*) ISET,INONPERT,IFLAG_C3,C1,C2,C3,C4
      WRITE(NOUT,*) ' MT,MW,MZ,MH,MA'
      WRITE(NOUT,*) MT,MW,MZ,MH,MA
      WRITE(NOUT,*) ' TYPE_V'
      WRITE(NOUT,'(2X,A2)') TYPE_V
      WRITE(NOUT,*)' S_Q,S_LAM,S_BETA1,B0,A1,A2,B1,B2,NF_EFF,NORDER'
      WRITE(NOUT,*)S_Q,S_LAM,S_BETA1,B0,A1,A2,B1,B2,NF_EFF,NORDER
      WRITE(NOUT,*) ' Q_V'
      WRITE(NOUT,*) Q_V
      WRITE(NOUT,*) ' '
      endif

      nrun=1
      DO 900 IRUN=1,NRUN

c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c enters the values of the parameters for the fit AFTER the setup routine
	g1=gees(1)
	g2=gees(2)
	g3=gees(3)

	qt_v=vpt
	y_v=rapin

        q0=q0in
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      RESUM=0.D0
      RESUM_1=0.D0
      RESUM_2=0.D0
      I=0
      CALL YMAXIMUM
      IF(DABS(Y_V).GT.YMAX) THEN
C       WRITE(NWRT,*) 'ABS(Y_V)  > YMAX'
C       WRITE(NWRT,*)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
CSM        WRITE(NOUT,*) 'ABS(Y_V)  > YMAX'
CSM        WRITE(NOUT,910)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
       ier_res=1
       GOTO 900
      elseIF(abs(abs(ymax)-ABS(Y_V)).LT.YMAX/100) THEN
c if rapidity is close to YMAX, then we are in a region where the program
c   is not accurate; we flag this situation
       ier_res=2
      ENDIF

C SET UP FLAG TO TURN ON HIGH SPEED
C      SPEED=.TRUE.
      SPEED=.FALSE.
      IF(SPEED) THEN
C TO SPEED UP CODE
C WE HAVE CHECKED THIS FOR W+ PRODUCTION, IT WORKS VERY WELL
       ACC_RES1=0.05
       YNPT_CUT=15.0
       NRES1=2
       IF(QT_V.LT.25.0) THEN
        ACC_RERR=1.D-2
        FIT3=.FALSE.
       ELSE
        ACC_RERR=1.D-3
        FIT3=.FALSE.
       ENDIF
       IF(QT_V.LT.70.0 .AND. QT_V.GT.50.0) THEN
        NRES1=7
       ELSE
        NRES1=4
       ENDIF
       IF(QT_V.LT.1.0) THEN
        YNPT_CUT=35.0
       ELSEIF(QT_V.LT.4.0) THEN
        YNPT_CUT=25.0
       ELSE
        YNPT_CUT=15.0
       ENDIF
      ENDIF

      FIRST_YNPT=.TRUE.

      IF(DEBUG) THEN
       WRITE(NWRT,*)'YNPT_CUT,ACC_RERR,ACC_SUDNPT,ACC_RES1,NRES1'
       WRITE(NWRT,*)YNPT_CUT,ACC_RERR,ACC_SUDNPT,ACC_RES1,NRES1
       IF(FIT3) THEN
        WRITE(NWRT,*) 'FIT3 = TRUE'
       ELSE
        WRITE(NWRT,*) 'FIT3 = FALSE'
       ENDIF
      ENDIF

CsB___Definition of momentum fractions x1 and x2
      If (KinCorr.Eq.0) then
        X_A=Q_V/ECM*DEXP(Y_V)
        X_B=Q_V/ECM/DEXP(Y_V)
      Else If (KinCorr.Eq.1) then
        X_A=Sqrt(Q_V**2+QT_V**2)/ECM*DEXP(Y_V)
        X_B=Sqrt(Q_V**2+QT_V**2)/ECM/DEXP(Y_V)
      End If

cpn  
      if (x_a.ge.1.or.x_b.ge.1) then
        fresum = 0.0
        return
      endif 

CBY, May,2019
      if(QT_V .ge. 2.0*Q_V) then
        fresum=0.0
        return
      endif

      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2

      FIT3_TEMP=FIT3
      FIT3=.FALSE.

      IF(DEBUG) THEN
       WRITE(NWRT,*)'I,RES_1(I),RESUM_1'
      ENDIF

      I=1
111   DNLIM=ZEROJ0(I-1)
      UPLIM=ZEROJ0(I)

       RES_1(I)=ADZ2NT(ZJSUDNPT,DNLIM,UPLIM,AERR,RERR,ERREST,
     > IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN SUDNPT'
CCPY        CALL QUIT
          RES_1(I)=0.D0
       ENDIF

       RESUM_1=RESUM_1+RES_1(I)

       IF(DEBUG) THEN
        WRITE(NWRT,*)I,RES_1(I),RESUM_1
       ENDIF

       IF(I.EQ.1) THEN
        I=I+1
        GOTO 111
       ENDIF

       IF(RESUM_1.EQ.0.0) THEN
c        WRITE(NWRT,*) ' RESUM_1=0 IN RESUM.FOR'
        WRITE(NWRT,*)' WARNING: RESUM_1=0 IN RESUM.FOR; SET RESUM=0.'
        WRITE(NWRT,*) ' I,RES_1(I),RESUM_1'
        WRITE(NWRT,*) I,RES_1(I),RESUM_1
c        CALL QUIT
        resum=0.d0
        ier_res=3
        goto 900
       ENDIF

CCPY: OCT 2002
C         IF(DABS(RES_1(I)/RESUM_1).LT.ACC_SUDNPT) GOTO 210
       IF(DABS(RES_1(I)/RESUM_1).LT.ACC_SUDNPT) THEN
         IF(I.LT.2) THEN
           GOTO 210
         ELSE
           GOTO 200
         ENDIF 
       ENDIF
       
CCPY July 2019:       IF(I.EQ.25) THEN
CCPY Oct 2019: Changing the value of NRES1 could change the result of calculation
       IF(I .EQ. NRES1) THEN
C        WRITE(NWRT,*)' ACC_SUDNPT IS TOO SMALL: '
C        WRITE(NWRT,*)' RES_1(I),RESUM_1,ACC_SUDNPT'
C        WRITE(NWRT,*) RES_1(I),RESUM_1,ACC_SUDNPT
C        WRITE(NWRT,*) 'I,NRES1,QT_V,Q_V,y_V'
C        WRITE(NWRT,*) I,NRES1,QT_V,Q_V,y_V
        GOTO 200
       ENDIF

        I=I+1
        GOTO 111

200   CONTINUE

      FIT3=FIT3_TEMP

      B_E=ZEROJ0(I)/QT_V
      CALL RESUM_ASY(B_E,RESUM_2)


210   CONTINUE

      RESUM_1=RESUM_1/(2.0*QT_V**2)
      RESUM=RESUM_1+RESUM_2
C THIS IS d(SIGMA)/d(Q^2)d(QT^2)d(Y)
C      print*, resum
      RESUM=1.D0/ECM2*CONST*HBARC2*RESUM/PREV_UNDER/PREV_UNDER
C CONVERT TO d(SIGMA)/d(Q^2)d(QT)d(Y)
      RESUM=RESUM*2.0*QT_V
      if(HASJET.EQ.1) RESUM = RESUM*X_A*X_B
CJI September 2016, add in the choice of scheme
      if(SCHEME.eq.'CFG') then
        If( (TYPE_V.EQ.'W+' .or. TYPE_V.EQ.'W-'
     >  .or. TYPE_V.EQ.'Z0' .or. TYPE_V.EQ.'A0')
     >  .or. (Type_V.Eq.'W+_RA_UDB' .or. Type_V.Eq.'W-_RA_DUB'
     >     .or. Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
CJI 2013: Added in ZU and ZD to match c++ code
     >     .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') )then
            S_ALPI = ALPI(Q_V)
          H2 = CF*CA*(59*ZETA3/18-1535/192d0+215*PI2/216d0-PI2**2/240d0)
     >       + CF**2/4d0*(-15*ZETA3+511/16d0-67*PI2/12d0+17*PI2**2/45d0)
     >       + CF*NF/864d0*(192*ZETA3+1143-152*PI2)
            RESUM=RESUM*(1+0.5*CF*(Pi2-8)*S_ALPI+H2*S_ALPI**2)
        endif
      endif
CsB
C      WRITE(NWRT,*)'QT_V,Y_V,RESUM,RESUM_1,N_1,RESUM_2,X_A,X_B,YMAX'
C      WRITE(NWRT,*)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2,X_A,X_B,YMAX
C      WRITE(NWRT,*)' '

C      WRITE(NOUT,910)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
910   FORMAT(F7.2,2X,F7.2,2X,G14.7,2X,D16.4,2X,I3,2X,D16.4)

C      IF(ASK('MORE')) GOTO 10

900   CONTINUE

CsB ------------------------------------
cgal: send back d(SIGMA)/d(Q^2)d(QT)d(Y)
CsB ------------------------------------
	fresum=resum
CZL:
	call flush()

      END


C ---------------------------------------------------------------------------
      FUNCTION ZJSUDNPT(Z)
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 ZJSUDNPT
      REAL*8 Z,SUD,YNONPERT,BESJ0,S_B,CFCF
      INTEGER INPTCUT
      logical DUMB

      Integer KinCorr
      Common / CorrectKinematics / KinCorr

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

c for the fast g1g2g3
      character*40 type_hold
      integer index_b,nbloop,polyo,npoly,dipoly,iarray,ibeam_hold
      integer lepasy_hold
      parameter (nbloop=100,polyo=3)
      real divvec(4),divary(4),rarvar,divdif

C===========
CCPY SEPT 2009: This does not work.
C      real*8 bdep(nbloop),bdepfac(nbloop),b_high,b_low,s_b_in,deltab
C
C It has to be split as follows.
      REAL*8 bdep,bdepfac,b_high,b_low,s_b_in,deltab,s_b_str
      COMMON/BGRID/ bdep(nbloop),bdepfac(nbloop)
C===========

      
      real*8 arvar, cxfcxflepasy
      integer nb_grid,ifrstrn
      real*8 q_hold,qt_hold,y_hold,ecm_hold

      common/holdit/q_hold,qt_hold,y_hold,ecm_hold,ibeam_hold,type_hold
      external divdif

CCPY Sept 2009      data qt_hold/-9d9/
      data q_hold,qt_hold,y_hold,ecm_hold,ibeam_hold,type_hold
     >    /-9d9,-9d9,-9d9,-9d9,-99999,'00'/
cpn                    
      external besj0

      INTEGER I,J, lQ, lY, lB

      S_B=Z/QT_V

      CALL NONPERT(S_B,YNONPERT,INPTCUT)
      IF(INPTCUT.EQ.1) THEN
       ZJSUDNPT=0.0
       RETURN
      ENDIF

c****************************************************************************
c Purpose of this conditional:  To toggle between modifications specifically
c        designed to increase the speed of cross section tabulations when
c        different values of nonperturbative parameters are to be tabulated
c	 and the old way.
c****************************************************************************
CCPY OCT 2007, turn on Glenn SPEED-UP version when Kincorr=1. 
C OTHERWISE, use Pavel's version which only works for 
C KINCORR NOT EQUAL TO 1, BUT IFSAT CAN BE EITHER 1 OR 0. 

      Z = S_B*QT_V

CCPY June 2019
C      if(KINCORR.EQ.0) then
C        CALL SUDAKOV(S_B,SUD)
C        CALL CXFCXF(S_B,CFCF)
C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'S_B,SUD,YNONPERT,CFCF,QT_V,Z'
C      ENDIF

C        ZJSUDNPT=Z*BESJ0(Z)*SUD*YNONPERT*CFCF

C      else

c this is the g1g2g3 short-cut
ccpy Oct 2007, correct the code for implementing LEPASY=1 OPTION  
        if(q_hold.ne.q_v.or.qt_hold.ne.qt_v.or.ibeam.ne.ibeam_hold.or.
     >     type_v.ne.type_hold.or.y_hold.ne.y_v.or.ecm_hold.ne.ecm.or.
     >     lepasy_hold.ne.lepasy.or.I_RESET.NE.0) then
CCPY SEPT 2009:    >     lepasy_hold.ne.lepasy) then
          q_hold=q_V
          qt_hold=qt_V
          y_hold=y_V
          ecm_hold=ecm
          ibeam_hold=ibeam
          type_hold=type_v
          lepasy_hold=lepasy
          ifrstrn=1
          I_RESET=0
        endif
c in 1st run, create array in
c atan(nb_grid*b)/(pi/2) which has interval (0,1)
CCPY June 2019
        nb_grid=20
C        nb_grid=30
        if(ifrstrn.eq.1) then
CCPY Sept 2009: Initialization
          DO I=1,nbloop
            bdep(I)=0.d0
            bdepfac(I)=0.d0      
          ENDDO

          b_high=1.d0
          b_low =0.d0
CCPY Sept 2009:          deltab=(b_high-b_low)/(nbloop)
          deltab=(b_high-b_low)/(real(nbloop))
          arvar=-deltab/2
          do index_b=1,nbloop
            arvar=arvar+deltab
            s_b_in=tan(pi*arvar/2)/nb_grid

            bdep(index_b)=arvar
CCPY!cdump
CCPY!      print*,'index_b,arvar,bdep(index_b) =',
CCPY!     >index_b,arvar,bdep(index_b)  
CCPY!
C            S_B_str=S_B_in/(sqrt(1+S_B_in**2/bmax**2))
            CALL SUDAKOV(S_B_in,SUD)
            CALL CXFCXF(S_B_in,CFCF)
C            CALL NONPERT(S_B_in,YNONPERT,INPTCUT)
            bdepfac(index_b)=sud*cfcf
          enddo
          ifrstrn=-999

        endif
c procedure: (1) for this call make subarray about s_b value (2) interpolate
c            using something like divdif to get appropriate sud*cfcf value;
c            this procedure assumed  0<arvar<1
c            polyo=order of polynomial interpolation
c            npoly=#grid points needed above and below s_b for polyo
        rarvar=real( atan(nb_grid*s_b)/(pi/2.d0) )
        if(rarvar.lt.0.or.rarvar.gt.1.0) write(6,*) 'RARVAR=',rarvar
        index_b=int(rarvar*float(nbloop))+1
        npoly=polyo+1-int(float(polyo+1)/2.0)
        if(index_b+npoly.gt.nbloop) then
          index_b=nbloop-npoly
        elseif(index_b.lt.npoly) then
          index_b=npoly
        endif
        do iarray=1,npoly
         divary(iarray)=real(bdepfac(index_b-iarray+1))
         divary(npoly+iarray)=real(bdepfac(index_b+iarray))
         divvec(iarray)=real(bdep(index_b-iarray+1))
         divvec(npoly+iarray)=real(bdep(index_b+iarray))
        enddo
        dipoly=2*npoly
c
        ZJSUDNPT=Z*BESJ0(Z)*YNONPERT*
     >                        divdif(divary,divvec,dipoly,rarvar,polyo)

CCPY!cdump
CCPY!       WRITE(NWRT,*)'bdep(92),bdep(93),bdep(94),bdep(95),bdep(96) ='
CCPY!       WRITE(NWRT,*)bdep(92),bdep(93),bdep(94),bdep(95),bdep(96)
CCPY!       WRITE(NWRT,*)'index_b,npoly,nbloop ='
CCPY!       WRITE(NWRT,*)index_b,npoly,nbloop
CCPY!       WRITE(NWRT,*)'divary=',divary
CCPY!       WRITE(NWRT,*)'divvec=',divvec
CCPY!       WRITE(NWRT,*)'dipoly=',dipoly
CCPY!       WRITE(NWRT,*)'rarvar=',rarvar
CCPY!       WRITE(NWRT,*)'polyo=',polyo
CCPY!
CCPY!       WRITE(NWRT,*)'divdif(divary,divvec,dipoly,rarvar,polyo)'
CCPY!       WRITE(NWRT,*)divdif(divary,divvec,dipoly,rarvar,polyo)
CCPY!       
        CLOSE(999)

CCPU June 2019      endif  ! if(KINCORR.EQ.0)
c****************************************************************************
c end of toggle
c****************************************************************************

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'ZJSUDNPT = ',ZJSUDNPT
C      ENDIF
C        Print*, Z, ZJSUDNPT

      RETURN
      END ! ZJSUDNPT

c
c I'll keep the old version of ZJSUDNPT around for a while, but it may be
c    deleted at any time.
c
C ---------------------------------------------------------------------------
      real*8 FUNCTION ZJSUDNPT_OLD(Z)
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 ZJSUDNPT
      REAL*8 Z,SUD,YNONPERT,BESJ0,S_B,CFCF
      INTEGER INPTCUT
cpn                    
      external besj0

      S_B=Z/QT_V

      ZJSUDNPT_OLD=0.D0
      CALL NONPERT(S_B,YNONPERT,INPTCUT)
      IF(INPTCUT.EQ.1) THEN
       ZJSUDNPT=0.0
       RETURN
      ENDIF
      CALL SUDAKOV(S_B,SUD)
      CALL CXFCXF(S_B,CFCF)

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'S_B,SUD,YNONPERT,CFCF,QT_V,Z'
C       WRITE(NWRT,*)S_B,SUD,YNONPERT,CFCF,QT_V,Z
C      ENDIF

      ZJSUDNPT_OLD=Z*BESJ0(Z)*SUD*YNONPERT*CFCF

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'ZJSUDNPT_OLD = ',ZJSUDNPT_OLD
C      ENDIF

      RETURN
      END ! ZJSUDNPT_OLD(Z)

C ---------------------------------------------------------------------------
      SUBROUTINE SUDAKOV(B_SUD,SUD)
C ---------------------------------------------------------------------------
C THIS IS FOR   SUDAKOV FORM FACTOR

      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 SUD,B_SUD,S_BSTAR
      REAL*8 SUD_INT
      EXTERNAL SUD_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZINT
      Integer iStore
      Common / iS / iStore

      S_BSTAR=B_SUD/(DSQRT(1.0+(B_SUD/BMAX)**2))

      UPLIM=2.0D0*DLOG(C2*S_Q)
      DNLIM=2.0D0*DLOG(C1/S_BSTAR)

CsB To study the difference between the CSS and ERV formalism we force the lower
C   limit to be always smaller than the upper limit:
CC      iStore = iStore + 1
C      If (DnLim.Ge.UpLim) then
CC        Print*, iStore, DNLIM, UPLIM
C        Sud = 1.d0
C        Return
C      End If

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'B_SUD,S_BSTAR,DNLIM,UPLIM,BMAX'
C       WRITE(NWRT,*)B_SUD,S_BSTAR,DNLIM,UPLIM,BMAX
C      ENDIF

      IF(DNLIM.LT.0.0) THEN
       WRITE(NWRT,*)' DNLIM < 0'
      ENDIF

      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2

cpn Aug 19, 1999 In the new version of AdzInt in UtlPac, it is required that
cpn              uplim >=dnlim
      if (uplim.ge.dnlim) then
        SUD=ADZINT(SUD_INT,DNLIM,UPLIM,AERR,RERR,ERREST,
     >    IER,IACTA,IACTB)
      else
        SUD= -ADZINT(SUD_INT,UPLIM,DNLIM,AERR,RERR,ERREST,
     >    IER,IACTA,IACTB)
      endif  

       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN SUDAKOV'
        CALL QUIT
       ENDIF

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'B_SUD,SUD,ERREST,IER'
C       WRITE(NWRT,*)B_SUD,SUD,ERREST,IER
C      ENDIF

      SUD=DEXP(-SUD)

      RETURN
      END

C ---------------------------------------------------------------------------
      FUNCTION SUD_INT(DLMU2)
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 SUD_INT
      REAL*8 AMU,DLMU2,S_ALPI
CSM
      Integer HasJet
      REAL*8 XMT2Q2,LOGXB,BETA
      REAL*8 SS_Q
      Real*8 D1s, R,t,ptj, u
      Common /Jet/ u,D1s, R,t,ptj, HasJet

      INTEGER NF_EFF

CCPY Dec 2014:
      AMU=DSQRT(DEXP(DLMU2))
      S_ALPI=ALPI(AMU)

C=========================================
CCPY Dec 2014:
      NF_EFF = NFL(amu)

c this setup depends on NF_EFF and thus S_Q, and so must be done repetively
      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR.
     .   TYPE_V.EQ.'ZG' .OR. TYPE_V.EQ.'GG') THEN ! g g initial state
C       CALL HSETUP(NF_EFF)
      ELSE                                        ! q Q initial state
       CALL WSETUP(NF_EFF)
      ENDIF

      S_BETA1=(33.0-2.0*NF_EFF)/12.0
C=========================================

CSM
C.....FSR same for 'GL' and 'GG'
CCPY
      IF(I_FSR.EQ.1) THEN
        IF(TYPE_V.EQ.'GG') THEN
         B1=-2.*BETA1
        ELSE
         B1=-3.0D0/2.0D0*CF
        ENDIF
        XMT2Q2=(MT/Q_V)**2
        BETA=SQRT(1.0D0-4.0D0*XMT2Q2)
        IF(BETA.GT.1.D-6) THEN
           LOGXB=LOG((1.0D0-BETA)/(1.0D0+BETA))/BETA
        ELSE
           LOGXB=-2.D0
        ENDIF
        B1=B1+CF*(1.0D0+LOGXB*(1.0D0-2.D0*XMT2Q2))
      ENDIF

CCPY OCT. 17, 1995
C      IF(NORDER.EQ.1) THEN
C       SUD_INT=S_ALPI*((DLOG(S_Q**2)-DLMU2)*A1+B1)
C      ELSE IF(NORDER.EQ.2) THEN
C       SUD_INT=S_ALPI*((DLOG(S_Q**2)-DLMU2)*
C     >         (A1+S_ALPI*A2)+(B1+S_ALPI*B2))
C      ENDIF

CSM...APR 15 1996
CJI Sept 2013: Seperated the A and B Sudakov pieces to match new desired form
      SS_Q=C2*S_Q
      IF(N_SUD_A.EQ.1 .and. N_SUD_B.EQ.1) THEN
        SUD_INT=S_ALPI*((DLOG(SS_Q**2)-DLMU2)*A1+B1)
      ELSE IF(N_SUD_A.EQ.2 .and. N_SUD_B.EQ.2) THEN
        SUD_INT=S_ALPI*((DLOG(SS_Q**2)-DLMU2)*
     >          (A1+S_ALPI*A2)+(B1+S_ALPI*B2))
C        Print*, ' A2 =', A2
C        Stop
C        print*, a3, a2, a1, b2, b1
      ELSE IF(N_SUD_A.EQ.3 .and. N_SUD_B.EQ.2) THEN
        SUD_INT=S_ALPI*((DLOG(SS_Q**2)-DLMU2)*
     >          (A1+S_ALPI*A2+S_ALPI**2*A3)+(B1+S_ALPI*B2))
C        Print*, ' A3 =', A3
C        Stop
      ELSE IF(N_SUD_A.EQ.3 .and. N_SUD_B.EQ.2) THEN
        SUD_INT=S_ALPI*((DLOG(SS_Q**2)-DLMU2)*
     >          (A1+S_ALPI*A2+S_ALPI**2*A3)
     >         +(B1+S_ALPI*B2))
      ELSE IF(N_SUD_A.EQ.3 .and. N_SUD_B.EQ.3) THEN
        SUD_INT=S_ALPI*((DLOG(SS_Q**2)-DLMU2)*
     >          (A1+S_ALPI*A2+S_ALPI**2*A3)
     >         +(B1+S_ALPI*B2+S_ALPI**2*B3))
      ELSE IF(N_SUD_A.EQ.4 .and. N_SUD_B.EQ.3) THEN
          SUD_INT=S_ALPI*((DLOG(SS_Q**2)-DLMU2)*
     >          (A1+S_ALPI*A2+S_ALPI**2*A3+S_ALPI**3*A4)
     >         +(B1+S_ALPI*B2+S_ALPI**2*B3))
      ELSE
        WRITE(*,*) ' NO SUCH N_SUD_A AND N_SUD_B COMBINATION '
        CALL QUIT
      ENDIF

      if(HasJet.EQ.1) then
          SUD_INT = SUD_INT + D1s*dlog(1.0/R)*s_alpi
      ENDIF

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE NONPERT(B_NPT,YNONPERT,INPTCUT)
C --------------------------------------------------------------------------
C     THIS IS TESTING THE NON-PERTURBATIVE PIECE

      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 YNONPERT,EXPONENT,B_NPT,BSTAR,TERM
      REAL*8 X_0, Q_0, LAMBDA
      INTEGER INPTCUT
      INTEGER IKNT
      DATA IKNT/0/
      SAVE IKNT

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'S_Q,G1,G2,Q0,G3'
C       WRITE(NWRT,*)S_Q,G1,G2,Q0,G3
C      ENDIF

      IF(INONPERT.EQ.1) THEN
       EXPONENT=B_NPT**2*(G1+G2*DLOG(S_Q/2.0/Q0))
      ELSEIF(INONPERT.EQ.2) THEN
       BSTAR=B_NPT/(DSQRT(1.0+(B_NPT/BMAX)**2))
       EXPONENT=G2*4.0*A1*DLOG(C2*S_Q*BMAX/C1)*DLOG(B_NPT/BSTAR)
CSM
       IF(TYPE_V.EQ.'H0'.OR.TYPE_V.EQ.'GG') THEN
        IF(IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     >   IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) THEN
C CANONICAL CHOICE FOR H0
         TERM=BSTAR
CCPY        Else If (IFLAG_C3.EQ.5) then
        Else 
          TERM=BSTAR
        ENDIF
       ELSE
        IF(IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     $     IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) THEN
C CANONICAL CHOICE FOR W-BOSON
         TERM=BSTAR
CCPY        Else If (IFLAG_C3.EQ.5) then
        Else 
          TERM=BSTAR
        ENDIF
       ENDIF
       EXPONENT=EXPONENT*ALPI(C3/TERM)+G1*B_NPT
      ELSEIF(INONPERT.EQ.3) THEN
CCPY ADDED SEPT 1993
       EXPONENT=B_NPT*(G1*(B_NPT+G3*DLOG(100.0D0*X_A*X_B))+
     >                      G2*B_NPT*DLOG(S_Q/2.0/Q0))
      ELSEIF(INONPERT.EQ.4) THEN
CJI Modified to make g3(new) = g1*g3(old)
       EXPONENT=B_NPT**2*(G1*1.0D0+G3*DLOG(100.0D0*X_A*X_B)+
     >                      G2*DLOG(S_Q/2.0/Q0))
      ELSEIF(INONPERT.EQ.5) THEN
CCPY ADDED Feb 2004: BLNY
       EXPONENT=B_NPT**2*(G1*(1.0D0+G3*DLOG(100.0D0*X_A*X_B))
     >                      +G2*DLOG(S_Q/2.0/Q0))
      ELSEIF(INONPERT.EQ.6) THEN
CJI ADDED Jan 2014: TMD-BLNY
       BSTAR=B_NPT/(DSQRT(1.0+(B_NPT/BMAX)**2))
       X_0=0.01
       LAMBDA=0.2
       EXPONENT=G1*B_NPT**2+G2*DLOG(B_NPT/BSTAR)*DLOG(S_Q/Q0)
     >         +G3*B_NPT**2*((X_0/X_A)**LAMBDA+(X_0/X_B)**LAMBDA)
      ELSEIF(INONPERT.EQ.7) THEN
        EXPONENT=B_NPT**2*G1
      ELSEIF(INONPERT.EQ.16) THEN
CZL only prepared for W+
       EXPONENT=B_NPT**2*(G1*(1+G3*log(100d0*X_A*X_B)) 
     &   + G3*log(S_Q/3.2)
     &   + G2*(sqrt(1.0/X_A/X_A+1.0/G1/G1)-1/G1
     &        +sqrt(1.0/X_B/X_B+1.0/G1/G1)-1/G1) )
      ENDIF
      IF(EXPONENT.GT.YNPT_CUT) THEN
       IF(IKNT.EQ.0) THEN
        WRITE(NWRT,*)' SET YNONPERT = 0     FOR:'
        WRITE(NWRT,*)' QT_V,Y_V,B_NPT,EXPONENT'
        WRITE(NWRT,*)QT_V,Y_V,B_NPT,EXPONENT
        IKNT=1
       ENDIF
       INPTCUT=1
       YNONPERT=0.0
      ELSE
       INPTCUT=0
       YNONPERT=DEXP(-EXPONENT)
      ENDIF

C      IF(DEBUG) THEN
C       WRITE(*,*)'g1,g2,g3 =',g1,g2,g3
C       WRITE(NWRT,*)'B_NPT,YNONPERT,EXPONENT'
C       WRITE(NWRT,*)B_NPT,YNONPERT,EXPONENT
C      ENDIF

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE RESUM_ASY(B_E,RESUM_2)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 B_E,DBE,FIC0,FIC1,FIC2,FIS0,FIS1,FIS2
      REAL*8 T1,T2,T3,T4,Y_IC,Y_IS,RESUM_2

      REAL*8 F_IC,F_IS
      EXTERNAL F_IC,F_IS

      RESUM_2=0.D0
C      IF(DEBUG) THEN
       WRITE(NWRT,*) ' RESUM_ASY IS CALLED FOR B_E, QT_V = ',
     >B_E,QT_V
C      ENDIF
      DBE=1.0/(B_E*QT_V**2)
      IF((B_E-DBE/2.0) .LT. 0.0 .OR. (B_E-DBE) .LT. 0.0) THEN
       WRITE(NWRT,*)' ERROR IN DBE'
       CALL QUIT
      ENDIF
      FIC0=F_IC(B_E)
      FIC1=(F_IC(B_E+DBE/2.0)-F_IC(B_E-DBE/2.0))/DBE
      FIC2=(F_IC(B_E+DBE)-2.0*F_IC(B_E)+F_IC(B_E-DBE))/DBE**2
      FIS0=F_IS(B_E)
      FIS1=(F_IS(B_E+DBE/2.0)-F_IS(B_E-DBE/2.0))/DBE
      FIS2=(F_IS(B_E+DBE)-2.0*F_IS(B_E)+F_IS(B_E-DBE))/DBE**2
      T1=B_E*QT_V-PI/4.0
      T2=DCOS(T1)
      T3=DSIN(T1)
      T4=1.0/QT_V
      Y_IC=T2*(-T4**2*FIC1)+T3*(-T4*FIC0+T4**3*FIC2)
      Y_IS=T2*(T4*FIS0-T4**3*FIS2)+T3*(-T4**2*FIS1)
      RESUM_2=0.5*(Y_IC-Y_IS)

      RETURN
      END ! RESUM_ASY

C --------------------------------------------------------------------------
      SUBROUTINE F_ICIS(B_E)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 TERM1,TERM2,YNPT,YSUD,CFCF
      COMMON/FICIS/ TERM1,TERM2,YNPT,YSUD,CFCF
      REAL*8 B_E
      INTEGER INPTCUT

      CALL NONPERT(B_E,YNPT,INPTCUT)
      IF(INPTCUT.EQ.1) THEN
       YNPT=0.0
       TERM1=0.0
       TERM2=0.0
       YSUD=0.0
       CFCF=0.0
       RETURN
      ENDIF
      TERM1=1.0/(B_E*QT_V)
      TERM2=DSQRT(2.0*B_E/PI/QT_V)
      CALL SUDAKOV(B_E,YSUD)
      CALL CXFCXF(B_E,CFCF)

      RETURN
      END ! F_ICIS

C --------------------------------------------------------------------------
      FUNCTION F_IC(B_E)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 F_IC
      REAL*8 TERM1,TERM2,YNPT,YSUD,CFCF
      COMMON/FICIS/ TERM1,TERM2,YNPT,YSUD,CFCF
      REAL*8 B_E,T3

      CALL F_ICIS(B_E)
      T3=P00+P02*TERM1**2+P04*TERM1**4
      F_IC=TERM2*T3*YNPT*YSUD*CFCF

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'F_IC:   B_E,F_IC,T3,TERM1,TERM2,YNPT,YSUD,CFCF'
C       WRITE(NWRT,*)B_E,F_IC,T3,TERM1,TERM2,YNPT,YSUD,CFCF
C      ENDIF

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION F_IS(B_E)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 F_IS
      REAL*8 TERM1,TERM2,YNPT,YSUD,CFCF
      COMMON/FICIS/ TERM1,TERM2,YNPT,YSUD,CFCF
      REAL*8 B_E,T3

      CALL F_ICIS(B_E)
      T3=Q01*TERM1+Q03*TERM1**3
      F_IS=TERM2*T3*YNPT*YSUD*CFCF

C      IF(DEBUG) THEN
C       WRITE(NWRT,*)'F_IS:   B_E,F_IC,T3,TERM1,TERM2,YNPT,YSUD,CFCF'
C       WRITE(NWRT,*)B_E,F_IS,T3,TERM1,TERM2,YNPT,YSUD,CFCF
C      ENDIF

      RETURN
      END

C-------------------------------------------------------------------------
      subroutine GetCfCf(amu,cfcf)
C-------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 AMU,CONV,CFCF
      REAL*8 U_DB,DB_U,D_UB,UB_D,U_UB,D_DB,UB_U,DB_D
      REAL*8 CONVP1,CONVP2,CONVP3,CONVP4,CONVP5,
     >CONVM1,CONVM2,CONVM3,CONVM4,CONVM5
      REAL*8 CFCF_P,CFCF_N

      integer icnvx
      real*8 cnvxa(-6:6),cnvxb(-6:6)

CCPY      INTEGER I_CXFF
      REAL*8 PDF1,PDF2,PDF3,PDF4,PDF5
      REAL*8 PDFM1,PDFM2,PDFM3,PDFM4,PDFM5

      Real*8 C_BB,BB_C,B_BB,BB_B

      INTEGER NF_EFF

CPPY Feb 2016: To have the correct NF_EFF in SUDAKOV and CXFCXF
      S_LAM=ALAMBD(amu)
      NF_EFF=NFL(amu)
c this setup depends on NF_EFF and thus S_Q, and so must be done repetively
      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR.
     .     TYPE_V.EQ.'ZG' .OR. TYPE_V.EQ.'GG') THEN ! g g initial state
C        CALL HSETUP(NF_EFF)
      ELSE                      ! q Q inital state
        CALL WSETUP(NF_EFF)
      ENDIF
      
CCPY Error: April 14, 2004 (correct ALEPASY for A0, etc)
      CALL SET_ZFF2

      IF(IBEAM.NE.-1) GOTO 999
C THIS IS FOR IBEAM=-1
C FOR P PBAR

      IF(TYPE_V.EQ.'W+'.OR. TYPE_V.EQ.'HP'
     >  .OR. Type_V.Eq.'W+_RA_UDB') THEN

CCPY IF I_CXFF=1 THEN DO
C       CXF(X_A)*F(X_B) + F(X_A)*CXF(X_B)
CsB This is to return only O(alpha_s) contributions NOT O(alpha_s^2) ones.
C   To check the total cross section, LCA, ect. against the NLO codes.
CCPY       I_CXFF=0
       IF(I_CXFF.NE.1) THEN
CCPY THIS IS FOR NORMAL RUNS: DO    CXF(X_A)*CXF(X_B)

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(2,X_B,AMU)
       CONVP3=CONV(3,X_B,AMU)
       CONVP5=CONV(5,X_B,AMU)

       U_DB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-1,X_B,AMU)
       CONVM4=CONV(-4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       DB_U=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       ELSEIF(I_CXFF.EQ.1) THEN
C FOR CHECKING:  DO      CXF(X_A)*F(X_B) + F(X_A)*CXF(X_B)

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(2,X_B,AMU)
       CONVP3=CONV(3,X_B,AMU)
       CONVP5=CONV(5,X_B,AMU)

       PDF1=APDF(1,X_A,AMU)
       PDF4=APDF(4,X_A,AMU)
       PDF2=APDF(2,X_B,AMU)
       PDF3=APDF(3,X_B,AMU)
       PDF5=APDF(5,X_B,AMU)

       U_DB=VKM(1,2)*(CONVP1*PDF2+PDF1*CONVP2-PDF1*PDF2)
     >     +VKM(1,3)*(CONVP1*PDF3+PDF1*CONVP3-PDF1*PDF3)
     >     +VKM(1,5)*(CONVP1*PDF5+PDF1*CONVP5-PDF1*PDF5)
     >     +VKM(4,2)*(CONVP4*PDF2+PDF4*CONVP2-PDF4*PDF2)
     >     +VKM(4,3)*(CONVP4*PDF3+PDF4*CONVP3-PDF4*PDF3)
     >     +VKM(4,5)*(CONVP4*PDF5+PDF4*CONVP5-PDF4*PDF5)

       CONVM1=CONV(-1,X_B,AMU)
       CONVM4=CONV(-4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       PDFM1=APDF(-1,X_B,AMU)
       PDFM4=APDF(-4,X_B,AMU)
       PDFM2=APDF(-2,X_A,AMU)
       PDFM3=APDF(-3,X_A,AMU)
       PDFM5=APDF(-5,X_A,AMU)

       DB_U=VKM(1,2)*(CONVM1*PDFM2+PDFM1*CONVM2-PDFM1*PDFM2)
     >     +VKM(1,3)*(CONVM1*PDFM3+PDFM1*CONVM3-PDFM1*PDFM3)
     >     +VKM(1,5)*(CONVM1*PDFM5+PDFM1*CONVM5-PDFM1*PDFM5)
     >     +VKM(4,2)*(CONVM4*PDFM2+PDFM4*CONVM2-PDFM4*PDFM2)
     >     +VKM(4,3)*(CONVM4*PDFM3+PDFM4*CONVM3-PDFM4*PDFM3)
     >     +VKM(4,5)*(CONVM4*PDFM5+PDFM4*CONVM5-PDFM4*PDFM5)

       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=U_DB-DB_U
       ELSE
         CFCF=U_DB+DB_U
       ENDIF

      ELSE IF(TYPE_V.EQ.'W-'.OR. TYPE_V.EQ.'HM'
     > .OR. Type_V.Eq.'W-_RA_DUB') THEN

       CONVP1=CONV(1,X_B,AMU)
       CONVP4=CONV(4,X_B,AMU)
       CONVP2=CONV(2,X_A,AMU)
       CONVP3=CONV(3,X_A,AMU)
       CONVP5=CONV(5,X_A,AMU)

       D_UB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-1,X_A,AMU)
       CONVM4=CONV(-4,X_A,AMU)
       CONVM2=CONV(-2,X_B,AMU)
       CONVM3=CONV(-3,X_B,AMU)
       CONVM5=CONV(-5,X_B,AMU)

       UB_D=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=D_UB-UB_D
       ELSE
         CFCF=D_UB+UB_D
       ENDIF

      Else If (TYPE_V.EQ.'H+') then

CCPY       I_CXFF=0
       IF(I_CXFF.NE.1) THEN

         CONVP4=CONV(4,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         C_BB=CONVP4*CONVP5

         CONVM4=CONV(-4,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         BB_C=CONVM4*CONVM5

       ELSEIF(I_CXFF.EQ.1) THEN

         CONVP4=CONV(4,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         PDF4=APDF(4,X_A,AMU)
         PDF5=APDF(5,X_B,AMU)

         C_BB=CONVP4*PDF5+PDF4*CONVP5-PDF4*PDF5

         CONVM4=CONV(-4,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         PDFM4=APDF(-4,X_B,AMU)
         PDFM5=APDF(-5,X_A,AMU)

         BB_C=CONVM4*PDFM5+PDFM4*CONVM5-PDFM4*PDFM5

       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=C_BB-BB_C
       ELSE
         CFCF=C_BB+BB_C
       ENDIF

      Else If (TYPE_V.EQ.'HB') then

CCPY       I_CXFF=0
       IF(I_CXFF.NE.1) THEN

         CONVP4=CONV(5,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         B_BB=CONVP4*CONVP5

         CONVM4=CONV(-5,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         BB_B=CONVM4*CONVM5

       ELSEIF(I_CXFF.EQ.1) THEN

         CONVP4=CONV(5,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         PDF4=APDF(5,X_A,AMU)
         PDF5=APDF(5,X_B,AMU)

         B_BB=CONVP4*PDF5+PDF4*CONVP5-PDF4*PDF5

         CONVM4=CONV(-5,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         PDFM4=APDF(-5,X_B,AMU)
         PDFM5=APDF(-5,X_A,AMU)

         BB_B=CONVM4*PDFM5+PDFM4*CONVM5-PDFM4*PDFM5

       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=B_BB-BB_B
       ELSE
         CFCF=B_BB+BB_B
       ENDIF

CsB___The pP parton luminosity is the same for the following:
      ELSE IF(TYPE_V.EQ.'Z0'.OR. TYPE_V.EQ.'A0' .OR.
     > TYPE_V.EQ.'AA' .Or. TYPE_V.EQ.'ZZ' .Or.
     > TYPE_V.EQ.'HZ' .OR. TYPE_V.EQ.'GL' ) THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       u_ub=cnvxa(1)*cnvxb(1)
     >     +cnvxa(4)*cnvxb(4)

       d_db=cnvxa(2)*cnvxb(2)
     >     +cnvxa(3)*cnvxb(3)
     >     +cnvxa(5)*cnvxb(5)

       ub_u=cnvxb(-1)*cnvxa(-1)
     >     +cnvxb(-4)*cnvxa(-4)

       db_d=cnvxb(-2)*cnvxa(-2)
     >     +cnvxb(-3)*cnvxa(-3)
     >     +cnvxb(-5)*cnvxa(-5)

c
c       U_UB=CONV(1,X_A,AMU)*CONV(1,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(2,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(5,X_B,AMU)
c
c       UB_U=CONV(-1,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(-4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(-2,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(-3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(-5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)
      
       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

!ZL
      ELSE IF(Type_V.Eq.'WW_UUB' .or. Type_V.Eq.'WW_DDB') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       IF(Type_V.Eq.'WW_UUB') THEN
       u_ub=cnvxa(1)*cnvxb(1)
     >     +cnvxa(4)*cnvxb(4)

       ub_u=cnvxb(-1)*cnvxa(-1)
     >     +cnvxb(-4)*cnvxa(-4)

       d_db=0.d0
       db_d=0.d0

       ELSEIF(Type_V.Eq.'WW_DDB') THEN
       d_db=cnvxa(2)*cnvxb(2)
     >     +cnvxa(3)*cnvxb(3)
     >     +cnvxa(5)*cnvxb(5)

       db_d=cnvxb(-2)*cnvxa(-2)
     >     +cnvxb(-3)*cnvxa(-3)
     >     +cnvxb(-5)*cnvxa(-5)

       u_ub=0.d0
       ub_u=0.d0

       ENDIF

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)
      
       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF


CCPY
CJI Sept 2013: Added in ZU and ZD to match c++ code
      ELSE IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     >       .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') THEN

        do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
        enddo

        IF( Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'ZU') THEN
       u_ub=cnvxa(1)*cnvxb(1)
     >     +cnvxa(4)*cnvxb(4)

       ub_u=cnvxb(-1)*cnvxa(-1)
     >     +cnvxb(-4)*cnvxa(-4)
       
       d_db=0.d0
       db_d=0.d0
       
        ELSEIF( Type_V.Eq.'Z0_RA_DDB' .or. Type_V.Eq.'ZD') THEN

       d_db=cnvxa(2)*cnvxb(2)
     >     +cnvxa(3)*cnvxb(3)
     >     +cnvxa(5)*cnvxb(5)

       db_d=cnvxb(-2)*cnvxa(-2)
     >     +cnvxb(-3)*cnvxa(-3)
     >     +cnvxb(-5)*cnvxa(-5)

       u_ub=0.d0
       ub_u=0.d0
       
        ENDIF
c
       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF


C------------------------------
      ENDIF

      GOTO 888

999   CONTINUE
      IF(IBEAM.NE.1) GOTO 997
C THIS IS FOR IBEAM=1
C FOR P P

      IF(TYPE_V.EQ.'W+'.OR. TYPE_V.EQ.'HP'
     > .OR. Type_V.Eq.'W+_RA_UDB') THEN

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(-2,X_B,AMU)
       CONVP3=CONV(-3,X_B,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       U_DB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(1,X_B,AMU)
       CONVM4=CONV(4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       DB_U=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=U_DB-DB_U
       ELSE
         CFCF=U_DB+DB_U
       ENDIF

      ELSE IF(TYPE_V.EQ.'W-'.OR. TYPE_V.EQ.'HM'
     > .OR. Type_V.Eq.'W-_RA_DUB') THEN

       CONVP1=CONV(-1,X_B,AMU)
       CONVP4=CONV(-4,X_B,AMU)
       CONVP2=CONV(2,X_A,AMU)
       CONVP3=CONV(3,X_A,AMU)
       CONVP5=CONV(5,X_A,AMU)

       D_UB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-1,X_A,AMU)
       CONVM4=CONV(-4,X_A,AMU)
       CONVM2=CONV(2,X_B,AMU)
       CONVM3=CONV(3,X_B,AMU)
       CONVM5=CONV(5,X_B,AMU)

       UB_D=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=D_UB-UB_D
       ELSE
         CFCF=D_UB+UB_D
       ENDIF

      Else If (TYPE_V.EQ.'H+') then

       CONVP4=CONV(4,X_A,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       C_BB=CONVP4*CONVP5

       CONVM4=CONV(4,X_B,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       BB_C=CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=C_BB-BB_C
       ELSE
         CFCF=C_BB+BB_C
       ENDIF

      Else If (TYPE_V.EQ.'HB') then

       CONVP4=CONV(5,X_A,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       B_BB=CONVP4*CONVP5

       CONVM4=CONV(5,X_B,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       BB_B=CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=B_BB-BB_B
       ELSE
         CFCF=B_BB+BB_B
       ENDIF

CsB___The pp parton luminosity is the same for the following:
      ELSE IF(TYPE_V.EQ.'Z0'.OR. TYPE_V.EQ.'A0' .OR.
     > TYPE_V.EQ.'AA' .Or. TYPE_V.EQ.'ZZ' .Or.
     > TYPE_V.EQ.'HZ' .OR. TYPE_V.EQ.'GL' ) THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

c       U_UB=CONV(1,X_A,AMU)*CONV(-1,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(-4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(-2,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(-3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(-5,X_B,AMU)
c
c       UB_U=CONV(1,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(2,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

!ZL
      ELSE IF(Type_V.Eq.'WW_UUB' .or. Type_V.Eq.'WW_DDB') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       IF(Type_V.Eq.'WW_UUB') THEN
       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       d_db=0.d0
       db_d=0.d0

       ELSEIF(Type_V.Eq.'WW_DDB') THEN
       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

       u_ub=0.d0
       ub_u=0.d0

       ENDIF


       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

CCPY
CJI Sept 2013: Added ZU and ZD to match c++ code
      ELSE IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     >   .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
       enddo

        IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'ZU') THEN

       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)
        
       d_db=0.d0
       db_d=0.d0
       
        ELSEIF(Type_V.Eq.'Z0_RA_DDB' .or. Type_V.Eq.'ZD') THEN 

       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)
     
       u_ub=0.d0
       ub_u=0.d0

        ENDIF

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF


c------------------------------------
      ENDIF

      GOTO 888

997   CONTINUE

C THIS IS FOR IBEAM=0 OR IBEAM=-2
C FOR PROTON-NUCLEUS SCATTERING
C FOR (P P)

       do icnvx=1,5

C IHADRON=1 FOR PROTON AND -2 FOR PION_MINUS
         IF(IBEAM.EQ.-2)THEN
C xa is alwasy for PION, xb for proton or neutron
          IHADRON=-2
         ELSE
C for proton
          IHADRON=1
         ENDIF

         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
         endif

C This is for proton
         IHADRON=1
         cnvxb(icnvx)=conv(icnvx,x_b,amu)

         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

c       U_UB=CONV(1,X_A,AMU)*CONV(-1,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(-4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(-2,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(-3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(-5,X_B,AMU)
c
c       UB_U=CONV(1,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(2,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       CFCF_P=U_UB+D_DB+UB_U+DB_D

       IF(FRACT_N.LT.1.0D-4) GOTO 886

C FOR (P N)
CsB This process is set up for A0, AA, AG, Z0 and GG cases.
C   (We do not produce W's at fixed target and we do not calculate p N -> W at a
C   collider.)
CCPY MAY 1996
C THESE LINES PRESENT IN THE OLD VERSION WHICH GAVE RESULTS IN THE
C PAPER, PRD 50 (1994) 4415, WITH LADINSKY WERE WRONG.
C       u_ub=cnvxa(2)*cnvxb(-2)
C     >     +cnvxa(4)*cnvxb(-4)
C
C       d_db=cnvxa(1)*cnvxb(-1)
C     >     +cnvxa(3)*cnvxb(-3)
C     >     +cnvxa(5)*cnvxb(-5)
C
C       ub_u=cnvxb(2)*cnvxa(-2)
C     >     +cnvxb(4)*cnvxa(-4)
C
C       db_d=cnvxb(1)*cnvxa(-1)
C     >     +cnvxb(3)*cnvxa(-3)
C     >     +cnvxb(5)*cnvxa(-5)
C
CCPY THE FOLLOWING ARE THE CORRECTED CODE
C xa is for P, and xb is for N

       u_ub=cnvxa(1)*cnvxb(-2)
     >     +cnvxa(4)*cnvxb(-4)

       d_db=cnvxa(2)*cnvxb(-1)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       ub_u=cnvxb(2)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       db_d=cnvxb(1)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

CCPY THE FOLLOWING LINES WERE CORRECTED IN MAY 1996.
c       U_UB=CONV(1,X_A,AMU)*CONV(-2,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(-4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(-1,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(-3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(-5,X_B,AMU)
c
c       UB_U=CONV(2,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(1,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       CFCF_N=U_UB+D_DB+UB_U+DB_D

  886  CONTINUE

       CFCF=(1.0-FRACT_N)*CFCF_P+FRACT_N*CFCF_N

C       print*, "x_a, x_b, amu: ", x_a, x_b, amu
C       print*, "cfcf: ",  cfcf

888   CONTINUE
      RETURN

      End!GetCFCF

C --------------------------------------------------------------------------
      subroutine SetCf(do_asym)
C --------------------------------------------------------------------------
CPN Sets up the grids ConvGrdS and ConvGrdA with the convolutions
C   of the C-functions for the resummed piece.
      implicit NONE
      include 'common.for'
      real*8 qg(200),pt(200),y(200)
      integer IPTMIN,IPTMAX, IPTSTP,IYMIN,IYMAX, IYSTP,
     >     IQMIN,IQMAX, IQSTP
      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &     IQMIN, IQMAX, IQSTP

      REAL*8 amu,b,tem,tiny
      integer ib,j,nf_eff
      logical do_asym
      data tiny/1d-8/
      
      
      print *,'Creating the grids for the convolutions',
     >        ' of C-functions...'

CPN Setup some parameters
      S_Q=QG(iqin)

CCPY Feb 2016: To have the correct NF_EFF in SUDAKOV and CXFCXF
C      S_LAM=ALAMBD(S_Q)
C      NF_EFF=NFL(S_Q)
Cc this setup depends on NF_EFF and thus S_Q, and so must be done repetively
C      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR.
C     .     TYPE_V.EQ.'ZG' .OR. TYPE_V.EQ.'GG') THEN ! g g initial state
C        CALL HSETUP(NF_EFF)
C      ELSE                      ! q Q initial state
C        CALL WSETUP(NF_EFF)
C      ENDIF


cpn The smallest value of amu in the grid for the space-like processes is
cpn given by C3/ecm
      sml_b=C3/ecm +tiny

      DO J = IYMIN, IYMAX, IYSTP !loop  over Y
CsB____Definition of momentum fractions x1 and x2
        x_a=DEXP(Y(J))*QG(IQIN)/ECM
        x_b=DEXP(-Y(J))*QG(IQIN)/ECM

CPN                                           For AmuGrd(0), we use the value
C                                             of the cutoff at small b

        AmuGrd(Nconv)=C3/sml_b
        amu=AmuGrd(Nconv)

        if (x_a.gt.1.or.x_b.gt.1) then
          do ib = 0, NConv
            ConvGrdS(j,ib) =0.0 
            ConvGrdA(j,ib) =0.0
          enddo
        else
cpn                                                   Get CfCf at b=sml_b
          lepasy=0              !Symmetric piece
          Call GetCfCf(amu,tem)
          ConvGrdS(j,Nconv) = tem
          if (do_asym) then
            lepasy=1            !Antisymmetric piece
            Call GetCfCf(amu,tem)
            ConvGrdA(j,Nconv) = tem
          endif
          

          do ib=1,Nconv         !loop over b
            b=bmax*ib/Nconv
            AmuGrd(Nconv-ib)=C3/b
            amu=AmuGrd(Nconv-ib)
            
CCPY: IMHQMASS
            B_INT=B

            lepasy=0            !Symmetric piece
            Call GetCfCf(amu,tem) !Get CfCf at this b
            ConvGrdS(j,Nconv-ib) = tem
            if (do_asym) then
              lepasy=1          !Antisymmetric piece
              Call GetCfCf(amu,tem) !Get CfCf at this b
              ConvGrdA(j,Nconv-ib) = tem
            endif 

          enddo                 !loop over b
        endif    !x_a < 1 and x_b <1
      enddo                     !loop over Y

      print *,'... created successfully'
      RETURN
      End!SetCf


C --------------------------------------------------------------------------
      SUBROUTINE CXFCXF(B_SUD,CFCF)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 B_SUD,CFCF
      Integer HasJet
      REAL*8 D1s, R,t,ptj,u
      COMMON /Jet/ u,D1s, R,t,ptj, HasJet

      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR.
     .   TYPE_V.EQ.'ZG' .OR. TYPE_V.EQ.'GG') THEN
C FOR HIGGS
       CALL HCXFCXF(B_SUD,CFCF)
       ELSE IF(TYPE_V.EQ.'HJ') THEN
          CALL HJCXFCXF(B_SUD,CFCF)
      ELSE
C FOR W+, W-, Z0, A0 AND GL
       if (IFast.ge.1) then
        call WCXFCXF_PN(B_SUD,CFCF)
       else
        CALL WCXFCXF(B_SUD,CFCF)
       Endif

      ENDIF

      RETURN
      END!cxfcxf

C --------------------------------------------------------------------------
      SUBROUTINE WCXFCXF_PN(B_SUD,CFCF)
C --------------------------------------------------------------------------
cpn Starting from version 0.3, the purpose of this function is
C to interpolate the convolutions of the C-function with the  PDFs of
C the incoming particles. The convolutions are supplied through the array
C ConvGrd.

      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONV,X,AMU,tem,Dtem,b_sud,s_bstar,CfCf

CPN   New stuff
      integer nmu,Lmu,k
      parameter (nmu=4) !The (even) # of points used for the interpolation
      real*8 Ymu(nmu)
      INTEGER NF_EFF
      logical Testing

      REAL*8 Q0_PDF,q0pdf

      data Testing /.False./

      S_BSTAR=B_SUD/(DSQRT(1.0+(B_SUD/BMAX)**2))

CCPY: IMHQMASS
      B_INT=B_SUD

cpn   To prevent the accidents at small b
cpn      if (s_bstar.lt.sml_b) s_bstar=sml_b

CsB___Set the scale for C functions
      IF (IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     .    IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) THEN
C     CANONICAL CHOICE FOR W-BOSON
        AMU=C3/(S_BSTAR)
CCPY        Else If (IFLAG_C3.EQ.5) then
      Else 
        AMU=C3/S_BSTAR
      ENDIF
ccpy dec 1993: for very small b_sud
      if(amu.gt.ecm) then
        cfcf=0.d0
        return
      endif


CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
CBY April,2019
C       Q0_PDF=1.3
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF
CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
C      if(initlhapdf==1) then
C       Q0_PDF=q0pdf()
C      else
C       Q0_PDF=1.3D0
C      endif
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF

CPN                                              Find out the location of AMU
      call LOCATE(AmuGrd, NConv, AMU, LMU)
CPN *                                                      If not on the grid
C      If ((AmuGrd(NConv)-Amu)*(Amu-AmuGrd(0)).LT.0.e-9 ) then
C        Print*, ' Warning from Conv: Extrapolation used'
C        Print*, ' Amu,Lmu,NConv= ', Amu,Lmu,NConv
C      EndIf

CsB *                                                    If close to the ends
CsB                                      ( This works for even nmu for sure )
      If ( Lmu.LT.nmu/2)      Lmu = nmu/2
      If ( Lmu.GT.NConv-nmu/2)  Lmu = NConv-nmu/2+1
      If (Testing) Print *, ' LMu = ', LMu

CsB * Fill the dummy array(s) you want to interpolate
C     ( This works for even nmu for sure )


      if (lepasy.eq.0) then                 !interpolate the symmetric piece
       Do k = 1,nmu
         Ymu(k) = ConvGrdS(iyin,Lmu-nmu/2+k-1)
        End Do
      elseif (lepasy.eq.1) then         !interpolate the anti-symmetric piece
       Do k = 1,nmu
         Ymu(k) = ConvGrdA(iyin,Lmu-nmu/2+k-1)
        End Do
      endif !lepasy

CsB * Interpolate
      If (Testing) then
        Print*, ' AmuGrd(lmu-nmu/2),amu =',AmuGrd(lmu-nmu/2),amu
      EndIf

cpn      Print*, AmuGrd(lmu-nmu/2)
c      Print*, '******'
c      Print*, Ymu
c      Print*, '******'
c      Print*, 'nmu,amu,tem,Dtem'
c      Print*, nmu,amu,tem,Dtem

      Call PolInt( AmuGrd(lmu-nmu/2),Ymu,nmu,amu,tem,Dtem)

      CfCf=tem

      End!WCxfCxf_PN

C --------------------------------------------------------------------------
      SUBROUTINE WCXFCXF(B_SUD,CFCF)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 AMU,CONV,B_SUD,CFCF,S_BSTAR
      REAL*8 U_DB,DB_U,D_UB,UB_D,U_UB,D_DB,UB_U,DB_D
C_yfu
      REAL*8 S_SB,SB_S,C_CB,CB_C
      REAL*8 CONVP1,CONVP2,CONVP3,CONVP4,CONVP5,
     >CONVM1,CONVM2,CONVM3,CONVM4,CONVM5
      REAL*8 CFCF_P,CFCF_N

      integer icnvx
      real*8 cnvxa(-6:6),cnvxb(-6:6)

CCPY      INTEGER I_CXFF
      REAL*8 PDF1,PDF2,PDF3,PDF4,PDF5
      REAL*8 PDFM1,PDFM2,PDFM3,PDFM4,PDFM5

      Real*8 C_BB,BB_C,B_BB,BB_B
      REAL*8 Q0_PDF, q0pdf

      INTEGER NF_EFF
C yfu
      CHARACTER*40 envZUZDFLAG
      INTEGER iZUZDFLAG

CsB___Definition of b*
      S_BSTAR=B_SUD/(DSQRT(1.0+(B_SUD/BMAX)**2))

CCPY: IMHQMASS
      B_INT=B_SUD

CsB___Set the scale for C functions
      IF (IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     .    IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) THEN
C     CANONICAL CHOICE FOR W-BOSON
        AMU=C3/S_BSTAR
CCPY      Else If (IFLAG_C3.EQ.5 .or. iflag_c3.eq.6) then
      Else
        AMU=C3/S_BSTAR
      ENDIF
ccpy dec 1993: for very small b_sud
      if(amu.gt.ecm) then
        cfcf=0.d0
        return
      endif

CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
CBY April,2019
C       Q0_PDF=1.3
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF
CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
C      if(initlhapdf==1) then
C          call getq2min(1, q0_pdf)
C          q0_pdf = dsqrt(q0_pdf)
C      else
C       Q0_PDF=1.3D0
C      endif
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF

CJI September 2014: Add running of NF
      NF_EFF = NFL(amu)
CCPY Dec 2014:
      CALL WSETUP(NF_EFF)

      IF(IBEAM.NE.-1) GOTO 999
C THIS IS FOR IBEAM=-1
C FOR P PBAR

      IF(TYPE_V.EQ.'W+'.OR. TYPE_V.EQ.'HP'
     > .OR. Type_V.Eq.'W+_RA_UDB') THEN

CCPY IF I_CXFF=1 THEN DO
C       CXF(X_A)*F(X_B) + F(X_A)*CXF(X_B)
CsB This is to return only O(alpha_s) contributions NOT O(alpha_s^2) ones.
C   To check the total cross section, LCA, ect. against the NLO codes.
CCPY       I_CXFF=0
       IF(I_CXFF.NE.1) THEN
CCPY THIS IS FOR NORMAL RUNS: DO    CXF(X_A)*CXF(X_B)

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(2,X_B,AMU)
       CONVP3=CONV(3,X_B,AMU)
       CONVP5=CONV(5,X_B,AMU)

       U_DB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-1,X_B,AMU)
       CONVM4=CONV(-4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       DB_U=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       ELSEIF(I_CXFF.EQ.1) THEN
C FOR CHECKING:  DO      CXF(X_A)*F(X_B) + F(X_A)*CXF(X_B)

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(2,X_B,AMU)
       CONVP3=CONV(3,X_B,AMU)
       CONVP5=CONV(5,X_B,AMU)

       PDF1=APDF(1,X_A,AMU)
       PDF4=APDF(4,X_A,AMU)
       PDF2=APDF(2,X_B,AMU)
       PDF3=APDF(3,X_B,AMU)
       PDF5=APDF(5,X_B,AMU)

       U_DB=VKM(1,2)*(CONVP1*PDF2+PDF1*CONVP2-PDF1*PDF2)
     >     +VKM(1,3)*(CONVP1*PDF3+PDF1*CONVP3-PDF1*PDF3)
     >     +VKM(1,5)*(CONVP1*PDF5+PDF1*CONVP5-PDF1*PDF5)
     >     +VKM(4,2)*(CONVP4*PDF2+PDF4*CONVP2-PDF4*PDF2)
     >     +VKM(4,3)*(CONVP4*PDF3+PDF4*CONVP3-PDF4*PDF3)
     >     +VKM(4,5)*(CONVP4*PDF5+PDF4*CONVP5-PDF4*PDF5)

       CONVM1=CONV(-1,X_B,AMU)
       CONVM4=CONV(-4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       PDFM1=APDF(-1,X_B,AMU)
       PDFM4=APDF(-4,X_B,AMU)
       PDFM2=APDF(-2,X_A,AMU)
       PDFM3=APDF(-3,X_A,AMU)
       PDFM5=APDF(-5,X_A,AMU)

       DB_U=VKM(1,2)*(CONVM1*PDFM2+PDFM1*CONVM2-PDFM1*PDFM2)
     >     +VKM(1,3)*(CONVM1*PDFM3+PDFM1*CONVM3-PDFM1*PDFM3)
     >     +VKM(1,5)*(CONVM1*PDFM5+PDFM1*CONVM5-PDFM1*PDFM5)
     >     +VKM(4,2)*(CONVM4*PDFM2+PDFM4*CONVM2-PDFM4*PDFM2)
     >     +VKM(4,3)*(CONVM4*PDFM3+PDFM4*CONVM3-PDFM4*PDFM3)
     >     +VKM(4,5)*(CONVM4*PDFM5+PDFM4*CONVM5-PDFM4*PDFM5)

       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=U_DB-DB_U
       ELSE
         CFCF=U_DB+DB_U
       ENDIF

      ELSE IF(TYPE_V.EQ.'W-'.OR. TYPE_V.EQ.'HM'
     > .or. Type_V.Eq.'W-_RA_DUB') THEN

       CONVP1=CONV(1,X_B,AMU)
       CONVP4=CONV(4,X_B,AMU)
       CONVP2=CONV(2,X_A,AMU)
       CONVP3=CONV(3,X_A,AMU)
       CONVP5=CONV(5,X_A,AMU)

       D_UB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-1,X_A,AMU)
       CONVM4=CONV(-4,X_A,AMU)
       CONVM2=CONV(-2,X_B,AMU)
       CONVM3=CONV(-3,X_B,AMU)
       CONVM5=CONV(-5,X_B,AMU)

       UB_D=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=D_UB-UB_D
       ELSE
         CFCF=D_UB+UB_D
       ENDIF

      Else If (TYPE_V.EQ.'H+') then

CCPY       I_CXFF=0
       IF(I_CXFF.NE.1) THEN

         CONVP4=CONV(4,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         C_BB=CONVP4*CONVP5

         CONVM4=CONV(-4,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         BB_C=CONVM4*CONVM5

       ELSEIF(I_CXFF.EQ.1) THEN

         CONVP4=CONV(4,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         PDF4=APDF(4,X_A,AMU)
         PDF5=APDF(5,X_B,AMU)

         C_BB=CONVP4*PDF5+PDF4*CONVP5-PDF4*PDF5

         CONVM4=CONV(-4,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         PDFM4=APDF(-4,X_B,AMU)
         PDFM5=APDF(-5,X_A,AMU)

         BB_C=CONVM4*PDFM5+PDFM4*CONVM5-PDFM4*PDFM5

       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=C_BB-BB_C
       ELSE
         CFCF=C_BB+BB_C
       ENDIF

      Else If (TYPE_V.EQ.'HB') then

CCPY       I_CXFF=0
       IF(I_CXFF.NE.1) THEN

         CONVP4=CONV(5,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         B_BB=CONVP4*CONVP5

         CONVM4=CONV(-5,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         BB_B=CONVM4*CONVM5

       ELSEIF(I_CXFF.EQ.1) THEN

         CONVP4=CONV(5,X_A,AMU)
         CONVP5=CONV(5,X_B,AMU)

         PDF4=APDF(5,X_A,AMU)
         PDF5=APDF(5,X_B,AMU)

         B_BB=CONVP4*PDF5+PDF4*CONVP5-PDF4*PDF5

         CONVM4=CONV(-5,X_B,AMU)
         CONVM5=CONV(-5,X_A,AMU)

         PDFM4=APDF(-5,X_B,AMU)
         PDFM5=APDF(-5,X_A,AMU)

         BB_B=CONVM4*PDFM5+PDFM4*CONVM5-PDFM4*PDFM5

       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=B_BB-BB_B
       ELSE
         CFCF=B_BB+BB_B
       ENDIF

CsB___The pP parton luminosity is the same for the following:
      ELSE IF(TYPE_V.EQ.'Z0'.OR. TYPE_V.EQ.'A0' .OR.
     > TYPE_V.EQ.'AA' .Or. TYPE_V.EQ.'ZZ' .Or.
     > TYPE_V.EQ.'HZ' .OR. TYPE_V.EQ.'GL' ) THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       u_ub=cnvxa(1)*cnvxb(1)
     >     +cnvxa(4)*cnvxb(4)

       d_db=cnvxa(2)*cnvxb(2)
     >     +cnvxa(3)*cnvxb(3)
     >     +cnvxa(5)*cnvxb(5)

       ub_u=cnvxb(-1)*cnvxa(-1)
     >     +cnvxb(-4)*cnvxa(-4)

       db_d=cnvxb(-2)*cnvxa(-2)
     >     +cnvxb(-3)*cnvxa(-3)
     >     +cnvxb(-5)*cnvxa(-5)

c
c       U_UB=CONV(1,X_A,AMU)*CONV(1,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(2,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(5,X_B,AMU)
c
c       UB_U=CONV(-1,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(-4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(-2,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(-3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(-5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

!ZL
      ELSE IF(Type_V.Eq.'WW_UUB' .or. Type_V.Eq.'WW_DDB') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       IF(Type_V.Eq.'WW_UUB') THEN
       u_ub=cnvxa(1)*cnvxb(1)
     >     +cnvxa(4)*cnvxb(4)

       ub_u=cnvxb(-1)*cnvxa(-1)
     >     +cnvxb(-4)*cnvxa(-4)

       d_db=0.d0
       db_d=0.d0

       ELSEIF(Type_V.Eq.'WW_DDB') THEN

       d_db=cnvxa(2)*cnvxb(2)
     >     +cnvxa(3)*cnvxb(3)
     >     +cnvxa(5)*cnvxb(5)

       db_d=cnvxb(-2)*cnvxa(-2)
     >     +cnvxb(-3)*cnvxa(-3)
     >     +cnvxb(-5)*cnvxa(-5)

       u_ub=0.d0
       ub_u=0.d0

       ENDIF

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

CCPY
CJI Sept 2013: Added in ZU and ZD to match c++ code
      ELSE IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     >     .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
       enddo

        IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'ZU') THEN
        
       u_ub=cnvxa(1)*cnvxb(1)
     >     +cnvxa(4)*cnvxb(4)

       ub_u=cnvxb(-1)*cnvxa(-1)
     >     +cnvxb(-4)*cnvxa(-4)

       d_db=0.d0
       db_d=0.d0
       
        ELSEIF(Type_V.Eq.'Z0_RA_DDB' .or. Type_V.Eq.'ZD') THEN 
        
       d_db=cnvxa(2)*cnvxb(2)
     >     +cnvxa(3)*cnvxb(3)
     >     +cnvxa(5)*cnvxb(5)

       db_d=cnvxb(-2)*cnvxa(-2)
     >     +cnvxb(-3)*cnvxa(-3)
     >     +cnvxb(-5)*cnvxa(-5)

       u_ub=0.d0
       ub_u=0.d0

        ENDIF
        
       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF


C--------------------------
      ENDIF

      GOTO 888

999   CONTINUE
      IF(IBEAM.NE.1) GOTO 997
C THIS IS FOR IBEAM=1
C FOR P P

C_yfu
      call getenv('ZUZDFLAG',envZUZDFLAG)

      IF(TYPE_V.EQ.'W+'.OR. TYPE_V.EQ.'HP'
     >  .OR. Type_V.Eq.'W+_RA_UDB') THEN

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(-2,X_B,AMU)
       CONVP3=CONV(-3,X_B,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       CONVM1=CONV(1,X_B,AMU)
       CONVM4=CONV(4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

C_yfu seperate five quarks
       IF ( envZUZDFLAG(1:1).eq.' ') then
       U_DB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       DB_U=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5
       ELSE
         READ (envZUZDFLAG,*) iZUZDFLAG
         IF(iZUZDFLAG.EQ.12) THEN
       U_DB=VKM(1,2)*CONVP1*CONVP2
       DB_U=0.d0
         ELSE IF(iZUZDFLAG.EQ.13) THEN
       U_DB=VKM(1,3)*CONVP1*CONVP3
       DB_U=0.d0
         ELSE IF(iZUZDFLAG.EQ.15) THEN
       U_DB=VKM(1,5)*CONVP1*CONVP5
       DB_U=0.d0
         ELSE IF(iZUZDFLAG.EQ.42) THEN
       U_DB=VKM(4,2)*CONVP4*CONVP2
       DB_U=0.d0
         ELSE IF(iZUZDFLAG.EQ.43) THEN
       U_DB=VKM(4,3)*CONVP4*CONVP3
       DB_U=0.d0
         ELSE IF(iZUZDFLAG.EQ.45) THEN
       U_DB=VKM(4,5)*CONVP4*CONVP5
       DB_U=0.d0
         ELSE IF(iZUZDFLAG.EQ.-12) THEN
       U_DB=0.d0
       DB_U=VKM(1,2)*CONVM1*CONVM2
         ELSE IF(iZUZDFLAG.EQ.-13) THEN
       U_DB=0.d0
       DB_U=VKM(1,3)*CONVM1*CONVM3 
         ELSE IF(iZUZDFLAG.EQ.-15) THEN
       U_DB=0.d0
       DB_U=VKM(1,5)*CONVM1*CONVM5 
         ELSE IF(iZUZDFLAG.EQ.-42) THEN
       U_DB=0.d0
       DB_U=VKM(4,2)*CONVM4*CONVM2 
         ELSE IF(iZUZDFLAG.EQ.-43) THEN
       U_DB=0.d0
       DB_U=VKM(4,3)*CONVM4*CONVM3 
         ELSE IF(iZUZDFLAG.EQ.-45) THEN
       U_DB=0.d0
       DB_U=VKM(4,5)*CONVM4*CONVM5 
         ENDIF
       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=U_DB-DB_U
       ELSE
         CFCF=U_DB+DB_U
       ENDIF

      ELSE IF(TYPE_V.EQ.'W-'.OR. TYPE_V.EQ.'HM'
     > .or. Type_V.Eq.'W-_RA_DUB') THEN

       CONVP1=CONV(-1,X_B,AMU)
       CONVP4=CONV(-4,X_B,AMU)
       CONVP2=CONV(2,X_A,AMU)
       CONVP3=CONV(3,X_A,AMU)
       CONVP5=CONV(5,X_A,AMU)

       CONVM1=CONV(-1,X_A,AMU)
       CONVM4=CONV(-4,X_A,AMU)
       CONVM2=CONV(2,X_B,AMU)
       CONVM3=CONV(3,X_B,AMU)
       CONVM5=CONV(5,X_B,AMU)

C_yfu seperate five quarks
       IF ( envZUZDFLAG(1:1).eq.' ') then
       D_UB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       UB_D=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5
       ELSE
         READ (envZUZDFLAG,*) iZUZDFLAG
         IF(iZUZDFLAG.EQ.12) THEN
       D_UB=VKM(1,2)*CONVP1*CONVP2
       UB_D=0.d0
         ELSE IF(iZUZDFLAG.EQ.13) THEN
       D_UB=VKM(1,3)*CONVP1*CONVP3
       UB_D=0.d0
         ELSE IF(iZUZDFLAG.EQ.15) THEN
       D_UB=VKM(1,5)*CONVP1*CONVP5
       UB_D=0.d0
         ELSE IF(iZUZDFLAG.EQ.42) THEN
       D_UB=VKM(4,2)*CONVP4*CONVP2
       UB_D=0.d0
         ELSE IF(iZUZDFLAG.EQ.43) THEN
       D_UB=VKM(4,3)*CONVP4*CONVP3
       UB_D=0.d0
         ELSE IF(iZUZDFLAG.EQ.45) THEN
       D_UB=VKM(4,5)*CONVP4*CONVP5
       UB_D=0.d0
         ELSE IF(iZUZDFLAG.EQ.-12) THEN
       D_UB=0.d0
       UB_D=VKM(1,2)*CONVM1*CONVM2
         ELSE IF(iZUZDFLAG.EQ.-13) THEN
       D_UB=0.d0
       UB_D=VKM(1,3)*CONVM1*CONVM3
         ELSE IF(iZUZDFLAG.EQ.-15) THEN
       D_UB=0.d0
       UB_D=VKM(1,5)*CONVM1*CONVM5
         ELSE IF(iZUZDFLAG.EQ.-42) THEN
       D_UB=0.d0
       UB_D=VKM(4,2)*CONVM4*CONVM2
         ELSE IF(iZUZDFLAG.EQ.-43) THEN
       D_UB=0.d0
       UB_D=VKM(4,3)*CONVM4*CONVM3
         ELSE IF(iZUZDFLAG.EQ.-45) THEN
       D_UB=0.d0
       UB_D=VKM(4,5)*CONVM4*CONVM5
         ENDIF
       ENDIF

       IF(LEPASY.EQ.1) THEN
         CFCF=D_UB-UB_D
       ELSE
         CFCF=D_UB+UB_D
       ENDIF

      Else If (TYPE_V.EQ.'H+') then

       CONVP4=CONV(4,X_A,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       C_BB=CONVP4*CONVP5

       CONVM4=CONV(4,X_B,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       BB_C=CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=C_BB-BB_C
       ELSE
         CFCF=C_BB+BB_C
       ENDIF

      Else If (TYPE_V.EQ.'HB') then

       CONVP4=CONV(5,X_A,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       B_BB=CONVP4*CONVP5

       CONVM4=CONV(5,X_B,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       BB_B=CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF=B_BB-BB_B
       ELSE
         CFCF=B_BB+BB_B
       ENDIF

CsB___The pp parton luminosity is the same for the following:
      ELSE IF(TYPE_V.EQ.'Z0'.OR. TYPE_V.EQ.'A0' .OR.
     > TYPE_V.EQ.'AA' .Or. TYPE_V.EQ.'ZZ' .Or.
     > TYPE_V.EQ.'HZ' .OR. TYPE_V.EQ.'GL' ) THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

c       U_UB=CONV(1,X_A,AMU)*CONV(-1,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(-4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(-2,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(-3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(-5,X_B,AMU)
c
c       UB_U=CONV(1,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(2,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

!ZL
      ELSE IF(Type_V.Eq.'WW_UUB' .or. Type_V.Eq.'WW_DDB') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

      IF(Type_V.Eq.'WW_UUB') THEN
       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

      d_db=0.d0
      db_d=0.d0

      ELSEIF(Type_V.Eq.'WW_DDB') THEN

       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

      u_ub=0.d0
      ub_u=0.d0

      ENDIF

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
         CFCF=U_UB+D_DB-UB_U-DB_D
       ELSE
         CFCF=U_UB+D_DB+UB_U+DB_D
       ENDIF

CCPY
CJI Sept 2013: Added in ZU and ZD to match c++ code
      ELSE IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     >   .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD') THEN

       do icnvx=1,5
         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         cnvxb(icnvx)=conv(icnvx,x_b,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
       enddo

        IF(Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'ZU')THEN
        
C       u_ub=cnvxa(1)*cnvxb(-1)
C     >     +cnvxa(4)*cnvxb(-4)

C       ub_u=cnvxb(1)*cnvxa(-1)
C     >     +cnvxb(4)*cnvxa(-4)

C_yfu seperate five quarks
       IF ( envZUZDFLAG(1:1).eq.' ') then
         u_ub=cnvxa(1)*cnvxb(-1)
         c_cb=cnvxa(4)*cnvxb(-4)
         ub_u=cnvxb(1)*cnvxa(-1)
         cb_c=cnvxb(4)*cnvxa(-4)
       ELSE
         READ (envZUZDFLAG,*) iZUZDFLAG
         IF(iZUZDFLAG.EQ.1) THEN
         u_ub=cnvxa(1)*cnvxb(-1)
         c_cb=0.d0
         ub_u=0.d0
         cb_c=0.d0
         ELSE IF(iZUZDFLAG.EQ.4) THEN
         u_ub=0.d0
         c_cb=cnvxa(4)*cnvxb(-4)
         ub_u=0.d0
         cb_c=0.d0
         ELSE IF(iZUZDFLAG.EQ.-1) THEN
         u_ub=0.d0
         c_cb=0.d0
         ub_u=cnvxb(1)*cnvxa(-1)
         cb_c=0.d0
         ELSE IF(iZUZDFLAG.EQ.-4) THEN
         u_ub=0.d0
         c_cb=0.d0
         ub_u=0.d0
         cb_c=cnvxb(4)*cnvxa(-4)
         ENDIF
       ENDIF

       d_db=0.d0
       db_d=0.d0
       
C_yfu
       s_sb=0.d0
       sb_s=0.d0
       b_bb=0.d0
       bb_b=0.d0

        ELSEIF(Type_V.Eq.'Z0_RA_DDB' .or. Type_V.Eq.'ZD')THEN 
        
C       d_db=cnvxa(2)*cnvxb(-2)
C     >     +cnvxa(3)*cnvxb(-3)
C     >     +cnvxa(5)*cnvxb(-5)

C       db_d=cnvxb(2)*cnvxa(-2)
C     >     +cnvxb(3)*cnvxa(-3)
C     >     +cnvxb(5)*cnvxa(-5)
      
C_yfu seperate five quarks
       IF ( envZUZDFLAG(1:1).eq.' ') then
         d_db=cnvxa(2)*cnvxb(-2)
         s_sb=cnvxa(3)*cnvxb(-3)
         b_bb=cnvxa(5)*cnvxb(-5)
         db_d=cnvxb(2)*cnvxa(-2)
         sb_s=cnvxb(3)*cnvxa(-3)
         bb_b=cnvxb(5)*cnvxa(-5)
       ELSE
         READ (envZUZDFLAG,*) iZUZDFLAG
         IF(iZUZDFLAG.EQ.2) THEN
         d_db=cnvxa(2)*cnvxb(-2)
         s_sb=0.d0
         b_bb=0.d0
         db_d=0.d0
         sb_s=0.d0
         bb_b=0.d0
         ELSE IF(iZUZDFLAG.EQ.3) THEN
         d_db=0.d0
         s_sb=cnvxa(3)*cnvxb(-3)
         b_bb=0.d0
         db_d=0.d0
         sb_s=0.d0
         bb_b=0.d0
         ELSE IF(iZUZDFLAG.EQ.5) THEN
         d_db=0.d0
         s_sb=0.d0
         b_bb=cnvxa(5)*cnvxb(-5)
         db_d=0.d0
         sb_s=0.d0
         bb_b=0.d0
         ELSE IF(iZUZDFLAG.EQ.-2) THEN
         d_db=0.d0
         s_sb=0.d0
         b_bb=0.d0
         db_d=cnvxb(2)*cnvxa(-2)
         sb_s=0.d0
         bb_b=0.d0
         ELSE IF(iZUZDFLAG.EQ.-3) THEN
         d_db=0.d0
         s_sb=0.d0
         b_bb=0.d0
         db_d=0.d0
         sb_s=cnvxb(3)*cnvxa(-3)
         bb_b=0.d0
         ELSE IF(iZUZDFLAG.EQ.-5) THEN
         d_db=0.d0
         s_sb=0.d0
         b_bb=0.d0
         db_d=0.d0
         sb_s=0.d0
         bb_b=cnvxb(5)*cnvxa(-5)
         ENDIF
       ENDIF

       u_ub=0.d0
       ub_u=0.d0

C_yfu
       c_cb=0.d0
       cb_c=0.d0

        ENDIF
        
       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

C_yfu seperate five quarks
       C_CB=C_CB*ZFF2(1)
       S_SB=S_SB*ZFF2(2)
       B_BB=B_BB*ZFF2(2)
       CB_C=CB_C*ZFF2(1)
       SB_S=SB_S*ZFF2(2)
       BB_B=BB_B*ZFF2(2)

       IF(LEPASY.EQ.1) THEN
C yfu
C         CFCF=U_UB+D_DB-UB_U-DB_D
           CFCF=U_UB+D_DB+C_CB+S_SB+B_BB
     >         -UB_U-DB_D-CB_C-SB_S-BB_B
       ELSE
C yfu add a flag for separate u_ub and ub_u, d_db and db_d
C         CFCF=U_UB+D_DB+UB_U+DB_D
           CFCF=U_UB+D_DB+C_CB+S_SB+B_BB
     >         +UB_U+DB_D+CB_C+SB_S+BB_B
       ENDIF


c----------------------------------------------
      ENDIF

      GOTO 888

997   CONTINUE

C THIS IS FOR IBEAM=0 OR IBEAM=-2
C FOR PROTON-NUCLEUS SCATTERING
C FOR (P P)
      IF(TYPE_V.EQ.'W+'.OR. TYPE_V.EQ.'HP'
     > .OR. Type_V.Eq.'W+_RA_UDB') THEN

       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(-2,X_B,AMU)
       CONVP3=CONV(-3,X_B,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       U_DB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(1,X_B,AMU)
       CONVM4=CONV(4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       DB_U=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF_P=U_DB-DB_U
       ELSE
         CFCF_P=U_DB+DB_U
       ENDIF
       IF(FRACT_N.LT.1.0D-4) GOTO 886
       CONVP1=CONV(1,X_A,AMU)
       CONVP4=CONV(4,X_A,AMU)
       CONVP2=CONV(-1,X_B,AMU)
       CONVP3=CONV(-3,X_B,AMU)
       CONVP5=CONV(-5,X_B,AMU)

       U_DB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(2,X_B,AMU)
       CONVM4=CONV(4,X_B,AMU)
       CONVM2=CONV(-2,X_A,AMU)
       CONVM3=CONV(-3,X_A,AMU)
       CONVM5=CONV(-5,X_A,AMU)

       DB_U=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF_N=U_DB-DB_U
       ELSE
         CFCF_N=U_DB+DB_U
       ENDIF

      ELSE IF(TYPE_V.EQ.'W-'.OR. TYPE_V.EQ.'HM'
     > .or. Type_V.Eq.'W-_RA_DUB') THEN

       CONVP1=CONV(-1,X_B,AMU)
       CONVP4=CONV(-4,X_B,AMU)
       CONVP2=CONV(2,X_A,AMU)
       CONVP3=CONV(3,X_A,AMU)
       CONVP5=CONV(5,X_A,AMU)

       D_UB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-1,X_A,AMU)
       CONVM4=CONV(-4,X_A,AMU)
       CONVM2=CONV(2,X_B,AMU)
       CONVM3=CONV(3,X_B,AMU)
       CONVM5=CONV(5,X_B,AMU)

       UB_D=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF_P=D_UB-UB_D
       ELSE
         CFCF_P=D_UB+UB_D
       ENDIF
       IF(FRACT_N.LT.1.0D-4) GOTO 886
       CONVP1=CONV(-2,X_B,AMU)
       CONVP4=CONV(-4,X_B,AMU)
       CONVP2=CONV(2,X_A,AMU)
       CONVP3=CONV(3,X_A,AMU)
       CONVP5=CONV(5,X_A,AMU)

       D_UB=VKM(1,2)*CONVP1*CONVP2
     >     +VKM(1,3)*CONVP1*CONVP3
     >     +VKM(1,5)*CONVP1*CONVP5
     >     +VKM(4,2)*CONVP4*CONVP2
     >     +VKM(4,3)*CONVP4*CONVP3
     >     +VKM(4,5)*CONVP4*CONVP5

       CONVM1=CONV(-2,X_A,AMU)
       CONVM4=CONV(-4,X_A,AMU)
       CONVM2=CONV(2,X_B,AMU)
       CONVM3=CONV(3,X_B,AMU)
       CONVM5=CONV(5,X_B,AMU)

       UB_D=VKM(1,2)*CONVM1*CONVM2
     >     +VKM(1,3)*CONVM1*CONVM3
     >     +VKM(1,5)*CONVM1*CONVM5
     >     +VKM(4,2)*CONVM4*CONVM2
     >     +VKM(4,3)*CONVM4*CONVM3
     >     +VKM(4,5)*CONVM4*CONVM5

       IF(LEPASY.EQ.1) THEN
         CFCF_N=D_UB-UB_D
       ELSE
         CFCF_N=D_UB+UB_D
       ENDIF

       else 
       do icnvx=1,5

C IHADRON=1 FOR PROTON AND -2 FOR PION_MINUS
         IF(IBEAM.EQ.-2)THEN
C xa is alwasy for PION, xb for proton or neutron
          IHADRON=-2
         ELSE
C for proton
          IHADRON=1
         ENDIF

         cnvxa(icnvx)=conv(icnvx,x_a,amu)
         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxa(-icnvx)=conv(-icnvx,x_a,amu)
         else
           cnvxa(-icnvx)=cnvxa(icnvx)
         endif

C This is for proton
         IHADRON=1
         cnvxb(icnvx)=conv(icnvx,x_b,amu)

         if(icnvx.eq.1.or.icnvx.eq.2) then
           cnvxb(-icnvx)=conv(-icnvx,x_b,amu)
         else
           cnvxb(-icnvx)=cnvxb(icnvx)
         endif
      enddo

       u_ub=cnvxa(1)*cnvxb(-1)
     >     +cnvxa(4)*cnvxb(-4)

       d_db=cnvxa(2)*cnvxb(-2)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       ub_u=cnvxb(1)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       db_d=cnvxb(2)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

c       U_UB=CONV(1,X_A,AMU)*CONV(-1,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(-4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(-2,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(-3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(-5,X_B,AMU)
c
c       UB_U=CONV(1,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(2,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       CFCF_P=U_UB+D_DB+UB_U+DB_D

       IF(FRACT_N.LT.1.0D-4) GOTO 886

C FOR (P N)
CsB This process is set up for A0, AA, AG, Z0 and GG cases.
C   (We do not produce W's at fixed target and we do not calculate p N -> W at a
C   collider.)
CCPY MAY 1996
C THESE LINES PRESENT IN THE OLD VERSION WHICH GAVE RESULTS IN THE
C PAPER, PRD 50 (1994) 4415, WITH LADINSKY WERE WRONG.
C       u_ub=cnvxa(2)*cnvxb(-2)
C     >     +cnvxa(4)*cnvxb(-4)
C
C       d_db=cnvxa(1)*cnvxb(-1)
C     >     +cnvxa(3)*cnvxb(-3)
C     >     +cnvxa(5)*cnvxb(-5)
C
C       ub_u=cnvxb(2)*cnvxa(-2)
C     >     +cnvxb(4)*cnvxa(-4)
C
C       db_d=cnvxb(1)*cnvxa(-1)
C     >     +cnvxb(3)*cnvxa(-3)
C     >     +cnvxb(5)*cnvxa(-5)
C
CCPY THE FOLLOWING ARE THE CORRECTED CODE
C xa is for P, and xb is for N

       u_ub=cnvxa(1)*cnvxb(-2)
     >     +cnvxa(4)*cnvxb(-4)

       d_db=cnvxa(2)*cnvxb(-1)
     >     +cnvxa(3)*cnvxb(-3)
     >     +cnvxa(5)*cnvxb(-5)

       ub_u=cnvxb(2)*cnvxa(-1)
     >     +cnvxb(4)*cnvxa(-4)

       db_d=cnvxb(1)*cnvxa(-2)
     >     +cnvxb(3)*cnvxa(-3)
     >     +cnvxb(5)*cnvxa(-5)

CCPY THE FOLLOWING LINES WERE CORRECTED IN MAY 1996.
c       U_UB=CONV(1,X_A,AMU)*CONV(-2,X_B,AMU)
c     >     +CONV(4,X_A,AMU)*CONV(-4,X_B,AMU)
c
c       D_DB=CONV(2,X_A,AMU)*CONV(-1,X_B,AMU)
c     >     +CONV(3,X_A,AMU)*CONV(-3,X_B,AMU)
c     >     +CONV(5,X_A,AMU)*CONV(-5,X_B,AMU)
c
c       UB_U=CONV(2,X_B,AMU)*CONV(-1,X_A,AMU)
c     >     +CONV(4,X_B,AMU)*CONV(-4,X_A,AMU)
c
c       DB_D=CONV(1,X_B,AMU)*CONV(-2,X_A,AMU)
c     >     +CONV(3,X_B,AMU)*CONV(-3,X_A,AMU)
c     >     +CONV(5,X_B,AMU)*CONV(-5,X_A,AMU)

       U_UB=U_UB*ZFF2(1)
       D_DB=D_DB*ZFF2(2)
       UB_U=UB_U*ZFF2(1)
       DB_D=DB_D*ZFF2(2)

       CFCF_N=U_UB+D_DB+UB_U+DB_D
       endif

  886  CONTINUE

       CFCF=(1.0-FRACT_N)*CFCF_P+FRACT_N*CFCF_N


888   CONTINUE
      RETURN
      END!Wcxfcxf


C --------------------------------------------------------------------------
      FUNCTION CONV(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONV,X,AMU
      INTEGER IPDF, i
      REAL*8 CONVQ,CONVG, DCONVQQ,CONVQQ, DCONVQQP, CONVQQP, CONVGG
      EXTERNAL CONVQ,CONVG, CONVQQ, CONVQQP, CONVGG

CsB   I_Proc = 1 Calculate only QI QJ -->  G V process
C     I_Proc = 2 Calculate only  G QI --> QJ V process

      CONV = 0
      DCONVQQP = 0

      If (i_Proc.Eq.0) then
CCPY THIS IS CORRECTED ON SEPT 1996.
C IT HAS NO EFFECT ON PREVIOUS NUMERICAL RESULTS
C BECAUSE OF THE STRUCTURE OF FUNCTION CJXFG.
C       CONV = CONVQ(IPDF,X,AMU) + CONVG(IPDF,X,AMU)
        CONV = CONVQ(IPDF,X,AMU) + CONVG(0,X,AMU)
     >         + CONVQQ(-ipdf,x,amu)
        do i=-5,5
            if(i.ne.0 .and. i.ne.ipdf .and. i.ne.-ipdf) then
                DCONVQQP = CONVQQP(i,x,amu)
                CONV = CONV + DCONVQQP
            endif
        enddo

      Else If (i_Proc.Eq.1) then
        CONV = CONVQ(IPDF,X,AMU)
      Else If (i_Proc.Eq.2) then
CsB        CONV = CONVG(IPDF,X,AMU)
        CONV = CONVG(0,X,AMU)
CJI Jan 2015: Add in qqbar <- qq(bar)'
      Else If (i_Proc.Eq.3) then
        DO i=-5,5 
            if(i.ne.0 .and. i.ne.ipdf.and.i.ne.-ipdf) then
                DCONVQQP = CONVQQP(i,x,amu)
                CONV = CONV + DCONVQQP
            endif
        enddo
CJI Jan 2015: Add in qqbar <- qq (qbarqbar)
      Else If (i_Proc.Eq.4) then
            CONV = CONVQQ(-ipdf,x,amu)
!        enddo
      Else
        Print*, ' Not implemented ! '
        Stop
      End if
        
      RETURN
      END



C --------------------------------------------------------------------------
      FUNCTION CONVQ(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONVQ,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT
      REAL*8 CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ_3,CJJ1XFJ,
     >S_ALPI, CJJ2XFJ,CJJ2XFJ_1,CJJ2XFJ_2,CJJ2XFJ_3,CJJ2XFJ_4,
     >CJJ2XFJ_5,CJJ2XFJ_6
      REAL*8 CJXFJ,SPCJFJ,CJXFJ2,SPCJFJ2,CJXFJ3
      Real*8 Omega
      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC
      EXTERNAL CJXFJ,SPCJFJ,CJXFJ2,SPCJFJ2,CJXFJ3

      REAL*8 H1_C2

      real*8 a1c,a2c,a3c,b1c,b2c,b3c
      common /connonical/ a1c,a2c,a3c,b1c,b2c,b3c

      real*8 splithoppet,c1p1xf,c1xf
      external splithoppet,c1p1xf,c1xf

CJI September 2016, add in the choice of scheme
      CHARACTER*3 SCHEME
      COMMON / RESUMSCHEME / SCHEME

      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU
      NF=NFL(AMU)

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
      ELSE
        FPDF=APDF(IPDF,X,AMU)
      ENDIF

C      print*, fpdf

      IF(GOAROUND) THEN
       CJJ0XFJ=FPDF
       CONVQ=CJJ0XFJ
       RETURN
      ENDIF

      CJJ0XFJ=0.
      CJJ1XFJ_1=0.
      CJJ1XFJ_2=0.
      CJJ1XFJ_3=0.
      CJJ1XFJ=0.
      CJJ2XFJ=0.
      CJJ2XFJ_1=0.
      CJJ2XFJ_2=0.
      CJJ2XFJ_3=0.
      CJJ2XFJ_4=0.
      CJJ2XFJ_5=0.
      CJJ2XFJ_6=0.
      S_ALPI=ALPI(AMU)

C      print*,"In convq: amu, s_alpi", amu, s_alpi

C OCT. 17, 1995
C      IF(NORDER.EQ.1) THEN
C       CJJ0XFJ=FPDF
C       CONVQ=CJJ0XFJ
C
C      ELSE IF(NORDER.EQ.2) THEN
CCPY
C      IF(NORDER.EQ.1 .OR. NORDER.EQ.2) THEN

      IF(N_WIL_C.EQ.0) THEN
       CJJ0XFJ=FPDF
       CONVQ=CJJ0XFJ
       RETURN

      ELSEIF(N_WIL_C.EQ.1) THEN

C DO \delta(1-Z)
        CJJ0XFJ=FPDF

CJI September 2016, add in the choice of scheme
      if(SCHEME.eq.'CSS') then
CsB iCF1 separates the alpha_s delta(1-Z) and LO + alpha_s (1-Z)+gluon
C   contributions.
C   Settings: iCF1 = 1 alpha_s delta(1-Z) qqB contribution returned
C             iCF1 = 2 LO + alpha_s (1-Z) qqB + Gq contributions
C DO \delta(1-Z)*CF*(PI^2-8)/4
      If (iCF1.NE.2) then 
CCPY THIS IS FOR USING RUNNING MASS, OR RUNNING YUKAWA COUPLING.
C IN THIS CASE, OMEGA=0
        If (Type_V.Eq.'H+') then ! C(1) for H+ production
CCPY Error: March 1, 2005         CJJ1XFJ_1 = 0.25d0*CF*(PI**2)*FPDF
          CJJ1XFJ_1 = 0.25d0*CF*(PI**2-2.0)*FPDF
        Else If (Type_V.Eq.'HB') then ! C(1) for b B -> A0
CCPY Error: March 1, 2005         CJJ1XFJ_1 = 0.25d0*CF*(PI**2)*FPDF
          CJJ1XFJ_1 = 0.25d0*CF*(PI**2-2.0)*FPDF
        else If (Type_V.Eq.'AA'.OR.Type_V.Eq.'ZZ'
     &    .or.Type_V.Eq.'WW_UUB'.or.Type_V.Eq.'WW_DDB') then 
CCPY PUT IN AN "EFFECTIVE" C(1) TERM FOR 'AA' PRODUCTION
          H1_C2=0.6D0
          CJJ1XFJ_1 = H1_C2*FPDF 
	else
C for the other qqB initiated process
          CJJ1XFJ_1 = 0.25d0*CF*(PI**2-8.d0)*FPDF

        ENDIF

        CJJ1XFJ_1 = CJJ1XFJ_1*S_ALPI

      ENDIF    ! iCF1.NE.2
      endif
      
C DO CF*(1-Z)/2
        UPLIM=1.0
        DNLIM=X
        AERR=ACC_AERR
        RERR=ACC_RERR
        IACTA=2
        IACTB=2

          CJJ1XFJ_2=ADZ3NT(CJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >    IER,IACTA,IACTB)
    

        IF(IER.NE.0) THEN
          WRITE(NWRT,*)' ERROR IN CJXFJ_2', CJJ1XFJ_2
CZL
CCPY	  CJJ1XFJ_2=0
!          CALL QUIT
        ENDIF

  100   Continue

CCPY April 16, 1996
C This is only correct for changing C3, but not for changing C2.
CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
       IF(IFLAG_C3.NE.1) THEN
C       IF(IFLAG_C3.GT.4) THEN
C ADD EXTRA TERMS
CsB Add CF*(9/16 - (ln(C1/(2 C2))-3/4+gamma_E)^2) (See Cfns.nb)
        CJJ1XFJ_1=CJJ1XFJ_1+CF*(9.D0/16.D0-
     >           (DLOG(C1/2.D0/C2)-3.D0/4.D0+EULER)**2)*FPDF*s_alpi
C DO -CF*(LN(C1/2)+EULER)*((1+Z^2)/(1-Z))_{+}
        CJJ1XFJ_3=-(DLOG(C3/2.D0)+EULER)*SPCJFJ(IPDF,X,AMU)*s_alpi
       ENDIF

       CJJ1XFJ=CJJ1XFJ_1+CJJ1XFJ_2+CJJ1XFJ_3


       If (iCF1.EQ.1) then
C Used only for double check
         ConvQ = Sqrt(CJJ0XFJ*S_ALPI*2.d0*(CJJ1XFJ_1+CJJ1XFJ_3))
       Else
         CONVQ=CJJ0XFJ+CJJ1XFJ!*S_ALPI
       End If

C       IF(DEBUG) THEN
C        Print*,'X,AMU,S_ALPI'
C        Print*,X,AMU,S_ALPI
C        Print*,'CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ,CONVQ,IPDF'
C        Print*,CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ,CONVQ,IPDF
C        Print*,'ICF1 =',ICF1
C       ENDIF


CJI: 2014 add C2 coefficient and improve organization to optimize the
C       calculation time
      ELSE IF(N_WIL_C.EQ.2) THEN
       
       S_ALPI=ALPI(AMU)

C DO \delta(1-Z)
        CJJ0XFJ=FPDF

      if(SCHEME.eq."CSS") then
CsB iCF1 separates the alpha_s delta(1-Z) and LO + alpha_s (1-Z)+gluon
C   contributions.
C   Settings: iCF1 = 1 alpha_s delta(1-Z) qqB contribution returned
C             iCF1 = 2 LO + alpha_s (1-Z) qqB + Gq contributions
C DO \delta(1-Z)*CF*(PI^2-8)/4
      If (iCF1.NE.2) then 
CCPY THIS IS FOR USING RUNNING MASS, OR RUNNING YUKAWA COUPLING.
C IN THIS CASE, OMEGA=0
        If (Type_V.Eq.'H+') then ! C(1) for H+ production
CCPY Error: March 1, 2005         CJJ1XFJ_1 = 0.25d0*CF*(PI**2)*FPDF
          CJJ1XFJ_1 = 0.25d0*CF*(PI**2-2.0)*FPDF*S_ALPI
        Else If (Type_V.Eq.'HB') then ! C(1) for b B -> A0
CCPY Error: March 1, 2005         CJJ1XFJ_1 = 0.25d0*CF*(PI**2)*FPDF
          CJJ1XFJ_1 = 0.25d0*CF*(PI**2-2.0)*FPDF*S_ALPI
        else If (Type_V.Eq.'AA'.OR.Type_V.Eq.'ZZ'
     &    .or.Type_V.Eq.'WW_UUB'.or.Type_V.Eq.'WW_DDB') then 
CCPY PUT IN AN "EFFECTIVE" C(1) TERM FOR 'AA' PRODUCTION
          H1_C2=0.6D0
          CJJ1XFJ_1 = H1_C2*FPDF*S_ALPI
	else
C for the other qqB initiated process
CJI January 2014: Added complete CJJ2XFJ term (Eur. Phys. J. C (2012) 72:2195)
C delta(1-z) piece
C C(1) Contribution
          CJJ1XFJ_1 = (0.25d0*CF*(PI**2-8.d0)*S_ALPI
C C(2) Contribution
     >      +(CA*CF*(59.D0/18.D0*ZETA3-1535.D0/192.D0
     >          +215.D0/216.D0*PI2-PI2**2/240.D0)
     >          +0.25*CF**2*(-15.D0*ZETA3+511.D0/16.D0-67.D0*PI2/12.D0
     >          +17.D0/45.D0*PI2**2)
     >          +1.D0/864.D0*CF*NF*(192.D0*ZETA3+1143.D0-152.D0*PI2)
     >          -1.D0/4.D0*((PI2/2.0-4.0)*CF)**2)
     >          *0.5*S_ALPI**2)*FPDF

        ENDIF

      ENDIF    ! iCF1.NE.2
      endif
      
C Z Dependent piece minus the plus-functions
        UPLIM=1.0
        DNLIM=X
CJI Lower accuracy
        AERR=ACC_AERR
        RERR=ACC_RERR
C       AERR=1E-8
C       RERR=1E-1
        IACTA=2
        IACTB=2

          CJJ1XFJ_2=ADZ3NT(CJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >    IER,IACTA,IACTB)
    

        IF(IER.NE.0) THEN
          WRITE(NWRT,*)' ERROR IN CJXFJ_2', CJJ1XFJ_2
CZL
C	  CJJ1XFJ_2=0
!          CALL QUIT
        ENDIF

 2100   Continue

CCPY April 16, 1996
C This is only correct for changing C3, but not for changing C2.
CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
       IF(IFLAG_C3.NE.1) THEN
C       IF(IFLAG_C3.GT.4) THEN
C ADD EXTRA TERMS
CsB Add CF*(9/16 - (ln(C1/(2 C2))-3/4+gamma_E)^2) (See Cfns.nb)
        CJJ1XFJ_1=CJJ1XFJ_1+CF*(9.D0/16.D0-
     >  (DLOG(C1/2.D0/C2)-3.D0/4.D0+EULER)**2)*FPDF*S_ALPI
C DO -CF*(LN(C1/2)+EULER)*((1+Z^2)/(1-Z))_{+}
        CJJ1XFJ_3=-(DLOG(C3/2.D0)+EULER)*SPCJFJ(IPDF,X,AMU)*S_ALPI
        CJJ2xFJ_3=(-2*beta1*A1c*log(b0*C2/C1)**2*log(C3/b0)
     >  +2*beta1*B1c*log(b0*C2/C1)*log(C3/b0)
     >  +1.0/32.0*A1c**3*log(b0**2*C2**2/C1**2)**4
     >  -1.0/12.0*beta1*A1c*log(b0**2*C2**2/C1**2)**3
     >  -1.0/8.0*A1c*B1c*log(b0**2*C2**2/C1**2)**3
     >  +B2c*log(b0*C2/C1)-A2c*log(b0*C2/C1)**2
     >  +1.0/2.0*(B1c*log(b0*C2/C1))**2
     >  +beta1*B1c*log(b0*C2/C1)**2)*FPDF
     >  +C1xf(IPDF,x,amu)*(B1c*log(b0*C2/C1)
     >     -A1c*log(b0*C2/C1)**2
     >     +beta1*log(C3**2/b0**2))
     >  -0.5*C1P1xF(IPDF,x,amu)*log(C3**2/b0**2)
     >  -0.25*splithoppet(2,ipdf,x,amu)*log(C3**2/b0**2)
     >  +1.0/8.0*A1c*splithoppet(11,ipdf,x,amu)*log(C3**2/b0**2)**2+
     >  splithoppet(1,ipdf,x,amu)*(
     >     -B1c*log(b0*C2/C1)*log(C3/b0)
     >     +A1c*log(b0*C2/C1)**2*log(C3/b0)
     >     -beta1*log(c3/b0)**2)
        CJJ2XFJ_3=CJJ2XFJ_3*s_ALPI**2
       ENDIF

       CJJ1XFJ=CJJ1XFJ_1+CJJ1XFJ_2+CJJ1XFJ_3

C (1/(1-z))_{+} piece
       CJJ2XFJ_2=(CA*CF*(7.D0*ZETA3/2.D0-101.D0/27.D0)
     > +CF*NF*(14.D0/27.D0))*SPCJFJ2(IPDF,X,AMU)*0.5*S_ALPI**2

       If (iCF1.EQ.1) then
C Used only for double check
         ConvQ = Sqrt(CJJ0XFJ*S_ALPI*2.d0*(CJJ1XFJ_1+CJJ1XFJ_3))
       Else
         CONVQ=CJJ0XFJ+CJJ2XFJ_2+CJJ1XFJ+CJJ2XFJ_3
       End If


C       IF(DEBUG) THEN
c        WRITE(NWRT,*)'X,AMU,S_ALPI'
c        WRITE(NWRT,*)X,AMU,S_ALPI
c        WRITE(NWRT,*)'CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ,CONVQ,IPDF'
c        WRITE(NWRT,*)CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ,CONVQ,IPDF
c        WRITE(NWRT,*)'ICF1 =',ICF1
C       ENDIF

      ENDIF
      END


C --------------------------------------------------------------------------
      FUNCTION CONVG(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONVG,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ4NT
      REAL*8 CJG1XFG,S_ALPI, CJG2XFG
      REAL*8 CJXFG, CJXFGH1
      EXTERNAL CJXFG, CJXFGH1

      ConvG = 0.d0
CsB iCF1 separates the alpha_s delta(1-Z) and LO + alpha_s (1-Z)+gluon
C   contributions.
C   Settings: iCF1 = 0 -> full contribution is returned in convq(...)
C             iCF1 = 1 -> alpha_s delta(1-Z) qqB contribution returned
C             iCF1 = 2 -> LO + alpha_s (1-Z) qqB + Gq contributions

CCPY Sept 2009 
      If (Type_V.Eq.'AA' .or. Type_V.Eq.'ZZ'
     &    .or.Type_V.Eq.'WW_UUB'.or.Type_V.Eq.'WW_DDB') then 
        If (iCF1.Eq.1) Return
      End If

      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

CsB 9/26/96 This is not needed here. ConvG should not know about PDF(IPDF)
C   other than IPDF = 0.
C      IF(IHADRON.EQ.-2) THEN
CC FOR PION_MINUS-NUCLEUS SCATTERING
C        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
C      ELSE
C        FPDF=APDF(IPDF,X,AMU)
C      ENDIF

      IF(GOAROUND) THEN
       CONVG=0.D0
       RETURN
      ENDIF

      CJG1XFG=0.
      CJG2XFG=0.

C OCT. 17, 1995
C      IF(NORDER.EQ.1) THEN
C       CONVG=0.D0
C
C      ELSE IF(NORDER.EQ.2) THEN
CCPY
C      IF(NORDER.EQ.1 .OR. NORDER.EQ.2) THEN

      IF(N_WIL_C.EQ.0) THEN
       CONVG=0.D0
       RETURN

      ELSEIF(N_WIL_C.GE.1) THEN

       UPLIM=1.0
       DNLIM=X
CJI Lower accuracy
       AERR=ACC_AERR
       RERR=ACC_RERR
C       AERR=1E-8
C       RERR=1E-1
       IACTA=2
       IACTB=2

C FOR GLUONS
C       IF(LEPASY.EQ.1 .AND. NOT(TESTASY)) THEN
C        CJG1XFG=0.D0
C       ELSE
C DO Z*(1-Z)/2
       CJG1XFG=ADZ4NT(CJXFG,DNLIM,UPLIM,AERR,RERR,ERREST,
     >      IER,IACTA,IACTB)

        IF(IER.NE.0) THEN
         WRITE(NWRT,*)' ERROR IN CJXFG'
CCPY         CALL QUIT
        ENDIF
C       ENDIF

       S_ALPI=ALPI(AMU)
       CONVG=CJG1XFG!*S_ALPI!-0.5*CJG2XFG*S_ALPI
C        print*, CJG1XFG, -CJG2XFG*0.5*S_ALPI,CONVG

C       IF(DEBUG) THEN
C        Print*,'X,AMU,S_ALPI'
C        Print*,X,AMU,S_ALPI
C        Print*,'CJG1XFG,CONVG,IPDF'
C        Print*, CJG1XFG,CONVG,IPDF
C       ENDIF
C        IF(N_WIL_C.GE.2) THEN
C           CJG2XFG_1=ADZ4NT(CJXFG2,DNLIM,UPLIM,AERR,RERR,ERREST,
C     >          IER,IACTA,IACTB)
C
C           CJG2XFG=CJG2XFG_1
C
C           CONVG=CONVG+S_ALPI**2*CJG2XFG
C        ENDIF
      ENDIF
      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CONVQQP(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONVQQP,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ4NT
      REAL*8 CJQQP1XFQ,S_ALPI
      REAL*8 CJXFQQP
      EXTERNAL CJXFQQP

      ConvQQP = 0.d0
      
      If (Type_V.Eq.'AA' .or. Type_V.Eq.'ZZ'
     &    .or.Type_V.Eq.'WW_UUB'.or.Type_V.Eq.'WW_DDB') then 
        If (iCF1.Eq.1) Return
      End If

      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      IF(GOAROUND) THEN
       CONVQQP=0.D0
       RETURN
      ENDIF

      CJQQP1XFQ=0.

      IF(N_WIL_C.lt.2) THEN
       CONVQQP=0.D0
       RETURN

      ELSEIF(N_WIL_C.Eq.2) THEN

       UPLIM=1.0
       DNLIM=X
CJI Lower accuracy
       AERR=ACC_AERR
       RERR=ACC_RERR
C        AERR=1E-1
C        RERR=1E-1
       IACTA=2
       IACTB=2

       CJQQP1XFQ=ADZ4NT(CJXFQQP,DNLIM,UPLIM,AERR,RERR,ERREST,
     >      IER,IACTA,IACTB)

        IF(IER.NE.0) THEN
         WRITE(NWRT,*)' ERROR IN CJXFG'
CCPY         CALL QUIT
        ENDIF

       S_ALPI=ALPI(AMU)
       CONVQQP=CJQQP1XFQ*S_ALPI**2
      ENDIF
      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CONVQQ(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONVQQ,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ4NT
      REAL*8 CJQQ1XFQ,S_ALPI
      REAL*8 CJXFQQ
      EXTERNAL CJXFQQ

      ConvQQ = 0.d0
      
      If (Type_V.Eq.'AA' .or. Type_V.Eq.'ZZ'
     &    .or.Type_V.Eq.'WW_UUB'.or.Type_V.Eq.'WW_DDB') then 
        If (iCF1.Eq.1) Return
      End If

      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      IF(GOAROUND) THEN
       CONVQQ=0.D0
       RETURN
      ENDIF

      CJQQ1XFQ=0.

      IF(N_WIL_C.lt.2) THEN
       CONVQQ=0.D0
       RETURN

      ELSEIF(N_WIL_C.Eq.2) THEN

       UPLIM=1.0
       DNLIM=X
CJI Lower accuracy
       AERR=ACC_AERR
       RERR=ACC_RERR
C        AERR=1E-1
C        RERR=1E-1
       IACTA=2
       IACTB=2

       CJQQ1XFQ=ADZ4NT(CJXFQQ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >      IER,IACTA,IACTB)

        IF(IER.NE.0) THEN
         WRITE(NWRT,*)' ERROR IN CJXFG'
CCPY         CALL QUIT
        ENDIF

       S_ALPI=ALPI(AMU)
       CONVQQ=CJQQ1XFQ*S_ALPI**2
      ENDIF
      RETURN
      END


C --------------------------------------------------------------------------
      FUNCTION CONVGG(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONVGG,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ4NT
      REAL*8 CJGG1XFG,S_ALPI
      REAL*8 CJXFGG
      EXTERNAL CJXFGG

      ConvGG = 0.d0
      
      If (Type_V.Eq.'AA' .or. Type_V.Eq.'ZZ'
     &    .or.Type_V.Eq.'WW_UUB'.or.Type_V.Eq.'WW_DDB') then 
        If (iCF1.Eq.1) Return
      End If

      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      IF(GOAROUND) THEN
       CONVGG=0.D0
       RETURN
      ENDIF

      CJGG1XFG=0.

      IF(N_WIL_C.lt.2) THEN
       CONVGG=0.D0
       RETURN

      ELSEIF(N_WIL_C.Eq.2) THEN

       UPLIM=1.0
       DNLIM=X
CJI Lower accuracy
       AERR=ACC_AERR
       RERR=ACC_RERR
C        AERR=1E-1
C        RERR=1E-1
       IACTA=2
       IACTB=2

       CJGG1XFG=ADZ4NT(CJXFGG,DNLIM,UPLIM,AERR,RERR,ERREST,
     >      IER,IACTA,IACTB)

        IF(IER.NE.0) THEN
         WRITE(NWRT,*)' ERROR IN CJXFG'
CCPY         CALL QUIT
        ENDIF

       S_ALPI=ALPI(AMU)
       CONVGG=CJGG1XFG*S_ALPI**2
      ENDIF
      RETURN
      END

CsB This function is not active as 3/96.
C   The function CONVQ(IPDF,X,AMU) replaces it.
C --------------------------------------------------------------------------
      FUNCTION CONV_SAVE(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CONV_SAVE,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT
      REAL*8 CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ_3,CJJ1XFJ,CJG1XFG,
     >S_ALPI
      REAL*8 CJXFJ,CJXFG,SPCJFJ
      EXTERNAL CJXFJ,CJXFG,SPCJFJ

      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
      ELSE
        FPDF=APDF(IPDF,X,AMU)
      ENDIF


      IF(GOAROUND) THEN
       CJJ0XFJ=FPDF
       CONV_SAVE=CJJ0XFJ
       RETURN
      ENDIF

      CJJ0XFJ=0.d0
      CJJ1XFJ_1=0.d0
      CJJ1XFJ_2=0.d0
      CJJ1XFJ_3=0.d0
      CJJ1XFJ=0.d0
      CJG1XFG=0.d0

C OCT. 17, 1995
C      IF(NORDER.EQ.1) THEN
C       CJJ0XFJ=FPDF
C       CONV_SAVE=CJJ0XFJ
C
C      ELSE IF(NORDER.EQ.2) THEN
CCPY
C      IF(NORDER.EQ.1 .OR. NORDER.EQ.2) THEN

      IF(N_WIL_C.EQ.0) THEN
       CJJ0XFJ=FPDF
       CONV_SAVE=CJJ0XFJ

      ELSEIF(N_WIL_C.EQ.1) THEN

C DO \delta(1-Z)
       CJJ0XFJ=FPDF

C DO \delta(1-Z)*CF*(PI^2-8)/4
        CJJ1XFJ_1=0.25*CF*(PI**2-8.0)*FPDF

C DO CF*(1-Z)/2
        UPLIM=1.0
        DNLIM=X
        AERR=ACC_AERR
        RERR=ACC_RERR
        IACTA=2
        IACTB=2

        CJJ1XFJ_2=ADZ3NT(CJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >  IER,IACTA,IACTB)
        IF(IER.NE.0) THEN
          WRITE(NWRT,*)' ERROR IN CJXFJ_2', CJJ1XFJ_2
CCPY          CALL QUIT
        ENDIF

CCPY April 16, 1996
C This is consistent with our PRD paper. 9/15/97
CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
C       IF(IFLAG_C3.NE.1) THEN
       IF(IFLAG_C3.GT.4) THEN
C ADD EXTRA TERMS
CsB Add CF*(9/16 - (ln(C1/(2 C2))-3/4+gamma_E)^2) (See Cfns.nb)
        CJJ1XFJ_1=CJJ1XFJ_1+CF*(9.D0/16.D0-
     >           (DLOG(C1/2.D0/C2)-3.D0/4.D0+EULER)**2)*FPDF
C DO -CF*(LN(C1/2)+EULER)*((1+Z^2)/(1-Z))_{+}
        CJJ1XFJ_3=-(DLOG(C3/2.D0)+EULER)*SPCJFJ(IPDF,X,AMU)
       ENDIF

       CJJ1XFJ=CJJ1XFJ_1+CJJ1XFJ_2+CJJ1XFJ_3

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C FOR GLUONS
C DO Z*(1-Z)/2
       CJG1XFG=ADZ4NT(CJXFG,DNLIM,UPLIM,AERR,RERR,ERREST,
     > IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CJXFG'
CCPY        CALL QUIT
       ENDIF

       S_ALPI=ALPI(AMU)
       CONV_SAVE=CJJ0XFJ+S_ALPI*(CJJ1XFJ+CJG1XFG)

C       IF(DEBUG) THEN
C        WRITE(NWRT,*)'X,AMU,S_ALPI'
C        WRITE(NWRT,*)X,AMU,S_ALPI
C        WRITE(NWRT,*)'CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ,CJG1XFG,CONV,IPDF'
C        WRITE(NWRT,*)CJJ0XFJ,CJJ1XFJ_1,CJJ1XFJ_2,CJJ1XFJ,CJG1XFG,CONV,IPDF
C       ENDIF

      ELSE
       WRITE(*,*) ' NO SUCH N_WIL_C VALUE '
       CALL QUIT
      ENDIF
      RETURN
      END


C --------------------------------------------------------------------------
      FUNCTION CJXFJ(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CJXFJ,Z
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2,T3,Y
      DOUBLE PRECISION xLi, S_ALPI
      EXTERNAL xLi

C DO CF*(1-Z)/2
      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF_INT,Z,AMU_INT)
      ELSE
        FPDF=APDF(IPDF_INT,Z,AMU_INT)
      ENDIF

      S_ALPI=ALPI(AMU_INT)

      if(x_int .eq. z) then
          cjxfj=0
          return
      endif

      T1=1.0/Z
CC C(1) Piece
      T2=0.5*CF*(1.0-X_INT/Z)*S_ALPI
        If(N_WIL_C.EQ.2) THEN
          Y=X_INT/Z
          T3=CA*CF*((1.0+Y**2)/(1-Y)*(-xLi(3,1-y)/2.0
     >+xLi(3,y)-xLi(2,y)*log(y)/2d0
     >-0.5*xLi(2,y)*log(1-y)-1.0/24d0*log(y)**3-0.5*log(1-y)**2*log(y)
     >+1.0/12.0*PI2*log(1-y)-pi2/8d0)
     >+1.0/(1-y)*(-0.25*(11-3*y*y)*ZETA3
     >-1/48d0*(-y*y+12*y+11)*log(y)**2
     >-1.0/36.0*(83*y*y-36*y+29)*log(y)+Pi2*y/4d0)
     >+(1-y)*(xLi(2,y)/2d0+0.5*log(1-y)*log(y))
     >+(y+100)/27d0+0.25*y*log(1-y))
     >+CF*NF*((1+y*y)/(72*(1-y))*log(y)*(3.0*log(y)+10)+(-19*y-37)/108)
     >+CF**2*((1+y*y)/(1-y)*(xLi(3,1-y)/2d0+xLi(2,y)*log(1-y)/2d0
     >+3*xLi(2,y)*log(y)/2d0+0.75*log(y)*log(1-y)**2
     >+log(y)**2*log(1-y)/4d0-1/12d0*pi2*log(1-y))
     >+(1-y)*(-xLi(2,y)-1.5*log(1-y)*log(y)+2*pi2/3d0-29/4d0)
     >+1/24d0*(1+y)*log(y)**3
     >+1/(1-y)*(1/8d0*(-2*y*y+2*y+3)*log(y)**2
     >+0.25*(17*y*y-13*y+4)*log(y))-y/4d0*log(1-y))
     >+CF*CF*((1+y*y)/(1-y)*(5*ZETA3
     >-5*xLi(3,y))/2d0)
     >-CF*CF/4d0*((2*pi2-18)*(1-y)-(1+y)*log(y))     
     >+CF*(1/y*(1-y)*(2*y*y-y+2)*(xLi(2,y)/6+1/6d0*log(1-y)*log(y)
     >-pi2/36)
     >+1/(216*y)*(1-y)*(136*y*y-143*y+172)
     >-1/48d0*(8*y*y+3*y+3)*log(y)**2
     >+1d0/36d0*(32*y*y-30*y+21)*log(y)
     >+1/24d0*(1+y)*log(y)**3 )
       T3 =T3/2d0
        T3 = T3 + 0.25*CF*(1-y)*CF*(Pi2/2d0-4)
       T2 = T2 + T3*S_ALPI*S_ALPI
        ENDIF
      CJXFJ=T1*T2*FPDF

      if(CJXFJ .ne. cjxfj) then
          print*, "here", t1, t2, z, x_int
          call quit()
      endif

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFJCFG(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CJXFJCFG,Z
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT,X
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2,T3,TTOT, y
      DOUBLE PRECISION xLi, S_ALPI
      EXTERNAL xLi

C DO CF*(1-Z)/2
      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF_INT,Z,AMU_INT)
      ELSE
        FPDF=APDF(IPDF_INT,Z,AMU_INT)
      ENDIF

      X = X_INT/Z
C      if(x.eq.1) x=x-0.000001
      if(X.eq.1) then
          cjxfjcfg=0d0
          return
      endif
    

      S_ALPI=ALPI(AMU_INT)
!      print*, "In cjxfj: amu, s_alpi", amu_int, s_alpi

C      if(x_int .eq. z) then
C          cjxfj_CFG=0
C          return
C      endif

      T1=0
      T2=0
      T3=0

      T1=1.0/Z
CC C(1) Piece
      T2=0.5*CF*(1.0-X)
        If(N_WIL_C.EQ.2) THEN
CC C(2) Piece
          Y=X_INT/Z
          T3=CA*CF*((1.0-Y**2)/(1-Y)
     >*(-xLi(3,1-y)/2.0
     >+xLi(3,y)-xLi(2,y)*log(y)/2d0
     >-0.5*xLi(2,y)*log(1-y)-1.0/24d0*log(y)**3
     >-1.0/2.0*log(1-y)**2*log(y)
     >+1.0/12.0*PI2*log(1-y)-pi2/8d0)
     >+1.0/(1-y)*(-0.25*(11-3*y*y)*ZETA3
     >-1/48d0*(-y*y+12*y+11)*log(y)**2
     >-1.0/36.0*(83*y*y-36*y+29)*log(y)+Pi2*y/4d0)
     >+(1-y)*(xLi(2,y)/2.0
     >+log(1-y)*log(y)/2.0)
     >+(y+100d0)/27d0+1.0/4.0*y*log(1-y))
     >+CF*NF*((1+y*y)/(72*(1-y))*log(y)*(3.0*log(y)+10.0)
     >+1/108d0*(-19*y-37.0))
     >+CF**2*((1+y*y)/(1-y)*(xLi(3,1-y)/2d0+xLi(2,y)*log(1-y)/2d0
     >+3*xLi(2,y)*log(y)/2d0+0.75*log(y)*log(1-y)**2
     >+log(y)**2*log(1-y)/4d0-1/12d0*pi2*log(1-y))
     >+(1-y)*(-xLi(2,y)-1.5*log(1-y)*log(y)+2*pi2/3d0-29/4d0)
     >+1/24d0*(1+y)*log(y)**3
     >+1/(1-y)*(1/8d0*(-2*y*y+2*y+3)*log(y)**2
     >+0.25*(17*y*y-13*y+4)*log(y))-y/4d0*log(1-y))
     >+CF*CF*((1+y*y)/(1-y)*(5*ZETA3
     >-5*xLi(3,y))/2d0)
     >-CF*CF/4d0*((2*pi2-18)*(1-y)-(1+y)*log(y))
     >+CF*(1/y*(1-y)*(2*y*y-y+2)*(xLi(2,y)/6+1/6d0*log(1-y)*log(y)
     >-pi2/36)
     >+1/(216*y)*(1-y)*(136*y*y-143*y+172)
     >-1/48d0*(8*y*y+3*y+3)*log(y)**2
     >+1d0/36d0*(32*y*y-30*y+21)*log(y)
     >+1/24d0*(1+y)*log(y)**3 )
       T3 =T3/2d0
       T2 = T2 + T3*S_ALPI*S_ALPI
        ENDIF
        TTOT = T2*S_ALPI+T3*S_ALPI**2
      CJXFJCFG=T1*TTOT*FPDF

      if(CJXFJCFG .ne. cjxfjCFG) then
          print*, "here", t1, t2, z, x_int, t3,x
          stop
      endif

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFG(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      DOUBLE PRECISION CJXFG,Z
      INTEGER IPDF_INT
      DOUBLE PRECISION X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      DOUBLE PRECISION FPDF
      INTEGER IPDFG
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION CA_TERM, CF_TERM
      double precision xli, s_alpi
      external xli

cpn                    For heavy quarks
      double precision fc0, fc1, bM, bessk0
      external fc0,fc1, bessk0

cpn 
      double precision tmp0, tmp1, tmp2,t3,y
C DO Z*(1-Z)/2
      IPDFG=0

      S_ALPI=ALPI(AMU_INT)

      IF(IHADRON.EQ.2) THEN
cpn This function cannot be used for the outgoing hadron
        print *,'STOP: CJXFG cannot be used with ihadron = 2'
        call quit
      ELSE
        FPDF=APDF(IPDFG,x_int/Z,AMU_INT)
      ENDIF

C      print*, x_int/z, fpdf

      T1=1.0d0/Z
      T2=0.5d0*Z*(1.0d0-Z)


cpn                    Add a term that is proportional to P_{qG} and
cpn                    vanishes when C3 = b_0
CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
       IF(IFLAG_C3.NE.1)
C       IF(IFLAG_C3.GT.4)
     >   T2 = T2  - 0.5d0*(DLOG(C3/2.D0)+EULER) *(Z**2+(1.0d0-Z)**2)

cpn                    For heavy quarks, add a mass-dependent part
cpn                    If the massive heavy quark resummation is
cpn                     performed, and we are below the mass threshold
cpn                     of the heavy quark, the heavy quark PDF is zero.
cpn                     Here we put back the \log(mu/M_Q) piece that
cpn                     is subtracted from the parton cross section 
cpn                     in the process of finding the hard part above
cpn                     the threshold.
cpn                     The lines below restore the \log(mu/M_Q) piece 
cpn                     associated with the NLO quark C-functions. This
cpn                     uncompensated log term grows at large b, so that the
cpn                     calculation with the NLO quark C-functions
cpn                     should also include a NP Sudakov contribution to
cpn                     suppress this log term. 
cpn 
cpn                     The \log(mu/M_Q) piece  associated with the LO 
cpn                     quark C-functions is restored in the function
cpn                     CJXFG. Since the NLO gluon C-function is
cpn                     complete, the log term is properly cancelled at
cpn                     large b. Hence no NP Sudakov is needed. 
cpn
cpn                    Note: the scale in the PDF here is the "real
cpn                     scale" amu = b_0/b, not b_0/b_star

C===============================
CCPY HARD-WIRED IN main.for FOR HEAVY QUARK MASS EFFECT: 
C This is for 'HB" production, from b+bbar fusion.
c      IHQMASS=1
c      IHQPDF=5
ccpy: To avoid oscillating CSS output in large QT, we have to use 
C     the consistent value for the mass of bottom quark used in 
C     deriving the PDF, which is TMPVAR(8), not qmas_b.
C     Note that qmas_b is used to calculate the Yukawa coupling.
C             
C        XMHQ=TMPVAR(8)
C===============================

C      print*,'ihqmass,ihqpdf,xmhq,ipdf_int,b_int,amu_int ='
C      print*,ihqmass,ihqpdf,xmhq,ipdf_int,b_int,amu_int 
      
      if (ihqmass.eq.1) then
cpn      if (ihqmass.eq.1.and.abs(ipdf_int).eq.ihqpdf) then
cpn                    Note: bM is calculated using b, not b_star
        bM = b_int*xmhq

c        print*,' b_int,xmhq =',b_int,xmhq
        
cpn45                   
cpn                    If the massive heavy quark resummation is
cpn                     performed, and we are below the mass threshold
cpn                     of the heavy quark, the heavy quark PDF is zero.
cpn                     Here we put back the leading order 
cpn                     \log(mu/M_Q) piece that
cpn                     is subtracted from the parton cross section
cpn                     above the threshold  in the process of finding
cpn                     the hard part. 
cpn                     Note: the scale in the PDF here is the "real
cpn                     scale" amu = b_0/b, not b_0/b_star

        if (amu_int.ge.xmhq) then
cpn2005          tmp0 = fc0(bM)
          tmp0 = fc0(bM) - log(amu_int*b_int/b0)
        else
          tmp0 = bessk0(bM)
        endif

        tmp1 = fc1(bM) 

c        print*,' tmp0,tmp1 =',tmp0,tmp1

        T2 = T2 + 0.5d0*Z*(1.0d0-Z)*tmp1
     >    + 0.5d0*(Z**2+(1.0d0-Z)**2)*tmp0

      endif 

      t2=t2*S_ALPI

      if(N_WIL_C.GE.2) then
            y=z
      T3=CA*(-1.0/(12*y)*(1-y)*(11*y*y-y+2)*xLi(2,1-y)
     >+(2*y*y-2*y+1)*(xLi(3,1-y)/8-xLi(2,1-y)*dlog(1-y)/8
     >+dlog(1-y)**3/48d0)
     >+(2*y*y+2*y+1)*(3*xLi(3,-y)/8+xLi(3,1/(1+y))/4d0
     >-xLi(2,-y)*dlog(y)/8-dlog(1+y)**3/24+dlog(y)**2*dlog(1+y)/16
     >+pi2/48*dlog(1+y))
     >+0.25*y*(1+y)*xLi(2,-y)+y*xLi(3,y)-0.5*y*xLi(2,1-y)*dlog(y)
     >-y*xLi(2,y)*dlog(y)-3.0/8.0*(2*y*y+1)*ZETA3-149*y*y/216
     >-(44*y*y-12*y+3)/96*dlog(y)**2
     >+(68*y*y+6*pi2*y-30*y+21)/72*dlog(y)+pi2*y/24+43*y/48+43/(108*y)
     >+(2*y+1)/48*dlog(y)**3-0.5*y*dlog(1-y)*dlog(y)**2
     >-(1-y)/8*y*dlog(1-y)**2
     >+0.25*y*(1+y)*dlog(1+y)*dlog(y)+(3-4*y)/16*y*dlog(1-y)-35/48d0)
     >+CF*((2*y*y-2*y+1)*(ZETA3-xLi(3,1-y)/8-xLi(3,y)/8
     >+xLi(2,1-y)*dlog(1-y)/8+xLi(2,y)*dlog(y)/8-dlog(1-y)**3/48
     >+dlog(y)*dlog(1-y)**2/16+dlog(y)**2*dlog(1-y)/16)-3*y*y/8
     >-(4*y*y-2*y+1)/96*dlog(y)**3+(-8*y*y+12*y+1)/64*dlog(y)**2
     >+(-8*y*y+23*y+8)/32*dlog(y)+5*pi2/24*(1-y)*y+11*y/32
     >+(1-y)*y/8*dlog(1-y)**2-0.25*(1-y)*y*dlog(1-y)*dlog(y)
     >-(3-4*y)*y/16*dlog(1-y)-9/32d0)
     >-CF/4.0*(y*dlog(y)+0.5*(1-y*y)+(pi2-8)*y*(1-y))
     >+0.25*y*(1-y)*CF*(pi2/2d0-4)

      t2 = t2+t3*s_alpi**2
      endif

      tmp2 = t1*t2
      if (abs(tmp2).le.1e-16) tmp2 = 0d0
      CJXFG=tmp2*FPDF

C      print*, "x, z, vpdf"
C      print*, x_int, z, fpdf
C      print*, "t1, t2"
C      print*, t1, t2

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFGCFG(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      DOUBLE PRECISION CJXFGCFG,Z
      INTEGER IPDF_INT
      DOUBLE PRECISION X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      DOUBLE PRECISION FPDF
      INTEGER IPDFG
      DOUBLE PRECISION T1,T2,T3,TTOT
      DOUBLE PRECISION CA_TERM, CF_TERM
      double precision xli, s_alpi,y
      external xli

cpn                    For heavy quarks
      double precision fc0, fc1, bM, bessk0
      external fc0,fc1, bessk0

cpn 
      double precision tmp0, tmp1, tmp2
C DO Z*(1-Z)/2
      IPDFG=0
      T1=0
      T2=0
      T3=0
      TTOT=0

      S_ALPI=ALPI(AMU_INT)

      IF(IHADRON.EQ.2) THEN
cpn This function cannot be used for the outgoing hadron
        print *,'STOP: CJXFGCFG cannot be used with ihadron = 2'
        call quit
      ELSE
        FPDF=APDF(IPDFG,x_int/Z,AMU_INT)
      ENDIF

C      print*, x_int/z, fpdf

      T1=1.0d0/Z
      T2=0.5d0*Z*(1.0d0-Z)


cpn                    Add a term that is proportional to P_{qG} and
cpn                    vanishes when C3 = b_0
CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
C       IF(IFLAG_C3.NE.1) THEN
       IF(IFLAG_C3.GT.4)
     >   T2 = T2  - 0.5d0*(DLOG(C3/2.D0)+EULER) *(Z**2+(1.0d0-Z)**2)

cpn                    For heavy quarks, add a mass-dependent part
cpn                    If the massive heavy quark resummation is
cpn                     performed, and we are below the mass threshold
cpn                     of the heavy quark, the heavy quark PDF is zero.
cpn                     Here we put back the \log(mu/M_Q) piece that
cpn                     is subtracted from the parton cross section 
cpn                     in the process of finding the hard part above
cpn                     the threshold.
cpn                     The lines below restore the \log(mu/M_Q) piece 
cpn                     associated with the NLO quark C-functions. This
cpn                     uncompensated log term grows at large b, so that the
cpn                     calculation with the NLO quark C-functions
cpn                     should also include a NP Sudakov contribution to
cpn                     suppress this log term. 
cpn 
cpn                     The \log(mu/M_Q) piece  associated with the LO 
cpn                     quark C-functions is restored in the function
cpn                     CJXFGCFG. Since the NLO gluon C-function is
cpn                     complete, the log term is properly cancelled at
cpn                     large b. Hence no NP Sudakov is needed. 
cpn
cpn                    Note: the scale in the PDF here is the "real
cpn                     scale" amu = b_0/b, not b_0/b_star

C===============================
CCPY HARD-WIRED IN main.for FOR HEAVY QUARK MASS EFFECT: 
C This is for 'HB" production, from b+bbar fusion.
c      IHQMASS=1
c      IHQPDF=5
ccpy: To avoid oscillating CSS output in large QT, we have to use 
C     the consistent value for the mass of bottom quark used in 
C     deriving the PDF, which is TMPVAR(8), not qmas_b.
C     Note that qmas_b is used to calculate the Yukawa coupling.
C             
C        XMHQ=TMPVAR(8)
C===============================

C      print*,'ihqmass,ihqpdf,xmhq,ipdf_int,b_int,amu_int ='
C      print*,ihqmass,ihqpdf,xmhq,ipdf_int,b_int,amu_int 
      
      if (ihqmass.eq.1) then
cpn      if (ihqmass.eq.1.and.abs(ipdf_int).eq.ihqpdf) then
cpn                    Note: bM is calculated using b, not b_star
        bM = b_int*xmhq

c        print*,' b_int,xmhq =',b_int,xmhq
        
cpn45                   
cpn                    If the massive heavy quark resummation is
cpn                     performed, and we are below the mass threshold
cpn                     of the heavy quark, the heavy quark PDF is zero.
cpn                     Here we put back the leading order 
cpn                     \log(mu/M_Q) piece that
cpn                     is subtracted from the parton cross section
cpn                     above the threshold  in the process of finding
cpn                     the hard part. 
cpn                     Note: the scale in the PDF here is the "real
cpn                     scale" amu = b_0/b, not b_0/b_star

        if (amu_int.ge.xmhq) then
cpn2005          tmp0 = fc0(bM)
          tmp0 = fc0(bM) - log(amu_int*b_int/b0)
        else
          tmp0 = bessk0(bM)
        endif

        tmp1 = fc1(bM) 

c        print*,' tmp0,tmp1 =',tmp0,tmp1

        T2 = T2 + 0.5d0*Z*(1.0d0-Z)*tmp1
     >    + 0.5d0*(Z**2+(1.0d0-Z)**2)*tmp0

      endif 

      if(N_WIL_C.GE.2) then
          y=z
      T3=CA*(-1/(12*y)*(1-y)*(11*y*y-y+2)*xLi(2,1-y)
     >+(2*y*y-2*y+1)*(xLi(3,1-y)/8-xLi(2,1-y)*log(1-y)/8
     >+log(1-y)**3/48d0)
     >+(2*y*y+2*y+1)*(3*xLi(3,-y)/8+xLi(3,1/(1+y))/4d0
     >-xLi(2,-y)*log(y)/8-log(1+y)**3/24+log(y)**2*log(1+y)/16
     >+pi2/48*log(1+y))
     >+0.25*y*(1+y)*xLi(2,-y)+y*xLi(3,y)-0.5*y*xLi(2,1-y)*log(y)
     >-y*xLi(2,y)*log(y)-3.0/8.0*(2*y*y+1)*ZETA3-149*y*y/216
     >-(44*y*y-12*y+3)/96*log(y)**2
     >+(68*y*y+6*pi2*y-30*y+21)/72*log(y)+pi2*y/24+43*y/48+43/(108*y)
     >+(2*y+1)/48*log(y)**3-0.5*y*log(1-y)*log(y)**2
     >-(1-y)/8*y*log(1-y)**2
     >+0.25*y*(1+y)*log(1+y)*log(y)+(3-4*y)/16*y*log(1-y)-35/48d0)
     
     >+CF*((2*y*y-2*y+1)*(ZETA3-xLi(3,1-y)/8-xLi(3,y)/8
     >+xLi(2,1-y)*log(1-y)/8+xLi(2,y)*log(y)/8-log(1-y)**3/48
     >+log(y)*log(1-y)**2/16+log(y)**2*log(1-y)/16)-3*y*y/8
     >-(4*y*y-2*y+1)/96*log(y)**3+(-8*y*y+12*y+1)/64*log(y)**2
     >+(-8*y*y+23*y+8)/32*log(y)+5*pi2/24*(1-y)*y+11*y/32
     >+(1-y)*y/8*log(1-y)**2-0.25*(1-y)*y*log(1-y)*log(y)
     >-(3-4*y)*y/16*log(1-y)-9/32d0)
     >-CF/4.0*(y*log(y)+0.5*(1-y*y)+(pi2-8)*y*(1-y))
      endif

      ttot=t2*S_ALPI+t3*S_ALPI**2

      tmp2 = t1*ttot
      if (abs(tmp2).le.1e-16) tmp2 = 0d0
      CJXFGCFG=tmp2*FPDF

C      print*, "x, z, vpdf"
C      print*, x_int, z, fpdf
C      print*, "t1, t2"
C      print*, t1, t2

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFQQP(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      DOUBLE PRECISION CJXFQQP,Z
      INTEGER IPDF_INT
      DOUBLE PRECISION X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      DOUBLE PRECISION FPDF
      INTEGER IPDFG
      DOUBLE PRECISION T1,T2
      double precision xli, s_alpi
      external xli

cpn 
      double precision tmp0, tmp1, tmp2

      IF(IHADRON.EQ.2) THEN
cpn This function cannot be used for the outgoing hadron
        print *,'STOP: CJXFG cannot be used with ihadron = 2'
        call quit
      ELSE
        FPDF=APDF(IPDF_INT,x_int/Z,AMU_INT)
      ENDIF

C      print*, x_int/z, fpdf

      T1=1.0d0/Z
      T2=CF*(1d0/(12d0*Z)*(1d0-Z)*(2d0*Z**2-Z+2d0)
     >  *(xLi(2,Z)+log(1d0-Z)*log(Z)-pi2/6d0)
     >  +(1d0/(432d0*Z))*(1d0-Z)*(136*Z**2-143*Z+172)
     >  +(1d0/48d0)*(1d0+Z)*log(Z)**3
     >  -(1d0/96d0)*(8d0*Z**2+3*Z+3)*log(z)**2
     >  +(1d0/72d0)*(32d0*Z**2-30*Z+21)*log(Z))

      tmp2 = t1*t2
      if (abs(tmp2).le.1e-16) tmp2 = 0d0
      CJXFQQP=tmp2*FPDF

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFQQ(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      DOUBLE PRECISION CJXFQQ,Z
      INTEGER IPDF_INT
      DOUBLE PRECISION X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      DOUBLE PRECISION FPDF
      INTEGER IPDFG
      DOUBLE PRECISION T1,T2
      double precision xli, s_alpi
      external xli

cpn 
      double precision tmp0, tmp1, tmp2

      IF(IHADRON.EQ.2) THEN
cpn This function cannot be used for the outgoing hadron
        print *,'STOP: CJXFG cannot be used with ihadron = 2'
        call quit
      ELSE
        FPDF=APDF(IPDF_INT,x_int/Z,AMU_INT)
      ENDIF

C      print*, x_int/z, fpdf

      T1=1.0d0/Z
      T2=CF*(CF-1d0/2d0*CA)*((1d0+Z**2)/(1d0+Z)*(3d0*xLi(3,-Z)/2d0
     >  +xLi(3,Z)+xLi(3,1d0/(1d0+Z))-xLi(2,-Z)*log(Z)/2d0
     >  -xLi(2,Z)*log(Z)/2d0-1d0/24d0*log(Z)**3-1d0/6d0*log(1d0+Z)**3
     >  +1d0/4d0*log(1d0+Z)*log(Z)**2+pi2/12d0*log(1d0+Z)-3d0*Zeta3/4d0)
     >  +(1d0-Z)*(xLi(2,Z)/2d0+1d0/2d0*log(1d0-Z)*log(Z)+15d0/8d0)
     >  -1d0/2d0*(1d0+Z)*(xLi(2,-Z)+log(Z)*log(1d0+Z))
     >  +pi2/24d0*(Z-3d0)+1d0/8d0*(11d0*Z+3d0)*log(Z))
     >  +CF*(1d0/(12d0*Z)*(1d0-Z)*(2d0*Z**2-Z+2d0)
     >  *(xLi(2,Z)+log(1d0-Z)*log(Z)-pi2/6d0)
     >  +1d0/(432d0*Z)*(1d0-Z)*(136*Z**2-143d0*Z+172D0)
     >  -1d0/96d0*(8d0*Z**2+3d0*Z+3d0)*log(z)**2
     >  +1d0/72d0*(32d0*z**2-30d0*z+21d0)*log(Z)
     >  +1d0/48d0*(1d0+z)*log(Z)**3)

      tmp2 = t1*t2
      if (abs(tmp2).le.1e-16) tmp2 = 0d0
      CJXFQQ=tmp2*FPDF

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFGG(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      DOUBLE PRECISION CJXFGG,Z
      INTEGER IPDF_INT
      DOUBLE PRECISION X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      DOUBLE PRECISION FPDF
      INTEGER IPDFG
      DOUBLE PRECISION T1,T2

cpn 
      double precision tmp0, tmp1, tmp2
C DO Z*(1-Z)/2
      IPDFG=0

      IF(IHADRON.EQ.2) THEN
cpn This function cannot be used for the outgoing hadron
        print *,'STOP: CJXFG cannot be used with ihadron = 2'
        call quit
      ELSE
        FPDF=APDF(IPDFG,x_int/Z,AMU_INT)
      ENDIF

C      print*, x_int/z, fpdf

      T1=1.0d0/Z
      T2=-Z/2d0*(1d0-Z+1d0/2d0*(1d0+Z)*Log(Z))

      tmp2 = t1*t2
      if (abs(tmp2).le.1e-16) tmp2 = 0d0
      CJXFGG=tmp2*FPDF

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION CJXFG_SAVE(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CJXFG_SAVE,Z
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      INTEGER IPDFG
      REAL*8 T1,T2

C DO Z*(1-Z)/2
      IPDFG=0

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDFG,Z,AMU_INT)
      ELSE
        FPDF=APDF(IPDFG,Z,AMU_INT)
      ENDIF

      T1=1.0/Z
      T2=0.5*(X_INT/Z)*(1.0-X_INT/Z)
      CJXFG_SAVE=T1*T2*FPDF

C This is consistent with our PRD paper. 9/15/97
CCPY JUNE 2013: 
C FOR IFALG_C3=1,2,3,4; WE HAVE C1=C2*B0 AND C3=B0=2*EXP(-EULER)
C       IF(IFLAG_C3.NE.1) THEN
       IF(IFLAG_C3.GT.4) THEN
       T2=-0.5*(DLOG(C3/2.D0)+EULER)*((X_INT/Z)**2+(1.0-X_INT/Z)**2)
       CJXFG_SAVE=CJXFG_SAVE+T1*T2*FPDF
      ENDIF

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION SPCJFJ(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 SPCJFJ,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/SPLIT_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ5NT
      REAL*8 PQQFQ_1,PQQFQ_2
      REAL*8 PQQXFQ_1,PQQXFQ_2
      EXTERNAL PQQXFQ_1,PQQXFQ_2

      SPCJFJ=0.D0
      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2

      PQQFQ_1=0.
      PQQFQ_2=0.

C DO THE SPLITTING OF q --> g q  OR  qbar --> g qbar
C FIRST PART:
       PQQFQ_1=ADZ5NT(PQQXFQ_1,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN PQQXFQ_1',IER,IER
CCPY        CALL QUIT
       ENDIF
C SECOND PART
      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
      ELSE
        FPDF=APDF(IPDF,X,AMU)
      ENDIF

      PQQFQ_2=FPDF*PQQXFQ_2(X)
C FINAL RESULT OF q --> g q  OR  qbar --> g qbar
      SPCJFJ=PQQFQ_1-PQQFQ_2

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION SPCJFJ2(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 SPCJFJ2,X,AMU
      INTEGER IPDF
      REAL*8 FPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/SPLIT_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ5NT
      REAL*8 PQQFQ_1,PQQFQ_2
      REAL*8 PQQXFQ_1,PQQXFQ_2,PQQXFQ_3,PQQXFQ_4
      EXTERNAL PQQXFQ_1,PQQXFQ_2,PQQXFQ_3,PQQXFQ_4

      SPCJFJ2=0.D0
      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2

      PQQFQ_1=0.
      PQQFQ_2=0.

C DO THE SPLITTING OF q --> g q  OR  qbar --> g qbar
C FIRST PART:
       PQQFQ_1=ADZ5NT(PQQXFQ_3,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN PQQXFQ_3',IER,IER
CCPY        CALL QUIT
       ENDIF
C SECOND PART
      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
      ELSE
        FPDF=APDF(IPDF,X,AMU)
      ENDIF

      PQQFQ_2=FPDF*PQQXFQ_4(X)
C FINAL RESULT OF q --> g q  OR  qbar --> g qbar
      SPCJFJ2=PQQFQ_1-PQQFQ_2

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION PQQXFQ_1(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PQQXFQ_1,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT
      COMMON/SPLIT_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF1,FPDF2
      REAL*8 T1,T2
      INTEGER IPDF

C DO THE FIRST PART OF SPLITTING q --> g q  OR  qbar --> g qbar
      IPDF=IPDF_INT
      IF(IPDF.EQ.0) THEN
       WRITE(NWRT,*) ' WRONG IPDF IN PQQXFQ_1'
       CALL QUIT
      ENDIF

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF1=APDF_PION_MINUS(IPDF,X_INT/Z,AMU_INT)
        FPDF2=APDF_PION_MINUS(IPDF,X_INT,AMU_INT)
      ELSE
        FPDF1=APDF(IPDF,X_INT/Z,AMU_INT)
        FPDF2=APDF(IPDF,X_INT,AMU_INT)
      ENDIF

      T1=FPDF1/Z-FPDF2
      T2=CF*(1.0+Z**2)/(1.0-Z)
      PQQXFQ_1=T1*T2

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION PQQXFQ_2(X)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PQQXFQ_2,X

C DO THE SECOND PART OF SPLITTING  q --> g q  OR  qbar --> g qbar
      PQQXFQ_2=CF*(-2.0*DLOG(1.0-X)-X-0.5*X**2)

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION PQQXFQ_3(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PQQXFQ_3,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT
      COMMON/SPLIT_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF1,FPDF2
      REAL*8 T1,T2
      INTEGER IPDF

C DO THE FIRST PART OF SPLITTING q --> g q  OR  qbar --> g qbar
      IPDF=IPDF_INT
      IF(IPDF.EQ.0) THEN
       WRITE(NWRT,*) ' WRONG IPDF IN PQQXFQ_1'
       CALL QUIT
      ENDIF

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF1=APDF_PION_MINUS(IPDF,X_INT/Z,AMU_INT)
        FPDF2=APDF_PION_MINUS(IPDF,X_INT,AMU_INT)
      ELSE
        FPDF1=APDF(IPDF,X_INT/Z,AMU_INT)
        FPDF2=APDF(IPDF,X_INT,AMU_INT)
      ENDIF

C      if(Z.ge.0.9999) then
C          PQQXFQ_3 = 0
C          return
C      endif

      T1=FPDF1/Z-FPDF2
      T2=(1.0)/(1.0-Z)
      PQQXFQ_3=T1*T2

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION PQQXFQ_4(X)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PQQXFQ_4,X

C DO THE SECOND PART OF SPLITTING  q --> g q  OR  qbar --> g qbar
      PQQXFQ_4=-log(1-x)

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE YMAXIMUM
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 X,ACOSH,T1,T2,T3

      ACOSH(X)=DLOG(X+DSQRT(X**2-1.0))

      T1=ECM2+Q_V**2
      T2=DSQRT(QT_V**2+Q_V**2)
      T3=2.0*ECM*T2
      YMAX=ACOSH(T1/T3)
CJI 2014: To handle high pt region for resummation
      IF(ISNAN(YMAX)) THEN
        YMAX=-9999
      ENDIF
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE HCXFCXF(B_SUD,CFCF)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 AMU,B_SUD,CFCF,S_BSTAR,TERM
      REAL*8 HCONV,HCONVG,HCONVGprime, HCONVQfuncG, HCONVGfuncG
      EXTERNAL HCONV,HCONVG,HCONVGprime, HCONVQfuncG, HCONVGfuncG

      REAL*8 Q0_PDF, q0pdf

      S_BSTAR=B_SUD/(DSQRT(1.0+(B_SUD/BMAX)**2))
      IF(IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     $     IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) THEN
C CANONICAL CHOICE FOR H0
        TERM=S_BSTAR
CCPY      Else If (IFLAG_C3.EQ.5) then
      Else
        TERM=S_BSTAR
      ENDIF
      AMU=C3/TERM


CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
CBY April, 2019
C       Q0_PDF=1.3
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF
CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
C      if(initlhapdf==1) then
C       Q0_PDF=q0pdf()
C      else
C       Q0_PDF=1.3D0
C      endif
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF

      IF(TYPE_V.EQ.'AG'.AND.LEPASY.EQ.1.AND.ICF1.EQ.0)THEN
        CFCF=HCONVGprime(X_A,AMU)*HCONVG(X_B,AMU)
     &      +HCONVG(X_A,AMU)*HCONVGprime(X_B,AMU)
      ELSE IF(TYPE_V.EQ.'AG'.AND.LEPASY.EQ.1.AND.ICF1.EQ.2)THEN
        CFCF=HCONVGprime(X_A,AMU)*HCONVGprime(X_B,AMU)
      ELSE
        CFCF=HCONV(0,X_A,AMU,S_BSTAR)*HCONV(0,X_B,AMU,S_BSTAR)
      ENDIF

      IF(TYPE_V.EQ.'H0')THEN
CZLdebug
	CFCF=CFCF+( HCONVQfuncG(X_A,AMU)+HCONVGfuncG(X_A,AMU) )
     &	         *( HCONVQfuncG(X_B,AMU)+HCONVGfuncG(X_B,AMU) )
      ENDIF

      RETURN
      END ! HCXFCXF

C --------------------------------------------------------------------------
      SUBROUTINE HJCXFCXF(B_SUD,CFCF)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'

      REAL*8 AMU,B_SUD,CFCF,S_BSTAR,TERM
      REAL*8 CFCF_1, CFCF_2
      REAL*8 S_ALPI, FPDF1, FPDF2
      REAL*8 H1, HardPart1, HJCONV
      REAL*8 D1s, R, t, ptj
      Common /HJ/ D1s, R, t, ptj
      External HardPart1, HJCONV
      integer iflavor, i
      common/flavor/iflavor

      REAL*8 Q0_PDF, q0pdf
      
      S_BSTAR=B_SUD/(DSQRT(1.0+(B_SUD/BMAX)**2))
      IF(IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     $     IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) THEN
C CANONICAL CHOICE FOR H0
        TERM=S_BSTAR
CCPY      Else If (IFLAG_C3.EQ.5) then
      Else
        TERM=S_BSTAR
      ENDIF
      AMU=C3/TERM


CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
CBY April,2019
C       Q0_PDF=1.3
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF
CCPY March 2019: Force AMU to be larger than the initial PDF scale Q0_PDF
C      if(initlhapdf==1) then
C       Q0_PDF=q0pdf()
C      else
C       Q0_PDF=1.3D0
C      endif
C       IF (AMU.LT.Q0_PDF) then
C         AMU=C3*DSQRT(1.0+(B_SUD*Q0_PDF/C3)**2)/B_SUD
C       ENDIF


      if(iflavor.eq.0) then

        CFCF_1 = HJCONV(0,X_A,AMU)*HJCONV(0,X_B,AMU)

      else
        do i=-5,5
          if(i.ne.0) then

            CFCF_1 = CFCF_1+HJCONV(0,X_A,AMU)*HJCONV(i,X_B,AMU)

          endif
        enddo
      endif

      CFCF = CFCF_1*(1-X_A)*(1-X_B)

C      S_ALPI = ALPI(AMU)

CJI: The b-dependent pieces
C      IF(LEPASY.EQ.0) THEN
C        CFCF_1 = FPDF1*FPDF2*(
C     >      (-(dlog(S_Q**2*B_SUD**2/C1**2))**2)+
C     >      (4d0*BETA1/CA-dlog(1d0/R**2))*dlog(S_Q**2*B_SUD**2/C1**2))
C
C        CFCF_2 = HJCONV(0,X_A,AMU)*FPDF2 
C     >         + HJCONV(0,X_B,AMU)*FPDF1
C        CFCF_2 = CFCF_2*dlog(C1**2/(B_SUD**2*AMU**2))
C
C        CFCF = (CFCF_1 + CFCF_2)
CCprint*, CFCF_1, CFCF_2
CCJI: The b-independent pieces
C      ELSE
C        CFCF = FPDF1*FPDF2
C      ENDIF
C
C      CFCF = CFCF*S_ALPI*CA/(2d0*pi)

C      CFCF = FPDF1*FPDF2

      RETURN
      END ! HCXFCXF

C --------------------------------------------------------------------------
      FUNCTION HJCONV(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HJCONV,X,AMU
      INTEGER IPDF, JPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ, HJCXF
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ, HJCXF
      REAL*8 SMALL, T1
      DATA SMALL/1.D-6/


      REAL*8 H1_C2
      
      HJCONV=0.D0
      X_INT=X
      AMU_INT=AMU

      S_ALPI = ALPI(AMU)

      T1 = 0
      IPDF_INT=IPDF
      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
C      AERR=1E-6
C      RERR=ACC_RERR
      IACTA=2
      IACTB=2

      FPDF=APDF(IPDF,X,AMU)

      T1 = T1 + FPDF/(1-X)
      T1=T1+S_ALPI*ADZ4NT(HJCXF,DNLIM,UPLIM,AERR,RERR,ERREST,
     >                    IER,IACTA,IACTB)

      HJCONV=T1
      !print*, hjconv

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION HJCXF(Z)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HJCXF,Z
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2
      DOUBLE PRECISION xLi, S_ALPI
      EXTERNAL xLi
      integer iflavor, jpdf
      common/flavor/iflavor

      HJCXF = 0

      if(ipdf_int.eq.0) then
          do jpdf = -3, 3
            if(jpdf .ne. 0) then 
                FPDF=APDF(IPDF_INT,Z,AMU_INT)
                HJCXF=HJCXF+FPDF/Z*(0.5*CF*X_INT/Z)
            endif
          enddo
      else
          FPDF = APDF(IPDF_INT,Z,AMU_INT)
          HJCXF = HJCXF+FPDF/Z*(2.0/3.0*(1-X_INT/Z))
          FPDF = APDF(0, Z, AMU_INT)
          HJCXF = HJCXF + FPDF/Z*(0.5*X_INT/Z*(1-X_INT/Z))
      endif

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION H0Approx(VMAS,vpt,rapin)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 H0Approx, VMAS,vpt,rapin,s,t,u,pj

      pj=VMAS

C       s=(2d0*(mh**2+pj**2)**(3d0/2d0)*(mh**2+2d0*pj**2)
C     >    +pj*(4d0*(mh**2+pj**2)**2+(2d0*mh**2+pj**2)*vpt**2)
C     >    *cosh(rapin))/(2d0*(mh**2+pj**2)**(3d0/2d0))
C
C
CC      s=mh**2+2d0*vmas*dsqrt(mh**2+vmas**2+vpt**2)*cosh(rapin)
CC     > +2d0*vmas**2+2d0*vpt**2
CC      t=-vmas*dsqrt(mh**2+vmas**2+vpt**2)*cosh(rapin)-vmas**2
CC      t=(mh**2-s)/2d0
C      t=-(pj*(4d0*pj*(mh**2+pj**2)**(3d0/2d0)*(mh**4-
C     >  (mh**2-2d0*pj**2)*vpt**2)+Exp(rapin)*(4d0*mh**4*
C     >  (mh**2+pj**2)**2+(-2d0*mh**6-3d0*mh**4*pj**2+8*mh**2*pj*4
C     >  +8*pj**6)*vpt**2)))/(4d0*mh**4*(mh**2+pj**2)**(3d0/2d0))
C      u=-exp(-rapin)*pj*(4d0*mh**8+8d0*pj**5*(pj+exp(rapin)*
C     >  dsqrt(mh**2+pj**2))*vpt**2+4d0*mh**2*pj**3*(2d0*pj
C     >  +Exp(rapin)*dsqrt(mh**2+pj**2))*vpt**2+mh**6*(8d0*pj**2
C     >  +4d0*exp(rapin)*pj*dsqrt(mh**2+pj**2)-2d0*vpt**2)
C     >  +mh**4*pj*(4d0*pj**3+4d0*exp(rapin)*pj**2*dsqrt(mh**2+pj**2)
C     >  -3d0*pj*vpt**2-4d0*exp(rapin)*dsqrt(mh**2+pj**2)*vpt**2))/
C     >  (4d0*mh**4*(mh**2+pj**2)**(3d0/2d0))
CC      print*, vmas, vpt, rapin, mh, x_a, x_b
CC      print*, s, t, u, dsqrt(s+t+u)
CC      if(sqrt(s) .le. sqrt(2d0)*mh) then
CC          H0Approx = 0
CC          Return
CC      endif

      s =mh**2+2d0*vmas**2
     >  +2d0*vmas*dsqrt(mh**2+vmas**2)*Cosh(rapin)

      t = -vmas*(vmas+dexp(0d0)*dsqrt(mh**2+vmas**2))

      u = -dexp(-0d0)*vmas*(dexp(0d0)*vmas+sqrt(mh**2+vmas**2))
C      print*, t, u

      H0Approx = (s**4+t**4+u**4+mh**8)/(s*t*u)*NC*(NC**2-1)

      if(H0Approx.lt.0) H0Approx=0

      Return
      END

C --------------------------------------------------------------------------
      FUNCTION H1Approx(VMAS,vpt,rapin)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 H1Approx, VMAS, xLi,vpt,rapin,s,u
      INTEGER NF_EFF
      REAL*8 D1s, R, xLi1, xLi2, xLi3,z,t,ptj,pj
      REAL*8 ttil, util
      REAL*8 S_ALPI, H0Approx
      External xLi, H0Approx
      Common /HJ/ D1s, R,t,ptj
      pj=VMAS

C      s=mh**2+2d0*vmas*dsqrt(mh**2+vmas**2+vpt**2)*cosh(rapin)
C     > +2d0*vmas**2+2d0*vpt**2
C       s=(2d0*(mh**2+pj**2)**(3d0/2d0)*(mh**2+2d0*pj**2)
C     >    +pj*(4d0*(mh**2+pj**2)**2+(2d0*mh**2+pj**2)*vpt**2)
C     >    *cosh(rapin))/(2d0*(mh**2+pj**2)**(3d0/2d0))
C
CC      t=-vmas*dsqrt(mh**2+vmas**2+vpt**2)*cosh(rapin)-vmas**2
CC      t=(mh**2-s)/2d0
CC      t=(mh**2-x_a*ECM*vmas-x_b*ECM*sqrt(mh**2+vmas**2+vpt**2))/2d0
CC      u=(mh**2-x_a*ECM*sqrt(mh**2+vmas**2+vpt**2)-x_b*ECM*vmas)/2d0
C      t=-(pj*(4d0*pj*(mh**2+pj**2)**(3d0/2d0)*(mh**4-
C     >  (mh**2-2d0*pj**2)*vpt**2)+Exp(rapin)*(4d0*mh**4*
C     >  (mh**2+pj**2)**2+(-2d0*mh**6-3d0*mh**4*pj**2+8*mh**2*pj*4
C     >  +8*pj**6)*vpt**2)))/(4d0*mh**4*(mh**2+pj**2)**(3d0/2d0))
C      u=-exp(-rapin)*pj*(4d0*mh**8+8d0*pj**5*(pj+exp(rapin)*
C     >  dsqrt(mh**2+pj**2))*vpt**2+4d0*mh**2*pj**3*(2d0*pj
C     >  +Exp(rapin)*dsqrt(mh**2+pj**2))*vpt**2+mh**6*(8d0*pj**2
C     >  +4d0*exp(rapin)*pj*dsqrt(mh**2+pj**2)-2d0*vpt**2)
C     >  +mh**4*pj*(4d0*pj**3+4d0*exp(rapin)*pj**2*dsqrt(mh**2+pj**2)
C     >  -3d0*pj*vpt**2-4d0*exp(rapin)*dsqrt(mh**2+pj**2)*vpt**2))/
C     >  (4d0*mh**4*(mh**2+pj**2)**(3d0/2d0))

      s =mh**2+2d0*vmas**2
     >  +2d0*vmas*dsqrt(mh**2+vmas**2)*Cosh(rapin)

      t = -vmas*(vmas+dexp(0d0)*dsqrt(mh**2+vmas**2))

      u = -dexp(-0d0)*vmas*(dexp(0d0)*vmas+sqrt(mh**2+vmas**2))

      S_ALPI = ALPI(dsqrt(s))

C      s=mh**2+2d0*vmas*sqrt(mh**2+vmas**2)+2d0*vmas**2
      NF_EFF= NFL(dsqrt(s))
C      t=-vmas*sqrt(mh**2+vmas**2)*cosh(rapin)-vmas**2
C      u=t
      ttil=mh**2-t
      util=mh**2-u

C      if(sqrt(s) .le. sqrt(2d0)*mh) then
C          H1Approx = 0
C          Return
C      endif

C      print*, s,t,u

      z=1d0-mh**2/s
      if(dabs(z).gt.1d0) then
          xLi1=-xLi(2,1d0/z)-pi2/6d0-0.5d0*dlog(-z)**2
      else
          xLi1=xLi(2,z)
      endif

      z=t/mh**2
      if(dabs(z).gt.1d0) then
          xLi2=-xLi(2,1d0/z)-pi2/6d0-0.5d0*dlog(-z)**2
      else
          xLi2=xLi(2,z)
      endif

      z=u/mh**2
      if(dabs(z).gt.1d0) then
          xLi3=-xLi(2,1d0/z)-pi2/6d0-0.5d0*dlog(-z)**2
      else
          xLi3=xLi(2,z)
      endif

C      print*, xLi1, Xli2,xLi3,nf_eff,R

C      H1Approx =16.3071-0.579959*NF_EFF+1.38629*dlog(0.5-mh**2/VMAS**2)
C     > +dlog(VMAS**2/(2d0*mh**2))**2
C     > -dlog(VMAS**2/(VMAS**2-2d0*mh**2))**2
C     > +dlog(1d0+(0.5*VMAS**2)/(mh**2))**2
C     > +2d0*Beta1*dlog(1d0/R**2)
C     > +2d0*xLi1+2d0*xLi2+2d0*xLi3

      H1Approx = ((dlog(s/vmas**2))**2+2d0*(11d0-2d0*NF_EFF/3d0)/12d0
     >   *dlog(s/(R**2*vmas**2))
     >   +dlog(1d0/R**2)*dlog(s/vmas**2)-2d0*dlog(-t/s)*dlog(-u/s)
     >   +dlog(ttil/mh**2)**2-dlog(-ttil/t)**2+dlog(util/mh**2)
     >   -dlog(-util/u)**2+2d0*xLi1+2d0*xLi2+2d0*xLi3+67d0/9d0
     >   +Pi2/2d0-23d0/54d0*NF_EFF)*S_ALPI*CA*0.5

      H1Approx = H1Approx +
     >  (NC-NF_EFF)*(1d0/4d0)*(2d0/3d0)*NC*(NC**2-1)*mh**2
     > /(s*t*u)*(s*t*u+mh**2*(s*t+s*u+t*u))
     > /H0Approx(VMAS,vpt,rapin)*S_ALPI


      if(H1Approx.ne.H1Approx) H1Approx=0


      Return
      END

C --------------------------------------------------------------------------
      FUNCTION HCONV_OLD(IPDF,X,AMU,BSTAR)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONV_OLD,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      REAL*8 SMALL
      DATA SMALL/1.D-6/


      REAL*8 H1_C2
      
      HCONV_OLD=0.D0
      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      IF(IPDF.NE.0) THEN
       WRITE(NWRT,*) ' IPDF IS NOT ZERO IN HCONV'
       CALL QUIT
      ENDIF

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
      ELSE
        FPDF=APDF(IPDF,X,AMU)
      ENDIF


      IF(GOAROUND) THEN
       CGG0XFG=FPDF
       HCONV_OLD=CGG0XFG
       RETURN
      ENDIF

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
CSM
      IF(N_WIL_C.EQ.0) THEN
       CGG0XFG=FPDF
       HCONV_OLD=CGG0XFG
       RETURN
      ENDIF

      CGG1XFG_1 = 0.d0
      CGG1XFG = 0.d0
      
CCPY April 16, 1996
C This is only correct for changing C3, but not for changing C2.
C      IF(IFLAG_C3.EQ.1) GOTO 200
        IF(IFLAG_C3.EQ.1 .OR. IFLAG_C3.EQ.2 .OR.
     $     IFLAG_C3.EQ.3 .OR. IFLAG_C3.EQ.4 ) GOTO 200

C THIS IS FOR IFLAG_C3.NE.1
C CHOOSE THE CANONICAL C3 WITH FULL EXPRESSION

C DO C_{g/g}
C DO \delta(1-Z) FOR ZEROTH ORDER
      CGG0XFG=FPDF

      TERM1=LOG(BSTAR*AMU/2.0)+EULER
      TERM2=LOG(C1/C2/2.0)+EULER

C DO \delta(1-Z)*(11+3*PI^2)/4 FOR FIRST ORDER
        CGG1XFG_1=FPDF*( (11.0D0+3.0D0*PI**2)/4.0D0
     >-3.0*TERM2**2+3.0*TERM2+(2.0*BETA1-3.0)*TERM1 )

C DO \log(..)*P_{g/g}
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DO THE SPLITTING OF g --> g g
C FIRST PART:
C FIRST PART OF 1/(1-Z)_+
      PGGFG_1=ADZ3NT(PGGXFG_1,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN PGGXFG_1'
CCPY        CALL QUIT
       ENDIF
C SECOND PART OF 1/(1-Z)_+
      PGGFG_1=6.0D0*(PGGFG_1+FPDF*DLOG(1.0D0-X))
C SECOND PART:
      PGGFG_2=ADZ4NT(PGGXFG_2,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN PGGXFG_2',IER
CCPY        CALL QUIT
       ENDIF
      PGGFG_2=PGGFG_2+2.0*BETA1*FPDF
C FINAL RESULT FOR g --> g g
      PGGFG=PGGFG_1+PGGFG_2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      CGG1XFG_2=-TERM1*PGGFG
      CGG1XFG=CGG1XFG_1+CGG1XFG_2

C DO C_{g/q}
C DO \log(...)*P_{g/q}

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DO THE SPLITTING OF (q, qbar) --> (q, qbar) g
C FOR ALL QUARKS AND ANTI-QUARKS
      PGJFJ=ADZ5NT(PGJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN PGJXFJ'
        CALL QUIT
       ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C FIRST PART:
      CGJ1XFJ_1=-TERM1*PGJFJ
C SECOND PART:
      CGJ1XFJ_2=ADZ6NT(CGJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ_2'
CCPY        CALL QUIT
       ENDIF

C FINAL RESULT FOR (q, qbar) --> (q, qbar) g
      CGJ1XFJ=CGJ1XFJ_1+CGJ1XFJ_2

      S_ALPI=ALPI(AMU)
      HCONV_OLD=CGG0XFG+S_ALPI*(CGG1XFG+CGJ1XFJ)

      IF(DEBUG) THEN
       WRITE(NWRT,*) ' CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG'
       WRITE(NWRT,*) ' CGG1XFG_2,CGG1XFG,PGJFJ,'
     $ ,'CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ'
       WRITE(NWRT,*) CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,
     >CGG1XFG_2,CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ
       WRITE(NWRT,*) ' HCONV_OLD,S_ALPI,AMU'
       WRITE(NWRT,*) HCONV_OLD,S_ALPI,AMU
       IF(HCONV_OLD.LT.0.D0) THEN
        WRITE(NWRT,*) ' HCONV_OLD < 0'
        CALL QUIT
       ENDIF
      ENDIF

      GOTO 999

200   CONTINUE

C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3

C DO C_{g/g}
C DO \delta(1-Z) FOR ZEROTH ORDER
      CGG0XFG=FPDF

      If (iCF1.NE.2) then 
CCPY        If (Type_V.Eq.'AG') then 
        If (Type_V.Eq.'AG'.OR.Type_V.Eq.'ZG') then 
CCPY PUT IN AN "EFFECTIVE" C(1) TERM FOR 'AA' PRODUCTION
          H1_C2=6.64D0
          CGG1XFG = H1_C2*FPDF 
	else
C DO \delta(1-Z)*(11+3*PI^2)/4 FOR FIRST ORDER
          CGG1XFG=FPDF*(11.0D0+3.0D0*PI**2)/4.0D0	
        ENDIF

      ENDIF    ! iCF1.NE.2
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CGJ1XFJ=ADZ6NT(CGJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ'
CCPY        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
      HCONV_OLD=CGG0XFG+S_ALPI*(CGG1XFG+CGJ1XFJ)

999   CONTINUE

      RETURN
      END

C --------------------------------------------------------------------------
      FUNCTION HCONV(IPDF,X,AMU,BSTAR)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONV,X,AMU,BSTAR
      INTEGER IPDF
      REAL*8 HCONVQ,HCONVG, HCONVQ2, HCONVG2, HCONVQG2
      EXTERNAL HCONVQ,HCONVG, HCONVQ2, HCONVG2, HCONVQG2

      IF(IPDF.NE.0) THEN
       WRITE(NWRT,*) ' IPDF IS NOT ZERO IN HCONV'
       CALL QUIT
      ENDIF

        HCONV = HCONVQ(X,AMU) + HCONVG(X,AMU)
CZL--------------------------------
CZL add C^(2) function for Higgs production
      if(Type_V.eq.'H0') then
CZLdebug
        HCONV = HCONV + HCONVQG2(X,AMU) 
      endif
CZL--------------------------------

      RETURN
      END


C --------------------------------------------------------------------------
      FUNCTION HCONVG(X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVG,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ
	double precision CGG2XFG
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      REAL*8 H1_C2
	double precision Lt
      
CCPY NOV 2009:
      IPDF=0

      HCONVG=0.D0
      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(IPDF,X,AMU)
      ELSE
        FPDF=APDF(IPDF,X,AMU)
      ENDIF


      IF(GOAROUND) THEN
       CGG0XFG=FPDF
       HCONVG=CGG0XFG
       RETURN
      ENDIF

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
CSM
      IF(N_WIL_C.EQ.0) THEN
       CGG0XFG=FPDF
       HCONVG=CGG0XFG
       RETURN
      ENDIF

      CGG1XFG_1 = 0.d0
      CGG1XFG = 0.d0
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3

C DO C_{g/g}
C DO \delta(1-Z) FOR ZEROTH ORDER
      CGG0XFG=FPDF

      If (iCF1.NE.2) then 
CCPY SEPT 2009        If (Type_V.Eq.'AG') then 
        If (Type_V.Eq.'AG'.OR.Type_V.Eq.'ZG') then 
CCPY PUT IN AN "EFFECTIVE" C(1) TERM FOR 'AA' PRODUCTION
          H1_C2=6.64D0
          CGG1XFG = H1_C2*FPDF 
	else if(Type_V.Eq.'H0') then
C DO \delta(1-Z)*(11+3*PI^2)/4 FOR FIRST ORDER
          CGG1XFG=FPDF*(11.0D0+3.0D0*PI**2)/4.0D0	
CZL add C^(2) for the second order C_gg function in Higgs production, here is the only delta function part
	  Lt=2*Log(MH/MT)
	  CGG2XFG=FPDF/2d0*(124.41658709962584 + 5.708333333333334*Lt)!(
!     &      (-5*CA)/96. - CF/12. + (9*CF**2)/4. + 
!     &      CA*CF*(-6.041666666666667 - (11*Lt)/8. - 
!     &      (3*Pi2)/4.) - (-3*CF + CA*(5 + Pi2))**2/16. + 
!     &      CA**2*(11.065972222222221 + (7*Lt)/8. + 
!     &      (157*Pi2)/72. + (13*Pi**4)/144. - (55*Zeta3)/18.) - 
!     &      CA*nf*(1.9930555555555556 + (5*Pi2)/36. + 
!     &      (4*Zeta3)/9.) + CF*nf*(-1.7083333333333333 + 
!     &      Lt/2. + Zeta3)
!     &      )
	else
C DO \delta(1-Z)*(11+3*PI^2)/4 FOR FIRST ORDER
          CGG1XFG=FPDF*(11.0D0+3.0D0*PI**2)/4.0D0	
        ENDIF

      ENDIF    ! iCF1.NE.2
      
      S_ALPI=ALPI(AMU)
CZLdebug
      HCONVG=CGG0XFG+S_ALPI*CGG1XFG+S_ALPI**2*CGG2XFG

999   CONTINUE

      RETURN
      END


C --------------------------------------------------------------------------
      FUNCTION HCONVQ(X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVQ,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      
      HCONVQ=0.D0
      
CCPY Feb 22, 2007: Adding the flag for N_WIL_C
      IF(N_WIL_C.EQ.0) RETURN

      X_INT=X
CCPY NOV 2009: IPDF IS NOT NEEDED FOR THIS FUNCTION. (SEE FUNCTION CGJXFJ(Z))
C      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CGJ1XFJ=ADZ6NT(CGJXFJ,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ'
!        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
      HCONVQ=S_ALPI*CGJ1XFJ

999   CONTINUE

      RETURN
      END

C --------------------------------------------------------------------------
CZL: G function for quark sum in Eq.(8) in 1106.4652
      FUNCTION HCONVQfuncG(X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVQfuncG,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJQfuncG
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJQfuncG
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      
      HCONVQfuncG=0.D0
      
CCPY Feb 22, 2007: Adding the flag for N_WIL_C
      IF(N_WIL_C.EQ.0) RETURN

      X_INT=X
CCPY NOV 2009: IPDF IS NOT NEEDED FOR THIS FUNCTION. (SEE FUNCTION CGJXFJ(Z))
C      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CGJ1XFJ=ADZ6NT(CGJXFJQfuncG,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ'
!        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
      HCONVQfuncG=S_ALPI*CGJ1XFJ

999   CONTINUE

      RETURN
      END

C --------------------------------------------------------------------------
CZL: G function for gluon flavor in Eq.(8) in 1106.4652
      FUNCTION HCONVGfuncG(X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVGfuncG,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJGfuncG
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJGfuncG
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      
      HCONVGfuncG=0.D0
      
CCPY Feb 22, 2007: Adding the flag for N_WIL_C
      IF(N_WIL_C.EQ.0) RETURN

      X_INT=X
CCPY NOV 2009: IPDF IS NOT NEEDED FOR THIS FUNCTION. (SEE FUNCTION CGJXFJ(Z))
C      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CGJ1XFJ=ADZ6NT(CGJXFJGfuncG,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ'
!        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
      HCONVGfuncG=S_ALPI*CGJ1XFJ

999   CONTINUE

      RETURN
      END

C --------------------------------------------------------------------------
CZL: P'xF
      FUNCTION HCONVGprime(X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVGprime,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,PGGFG,CGG1XFG,PGJFJ,CGJ1XFJ
      REAL*8 CGJXFJprime
      EXTERNAL CGJXFJprime
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      
      HCONVGprime=0.D0
      
CCPY Feb 22, 2007: Adding the flag for N_WIL_C
      IF(N_WIL_C.EQ.0) RETURN

      X_INT=X
CCPY NOV 2009: IPDF IS NOT NEEDED FOR THIS FUNCTION. (SEE FUNCTION CGJXFJ(Z))
C      IPDF_INT=IPDF
CZL      AMU_INT=C2*S_Q!AMU
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CGJ1XFJ=ADZ6NT(CGJXFJprime,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJprime'
CCPY        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
CZL      HCONVGprime=0.5*Log(AMU**2/(C2*S_Q)**2)*S_ALPI*CGJ1XFJ
      HCONVGprime=S_ALPI*CGJ1XFJ

999   CONTINUE

      RETURN
      END


C --------------------------------------------------------------------------
CZL HCONVQ2 copied from HCONVQ but for second order C_gq function C^(2)
      FUNCTION HCONVQG2(X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVQG2,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CQGG1XFJ2
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CQGGXFJ2
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CQGGXFJ2
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      
      HCONVQG2=0.D0
      
CCPY Feb 22, 2007: Adding the flag for N_WIL_C
      IF(N_WIL_C.EQ.0) RETURN

      X_INT=X
CCPY NOV 2009: IPDF IS NOT NEEDED FOR THIS FUNCTION. (SEE FUNCTION CGJXFJ(Z))
C      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CQGG1XFJ2=ADZ6NT(CQGGXFJ2,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ2'
	CQGG1XFJ2=0
CZL        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
      HCONVQG2=S_ALPI**2*CQGG1XFJ2

999   CONTINUE

      RETURN
      END



C --------------------------------------------------------------------------
      FUNCTION HCONVQ_ONE(IPDF,X,AMU)
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 HCONVQ_ONE,X,AMU,BSTAR
      INTEGER IPDF
      INTEGER IPDF_INT
      REAL*8 X_INT,AMU_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 DNLIM,UPLIM,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      REAL*8 ADZ3NT,ADZ4NT,ADZ5NT,ADZ6NT
      REAL*8 FPDF,S_ALPI
      REAL*8 CGG0XFG,CGG1XFG_1,PGGFG_1,PGGFG_2,PGGFG,CGG1XFG_2,
     >CGG1XFG,PGJFJ,CGJ1XFJ_1,CGJ1XFJ_2,CGJ1XFJ,CGJ1XFJ_ONE
      REAL*8 TERM1,TERM2
      REAL*8 PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ,CGJXFJ_ONE
      EXTERNAL PGGXFG_1,PGGXFG_2,PGJXFJ,CGJXFJ_2,CGJXFJ,CGJXFJ_ONE
      REAL*8 SMALL
      DATA SMALL/1.D-6/

      
      HCONVQ_ONE=0.D0
      X_INT=X
      IPDF_INT=IPDF
      AMU_INT=AMU

      UPLIM=1.0
      DNLIM=X
      AERR=ACC_AERR
      RERR=ACC_RERR
      IACTA=2
      IACTB=2
      
C THIS IS FOR IFLAG_C3=1
C THE CANONICAL C1, C2 AND C3
      
C DO C_{g/q}
C DO \log(...)*P_{g/q}
      CGJ1XFJ_ONE=ADZ6NT(CGJXFJ_ONE,DNLIM,UPLIM,AERR,RERR,ERREST,
     >IER,IACTA,IACTB)
       IF(IER.NE.0) THEN
        WRITE(NWRT,*)' ERROR IN CGJXFJ'
CCPY        CALL QUIT
       ENDIF

      S_ALPI=ALPI(AMU)
      HCONVQ_ONE=S_ALPI*CGJ1XFJ_ONE

999   CONTINUE

      RETURN
      END


      FUNCTION PGGXFG_1(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PGGXFG_1,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF1,FPDF2
      REAL*8 T1,T2
      INTEGER IPDF

C DO THE FIRST PART OF SPLITTING g --> g g
      IPDF=IPDF_INT
      IF(IPDF.NE.0) THEN
       WRITE(NWRT,*) ' WRONG IPDF IN PGGXFG_1'
       CALL QUIT
      ENDIF

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF1=APDF_PION_MINUS(IPDF,X_INT/Z,AMU_INT)
        FPDF2=APDF_PION_MINUS(IPDF,X_INT,AMU_INT)
      ELSE
        FPDF1=APDF(IPDF,X_INT/Z,AMU_INT)
        FPDF2=APDF(IPDF,X_INT,AMU_INT)
      ENDIF

      T1=FPDF1/Z-FPDF2
      T2=1.0D0/(1.0D0-Z)
      PGGXFG_1=T1*T2

      RETURN
      END

      FUNCTION PGGXFG_2(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PGGXFG_2,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2
      INTEGER IPDF

C DO THE SECOND PART OF SPLITTING g --> g g
      IPDF=IPDF_INT
      IF(IPDF.NE.0) THEN
       WRITE(NWRT,*) ' WRONG IPDF IN PGGXFG_2'
       CALL QUIT
      ENDIF
      FPDF=APDF(IPDF,X_INT/Z,AMU_INT)
      T1=FPDF/Z
      T2=6.0*(1.0/Z-2.0+Z*(1.0-Z))
      PGGXFG_2=T1*T2

      RETURN
      END


      FUNCTION PGJXFJ(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 PGJXFJ,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2
      INTEGER k

      PGJXFJ=0.D0
C DO THE SPLITTING (q, qbar) --> (q, qbar) g
      DO 10 K=-5,5
       IF(K.EQ.0) GOTO 10

       IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
         FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
       ELSE
         FPDF=APDF(K,X_INT/Z,AMU_INT)
       ENDIF

       T1=FPDF/Z
       T2=(4.0/3.0)*(1.0+(1.0-Z)**2)/Z
       PGJXFJ=PGJXFJ+T1*T2
10    CONTINUE

      RETURN
      END

      FUNCTION CGJXFJ_2(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CGJXFJ_2,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2

      CGJXFJ_2=0.D0
C DO THE SECOND PART OF C_{g/q} FOR (q, qbar) --> (q, qbar) g
      DO 10 K=-5,5
       IF(K.EQ.0) GOTO 10

       IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
         FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
       ELSE
         FPDF=APDF(K,X_INT/Z,AMU_INT)
       ENDIF

       T1=FPDF/Z
       T2=CF*(1.0-1.0/Z)
       CGJXFJ_2=CGJXFJ_2+T1*T2
10    CONTINUE

      RETURN
      END

      FUNCTION CGJXFJ(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CGJXFJ,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2

      CGJXFJ=0.D0
C DO C_{g/q} FOR (q, qbar) --> (q, qbar) g
C FOR CANONICAL C3
      DO 10 K=-5,5
       IF(K.EQ.0) GOTO 10

       IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
         FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
       ELSE
         FPDF=APDF(K,X_INT/Z,AMU_INT)
       ENDIF

       T1=FPDF/Z
       T2=CF*Z/2.0
       CGJXFJ=CGJXFJ+T1*T2
10    CONTINUE

      RETURN
      END

CZL: G function for quark in 1106.4652
      FUNCTION CGJXFJQfuncG(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CGJXFJQfuncG,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2

      CGJXFJQfuncG=0.D0
C DO C_{g/q} FOR (q, qbar) --> (q, qbar) g
C FOR CANONICAL C3
      DO 10 K=-5,5
       IF(K.EQ.0) GOTO 10

       IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
         FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
       ELSE
         FPDF=APDF(K,X_INT/Z,AMU_INT)
       ENDIF

       T1=FPDF/Z
       T2=CF*(1-Z)/Z
       CGJXFJQfuncG=CGJXFJQfuncG+T1*T2
10    CONTINUE

      RETURN
      END

CZL: G function for gluon in 1106.4652
      FUNCTION CGJXFJGfuncG(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CGJXFJGfuncG,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2

      CGJXFJGfuncG=0.D0
C DO C_{g/q} FOR (q, qbar) --> (q, qbar) g
C FOR CANONICAL C3
      DO 10 K=-5,5
       IF(K.NE.0) GOTO 10

       IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
         FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
       ELSE
         FPDF=APDF(K,X_INT/Z,AMU_INT)
       ENDIF

       T1=FPDF/Z
       T2=CA*(1-Z)/Z
       CGJXFJGfuncG=CGJXFJGfuncG+T1*T2
10    CONTINUE

      RETURN
      END

CZL: integrand kernel for P'xF
      FUNCTION CGJXFJprime(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CGJXFJprime,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2

      CGJXFJprime=0.D0
C DO C_{g/q} FOR (q, qbar) --> (q, qbar) g
C FOR CANONICAL C3

      IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
        FPDF=APDF_PION_MINUS(0,X_INT/Z,AMU_INT)
      ELSE
        FPDF=APDF(0,X_INT/Z,AMU_INT)
      ENDIF

      T1=FPDF/Z
      T2=2*CA*(1-Z)/Z
      CGJXFJprime=T1*T2

      RETURN
      END


CZL----------------------------------------------------
CZL copied from CGJXFJ2 for second order C_gg and C_qg function of non-delta function part
      FUNCTION CQGGXFJ2(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CQGGXFJ2,Z
        real*8 CQGXFJ2, CGGXFJ2
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF, FPDF0
      REAL*8 T1,T2, T10, T20
	external xLi
	double precision xLi
        double precision z2, z3, logz, log1mz, log1pz
        double precision logz2, log1mz2, log1pz2
        double precision logz3, log1mz3, log1pz3
        double precision poly1, poly2, poly3, poly4
        double precision poly5, poly6, poly7, poly8
        double precision temp

      CQGGXFJ2=0.D0
      CQGXFJ2=0.D0
      CGGXFJ2=0.D0

      z2=z**2
      z3=z**3
      logz=log(z)
      log1mz=log(1-z)
      log1pz=log(1+z)
      log1mz3=log1mz**3  
      logz2=logz**2  
      logz3=logz**3  
      log1pz3=log1pz**3  
      log1mz2=log1mz**2  
      log1pz2=log1pz**2 
      poly1=xLi(3,1/(1 + z))
      poly2=xLi(3,z)
      poly3=xLi(3,z2)
      poly4=xLi(3,z/(1 + z))
      poly5=xLi(3,-z)
      poly6=xLi(2,-z)
      poly7=xLi(2,z)
      poly8=xLi(2,z2)

C DO C_{g/q} FOR (q, qbar) --> (q, qbar) g
C FOR CANONICAL C3
      DO 10 K=-5,5
       IF(K.NE.0) THEN 

         IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
           FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
         ELSE
           FPDF=APDF(K,X_INT/Z,AMU_INT)
         ENDIF

         T1=FPDF/Z
         T2=(0.0030864197530864196*(-4451.636657168614 - 198.*log1mz2 - 
     &    60.*log1mz3 - 1065.9172753176506*log1pz + 216.*log1pz3 + 
     &      324.*log1mz2*logz - 324.*log1pz*logz2 + 4468.330898729398*z 
     &   +198.*log1mz2*z + 60.*log1mz3*z - 1065.9172753176506*log1pz*z+ 
     &      216.*log1pz3*z - 5022.*logz*z - 324.*log1mz2*logz*z + 
     &    900.*logz2*z - 324.*log1pz*logz2*z - 84.*logz3*z - 
     &      648.*poly2*(2. + (-4. + z)*z) - 648.*poly1*(2. + z*(2. + z))
     &   -243.*poly3*(2. + z*(2. + z)) + 
     &      778.9328732474171*(8. + z*(-5. + 4.*z)) + 
     &    266.7931882941266*z2 + 126.*log1mz2*z2 - 30.*log1mz3*z2 - 
     &      532.9586376588253*log1pz*z2 + 108.*log1pz3*z2 + 
     &    288.*logz*z2 + 162.*log1mz2*logz*z2 + 189.*logz2*z2 + 
     &    162.*log1pz*logz2*z2 - 
     &      66.*logz3*z2 - 3.*(45.*log1mz2 + 108.*log1pz*(-1. + logz)*
     &    logz + 8.*(-46.39118679673193 + 66.*logz - 9.*logz2)*z)*z2 + 
     &      162.*poly8*(logz*(2. + z*(2. + z)) + z2) + 
     &    108.*poly7*(-22. + 6.*logz*(2. + (-2. + z)*z) - 
     &    9.*z2 + 4.*z*(6. + z2)) + 
     &    6.*log1mz*(-164. + 4.*z*(41. + 20.*z) + 
     &    27.*logz2*(2. + (-2. + z)*z) - 129.*z2 + 
     &    18.*logz*(-22. - 9.*z2 + 4.*z*(6. + z2)))))/z

         CQGXFJ2=CQGXFJ2+T1*T2

       ELSE

         IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
           FPDF0=APDF_PION_MINUS(K,X_INT,AMU_INT)
           FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
         ELSE
           FPDF0=APDF(K,X_INT,AMU_INT)
           FPDF=APDF(K,X_INT/Z,AMU_INT)
         ENDIF

         T10=FPDF0
         T1=FPDF/Z
         T20=11.97590356063833!(14*CA*nf)/27.+CA**2*(-3.740740740740741+(7*Zeta3)/2.)
         T2=0.013888888888888888*(11372. + 1758.*logz + 270.*logz2 + 
     &    40.*logz3 - 8350./z + (648.*logz)/z - 9860.*z - 36.*log1mz*z+ 
     &    1020.*logz*z + 150.*logz2*z + 40.*logz3*z + (108.*logz3*(1.+ 
     &    z - z2)**2)/(-1. + z2) + 8710.*z2 - (80.*z2)/z + 
     &    (1065.9172753176506*log1pz*(1. + z + z2)**2)/(z*(1. + z)) - 
     &    (216.*log1pz3*(1. + z + z2)**2)/(z*(1. + z)) + 
     &    (1296.*poly4*(1. + z + z2)**2)/(z*(1. + z)) - (648.*poly5*
     &    (1. + z + z2)**2)/(z*(1. + z)) + 
     &    (648.*logz*poly6*(1. + z + z2)**2)/(z*(1. + z)) + 
     &    (355.3057584392169*(-1. + z)*(z - 11.*(1. + z2)))/z - 
     &    (216.*log1mz*logz*(-1. + z)*(z - 11.*(1. + z2)))/z + 
     &    27.*logz2*(25. - 11.*z + 44.*z2 - (12.*log1mz*(1. - z + 
     &    z2)**2)/((-1. + z)*z) - (12.*log1pz*(1. + z + z2)**2)/
     &    (z*(1. + z))) - 
     &    (9.*logz*(72. + 773.*z + 149.*z2 + (36.*log1mz2*
     &    (1. - z + z2)**2)/(-1. + z) - (72.*log1pz2*
     &    (1. + z + z2)**2)/(1. + z) + 
     &    536.*z3))/z - (648.*poly2*(-5. + z + 5.*z2**2 - 
     &    z3 - z2*(5. + z3)))/(z*(-1. + z2)) - 
     &    (389.46643662370855*(17.*z2 + z*(-21. + 
     &    2.*z2*(11. + 6.*z2) + 10.*z3) + 2.*(6. - 4.*z3 + 
     &    z2*(-6. - 11.*z2 + z3))))/
     &     (z*(-1. + z2)) - (216.*poly7*((-1. + z)**2*(1. + z)*
     &    (z - 11.*(1. + z2)) + 3.*logz*(3. - z + z3 + 
     &    z2*(3. - 3.*z2 + z3))))/
     &     (z*(-1. + z2)))

         CGGXFJ2 = CGGXFJ2 + T1*T2/2d0 + (T1-T10)*T20/2d0/(1-z)
       ENDIF
10    CONTINUE

        CQGGXFJ2 = CQGXFJ2 + CGGXFJ2

      RETURN
      END
CZL----------------------------------------------------

      FUNCTION CGJXFJ_ONE(Z)
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 CGJXFJ_ONE,Z
      REAL*8 X_INT,AMU_INT
      INTEGER IPDF_INT,k
      COMMON/CONV_INT/ X_INT,AMU_INT,IPDF_INT
      REAL*8 FPDF
      REAL*8 T1,T2

      CGJXFJ_ONE=0.D0
C DO C_{g/q} FOR (q, qbar) --> (q, qbar) g
C FOR CANONICAL C3

       K=IPDF_INT
      
       IF(IHADRON.EQ.-2) THEN
C FOR PION_MINUS-NUCLEUS SCATTERING
         FPDF=APDF_PION_MINUS(K,X_INT/Z,AMU_INT)
       ELSE
         FPDF=APDF(K,X_INT/Z,AMU_INT)
       ENDIF

       T1=FPDF/Z
       T2=CF*Z/2.0
       CGJXFJ_ONE=T1*T2

      RETURN
      END


C THIS IS C312.FOR

      FUNCTION BESJY(X)

      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL L

      ENTRY BESJ0(X)

      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 4
    8 F=0.0625*X**2-2.0
      A =           - 0.00000 00000 000008
      B = F * A     + 0.00000 00000 000413
      A = F * B - A - 0.00000 00000 019438
      B = F * A - B + 0.00000 00000 784870
      A = F * B - A - 0.00000 00026 792535
      B = F * A - B + 0.00000 00760 816359
      A = F * B - A - 0.00000 17619 469078
      B = F * A - B + 0.00003 24603 288210
      A = F * B - A - 0.00046 06261 662063
      B = F * A - B + 0.00481 91800 694676
      A = F * B - A - 0.03489 37694 114089
      B = F * A - B + 0.15806 71023 320973
      A = F * B - A - 0.37009 49938 726498
      B = F * A - B + 0.26517 86132 033368
      A = F * B - A - 0.00872 34423 528522
      A = F * A - B + 0.31545 59429 497802
      BESJY=0.5*(A-B)
      IF(L) RETURN

      A =           + 0.00000 00000 000016
      B = F * A     - 0.00000 00000 000875
      A = F * B - A + 0.00000 00000 040263
      B = F * A - B - 0.00000 00001 583755
      A = F * B - A + 0.00000 00052 487948
      B = F * A - B - 0.00000 01440 723327
      A = F * B - A + 0.00000 32065 325377
      B = F * A - B - 0.00005 63207 914106
      A = F * B - A + 0.00075 31135 932578
      B = F * A - B - 0.00728 79624 795521
      A = F * B - A + 0.04719 66895 957634
      B = F * A - B - 0.17730 20127 811436
      A = F * B - A + 0.26156 73462 550466
      B = F * A - B + 0.17903 43140 771827
      A = F * B - A - 0.27447 43055 297453
      A = F * A - B - 0.06629 22264 065699
      BESJY=0.636619772367581*LOG(X)*BESJY+0.5*(A-B)
      RETURN

    4 continue
C      print *,'C312: x= ',x
      F=256.0/X**2-2.0
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000051
      B = F * A - B + 0.00000 00000 000433
      A = F * B - A - 0.00000 00000 004305
      B = F * A - B + 0.00000 00000 051683
      A = F * B - A - 0.00000 00000 786409
      B = F * A - B + 0.00000 00016 306465
      A = F * B - A - 0.00000 00517 059454
      B = F * A - B + 0.00000 30751 847875
      A = F * B - A - 0.00053 65220 468132
      A = F * A - B + 1.99892 06986 950373
      P=A-B
      B =           - 0.00000 00000 000006
      A = F * B     + 0.00000 00000 000043
      B = F * A - B - 0.00000 00000 000334
      A = F * B - A + 0.00000 00000 003006
      B = F * A - B - 0.00000 00000 032067
      A = F * B - A + 0.00000 00000 422012
      B = F * A - B - 0.00000 00007 271916
      A = F * B - A + 0.00000 00179 724572
      B = F * A - B - 0.00000 07414 498411
      A = F * B - A + 0.00006 83851 994261
      A = F * A - B - 0.03111 17092 106740
      Q=8.0*(A-B)/V
      F=V-0.785398163397448
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 6
      BESJY=F*(Q*A+P*B)
      RETURN
    6 BESJY=F*(P*A-Q*B)
      RETURN

      ENTRY BESJ1(X)

      L=.TRUE.
      V=ABS(X)
      IF(V .GE. 8.0) GO TO 5
    3 F=0.0625*X**2-2.0
      B =           + 0.00000 00000 000114
      A = F * B     - 0.00000 00000 005777
      B = F * A - B + 0.00000 00000 252812
      A = F * B - A - 0.00000 00009 424213
      B = F * A - B + 0.00000 00294 970701
      A = F * B - A - 0.00000 07617 587805
      B = F * A - B + 0.00001 58870 192399
      A = F * B - A - 0.00026 04443 893486
      B = F * A - B + 0.00324 02701 826839
      A = F * B - A - 0.02917 55248 061542
      B = F * A - B + 0.17770 91172 397283
      A = F * B - A - 0.66144 39341 345433
      B = F * A - B + 1.28799 40988 576776
      A = F * B - A - 1.19180 11605 412169
      A = F * A - B + 1.29671 75412 105298
      BESJY=0.0625*(A-B)*X
      IF(L) RETURN

      B =           - 0.00000 00000 000244
      A = F * B     + 0.00000 00000 012114
      B = F * A - B - 0.00000 00000 517212
      A = F * B - A + 0.00000 00018 754703
      B = F * A - B - 0.00000 00568 844004
      A = F * B - A + 0.00000 14166 243645
      B = F * A - B - 0.00002 83046 401495
      A = F * B - A + 0.00044 04786 298671
      B = F * A - B - 0.00513 16411 610611
      A = F * B - A + 0.04231 91803 533369
      B = F * A - B - 0.22662 49915 567549
      A = F * B - A + 0.67561 57807 721877
      B = F * A - B - 0.76729 63628 866459
      A = F * B - A - 0.12869 73843 813500
      A = F * A - B + 0.04060 82117 718685
      BESJY=0.636619772367581*LOG(X)*BESJY-0.636619772367581/X
     1     +0.0625*(A-B)*X
      RETURN

    5 F=256.0/X**2-2.0
      B =           - 0.00000 00000 000007
      A = F * B     + 0.00000 00000 000055
      B = F * A - B - 0.00000 00000 000468
      A = F * B - A + 0.00000 00000 004699
      B = F * A - B - 0.00000 00000 057049
      A = F * B - A + 0.00000 00000 881690
      B = F * A - B - 0.00000 00018 718907
      A = F * B - A + 0.00000 00617 763396
      B = F * A - B - 0.00000 39872 843005
      A = F * B - A + 0.00089 89898 330859
      A = F * A - B + 2.00180 60817 200274
      P=A-B
      B =           + 0.00000 00000 000007
      A = F * B     - 0.00000 00000 000046
      B = F * A - B + 0.00000 00000 000360
      A = F * B - A - 0.00000 00000 003264
      B = F * A - B + 0.00000 00000 035152
      A = F * B - A - 0.00000 00000 468636
      B = F * A - B + 0.00000 00008 229193
      A = F * B - A - 0.00000 00209 597814
      B = F * A - B + 0.00000 09138 615258
      A = F * B - A - 0.00009 62772 354916
      A = F * A - B + 0.09355 55741 390707
      Q=8.0*(A-B)/V
      F=V-2.356194490192345
      A=COS(F)
      B=SIN(F)
      F=0.398942280401432/SQRT(V)
      IF(L) GO TO 7
      BESJY=F*(Q*A+P*B)
      RETURN
    7 BESJY=F*(P*A-Q*B)
      IF(X .LT. 0.0) BESJY=-BESJY
      RETURN

      ENTRY BESY0(X)

      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 4
      GO TO 8

      ENTRY BESY1(X)

      IF(X .LE. 0.0) GO TO 9
      L=.FALSE.
      V=X
      IF(V .GE. 8.0) GO TO 5
      GO TO 3

    9 BESJY=0.
      PRINT 100,X
      RETURN
  100 FORMAT(1X,36HBESJY ... NON-POSITIVE ARGUMENT X = ,E15.4)

      END
C THIS IS INTWKT.FOR

CCPY                         6TH copy of ADZINT to be used for 6 integrals
C List of GLOBAL Symbols

C     FUNCTION   ADZ6NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ6PL (F, I, IER)
C     SUBROUTINE ADZ6AL (F,I)
C     SUBROUTINE SGL6NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT6LZ
C     FUNCTION   INT6SZ ()
C
C     COMMON / ADZ6RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
C     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
C     > FA, FB, NUMINT, ICTA, ICTB, IB

C                        ****************************

      FUNCTION ADZ6NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 2000)

	double precision ERS_OLD
C
C                   Work space:
      COMMON / ADZ6RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ6RK /

      IER = 0
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ6NT call',
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ6NT call',
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB

      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZ6AL(F,I)
   10     CONTINUE
       CALL TOT6LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
 	  ERS_OLD = ERS
          IF (ERS .GT. TARGET)  THEN
              OLDINT = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*OLDINT .GT. TARGET) CALL ADZ6PL(F,I,IER)
   40         CONTINUE
CZLdebug
              IF (IER .EQ. 0.AND.ERS_OLD.NE.ERS)  GOTO 30
          ENDIF
      ADZ6NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ6PL (F, I, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ6RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ6RK /
      DATA TINY / 1.D-20 /

      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)

      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)

      CALL ADZ6AL (F, I)
      CALL ADZ6AL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES)
C                                                             Update new global
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR
      ELSE
         FAC = 1.
      ENDIF

      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC

      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ6AL (F,I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E25)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ6RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ6RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.

      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL6NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL6NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF

      RESULT(I) = TEM
      ERR   (I) = ABS (ER)

      RETURN
C                        ****************************
      END

      SUBROUTINE SGL6NT (IACT, F1, F2, F3, DX, FINT, ESTER)

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E30 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3) .LT. T1**2/HUGE) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END

      FUNCTION INT6SZ ()
C                    Return number of intervals used in last call to ADZ6NT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ6RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ6RK /
      INT6SZ = NUMINT
      RETURN
C                        ****************************
      END
C
      SUBROUTINE TOT6LZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ6RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ6RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCPY                         5TH copy of ADZINT to be used for 5 integrals
C List of GLOBAL Symbols

C     FUNCTION   ADZ5NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ5PL (F, I, IER)
C     SUBROUTINE ADZ5AL (F,I)
C     SUBROUTINE SGL5NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT5LZ
C     FUNCTION   INT5SZ ()
C
C     COMMON / ADZ5RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
C     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
C     > FA, FB, NUMINT, ICTA, ICTB, IB

C                        ****************************

      FUNCTION ADZ5NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 5000)
C
C                   Work space:
      COMMON / ADZ5RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ5RK /

	double precision ERS_OLD

      IER = 0
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ5NT call',
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ5NT call',
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB

      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZ5AL(F,I)
   10     CONTINUE
       CALL TOT5LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
	  ERS_OLD=ERS
          IF (ERS .GT. TARGET)  THEN
              OLDINT = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*OLDINT .GT. TARGET) CALL ADZ5PL(F,I,IER)
   40         CONTINUE
CZLdebug
              IF (IER .EQ. 0.AND.ERS_OLD.NE.ERS)  GOTO 30
          ENDIF
      ADZ5NT = RES
      ERREST = ERS

      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ5PL (F, I, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 5000)
      COMMON / ADZ5RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ5RK /
      DATA TINY / 1.D-20 /

      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)

      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)

      CALL ADZ5AL (F, I)
      CALL ADZ5AL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES)
C                                                             Update new global
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR
      ELSE
         FAC = 1.
      ENDIF

      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC

      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ5AL (F,I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E25)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 5000)
      COMMON / ADZ5RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ5RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.

      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL5NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL5NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF

      RESULT(I) = TEM
      ERR   (I) = ABS (ER)

      RETURN
C                        ****************************
      END

      SUBROUTINE SGL5NT (IACT, F1, F2, F3, DX, FINT, ESTER)

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E30 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3) .LT. T1**2/HUGE) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END

      FUNCTION INT5SZ ()
C                    Return number of intervals used in last call to ADZ5NT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 5000)
      COMMON / ADZ5RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ5RK /
      INT5SZ = NUMINT
      RETURN
C                        ****************************
      END
C
      SUBROUTINE TOT5LZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 5000)
      COMMON / ADZ5RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ5RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCPY                         4TH copy of ADZINT to be used for 4 integrals
C List of GLOBAL Symbols

C     FUNCTION   ADZ4NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ4PL (F, I, IER)
C     SUBROUTINE ADZ4AL (F,I)
C     SUBROUTINE SGL4NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT4LZ
C     FUNCTION   INT4SZ ()
C
C     COMMON / ADZ4RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
C     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
C     > FA, FB, NUMINT, ICTA, ICTB, IB

C                        ****************************

      FUNCTION ADZ4NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 2000)
C
C                   Work space:
      COMMON / ADZ4RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ4RK /

CZL
      double precision ERS_OLD

      IER = 0
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ4NT call',
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ4NT call',
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB

      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZ4AL(F,I)
   10     CONTINUE
       CALL TOT4LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          ERS_OLD=ERS
          IF (ERS .GT. TARGET)  THEN
              OLDINT = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*OLDINT .GT. TARGET) CALL ADZ4PL(F,I,IER)
   40         CONTINUE
CZLdebug
              IF (IER .EQ. 0.AND.ERS_OLD.NE.ERS)  GOTO 30
          ENDIF
      ADZ4NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ4PL (F, I, IER)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ4RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ4RK /
      DATA TINY / 1.D-20 /

      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)

      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)

      CALL ADZ4AL (F, I)
      CALL ADZ4AL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES)
C                                                             Update new global
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR
      ELSE
         FAC = 1.
      ENDIF

      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC

      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ4AL (F,I)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E25)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ4RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ4RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.

      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL4NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL4NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF

      RESULT(I) = TEM
      ERR   (I) = ABS (ER)

      RETURN
C                        ****************************
      END

      SUBROUTINE SGL4NT (IACT, F1, F2, F3, DX, FINT, ESTER)

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E30 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3) .LT. T1**2/HUGE) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -D1)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END

      FUNCTION INT4SZ ()
C                    Return number of intervals used in last call to ADZ4NT
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ4RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ4RK /
      INT4SZ = NUMINT
      RETURN
C                        ****************************
      END
C
      SUBROUTINE TOT4LZ
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 2000)
      COMMON / ADZ4RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES,
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT),
     > FA, FB, NUMINT, ICTA, ICTB, IB

      SAVE / ADZ4RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

CDECK  ID>, INPTN
      SUBROUTINE INPTN

C     DUMMY SUBROUTINE TO SET UP NECESSAY COMMON BLOCKS FOR FINI - USER DEFINED

      RETURN
C                        ****************************
      END

cpn=======================================================================
c
c             Finite parts of modified Bessel functions 
c              (used in the heavy-quark resummation)
c                 FOR  IHQMASS=1
c 
cpn=======================================================================
      FUNCTION fc0(x) 
cpn Returns the finite part of the modified Bessel function K0(x)
      DOUBLE PRECISION fc0,x,tmpC 
      double precision euler

      DOUBLE PRECISION pp2,pp3,pp4,pp5,pp6,pp7
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1, 
     >                 q2,q3,q4,q5,q6,q7,y

      DATA pp2,pp3,pp4,pp5,pp6,pp7
     > /3.5156229d0,3.0899424d0,1.2067492d0,  
     > 0.2659732d0,0.360768d-1,0.45813d-2/ 
      
      DATA p2,p3,p4,p5,p6,p7
     >  /0.42278420d0,0.23069756d0, 
     >   0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/ 
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,
     >  -0.7832358d-1,0.2189568d-1, -0.1062446d-1,0.587872d-2,
     >  -0.251540d-2,0.53208d-3/ 

cpn                    Euler constant
      data euler /0.577215664901532860606512/

      SAVE pp2,pp3,pp4,pp5,pp6,pp7, 
     >  p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,
     >  euler

      if (x.le.2.0) then  
        y=(x/3.75)**2 
        tmpC=y*(pp2+y*(pp3+y*(pp4+y*(pp5+y*(pp6+y*pp7))))) 

        y=x*x/4.0
        fc0=-log(x/2.0)*tmpC +(y*(p2+y*(p3+ 
     >    y*(p4+y*(p5+y*(p6+y*p7))))))   
      else 
        y=(2.0/x) 
        fc0=(exp(-x)/sqrt(x))*(q1+y
     >  *(q2+y*(q3+ y*(q4+y*(q5+y*(q6+y*q7)))))) + log(x/2.0) + euler 
      endif 
      return 
      end !fc0

cpn --------------------------------------------------------------------
      FUNCTION fc1(x) 
c Returns the finite part of the modified Bessel function x*K1(x) 
      DOUBLE PRECISION fc1,x, tmpC 

      DATA pp1,pp2,pp3,pp4,pp5,pp6,pp7
     > /0.5d0, 0.87890594d0, 0.51498869d0,
     >  0.15084934d0, 0.2658733d-1, 0.301532d-2, 0.32411d-3/ 

      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1, 
     >                 q2,q3,q4,q5,q6,q7,y
c Accumulate polynomials in double precision. 
      DATA p1,p2,p3,p4,p5,p6,p7
     > /1.0d0, 0.15443144d0, -0.67278579d0,
     >  -0.18156897d0, -0.1919402d-1, -0.110404d-2, -0.4686d-4/ 
      DATA q1,q2,q3,q4,q5,q6,q7
     > /1.25331414d0, 0.23498619d0, -0.3655620d-1,
     >  0.1504268d-1, -0.780353d-2, 0.325614d-2, -0.68245d-3/

      SAVE pp1,pp2,pp3,pp4,pp5,pp6,pp7,
     >  p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 

      if (x.le.2.0) then 
c Polynomial fit. 
        y=(x/3.75)**2 
        tmpC=pp1+y*(pp2+y*(pp3+y*(pp4+y*(pp5+y*(pp6+y*pp7))))) 

        y=x*x/4.0
        fc1= x*x*log(x/2.0)*tmpC +(p1+y*(p2+y*(p3+ 
     >    y*(p4+y*(p5+y*(p6+y*p7)))))) -1.0d0 
      else 
        y=2.0/x 
        fc1= x * exp(-x)/sqrt(x) *(q1+y
     >  *(q2+y*(q3+ y*(q4+y*(q5+y*(q6+y*q7)))))) -1.0d0 
      endif 
      return 
      end !fc1

cpn-----------------------------------------------------------------------
      FUNCTION bessk0(x) 
C USES bessi0
c Returns the modified Bessel function K0(x) for positive double precision x. 
      DOUBLE PRECISION bessk0,x 
      DOUBLE PRECISION bessi0 
      external bessi0
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1, 
     >                 q2,q3,q4,q5,q6,q7,y
c Accumulate polynomials in double precision. 
      DATA p1,p2,p3,p4,p5,p6,p7
     >  /-0.57721566d0,0.42278420d0,0.23069756d0, 
     >   0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/ 
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,
     >  -0.7832358d-1,0.2189568d-1, -0.1062446d-1,0.587872d-2,
     >  -0.251540d-2,0.53208d-3/ 

      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7 

      if (x.le.2.0) then 
c Polynomial fit. 
        y=x*x/4.0
        bessk0=(-log(x/2.0)*bessi0(x)) +( p1+y*(p2+y*(p3+ 
     >    y*(p4+y*(p5+y*(p6+y*p7)))))) 
      else 
        y=(2.0/x) 
        bessk0=(exp(-x)/sqrt(x))*(q1+y
     >  *(q2+y*(q3+ y*(q4+y*(q5+y*(q6+y*q7)))))) 
      endif 
      return 
      end !bessk0
cpn-----------------------------------------------------------------------
cpn===================================================================
cpn Modified Bessel functions I_0(x), K_0(x), I_1(x), and K_1(x)
cpn                    Copied from the Numerical recipies in Fortran
      function  bessi0(x) 
c Returns the modified Bessel function I0(x) for any double precision x. 
      DOUBLE PRECISION bessi0,x 
      DOUBLE PRECISION ax 
      DOUBLE PRECISION pp1,pp2,pp3,pp4,pp5,pp6,pp7,
     >  qq1,qq2,qq3,qq4,qq5,qq6,qq7,qq8,qq9,y 
c Accumulate polynomials in double precision. 

      DATA pp1,pp2,pp3,pp4,pp5,pp6,pp7
     > /1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,  
     > 0.2659732d0,0.360768d-1,0.45813d-2/ 
      DATA qq1,qq2,qq3,qq4,qq5,qq6,qq7,qq8,qq9
     >  /0.39894228d0,0.1328592d-1, 0.225319d-2,-0.157565d-2,
     >   0.916281d-2,-0.2057706d-1, 0.2635537d-1,-0.1647633d-1,
     >   0.392377d-2/

      SAVE pp1,pp2,pp3,pp4,pp5,pp6,pp7,
     >  qq1,qq2,qq3,qq4,qq5,qq6,qq7,qq8,qq9 

      if (abs(x).lt.3.75) then 
        y=(x/3.75)**2 
        bessi0=pp1+y*(pp2+y*(pp3+y*(pp4+y*(pp5+y*(pp6+y*pp7))))) 
      else 
        ax=abs(x) 
        y=3.75/ax 
        bessi0=(exp(ax)/sqrt(ax))*(qq1+y*(qq2+y*(qq3+y*(qq4 
     >    +y*(qq5+y*(qq6+y*(qq7+y*(qq8+y*qq9 )))))))) 
      endif 
      return 
      end !bessi0


cpn-----------------------------------------------------------------------
