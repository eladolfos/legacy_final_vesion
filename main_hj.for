CCPY Dec 2005: Modify the output format for ResBos-A version 
C    (W+_RA_UDB, Z0_RA_UUB, etc)
C MArch 7, 2005: Add heavy quark amss effect in the CSS piece for 
c (b+bbar --> HB) process. The keyword is iHQMass.
C Error: March 1, 2005: Correct the error of C(1) for H+ and HB;
c    SET I_CXFF=1 FOR NLO CALCULATION 
C Feb 10, 2005, B(2) for bb-->HB process is added, and eliminate 
C "omega" in HB process which is not needed if the input
C B quark mass is the MS-bar running mass at M_B scale. 
C Error: April 14, 2004 (correct ALEPASY for A0, etc)
C Feb 2004, Add (INONPERT.EQ.5) 
C Change the input format of TYPE_V (also inside  common.for)
C I have checked that this fortran version agrees with the 
C Legacy++2.2.3 version with the input lines in legacy.in file as follows:
C Input file for C++ version of Legacy 
C 1,0,0                               > Fast flags for CxFCxF, CxF, Sud (0/1)
C Namely, insider this fortran version, the ifast flag is set for CXFCXF only.
C ==> subroutine SetCf
C
C
C-----------------------------------------------------------------------
C%%%%%%%%%%%%%%%%%%%%            MAIN.FOR          %%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%         VERSION 4.4.3        %%%%%%%%%%%%%%%%%%
C-----------------------------------------------------------------------
C
CCPY The last time this code is modified: Jan. 31, 1996
CSM  The last time this code is modified: Feb.  7, 1996
CsB  The last time this code is modified: May   2, 1998
CPN  The last time this code is modified: Jan. 15, 1999
C
C ----------------------------------------------------------------------
      PROGRAM MAIN
C ----------------------------------------------------------------------
C Open output file
cZL check the argument of program, if yes use it as the output file name, otherwise use "legacy.out" as output file name
cZL unit No. of output file is 22 (following the original setup)
	character*100 jobname
	integer length
	common /jobcommon/ jobname,length

      if(iargc().eq.0)then
        jobname="legacy"
      else
        call getarg(1,jobname)
      endif
	call trmstr(jobname,length)
      open(unit=22,file=jobname(1:length)//'.out',
     &  form='formatted',status="unknown")
c      OPEN(UNIT=22,FILE='legacy.out',FORM='FORMATTED'
c     >)
!    >,    Action='Write')
!    >,    STATUS='New' )
      CALL LEGACY
      CLOSE(22)
      WRITE(*,*) ' The Program has completed.  Long live the program!'
      STOP
      END

C ----------------------------------------------------------------------
      SUBROUTINE LEGACY
C ----------------------------------------------------------------------
C
C Purpose and Methods:  Compute resummed portion of
C                                  d(sigma)/d(P_T)/dy/dQ^2
C                       for a ... vector boson with a given
C                       rapidity, mass, transverse momentum, and
C                       nonperturbative parameters g1,g2,g3,Q0.
C                       This is performed by the call to fresum().
C                       In the case of on-shell W and Z bosons a delta
C                       function over Q^2 takes the result to
C                       d(sigma)/d(P_T)/dy.
C     N.B.  This just gives the "CSS" piece.  To this one
C           must add the "perturbative" contribution and
C           subtract the "asymptotic" contribution; however,
C           at low transverse momenta, say,
C                   P_T (of W or Z) < 30 GeV at the Tevatron,
C           the "CSS" piece can provide a decent approximation.
C
C Inputs  :  rap_input  -> desired rapidity of boson
C            pt_input   -> desired transverse momentum of boson
C            q0_input   -> value for nonpert. parameter Q0
C            fit_par(i) -> value for nonpert. parameter gi (i=1,2,3)
C
C            A file called legacy.in
C                   -------------------------------------------------
C                                 Format for legacy.in
C                   -------------------------------------------------
C                    p+p or p+pbar COLLIDER          p+N Collisions
C                   -------------------------------------------------
C                   IBEAM                           IBEAM, Fract_N
C                   ECM,LTO                         ECM,LTO
C                   MT,MW,MZ,MH,MA                  MT,MW,MZ,MH,MA
C                   TYPE_V                          TYPE_V
C                   ISET,iPionPDF                   ISET,iPionPDF
C                   INONPERT,IFLAG_C3               INONPERT,IFLAG_C3
C                   N_SUD_A,N_SUD_B,N_WIL_C,I_FSR   N_SUD_A,N_SUD_B,N_WIL_C,I_FSR
C                   PDF_FILE                        PDF_FILE
C                   MU,MD,MS,MC,MB                  MU,MD,MS,MC,MB
C                   BMAX                            BMAX
C                   Q_GRID.INP                      Q_GRID.INP
C                   QT_GRID.INP                     QT_GRID.INP
C                   Y_GRID.INP                      Y_GRID.INP
C                   iProc                           iProc
C                   IqTMn,IqTMx,IqTSt,IyMn,IyMx,IySt,IQMn,IQMx,IQSt
C                   -------------------------------------------------
C            IBEAM   -> 1(pp collider),-1(p pbar collider), 0(p+N collisions)
C                       -2(pion_minus+N scattering)
C            FRACT_N -> (IBEAM=0 only) fraction of neutrons in nucleus
C            ECM     -> center of mass energy for nucleon-nucleon collision
C            LEPASY  -> 0 for parity-conserving CSS contribution
C                       1 for parity-violating CSS contribution
C            LTO     ->  0 CSS, perturbative and asymptotic parts
C                        1 DELTA_SIGMA FROM QT=0 TO PT
C                       -1 LEADING ORDER
C                        2 ASYMPTOTIC PART without calculating CSS
C                        3 Y piece
C                        4 CSS pieces separating out alpha_s delta(1-z) in CFns
C            MT,MW,MZ,MH,MA    -> Masses(top,W,Z,Higgs,virtual Photon)
C            TYPE_V  -> 2 characters for boson choice (either W+,W-,H0,A0,Z0)
C  Iset assignments as of Aug 19, 1999:

C  Iset      Description       FortranFun needed     Data/Input File(s) needed
C ------------------------------------------------------------------------ 
C     0       Test                  --                    --
C  -----------------------------------------------------------------------
C             Evolved from
C  Isetev0,1    input fn       Fini.f (DummyArg)        (Block Data DatPdf)
C  (10, 11)
 
C             Evolved from
C  Isetin0   input para tbl       (Evlini.f)            pdf.ini   
C  (900)      (use default grid pts)
C  Isetin1   input para tbl       (Evlini.f)            xxx.ini   
C  (901)      (prompt for grid pts)
C  Isetin2   input para tbl       (Evlini.f)            pdf.ini   
C  (902)     (set grid pts by ParPdf or use defaults in Evlpac commons)

C             Read in table
C  Isettbl  fr evolved results    (Tblini.f)            xxx.tbl
C  (911)
C  -------------------------------------------------------------------------
C  
C  1101 - 5    CTEQ1 M,MS,ML,D,L       (Ctq1Pd.f)
C  1201 - 6    CTEQ2 M,MS,MF,ML,L,D    (Ctq2Pd.f)
C  1301 - 3    CTEQ3 M,L,D             (Ctq3Pd.f)
C  1401 - 3    CTEQ4 M,D,L   \
C  1404 - 9    CTEQ4 A1 - A5  |        (Ctq4Pdf.f)         cteq4xx.tbl
C  1409 -10    CTEQ4 HJ,LQ    |
C  1411 -12    CTEQ4 HQ,HQ1   |
C  1413 -14    CTEQ4 F3,F4   /
C  
C  2001 - 4    KMRS, MRSS0,D0,D-      Strc78 - 81.f       For078 - 81.Dat
C  2005 - 7      MRSS0',D0',D-'       Strc82 - 84.f       For097 - 99.Dat
C  2008 - 10      MRSA, A', J         Strc33,30,37.f      For033,30,37.Dat
C  2011 - 14       MRSR1 - 4            mrsr.f            mrsr1 - 4.Dat
C  2021 - 25     MRS98 1 - 5        Mrs98.f; Mrs98x.f     Ftxx.dat
C                             (all MRS's need PDMRS.f and MRSEB.f)
C
C  3001 - 3    GRV LO,NLO,DIS ('94')  PDGRV.f + GRV94LO,HO,DI.f  
C  
C  - - - - - - - - - - - - - - Some old PDF's - - - - - - - - - - - - - -
C    801         EHLQ1                  PSHK.f
C    802         Duke-Owens 1           PSHK.f
C    803         DFLM                   DFLM.f
C    804         MT90-S1                Pdxmt.f
C

C            iPionPDF -> PDF choice for pions
C                       1 = ABFKW (Phys. Lett. 233B (1989) 517)
C                       2 = GRV92 appears (Z Phys C53 (1992) 651)
C            INONPERT,IFLAG_C3 -> just set to 3,1
C                                 IFLAG_C3=1: canonical choice of C1,C2,C3,C4
C            N_SUD_A,N_SUD_B,N_WIL_C,(I_FSR) -> Order upto calc. Sudakov and C fns
C                                        I_FSR for final state radiation
C            PDF_FILE           -> if ISET>=900 (i.e., input of initial PDF
C                                 for evolution is requested), then this
C                                 is where you put the input file name
C                                 (e.g., le26.ini)
C            MU,MD,MS,MC,MB -> quark masses: up,down,strange,charm,bottom
C            BMAX -> b-space cutoff indicating nonperturbative regime(1/GeV)
C
C Outputs :  A file called output.dat is created which tabulates the
C            g1,g2,g3,Pt,and "CSS" contribution in those respective columns;
C            output is also dumped into units 2 and 13.
C Controls:  ????
C
C Program History:
C  Created:  Sometime 1990  C.-P. Yuan
C          4 November 1993  G.A.Ladinsky and C.-P.Yuan
C                            using QCD library of John Collins,
C                            Porter Johnson, Sijin Qian, and Wu-Ki Tung.
C          16 June 1994     G.A.Ladinsky
C                            updated for latest libraries of Wu-Ki Tung''s;
C                            major frontline modifications for flexibility.
C          9 December 1994  C.-P. Yuan
C                            updated for POLARIZED W and Z productions;
C                            On-shell scheme is implemented for electroweak
C                            parameters in STANDARD subroutine;
C                            Two erros in function APDF (for getting PDF)
C                            are corrected.
C          11 October 1995  C. Balazs and C.-P. Yuan
C                            confirm NO_SIGMA0=0 case for TYPE_V='A0','W+',
C                            'W-' and 'Z0', and for LEPASY=0 or 1.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NIN,NOUT,NWRT
      COMMON/SETUP1/ DEBUG,NIN,NOUT,NWRT
      LOGICAL DEBUG
C      DEBUG=.TRUE.
C a dummy routine at the start of the program, before anything happens
      CALL BIRTH
C set some unit #''s and read the input file
      CALL WHATIS
C setup the electroweak parameters
      CALL STANDARD
cpn      CALL STANDARD_old
C setup the PDF
      CALL SETTHEPDF
C the remaining setup requirements, including more QCD physics
      CALL FSTAGE
C compute the cross sections
      CALL XSECT
C close some input files (e.g., PDF inputs); perform any other clean-up
      CALL CLEANUP
C a dummy routine after everything has happened
      CALL DEATH
C
      RETURN
      END
C --------------------------------------------------------------------------
      SUBROUTINE BIRTH
C --------------------------------------------------------------------------
C a routine executed at the start of the program
      RETURN
      END
C --------------------------------------------------------------------------
      SUBROUTINE DEATH
C --------------------------------------------------------------------------
      RETURN
      END
C --------------------------------------------------------------------------
      SUBROUTINE WHATIS
C --------------------------------------------------------------------------
C Reads in the information from the file wp_res.in
C Details were given in the subroutine legacy.
C
      IMPLICIT NONE
      INCLUDE 'common.for'
      CHARACTER*78 TMPCHR
      REAL*8 TMPVAR(8)
      REAL*8 VMAS
      COMMON/BOSONMASS/ VMAS
      INTEGER IQTMN, IQTMX, IQTST, IYMN, IYMX, IYST,
     &        IQMN, IQMX, IQST
      COMMON/IMMS/ IQTMN, IQTMX, IQTST, IYMN, IYMX, IYST,
     &             IQMN, IQMX, IQST
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME
C Legacy input
cc      INTEGER I_PROC
cc      COMMON / PARTPROC / I_PROC
      CHARACTER*40 QGFN, QTGFN, YGFN, DUMMY
      REAL*8 QG(200),PT(200),Y(200)
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q,N_QT,N_Y, LTO
      COMMON / NGRID / N_Q,N_QT,N_Y, LTO
      REAL*8 GEES(3),Q0IN
      Common /NonPertC/ GEES,Q0IN
      Integer iPionPDF
      Common / PionPDF / iPionPDF
      Integer ngag
      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC
      Integer KinCorr
      Common / CorrectKinematics / KinCorr


      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

CCPY
      CHARACTER*40 RA_TYPE_V,RA_PDF_FILE
      INTEGER LEN_TYPE_V,LEN_PDF_FILE
      DATA LEN_TYPE_V,LEN_PDF_FILE /2,10/

CJI Jan 2015: Add in hj calculation
      REAL*8 D1s, R, t, ptj
      Common /HJ/ D1s, R, t, ptj

cZL
	character*100 jobname
	integer length
	common /jobcommon/ jobname,length

C UNIT: intended for simplicity in unit conversion, now needs work.
C       must be set equal to one until repaired
      UNIT=1.D0
C I/O UNITS
      NIN=23
      NOUT=2
      NWRT=13
C      nin=5
C      nout=6
C      nwrt=6
C
!      OPEN(UNIT=NIN,FILE='legacy.in',STATUS='old')
      OPEN(UNIT=NIN,FILE=jobname(1:length)//'.in',STATUS='old')
C
C            IBEAM -> 1(pp collider),-1(p pbar collider), 0(p+N collisions)
C                     -2(pion_minus-nucleus scattering)
C            KinCorr -> 0: x_{1,2} =   Q/Sqrt[S] e^{+/-y} in Asy and CSS 
C                       1: x_{1,2} = M_T/Sqrt[S] e^{+/-y} in Asy and CSS
C 
 
      READ(NIN,*) IBEAM, FRACT_N, KinCorr
C FRACT_N=0.6d0 is used for a Cu nucleus (cf. E288 paper)
C FRACT_N=0.5556d0 is used for E706 (Andre Maul)
      IF(IBEAM.EQ.0 .OR. IBEAM.EQ.-2) THEN
        IF(FRACT_N.GT.1.D0.OR.FRACT_N.LT.0.D0) THEN
          PRINT *,
     >      ' ERROR: Forget FRACT_N in input? (IBEAM,FRACT_N)> ',
     >      IBEAM,FRACT_N
          CALL QUIT
        ENDIF
      ELSE
        FRACT_N=0.0D0
      ENDIF
C C.M. ENERGY, ECM,  IN GeV
C
C FOR PARITY-CONSERVING PART: LEPASY=0
C FOR PARITY-VIOLATING PART: LEPASY=1
C
C SET NO_SIGMA0=0 FOR INCLUDING CONST=\sigma_0 IN THE RATES
C SET NO_SIGMA0=1 FOR NOT INCLUDING CONST=\sigma_0 IN THE RATES
      NO_SIGMA0=0
C
CCPY SET LTOPT=0 FOR CALCULATING THE ASYMPTOTIC PART
C    SET LTOPT=1 FOR CALCULATING THE DELTA_SIGMA FROM QT=0 TO PT
C    SET LTOPT=-1 FOR CALCULATING THE LEADING ORDER RESULT
C    Set LTOpt=2 for CALCULATING THE ASYMPTOTIC PART without calculating
C                the resummed one
C
      READ(NIN,*) ECM,LTO, iFast
      ECM=ECM*UNIT
      IF(ECM.LT.1.D0*UNIT) THEN
        WRITE(*,*)
     >     'ERROR(FRESUM): Forget to remove FRACT_N in input? (ECM)> ',
     >     ECM,FRACT_N
        CALL QUIT
      ENDIF

CCPY SET I_CXFF=1 FOR NLO CALCULATION (using LTOPT to control)
C HARD-WIRED FLAG
      I_CXFF=0
      IF(LTO.EQ.1) I_CXFF=1
      WRITE(*,*)' CHECK LTOPT, I_CXFF =',LTO, I_CXFF 

CCPY MARCH 2005
CJI Sept 2014: Force only KINCORR = 1 to not use ifast
      IF(KINCORR.EQ.1) THEN
        WRITE(*,*) ' *** KINCORR=1 ONLY WORKS WITH GLEEN SPEED-UP 
     >VERSION, IN ZJSUDNPT *** '
        IFAST=0
C      ELSE
C        WRITE(*,*) ' *** KINCORR=0 USE PAVEL IFAST     
C     >VERSION, IN ZJSUDNPT *** '
C        IFAST=1
      ENDIF       
      
CCPY Feb 2004
C      IF(LEPASY.GT.1) THEN
C        WRITE(*,*) ' ERROR: Forget LEPASY in input? (ECM,LEPASY)> ',
C     >    ECM,LEPASY
C        CALL QUIT
C      ENDIF
C CHECK INPUT
C      IF(IBEAM.NE.-1 .AND. LEPASY.EQ.1) THEN
C        WRITE(*,*) ' LEPASY=1, only for p-pbar collision '
C        CALL QUIT
C      ENDIF
C INPUT MASS IN GeV
C MT IS FOR TOP, MW FOR W BOSON, MZ FOR Z BOSON;
C MH IS HIGGS MASS, MA IS VIRTUAL PHOTON MASS.
      READ(NIN,*) MT,MW,MZ,MH,MA
      MT=MT*UNIT
      MW=MW*UNIT
      MZ=MZ*UNIT
      MH=MH*UNIT
      MA=MA*UNIT
C there is a second common block holding quark masses
      QMAS_T=MT
C______Type of produced particle(s)
C      TYPE_V='W+','W-','Z0','A0','H0','AA','AG','ZZ','ZG','HP','HM','HZ',
C             'H+','HB','WW_UUB','WW_DDB',"ZU","ZD"
C Add in for ResBos-A
C______TYPE_V = 'W+_RA_UDB','W-_RA_DUB','Z0_RA_UUB','Z0_RA_DDB'
CJI January 2015: Added in 'hj' for higgs + 1 Jet
CCPY      READ(NIN,'(A2)') TYPE_V
      
      READ(NIN,'(A)') RA_TYPE_V
      Call TrmStr(RA_TYPE_V, LEN_TYPE_V)      
      TYPE_V=RA_TYPE_V(1:LEN_TYPE_V) 
      
      print*,' RA_TYPE_V, LEN_TYPE_V ,  TYPE_V '
      print*,RA_TYPE_V,LEN_TYPE_V,TYPE_V 
        
      
C CHECK INPUT
CCPY      IF( (TYPE_V.NE.'W+' .AND. TYPE_V.NE.'W-' .AND. TYPE_V.NE.'Z0')
C     >  .AND. (LEPASY.EQ.1) )THEN
C        WRITE(*,*) ' LEPASY=1, only for TYPE_V=W+,W-,or Z0 '
C        CALL QUIT
C      ENDIF

C CONVERT TO UPPER CASE
      CALL UC(TYPE_V)
C CHECK ON INPUTS
        IF(IBEAM.EQ.-2 .AND. (TYPE_V.NE.'AA' .AND. TYPE_V.NE.'AG')) THEN
          PRINT*,' NO SUCH PROCESS FOR IBEAM = ',IBEAM
          CALL QUIT
        ENDIF

C PDF choice for proton and pion, and parameters for H+ production:
C i_RunMass = 0 ->    pole top-mass is used in the b h+ t coupling
C i_RunMass = 1 -> running top-mass is used in the b h+ t coupling
C I_MODEL.EQ.1 FOR TOPCOLOR, AND 2 FOR 2HDM
      READ(NIN,*) ISET,iPionPDF, i_RunMass,i_Model
      IF(TYPE_V.EQ.'HB' .AND. I_MODEL.NE.2) THEN
        WRITE(*,*) ' THIS MODEL IS NOT YET IMPLEMENTED '
        CALL QUIT
      ENDIF

C Choose nonperturbative parametrization with NONPERT (1,2 or 3);
C NONPERT IS DEFINED IN SUBROUTINES WSETUP AND HSETUP
C      INONPERT=3
C select the choice of C1,C2,C3 and C4 with IFLAG_C3: 1 for CANONICAL
C choice, 2 for ELSE
C       IFLAG_C3=1
      READ(NIN,*) INONPERT,IFLAG_C3
CCPY Feb 2004
       IF(LTO.EQ.1 .AND. IFLAG_C3.NE.1) THEN
         PRINT*,' THIS ONLY WORKS FOR IFLAG_C3 = 1 '
         CALL QUIT
       ENDIF
       IF(TYPE_V.EQ.'H0' .AND. IFLAG_C3.NE.1) THEN
         PRINT*,' THIS ONLY WORKS FOR IFLAG_C3 = 1 '
         CALL QUIT
       ENDIF
       
CCPY SELECT THE ORDER IN ALPHA_STRONG FOR CALCULATING A AND B
C FUNCTION IN SUDKOV FACTOR (N_SUD_AB = 1 OR 2); AND WILSON COEFFICIENT
C FUNCTIOPN C (N_WIL_C = 0 OR 1).
C      N_SUD_AB=2
C      N_WIL_C=1
CCPY NOV. 1995, S_Q IS NOW DEFINED TO BE THE HARD SCALE FOR
C SUDAKOV FACTOR
C ONLY FOR TOP QUARK PAIR PRODUCTION (I.E. FOR 'GL' or 'GG' PROCESS):
      READ(NIN,*) N_SUD_A,N_SUD_B,N_WIL_C,I_FSR
      IF(TYPE_V.NE.'GL' .AND. TYPE_V.NE.'GG') THEN
        I_FSR=0
      ENDIF
C if evolving from an initial PDF, indicate the input file
C      READ(NIN,'(A10)') PDF_FILE
      READ(NIN,'(A)') RA_PDF_FILE
      Call TrmStr(RA_PDF_FILE,LEN_PDF_FILE)
      PDF_FILE=RA_PDF_FILE(1:LEN_PDF_FILE)
      
      PDF_EVL_NAME=PDF_FILE
C quark masses other than for the top quark
      READ(NIN,*) QMAS_U,QMAS_D,QMAS_S,QMAS_C,XMBOT
      QMAS_U=QMAS_U*UNIT
      QMAS_D=QMAS_D*UNIT
      QMAS_S=QMAS_S*UNIT
      QMAS_C=QMAS_C*UNIT
      QMAS_B=XMBOT*UNIT
C consistency check on quark masses if using WKT evolution files *.ini
CCPY      IF(ISET.GE.900) THEN
      IF(ISET.EQ.902) THEN
        OPEN(UNIT=12,FILE=PDF_FILE,STATUS='old')
        READ(12,'(a78)') TMPCHR
        READ(12,'(a78)') TMPCHR
        READ(12,*) TMPVAR(1),TMPVAR(2),TMPVAR(3),TMPVAR(4),TMPVAR(5),
     >             TMPVAR(6),TMPVAR(7),TMPVAR(8)
        CLOSE(12)
        IF(TMPVAR(7).NE.QMAS_C.OR.TMPVAR(8).NE.QMAS_B) THEN
          WRITE(*,*)
     >            'WARNING: Quark masses inconsistent with ',PDF_FILE
CCPY          WRITE(*,*) '   Input Charm:  ',QMAS_C,' reset to ',TMPVAR(7)
C          WRITE(*,*) '   Input Bottom: ',QMAS_B,' reset to ',TMPVAR(8)
C          QMAS_C=TMPVAR(7)
C          QMAS_B=TMPVAR(8)

        WRITE(*,*) 
     >    'In the following, using qmas_c from the legacy.in file'
        WRITE(*,*) 'c masses from the .IN and .INI files are', 
     >    qmas_c, TMPVAR(7)
        WRITE(*,*) 
     >    'In the following, using qmas_b from the legacy.in file'
        WRITE(*,*) 'b masses from the .IN and .INI files are', 
     >    qmas_b, TMPVAR(8)

	ENDIF

CCPY Feb 2009
CCPY!      ELSEIF(ISET.EQ.903) THEN
CCPY!        OPEN(UNIT=12,FILE=PDF_FILE,STATUS='old')
CCPY!        READ(12,'(a78)') TMPCHR
CCPY!        READ(12,'(a78)') TMPCHR
CCPY!        READ(12,*) TMPVAR(1),TMPVAR(2),TMPVAR(3),TMPVAR(4),TMPVAR(5),
CCPY!     >             TMPVAR(6),TMPVAR(7),TMPVAR(8)
CCPY!        CLOSE(12)
CCPY!        IF(TMPVAR(7).NE.QMAS_C.OR.TMPVAR(8).NE.QMAS_B) THEN
CCPY!          WRITE(*,*)
CCPY!     >            'WARNING: Quark masses inconsistent with ',PDF_FILE
CCPY!CCPY          WRITE(*,*) '   Input Charm:  ',QMAS_C,' reset to ',TMPVAR(7)
CCPY!C          WRITE(*,*) '   Input Bottom: ',QMAS_B,' reset to ',TMPVAR(8)
CCPY!C          QMAS_C=TMPVAR(7)
CCPY!C          QMAS_B=TMPVAR(8)
CCPY!
CCPY!        WRITE(*,*) 
CCPY!     >    'In the following, using qmas_c from the legacy.in file'
CCPY!        WRITE(*,*) 'c masses from the .IN and .INI files are', 
CCPY!     >    qmas_c, TMPVAR(7)
CCPY!        WRITE(*,*) 
CCPY!     >    'In the following, using qmas_b from the legacy.in file'
CCPY!        WRITE(*,*) 'b masses from the .IN and .INI files are', 
CCPY!     >    qmas_b, TMPVAR(8)
CCPY!
CCPY!        ENDIF
CCPY!
      ELSE
CCPY
        WRITE(*,*) '   Use PDF Table file !!!!!!!!! '     
      ENDIF
C bmax, the coordinate space cutoff for entering the nonperturbative region
C      BMAX=1.0D0/2.0D0/UNIT
      READ(NIN,*) BMAX
      BMAX=BMAX/UNIT

C Grid Filenames
      READ(NIN,*) I_PROC
CsB   I_Proc = 1 Calculate only QI QJ -->  G V process
C     I_Proc = 2 Calculate only  G QI --> QJ V process

      READ(NIN, '(A40)') QGFN
      READ(NIN, '(A40)') QTGFN
      READ(NIN, '(A40)') YGFN
      IF (QGFN.NE.'-') THEN
        PRINT*, ' Reading  Q values from ', QGFN
        OPEN(UNIT=31,FILE=QGFN,STATUS='OLD')
        DO II = 1,200
          READ(31,*,END=98) QG(II)
c          Print*, '>>>', ii,  QG(ii)
        END DO
   98   N_Q = II - 1
        CLOSE(31)
      END IF
      IF (QTGFN.NE.'-') THEN
        PRINT*, ' Reading qT values from ', QTGFN
        OPEN(UNIT=31,FILE=QTGFN,STATUS='OLD')
        DO II = 1,200
          READ(31,*,END=99) PT(II)
C          Print*, '>>>', ii, pt(ii)
        END DO
   99   N_QT = II - 1
        CLOSE(31)
      END IF
      IF (YGFN.NE.'-') THEN
        PRINT*, ' Reading  y values from ', YGFN
        OPEN(UNIT=31,FILE=YGFN,STATUS='OLD')
        DO II = 1,200
          READ(31,*,END=100) Y(II)
C          Print*, '>>>', ii, y(ii)
        END DO
  100   N_Y = II - 1
        CLOSE(31)
      END IF

CsB * Input g s & Q0
      Read(nIn, *) Gees(1), Gees(2), Gees(3), Q0In, ngag

C FOR G G INITIATED PROCESSES A1 = CA (= 3   for QCD)
C FOR Q Q INITIATED PROCESSES A1 = CF (= 4/3 for QCD)
C For iNONPERT contribution
 
CCPY      IF(TYPE_V.EQ.'H0'.OR.TYPE_V.EQ.'GG'.OR. TYPE_V.EQ.'AG') THEN
      IF(TYPE_V.EQ.'H0'.OR.TYPE_V.EQ.'GG'.OR. 
     >TYPE_V.EQ.'AG'.OR.TYPE_V.EQ.'ZG'.or.TYPE_V.EQ.'HJ') THEN
        If (ngag.Eq.2) then
          Gees(2) = 9.d0/4.d0 * Gees(2)
        Else If (ngag.Eq.3) then
          Gees(1) = 9.d0/4.d0 * Gees(1)
          Gees(2) = 9.d0/4.d0 * Gees(2)
          Gees(3) = 9.d0/4.d0 * Gees(3)
        End If
      END IF

C IqTMn, IqTMx, IqTSt, IyMn, IyMx, IySt, IQMn, IQMx, IQSt =
C   1  ,  24  ,   1  ,  1  ,  13 ,  1  ,   1 ,  13 ,  1
C Middle (peak) point at
C   8  ,   8  ,   1  ,  7  ,  7  ,  1  ,   5 ,  5  ,  1
      READ(NIN, *) IQTMN, IQTMX, IQTST, IYMN, IYMX, IYST,
     &             IQMN, IQMX, IQST

CJI Jan 2015: Added in hj and need to read jet raidus
      READ(NIN, *) DUMMY
      READ(NIN, *) DUMMY
      if(Type_v.eq.'HJ') then
          read(nin,*) R
      endif

      CLOSE(NIN)
C
C===============================
CCPY HARD-WIRED HERE FOR HEAVY QUARK MASS EFFECT: 
C SET IHQMASS=1 FOR TURNING ON HEAVY QUARK MASS EFFECT
C This is for 'HB' production, from b+bbar fusion.
      IHQMASS=0
      IF(TYPE_V.EQ.'HB') THEN
        IHQMASS=1      
        IHQPDF=5
ccpy: To avoid oscillating CSS output in large QT, we have to use 
C     the consistent value for the mass of bottom quark used in 
C     deriving the PDF, which is TMPVAR(8), not qmas_b.
C     Note that qmas_b is used to calculate the Yukawa coupling.
C             
        XMHQ=TMPVAR(8)
        WRITE(*,*) ' IHQMASS,IHQPDF,XMHQ,qmas_b =',
     >  IHQMASS,IHQPDF,XMHQ,qmas_b
        
        IF(IFAST.EQ.1) THEN
          WRITE(*,*) ' IFAST=1 IS NOT WORKING FOR THIS OPTION '
          CALL EXIT
        ENDIF  
      ENDIF
      
C===============================
CCPY Initialization
      I_RESET=0

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE FSTAGE
C --------------------------------------------------------------------------
C The resummation routines themselves have some individual setup requirements;
C  there is also some remaining QCD setup to perform.
C  This is accomplished in this routine.
      IMPLICIT NONE
      INCLUDE 'common.for'
C especially, if the evolve PDF isnot used, we must set LAMBDA, etc.
CCPY NORDER IS NOT REALLY USED
      NORDER=2
      CALL QCD(NORDER)
C
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE XSECT
C --------------------------------------------------------------------------
C See, subroutine legacy.
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 VMAS, PTSTRACH
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP, NF_EFF
CCPY 
      Integer KinCorr
      Common / CorrectKinematics / KinCorr


      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME
      INTEGER I_GEES
      REAL*8 RDUMP
      LOGICAL FIRST_QT
CJI Jan 2015: Add in hj calculation
      REAL*8 D1s, R, t, ptj
      Common /HJ/ D1s, R, t, ptj

CJI October 2013: Variables added to match the header... What are they used for in the code???
      INTEGER iResScheme, fCxFCxF, fCxF, fSud, i_pma, iFract_n

C      TESTING = .TRUE.
      TESTING = .FALSE.

      If (N_Q_GD.LT.IQMAX) then
        Print*, ' Resetting IQMAX to N_Q_GD = ', N_Q_GD
        IQMAX = N_Q_GD
      End If
      If (N_qT_GD.LT.IpTMAX) then
        Print*, ' Resetting IpTMAX to N_qT_GD = ', N_qT_GD
        IpTMAX = N_qT_GD
      End If
      If (N_y_GD.LT.IyMAX) then
        Print*, ' Resetting IyMAX to N_y_GD = ', N_y_GD
        IyMAX = N_y_GD
      End If

CCPY      open(unit=22,file='output.dat',status='new')
CsB      OPEN(UNIT=22,FILE='LEGACY_W.OUT', STATUS='unknown')
      TEMP1(1,1)=-9999.9999D9

      Print*, ' Type_V : ', Type_V
      If (Type_V.Eq.'W+' .or. Type_V.Eq.'W-'
     >  .or. Type_V.Eq.'W+_RA_UDB' .or. Type_V.Eq.'W-_RA_DUB'
CJI September 2013 added in new names for W boson processes
     >  .or. Type_V.Eq.'WU' .or. Type_V.Eq.'WD' ) then
        VMAS = MW
        Print*, ' MW = ', MW
      Else If (Type_V.Eq.'Z0' .or. Type_V.Eq.'Z0_RA_UUB'
     >  .or. Type_V.Eq.'Z0_RA_DDB' 
CJI September 2013 added in new names for Z boson processes
     >  .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD' ) then
        VMAS = MZ
        Print*, ' MZ = ', MZ
      Else If (Type_V.Eq.'ZZ' .or. Type_V.Eq.'ZG') then
        VMAS = 2.d0*MZ
      Else If (Type_V.Eq.'WW_UUB' .or. Type_V.Eq.'WW_DDB') then
        VMAS = 2.d0*MW
      Else If (Type_V.Eq.'A0' .or. Type_V.Eq.'AA' .or.
     .    Type_V.Eq.'AG') then
CsB        VMAS = 0.d0
C        Print*, ' QG(1) = ', QG(1)
        VMas = QG(1)
      Else If (Type_V.Eq.'H0' .or. Type_V.Eq.'H+' .or.
     .         Type_V.Eq.'HB' .or. Type_V.EQ.'HJ' ) then
        VMAS = MH
      Else If (Type_V.Eq.'HP' .or. Type_V.Eq.'HM') then
        VMAS = MH + MW
      Else If (Type_V.Eq.'HZ') then
        VMAS = MH + MZ
      Else If (Type_V.Eq.'GL' .or. Type_V.Eq.'GG') then
        VMAS = 2.d0*MT
      Else
        Print*, ' Set VMas in XSect for boson: ', Type_V
        Stop
      End If

      NF_EFF=NFL(VMAS)

C SETUP THE DEFAULT VALUE FOR IHADRON
C IHADRON=1 IS FOR PP, PPBAR, PB COLLISIONS.
C THIS IS NEEDED FOR CONV...
      IHADRON=1

CCPY      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR. TYPE_V.EQ.'GG') THEN
      IF(TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'AG' .OR. 
     >TYPE_V.EQ.'GG'.OR. TYPE_V.EQ.'ZG') THEN
        CALL HSETUP(NF_EFF)
      ELSE IF(TYPE_V.EQ.'HJ') THEN
        CALL HJSETUP(NF_EFF)
      ELSE
        CALL WSETUP(NF_EFF)
      ENDIF

CJI Variables to be used in the header... Purpose in the grid generation???
      iResScheme=1
      fCxFCxF=0
      fCxF=0
      fSud=0
      i_pma=1
      iFract_n=INT(Fract_n)

C Write of the header of the output file
CCPY      WRITE(22,*) ' ECM,IBEAM,LEPASY,LTO,iProc,Fract_n,iFast'
CJI       WRITE(22,*) ' ECM,IBEAM,idummy,LTO,iProc,Fract_n,iFast'
CJI       WRITE(22,102) ECM,IBEAM,LEPASY,LTO,i_Proc,Fract_n,iFast
CJI October 2013: Modified to match the format of the c++ version
      WRITE(22,*) 'ecm,  ibeam, lto,iproc,iResScheme,fract_n,',
     > 'fCxFCxF,fCxF,fSud, i_pma'
      WRITE(22,106) ECM,IBEAM,LTO,i_Proc,iResScheme,iFract_n,fCxFCxF,
     >  fCxF,fSud,i_pma
CCPY      IF(ISET.GE.900) THEN
      IF(ISET.EQ.902) THEN
CCPY
        WRITE(22,*) ' ISET,.INI,INONPERT,IFLAG_C3,C1,C2,C3,C4,KinCorr '
        WRITE(22,*) ISET,'  ',PDF_EVL_NAME,INONPERT,IFLAG_C3,   
     >       C1,C2
        write(22,*) C3, C4,KinCorr
      ELSE
        WRITE(22,'(A43)') 'ISET,INONPERT,IFLAG_C3,C1,C2,C3,C4,KinCorr '
CCPY        WRITE(22,*) ISET,INONPERT,IFLAG_C3,C1,C2,C3,C4,KinCorr
        WRITE(22,110) ISET,INONPERT,IFLAG_C3,C1,C2
        WRITE(22,111) C3,C4,KinCorr
      ENDIF
      WRITE(22,'(A14)') 'MT,MW,MZ,MH,MA'
      PRINT*, MT, MW, MZ
      WRITE(22,107) MT,MW,MZ
      WRITE(22,107) MH,MA
      WRITE(22,'(A9,3X,A40)') ' TYPE_V:', TYPE_V
CJI      WRITE(22,101) '  g1,g2,g3,Q0:', GEES(1), GEES(2), GEES(3), Q0IN
CJI October 2013: Modified to match the format of the c++ version
      WRITE(22,105) 'bmax,g1,g2,g3,Q0in: ',BMAX,GEES(1),GEES(2),
     >     GEES(3),Q0IN
CJI      WRITE(22,*)'VMAS,B0,A1,A2,A3,B1,B2,NF_EFF,NORDER,N_SUD_A,N_SUD_B,',
CJI     >   'N_WIL_C, I_FSR'
CJI October 2013: Modified to match the format of the c++ version
      WRITE(22,*)'Mass of V,h1,A1,A2,A3,B1,B2,Nf,norder,n_sud_a,',
     >   'n_sud_b,n_wil_c,i_fsr,R'
      WRITE(22,107) VMAS,B0,A1
      WRITE(22,108) A2,A3,B1,B2
      WRITE(22,109) NF_EFF,NORDER,N_SUD_A,N_SUD_B,N_WIL_C, I_FSR, R

  101 FORMAT(1X,A14,2X,4(G10.4,2X))
  102 FORMAT(1X,G10.4,2X,4I3,2X,G10.4,2X,I3)
 103  FORMAT (1X,I3,2X,A15,2(2X,I3),2X,(2X,G10.4))
 104  FORMAT (1X, 2(G10.4,2X))
CJI FORMATS Added to get proper formatting to match c++ version header 
  105 FORMAT(A20,2x,5(G10.4,2X))
  106 FORMAT(F6.1,9(3X,I3))
  107 FORMAT(1X,G2.4,2(2X,G3.4))
  108 FORMAT(1X,4G1.4)
  109 FORMAT(1X,6(I1,9X),F3.1)
  110 FORMAT(1X,I5,2X,I1,2X,I1,3X,G2.4,3X,G3.2)
  111 FORMAT(G10.4,3X,G10.4,3X,I1)
 
CCPY To note this is generated from fortran version of Legacy code
CJI Removed posting of where the grid file is generated to match current header used
CJI      WRITE(22,*) 'This is generated from fortran version of 
CJI     > Legacy code'
C
C      If (Testing) Call PDFTest
C
CJI      WRITE(22,*) '  '
C
C      If(.False.) Call TotalLOXSect
C

C Initialize iCF1
      ICF1 = 0

      IF(TYPE_V.EQ.'W+'.OR.TYPE_V.EQ.'W-'.OR.TYPE_V.EQ.'Z0'.OR.
     >TYPE_V.EQ.'HP'.OR.TYPE_V.EQ.'HM'.OR.TYPE_V.EQ.'HZ'.Or.
     >Type_V.Eq.'A0'.or.TYPE_V.Eq.'H+'.or.TYPE_V.Eq.'HB') THEN
CCPY HP FOR H W^+; HM FOR H W^-; HZ FOR H Z  ASSOCIATED PRODUCTION
        Call VecBosXSec
        Return
CJI September 2013: Modified the Type_V to match those from the c++ version
      Else If (Type_V.Eq.'W+_RA_UDB' .or. Type_V.Eq.'W-_RA_DUB' 
     > .or. Type_V.Eq.'Z0_RA_UUB' .or. Type_V.Eq.'Z0_RA_DDB'
     > .or. Type_V.Eq.'WU' .or. Type_V.Eq.'WD'
     > .or. Type_V.Eq.'ZU' .or. Type_V.Eq.'ZD' ) then
         Call VecBosXSec
         Return
      Else If (TYPE_V.Eq.'AA' .or. TYPE_V.Eq.'ZZ') then
        Call AAXSect
        Return
      Else If (TYPE_V.Eq.'WW_UUB' .or. TYPE_V.Eq.'WW_DDB') then
        Call WWXSect
        Return
      Else If (TYPE_V.Eq.'AG' .or. TYPE_V.Eq.'ZG') then
        Call AGXSect
        Return
      Else If (TYPE_V.EQ.'H0' .OR. TYPE_V.EQ.'GG') THEN
        Call HXSect
        Return
      Else If (TYPE_V.EQ.'HJ') THEN
        Call HJXSect
        Return
      End If

      END

C --------------------------------------------------------------------------
      SUBROUTINE AAXSect
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 VMAS
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      INTEGER I_GEES
      REAL*8 RDUMP
      LOGICAL FIRST_QT
      Real*8 pT_Stretch

      pT_Stretch = 1.d0
CCPY Sep 2006: switch off the hard-wired pT_Stretch value
C      If (ECM.GT.5.d3) pT_Stretch = 1.2d0
C      If (TYPE_V.Eq.'AA') pT_Stretch = 1.d0
C      If (TYPE_V.Eq.'ZZ') pT_Stretch = 1.2d0

      IF (LTO.EQ.-1) THEN ! LO piece
C ----------------------------------------------
C LO piece: L0 and zero.
C ----------------------------------------------
        LTOPT = -1
        WRITE(22,*) ' Q,qT,y, Pert(LO), Zero'
        PRINT *,    ' Q,qT,y, Pert(LO), Zero'
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            FIRST_QT=.TRUE.
            DO 470 I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(FIRST_QT) THEN
                VPT=0.D0
                FIRST_QT=.FALSE.
                GOTO 530
              ELSE
                GOTO 470
              ENDIF
  530         CONTINUE
C Get the LO L0 term from Pert.for
              LepAsy = 0
CCPY April 14, 2004
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,0) = Asy
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                    Temp1(0,0),0.d0
              PRINT 102, QG(K),VPT,RAPIN,
     &                    Temp1(0,0),0.d0
  470       CONTINUE   ! DO QT
          ENDDO        ! DO Y
        ENDDO          ! DO Q
        Return

      Else IF (LTO.EQ.0) THEN
        LTOPT = 0
C -------------------------------------------------------------------
C Resummed L0 pieces: qqB_delta(1-z), alpha_s^0+alpha_s(qqB_(1-z)+qG)
C -------------------------------------------------------------------
        WRITE(22,*) ' Q,qT,y, CSS, CSS "qqB_delta(1-z)", Asymtotic '
        PRINT *,    ' Q,qT,y, CSS, CSS "qqB_delta(1-z)", Asymtotic '

        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          iqin=k
CsB_______Set up the grid files with the convolutions of C-functions
          if (iFast.eq.1) call SetCf(no_asym)

          lepasy=0
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            iyin=j
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch

CsB Resummed L0 piece for qqB ~ alpha_s delta(1-z)
C   For testing purposes: sum of options iCF1 = 1 and 2 gives the
C   result of iCF1 = 0.
CsB iCF1 separates the alpha_s delta(1-Z) and LO + alpha_s (1-Z)+gluon
C   contributions.
C   Settings: iCF1 = 1 alpha_s delta(1-Z) qqB contribution returned
C             iCF1 = 2 LO + alpha_s (1-Z) qqB + Gq contributions
C
C              iCF1 = 1
C              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
C     >                  IER_RES)
C              TEMP1(1,1)=FRESUM

CsB Resummed L0 pieces for qqB and Gq ~ delta(1-z) + alpha_s (...) (1-z)
              ICF1 = 0
CCPY SEPT 2009
              I_RESET=1
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,ICF1) = FRESUM

cdump
c      print*,'ICF1, TEMP1(1,ICF1) =',ICF1, TEMP1(1,ICF1) 


CsB Resummed LO + alpha_s (1-Z) qqB + Gq contributions
              ICF1 = 2
              I_RESET=1
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,ICF1) = FRESUM

cdump
c      print*,'ICF1, TEMP1(1,ICF1) =',ICF1, TEMP1(1,ICF1) 


CsB Perturbative and asymptotic L0 pieces
C              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
C     >                        ASY,IER_PERT)
C              TEMP1(2,1) = Asy
              temp1(2,1) = 0
C Output the results
C             -> Q, Q_T, y, CSS, CSS qqB_delta(1-z), Asymtotic:
              PRINT 102,    QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(1,0)-TEMP1(1,2),TEMP1(2,1)
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(1,0)-TEMP1(1,2),TEMP1(2,1)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.1) THEN ! NLO_Sig piece
C -----------------------------------------------
C QT < QT_SEP perturbative piece.
C -----------------------------------------------
        LTOPT = 1
        WRITE(22,*) ' Q,qT,y, NLO_Sig (L0,A3) '
        PRINT *,    ' Q,qT,y, NLO_Sig (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,3,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Singular L0, A3:
              PRINT 122,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,122) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.2) THEN ! Asymptotic piece
        LTO = 3
C ---------------------------------
C Asymptotic (singular, NLO) piece.
C ---------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Asymptotic (L0 and zero) '
        PRINT *, ' Q,qT,y, Asymptotic (L0 and zero) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0
C                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y,  Aymptotic L0, A3:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),0.d0 !TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),0.d0 !TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.3) THEN
C -------------------------------------------
C Y piece: perturbative and asymptotic parts.
C -------------------------------------------
C This fakes the Y piece for the di-photon production.
C It calculates the DY Y piece with the quark charges^4 (instead of ^2).
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Asymptotic(L0,A3), Pert.(L0,A3,A1,A2,A4) '
        PRINT *, ' Q,qT,y, Asymptotic (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0
C                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert

ccpy sept 2009: 
C FOR 'AA' AND 'ZZ' PROCESSES, THE EXACT PERT PART IS INCLUDED 
C BY A SEPERATE CALCULAITON USING NAPRT=21,22,23 (IN RESBOS).
C HENCE, PERT PART SHOULD NOT NOE INCLUDED IN THE Y-GRID FILE, IN 
C ORDER TO BE CONSISTENT WITH RESBOS WHICH CALCULATES W-ASYM AND 
C PERT PART IN DIFFERENT RUNS.
              IF(TYPE_V.EQ.'AA'.OR.TYPE_V.EQ.'ZZ') THEN                
                DO I_PERT = 0,0
                 TEMP1(3,I_PERT+1) = 0.D0
                END DO ! I_Pert
              ENDIF


C Output the results
C             -> Q, Q_T, y:
              PRINT 103,    QG(K),VPT,RAPIN
              WRITE(22,103) QG(K),VPT,RAPIN
C             -> Asymptotic L0, A3 and perturbative L0, A3:
              PRINT 107,    TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
              WRITE(22,107) TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
C             -> Perturbative A1, A2, A42:
              PRINT 105,    TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
              WRITE(22,105) TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)

            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN


      ENDIF

  100 FORMAT(1X,5(G10.4,2X),G14.7)
  101 FORMAT(1X,A14,2X,4(G10.4,2X))
  102 FORMAT(1X,3(f9.4,2X),3(G14.8,2X))
  103 FORMAT(1X,3(G10.4,2X))
  104 FORMAT(1X,2(G14.8,2X))
  105 FORMAT(1X,3(G18.12,2X))
  106 FORMAT(2X,3(G10.4,2X),2(G16.10,2X))
  107 FORMAT(1X,4(G18.12,2X))
  112 FORMAT(1X,3(f9.4,2X),9(G14.8,2X))
  122 FORMAT(1X,G11.5,2X,2(f9.4,2X),3(G14.8,2X))
CsB <-
  300 FORMAT(A3,8X,A8,6X,A4,12X,A4,13X,A9)
  400 FORMAT(2(G10.4,2X),3(G14.7,2X))
C
  600 CONTINUE
      CLOSE(22)

      Return
      END

C --------------------------------------------------------------------------
      SUBROUTINE WWXSect
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 VMAS
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      INTEGER I_GEES
      REAL*8 RDUMP
      LOGICAL FIRST_QT
      Real*8 pT_Stretch

      pT_Stretch = 1.d0
CCPY Sep 2006: switch off the hard-wired pT_Stretch value

      IF (LTO.EQ.-1) THEN ! LO piece
C ----------------------------------------------
C LO piece: L0 and A3.
C ----------------------------------------------
        LTOPT = -1
        WRITE(22,*) ' Q,qT,y, Pert(LO, A3)'
        PRINT *,    ' Q,qT,y, Pert(LO, A3)'
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            FIRST_QT=.TRUE.
            DO 470 I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(FIRST_QT) THEN
                VPT=0.D0
                FIRST_QT=.FALSE.
                GOTO 530
              ELSE
                GOTO 470
              ENDIF
  530         CONTINUE
C Get the LO L0 term from Pert.for
              LepAsy = 0
CCPY April 14, 2004
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,0) = Asy

C Get the LO A3 term from Pert.for
              LepAsy = 1
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,1) = Asy

              WRITE(22,102) QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
              PRINT 102, QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
  470       CONTINUE   ! DO QT
          ENDDO        ! DO Y
        ENDDO          ! DO Q
        Return

      Else IF (LTO.EQ.0) THEN
        LTOPT = 0
C -------------------------------------------------------------------
C Resummed L0, A3 pieces: qqB_delta(1-z), alpha_s^0+alpha_s(qqB_(1-z)+qG)
C -------------------------------------------------------------------
        WRITE(22,*) ' Q,qT,y, CSS (L0,A3), CSS qqB_delta(1-z) (L0,A3),'
        PRINT *,    ' Q,qT,y, CSS (L0,A3), CSS qqB_delta(1-z) (L0,A3),'

        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          iqin=k
CsB_______Set up the grid files with the convolutions of C-functions
          if (iFast.eq.1) call SetCf(no_asym)

          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            iyin=j
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch

CsB Resummed L0 piece for qqB ~ alpha_s delta(1-z)
C   For testing purposes: sum of options iCF1 = 1 and 2 gives the
C   result of iCF1 = 0.
CsB iCF1 separates the alpha_s delta(1-Z) and LO + alpha_s (1-Z)+gluon
C   contributions.
C   Settings: iCF1 = 1 alpha_s delta(1-Z) qqB contribution returned
C             iCF1 = 2 LO + alpha_s (1-Z) qqB + Gq contributions

CsB Resummed L0 pieces for qqB and Gq ~ delta(1-z) + alpha_s (...) (1-z)
              ICF1 = 0
              I_RESET=1
              lepasy=0
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,ICF1) = FRESUM

              ICF1 = 0
              I_RESET=1
              lepasy=1
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(2,ICF1) = FRESUM

CsB Resummed LO + alpha_s (1-Z) qqB + Gq contributions
              ICF1 = 2
              I_RESET=1
              lepasy=0
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,ICF1) = FRESUM

              ICF1 = 2
              I_RESET=1
              lepasy=1
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(2,ICF1) = FRESUM

c              temp1(2,1) = 0
C Output the results
C             -> Q, Q_T, y, CSS, CSS qqB_delta(1-z), Asymtotic:
              PRINT 102,    QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(2,0),
     &                      TEMP1(1,0)-TEMP1(1,2),TEMP1(2,0)-TEMP1(2,2)
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(2,0),
     &                      TEMP1(1,0)-TEMP1(1,2),TEMP1(2,0)-TEMP1(2,2)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.1) THEN ! NLO_Sig piece
C -----------------------------------------------
C QT < QT_SEP perturbative piece.
C -----------------------------------------------
        LTOPT = 1
        WRITE(22,*) ' Q,qT,y, NLO_Sig (L0,A3) '
        PRINT *,    ' Q,qT,y, NLO_Sig (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              I_PERT = 0
              DO LEPASY = 0,1 ! Ex symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,LEPASY+1) = ASY
c                TEMP1(3,LEPASY+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Singular L0, A3:
              PRINT 122,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,2)
              WRITE(22,122) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,2)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.2) THEN ! Asymptotic piece
        LTO = 3
C ---------------------------------
C Asymptotic (singular, NLO) piece.
C ---------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Asymptotic (L0 and A3) '
        PRINT *, ' Q,qT,y, Asymptotic (L0 and A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              I_PERT = 0
              DO LEPASY = 0,1 ! Ex symmetric piece
c                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,LEPASY+1) = ASY
c                TEMP1(3,LEPASY+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y,  Aymptotic L0, A3:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,2)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,2)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.3) THEN
C -------------------------------------------
C Y piece: perturbative and asymptotic parts.
C -------------------------------------------
C This fakes the Y piece for the di-photon production.
C It calculates the DY Y piece with the quark charges^4 (instead of ^2).
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Asymptotic(L0,A3), Pert.(L0,A3,A1,A2,A4) '
        PRINT *, ' Q,qT,y, Asymptotic (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,3,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,LEPASY+1) = ASY
                TEMP1(3,LEPASY+1) = PERT
              END DO ! I_Pert

ccpy sept 2009: 
C FOR 'AA' AND 'ZZ' PROCESSES, THE EXACT PERT PART IS INCLUDED 
C BY A SEPERATE CALCULAITON USING NAPRT=21,22,23 (IN RESBOS).
C HENCE, PERT PART SHOULD NOT NOE INCLUDED IN THE Y-GRID FILE, IN 
C ORDER TO BE CONSISTENT WITH RESBOS WHICH CALCULATES W-ASYM AND 
C PERT PART IN DIFFERENT RUNS.
CZL              IF(TYPE_V.EQ.'AA'.OR.TYPE_V.EQ.'ZZ') THEN
              IF(TYPE_V.EQ.'AA'.OR.TYPE_V.EQ.'ZZ'
     &          .or.TYPE_V.EQ.'WW_UUB'.or.TYPE_V.EQ.'WW_DDB') THEN   
                DO LEPASY = 0,1
                 TEMP1(3, LEPASY+1) = 0.D0
                END DO ! I_Pert
              ENDIF


C Output the results
C             -> Q, Q_T, y:
              PRINT 103,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,2),
     & TEMP1(3,1),TEMP1(3,2),0d0,0d0,0d0,0d0,0d0,0d0
              WRITE(22,103) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,2),
     & TEMP1(3,1),TEMP1(3,2),0d0,0d0,0d0,0d0,0d0,0d0
C             -> Asymptotic L0, A3 and perturbative L0, A3:
c              PRINT 107,    TEMP1(2,1),TEMP1(2,2), TEMP1(3,1),TEMP1(3,4)
c              WRITE(22,107) TEMP1(2,1),TEMP1(2,2), TEMP1(3,1),TEMP1(3,4)
C             -> Perturbative A1, A2, A42:
c              PRINT 105,    TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
c              WRITE(22,105) TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)

            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN


      ENDIF

  100 FORMAT(1X,5(G10.4,2X),G14.7)
  101 FORMAT(1X,A14,2X,4(G10.4,2X))
  102 FORMAT(1X,3(f9.4,2X),3(G14.8,2X))
  103 FORMAT(1X,3(G10.4,2X))
  104 FORMAT(1X,2(G14.8,2X))
  105 FORMAT(1X,3(G18.12,2X))
  106 FORMAT(2X,3(G10.4,2X),2(G16.10,2X))
  107 FORMAT(1X,4(G18.12,2X))
  112 FORMAT(1X,3(f9.4,2X),9(G14.8,2X))
  122 FORMAT(1X,G11.5,2X,2(f9.4,2X),3(G14.8,2X))
CsB <-
  300 FORMAT(A3,8X,A8,6X,A4,12X,A4,13X,A9)
  400 FORMAT(2(G10.4,2X),3(G14.7,2X))
C
  600 CONTINUE
      CLOSE(22)

      Return
      END

C --------------------------------------------------------------------------
      SUBROUTINE AGXSect
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 VMAS
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      INTEGER I_GEES
      REAL*8 RDUMP
      LOGICAL FIRST_QT
      Real*8 pT_Stretch

      pT_Stretch = 1.d0
CCPY Sep 2006: switch off the hard-wired pT_Stretch value
C      If (ECM.GT.5.d3) pT_Stretch = 1.2d0
C        pT_Stretch = 1.d0*(Log(ECM)/Log(1.8d3))**3 *
C     >                   Log(QG(5))/Log(80.d0)

      IF (LTO.EQ.-1) THEN ! LO piece
C ----------------------------------------------
C LO pieces: L0 and A3.
C ----------------------------------------------
        LTOPT = -1
        WRITE(22,*) ' Q,qT,y, LO L0 & A3 terms'
        PRINT *,    ' Q,qT,y, LO L0 & A3 terms'
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            FIRST_QT=.TRUE.
            DO 470 I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(FIRST_QT) THEN
                VPT=0.D0
                FIRST_QT=.FALSE.
                GOTO 530
              ELSE
                GOTO 470
              ENDIF
  530         CONTINUE
C Get the LO L0 term from Pert.for
              LepAsy = 0
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,0) = Asy
C Get the LO A3 term from Pert.for = 0
c              LepAsy = 1
c              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
c     >                       PERT,ASY,IER_PERT)
c              Temp1(0,1) = Asy
              Temp1(0,1) = 0.d0
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
              PRINT 102, QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
  470       CONTINUE   ! DO QT
          ENDDO        ! DO Y
        ENDDO          ! DO Q
        Return

      Else IF (LTO.EQ.0) THEN
        LTOPT = 0
C -----------------
C Resummed L0 piece
C -----------------
c        WRITE(22,*) ' Q,qT,y, CSS, Perturbative-Asymtotic '
c        PRINT *,    ' Q,qT,y, CSS, Perturbative-Asymtotic '

	IF(TYPE_V.EQ.'AG') THEN
          WRITE(22,*) ' Q,qT,y, CSS, CSS gg_delta(1-z), 
     &                  CSS prime, CSS prime prime, Asymtotic: '
          PRINT *,    ' Q,qT,y, CSS, CSS gg_delta(1-z), 
     &                  CSS prime, CSS prime prime, Asymtotic: '
        ELSE
          WRITE(22,*) ' Q,qT,y, CSS, CSS gg_delta(1-z), Asymtotic: '
          PRINT *,    ' Q,qT,y, CSS, CSS gg_delta(1-z), Asymtotic: '
        ENDIF

        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
CsB_______Set up the grid files with the convolutions of C-functions
          iFast=0 ! Not implemented in HCXFCXF !
          iqin=k
          if (iFast.eq.1) call SetCf(no_asym)
          lepasy=0
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
CsB Resummed L0 piece GG -> di-photon ~ delta(1-z) + alpha_s (...) (1-z)
              ICF1 = 0
              I_RESET=1
	      LEPASY=0
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,0) = FRESUM
CsB Resummed L0 piece GG -> di-photon ~ delta(1-z)
              ICF1 = 2
              I_RESET=1
	      LEPASY=0
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,2) = FRESUM
	      IF (TYPE_V.EQ.'AG') THEN
CZL: compute spin-flip prime term
                ICF1 = 0
                I_RESET=1
	        LEPASY=1
                CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                    PERT,ASY,IER_RES)
                TEMP1(1,0+LEPASY) = FRESUM
CZL: compute spin-flip prime prime term
                ICF1 = 2
                I_RESET=1
	        LEPASY=1
                CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                    PERT,ASY,IER_RES)
                TEMP1(1,2+LEPASY) = FRESUM
              ENDIF
CsB Perturbative and asymptotic L0 pieces
cpn              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
cpn     >                        ASY,IER_PERT)
cpn              TEMP1(2,0) = Pert-Asy
c              temp1(2,0) = 0d0
cpn              TEMP1(2,1) = Asy
              temp1(2,1) = 0d0
C Output the results
C             -> Q, Q_T, y,CSS, CSS gg_delta(1-z), Asymtotic:
              IF (TYPE_V.EQ.'AG') THEN
                PRINT 102,    QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(1,0)-TEMP1(1,2),
     &                      TEMP1(1,1),TEMP1(1,3),TEMP1(2,1)
                WRITE(22,102) QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(1,0)-TEMP1(1,2),
     &                      TEMP1(1,1),TEMP1(1,3),TEMP1(2,1)
              ELSE
                PRINT 102,    QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(1,0)-TEMP1(1,2),TEMP1(2,1)
                WRITE(22,102) QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(1,0)-TEMP1(1,2),TEMP1(2,1)
              ENDIF
C Output the results
C             -> Q, Q_T, y, CSS, Perturbative-Asymtotic:
c              PRINT 102,    QG(K),VPT,RAPIN,
c     &                      TEMP1(1,0),TEMP1(2,0)
c              WRITE(22,102) QG(K),VPT,RAPIN,
c     &                      TEMP1(1,0),TEMP1(2,0)

            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.1) THEN ! NLO_Sig piece
C -----------------------------------------------
C For QT < QT_SEP perturbative piece.
C -----------------------------------------------
        LTOPT = 1
        WRITE(22,*) ' Q,qT,y, NLO_Sig (L0,A3) '
        PRINT *,    ' Q,qT,y, NLO_Sig (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,3,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Singular L0, A3:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.3) THEN
C ------------------------------------------
C Y piece: perturbative and asymptotic parts
C ------------------------------------------
C This fakes the Y piece for the gluon-gluon initiated di-photon production.
C Calculates the same factor as for GG -> hG.
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        PRINT *, ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              I_PERT = 0
C Call the computation
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
              TEMP1(2,I_PERT+1) = ASY
              TEMP1(3,I_PERT+1) = PERT

CZL: trigger to compute spin-flip term for gg->\gamma\gamma subprocess
              IF (TYPE_V.EQ.'AG') THEN
                I_PERT=3
                LEPASY=1 
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              ENDIF

ccpy sept 2009: 
C FOR 'AG' (BUT NOT FOR 'ZG') PROCESS, THE EXACT PERT PART IS INCLUDED 
C BY A SEPERATE CALCULAITON USING NAPRT=21,22,23 (IN RESBOS).
C HENCE, PERT PART SHOULD NOT NOE INCLUDED IN THE Y-GRID FILE, IN 
C ORDER TO BE CONSISTENT WITH RESBOS WHICH CALCULATES W-ASYM AND 
C PERT PART IN DIFFERENT RUNS.
              IF(TYPE_V.EQ.'AG') THEN                
                DO I_PERT = 0,4
                 TEMP1(3,I_PERT+1) = 0.D0
                END DO ! I_Pert
              ENDIF


C Output the results
C             -> Q, Q_T, y:
              PRINT 103,    QG(K),VPT,RAPIN
              WRITE(22,103) QG(K),VPT,RAPIN
C             -> Singular (asymptotic) L0, A3 and perturbative L0, A3:
              PRINT 107,    TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
              WRITE(22,107) TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
C             -> Perturbative A1, A2, A42:
              PRINT 105,    TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
              WRITE(22,105) TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.2) THEN ! Asymptotic piece
        LTO = 3
C ---------------------------------
C Asymptotic (singular, NLO) piece.
C ---------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Asymptotic (L0 and zero) '
        PRINT *, ' Q,qT,y, Asymptotic (L0 and zero) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0
C                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y,  Aymptotic L0, A3:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),0.d0 !TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),0.d0 !TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ENDIF

  100 FORMAT(1X,5(G10.4,2X),G14.7)
  101 FORMAT(1X,A14,2X,4(G10.4,2X))
  102 FORMAT(1X,3(f9.4,2X),5(G14.8,2X))
  103 FORMAT(1X,3(G10.4,2X))
  104 FORMAT(1X,2(G14.8,2X))
  105 FORMAT(1X,3(G18.12,2X))
  106 FORMAT(2X,3(G10.4,2X),2(G16.10,2X))
  107 FORMAT(1X,4(G18.12,2X))
  112 FORMAT(1X,3(f9.4,2X),9(G14.8,2X))
CsB <-
  300 FORMAT(A3,8X,A8,6X,A4,12X,A4,13X,A9)
  400 FORMAT(2(G10.4,2X),3(G14.7,2X))
C
  600 CONTINUE
      CLOSE(22)

      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE QUIT
C --------------------------------------------------------------------------
      IMPLICIT NONE
CsB       Dummy routine
      PRINT *, ' Stopping in QUIT.'
      STOP
      END

C --------------------------------------------------------------------------
      SUBROUTINE PDFTEST
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IPDF,I,IER_RES
      REAL*8 VMAS,TMP,X,AMU,VPT,RAPIN,GEES(3),Q0IN
      REAL*8 FRESUM,PMA,PERT,ASY
      INCLUDE 'common.for'
      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      VMAS  = 80.D0
      VPT   = .4D0
      RAPIN = 0.D0
      GEES(1)=0.11D0
      GEES(2)=0.580D0
      GEES(3)=-1.5D0
      Q0IN=1.6D0
      I_RESET=1
      CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,IER_RES)

      IPDF = 1
      X    = .1D0
      AMU  = 5.D0
      PRINT*, ' IPDF,x,aMu,APDF:'
C      Do IPDF = -6,6
C       Do i = 2,10
C         aMu  = i*1.d0
      DO I = 1,9
        X    = .01D0 + I*.01D0
        TMP  = APDF(IPDF,X,AMU)
        PRINT*, IPDF,X,AMU,TMP
      ENDDO

      IF (.TRUE.) STOP
      END

CsB ========================================================================
      SUBROUTINE TOTALLOXSECT
CsB ========================================================================
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 TMP, DSDQDY,DSDYDQ, DSDQ,DSDY, VMAS,VRAP
      REAL*8 ADZ2NT,YMIN,QMIN,QMAX,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      COMMON / BOSONMASS / VMAS
      COMMON / BOSONRAP / VRAP
      COMMON / ERREST / ERREST
C      External AdZInt

C      Goto 100
CC Check dsdQdy & dsdydQ
C      VRap = 2.d0
C      VMas = 70.d0
C      tmp = dsdQdy (VRap)
C      Print*, '[TLOXS] VRap, dsdQdy(VRap): ', VRap, tmp
C      tmp = dsdydQ (VMas)
C      Print*, '[TLOXS] VMas, dsdQdy(VMas): ', VMas, tmp
C
CC Check dsdQ & dsdy
C      VRap = .33d0
C      tmp = dsdy (VRap)
C      Print*, '[TLOXS] VRap, dsdQ(VRap): ', VRap, tmp, '  +-', ErrEst
C      VMas = 110.d0
C      tmp = dsdQ (VMas)
C      Print*, '[TLOXS] VMas, dsdy(VMas): ', VMas, tmp, '  +-', ErrEst

  100 CONTINUE
C Calculate the total cross section
      PRINT*, '[TLOXS] VRap, dsdQdy(VRap): ', 2.D0, DSDQDY (2.D0)

      QT_V = 0.D0
      CALL YMAXIMUM
      YMIN = - YMAX
      AERR = 0.D0
      RERR = 1.D-2
      TMP = ADZ2NT (DSDY, YMIN,YMAX, AERR, RERR,
     >              ERREST, IER,IACTA,IACTB)
      PRINT*, '[TLOXS] TotLOXSec: ', TMP, '  +-', ERREST

      QMIN = 70.D0
      QMAX = 90.D0
      AERR = 0.D0
      RERR = 1.D-2
      TMP = ADZ2NT (DSDQ, QMIN,QMAX, AERR, RERR,
     >              ERREST, IER,IACTA,IACTB)
      PRINT*, '[TLOXS] TotLOXSec: ', TMP, '  +-', ERREST

      STOP
      END

CsB ************
      FUNCTION DSDQ (Q)
CsB ************
      IMPLICIT NONE
      INCLUDE 'common.for'
      REAL*8 TMP, DSDQ,DSDY,Q,Y, VMAS,VRAP, YMIN,QMIN,QMAX
      REAL*8 ADZINT,AERR,RERR,ERREST
      INTEGER IER,IACTA,IACTB
      COMMON / BOSONMASS / VMAS
      COMMON / BOSONRAP / VRAP
      COMMON / ERREST / ERREST
      EXTERNAL ADZINT, DSDQDY, DSDYDQ

      VMAS = Q
      Q_V  = Q
      QT_V = 0.D0
      CALL YMAXIMUM
      YMIN = - YMAX
      AERR = 0.D0
      RERR = 1.D-2
      TMP = ADZINT (DSDQDY, YMIN,YMAX, AERR, RERR,
     >              ERREST, IER,IACTA,IACTB)
C      Print*, '[dsdQ] VMas, AdZInt:', VMas, tmp, '  +-', ErrEst
      DSDQ = TMP
      RETURN
CsB ************
      ENTRY DSDY (Y)
C QMin & QMax are hard-wired for W here. Get them from input!!
      VRAP = Y
      QMIN = 70.D0
      QMAX = 90.D0
      AERR = 0.D0
      RERR = 1.D-2
      TMP = ADZINT (DSDYDQ, QMIN,QMAX, AERR, RERR,
     >              ERREST, IER,IACTA,IACTB)
C      Print*, '[dsdQ] VRap, AdZInt:', VRap, tmp, '  +-', ErrEst
      DSDY = TMP

  100 CONTINUE
      END

CsB ************
      FUNCTION DSDQDY (Y)
CsB ************
C The idea is that y is passed in the argument and VMas in the common block
C so that we can use AdZInt.
      IMPLICIT NONE
      REAL*8 DSDQDY,Y, DSDYDQ,Q,
     &  VPT, TMP, RAPIN,Q0IN, FGETPERT, GEES(3),
     &  FRESUM,PMA,PERT,ASY, VMAS,VRAP,
     &     MT,MW,MZ,MH,MA,MW2,MZ2,WCOUPL,ZCOUPL,HCOUPL,ACOUPL,HpCoupl,
     >     HBCoupl,FL,FR, EW_ALFA,GAMW,GAMZ,GAMW2,GAMZ2,GFermi
      INTEGER IER_PERT
      EXTERNAL FGETPERT
      COMMON / BOSONMASS / VMAS
      COMMON / BOSONRAP / VRAP
      COMMON/STAND1/ MT,MW,MZ,MH,MA,MW2,MZ2,
     >WCOUPL,ZCOUPL,HCOUPL,ACOUPL,HpCOUPL,HbCOUPL
      COMMON/STAND3/ EW_ALFA,GFermi,GAMW,GAMZ,GAMW2,GAMZ2

      RAPIN = Y
      GOTO 100
      ENTRY DSDYDQ (Q)
CsB ************
      VMAS = Q
C      Q_V  = VMas
      RAPIN = VRAP
      GOTO 100

  100 CONTINUE
      VPT   = 0.D0
      Q0IN  = 1.6D0

      TMP=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_PERT)

      FL = 1.D0
      FR = 0.D0

      ASY = ASY *
C Conversion to d(sigma)/d(Q), extra factor of 2, factor of Q^2
     *      4.D0*VMAS**4.D0 *
C V propagator: !! W mass hard-wired
     *      1.D0/((VMAS**2-MW2)**2+(MW2*GAMW)**2) *
C V-l-l couplings:
     &      (FL**2 + FR**2)*
C Angular integral over (1.d0 + Cos(the_sta)**2.d0)
     *      8.D0/3.D0

      DSDQDY = ASY
      DSDYDQ = ASY
C      Print*, VMas,VPT,RAPIN,Q0IN, dsdydQ, dsdQdy
      RETURN
      END

C --------------------------------------------------------------------------
      SUBROUTINE VecBosXSec
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 VMAS, PTSTRACH
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP, NF_EFF

      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      INTEGER I_GEES
      REAL*8 RDUMP
      Real*8 L0A3(2)
      Real*8 pT_Stretch
      LOGICAL FIRST_QT
cc      INTEGER I_PROC
cc      COMMON / PARTPROC / I_PROC

CsB   Choosing different VPT grids for different ECM by scaling
C      the VPT grids at 1.8 TeV
C      pT_Stretch = 1.d0*(Log(ECM)/Log(1.8d3))**3*
C     >                 Log(QGrid(5))/Log(80.d0)
      pT_Stretch = 1.0d0
CCPY      If (ECM.GT.5.d3) pT_Stretch = 1.2d0

CCPY         IF(TYPE_V.EQ.'GL') THEN
C STRETCH QT GRID GIVEN BY QTGRID
C              pT_Stretch = 1.5D0
C              IF(VMAS.GT.2.0*MT+100.D0)pTStrach = 2.0D0
C              IF(VMAS.GT.2.0*MT+500.D0)pTStrach = 2.5D0
C            ELSE
C      pT_Stretch = 1.0D0
C            ENDIF

      IF (LTO.EQ.-1) THEN ! LO piece
C ----------------------------------------------
C LO pieces: L0 and A3.
C ----------------------------------------------
        LTOPT = -1
        WRITE(22,*) ' Q,qT,y, LO L0 & A3 terms'
        PRINT *,    ' Q,qT,y, LO L0 & A3 terms'
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            FIRST_QT=.TRUE.
            DO 470 I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(FIRST_QT) THEN
                VPT=0.D0
                FIRST_QT=.FALSE.
                GOTO 530
              ELSE
                GOTO 470
              ENDIF
  530         CONTINUE
C Get the LO L0 term from Pert.for
              LepAsy = 0
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,0) = Asy
C Get the LO A3 term from Pert.for
              LepAsy = 1
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,1) = Asy
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
              PRINT 102, QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
  470       CONTINUE   ! DO QT
          ENDDO        ! DO Y
        ENDDO          ! DO Q
        Return
      ELSE IF (LTO.EQ.0) THEN
C ----------------------------------------------
C Resummed pieces: L0 and A3 and pert-asympt.
C ----------------------------------------------
        LTOPT = 0
CCPY Dec 2005: Modify the output format for ResBos-A version 
C    (W+_RA_UDB, Z0_RA_UUB, etc)
c        WRITE(22,*) ' Q,qT,y, CSS L0 & A3, Pert-Sing L0 '
c        PRINT *,    ' Q,qT,y, CSS L0 & A3, Pert-Sing L0 '

        WRITE(22,*) ' Q,qT,y, CSS L0, CSS A3, ASY L0, ASY A3 '


        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
cpn_______Calculate convolutions of C-functions
           iqin=k
           if (iFast.eq.1) call SetCf(yes_asym)
            

          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            iyin=j
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch

ccpy checking:
c	vmas=81.2504514 
c	vpt=0.986546309
c	rapin=-0.288222327

CsB Resummed L0 (=1+cos^2) piece
              LEPASY = 0
              I_RESET=1
C              CALL SRESUMCsB(VMAS,VPT,RAPIN,GEES,Q0IN,L0A3,PMA,PERT,ASY,
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
C              TEMP1(1,1) = L0A3(1)
C              TEMP1(1,3) = L0A3(2)
              TEMP1(1,1)=FRESUM
              If (Type_V.Eq.'HP' .or. Type_V.Eq.'HM' .or.
     .            Type_V.Eq.'HZ') Goto 360
CsB Resummed A3 (=2*cos) piece
              LEPASY = 1
              I_RESET=1
              CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,ASY,
     >                  IER_RES)
              TEMP1(1,3)=FRESUM
  360         Continue
CCPY Dec 2005: Modify the output format for ResBos-A version 
C    (W+_RA_UDB, Z0_RA_UUB, etc)
CsB Perturbative and asymptotic L0 pieces
              LEPASY = 0
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                        ASY,IER_PERT)
              TEMP1(2,1)=ASY
              If (Type_V.Eq.'HP' .or. Type_V.Eq.'HM' .or.
     .            Type_V.Eq.'HZ') Goto 362
              LEPASY = 1
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                        ASY,IER_PERT)
              TEMP1(2,3)=ASY
  362         Continue

C Output the results
C             -> Q, Q_T, y, CSS L0 & A3, Pert-Asymp L0:
              PRINT 132,   QG(K),VPT,RAPIN,TEMP1(1,1),TEMP1(1,3),
     &                      TEMP1(2,1),TEMP1(2,3)
              WRITE(22,132) QG(K),VPT,RAPIN,TEMP1(1,1),TEMP1(1,3),
     &                       TEMP1(2,1),TEMP1(2,3)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.1) THEN ! NLO_Sig piece
C -----------------------------------------------
C QT < QT_SEP perturbative piece.
C -----------------------------------------------
        LTOPT = 1
        WRITE(22,*) ' Q,qT,y, NLO_Sig (L0,A3) '
        PRINT *,    ' Q,qT,y, NLO_Sig (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,3,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Singular L0, A3:
              PRINT 122,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,122) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.2) THEN ! Singular piece
C ---------------------------------
C Asymptotic (singular, NLO) piece.
C ---------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3) '
        PRINT *, ' Q,qT,y, Singular (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,3,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Singular (asymptotic) L0, A3:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.3) THEN
C -----------------------------------------------
C For Y piece: perturbative and asymptotic parts.
C -----------------------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        PRINT *, ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
c        Print*, ' IQMIN, IQMAX, IYMIN, IYMAX, IPTMIN, IPTMAX '
c        Print*,   IQMIN, IQMAX, IYMIN, IYMAX, IPTMIN, IPTMAX
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,4
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                if (abs(asy).lt.1e-15) then
                    TEMP1(2,I_PERT+1) = 0.d0
                else   
                  TEMP1(2,I_PERT+1)=ASY
                endif 
                if (abs(pert).lt.1e-15) then
                    TEMP1(3,I_PERT+1) = 0.d0
                else    
                  TEMP1(3,I_PERT+1) = PERT
                endif
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y:
              PRINT 103,    QG(K),VPT,RAPIN
              WRITE(22,103) QG(K),VPT,RAPIN
C             -> Singular (asymptotic) L0, A3 and perturbative L0, A3:
              PRINT 107,    TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
              WRITE(22,107) TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
C             -> Perturbative A1, A2, A4:
              PRINT 105,    TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
              WRITE(22,105) TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN
      ENDIF

  100 FORMAT(1X,5(G10.4,2X),G14.7)
  101 FORMAT(1X,A14,2X,4(G10.4,2X))
  102 FORMAT(1X,3(f9.4,2X),4(G14.8,2X))
  103 FORMAT(1X,3(G10.4,2X))
  104 FORMAT(1X,2(G14.8,2X))
  105 FORMAT(1X,3(G18.12,2X))
  106 FORMAT(2X,3(G10.4,2X),2(G16.10,2X))
  107 FORMAT(1X,4(G18.12,2X))
  112 FORMAT(1X,3(f9.4,2X),9(G14.8,2X))
  122 FORMAT(1X,G11.5,2X,2(f9.4,2X),3(G14.8,2X))
  132 FORMAT(1X,3(f9.4,2X),4(G14.8,2X))
CsB <-
  300 FORMAT(A3,8X,A8,6X,A4,12X,A4,13X,A9)
  400 FORMAT(2(G10.4,2X),3(G14.7,2X))
C
  600 CONTINUE
      CLOSE(22)

      END

C --------------------------------------------------------------------------
      SUBROUTINE HXSect
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 VMAS
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      INTEGER I_GEES
      REAL*8 RDUMP
      LOGICAL FIRST_QT
cc      INTEGER I_PROC
cc      COMMON / PARTPROC / I_PROC
      Real*8 pT_Stretch

      pT_Stretch = 1.d0
CCPY Sep 2006: switch off the hard-wired pT_Stretch value
C      If (ECM.GT.5.d3) pT_Stretch = 1.2d0
C      pT_Stretch = 2.d0

      IF (LTO.EQ.-1) THEN ! LO piece
C ----------------------------------------------
C LO pieces: LO L0 & A3 terms
C ----------------------------------------------
        LTOPT = -1
        WRITE(22,*) ' Q,qT,y, LO L0 & A3 terms'
        PRINT *,    ' Q,qT,y, LO L0 & A3 terms'
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            FIRST_QT=.TRUE.
            DO 470 I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(FIRST_QT) THEN
                VPT=0.D0
                FIRST_QT=.FALSE.
                GOTO 530
              ELSE
                GOTO 470
              ENDIF
  530         CONTINUE
C Get the LO L0 term from Pert.for
              LepAsy = 0
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,0) = Asy
C Get the LO A3 term from Pert.for = 0
c              LepAsy = 1
c              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
c     >                       PERT,ASY,IER_PERT)
c              Temp1(0,1) = Asy
              Temp1(0,1) = 0.d0
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
              PRINT 102, QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
  470       CONTINUE   ! DO QT
          ENDDO        ! DO Y
        ENDDO          ! DO Q

      Else IF (LTO.EQ.0) THEN
        LTOPT = 0
C -----------------
C Resummed L0 piece
C -----------------
        WRITE(22,*) ' Q,qT,y, CSS, Perturbative-Asymtotic, Pert '
        PRINT *,    ' Q,qT,y, CSS, Perturbative-Asymtotic, Pert '

        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
CsB_______Set up the grid files with the convolutions of C-functions
          iqin=k
          if (iFast.eq.1) call SetCf(no_asym)
          lepasy=0
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(LTO.EQ.0) THEN
CsB Resummed L0 piece GG -> H0
                I_RESET=1
                CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                  PERT,ASY,IER_RES)
                TEMP1(1,0) = FRESUM
              ELSE
                TEMP1(1,0) = 0.D0
              ENDIF
CsB Perturbative and asymptotic L0 pieces
cpn              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
cpn     >                        ASY,IER_PERT)
cpn              TEMP1(2,0) = Pert-Asy
              Pert = 0d0
C Output the results
C             -> Q, Q_T, y, CSS, Perturbative-Asymtotic:
              PRINT 104,    QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(2,0),Pert
              WRITE(22,104) QG(K),VPT,RAPIN,
     &                      TEMP1(1,0),TEMP1(2,0),Pert
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)

CCPY Sept 2006: add the option LTO.eq.1
C This corresponds to  LTOPT.eq.1 option in PERT.FOR
      ELSE IF (LTO.EQ.1) THEN ! NLO_Sig piece
C -----------------------------------------------
C QT < QT_SEP perturbative piece.
C -----------------------------------------------
        LTOPT = 1
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3) '
        PRINT *, ' Q,qT,y, Singular (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(2,       4) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Asymptotic (Singular) L0, Perturbative:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.2) THEN ! Singular piece
C ---------------------------------
C Asymptotic (singular, NLO) piece.
C ---------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3) '
        PRINT *, ' Q,qT,y, Singular (L0,A3 '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT                
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Asymptotic (Singular) L0, Perturbative:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.3) THEN
C -----------------------------------------------
C For Y piece: perturbative and asymptotic parts.
C -----------------------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        PRINT *, ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
c        Print*, ' IQMIN, IQMAX, IYMIN, IYMAX, IPTMIN, IPTMAX '
c        Print*,   IQMIN, IQMAX, IYMIN, IYMAX, IPTMIN, IPTMAX
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y:
              PRINT 103,    QG(K),VPT,RAPIN
              WRITE(22,103) QG(K),VPT,RAPIN
C             -> Singular (asymptotic) L0, A3 and perturbative L0, A3:
              PRINT 107,    TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
              WRITE(22,107) TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
C             -> Perturbative A1, A2, A42:
              PRINT 105,    TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
              WRITE(22,105) TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ENDIF

  102 FORMAT(1X,3(f9.4,2X),2(G14.8,2X))
  104 FORMAT(1X,3(f9.4,2X),3(G14.8,2X))

  103 FORMAT(1X,3(G10.4,2X))
  105 FORMAT(1X,3(G18.12,2X))
  107 FORMAT(1X,4(G18.12,2X))

      CLOSE(22)
      RETURN
      END

CJI Jan 2015: Add H+1Jet
C --------------------------------------------------------------------------
      SUBROUTINE HJXSect
C --------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'common.for'
      INTEGER I,J,K,IER_RES,IER_PERT
      REAL*8 VPT,GEES(3),TEMP1(0:30,0:30),TEMP2,RAPIN,Q0IN
      Common /NonPertC/ GEES,Q0IN
      REAL*8 FRESUM,FGETPERT,ASY,PERT,PMA
      REAL*8 H0APPROX, H1APPROX
      REAL*8 VMAS
      COMMON/BOSONMASS/ VMAS
      INTEGER IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     >        IQMIN, IQMAX, IQSTP

      INTEGER N_Q,N_QT,N_Y
CCPY THE MAXIMAUM VALUES OF N_Q,N_QT AND N_Y ARE  30.
      PARAMETER (N_Q=9,N_QT=24,N_Y=13)

      CHARACTER*40 QGFN, QTGFN, YGFN
      COMMON / GRIDFILE / QG,PT,Y, QGFN,QTGFN,YGFN
      INTEGER II, N_Q_GD,N_QT_GD,N_Y_GD, LTO
      COMMON / NGRID / N_Q_GD,N_QT_GD,N_Y_GD, LTO

      REAL*8 QG(200),PT(200),Y(200)
      INTEGER I_Q,I_QT,I_Y
      REAL*8 A_QGRID(N_Q),W_QGRID(N_Q),Z_QGRID(N_Q),G_QGRID(N_Q)
      REAL*8 A_QTGRID(N_QT),W_QTGRID(N_QT),Z_QTGRID(N_QT),G_QTGRID(N_QT)
      REAL*8 A_YGRID(N_Y)
      REAL*8 QGRID(N_Q),YGRID(N_Y)

      EXTERNAL FGETPERT, H0APPROX, H1APPROX
      LOGICAL TESTING
      COMMON/IMMS/ IPTMIN, IPTMAX, IPTSTP, IYMIN, IYMAX, IYSTP,
     &             IQMIN, IQMAX, IQSTP
      CHARACTER*10 PDF_EVL_NAME
      COMMON/FILE_EVL/PDF_EVL_NAME

      INTEGER I_RESET
      COMMON/I_DIVDIF/I_RESET

      INTEGER I_GEES
      REAL*8 RDUMP
      LOGICAL FIRST_QT
cc      INTEGER I_PROC
cc      COMMON / PARTPROC / I_PROC
      Real*8 pT_Stretch

      pT_Stretch = 1.d0
CCPY Sep 2006: switch off the hard-wired pT_Stretch value
C      If (ECM.GT.5.d3) pT_Stretch = 1.2d0
C      pT_Stretch = 2.d0

      IF (LTO.EQ.-1) THEN ! LO piece
C ----------------------------------------------
C LO pieces: LO L0 & A3 terms
C ----------------------------------------------
        LTOPT = -1
        WRITE(22,*) ' Q,qT,y, LO L0 & A3 terms'
        PRINT *,    ' Q,qT,y, LO L0 & A3 terms'
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            FIRST_QT=.TRUE.
            DO 470 I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(FIRST_QT) THEN
                VPT=0.D0
                FIRST_QT=.FALSE.
                GOTO 530
              ELSE
                GOTO 470
              ENDIF
  530         CONTINUE
C Get the LO L0 term from Pert.for
              LepAsy = 0
              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                       PERT,ASY,IER_PERT)
              Temp1(0,0) = Asy
C Get the LO A3 term from Pert.for = 0
c              LepAsy = 1
c              TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
c     >                       PERT,ASY,IER_PERT)
c              Temp1(0,1) = Asy
              Temp1(0,1) = 0.d0
              WRITE(22,102) QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
              PRINT 102, QG(K),VPT,RAPIN,
     &                    Temp1(0,0),Temp1(0,1)
  470       CONTINUE   ! DO QT
          ENDDO        ! DO Y
        ENDDO          ! DO Q

      Else IF (LTO.EQ.0) THEN
        LTOPT = 0
C -----------------
C Resummed L0 piece
C -----------------
        WRITE(22,*) ' Q,qT,y, W~ '
        PRINT *,    ' Q,qT,y, W~ '

        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
CsB_______Set up the grid files with the convolutions of C-functions
          iqin=k
          if (iFast.eq.1) call SetCf(no_asym)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              IF(LTO.EQ.0) THEN
                lepasy=0
                I_RESET=1
                CALL SRESUM(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,
     >                  PERT,ASY,IER_RES)
                TEMP1(1,0) = FRESUM
              ELSE
                TEMP1(1,0) = 0.D0
              ENDIF
C Output the results
C             -> Q, Q_T, y, CSS b-dep, b-inp, h0 app, h1 app:
              PRINT 104,    QG(K),VPT,RAPIN,
     &              TEMP1(1,0),TEMP1(2,0),TEMP1(3,0),
     &              TEMP1(4,0)
              WRITE(22,104) QG(K),VPT,RAPIN,
     &              TEMP1(1,0),TEMP1(2,0),TEMP1(3,0),TEMP1(4,0)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)

CCPY Sept 2006: add the option LTO.eq.1
C This corresponds to  LTOPT.eq.1 option in PERT.FOR
      ELSE IF (LTO.EQ.1) THEN ! NLO_Sig piece
C -----------------------------------------------
C QT < QT_SEP perturbative piece.
C -----------------------------------------------
        LTOPT = 1
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3) '
        PRINT *, ' Q,qT,y, Singular (L0,A3) '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(2,       4) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Asymptotic (Singular) L0, Perturbative:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.2) THEN ! Singular piece
C ---------------------------------
C Asymptotic (singular, NLO) piece.
C ---------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3) '
        PRINT *, ' Q,qT,y, Singular (L0,A3 '
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0,3
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(2,I_PERT+1) = ASY
                TEMP1(3,I_PERT+1) = PERT                
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y, Asymptotic (Singular) L0, Perturbative:
              PRINT 102,    QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
              WRITE(22,102) QG(K),VPT,RAPIN,TEMP1(2,1),TEMP1(2,4)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ELSE IF (LTO.EQ.3) THEN
C -----------------------------------------------
C For Y piece: perturbative and asymptotic parts.
C -----------------------------------------------
        LTOPT = 0
        WRITE(22,*) ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
        PRINT *, ' Q,qT,y, Singular (L0,A3), Pert. (L0,A3,A1,A2,A4) '
c        Print*, ' IQMIN, IQMAX, IYMIN, IYMAX, IPTMIN, IPTMAX '
c        Print*,   IQMIN, IQMAX, IYMIN, IYMAX, IPTMIN, IPTMAX
        DO K = IQMIN, IQMAX, IQSTP
          VMAS = QG(K)
          !TEMP1(3,0) = H0Approx(VMAS)
          !TEMP1(4,0) = H1Approx(VMAS)
          DO J = IYMIN, IYMAX, IYSTP
            RAPIN = Y(J)
            DO I = IPTMIN, IPTMAX, IPTSTP
              VPT=PT(I)*pT_Stretch
              LEPASY = 0 ! Ex symmetric piece
              DO I_PERT = 0,0
                IF (I_PERT.EQ.3) LEPASY = 1 ! Ex anti-symmetric piece
C Call the computation
                TEMP2=FGETPERT(VMAS,VPT,RAPIN,GEES,Q0IN,FRESUM,PMA,PERT,
     >                         ASY,IER_PERT)
                TEMP1(3,0) = H0Approx(VMAS,vpt,rapin)
                TEMP1(2,I_PERT+1) = ASY*TEMP1(3,0)
                TEMP1(3,I_PERT+1) = PERT
              END DO ! I_Pert
C Output the results
C             -> Q, Q_T, y:
              PRINT 103,    QG(K),VPT,RAPIN
              WRITE(22,103) QG(K),VPT,RAPIN
C             -> Singular (asymptotic) L0, A3 and perturbative L0, A3:
              PRINT 107,    TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
              WRITE(22,107) TEMP1(2,1),TEMP1(2,4), TEMP1(3,1),TEMP1(3,4)
C             -> Perturbative A1, A2, A42:
              PRINT 105,    TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
              WRITE(22,105) TEMP1(3,2), TEMP1(3,3), TEMP1(3,5)
            END DO    ! I (QT)
          ENDDO       ! J (Y)
        ENDDO         ! K (Q)
        CLOSE(22)
        RETURN

      ENDIF

  102 FORMAT(1X,3(f9.4,2X),2(G14.8,2X))
  104 FORMAT(1X,3(f9.4,2X),4(G14.8,2X))

  103 FORMAT(1X,3(G10.4,2X))
  105 FORMAT(1X,3(G18.12,2X))
  107 FORMAT(1X,4(G18.12,2X))

      CLOSE(22)
      RETURN
      END

CsB ------------------------------------------------------------------------

CsB___Running alpha_em from PYTHIA
C --------------------------------------------------------------------------
      FUNCTION PYALEM(Q2)
C --------------------------------------------------------------------------
C...Double precision and integer declarations.
C      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      Real*8 PYALEM,Q2
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      REAL PARU(200)
      INTEGER MSTU(200)
C      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C      SAVE /PYDAT1/

C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.
CsB___See also Phys.Lett.B356:398-403,1995

      PARU(1)=3.1415927
      PARU(101)=1./137.04
      MSTU(101)=1
      PARU(104)=1.
      PARU(103)=1./128.8
      AEMPI=PARU(101)/(3D0*PARU(1))
      IF(MSTU(101).LE.0.OR.Q2.LT.2D-6) THEN
        RPIGG=0D0
      ELSEIF(MSTU(101).EQ.2.AND.Q2.LT.PARU(104)) THEN
        RPIGG=0D0
      ELSEIF(MSTU(101).EQ.2) THEN
        RPIGG=1D0-PARU(101)/PARU(103)
      ELSEIF(Q2.LT.0.09D0) THEN
        RPIGG=AEMPI*(13.4916D0+LOG(Q2))+0.00835D0*LOG(1D0+Q2)
      ELSEIF(Q2.LT.9D0) THEN
        RPIGG=AEMPI*(16.3200D0+2D0*LOG(Q2))+
     &  0.00238D0*LOG(1D0+3.927D0*Q2)
      ELSEIF(Q2.LT.1D4) THEN
        RPIGG=AEMPI*(13.4955D0+3D0*LOG(Q2))+0.00165D0+
     &  0.00299D0*LOG(1D0+Q2)
      ELSE
        RPIGG=AEMPI*(13.4955D0+3D0*LOG(Q2))+0.00221D0+
     &  0.00293D0*LOG(1D0+Q2)
      ENDIF

C...Calculate running alpha_em.
      PYALEM=PARU(101)/(1D0-RPIGG)
      PARU(108)=PYALEM

      RETURN
      END






