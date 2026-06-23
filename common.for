CCPY MARCH 2005
      LOGICAL DEBUG,GOAROUND,TESTASY
      INTEGER NIN,NOUT,NWRT
      REAL*8 UNIT,HBARC2,CONST
      INTEGER NC,NF,G
      REAL*8 PI,PI2,CF,CA,TR,BETA1,EULER,ZETA3
      REAL*8 A1,A2,A3,A4,B1,B2,B3,C1,C2,C3,C4,B0
      INTEGER IFLAG_C3, IFLAG_mu, FLAG_mu
      REAL*8 BMAX,Q0,G1,G2,G3
      REAL*8 P00,P02,P04,Q01,Q03
      REAL*8 ZEROJ0,ACC_SUDNPT,ACC_AERR,ACC_RERR,ACC_RES1
      INTEGER NRES1
      REAL*8 MT,MW,MZ,MH,MA,MW2,MZ2,
     ,WCOUPL,ZCOUPL,HCOUPL,ACOUPL,HpCOUPL,HbCOUPL
      REAL*8 VKM,ZFF_A,ZFF_V,ZFF2
      Real*8 EW_ALFA,GFermi,GAMW,GAMZ,GAMW2,GAMZ2
      Real*8 VEV,TanB
      INTEGER ISET,NORDER
      REAL*8 YNPT_CUT
      REAL*8 S_Q,S_LAM,S_BETA1
      REAL*8 QT_V,Y_V,Q_V,Y2_V
      REAL*8 X_A,X_B,ECM,ECM2,YMAX,GMU
      INTEGER IBEAM,LEPASY,LTOPT,i_Model,i_RunMass
      LOGICAL FIT3
      REAL*8 PREV_UNDER
      LOGICAL FIRST_YNPT
      INTEGER INONPERT
      INTEGER IERROR_XMIN,IERROR_XMAX,IERROR_QSQMIN,IERROR_QSQMAX
      REAL*8 APDF,APDF_PION_MINUS
      EXTERNAL APDF,APDF_PION_MINUS
      REAL*8 FRACT_N
      INTEGER IHADRON
      REAL*8 qmas_u,qmas_d,qmas_s,qmas_c,qmas_b,qmas_t
	Integer No_Sigma0

      INTEGER N_SUD_A,N_SUD_B,N_WIL_C,I_FSR
      REAL*8 ALPI,ALAMBD
      INTEGER NFL
CCPY
      real*8 CT14Alphas
      EXTERNAL CT14ALphas


c      LOGICAL ASK
CCPY      CHARACTER*2 TYPE_V
C      CHARACTER*10 PDF_FILE
      CHARACTER*40 TYPE_V
      CHARACTER*40 PDF_FILE
      EXTERNAL ALPI,NFL,ALAMBD
      Integer I_Pert, iCF1
      integer initlhapdf, iset_save
      COMMON/LHAPDF/ initlhapdf, iset_save

CsB   I_Proc = 1 Calculate only QI QJ -->  G V process
C     I_Proc = 2 Calculate only  G QI --> QJ V process
      INTEGER I_PROC
      COMMON / PARTPROC / I_PROC

      COMMON/SETUP1/ DEBUG,NIN,NOUT,NWRT
      COMMON/SETUP2/ UNIT,HBARC2,CONST,GMU
      COMMON/SETUP3/ PI,PI2,CF,CA,TR,BETA1,EULER,ZETA3,NC,NF
      COMMON/SETUP4/ A1,A2,A3,A4,B1,B2,B3,C1,C2,C3,
     >C4,B0,IFLAG_C3, IFLAG_mu,FLAG_mu
      COMMON/SETUP5/ BMAX,Q0,G1,G2,G3
      COMMON/SETUP6/ P00,P02,P04,Q01,Q03
      COMMON/SETUP7/ ZEROJ0(0:29),ACC_SUDNPT,ACC_AERR,
     >ACC_RERR,ACC_RES1,NRES1
      COMMON/SETUP8/ FIT3
      COMMON/SETUP9/ PREV_UNDER
      COMMON/STAND1/ MT,MW,MZ,MH,MA,MW2,MZ2,
     >WCOUPL,ZCOUPL,HCOUPL,ACOUPL,HpCOUPL,HbCOUPL
      COMMON/STAND2/ VKM(6,6),ZFF_A(6),ZFF_V(6),ZFF2(6)
      COMMON/STAND3/ EW_ALFA,GFermi,GAMW,GAMZ,GAMW2,GAMZ2
      Common/Stand4/ VEV,TanB
      COMMON/QCD1/ ISET,NORDER
      COMMON/NPERT/ YNPT_CUT,FIRST_YNPT,INONPERT
      COMMON/SUD1/ S_Q,S_LAM,S_BETA1
      COMMON/INPUT1/ QT_V,Y_V,Q_V,TYPE_V,Y2_V
      COMMON/INPUT2/ X_A,X_B,ECM,ECM2,YMAX,IBEAM,LEPASY,LTOPT,
     ,i_Model,i_RunMass
      COMMON/CFCF/ GOAROUND,TESTASY
      COMMON/MYHMRS/ IERROR_XMIN,IERROR_XMAX,IERROR_QSQMIN,IERROR_QSQMAX
      COMMON/DRELL/ FRACT_N,IHADRON
      COMMON/QMASES/ qmas_u,qmas_d,qmas_s,qmas_c,qmas_b,qmas_t
      COMMON/FILINP/ PDF_FILE 
      Common/Sigma0/ No_Sigma0
      COMMON/N_ABC/ N_SUD_A,N_SUD_B,N_WIL_C,I_FSR
C I_Pert Assignments. I_Pert:Piece -> 0:L0, 1:A1, 2:(A2+A0), 3:2*A3, 4:A4.
      Common/Pert_Piece/ I_Pert
	Common/ CFnSwitch / iCF1
CPN
      integer iFast
      common/Fast/iFast 
      integer Nconv,iqin,iyin,iQTin
      real*8 sml_b
      parameter(Nconv=200) 
      real*8 ConvGrdS(100,0:Nconv),ConvGrdA(100,0:Nconv),AmuGrd(0:Nconv)
      common/c_conv/ConvGrdS,ConvGrdA,AmuGrd,sml_b,iqin,iyin,iqtin
cpn                          Switches to determine if to calculate the
cpn                     asym. piece in SetCf
      logical yes_asym, no_asym
      data yes_asym, no_asym /.true.,.false./

ccpy
      REAL*8 SW2,SW2_EFF
      Common / WMA / SW2,SW2_EFF

      INTEGER I_CXFF
      COMMON/NLO_SIG/ I_CXFF

cpn                    For heavy quarks
      double precision xmhq,B_INT
      Integer iHQMass,ihqpdf
      common/heavyquark/xmhq,iHQMass,ihqpdf
      COMMON/CONV_INT_HQ/ B_INT

CJI ADD Option to use CFG or CSS resummation
      Integer ResumType
      COMMON /RESUMFORM/ ResumType

      real*8 mur, muf
      common /delscale/ mur, muf
      INTEGER iscale
      common/petscale/ iscale
