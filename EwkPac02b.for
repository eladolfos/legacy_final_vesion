      FUNCTION AlfaEm(Amu)
C                                                   -=-=- alfaem
C --------------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C                                                     Adapted from PYTHIA
C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.

      Pi = 3.1415927
      Alf0 = 1./137.04
      Alpi = Alf0/(3D0*Pi)
C                                                         AlfZ = 1./128.8
      Am2  =  Amu **2
C                                     approximate order alpha correction :
      IF(Amu .LT. 1.4D-3) Then
        Del = 0D0
      ElseIf(Amu .LT. 0.3D0) Then
        Del = Alpi*(13.4916D0 +LOG(Am2)) +0.00835D0*LOG(1D0 +Am2)
      ElseIf(Amu .LT. 3D0) Then
        Del = Alpi*(16.3200D0 + 2D0*LOG(Am2)) + 
     &  0.00238D0*LOG(1D0 + 3.927D0*Am2)
      ElseIf(Amu .LT. 1D2) Then
        Del = Alpi*(13.4955D0 + 3D0*LOG(Am2)) + 0.00165D0 + 
     &  0.00299D0*LOG(1D0 + Am2)
      Else
        Del = Alpi*(13.4955D0 + 3D0*LOG(Am2)) + 0.00221D0 + 
     &  0.00293D0*LOG(1D0 + Am2)
      EndIF
C...                                                     running alpha_em.
      AlfaEm = Alf0/(1D0-Del)

      RETURN
C                     *************************
      End
C                                                          =-=-= EwCpl1
      Function ALFEWK (IBOSON)
C                                                   -=-=- alfewk

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: EwCpl0
C Description: Function calls to extract the Ewk couplings + Running alfa_em
C ListOfFiles: alfewk alfaem
C ====================================================================
C Entry function not shown: (vbnmas swg2f gewlt gewlh gewqt gewqh) 

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/EwCpl0.f,v 6.1 98/08/16 17:21:35 wkt Exp $
C #Log:	EwCpl0.f,v $
c Revision 6.1  98/08/16  17:21:35  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 

C     +++++++++++++++++++++++++++++  Functions to extract EWK coefficients

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EWKPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / KMATRX / VKM (NGN, NGN)

C            Left, Right, Vector, and Axial-vector couplings of Leptons & Quarks

      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 

C                                 "Alpha QED" for all Electo-Weak Vector Bosons
      ALFEWK = ALFEW(IBOSON)

      Return
C                        ****************************

      Entry VBNMAS (IBOSON)

      IF     (IBOSON .EQ. 1) THEN
        BM = 0.0
      ELSEIF (IBOSON .EQ. 2 .OR. IBOSON .EQ. 3) THEN
        BM = WMS
      ELSEIF (IBOSON .EQ. 4) THEN
        BM = ZMS
      ELSE
        WRITE (NOUT, '(2A, I5)') ' Vector Boson Index out of range',
     >  ' in VBNMAS; IBOSON =', IBOSON
      ENDIF

      VBNMAS = BM

      Return
C                        ****************************

      Entry SWG2F ()

      SWG2F = SWG2

      Return
C                        ****************************

      Entry GEWLT (IT1, IBS, IBT, IT2)

C           (G) Electro-Weak coupling for Leptons in the Tensor base.
C            -  -       -                 -              -

C               To be used in the calculation of general E-WK matrix elements.

C           IT1, IT2 = 1,2 :  Weak Isospin (T3) of the two leptons
C           IBS =  1 -4    :  Boson label (see SETEWK)
C           IBT =             Boson polarization (tensor) label
C                  1       :  vector
C                 -1       :  axial-vector
C                  other   :  illegal

C                   For now, this still applies to fields rather then particles;
C                       hence there is no distinction between particle & anti-p.

C                 To distinquish between them, either write another module, or
C                       use the convention, IT1 = -1, -2 for the anti-part;
C                            and test the signs upon Entry to determine the 
C                            channel (i.e. scattering, decay, or production).

      IF ((IT1.LT.1 .OR. IT1.GT.NSP) .OR. (IT2.LT.1 .OR. IT2.GT.NSP)) 
     >  THEN
        WRITE (NOUT,*)' Lepton Label out of range in GEWLT; IT1,IT2 =',
     >  IT1, IT2
        STOP
      ENDIF

      IF     (IBT .EQ. 1) THEN
         TEM = GLV (IT1, IBS, IT2)
      ELSEIF (IBT .EQ.-1) THEN
         TEM = GLA (IT1, IBS, IT2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson Tensor index out of range in GEWLT; IBT = ', IBT
         STOP
      ENDIF

      GEWLT = TEM

      Return
C                        ****************************

      Entry GEWLH (IT1, IBS, IBH, IT2)

C           (G) Electro-Weak coupling for Leptons in the cHirality base.
C            -  -       -                 -               -

C               To be used in the calculation of general E-WK matrix elements.

C           IT1, IT2 = 1,2 :  Weak Isospin (T3) of the two leptons
C           IBS =  1 -4    :  Boson label (see SETEWK)
C           IBH =             Boson polarization (helicity) label
C                  1       :  right-handed
C                 -1       :  left -handed
C                  other   :  illegal

C                 For now, this still just apply to fiels rather then particles;
C                       hence there is no distinction between particle & anti-p.

C                 To distinquish between them, either write another module, or
C                       use the convention, IT1 = -1, -2 for the anti-part;
C                            and test the signs upon Entry to determine the 
C                            channel (i.e. scattering, decay, or production).

C         When make the field ---> particle transition, cHirality ---> Helicity

      IF ((IT1.LT.1 .OR. IT1.GT.NSP) .OR. (IT2.LT.1 .OR. IT2.GT.NSP)) 
     >  THEN
        WRITE (NOUT,*)' Lepton Label out of range in GEWLH; IT1,IT2 =',
     >  IT1, IT2
        STOP
      ENDIF

      IF     (IBH .EQ. 1) THEN
         TEM = GLR (IT1, IBS, IT2)
      ELSEIF (IBH .EQ.-1) THEN
         TEM = GLL (IT1, IBS, IT2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson cHirality index out of range in GEWLH; IBH = ', IBH
         STOP
      ENDIF

      GEWLH = TEM

      Return
C                        ****************************

      Entry GEWQT (IQ1, IBS, IBT, IQ2)

C           (G) Electro-Weak coupling for Quarks in the Tensor base.
C            -  -       -                 -             -

C               To be used in the calculation of general E-WK matrix elements.

C           IBT =  1     :  Vector
C                 -1     :  Axial-vector
C                  other :  illegal

      IF ((IQ1.LT.1 .OR. IQ1.GT.NFL) .OR. (IQ2.LT.1 .OR. IQ2.GT.NFL)) 
     >  THEN
        WRITE (NOUT,*)' Quark Label out of range in GEWQT; IQ1,IQ2 =',
     >  IQ1, IQ2
        STOP
      ENDIF

      IF     (IBT .EQ. 1) THEN
         TEM = HQV (IQ1, IBS, IQ2)
      ELSEIF (IBT .EQ.-1) THEN
         TEM = HQA (IQ1, IBS, IQ2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson Tensor index out of range in GEWQT; IBT = ', IBT
         STOP
      ENDIF

      GEWQT = TEM

      Return
C                        ****************************

      Entry GEWQH (IQ1, IBS, IBH, IQ2)

C           (G) Electro-Weak coupling for Quarks in the cHirality base.
C            -  -       -                 -             -

C               To be used in the calculation of general E-WK matrix elements.

C           IBH =  1     :  right-handed
C                 -1     :  left -handed
C                  other :  illegal

      IF ((IQ1.LT.1 .OR. IQ1.GT.NFL) .OR. (IQ2.LT.1 .OR. IQ2.GT.NFL)) 
     >  THEN
        WRITE (NOUT,*)' Quark Label out of range in GEWQH; IQ1,IQ2 =',
     >  IQ1, IQ2
        STOP
      ENDIF

      IF     (IBH .EQ. 1) THEN
         TEM = HQR (IQ1, IBS, IQ2)
      ELSEIF (IBH .EQ.-1) THEN
         TEM = HQL (IQ1, IBS, IQ2)
      ELSE
         WRITE (NOUT, *) 
     > ' IBoson cHirality index out of range in GEWQH; IBH = ', IBH
         STOP
      ENDIF

      GEWQH = TEM

      Return
C                        ****************************
      END

      Function EwCpl2An (JP1, JBN, JP2)
C                                                   -=-=- ewcplg2

C     Function call of the squared EW coupling constants for Drell-Yan
C     processes in the Annihilation channel.

C     EwCplC common block is set up by the SetEwCpl2 subroutine

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NFL = 6, NBN = 4)

      COMMON / EWCPLC / CPLANH(-NFL:NFL, NBN, -NFL:NFL), 
     >                  CPLSCT(-NFL:NFL, NBN)

      EwCpl2An = CPLANH (JP1, JBN, JP2)

      Return
C                        ****************************

      Entry EwCpl2Cn (JP, JBN)

C     Function call of the squared EW coupling constants for Drell-Yan
C     processes in the Compton Sc. channel.

      EwCpl2Cn = CPLSCT (JP, JBN)

      Return
C                        ****************************
      END

      Function EwCplSc (Iboson, Ihelic, Ipartn)
C                                                   -=-=- ewcplg

C                Electro-weak coupling for VecBos + Parton --> Parton2

C                Use helicity basis for the boson polarization index.

C                       Iboson (1:4)  \  Ihelic (-1:1)
C                                      \
C                                       \______________ 
C                                       /         (Parton2, redundant)
C                                      /
C                      Ipartn (-6:6)  /

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      COMMON / WKCPLG / QCPLN(4, -1:1, -6:6)

      SAVE

      If (  (Iboson.Lt. 1 .or. Iboson.Gt.4) 
     > .or. (Ihelic.Lt.-1 .or. Ihelic.Gt.1) 
     > .or. (Ipartn.Lt.-6 .or. Ipartn.Gt.6) ) Then
         Print '(/ A / A, 3I4 /)'
     >,'Iboson, Ihelicity, Iparton out of range in EwCplSc; '
     >,'incorrect input values : ', Iboson, Ihelic, Ipartn
         Stop
      EndIf

      tem = QCPLN (Iboson, Ihelic, Ipartn)

      EwCplSc = tem

      RETURN
C                        ****************************

      Entry EwCplAn (Iboson, Ihelic, Ipartn)

C           Electro-weak coupling for Parton + Parton2bar ---> VecBoson

C                      Ipartn (-6:6)  \ 
C                                      \
C                                       \______________Iboson (1:4) 
C                                       /              Ihelic (-1:1)
C                                      /
C                         Parton2bar  /

      Print '(/ A / A /)', 'EwCplAn has not yet been implemented.'
     > , 'Can be obtained by crossing from EwCplSc above.'

      EwCplAn = 0.

      END

C                                                          =-=-= EwCpl2
      Subroutine FILCPL (CH, GV, GA, GL, GR)
C                                                   -=-=- filcpl

C            Given the charge of the fermion, compute the ElectroWeak couplings 
C            for a given generation.  Used in SetEwk.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)

      COMMON / EWKPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)

      DIMENSION CH(NSP), GV(NSP,NBN,NSP), GA(NSP,NBN,NSP), 
     >          T3(NSP), GL(NSP,NBN,NSP), GR(NSP,NBN,NSP)

      DATA T3 / 0.5, -0.5 /

      DO 31 IS1 = 1, NSP
      DO 31 IS2 = 1, NSP
        IF     (IS1 .EQ. IS2) THEN
            GV (IS1, 1, IS2) = CH(IS1)
            GA (IS1, 1, IS2) = 0.0
            GV (IS1, 2, IS2) = 0.0
            GA (IS1, 2, IS2) = 0.0
            GV (IS1, 3, IS2) = 0.0
            GA (IS1, 3, IS2) = 0.0
            GV (IS1, 4, IS2) = T3(IS1) - 2.* CH(IS1) * SWG2
            GA (IS1, 4, IS2) =-T3(IS1)
        ELSEIF (IS1 .GT. IS2) THEN
            GV (IS1, 1, IS2) = 0.0
            GA (IS1, 1, IS2) = 0.0
            GV (IS1, 2, IS2) = 1.0
            GA (IS1, 2, IS2) =-1.0
            GV (IS1, 3, IS2) = 0.0
            GA (IS1, 3, IS2) = 0.0
            GV (IS1, 4, IS2) = 0.0
            GA (IS1, 4, IS2) = 0.0
        ELSEIF (IS1 .LT. IS2) THEN
            GV (IS1, 1, IS2) = 0.0
            GA (IS1, 1, IS2) = 0.0
            GV (IS1, 2, IS2) = 0.0
            GA (IS1, 2, IS2) = 0.0
            GV (IS1, 3, IS2) = 1.0
            GA (IS1, 3, IS2) =-1.0
            GV (IS1, 4, IS2) = 0.0
            GA (IS1, 4, IS2) = 0.0
        ENDIF

        DO 32 IBN = 1, NBN
            GL (IS1,IBN,IS2) = (GV(IS1,IBN,IS2) - GA(IS1,IBN,IS2)) / 2.
            GR (IS1,IBN,IS2) = (GV(IS1,IBN,IS2) + GA(IS1,IBN,IS2)) / 2.
   32   CONTINUE

   31 CONTINUE


      Return
C                        ****************************
      END

      Subroutine FILKMM 
C                                                   -=-=- filkmm
C                       Given the Mixing parameters for NGN generations, 
C                       calculate the KM matrix in the Cabbibo-KM-Maiani-
C                       -Wolfenstein-Chau-Keung...etc scheme.
C                       Used in SetEwk.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NGN = 3, NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)
 
      COMMON / KMATRX / VKM (NGN, NGN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)

C     (how to handle complex numbers due to phase factors efficiently,
C      bearing in mind that most applications do not need this info??)

C     For the moment, put the absolute values of the matrix elements in by hand

      DIMENSION AK (3, 3)
c CKM values from Particle data group
c W.-M. Yao et al., J. of Phys. G 33, 1 (2006)
      DATA (AK(I,1), I=1,3) / 0.97383, 0.2272, 0.00396 /
      DATA (AK(I,2), I=1,3) / 0.2271, 0.97296, 0.04221 /
      DATA (AK(I,3), I=1,3) / 0.00814, 0.04161, 0.9991 /

      DO 5 I = 1, 3
      DO 6 J = 1, 3
        VKM (I,J) = AK (I,J)
    6 CONTINUE
    5 CONTINUE

      Return
C                        ****************************
      END

C                                                          =-=-= EwCpl0
      Subroutine SetEwk
C                                                   -=-=- setewk

C 05-06-01 : Initialize the EW couplings to Zero first in Data Block.
C            Non-zero values will be filled by the program.
C            Can't assume all compilers will initialize the unset entries
C            to zero.
C 05-05-10 : trivial change: initialize Lprt ->0, so that the coupling
C            don't get printed out each time.  Bug: Lprt can't be set.

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: SetEwk
C Description: initial setup of the EwkPac.
C ListOfFiles: setewk filcpl filkmm
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/SetEwk.f,v 6.1 98/08/16 17:24:43 wkt Exp $
C #Log:	SetEwk.f,v $
c Revision 6.1  98/08/16  17:24:43  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 
C   In DisPac: QCPLN ---> EwCplSc
C   In VbpPac: SQCPAN --> EwCpl2An
C              SQCPCM --> EwCpl2Cn
c 
C Revision 6.0  96/12/25  14:38:28  wkt
C Synchronize with version 6 of all the other pac's

c Revision 1.1  94/02/22  11:10:59  lai
c Initial revision
c 
c Revision 2.2  92/03/08  15:21:24  wkt
c Setxxx.f uniformly added to force blockdata linking and 
c to perform other initiation Functions.
c 
C                 Setup the Ewkpac, which contains the following modules

C      Subroutine SetEwk
C      Subroutine FILCPL (CH, GV, GA, GL, GR)
C      Subroutine FILKMM 
C
C      Function ALFEwk (IBOSON)
C      Entry VBNMAS (IBOSON)
C      Entry SWG2F ()
C      Entry GEWLT (IT1, IBS, IBT, IT2)
C      Entry GEWLH (IT1, IBS, IBH, IT2)
C      Entry GEWQT (IQ1, IBS, IBT, IQ2)
C      Entry GEWQH (IQ1, IBS, IBH, IQ2)
C
C      Subroutine SetEwCpl2
C      Function SQCPAN (JP1, JBN, JP2)
C      Entry SQCPCM (JP, JBN)
C      Function EWCPL0 (IBSN, IP1, IP2, IRT)   (Inactive)

C     Set up the basic Electro-weak coupling matrices for the Boson-Fermion
C     Yukawa coupling term in the Effective Lagrangian

C     Boson label: (IBN)   1,   2,   3,   4
C                       gamma   W+   W-   Z    
  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN)
 
      DIMENSION IP(NSP, NGN), CHL(NSP), CHQ(NSP)
      DIMENSION GGQA(NSP,NGN,NSP,NGN), GGQV(NSP,NGN,NSP,NGN),
     >          GGQR(NSP,NGN,NSP,NGN), GGQL(NSP,NGN,NSP,NGN)
 
      External DatEwk

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EwkPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / KMATRX / VKM (NGN, NGN)

C            Left, Right, Vector, and Axial-vector couplings of Leptons & Quarks

      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 
      DATA IP, Lprt   / 1, 2, 4, 3, 6, 5, 0 /
      DATA CHL,CHQ / 0d0, -1d0, 0.666666667, -0.333333333 /
 
C     Adopt the convention:

C                |  IBSN      coupling       current           coupling**2/4pi
C                |
C   ------------ |------------------------------------------------------------
C      Photon    |    1          e             Gmu            ALFE = Alpha(QED) 
C                |
C      W+/W-     |   2/3      g/2Sqrt2      Gmu(1-G5)          ALFE / SWG2 / 8
C                |
C        Z       |    4        g/2CosW  Gmu((T3-2QSWG2)-T3G5)   ALFE / S2WG2
C                |
C   ------------ |-------------------------------------------------------------
C     where Gmu = Gamma-super-mu; G5 = Gamma-5; Cos W = Cos(Theta-Weinberg);
C     SWG2 = (Sin W)**2; & S2WG2 = (Sin 2W)**2.
C   ---------------------------------------------------------------------------
C                                       Fill the overall coupling**2 array ALFEW
CLai 11/14/06 replace the following by CP's note
c      ALFEW (1) = ALFE
c      ALFEW (2) = ALFE / SWG2 / 8.
c      ALFEW (3) = ALFEW (2)
c      ALFEW (4) = ALFE / SWG2 / (1.-SWG2) / 4.
C-----------------------------------------
C the following adopted from CP's note on EW couplings
C----------------------------------------- 
C Input parameters:

C_____W-BOSON MASS
      WMS = 80.41
C_____Z BOSON MASS      
      ZMS = 91.187
C_____FERMI CONSTANT      
      GMU = 1.16637D-5
       
C_____EFFECTIVE SIN^2(THETA_WEAK) AT Z-POLE
      SW2_EFF_MZ = 0.23143
      SWG2 = SW2_EFF_MZ

C_____RUNNING QED COUPLING AT Z-POLE
      ALFA_EM_MZ = 1.0/128.937

C_____RUNNING QED COUPLING AT ELECTRON MASS SCALE
      ALFA_EM_ME = 1.0D0/137.0359895D0

C_____PI
      PI = 3.1415927 
      
        
C------------------------------------ 
C Derived variables:

      R2 = SQRT(2.D0)
      XMW2 = WMS**2
      XMZ2 = ZMS**2
      
C_____Weak couplings: g_w
      GWEAK2 = 4.D0*R2*XMW2*GMU
      GWEAK = SQRT(GWEAK2)
      
      
C PHOTON COUPL IS (e**2)/2      
C_____RUNNING QED COUPLING AT THE SCALE AMU (CALCULATED IN 
C_____FUNCTION ALFAEM 
cLai  comment ALFE out since AMU depends on expt. 
c     Use the fine structure constant here from DatEwk
c      ALFE = ALFAEM(AMU) 
      
C WCOUPL IS (g**2)/8
      WCOUPL = GWEAK2/8.0D0
      
C WE IGNORE THE DIFFERENCE BETWEEN b-quark AND OTHER FERMIONS
C ZCOUPL IS (g/Cos_w)**2/4

C_____Sin^2 and Cos^2 of the weak angle - in the on shell scheme
      CWS = XMW2/XMZ2
      SWS = 1.D0-CWS
      ZCOUPL_MZ = PI*ALFA_EM_MZ/SWS/CWS
      ZCOUPL = ZCOUPL_MZ
      
C------------------------------------
C To a good approximation:

      ALFEW (1) = ALFE
      ALFEW (2) = WCOUPL/4.0/PI
      ALFEW (3) = ALFEW (2)
      ALFEW (4) = ZCOUPL/4.0/PI
  
C------------------------------------

C                                                   ----------------------------
C                                                       For Each Generation:
C                                     Fill the relative coupling in the currents
C                           GXL/R = (GXV -/+ GXA)/2.      (X = Lepton / Quark) ;
C                             hence the currents are:   GXL/R * Gmu * (1-/+G5)

      CALL FILCPL (CHL, GLV, GLA, GLL, GLR)
      CALL FILCPL (CHQ, GQV, GQA, GQL, GQR)
C                                                  ----------------------------
C                          For Quarks, put in Generations and generation-mixing

C                                 fill in the KM matrix array from data in the
C                           MIXing PARameter common block / MIXPAR / consisting
C                                 of Cosines of mixing angles and Phase angles.
      CALL FILKMM

C                     For W+ / W- , Multiply the KM matrix into the Down quarks
      DO 10 IG1 = 1, NGN
      DO 10 IG2 = 1, NGN
      DO 10 IS1 = 1, NSP
      DO 10 IS2 = 1, NSP
        GGQA (IS1, IG1, IS2, IG2) = GQA(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQV (IS1, IG1, IS2, IG2) = GQV(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQR (IS1, IG1, IS2, IG2) = GQR(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQL (IS1, IG1, IS2, IG2) = GQL(IS1, 2, IS2) * VKM(IG1, IG2)
   10 CONTINUE
C                                                  ----------------------------
C                         (Generation, ISPin) labels  ----> parton flavor label
      DO 20 IG1 = 1, NGN
      DO 20 IG2 = 1, NGN
      DO 20 IS1 = 1, NSP
      DO 20 IS2 = 1, NSP
        IP1 = IP(IS1, IG1)
        IP2 = IP(IS2, IG2)

        DO 21 IBN = 1, NBN
        IF (IBN .EQ. 1 .OR. IBN .EQ. 4) THEN
         IF (IG1 .NE. IG2) THEN
            HQA (IP1, IBN, IP2) = 0.0
            HQV (IP1, IBN, IP2) = 0.0
            HQL (IP1, IBN, IP2) = 0.0
            HQR (IP1, IBN, IP2) = 0.0
         ELSE
            HQA (IP1, IBN, IP2) = GQA(IS1, IBN, IS2)
            HQV (IP1, IBN, IP2) = GQV(IS1, IBN, IS2)
            HQL (IP1, IBN, IP2) = GQL(IS1, IBN, IS2)
            HQR (IP1, IBN, IP2) = GQR(IS1, IBN, IS2)
         ENDIF
        ENDIF
   21   CONTINUE

            HQA (IP1, 2, IP2) = GGQA(IS1, IG1, IS2, IG2)
            HQV (IP1, 2, IP2) = GGQV(IS1, IG1, IS2, IG2)
            HQL (IP1, 2, IP2) = GGQL(IS1, IG1, IS2, IG2)
            HQR (IP1, 2, IP2) = GGQR(IS1, IG1, IS2, IG2)

            HQA (IP1, 3, IP2) = GGQA(IS2, IG2, IS1, IG1)
            HQV (IP1, 3, IP2) = GGQV(IS2, IG2, IS1, IG1)
            HQL (IP1, 3, IP2) = GGQL(IS2, IG2, IS1, IG1)
            HQR (IP1, 3, IP2) = GGQR(IS2, IG2, IS1, IG1)

   20 CONTINUE
C                                                  ----------------------------
C                                                                      Finished
C                          Set Lprt>0 to print Quark couplings for confirmation
      If (Lprt .ge. 1) Then
        Print '(A)', 
     >  'Flavor-dependent part of Quark-VectorBoson Vertex:'
        Do 41 Ibn = 1, 4
	Print '(/A, I4)', 'Iboson = ', Ibn
	Print '(/A/(4F12.5))', 'Axial-v: Iflv=1-4:'
     >, ((HqA(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
	Print '(/A/(4F12.5))', 'Vector: Iflv=1-4:'
     >, ((HqV(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
	Print '(/A/(4F12.5))', 'Left-hand: Iflv=1-4:'
     >, ((HqL(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
	Print '(/A/(4F12.5))', 'Right-hand: Iflv=1-4:'
     >, ((HqR(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
 41     Continue
      Endif

      Return
C                        ****************************
      End

      BLOCK DATA DATEwk
C                       
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN, Ntt1=NSP*NBN*NSP, Ntt2=NFL*NBN*NFL)
      PARAMETER (NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)
 
      COMMON / EwkPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)
      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 
      DATA WMS, ZMS, SWG2, ALFE / 80.41, 91.187, 0.23124, 7.297353E-3 /
      Data (((GLL(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /
     >     (((GLR(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GLV(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GLA(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 / 
C
     >     (((GQL(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /
     >     (((GQR(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GQV(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GQA(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /   
C
     >     (((HQL(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQR(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQV(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQA(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
 
C                        ****************************
      END

      Subroutine StEwCpl2
C                                                   -=-=- stewcpl2

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: EwCpl2
C Description: Setup + Func. Ewk couplings^2 for DY: Scat. and Annih.
C ListOfFiles: stewcpl2 ewcplg2 
C ====================================================================
C Entry points (ewcpl2sc ewcpl2an)

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/EwCpl2.f,v 6.1 98/08/16 17:24:34 wkt Exp $
C #Log:	EwCpl2.f,v $
c Revision 6.1  98/08/16  17:24:34  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 

C          Setup EW SQuaRed coupling constants for Vector Boson Production
C          in the Annihilation and Compton channels; put in Common /EwCplC/. 

C                      --  Unpolarized - (Gv**2 + Ga**2)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NFL = 6, NBN = 4)

      COMMON / EWCPLC / CPLANH(-NFL:NFL, NBN, -NFL:NFL), 
     >                  CPLSCT(-NFL:NFL, NBN)

      Data Lprt / 0 /

      Call SetEwk

      NFLT = NFLTOT()

      If (Lprt .ge. 1)   PRINT '(A/)', 
     >  'Squared EW Couplings for Annih. and Compton processes:'
      DO 4 IBN = 1, NBN
         If (Lprt .ge. 1)   PRINT '(A/9A9)', ' ANNIH:'
     >, 'cb','cb','db','ub','G','u','d','s','c'
      DO 5 IP1 = -NFLT, NFLT
         IF (IP1 .EQ. 0) GOTO 5
C                                                                Scattering
C                                                        Final state summed
         CPLSCT (IP1, IBN) = 0.
         DO 6 IP2 = 1, NFLT
           IF  (IP1 .GT. 0) THEN
             CPLV = GEWQT (IP1, IBN, 1, IP2)
             CPLA = GEWQT (IP1, IBN,-1, IP2)
           ELSE
             CPLV = GEWQT (IP2, IBN, 1,-IP1)
             CPLA = GEWQT (IP2, IBN,-1,-IP1)
           ENDIF
           CPL2 = CPLV**2 + CPLA**2
           CPLSCT (IP1, IBN) = CPLSCT (IP1, IBN) + CPL2
    6 CONTINUE
C                                                                Annihilation
         DO 7 IP2 = -NFLT, NFLT
         IF (IP2 .EQ. 0) GOTO 7
C                                  To annihilate, must have quark - antiquark
         IF (SIGN(1,IP2) .EQ. SIGN(1,IP1)) THEN
            CPLANH (IP1, IBN, IP2) = 0.0
            GOTO 7
         ENDIF

         IF (IP1 .GT. 0) THEN
           CPLV = GEWQT (IP1, IBN, 1,-IP2)
           CPLA = GEWQT (IP1, IBN,-1,-IP2)
         ELSE
           CPLV = GEWQT (IP2, IBN, 1,-IP1)
           CPLA = GEWQT (IP2, IBN,-1,-IP1)
         ENDIF

         CPL2 = CPLV**2 + CPLA**2
         CPLANH (IP1, IBN, IP2) = CPL2
    7 CONTINUE

      If (Lprt .ge. 1)
     $  PRINT '(9F8.3)', (CPLANH(IP1, IBN, IP2), IP2=-4,4)
    5 CONTINUE
      If (Lprt .ge. 1)
     $  PRINT '(/A,9F7.3/)', ' Compt.:', (CPLSCT(IP1, IBN), IP1=-4,4)
    4 CONTINUE

      Return
C                        ****************************
      End

      SUBROUTINE StEwCpl
C                                                   -=-=- stewcpl

C 05-07-13 Mis-match of W+/W- corrected in the following lines.
C            QCPLN(2, -1, IPRTN) = CW * MOD(IS+1, 2)
C            QCPLN(3, -1, IPRTN) = CW * MOD(IS, 2)
C Mistake throughout these years! No consequences in CCFR experiments
C since always had W+ +/- W- 

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: EwCpl1
C Description: Function to extract the Ewk couplings for Scattering and Annih.
C ListOfFiles: stewcpl ewcplg 
C ====================================================================
C Entry points (ewcplsc ewcplan)

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/EwCpl1.f,v 6.1 98/08/16 17:24:17 wkt Exp $
C #Log:	EwCpl1.f,v $
c Revision 6.1  98/08/16  17:24:17  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 

C Set up the basic Electro-weak coupling matrices for the Boson-Fermion
C Yukawa coupling term in the Effective Lagrangian

C     8/16/98  wkt : migrated from DisPac over to EwkPac.

C    Really should be rewritten to use the GEwQx coefficents, 
C    especially if Z-exchanges are to be seriously incorported.
 
C    Leave the way it is for now for expediency!

C Use helicity basis for the boson polarization index in this section.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSPN = 2, NGN = 3, NBSN = 4, NPOL = 3)
      PARAMETER (MXQRK = 6,  CH1 = D1/D3, CH2 = D2/D3)

      DIMENSION IP(NSPN, NGN), T3(NSPN), CHQ(NSPN), JBSN(NBSN)

      COMMON / WKCPLG / QCPLN(NBSN, -1:1, -MXQRK:MXQRK)

c ip(1,1) = 1, ip(1,2) = 4, ip(1,3) = 6 - charge 2/3 quarks, u,c,t
c ip(2,1) = 2, ip(2,2) = 3, ip(2,3) = 5 - charge -1/3 quarks, d,s,b

      DATA IP / 1, 2, 4, 3, 6, 5 /
      DATA T3, CHQ / 1, -1, 0.66666667, -0.3333333 /
      
      DATA JBSN / 1, 3, 2, 4 /     ! Swap W+/W- when making quark-anti-quark
                                   ! transformation in Do 60.

      CW = SQRT (2.)
      Call SETEWK 
c Set up electroweak couplings
      DO 10 IPRTN = 1, MXQRK
C                     Longitudinal coupling vanishes for all bosons
         DO 20 IBSN = 1, 4
            QCPLN(IBSN, 0, IPRTN) = 0
            QCPLN(IBSN, 0,-IPRTN) = 0
 20      CONTINUE
C                     Right-handed coupling vanishes for Charged Bosons
         DO 30 IBSN = 2, 3
 30         QCPLN(IBSN, 1, IPRTN) = 0
 10   CONTINUE

      DO 40 IG = 1, NGN
         DO 50 IS = 1, NSPN
            IPRTN = IP(IS, IG)
C                                  Non-vanishing Photon couplings
            QCPLN(1, -1, IPRTN) = CHQ(IS)
            QCPLN(1,  1, IPRTN) = CHQ(IS)
C                                 Non-vanishing Charged boson couplings
            QCPLN(2, -1, IPRTN) = CW * MOD(IS+1, 2)
            QCPLN(3, -1, IPRTN) = CW * MOD(IS, 2)
C                                 Neutral Z couplings
            QCPLN(4, -1, IPRTN) = T3(IS) - 2.*CHQ(IS)*SWG2F()
            QCPLN(4,  1, IPRTN) = -2.*CHQ(IS)*SWG2F()
c swg2f returns sin^2(theta_weinberg)
C                                 Anti-quark couplings - by CP
            DO 60 ib = 1, 4
c iterate over ivb
               QCPLN(JBSN(IB), -1, -IPRTN) = QCPLN(IB,  1, IPRTN)
               QCPLN(JBSN(IB),  1, -IPRTN) = QCPLN(IB, -1, IPRTN)
 60            CONTINUE

 50      CONTINUE
 40   CONTINUE

      RETURN
C                        ****************************
      END

