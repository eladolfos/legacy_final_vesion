      FUNCTION AMOMEN (AM, NX, F)

C       Given the function F(x) defined in the interval (0,1] on the array
C       (1, NX), evenly distributed in the variable z where x = xfrmz (Z)
C       this routine returns its Amth moment defined as Int F(x) (x**Am) dx/x,
C
C       The value of F(x) as x --> 0 is obtained by extrapolation in the
C       variable z from the first three points F1, F2, and F3.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 204)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C


      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
C
      DIMENSION F(NX), HH(NDH, MXX), A(NDG), DI(0:3)
C
      IF (NX .LE. 10 .OR. NX .GT. MXX) THEN
        CALL WARNI(IWRN,NWRT,'NX out of range in AMOMEN subroutine ',
     >       'NX', NX, 10, MXX, 1)
      END IF
C                                                           Compute "A" matrix
      DO 30 I = NX-1, 1, -1
         DO 40 K = 1, NDG
           AK = AM - 2D0 + K
           IF (AK .EQ. 0D0) THEN
             A(K) = XA(I+1, 0) - XA(I, 0)
           Else
             A(K) = (XA(I+1,1)**AK - XA(I,1)**AK) / AK
           EndIf
   40    CONTINUE
C                                                            Compute "H" matrix
         DO 41 J = 1, NDH
           TEM = 0
           DO 43 L = 1, NDG
               TEM = TEM + A(L) * GB(L,J,I)
   43      CONTINUE
           HH(J,I) = TEM
   41    CONTINUE
   30  CONTINUE
C                                 Calculate integral from X= X(I=4) to 1.(I=Nx)
C                                                       Patterned after INTEGR
      TEM1 = 0D0
      DO 60 I = NX-1, 4, -1
         TEM1 = TEM1 + HH(1,I)*F(I-1) + HH(2,I)*F(I)
     >               + HH(3,I)*F(I+1) + HH(4,I)*F(I+2)
   60 CONTINUE
C                              Calculate contribution from the first 3 bins and
C                                 Use the first 3 bins to extrapolate to X = 0.
      TEM2 = 0
      DO 61 I = 1, 3
         DI(I) =  HH(1,I)*F(I-1) + HH(2,I)*F(I)
     >          + HH(3,I)*F(I+1) + HH(4,I)*F(I+2)
         TEM2 = TEM2 + DI(I)
   61 CONTINUE
C                                                 Quadratic extrapolation gives
      DI(0) = 3D0 * (DI(1) - DI(2)) + DI(3)
C
      TEM = TEM1 + TEM2 + DI(0)
C
      AMOMEN = TEM
      RETURN
C                        ****************************
      END

      FUNCTION C2QX (XX)
C
C Integrands for the convolution integral in DIS to MS-bar scheme
C  transformation

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1=1D0,D2=2D0,D3=3D0,D4=4D0,D5=5D0,D6=6D0,D10=1D1)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / MSTDS1 / X, Q, NF, JSET, JHDRN, JPRTN, NCNT

      DATA
     1 D0, DUM / 0.0, 0.0 /
     1 IW1, IW2, IW3 / 3 * 0 /
C                                                          Statement functions
      QRK(Y) = PDF (JSET, JHDRN, JPRTN, Y, Q, IRT)
      GLU(Y) = PDF (JSET, JHDRN,     0, Y, Q, IRT)
C
      C2Q(X) = (2./3.) * (1.+ X**2) / (1.- X) * LOG(X/(1.-X))
     >             + 1./ (1.- X) - 2.- (4./3.) * X
C
      C2G(X) = (1./4) * (X**2 + (1.- X)**2) * (LOG( X/(1.-X) ) + 1.)
     >             - (3./2.) * X * (1.-X)
C
      XV = XX
      IF (XV .LE. D0) THEN
         PRINT *, 'X < 0 out of range in C2QX, X = ', XV
         TEM = DUM
      ElseIF (XV .LT. D1) THEN
         TEM = C2Q(XV)
      Else
         PRINT *, 'X > 1 out of range in C2QX, X =', XV
         TEM = DUM
      EndIf

      C2QX = TEM
      RETURN
C                                                          --------------------
      ENTRY QQY (YY)
C                                                          Quark to Quark piece
      Y = YY
      IF (Y .GT. X .AND. Y .LT. D1) THEN
        TEM = C2Q (Y) * (QRK(X/Y) - QRK(X)*Y) / Y
      ElseIF (Y .EQ. X) THEN
        TEM = C2Q (Y) * (-QRK(X))
      Else
        TEM = 0.
        IF (Y .NE. D1)
     >  CALL WARNR(IW1, NWRT, 'Y out of range in QQY', 'Y',Y,X,D1,I1)
      EndIf
C
      QQY = TEM
C
      NCNT = NCNT + 1
      RETURN
C                         ----------------------------
      ENTRY GQY (YY)
C                                                      Singlet quark to gluon
      Y = YY
      IF (Y .GT. X .AND. Y .LT. D1) THEN
        TEM = C2Q (Y) * (SQRK(X/Y) - Y *SQRK(X)) / Y
      ElseIF (Y .EQ. X) THEN
        TEM = C2Q (Y) * (-SQRK(X))
      ElseIF (Y .EQ. D1) THEN
        TEM = 0.
      Else
        CALL WARNR (IW2, NWRT, 'Y out of range in GQY', 'Y',Y,X,D1,I1)
        TEM = 0.
      EndIf
C
      GQY = TEM
C
      NCNT = NCNT + 1
      RETURN
C                        ----------------------------
      ENTRY QGY (YY)
C
      Y = YY
      IF (Y .GT. X .AND. Y .LT. D1) THEN
        TEM = C2G (Y) * GLU(X/Y) / Y
      Else
        TEM = 0.
        IF (Y .LT. X .OR. Y .GT. D1)
     >  CALL WARNR(IW3, NWRT, 'Y out of range in QGY', 'Y',Y,X,D1,I1)
      EndIf
C
      QGY = TEM
C
      NCNT = NCNT + 1
      RETURN
C                        ****************************
      END

      SUBROUTINE DELDIS (XX, DFF2, DFG2, DGF2, DGG2, NFL, IRT)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      PARAMETER (D0 = 0.0, D1 = 1.0)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      
      external DelFF2_pn, DelFG2_pn, DelGF2_pn, DelGG2_pn

C
      DATA IWRN /0/
C
      IRT = 0
C
      X = XX
      IF (X .LE. 0. .OR. X .GE. 1.) THEN
        CALL WARNR(IWRN, NWRT, 'X out of range in DELDIS', 'X', X,
     >             D0, D1, 1)
        IRT = 1
        RETURN
      EndIf

cpn 2001
      if (iknl.eq.12) then
        DFF2 = DelFF2_pn(x)
        DFG2 = DelFG2_pn(x)
        DGF2 = DelGF2_pn(x)
        DGG2 = DelGG2_pn(x)
      else
        print *,'DelDIS cannot be called with iknl=',iknl
        stop
      endif                     !iknl
C
C                                       Take out singular factors at both ends
      XLG = (LOG(1./(1.-X)) + 1.)
      XG2 = XLG ** 2

      DFF2 = DFF2 * X * (1.- X)
      DFG2 = DFG2 * X / XG2
      DGF2 = DGF2 * X *(1.-x)
      DGG2 = DGG2 * X / XG2

      RETURN
C                        ****************************
      END!DelDIS


      double precision function beta_c2q(x)
cpn                    Returns the coefficient for transformation of the
cpn                    one-loop non-singlet kernels from the MS-bar
cpn                     scheme to DIs scheme
      implicit double precision (A-H,O-Z)
      double precision x
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv
      
      tem = (1.+x**2)/(1.-x)*(dlog((1.-x)/x)-0.75)+(9.+5.*x)/4.
      tem = 0.5*(11./3.*CfCg-4./3.*CfTr*Nflv)*tem
      
      beta_c2q = tem
      return
      end !beta_c2q

cpn---------------------------------------------------
      double precision function DelFF2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      DelFF2CFTR = 
     >        (-18.666666666666668 + Pi2/3. + 31*x + Pi2*x - 20*x2 - 
     >    (8*Pi2*x2)/3. + (16*x3)/3. + (4*Pi2*x3)/3. + (4*xi)/3. + 
     >    6*(-1 + x)*xLi2 + ((-11 - 45*x + 96*x2 - 44*x3)*xln)/3. + 
     >    (-7 + 19*x - 28*x2 + (44*x3)/3. + (8*xi)/3.)*xln1m + 
     >    4*(-1 + 3*x - 4*x2 + 2*x3)*xln*xln1m + 
     >  (2 - 6*x + 8*x2 - 4*x3)*xln1m2 + (-1 - 3*x + 8*x2 - 4*x3)*xln2)/
     >  (1 - x)

      DelFF2CFCG = 
     > -11*((9 + 5*x)/4. + (1 + x2)*(-0.75 + Log((1 - x)/x))/(1 - x))/6.

      DelFF2_pn = CfTr*Nflv*DelFF2CFTR + CfCg*DelFF2CFCG

      Return
      End!DelFF2_pn

cpn---------------------------------------------------
      double precision function DelFG2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      DelFG2CFTR = 
     >  -12 + (2*Pi2)/3. + 22*x - (4*Pi2*x)/3. - 16*x2 + (4*Pi2*x2)/3. + 
     >  (-4 + 8*x - 8*x2)*xLi2 + (-2 + 8*x - 20*x2)*xln + 
     >  (-2 - 20*x + 20*x2)*xln1m + (-2 + 4*x - 4*x2)*xln1m2

      DelFG2CGTR = 
     >  -14.333333333333334 - (242*x)/3. + 4*Pi2*x + (281*x2)/3. - 
     >  (4*Pi2*x2)/3. + (4*xi)/3. - 4*(1 + 4*x)*xLi2 + 
     >  (-2 - 64*x + (62*x2)/3.)*xln + 
     >  (2*(-3 + 72*x - 79*x2 + 4*xi)*xln1m)/3. + 
     >  (-4 + 8*x - 8*x2)*xln*xln1m + (4 - 8*x + 8*x2)*xln1m2 + 
     >  (-2 - 8*x)*xln2

      DelFG2TR2 = 
     >  22 - (2*Pi2)/3. + 96*x - (8*Pi2*x)/3. - 118*x2 - (8*Pi2*x2)/3. + 
     >  4*(1 + 2*x)**2*xLi2 + 8*(1 + 8*x + 5*x2)*xln + 
     >  8*(-1 - 2*x + 3*x2)*xln1m + 2*(1 + 2*x)**2*xln2

      DelFG2_pn = CfTr*Nflv*DelFG2CFTR + CgTr*Nflv*DelFG2CGTR + 
     >  (Tr*Nflv)**2*DelFG2TR2

      Return
      End!DelFG2_pn

cpn---------------------------------------------------
      double precision function DelGF2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      DelGF2CF2 = 
     >        (2.25 + Pi2/6. + 6*x + (Pi2*x)/2. + 3*x2 + (2*Pi2*x2)/3. + 
     >    (-5 + 3*x - 2*x2 + 4*xi)*xLi2 + ((4 + x - 5*x2)*xln)/2. + 
     >    (-15 + (19*x)/2. + (5*x2)/2. + 6*xi)*xln1m + 
     >    (-4 + 6*x + 2*x2 + 4*xi)*xln*xln1m + 
     >    (1 - 3*x - 2*x2 - 2*xi)*xln1m2 + ((1 - 3*x - 2*x2)*xln2)/2.)/
     >  (1 - x)

      DelGF2CFTR = 
     >    (12.666666666666666 - (2*Pi2)/3. - (40*x)/3. - (20*x2)/3. + 
     >    (2*Pi2*x2)/3. + (32*x3)/3. - (4*xi)/3. + (4 - 4*x2)*xLi2 + 
     >    (2*(5 + 12*x - 17*x2 + 4*x3)*xln)/3. - 
     >   (2*(1 + x)**2*(4 - 7*x + 4*x2)*xi*xln1m)/3. + (2 - 2*x2)*xln2)/
     >  (1 - x)

      DelGF2CFCG = 
     >     (3.5 - (50*x)/3. - Pi2*x + (20*x2)/3. + (Pi2*x2)/3. - 8*x3 - 
     >    (2*Pi2*x3)/3. + (6 - 6*x + 4*x2 - 4*xi)*xLi2 + 
     >    (-3.6666666666666665 + 3*x - (29*x2)/3. + 6*x3)*xln + 
     >    (20.666666666666668 - 21*x + (41*x2)/3. - 6*x3 - 6*xi)*xln1m + 
     >    (6 - 12*x + 6*x2 - 4*x3 - 4*xi)*xln*xln1m + 
     >  2*(-1 +3*x-x2 + x3 + xi)*xln1m2 +(-1 + 3*x - 2*x2 + 2*x3)*xln2)/
     >  (1 - x)

      DelGF2_pn = Cf2*DelGF2CF2 + CfTr*Nflv*DelGF2CFTR + CfCg*DelGF2CFCG

      Return
      End!DelGF2_pn

cpn---------------------------------------------------
      double precision function DelGG2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      DelGG2CFTR = 
     >  18.333333333333332 - Pi2/3. - 40*x/3. - 4*Pi2*x/3. + 16*x2/3. + 
     >  (4*Pi2*x2)/3. - (4*xi)/3. + 6*xLi2 + (3 + 18*x - 44*x2/3.)*xln + 
     >  (5 - 14*x + (44*x2)/3. - (8*xi)/3.)*xln1m + 
     >  (4 - 8*x + 8*x2)*xln*xln1m + (-2 + 4*x - 4*x2)*xln1m2 + 
     >  (1 + 4*x - 4*x2)*xln2

      DelGG2CGTR = 
     > (11*(-1 + 8*(1 - x)*x + ((1 - x)**2 + x2)*Log((1 - x)*xi)))/3.

      DelGG2TR2 = 
     >        1.3333333333333333 - (32*x)/3. + (32*x2)/3. - 
     >  (4*Log((1 - x)*xi))/3. + (8*x*Log((1 - x)*xi))/3. - 
     >  (8*x2*Log((1 - x)*xi))/3.

      DelGG2_pn = CfTr*Nflv*DelGG2CFTR + CgTr*Nflv*DelGG2CGTR + 
     >  (Tr*Nflv)**2*DelGG2TR2

      Return
      End!DelGG2_pn


cpn 2001 ---------------------------------------------------------
      double precision function RDFF2_pn(x)
      implicit NONE
      double precision DelFF2_pn, x
      
      RDFF2_pn = x* DelFF2_pn(x)
      return 
      End ! RDFF2_pn

cpn 2001 ---------------------------------------------------------
      double precision function RDGF2_pn(x)
      implicit NONE
      double precision DelGF2_pn, x
      
      RDGF2_pn = x* DelGF2_pn(x)
      return 
      End ! RDGF2_pn


cpn Exact expressions for some of the antiderivatives. They can be used
cpn to check the correctness of the numerical integration of these 
cpn antiderivatives in the subroutine STUPKL.
cpn---------------------------------------------------
      double precision function WFG2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv


      XI = 1./ X                !auxiliary definitions
      X2 = X ** 2
      X3=  x**3
      XLN = DLOG (X)
      XLN2 = XLN ** 2
      XLN1M = DLOG (1.- X)
         
      xLi2m=xLi(2,-x)
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1d0-x)
      xLi32=xLi(3,x2)
      xln1m2=xln1m*xln1m
      xln1p=dlog(1d0+x)
      x1m=1d0-x
      x1p=1d0+x
      x3m=3d0-x
      x3p=3d0+x
      
      wfgcft=
     >  (18 - 81*x + 6*Pi2*x + 123*x2 - 6*Pi2*x2 - 60*x3 +
     >  4*Pi2*x3 - 6*(-2 + 3*x - 3*x2 + 2*x3)*xln1m2 -33*x*xln +
     >  15*x2*xln - 24*x3*xln - 9*x*xln2 + 9*x2*xln2 -
     >  12*x3*xln2 - 12*x1m*xln1m*(-1 + 2*x2 + 2*xln - x*xln +
     >  2*x2*xln) - 24*xLi2)/9.
      
      wfgcgt=
     >  (2*(-67 + 2*Pi2 + x*(64 + x*(-91 + 3*Pi2 + 94*x)) +
     >  x1m*(7+x*(-5+16*x))*xln1m -3*x1m*(2+ x*(-1+2*x))*xln1m2 -
     >  20*xln - 3*x*xln*(13 + 16*x*x1p - 3*x1p*xln) +
     >  6*x1p*(2+x+2*x2)*xln*xln1p+6*x1p*(2+x+2*x2)*xLi2m))/9.
      
      WFG2_pn  = CfTR*Nflv*WFGCFT           + CgTr*Nflv * WFGCGT
      
      return 
      End !WFG2_pn
      
      double precision function WGF2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/ ! zeta(3.0)

      XI = 1./ X                !auxiliary definitions
      X2 = X ** 2
      X3=  x**3
      XLN = DLOG (X)
      XLN2 = XLN ** 2
      XLN1M = DLOG (1.- X)
         
      xLi2m=xLi(2,-x)
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1d0-x)
      xLi32=xLi(3,x2)
      xln1m2=xln1m*xln1m
      xln1p=dlog(1d0+x)
      x1m=1d0-x
      x1p=1d0+x
      x3m=3d0-x
      x3p=3d0+x
      
      wgfcft=
     >  (9 + 4*Pi2 - 22*x + 13*x2 + 6*(3 - 4*x + x2)*xln1m +
     >  40*xln - 24*xLi2)/9.

      wgfcf2=
     >  (6*(2*(-9 + Pi2) + 3*x*(5 + x)) +4*(3 +2*Pi2+3*x*(-3 + 2*x))*
     >  xln1m + 6*x3m*x1m*xln1m2 - 6*(x*(8 + 3*x) + 4*xln1m2)*
     >  xln - 3*(-4 + x)*x*xln2)/12 - 2*(3 + 2*xln1m)*xLi2 - 4*xLi31m

      wgfcfg=
     >  (3637-186*Pi2-x*(3198+72*Pi2+x*(231 + 208*x)))/108.- xln +
     >  (3*xln1m*(-33 - 4*Pi2 + (50 - 17*x)*x - 3*x3m*x1m*xln1m) +
     >  2*(x*(198 + x*(27+8*x))+9*xln1m*(3 - 4*x + x2 + 2*xln1m))*
     >  xln - 9*x*(4 + x)*xln2)/18- x1p*x3p*xln*xln1p-
     >  (x1p*x3p - 4*xln)*xLi2m + (31d0/3d0 +4*xln1m- 4*xln)*xLi2 +
     >  4*xLi31m + 12*xLi3 - 2*xLi32 - 10*zeta3


      
      WGF2_pn = Cf*Tr*WGFCFT +CF2* WGFCF2 + CfCg *WGFCFG     
      
      return 
      End !WFG2_pn
 

cpn---------------------------------------------------
      double precision function WDFF2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      WDFF2CFTR = 
     >  32*x/9. - 4*x**4/3. - (Pi2*x**4)/3. - (56*x2)/9. + (Pi2*x2)/6. + 
     >  (59*x3)/9. + (4*Pi2*x3)/9. + (2.6666666666666665 - 3*x2)*xLi2 + 
     >  (x*(24 - 46*x2 + 33*x3)*xln)/9. + 
     >  ((44 - 33*x**4 - 36*x2 + 34*x3)*xln1m)/9. + 
     >  (2.6666666666666665 - 2*x**4 - 2*x2 + (8*x3)/3.)*xln*xln1m + 
     >  (-1.3333333333333333 + x**4 + x2 - (4*x3)/3.)*xln1m2 + 
     >  (x2*(-3 - 8*x + 6*x2)*xln2)/6.

      WDFF2CFCG = 
     >        (-77*x)/18. - (55*x2)/18. - (11*x3)/9. - (11*xLi2)/3. - 
     >  (11*x*(12 + 3*x + 2*x2)*xln)/36. + 
     >  (11*(-26 + 12*x + 3*x2 + 2*x3)*xln1m)/36. - (11*xln*xln1m)/3. + 
     >  (11*xln1m2)/6.

      WDFF2_pn = CfTr*Nflv*WDFF2CFTR + CfCg*WDFF2CFCG

      Return
      End!WDFF2_pn

cpn---------------------------------------------------
      double precision function WDFG2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      WDFG2CFTR = 
     >   (14*x)/3. - (2*Pi2*x)/3. - 10*x2 + (2*Pi2*x2)/3. + (16*x3)/3. - 
     >  (4*Pi2*x3)/9. + (4*x*(3 - 3*x + 2*x2)*xLi2)/3. + 
     >  (2*x*(3 - 6*x + 10*x2)*xln)/3. - 
     >  (2*(8 - 5*x - 13*x2 + 10*x3)*xln1m)/3. + 
     >  (2*(-2 + 3*x - 3*x2 + 2*x3)*xln1m2)/3.

      WDFG2CGTR = 
     >  4.222222222222222 - (4*Pi2)/9. + (113*x)/9. + (238*x2)/9. - 
     >  2*Pi2*x2 - (313*x3)/9. + (4*Pi2*x3)/9. + 4*x*(1 + 2*x)*xLi2 - 
     >  (2*(6 + 21*x - 129*x2 + 35*x3)*xln)/9. + 
     >  (2*(-29 + 33*x - 87*x2 + 83*x3)*xln1m)/9. + 
     >  (4*(-1 + x)*(2 - x + 2*x2)*xln*xln1m)/3. + 
     >  (2.6666666666666665 - 4*x + 4*x2 - (8*x3)/3.)*xln1m2 + 
     >  2*x*(1 + 2*x)*xln2

      WDFG2TR2 = 
     >    0.8888888888888888 - (146*x)/9. + (2*Pi2*x)/3. - (280*x2)/9. + 
     >  4*Pi2*x2/3. + (418*x3)/9. + (8*Pi2*x3)/9. - 
     >  4*x*(3 + 6*x + 4*x2)*xLi2/3. - (4*x*(9 + 63*x + 26*x2)*xln)/9. + 
     >  4*(4 + 9*x + 9*x2 - 22*x3)*xln1m/9.-2*x*(3 + 6*x + 4*x2)*xln2/3.

      WDFG2_pn = CfTr*Nflv*WDFG2CFTR + CgTr*Nflv*WDFG2CGTR + 
     >  (Tr*Nflv)**2*WDFG2TR2

      Return
      End!WDFG2_pn

cpn---------------------------------------------------
      double precision function WDGF2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      WDGF2CF2 = 
     >   (-13*x)/3. - (4*Pi2*x)/3. - (17*x2)/12. - (7*Pi2*x2)/12. - x3 - 
     >  2*Pi2*x3/9. - 4*xLi3 - 8*xLi31m +x*(56 + 25*x + 10*x2)*xln/12. + 
     >  xLi2*(8.666666666666666 + 4*x - x2/2. + (2*x3)/3. + 4*xln - 
     >     8*xln1m) + ((4 - 40*x - 89*x2 - 10*x3)*xln1m)/12. - 
     >  (2*(-13 + 6*x + 6*x2 + x3)*xln*xln1m)/3. + 2*xln1m**3 + 
     >  (-8.666666666666666 + 4*x + (5*x2)/2. + 2*x3/3. -8*xln)*xln1m2 + 
     >  (2*x + (5*x2)/4. + x3/3. + 2*xln1m)*xln2 + 8*zeta3

      WDGF2CFTR = 
     >  -32*x/9. - (8*x**4)/3. + (32*x2)/9. - (Pi2*x2)/3. - (26*x3)/9. - 
     >  (2*Pi2*x3)/9. + (2*(-4 + 3*x2 + 2*x3)*xLi2)/3. - 
     >  (2*x*(12 + 3*x - 11*x2 + 3*x3)*xln)/9. + 
     >  (2*(-22 + 3*x**4 + 3*x2 + 7*x3)*xln1m)/9. - (8*xln*xln1m)/3. + 
     >  (4*xln1m2)/3. + (x2 + (2*x3)/3.)*xln2

      WDGF2CFCG = 
     >  155*x/18. + 4*Pi2*x/3. + 2*x**4 + (Pi2*x**4)/6. + (209*x2)/36. + 
     >  (2*Pi2*x2)/3. + (7*x3)/18. + (Pi2*x3)/9. + 4*xLi3 + 8*xLi31m - 
     >  (x*(18 + 15*x - 20*x2 + 27*x3)*xln)/18. + 
     >  ((137 - 6*x + 27*x**4 + 147*x2 - 44*x3)*xln1m)/18. + 
     >  (-5 + 4*x + x**4 + 5*x2 - (2*x3)/3.)*xln*xln1m - 2*xln1m**3 + 
     >  xLi2*(-5 - 4*x + x2 - (4*x3)/3. - 4*xln + 8*xln1m) + 
     >  (6.833333333333333 - 4*x - x**4/2. - 3*x2 + 8*xln)*xln1m2 + 
     >  ((-4*x - x**4 - 3*x2 - 4*xln1m)*xln2)/2. - 8*zeta3

      WDGF2_pn = Cf2*WDGF2CF2 + CfTr*Nflv*WDGF2CFTR + CfCg*WDGF2CFCG

      Return
      End!WDGF2_pn

cpn---------------------------------------------------
      double precision function WDGG2_pn(x)
      implicit double precision (A-H,O-Z)
      PARAMETER (PI=3.141592653589793,PI2=PI**2)
      external xLi
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg 
      common/nflav/nflv

      data zeta3/1.20205690315959/!zeta(3.0)

      xi=1./x
      x2=x**2
      x3=x**3
      xln=log (x)
      xln2=xln**2
      xln1m=log(1.-x)
      xln1m2=xln1m**2
      xLi2=xLi(2,x)
      xLi3=xLi(3,x)
      xLi31m=xLi(3,1-x)

      WDGG2CFTR = 
     >  2.7777777777777777 + (4*Pi2)/9. - (100*x)/9. + (Pi2*x)/3. + 
     >  91*x2/9. + 2*Pi2*x2/3. - (16*x3)/9. - (4*Pi2*x3)/9. - 6*x*xLi2 + 
     >  ((12 + 15*x - 69*x2 + 44*x3)*xln)/9. + 
     >  ((80 - 87*x + 51*x2 - 44*x3)*xln1m)/9. - 
     >  (4*(-1 + x)*(2 - x + 2*x2)*xln*xln1m)/3. + 
     >  2*(-2 + 3*x - 3*x2 + 2*x3)*xln1m2/3. + x*(-3-6*x + 4*x2)*xln2/3.

      WDGG2CGTR = 
     >  1.2222222222222223 + (22*x)/9. - (121*x2)/9. + (88*x3)/9. + 
     >  (11*x*(3 - 3*x + 2*x2)*xln)/9. - 
     >  (11*(-2 + 3*x - 3*x2 + 2*x3)*xln1m)/9.

      WDGG2TR2 = 
     >   -0.4444444444444444 - (8*x)/9. + (44*x2)/9. - (32*x3)/9. - 
     >  4*x*(3 - 3*x + 2*x2)*xln/9. + 4*(-2 + 3*x- 3*x2 + 2*x3)*xln1m/9.

      WDGG2_pn = CfTr*Nflv*WDGG2CFTR + CgTr*Nflv*WDGG2CGTR + 
     >    (Tr*Nflv)**2*WDGG2TR2

      Return
      End!WDGG2_pn
      function dFnsm_pn(x)
cpn                    Returns P_{NS,+}^{(1)} (conventional notations}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision dFnsm_pn, x
      double precision PCF2,PCFG,PCFT,PQQB,PQQ2,FFP,FFM,x2,xln1m,xln1m2,
     >   xln,xln2,spen2,spenc2,tem,xi,xm1i,xp1i

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PNSP: x=',x,' is out of range'

      XI = 1./ X
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      FFP = (1.+ X2) * XM1I
      FFM = (1.+ X2) * XP1I
 
      PCF2 = -2.* FFP *XLN*XLN1M - (3.*XM1I + 2.*X)*XLN
     >     - (1.+X)/2.*XLN2 - 5.*(1.-X)
      PCFG = FFP * (XLN2 + 11.*XLN/3.)+ (67./9.- PI**2 / 3.)*(-1.-x)
     >     + 2.*(1.+X) * XLN + 40.* (1.-X) / 3.
      PCFT = (-FFP *XLN + 5./3.*(1.-x) - 2.*(1.-X)) * 2./ 3.
      
      PQQB = 2.* FFM * SPEN2 + 2.*(1.+X)*XLN + 4.*(1.-X)
      PQQB = (CF2-CFCG/2.) * PQQB
      PQQ2 = CF2 * PCF2 + CFCG * PCFG / 2. + CFTR*NFlv * PCFT

      tem = PQQ2 + PQQB

      dFnsm_pn=tem

      Return
      End!dFnsm_pn
cpn---------------------------------------------------------------------------

      function dPFF1_pn(x)
cpn                                                Returns Delta P_{FF}^{(0)} 
      implicit NONE
      double precision dPFF1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PFF1: x=',x,' is out of range'

      dPFF1_pn=CF*(1.+x**2)/(1.-x)  
      Return
      End!PFF1_pn

cpn-------------------------------------------------------------------------


      function dPFF2_pn(x)
cpn                    Returns Delta P_{FF}^{(1)}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision dPFF2_pn, x
      double precision  FFCF2, FFCFG,FFCFT,FFP,FFM,x2,xln1m,xln1m2,xln,
     >   xln2,spen2,spenc2,tem,xi,xm1i,xp1i

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PFF2: x=',x,' is out of range'

      XI = 1./ X
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      FFP = (1.+ X2) * XM1I
      FFM = (1.+ X2) * XP1I
 
      FFCF2 = -2.* FFP *XLN*XLN1M - (3.*XM1I + 2.*X)*XLN
     >     - (1.+X)/2.*XLN2 - 5.*(1.-X)-
     >     (2.* FFM * SPEN2 + 2.*(1.+X)*XLN + 4.*(1.-X))
      FFCFG = FFP * (XLN2 + 11.*XLN/3.)- (67./9.- PI**2 / 3.)*(1+x)
     >     + 2.*(1.+X) * XLN + 40.* (1.-X) / 3.+
     >     (2.* FFM * SPEN2 + 2.*(1.+X)*XLN + 4.*(1.-X))
      FFCFG=0.5*FFCFG       !to account for extra 1/2 in the factor of NS piece
      
      FFCFT = (-FFP * XLN + 5./3.*(1+x) - 2.*(1.-X)) * 2./ 3.+
     >     2d0*(1d0-x-(1d0-3d0*x)*xln-(1d0+x)*xln2)
      
      tem = CFTR*NFlv * FFCFT + CF2 * FFCF2 + CFCG   * FFCFG

      dPFF2_pn=tem

      Return
      End!dPFF2_pn

cpn---------------------------------------------------------------------


      function dPFG1_pn(x)
cpn                      Returns Delta P_{FG}^{(0)}(conventional notations)
      implicit NONE
      double precision dPFG1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PFG1: x=',x,' is out of range'

      dPFG1_pn=2*Tr*Nflv *(2.* X -1.)

      Return
      End!dPGF1_pn

cpn-------------------------------------------------------------------------


      function dPFG2_pn(x)
cpn                    Returns P_{FG}^{(2)} (conventional notations}

      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision dPFG2_pn, x
      double precision  FGCFT,FGCGT,x2,xln1m,xln1m2,xln,xln2,spen2
     >     ,spenc2,tem

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGF2: x=',x,' is out of range'

      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

 
      FGCFT =-22.+27.*x -9.*xln+8*(1.-x)*xln1m
     >     +(2.*x-1.)*(2.*xln1m2-4.*xln*xln1m+xln2-2./3.*Pi2)
      
      FGCGT = 2.*(12.-11.*x) -8.*(1.-x)*xln1m+2.*(1.+8.*x)*xln
     >     -2.*(xln1m2-Pi2/6.)*(2.*x-1.)-
     >     (2.*Spen2-3.*xln2)*(-2.*x-1)

      tem = CFTR*NFlv * FGCFT               + CGTR*nflv * FGCGT
      dPFG2_pn=tem

      Return
      End!dPFG2_pn

cpn----------------------------------------------------------------------------
      function dPGF1_pn(x)
cpn                       Returns Delta P_{GF}^{(0)}(conventional notations)
      implicit NONE
      double precision dPGF1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg


      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGF1: x=',x,' is out of range'

      dPGF1_pn=CF*(2.- X)
      
      Return
      End!dPGF1_pn

cpn-------------------------------------------------------------------------

      function dPGF2_pn(x)
cpn                    Returns P_{GF}^{(1)} (conventional notations}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision dPGF2_pn, x
      double precision  GFCF2, GFCFG,GFCFT,x2,xln1m,xln1m2,xln,xln2,
     >  spen2,spenc2,tem

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGF2: x=',x,' is out of range'


      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)
     
      GFCF2 =-0.5-0.5*(4.-x)*xln-(2.+x)*xln1m
     >     + (-4.-xln1m2+0.5*xln2)*(2.-x)

      GFCFG =(4.-13.*x)*xln+(10.+x)*xln1m/3.+(41.+35.*x)/9.
     >     +0.5*(-2.*spen2+3.*xln2)*(2.+x)
     >     +(xln1m2-2*xln1m*xln-Pi2/6.)*(2.-x)

      GFCFT =-4./9.*(x+4.)-4./3.*(2.-x)*xln1m

      tem= CFTR*NFlv * GFCFT + CF2 * GFCF2 + CFCG   * GFCFG
      dPGF2_pn=tem

      Return
      End!dPGF2_pn

cpn----------------------------------------------------------------------------
      function dPGG1_pn(x)
cpn                                                Returns Delta P_{GG}^{(0)} 
      implicit NONE
      double precision dPGG1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGG1: x=',x,' is out of range'

      dPGG1_pn=2.*CG *(1./(1.-X) - 2./X +1.)

      Return
      End!dPGG1_pn

cpn---------------------------------------------------------------------------

   
      function dPGG2_pn(x)
cpn                    Returns P_{GG}^{(2)}(x) 
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision dPGG2_pn, x
      double precision  GGCFT,GGCGT,GGCG2,GGP,x2,xln1m,xln1m2,xln,xln2,
     >   spen2,spenc2,tem,xm1i,xp1i,xi

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGG2: x=',x,' is out of range'

      XI = 1./ X
      X2 = X ** 2
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      GGP=xm1i-2.*x+1

      GGCGT =-(4.*(1.-x)+4./3.*(1.+x)*xln+20./9.*(-2.*x+1.))
      
      GGCFT =-(10.*(1.-x)+2.*(5.-x)*xln+2.*(1.+x)*xln2)

      GGCG2 =(29.-67.*x)*xln/3.-19./2.*(1.-x)+4.*(1.+x)*xln2
     >     -2.*Spen2*(xp1i+2.*x+1) 
     >     +(-4.*xln1m*xln+xln2)*GGP +(67./9.-Pi2/3.)*(-2.*x+1.)
 
      tem = CFTR*NFlv * GGCFT + CG2 * GGCG2 + CGTR*NFlv * GGCGT

      dPGG2_pn=tem

      Return
      End!dPGG2_pn

cpn 2001 -------------------------------------------------------
cpn  Functions for the conversion of 2nd order MS-bar splitting
c    kernels into DIS splitting kernels

cpn 2001 A subroutine to calculate transformation coefficients for
cpn      2nd order splitting functions from the MS-bar to DIS scheme
      FUNCTION DXDZ (Z)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA HUGE, IWRN / 1.E20, 0 /
C
      ZZ = Z
      X = XFRMZ (ZZ)
C
      TEM = DZDX (X)
C
      IF     (TEM .NE. D0) THEN
        TMP = D1 / TEM
      Else
      CALL WARNR(IWRN, NWRT, 'DXDZ is singular in DXDZ; DXDZ set=HUGE',
     >             'Z', Z, D0, D1, 0)
        TMP = HUGE
      EndIf
C
      DXDZ = TMP
      RETURN
C                        ****************************
      END
      SUBROUTINE EVLGET (NAME, VALUE, IRET)
C                                                   -=-=- evlget
C
C                                           Gets VALUE of variable named NAME.
C
C               IRET =   0      variable not found.
C                        1      success.
C
C               NAME is assumed upper-case, and VALUE real.
C               If necessary, VALUE is converted to integer by NINT
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      CHARACTER*(*) NAME
      Character Flin*10
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
     > / PdfSwh / Iset, IpdMod, Iptn0, NuIni
C
      IRET = 1
C
      IF     (NAME .EQ. 'QINI')  THEN
          VALUE = QINI
      ElseIF (NAME .EQ. 'IPD0')  THEN
          VALUE = IPD0
      ElseIF (NAME .EQ. 'IHDN') THEN
          VALUE = IHDN
      ElseIF (NAME .EQ. 'QMAX')  THEN
          VALUE = QMAX
      ElseIF (NAME .EQ. 'IKNL') THEN
          VALUE = IKNL
      ElseIF (NAME .EQ. 'XCR') THEN
          VALUE = XCR
      ElseIF (NAME .EQ. 'XMIN') THEN
          VALUE = XMIN
      ElseIF (NAME .EQ. 'NX') THEN
          VALUE = NX
      ElseIF (NAME .EQ. 'NT') THEN
          VALUE = NT
      ElseIF (NAME .EQ. 'JT') THEN
          VALUE = JT
      ElseIF (NAME .EQ. 'NFMX') THEN
          VALUE = NfMx
      ElseIF (NAME .EQ. 'IPDMOD') THEN
          VALUE = IpdMod
      ElseIF (NAME .EQ. 'IPTN0') THEN
          VALUE = IPTN0
      ElseIF (NAME .EQ. 'NUINI') THEN
          VALUE = NuIni
      Else
          IRET = 0
      EndIf
C
      RETURN
C                       ____________________________
C
      ENTRY EVLOUT (NOUUT)
C
C                             Write current values of parameters to unit NOUUT
C
      WRITE (NOUUT, 131) QINI, IPD0, NuIni, QMAX, IKNL,
     >      XMIN, XCR, NX, NT, JT, NfMx
C
  131 FORMAT ( /
     >' Current parameters and values are: '//
     >' Initiation parameters: Qini, Ipd0, NuIni = ', F8.2, 2I8 //
     >' Maximum Q, Order of Alpha:   Qmax, IKNL = ', 1PE10.2, I6//
     >' X- mesh parameters   : Xmin, Xcr,  Nx   = ', 2(1PE10.2), I8//
     >' LnQ-mesh parameters  : Nt,   Jt         = ', 2I8       //
     >' # of parton flavors  : NfMx             = ',  I8       /)
C
      RETURN
C                        ****************************
      END

      SUBROUTINE EVLGT1 (NAME, VALUE, IRET)
C                                                   -=-=- evlgt1
C
C                                       COPY OF EVLGET for the alternate set
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL LSTX
      CHARACTER*(*) NAME
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / X1ARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / Q1RAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / E1LPAR / AL, IKNL, IPD0, IHDN, NfMx
C
      IRET = 1
C
      IF     (NAME .EQ. 'QINI')  THEN
          VALUE = QINI
      ElseIF (NAME .EQ. 'IPD0')  THEN
          VALUE = IPD0
      ElseIF (NAME .EQ. 'IHDN') THEN
          VALUE = IHDN
      ElseIF (NAME .EQ. 'QMAX')  THEN
          VALUE = QMAX
      ElseIF (NAME .EQ. 'IKNL') THEN
          VALUE = IKNL
      ElseIF (NAME .EQ. 'XCR') THEN
          VALUE = XCR
      ElseIF (NAME .EQ. 'XMIN') THEN
          VALUE = XMIN
      ElseIF (NAME .EQ. 'NX') THEN
          VALUE = NX
      ElseIF (NAME .EQ. 'NT') THEN
          VALUE = NT
      ElseIF (NAME .EQ. 'JT') THEN
          VALUE = JT
      ElseIF (NAME .EQ. 'NFMX') THEN
          VALUE = NfMx
      Else
          IRET = 0
      EndIf
C
      RETURN
C                       ____________________________
C
      ENTRY EVLOT1 (NOUUT)
C
C                             Write current values of parameters to unit NOUUT
C
      WRITE (NOUUT, 131) QINI, IPD0, IHDN, QMAX, IKNL,
     >      XMIN, XCR, NX, NT, JT, NfMx
C
  131 FORMAT ( /
     >' Current parameters and values are: '//
     >' Initiation parameters: Qini, Ipd0, Ihdn = ', F8.2, 2I8 //
     >' Maximum Q, Kernel ID : Qmax, IKNL       = ', 1PE10.2, I6//
     >' X- mesh parameters   : Xmin, Xcr,  Nx   = ', 2(1PE10.2), I8//
     >' LnQ-mesh parameters  : Nt,   Jt         = ', 2I8       //
     >' # of parton flavors  : NfMx             = ',  I8       /)
C
      RETURN
C                        ****************************
      END
C

      SUBROUTINE EVLIN
C
C                                       Solicits parameters in EVL calculations
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      COMMON / IOUNIT / NIN, NOUT, NWRT

      CALL EVLGET ('QINI', V1, IR1)
      CALL EVLGET ('IPD0', V2, IR2)
      CALL EVLGET ('NUINI', V3, IR3)

    1 WRITE  (NOUT, 101) V1, V2, V3
  101 FORMAT (/' Current values of parameters QINI, IPD0, NuIni: ',
     >        1PE10.2, 2I6 / '$Type new values: ' )
      READ(NIN, *, ERR=302) V1, V2, V3
C
      CALL EVLSET ('QINI', V1, IRET1)
      CALL EVLSET ('IPD0', V2, IRET2)
      CALL EVLSET ('NUINI', V3, IRET3)
      goto 301
 302  write(nout,120)
      goto 1
C
 301  IF ((IRET1.NE.1) .OR. (IRET2.NE.1) .OR. (IRET3.NE.1)) THEN
   20   WRITE (NOUT, 120)
        GOTO 1
      EndIf

      CALL EVLGET ('QMAX', V1, IR1)
      CALL EVLGET ('NT',   V2, IR2)
      CALL EVLGET ('JT',   V3, IR3)
C
    2 WRITE  (NOUT, 102)  V1, V2, V3
  102 FORMAT (/' Current values of parameters QMAX, NT, JT: ',
     >        1PE10.2, I6, I6 / '$Type new values: ' )
      READ(NIN, *, ERR=304) V1, V2, V3
C
      CALL EVLSET ('QMAX', V1, IRET1)
      CALL EVLSET ('NT',   V2, IRET2)
      CALL EVLSET ('JT',   V3, IRET3)
      goto 303
 304  write(nout,120)
      goto 2
C
 303  IF ((IRET1.NE.1) .OR. (IRET2.NE.1) .OR. (IRET3.NE.1)) THEN
   22   WRITE (NOUT, 120)
        GOTO 2
      EndIf
C
      CALL EVLGET ('XMIN', V1, IR1)
      CALL EVLGET ('XCR',  V2, IR2)
      CALL EVLGET ('NX',   V3, IR3)
      CALL EVLGET ('NFMX',   V4, IR4)
      CALL EVLGET ('IKNL', V5, IR5)

    3 WRITE  (NOUT, 103)  V1, V2, V3, V4, V5
  103 FORMAT(/' Current values of parameters XMIN, XCR, NX,NFMX,IKNL: ',
     >  2(1PE12.3), 3I6 / '$Type new values: ' )
      READ(NIN, *, ERR=306) V1, V2, V3, V4, V5
C
      CALL EVLSET ('XMIN', V1, IRET1)
      CALL EVLSET ('XCR',  V2, IRET2)
      CALL EVLSET ('NX',   V3, IRET3)
      CALL EVLSET ('NFMX',   V4, IRET4)
      CALL EVLSET ('IKNL', V5, IRET5)
      Goto 305

  306 Write (Nout, 120)
      Goto 3
C
  305 IF ( (IRET1 .NE. 1) .OR. (IRET2 .NE. 1) .OR. (IRET3 .NE. 1)
     > .OR. (IRET4 .NE. 1) .OR. (IRET5 .NE. 1) ) THEN
   23   WRITE (NOUT, 120)
        GOTO 3
      EndIf
C
  120 FORMAT(' Bad values, Try again!' /)
C
      RETURN
C                        ****************************
      END

      Subroutine Evlini (Nini, Xmin, Qmax, Nx, Nq)
C !! Warning:
C !! The input xxx.ini file must be open and associated with Unit# = Nini

C Modified 7/24/99

C Made true Evlini -- i.e. evolve according to the parameters in xxx.ini
C Got rid off the fixed # of flavors ( = 5 ) of the previous version.
C If a different number of flavors than the original .ini file is needed,
C generate a separate .ini file, and edit that number!

C Made more "logical".

C   Explanations:
 
C * Read in xxx.ini file (Unit # = Nini) for input evolution parameters
C                                                            from CTEQ
C * Do evolution in the range (Xmin, 1) and (Qini, Qmax)  
C                                           (Qini is given in .ini)
C * Nx and Nq are # of grids in (x,Q) - tradeoff between accuracy/efficiency

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (NF0 = 5, Nshp = 8,NEX = Nshp+2)
 
      External FitIni
C
      Common 
     >  /InpFin/ A(-NF0:NF0), B(-1:Nshp, -NF0:NF0), Ifun, Igfn(-NF0:NF0)
     >  /ShpNfl/ Mshp, Nfval, Nfsea 

      DATA Ahdn, Qm6 / 1.D0, 180.D0 /
C                                          ----------------------------
C                                             Initialize Qcd parameters
      Call SetQcd

      Print *, 'Reading input XXX.ini.....'
C                                              Nini is the input Xnnn.ini file
      Call iniread(Nini, qini, anfl, aknl, ordi, blam, Qm4, Qm5)

C                                          ----------------------------
C                                           Set up evolution parameters
      CALL ParPdf(1, 'IHDN', Ahdn, JR)
      Call ParQcd(1, 'M6', Qm6, jr)
C                                          From subroutine arguments
      Call ParPdf(1, 'QMAX', qmax, jr)
      Call ParPdf(1, 'XMIN', xmin, jr)
      Anx = Nx
      Ant = Nq
      Call ParPdf(1, 'NX', anx, jr)
      Call ParPdf(1, 'NT', ant, jr)
C                                           From xxx.ini
      CALL ParQcd(1, 'ORDR', ordi, JR)
      Call ParQcd(1, 'M4', Qm4, jr)
      Call ParQcd(1, 'M5', Qm5, jr)
      CALL ParQcd(1, 'NFL', Anfl, JR)

      CALL ParPdf(1, 'QINI', QINI, IR)
      CALL Evlpar(1, 'NFMX', Anfl, JR)
      CALL ParPdf(1, 'IKNL', aknl, JR)
C                                     Set Lambda for Anfl effective flavors
      Neff = Nint (Anfl)
      Iorda = Nint (Ordi)
      Call SetLam (Neff, Blam, Iorda)

       Print '(A/A)', 'Initial distribution parameters:'
     > ,'Igf    Mom     A0      A1      A2      A3      A4  ...'
       do 11 ifl = Nfval, Nfsea,-1
 11       print '(I3,10f7.3)', Igfn(ifl), (b(iex,ifl),iex=-1,Mshp)
C
      Print *, 'Start evolution ....'
C                                      ----------------------------
      Anout = 6D0
      Print *, 'Input evolution parameters are'
      Call ParPdf (4, 'DUM', ANout, Jr)
      Call Evolve (FitIni, Irt)
      Print *, 'Evolution Completed!'

      Return
C                        ****************************
      END
      SUBROUTINE EVLPAR (IACT, NAME, VALUE, IRET)
C
C               For STANDARD codes on Iact and Iret, see SUBROUTINE PARPDF
C
C               Additional options added for Version 7.2 and up of PDF package:
C
C          IACT =        5      find value of variable from the alternate set.
C                        6      type list of all values from the alternate set.
C
c               These options must be called from EVLPAR directly, not from the
C               front-end unit PARPDF because PARQCD cannot handle Iact > 4
C
C               NAME is assumed upper-case.
C               If necessary, VALUE is converted to integer by NINT
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      CHARACTER*(*) NAME
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      IRET = 1
      IF     (IACT .EQ. 0) THEN
              WRITE ( NINT (VALUE) , 101)
  101         FORMAT (/ ' Initiation parameters:   Qini, Ipd0, Ihdn ' /
     >                  ' Maximum Q, Order of Alpha:     Qmax, IKNL ' /
     >                  ' X- mesh parameters   :   Xmin, Xcr,   Nx  ' /
     >                  ' LnQ-mesh parameters  :         Nt,   Jt   ' /
     >                  ' # of parton flavors  :         NfMx       ' /)
              IRET = 4
      ElseIF (IACT .EQ. 1) THEN
              CALL EVLSET (NAME, VALUE, IRET)
      ElseIF (IACT .EQ. 2) THEN
              CALL EVLGET (NAME, VALUE, IRET)
      ElseIF (IACT .EQ. 5) THEN
              CALL EVLGT1 (NAME, VALUE, IRET)
      ElseIF (IACT .EQ. 3) THEN
              CALL EVLIN
              IRET = 4
      ElseIF (IACT .EQ. 4) THEN
              CALL EVLOUT ( NINT (VALUE) )
              IRET = 4
      ElseIF (IACT .EQ. 6) THEN
              CALL EVLOT1 ( NINT (VALUE) )
              IRET = 4
      Else
              IRET = 3
      EndIf
C
      RETURN
C                        ****************************
      END

      SUBROUTINE EVLRD1 (NU, HEADER, IRR)
C                                                   -=-=- evlrd1
C
C $Header: /home/wkt/common/flib/2evl/RCS/evlrd1.f,v 8.6 2003/06/30 00:41:49 wkt Exp $
C $Log: evlrd1.f,v $
C Revision 8.6  2003/06/30 00:41:49  wkt
C 1. Increase # of valence quarks to MxVal = 3 to allow s.ne.sbar;
C 2. Improve on the enforcement of positivity of the pdf's; warning only of pdf.lt.-1d-4
C
C Revision 8.5  2003/06/29 22:44:13  wkt

C *** Start another round of archiving ***
C
C                               READ DATA FROM FILE TO ALTERNATE COMMON BLOCKS
C
C       See comment lines of the Entry EVLRD section for detailed instructions.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      CHARACTER HEADER*78
C
      PARAMETER (MXX = 204, MXQ = 25, MxF = 6, MxVal = 3)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / X1ARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / Q1RAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / Q1RAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / E1LPAR / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / P1VLDT / UPD(MXPQX), KF, Nelmt
C
      Read  (Nu, '(A)') Header

      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, Am1, Am2, Am3, AM4, AM5, AM6

      Read  (Nu, '(A)') Line
      Read  (Nu, *) IPD0, IHDN, IKNL, NfMx, Nfval, KF, Nelmt

      Read  (Nu, '(A)') Line
      Read  (Nu, *) NX,  NT, JT,  NG, NTL(NG+1)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) (NTL(I), NTN(I), TLN(I), DTN(I), I =1, NG)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QV(I), TV(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, XCR,  (XV(I), I =1, NX)
C
C                  Since quark = anti-quark for nfl> Nfval at this stage,
C                  we Read  out only the non-redundent data points
C                  No of flavors = NfMx sea + 1 gluon + Nfval valence
      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+Nfval)
      Read  (Nu, '(A)') Line
      CALL RdUpd (UPD, Npts, NU, IRR)

      CLOSE (NU)
      IF (IRR .NE. 0) PRINT *, 'Read error in EVLRD1'

      If (NfMx .GT. Nfval) Then
C                                       Refill Upd for Nfval+1 -> NfMx
      Do 11 Nflv = Nfval+1, NfMx
         J0 = (-Nflv + NfMx) * Nblk
C offset = NfMx sea + 1 gluon + (Nflv-1) less-massive quarks = NfMx+Nflv
         J1 = ( Nflv + NfMx) * Nblk
         Do 12 I = 1, Nblk
            Upd (J1 + I) = Upd (J0 + I)
 12      Continue
 11   Continue
      EndIf
C
      PRINT '(/A/1X,A/)',' EVL Parameters for this set are:', HEADER
      CALL EVLPAR (6, 'ALL', DBLE(NOUT), IR)
      PRINT '(/ A /)',   ' QCD Parameters for this set are:'
      PRINT '(A, F7.3, 2F7.0)', ' ALAM, NFL, NRDR = ', AL, FL, DR
      PRINT '(A, 6F7.3)',' M1, M2, ..., M6 = ', AM1,AM2,AM3,AM4,AM5,AM6
C                                           Compare two sets of QCD parameters
      CALL PARPDF (2, 'ALAM', BL, IR)
      CALL PARPDF (2, 'ORDER',RDR, IR)
      CALL PARPDF (2, 'M1', BM1, IR)
      CALL PARPDF (2, 'M2', BM2, IR)
      CALL PARPDF (2, 'M3', BM3, IR)
      CALL PARPDF (2, 'M4', BM4, IR)
      CALL PARPDF (2, 'M5', BM5, IR)
      CALL PARPDF (2, 'M6', BM6, IR)
      CALL PARPDF (2, 'NFL', BFL, IR)
      LF = NINT(2.*BFL+2.)
C
      SML = 1.E-4
      DIF = MAX(ABS(AM1-BM1), ABS(AM2-BM2), ABS(AM3-BM3), ABS(AM4-BM4),
     >          ABS(AM5-BM5), ABS(AM6-BM6), ABS(RDR- DR), ABS(AL-BL))
      IF (DIF .GT. SML .OR. KF .NE. LF) THEN
       PRINT *, 'Warning!  Two PDF sets have different QCD parameters.'
       PRINT '(/A/)', ' Parameters for the regular set are:'
       CALL PARPDF(4, 'ALL', DBLE(NOUT), IR)
      EndIf
C
      RETURN
C                        ****************************
      END
      SUBROUTINE EVLSET (NAME, VALUE, IRET)
C                                                   -=-=- evlset
C
C                               Sets variable named NAME to VALUE.
C
C               IRET =   0      variable not found.
C                        1      success.
C                        2      variable found, but bad value.
C
C               NAME is assumed upper-case, and VALUE real.
C               If necessary, VALUE is converted to integer by NINT
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      CHARACTER*(*) NAME
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
     > / PdfSwh / Iset, IpdMod, Iptn0, NuIni
C
      IRET = 1
      IF     (NAME .EQ. 'QINI')  THEN
          IF (VALUE .LE. 0) GOTO 12
          QINI = VALUE
      ElseIF (NAME .EQ. 'IPD0')  THEN
          ITEM = NINT(VALUE)
          IF (Item .Eq. 10 .or. Item .Eq. 11) GOTO 12
          IPD0 = ITEM
C
      ElseIF (NAME .EQ. 'IHDN') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. -1 .OR. ITEM .GT. 5) GOTO 12
          IHDN = ITEM
      ElseIF (NAME .EQ. 'QMAX')  THEN
          IF (VALUE .LE. QINI) GOTO 12
          QMAX = VALUE
      ElseIF (NAME .EQ. 'IKNL') THEN
          ITMP = NINT(VALUE)
          ITEM = ABS(ITMP)
cpn 2001 iknl == 12 is now legal
cpn          IF (ITEM.NE.1.AND.ITEM.NE.2) GOTO 12
          IF (ITEM.NE.1.AND.ITEM.NE.2.AND.ITEM.NE.12) GOTO 12
          IKNL = ITMP
      ElseIF (NAME .EQ. 'XCR') THEN
          IF (VALUE .LT. XMIN .OR. VALUE .GT. 10.) GOTO 12
          XCR = VALUE
          LSTX = .FALSE.
      ElseIF (NAME .EQ. 'XMIN') THEN
          IF (VALUE .LT. 1D-7 .OR. VALUE .GT. 1D0) GOTO 12
          XMIN = VALUE
          LSTX = .FALSE.
      ElseIF (NAME .EQ. 'NX') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 10 .OR. ITEM .GT. MXX-1) GOTO 12
          NX = ITEM
          LSTX = .FALSE.
      ElseIF (NAME .EQ. 'NT') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 2 .OR. ITEM .GT. MXQ) GOTO 12
          NT = ITEM
      ElseIF (NAME .EQ. 'JT') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 1 .OR. ITEM .GT. 5) GOTO 12
          JT = ITEM
      ElseIF (NAME .EQ. 'NFMX') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LT. 1 .OR. ITEM .GT. MXPN) GOTO 12
          NfMx = ITEM
      ElseIF (NAME .EQ. 'IPDMOD') THEN
          ITEM = NINT(VALUE)
          IF (Abs(Item) .Gt. 1) GOTO 12
          IpdMod = ITEM
      ElseIF (NAME .EQ. 'IPTN0') THEN
          ITEM = NINT(VALUE)
          IF (ABS(ITEM) .GT. MXF) GOTO 12
          IPTN0 = ITEM
      ElseIF (NAME .EQ. 'NUINI') THEN
          ITEM = NINT(VALUE)
          IF (ITEM .LE. 0) GOTO 12
          NuIni = ITEM
      Else
          IRET = 0
      EndIf
C
      RETURN
C                                                                  Error exit:
   12 IRET = 2
C
      RETURN
C                        ****************************
      END
C                                                          =-=-= Pdfmis
      Subroutine EVLWT (NU, HEADER, IRR)

C $Header: /home/wkt/common/flib/2evl/RCS/evlwt.f,v 8.6 2003/06/30 00:41:49 wkt Exp $
C $Log: evlwt.f,v $
C Revision 8.6  2003/06/30 00:41:49  wkt
C 1. Increase # of valence quarks to MxVal = 3 to allow s.ne.sbar;
C 2. Improve on the enforcement of positivity of the pdf's; warning only of pdf.lt.-1d-4
C
C Revision 8.5  2003/06/29 22:44:13  wkt

C *** Start another round of archiving ***

C ====================================================================
C GroupName: Evlaux
C Description: Auxillary modules for evolution
C ListOfFiles: evlwt wtupd evlrd1 hqrk fpin fpinms spnint spence xli
C ====================================================================

C $Header: /home/wkt/common/flib/2evl/RCS/evlwt.f,v 8.6 2003/06/30 00:41:49 wkt Exp $
C $Log: evlwt.f,v $
C Revision 8.6  2003/06/30 00:41:49  wkt
C 1. Increase # of valence quarks to MxVal = 3 to allow s.ne.sbar;
C 2. Improve on the enforcement of positivity of the pdf's; warning only of pdf.lt.-1d-4
C
c Revision 7.1  99/08/20  23:14:37  wkt
c X1ARAY common block synchronized.
c 
c Revision 7.0  99/07/25  18:39:00  wkt
c Pavel's improved version: xli.f added.
c 
c Revision 6.3  98/08/11  12:36:09  wkt
c Gauss ==> GausInt
c
c Revision 6.2  97/11/15  19:03:38  wkt
c bug fix in EvlRd1: Read line missing before RdUpd
c
c Revision 6.1  97/11/15  17:17:45  wkt
c v.6.0 + HLL revisions
c
c Revision 5.96  96/10/28  00:20:44  wkt
c logistics
c
c Revision 5.93  96/10/27  15:20:40  wkt
c Cosmatics
c
c Revision 5.9  96/10/20  13:26:33  wkt
c 1. Consolidated with Liang's Cteq4 version;
c 2. Cleaned up numerous nuisances (variables not used, ..) which cause
c    compilation warnings
c 3. ZfrmX and XfrmZ and ancillaries replaced with new version based on
c    fractional power rather than log + linear transformation;
c 4. New Cteq, Mrs, Grv switches put in; but need their programs and
c    tables to run.
c
c Revision 5.8  96/10/15  23:43:18  wkt
c minor bug fixes
c
c Revision 5.5  96/06/02  23:03:16  wkt
c Write and Read to file converted to 'Formatted'.
c EvlWt and WtUpd and EvlRd1 modules moved to EvlAux package.
c Commons EVLPAC and PEVLDT modified.
c
c Revision 5.1  95/07/19  14:08:19  wkt
c comments only
c
C ====================================================================
C      Subroutine EVLWT (NU, HEADER, IRR)
C
C     Before calling this routine, the calling program must open an external
C              file for data transfer, provide a filename for the file, and
C              establish the equivalence of that file with unit=Nu, such as:
C
C       Open (Nu, File='FILENAME', Form='FORMATTED', status='NEW', ... )
C
C
C       The file is not automatically close because it is also used for memory
C       transfer (PEVLDT ---> PEVLD1) which required reading from this unit.

C       Input parameter:   Nu = unit number for the external file to be
C                               written to. (established in calling program)
C
C                      Header = informational header for the file
C
C       Output parameter: Irr
C
C                       ----------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      CHARACTER HEADER*78, Line*80
C
      PARAMETER (MXX = 204, MXQ = 25, MxF = 6, MxVal = 3)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX), KF, Nelmt
     >  /ShpNfl/ Mshp, Nfval, Nfsea 

      CALL PARPDF (2, 'ALAM', AL, IR)
      CALL PARPDF (2, 'NFL',  FL, IR)
      CALL PARPDF (2, 'ORDER', DR, IR)
      CALL PARPDF (2, 'M1', AM1, IR)
      CALL PARPDF (2, 'M2', AM2, IR)
      CALL PARPDF (2, 'M3', AM3, IR)
      CALL PARPDF (2, 'M4', AM4, IR)
      CALL PARPDF (2, 'M5', AM5, IR)
      CALL PARPDF (2, 'M6  ', AM6, IR)
C
      Write (Nu, '(A)') Header

      Write (Nu, '(2x, A, 8x, A, 9x, A)')
     $ 'Ordr, Nfl, lambda', 'Qmass 1,  2,  3,',  '4,  5,  6'
      Write (Nu, '(2F5.0, F8.4, 6F9.3)')
     $  Dr,   Fl,  Al,       Am1, Am2, Am3, AM4, AM5, AM6

      Write (Nu, '(3x, A)') 'IPD0, IHDN, IKNL, NfMx, Nfval,  KF, Nelmt'
      Write (Nu, '(6I6, I8)') IPD0, IHDN, IKNL, NfMx, Nfval, KF, Nelmt

      Write (Nu, '(3x,A)') 'NX,  NT,  JT,  NG, NTL(NG+1)'
      Write (Nu, '(5I5)') NX,  NT,  JT,  NG, NTL(NG+1)

      Write (Nu, '(A)') '(NTL(I), NTN(I), TLN(I), DTN(I), I =1, NG)'
      Write (Nu, '(2I5, 1pE13.5, E13.5)')
     $                   (NTL(I), NTN(I), TLN(I), DTN(I), I =1, NG)

      Write (Nu, '(A)') 'QINI, QMAX, (QV(I), TV(I), I =0, NT)'
      Write (Nu, '(1pE13.5, E13.5 / (2E13.5))')
     $                   QINI, QMAX, (QV(I), TV(I), I =0, NT)

      Write (Nu, '(A)') 'XMIN, XCR, (XV(I), I =1, NX)'
      Write (Nu, '(1pE13.5, E13.5 / (6E13.5))')
     $                   XMIN, XCR, (XV(I), I =1, NX)
C
C                  Since quark = anti-quark for nfl> Nfval at this stage,
C                  we write out only the non-redundent data points
C                  No of flavors = NfMx sea + 1 gluon + Nfval valence
      Write (Nu, '(A)') 'Parton Distribution Table:'
      Npts = (NX+1) * (NT+1) * (NfMx+1+Nfval)
      CALL WTUPD (UPD, Npts, NU, IR2)
C
      IRR = IR2
C
      IF (IRR .NE. 0) PRINT *, 'Write error in EVLWT'
C
      RETURN
C                                               ----------------------------
C                                                    Read external file
      ENTRY EVLRD (NU, HEADER, IRR)
C
C       See comment lines of the Entry EVLWT section for detailed instructions.
C
C       The call program must open the existing file to be read from with:
C
C       Open (Nu, File='FILENAME', Form='UNFORMATTED', status='OLD', Err=... )
C
C       where the file named 'FILENAME' must exist and it has to be written
C       originally by EVLWT.  (cf. above)
c
C                                     HEADER is an output character variable.
C
C                       ----------------------------
C
      Read  (Nu, '(A)') Header

      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, Am1, Am2, Am3, AM4, AM5, AM6

      Read  (Nu, '(A)') Line
      Read  (Nu, *) IPD0, IHDN, IKNL, NfMx, Nfval, KF, Nelmt

      Read  (Nu, '(A)') Line
      Read  (Nu, *) NX,  NT, JT,  NG, NTL(NG+1)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) (NTL(I), NTN(I), TLN(I), DTN(I), I =1, NG)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) QINI, QMAX, (QV(I), TV(I), I =0, NT)

      Read  (Nu, '(A)') Line
      Read  (Nu, *) XMIN, XCR,  (XV(I), I =1, NX)
C
C                  Since quark = anti-quark for nfl>2 at this stage,
C                  we Read  out only the non-redundent data points
C                  No of flavors = NfMx sea + 1 gluon + Nfval valence
      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+Nfval)
      Read  (Nu, '(A)') Line
      CALL RdUpd (UPD, Npts, NU, IRR)

      CLOSE (NU)
      IF (IRR .NE. 0) PRINT *, 'Read error in EVLRD'

      If (NfMx .GT. Nfval) Then
C                                       Refill Upd for Nfval+1 -> NfMx
      Do 11 Nflv = Nfval+1, NfMx
         J0 = (-Nflv + NfMx) * Nblk
C offset = NfMx sea + 1 gluon + (Nflv-1) less-massive quarks = NfMx+Nflv
         J1 = ( Nflv + NfMx) * Nblk
         Do 12 I = 1, Nblk
            Upd (J1 + I) = Upd (J0 + I)
 12      Continue
 11   Continue
      EndIf
C                             To check read-in result against saved file.
C      Print '(A/(1pE13.5, 5E13.5))', 'Upd', (UPD(I), I=1,Nelmt-1)
C
C                                       Set QCD parameters to current values
C      CALL PARPDF (1, 'ALAM', AL, IR)
C      CALL PARPDF (1, 'NFL',  FL, IR)
C      CALL PARPDF (1, 'ORDER', DR, IR)
C      CALL PARPDF (1, 'M1', AM1, IR)
C      CALL PARPDF (1, 'M2', AM2, IR)
C      CALL PARPDF (1, 'M3', AM3, IR)
C      CALL PARPDF (1, 'M4', AM4, IR)
C      CALL PARPDF (1, 'M5', AM5, IR)
C      CALL PARPDF (1, 'M6  ', AM6, IR)
C
C      PRINT '(/A/1X,A/)', ' Parameters from this data file are:',HEADER
C      CALL PARPDF(4, 'ALL', DBLE(NOUT), IR)
C                                               ------------------------------
      Return
C                          *****************
      End
      SUBROUTINE EVOLVE (FINI, IRET)
C                                                   -=-=- Evolve
C #Header: /Net/u52/wkt/h22/2evl/RCS/Evolve.f,v 7.0 99/07/25 18:42:40 wkt Exp $
C #Log: Evolve.f,v $

C ====================================================================
C GroupName: Evolve
C Description: Main evolution modules
C ListOfFiles: evolve nsevl nsrhsm nsrhsp snevl snrhs kernel pff1 stupkl
C ====================================================================
c Revision 7.0  99/07/25  18:42:40  wkt
c Pavel's improved version: Diagonal kernel calculation extensively
c modified, with subtraction at x = 1; 
c 
c Revision 6.2  98/07/21  16:11:03  wkt
c GAUSS --> GausInt to conform with change in UtlPac6
c
c Revision 6.1  97/11/15  17:19:39  wkt
c v.6.0 + HLL revisions
c

C ====================================================================
C      SUBROUTINE EVOLVE (FINI, IRET)
C
C               Input argument: FINI is a function
C
C                               FINI (LPARTN, X)
C
C                     where LPARTN = -6, ... 6 labels the parton flavor:
C                     t-bar(-6), b-bar(-5), ... gluon(0), u(1), ... t(6) res.
C
C               Output parameter:
C
C                       Iret = 0  :  normal execution
C
C                       ----------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      PARAMETER (MXX = 204, MXQ = 25, MxF = 6, MxVal = 3)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX), Kf, Nelmt
C
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
C
      DIMENSION QRKP(MXF), XxNv(Mxx), QrNv(Mxx), AqNv(Mxx)
      DIMENSION JI(-MXF : MXF+1)
C
      EXTERNAL NSRHSP, NSRHSM, FINI

      DATA D0, PdfMin, Mprt / 0.0, -1D-4, 0 /
      Data Xvd /0.3D0/  !where PDF positivity condition starts 
C
   11 IRET = 0
C                                             Set up number of "valence quarks"
C               3 allows s .ne. sbar; but must ensure that Int(s-sbar) = 0
C      MXVAL = 3						! now put in as a parameter statement
C                                    Set up X-mesh points and common parameters
      IF (.NOT. LSTX) CALL XARRAY
      DLX = 1D0 / NX
      J10 = NX / 10
CPN                                                       Set up QCD factors
       call StpQCD

C                                                       Set up Q-related arrays
      CALL PARPDF (2, 'ALAM', AL, IR)
C                         Nini and NfMx determined by Qini and Qmax in Qarray
      CALL QARRAY (NINI)

C                            Flavor # NFSN is for the singlet quark combination
      NFSN = NFMX + 1
C                  Kf is the range of the flavor index for the Upd data array
C                 = quark flavors + 1 for singlet combination and 1 for gluon
      KF = 2 * NFMX + 2
C                                 Total number data points to be stored in Upd
      Nelmt = Kf * (Nt+1) * (Nx+1)
C                                                       ----------------------
C                                                        Various initiations:
C
C      Offset, or Starting address -1, of Pdf(Iflv) in Upd for the first stage.
C                                                    Will be updated each turn
      DO 101 IFLV = -NFMX, NFMX+1
        JFL = NFMX + IFLV
        JI(IFLV) = JFL * (NT+1) * (NX+1)
  101 CONTINUE
C                               Define initial distributions for the evolution.
C                                        Input functions are defined on (1:Nx).
C
C                             Input gluon and quark distributions are inserted
C                          in the Upd array directly;   Output of each stage of
C                              evolution serve as the input for the next stage.
C
    3 DO 31 IZ = 1, NX
C                              Gluon:  (Position Ji(0)+1 is reserved for x = 0)
C
        UPD(JI(0)+IZ+1) = FINI (0, XV(IZ))
C                                                                 Singlet quark
C                        Input singlet quark distribution must be initiated for
C                                                        each stage separately.
C
        UPD(JI(NFSN)+IZ+1) = 0
C                                   If no quark, bypass filling of quark arrays
        IF (NFMX .EQ. 0) GOTO 31
C
        DO 331 IFLV = 1, NINI
          A = FINI ( IFLV, XV(IZ))
          B = FINI (-IFLV, XV(IZ))
          QRKP (IFLV) = A + B
C                                               Acculumate singlet distribution
          UPD(JI(NFSN)+IZ+1) = UPD(JI(NFSN)+IZ+1) + QRKP (IFLV)
C                                       Initialize the "minus" non-singlet com-
C                                        bination (Q - Q-bar) in area for Q-bar
          UPD(JI(-IFLV)+IZ+1) = A - B
  331   CONTINUE
C                            The "plus" non-singlet combination is initialized
C                           in the array area for the Quark of the same flavor
        DO 332 IFLV = 1, NINI
           UPD(JI( IFLV)+IZ+1) = QRKP(IFLV) - UPD(JI(NFSN)+IZ+1)/NINI
  332   CONTINUE
C
   31 CONTINUE
C                                                       -----------------------
C                                                   Start of the Q2- Evolution:
C
C                                 Outer loop is by the effective number of flvr
      DO 21 NEFF = NINI, NFMX
C                                    Set up 1st and 2nd order kernel functions
 
cpn 2001        IF (abs(IKNL) .EQ. 2) CALL STUPKL (NEFF)
        IF (abs(IKNL) .GE. 2) CALL STUPKL (NEFF)
C                                                           Singlet Calculation
        ICNT = NEFF - NINI + 1
C                                                 Skip if new quark mass = old
        IF (NTN(ICNT) .EQ. 0) GOTO 21
C                                       Otherwise, recall iteration parameters
        NITR = NTN (ICNT)
        DT   = DTN (ICNT)
        TIN  = TLN (ICNT)
C                                                            Perform evolution
        CALL SNEVL (IKNL, NX, NITR, JT, DT, TIN, NEFF,
     >   UPD(JI(NFSN)+2), UPD(JI(0)+2), UPD(JI(NFSN)+1), UPD(JI(0)+1
     >       ))
C
C                                        Non-singlet sector, one flavor a time.
C
C                                        Skip this section if quark flavor = 0.
        IF (NEFF .EQ. 0) GOTO 88
C
 5      DO 333 IFLV = 1, NEFF
C                             First evolve the (Q+Q-bar) "PLUS NS" part
          CALL NSEVL (NSRHSP, IKNL, NX, NITR, JT, DT, TIN, NEFF,
     >         UPD(JI( IFLV)+2), UPD(JI( IFLV)+1))
C
C                                  The (Q-Q-bar) "MINUS NS" evolution is needed
C                                           only for those flavors with valence
          IF (IFLV .LE. MXVAL)
     >         CALL NSEVL (NSRHSM, IKNL, NX, NITR, JT, DT, TIN, NEFF,
     >         UPD(JI(-IFLV)+2), UPD(JI(-IFLV)+1))
C
C                            To obtain  the real  quark distribution functions,
C                      combine the singlet piece to the two non-singlet pieces.
C                            Enforce positivity conditions also at this stage.
cpn 2001     The positivity condition is partially disabled for X <= Xvd,
c            where X>Xvd is considered to be the region dominated by   
c            the valence quarks. (let Xvd=0.3)                            
C
          DO 55 IS = 0, NITR
	      Nnv = 1 
            DO 56 IX = 0, NX
              TP = UPD (IS*(NX+1) + IX + 1 + JI( IFLV))
              TS = UPD (IS*(NX+1) + IX + 1 + JI( NFSN)) / NEFF
              TP = TP + TS
C              IF (IKNL.GT.0 .and. XV(IX).gt.Xvd) TP=MAX(TP,D0)
c                                                                IF (IKNL.GT.0) TP=MAX(TP,D0)
              
              IF (IFLV .LE. MXVAL) THEN
                TM = UPD (IS*(NX+1) + IX + 1 + JI(-IFLV))
C                IF (IKNL.GT.0 .and. XV(IX).gt.Xvd) THEN
c                                                                IF (IKNL.GT.0) THEN
C                  TM = MAX (TM, D0)
C                  TP = MAX (TP, TM)
C                EndIf
              Else
                TM = 0.
              EndIf

              Quark = (TP + TM)/2.
C											      Enforce positivity and control warning switch
              If (Quark.LT.D0 .and. IKNL.GE.0) Then
	           If (Quark .LT. PdfMin) Mprt = 1
			   Quark1 = Quark
			   Quark  = 0
			EndIf
              UPD (JI( IFLV) + IS*(NX+1) + IX + 1) = Quark

	        AnQrk = (TP - TM)/2.
              If (AnQrk.LT.D0 .and. IKNL.GE.0) Then
	           If (AnQrk .LT. PdfMin) Mprt = 1
			   AnQrk1 = AnQrk
			   AnQrk  = 0
			EndIf
              UPD (JI(-IFLV) + IS*(NX+1) + IX + 1) = AnQrk
C													Print a warning about negative Pdfs
	If (Mprt .EQ. 1) Then
	   XxNv(Nnv) = XV(Ix)
	   QrNv(Nnv) = Quark1
	   AqNv(Nnv) = AnQrk1
	   Nnv = Nnv + 1
	   Mprt = 0
	EndIf	   

 56         CONTINUE
      If (Nnv .GT. 1) Then
          If (Iflv0 .NE. Iflv) Then
            Write (Nwrt, '(A, I5)') 'Flavor number in evolve =', Iflv
	      Iflv0 = Iflv
	    EndIf
	    QQ = QV(NTL(ICNT) + IS)
	    Nprt = Min(Nnv-1, 10)                    ! We shall print a max of 10 values
	    Write (Nwrt, '(3x, A, 1pE12.4, A, I6)') 
     >          'Q = ', QQ, ' # of Negative Pdfs = ', Nnv-1
	    Write (Nwrt, '(5x, A, 1pE11.3, 9E11.3)') 
     >          '  X  = ', (XxNv(I), I=1,Nprt)
	    Write (Nwrt, '(5x, A, 1pE11.3, 9E11.3)') 
     >          'Qrk  = ', (QrNv(I), I=1,Nprt)
	    Write (Nwrt, '(5x, A, 1pE11.3, 9E11.3)') 
     >          'AnQr = ', (AqNv(I), I=1,Nprt)
	EndIf
 55       CONTINUE
 333    CONTINUE
C                                                                Check results
C      DO 533 IFLV = 0, NEFF, 2
C         WRITE (NOUT,'(/A, 2I5)') ' Neff, Iparton =', NEFF, IFLV
C         WRITE (44,  '(/A, 2I5)') ' Neff, Iparton =', NEFF, IFLV
C      DO 533 IS = 1, NITR, 3
C         QQ = QV(NTL(ICNT) + IS) - 0.00001
C         WRITE (NOUT,'(/A, I5, 1PE15.3)') ' IQ, Q =', IS, QQ
C         WRITE (44,  '(/A, I5, 1PE15.3)') ' IQ, Q =', IS, QQ
C      DO 533 IX = 1, NX, J10
C          XX = XV(IX) + 0.00001
C          TM = UPD (JI( IFLV) + IS*(NX+1) + IX + 1)
C          TL = DPDF(10, IHDN, IFLV, XX, QQ, JR)
C          TN = DPDF(IPD0, IHDN, IFLV, XX, QQ, JR)
C          WRITE (NOUT,'(I5, 4(1PE15.3))') IX, XX, TM, TL, TN
C          WRITE (44,  '(I5, 4(1PE15.3))') IX, XX, TM, TL, TN
C  533 CONTINUE
C                   Heavy quarks above current threshold have zero distribution
C
 
        DO 334 IFLV = NEFF + 1, NFMX
          DO 57 IS = 0, NITR
            DO 58 IX = 0, NX
              UPD(JI( IFLV) + IS*(NX+1) + IX + 1) = 0
              UPD(JI(-IFLV) + IS*(NX+1) + IX + 1) = 0
 58         CONTINUE
 57       CONTINUE
 334    CONTINUE
 88     CONTINUE
C                                        ------------------------------
C                                      Define initial parameters for next stage
        IF (NFMX .EQ. NEFF) GOTO 21
C                                                                   New Offsets
        DO 335 IFLV = -NFMX, NFMX+1
          JI(IFLV) = JI(IFLV) + NITR * (NX+1)
  335   CONTINUE
C                                       New distributions:
C                                           gluon input functions are in place;
C                       only non-singlet and singlet quark needs re-initiation
C
C                             Calculate initial heavy quark distribution due to
C                           change of renormalization scheme across threshold:
        CALL HQRK (NX, TT, NEFF+1, UPD(JI(0)+2), UPD(JI(NEFF+1)+2))
C
        DO 32 IZ = 1, NX
          QRKP (NEFF+1) = 2. * UPD(JI( NEFF+1) + IZ + 1)
C                                                 New Singlet piece
          UPD (JI(NFSN)+IZ+1) = UPD (JI(NFSN)+IZ+1)  + QRKP (NEFF+1)
          VS00 =  UPD (JI(NFSN)+IZ+1) / (NEFF+1)
C
C                                        "plus" non-singlet for the new flavor
          UPD(JI( NEFF+1) + IZ + 1) = QRKP(NEFF+1) - VS00
C
C                 Calculate the non-singulet parts of the other quark distr.
C
C                               Change from the output of last stage of calcu-
C                               lation due to two sources of change in Vs/Neff:
C                               change of Neff and addition of new quark distr.
          DO 321 IFL = 1, NEFF
            A = UPD(JI( IFL)+IZ+1)
            B = UPD(JI(-IFL)+IZ+1)
            QRKP(IFL) = A + B
C                                             "plus" non-singlet for flavor IFL
            UPD(JI( IFL)+IZ+1) = QRKP(IFL) - VS00
C                          "minus" non-singlet for flavors with valence
            IF (IFL .LE. MXVAL)  UPD(JI(-IFL)+IZ+1) = A - B
 321      CONTINUE
C
 32     CONTINUE
C                               Return of Q-2 evolution loop to the next stage
 21   CONTINUE
C                                            Conclusion of the full calculation
      Return
C                   **********************
      End
C==============================================================================
C Fini is "user defined".  Here it is just a dummy place-holder.
C==============================================================================
C      FUNCTION FINI (IPARTN, XX)
C
C      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C      FINI = 0.
C
C      RETURN
C                          ****************
C      END
C==============================================================================
C
c -----------------------------------------------------------------------------
      FUNCTION FitIni (LP, XX)
C                                        Initial Parton distribution at Q = Q0
C ----------------------------------------------------------------
C
C     Re-written 04-07-17  (along with inptn.f inptsup.f and gfun)  --- WKT

C     This is still an interim version: the most natural and economic way of
C     writing this routine is to have it return the full Bfun() array, cf. below,
C     to a parent subroutine, which produces the full initial 2-dim PDF array
C     in (Lp, X) at Q=Q0 for Evolve.f .
C
C     The current version demonstrate that the above can easily be done, but
C     preserves the old interface to Evolve.f for the time being,

C ----------------------------------------------------------------
C     Note:  Since our flavor assignments are:
C
C            svl dvl uvl glu ubr dbr sbr chm btm;
C
C and since uvl=u-ubr, dvl=d-dbr, it is convenient to define SEA quarks as:
C
C                Ifl =  -1       -2           -3       -4
C                -------------------------------------------------
C              / 1     2*ubr   2*dbr          s+sbr    c+cbr ....
C Ifun schemes|  2    dbr-ubr  2*(dbr+ubr)    s+sbr    c+cbr ....
C              \ 3    dbr/ubr  2*(dbr+ubr)    s+sbr    c+cbr ....
C
C so that the sum of momenta of (val+glu+sea) is 1.0 (excluding Ifl=-1 for Ifun=2,3);
C the factors of 2 for Ifl=-1, -2 compensate for the ubr,dbr taken out in the val def.
C ------------------------------------------------------------------

      Implicit Double Precision (A-H, O-Z)

      Parameter (NF0 = 5, Nshp = 8)

      Common
     >  /InpFin/ A(-NF0:NF0), B(-1:Nshp, -NF0:NF0), Ifun, Igfn(-NF0:NF0)
     >  /ShpNfl/ Mshp, Nfval, Nfsea

      dimension Bfun(-NF0:NF0)

      DATA D1M, SML / 0.9999999999d0, 1.0D-10 /

      IF (LP .LT. -NF0 .OR. LP .GT. NF0) Then
        FitIni = 0.
        Return
      EndIf

      X = XX
      IF ((x .LE. SML) .OR. (x .GE. D1M)) Then
        FitIni = 0.
        Return
      EndIf

      Do Ibp = Nfsea, Nfval
        Bfun(Ibp) = A(Ibp) * gfun(Igfn(Ibp), x, B(1,Ibp))
      EndDo
C                                                     Depending on Ifun, redefine Ibp=-1/-2
C                                                       to ubar/dbar
      If (Ifun .Eq. 2) Then
C                                                     Bfun(-2) =  2 *(ub+db)
C                                                     Bfun(-1) =  tan[2atan(1)*(db-ub)/(db+ub)]
        dpu = Bfun(-2) /2D0
        dmu = dpu * atan(Bfun(-1))/atan(1d0)/2d0

        Bfun(-2) = dpu + dmu                     ! 2 * dbar
        Bfun(-1) = dpu - dmu                     ! 2 * ubar
      ElseIf (Ifun .Eq. 3) Then
        dou  = Bfun(-1)
        dpu2 = Bfun(-2)
        Bfun(-2) = dou * dpu2 /(dou +1D0)        ! 2 * dbar
        Bfun(-1) = Bfun(-2) / dou                ! 2 * ubar
      Else
        Print *, 'Ifun must be 1,2,3 in fitini; you have Ifun = ', Ifun
        Stop
      EndIf
C                                                     If s_valence != 0, calculate s-sbar
      If (Igfn(3) .Ne. 0) Then
        Bfun(3) = Bfun(-3) * atan(Bfun(3))/atan(1d0)/2d0
      EndIf

      If     (Lp .Eq. 3) Then
        Tem = (Bfun(-3) + Bfun(3)) /2D0             ! s
      ElseIf (Lp .Eq. 2) Then
        Tem = Bfun(2) + Bfun(-2) /2D0               ! d
      ElseIf (Lp .Eq. 1) Then
        Tem = Bfun(1) + Bfun(-1) /2D0               ! u
      ElseIf (Lp .Eq. 0) Then
        Tem = Bfun(0)                               ! g
      ElseIf (Lp .Eq.-1) Then
        Tem = Bfun(-1) /2D0                         ! ub
      ElseIf (Lp .Eq.-2) Then
        Tem = Bfun(-2) /2D0                         ! db
      ElseIf (Lp .Eq.-3) Then
        Tem =(Bfun(-3) - Bfun(3)) /2D0              ! sb
      Else
        Tem = Bfun(-Abs(Lp)) /2D0                   ! c, cb, b, bb, ...
      EndIf

      Fitini = Tem

      Return
C                        ****************************
      End
      function Fnsm_pn(x)
cpn                    Returns P_{NS,-}^{(2)} (conventional notations}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision Fnsm_pn, x
      double precision  PCF2,PCFG,PCFT,PQQB,PQQ2, FFP,FFM,x2,xln1m,
     >     xln1m2,xln,xln2, spen2,spenc2,tem,xi,xm1i,xp1i

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

cpn                    For the DIS scheme
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      integer iknl, ipd0, ihdn, nfmx
      double precision al, beta_c2q
      external beta_c2q

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PNSM: x=',x,' is out of range'

      XI = 1./ X
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      FFP = (1.+ X2) * XM1I
      FFM = (1.+ X2) * XP1I
 
      PCF2 = -2.* FFP *XLN*XLN1M - (3.*XM1I + 2.*X)*XLN
     >     - (1.+X)/2.*XLN2 - 5.*(1.-X)
      PCFG = FFP * (XLN2 + 11.*XLN/3.+ 67./9.- PI**2 / 3.)
     >     + 2.*(1.+X) * XLN + 40.* (1.-X) / 3.
      PCFT = (FFP * (- XLN - 5./3.) - 2.*(1.-X)) * 2./ 3.
      
      PQQB = 2.* FFM * SPEN2 + 2.*(1.+X)*XLN + 4.*(1.-X)
      PQQB = (CF2-CFCG/2.) * PQQB
      PQQ2 = CF2 * PCF2 + CFCG * PCFG / 2. + CFTR*NFlv * PCFT

      tem = PQQ2 - PQQB 

      if (iknl.eq.12) then             !for the DIS scheme: subtract the
                                       !transformation terms
        tem = tem - beta_c2q(x)
      endif 

      Fnsm_pn=tem

      Return
      End!Fnsm_pn


cpn---------------------------------------------------------------------------
      function Fnsp_pn(x)
cpn                    Returns P_{NS,+}^{(2)} (conventional notations}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision Fnsp_pn, x
      double precision PCF2,PCFG,PCFT,PQQB,PQQ2,FFP,FFM,x2,xln1m,xln1m2,
     >  xln,xln2,spen2,spenc2,tem,xi,xm1i,xp1i

      external Spenc2
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

cpn                    For the DIS scheme
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      integer iknl, ipd0, ihdn, nfmx
      double precision al, beta_c2q
      external beta_c2q

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PNSP: x=',x,' is out of range'

      XI = 1./ X
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      FFP = (1.+ X2) * XM1I
      FFM = (1.+ X2) * XP1I
 
      PCF2 = -2.* FFP *XLN*XLN1M - (3.*XM1I + 2.*X)*XLN
     >     - (1.+X)/2.*XLN2 - 5.*(1.-X)
      PCFG = FFP * (XLN2 + 11.*XLN/3.+ 67./9.- PI**2 / 3.)
     >     + 2.*(1.+X) * XLN + 40.* (1.-X) / 3.
      PCFT = (FFP * (- XLN - 5./3.) - 2.*(1.-X)) * 2./ 3.
      
      PQQB = 2.* FFM * SPEN2 + 2.*(1.+X)*XLN + 4.*(1.-X)
      PQQB = (CF2-CFCG/2.) * PQQB
      PQQ2 = CF2 * PCF2 + CFCG * PCFG / 2. + CFTR*NFlv * PCFT

      tem = PQQ2 + PQQB

      if (iknl.eq.12) then             !for the DIS scheme: subtract the
                                       !transformation terms
        tem = tem - beta_c2q(x)
      endif 
      

      Fnsp_pn=tem

      Return
      End!Fnsp_pn

cpn---------------------------------------------------------------------------
      FUNCTION FPIN (LPARTN, X)
C                                       Initial parton distribution at Q = Qini
C                                       from published PDF parametrizations
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      IR  = 0
      IR1 = 0
      IR2 = 0
      IR3 = 0
C
      CALL PARPDF (2, 'IPD0', SET, IR1)
      CALL PARPDF (2, 'IHDN', HDRN, IR2)
      CALL PARPDF (2, 'QINI', QINI, IR3)
C
      IF (IR1 .NE. 1) PRINT *, 'Ir1 = ', IR1, ' in FPIN'
      IF (IR2 .NE. 1) PRINT *, 'Ir2 = ', IR2, ' in FPIN'
      IF (IR3 .NE. 1) PRINT *, 'Ir3 = ', IR3, ' in FPIN'
C
      IH = HDRN
      IS = SET
C
      FPIN = PDF (IS, IH, LPARTN, X, QINI, IR)
C
      IF (IR .NE. 0) PRINT *, 'Ir = ', IR, ' in FPIN'
      RETURN
C                        ****************************
      END

      FUNCTION FPINMS (LPARTN, X)
C                                                   -=-=- fpinms
C                                       Initial parton distribution at Q = Qini
C                                       from published PDF parametrizations
C
C                                    Conversion from DIS to MS-bar scheme is
C                                  performed before the function is passed on.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      IR  = 0
      IR1 = 0
      IR2 = 0
      IR3 = 0
C
      CALL PARPDF (2, 'IPD0', SET, IR1)
      CALL PARPDF (2, 'IHDN', HDRN, IR2)
      CALL PARPDF (2, 'QINI', QINI, IR3)
C
      IF (IR1 .NE. 1) PRINT *, 'Ir1 = ', IR1, ' in FPINMS'
      IF (IR2 .NE. 1) PRINT *, 'Ir2 = ', IR2, ' in FPINMS'
      IF (IR3 .NE. 1) PRINT *, 'Ir3 = ', IR3, ' in FPINMS'
C
      IH = HDRN
      IS = SET
C              MRS sets are in MS-Bar, the others in DIS scheme to begin with.
      IF (IS .GE. 5 .AND. IS .LE. 7) THEN
        IACT = -1
      Else
        IACT = 1
      EndIf

      TEM = PMSDIS (IS, IH, LPARTN, X, QINI, IACT)
C
      FPINMS = TEM
C
      IF (IR .NE. 0) PRINT *, 'Ir = ', IR, ' in FPINMS'
      RETURN
C                        ****************************
      END

c --------------------------------------------------------------------
      function gfun(iform,x,aaa)

C wkt: 2006-03-08:
C  * Simplification of the general structure by using the Fortran 95 Case construct;
C  * Added Case (21) - two-component form designed to probe heavy flavors
C
C Revision 8.7  2003/07/23 02:34:30  wkt
C Implement the "General Parametrization" of s-sbar according to the description
C of the Dimuon analysis paper.
C
C Revision 8.6  2003/06/30 00:45:34  wkt
C Implement s.ne.sbar option, Issbr=17 --- interim version
C
C Revision 8.5  2003/06/29 23:07:26  wkt

C *** Start another round of archiving: v8.5. ***
C
C Revision 80.4  2003/05/04 17:46:49  wkt
C Replaced Jon's version, which does not satisfy positive definiteness, with a new
C version of mine (Ifun=17) which does.  The implementation is non-trivial. Cf.
C notes ParSsbr.tex (.pdf) for details.
C
C Revision 80.3  2003/05/03 00:29:51  wkt
C Jon Pumplin's version with isea=16 option
C
c functional forms used by fitini.f

      implicit double precision (a-h,o-z)

      dimension aaa(*)

      data one, two  / 1d0, 2d0 /

      xx = x

      a1 = aaa(1)
      a2 = aaa(2)
      a3 = aaa(3)
      a4 = aaa(4)
      a5 = aaa(5)
      a6 = aaa(6)
      a7 = aaa(7)
      a8 = aaa(8)

      std = xx **(a1-one) * (one-xx) **a2

      Select Case (iform)
c Null -- used for default differences (dmu, s-sbr, ...)
      Case (0)
         gfun = 0D0

c constant -- used for default ratios (dpu/s+sbr, ...)
      Case (1)
         gfun = a1

c traditional CTEQ
      Case (2)
         gfun = std *
     >      (one + (exp(a3) - one)*xx**a4 +
     >           sinh(a5)*sqrt(xx) + sinh(a6)*xx)

c CTEQ5 style for dbar/ubar...
      Case (3)
           gfun = exp(a1)*xx**(a2-one)*(one-xx)**a3 +
     >                  (one + a4*xx)*(one - xx)**a5

c CTEQ5 style for dbar/ubar modified to approach zero smoothly
c instead of going negative...
      Case (4)

           dou = exp(a1)*xx**(a2-one)*(one-xx)**a3 +
     >                  (one + a4*xx)*(one - xx)**a5

c (formula is chosen to preserve dou exactly at dou=1, and
c  for dou .ge. bconst)
           bconst = 10.d0
           if(dou .gt. bconst) then
              fac = dou
           elseif(dou .lt. -bconst) then
              fac = 0.d0
           else
              tmp = 1.d0 + exp(-bconst*dou) - exp(-bconst)
              fac = dou + (1.d0/bconst)*log(tmp)
           endif
           gfun = fac

      Case (5)
         tem = (a1-one)*log(xx) + a2*log(one-xx) + a3*xx +
     >               a5*log(one + xx*exp(a4))
         gfun = exp(tem)

      Case (6)               ! 06.12.21: New Standard for general positive-def distributions
         tem = a3*xx + a4*(one-xx)**2 + a5*xx**2
         gfun = std * exp(tem)

      Case (8)              ! 06.12.21: New Standard for "a-/a+ =tanh(gfun)" type of distributions
         gfun = xx**a1 *(one-xx/a2) *sinh(a3+ a4*xx+ a5*xx**2)

      Case (18)
C     "General parametrization" of s-sbar of 03-07-02 used for Dimuon/NuTeV anomaly study
         a3 = (1D0 + tanh(a3)) / 2D0
         gfun = a5 *
     >   xx**a1 *(one -xx/a3) *sinh(one + a2*xx +a4*xx**2)
C
      Case (28)             ! used as the argument for the tanh function (for dbar-ubar)
         gfun =             ! a4,a5, if positive, are potential crossing points
     >   xx**a1 *(one-xx)**a2 *a3 *(one -xx/a4) *exp(xx*a5)

      Case (21)
C                                two-component form, intended to probe heavy flavors
C
         tem1 = (a1-one)*log(xx) + a2*log(one-xx) + a3*xx +
     >               a5*log(one + xx*exp(a4))
         tem2 = (a7-one)*log(xx) + a8*log(one-xx)
         gfun = exp(tem1) + exp(a6 + tem2)

      Case (45)
	   tmp = a1 + (a2-one)*log(xx) + a3*log(one-xx)
           dou1 = exp(tmp)
	   dou2 = (one + a4*xx)*(one-xx)**a5

	   gfun = dou1 + dou2


      Case (55)             ! Superceded by (6), -a3, a4, a5 --> a4, a5, a3
         tem = (a1-one)*log(xx) + a2*log(one-xx)
     >       - a3*(one-xx)**2   + a4*xx**2  + a5*xx
         gfun = exp(tem)

      Case (60)
c                                     BHPS (Brodsky et al.) form for intrinsic charm (note, contains no free parameters)
           brak1 = (1.d0/3.d0)*(1.d0-x)*(1.d0+10.d0*x+x**2)
           brak2 = 2.d0*x*(1.d0+x)
           gfun = x**2*(brak1 + brak2*log(x))

      Case (82)                      ! For tanh();  obsolete: identical to (8).
         gfun = xx**a1 *(one-xx/a2)*sinh(a3+a4*xx+a5*xx*xx)

      Case (83)
         gfun = xx**a1 *(one-xx/a2)*sinh(a3+a4*xx*(1d0-xx))

      Case (84)
         gfun = xx**a1 *(one-xx/a2)*a3*exp(a4*xx+a5*xx*xx)


      Case (85)
         gfun = xx**a1 *(one-xx/a2)
     >        * (sinh(a3) + sinh(a4) *sqrt(xx) + sinh(a5) *xx)

      Case (91)               ! 07.02.18 : better parametrization to ensure positivity
         reg = one + a3 * xx + a4 * xx**2 + a5 * xx**3
         gfun = std * log(one + exp(reg))

      Case (92)               ! 07.02.19 : mixed exp and log exp
         reg = one + a4 * xx + a5 * xx**2
         gfun = std * exp(a3*xx) * log(one + exp(reg))

      Case (93)               ! 07.02.19 : restore old CTEQ parametrization in new setting
         reg = one + (exp(a3)-0.9)*xx**a4 + a5*sqrt(xx)
         gfun = std * log(one + exp(reg))

      Case (94)               ! 07.02.21 :  refined 91
         reg = (one +a3*1d1*xx**(1d0/3)) 
     > *    (one +a4*1d1*xx) *(one +a5*3d1*xx**3)
         gfun = std * log(one + exp(reg))

      Case (95)               ! 07.02.21 :  refined 91
         reg = sinh(a3)*xx**(1d0/3) +sinh(a4)*xx +sinh(a5)*xx**3
         gfun = std * log(one + exp(reg))

      Case (96)                      ! For atan .
         std =  xx**a1 *(one-xx/a2)
         gfun = std * (a3 + a4 *xx + a5 *xx**2)

      Case Default
c                                                        undefined below here...
         print *,'gfun fatal error: undefined iform=',iform
         stop

      End Select
C                         ************************
      end
      SUBROUTINE HINTEG (NX, F, H)
C                                                   -=-=- hinteg
C
C       Computes the integral [yF(y)-xF(x)]/(y-x) * dy/y over the range [x, 1];
C       then add F(x) * Ln (1-x) to get Int 1/(1-x/y)(sub+)F(y)dy/y.
C
C       The input function F must be specified on an array of size NX over
C       the range (0, 1] of the x variable.  In order to allow a possible
c       singularity at x = 0, the first mesh-point is at x = 1/nx, not 0.
C       The output function H is given on the same array as above.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / HINTEC / GH(NDG, MXX)
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
C
      DIMENSION F(NX), H(NX), G(MXX)
C
      DZ = 1D0 / (NX-1)
C                                       Loop to calculate the required integral
C                                               |--|--|-----------|--|--|
C                                         Iy:   I I+1 ...              nx
C                                         Kz:   1  2   ...              np
C
C                                                ------------------------------
C                                                     Do each integral in turn.
      DO 20 I = 1, NX-2
C                                         Number of points in the I-th integral
         NP = NX - I + 1
C                                  Evaluate the first two bins of the integrals

         TEM = GH(1,I)*F(I) + GH(2,I)*F(I+1) + GH(3,I)*F(I+2)

C                                  Evaluate the integrand for the I-th integral
         DO 30 KZ = 3, NP
            IY = I + KZ - 1
C                       DXDZ is the Jacobian due to the change of variable X->Z
C
            W = XA(I,1) / XA(IY,1)
            G(KZ) = DXTZ(IY)*(F(IY)-W*F(I))/(1.-W)
C
   30    CONTINUE
C
         HTEM = SMPSNA (NP-2, DZ, G(3), ERR)
         TEM1 = F(I) * ELY(I)
         H(I) = TEM + HTEM + TEM1
C
   20 CONTINUE
C
      H(NX-1) = F(NX) - F(NX-1) + F(NX-1) * (ELY(NX-1) - XA(NX-1,0))
      H(NX)   = 0
C
      RETURN
C                        ****************************
      END

      SUBROUTINE HQRK (NX, TT, NQRK, Y, F)
C                                                   -=-=- hqrk
C
C       Subroutine to compute the (heavy) quark distribution from (given)
C       gluon distribution, as the result of a change in the renormalization
C       scheme (from MS-bar to BPH) as the threshold for quark flavor Nqrk
C       is crossed.
C
C       Nx is the number of mesh-points, Tt is the Log Q variable.
C       Y is the input g-distribution function defined on the mesh points
C       F is the outpur Qrk-distribution function defined on the same pts.
C
C       If the threshold is chosen at the 'natural boundary' Mu = Mass(qrk),
C       then there is no renormalization and this routine returns zero.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
C
      DIMENSION Y(NX), F(NX)
C      DIMENSION W0(MXX), W1(MXX), W2(MXX), WH(MXX), WM(MXX)
C
C                                    Returns zero, assuming 'natural boundary'.
      IF (NX .GT. 1) GOTO 11

C      Q = EXP(TT)
C      AL = ALPI(Q)
C      AMS = AMASS(NQRK)
C      AMH = AMHATF(NQRK)
C      FAC = 2.* LOG (AMS / AMH)
C      FAB = AL / 4.
C
C      CALL INTEGR (NX, 0, Y, W0, IR1)
C      CALL INTEGR (NX, 1, Y, W1, IR2)
C      CALL INTEGR (NX, 2, Y, W2, IR3)
C
   11 CONTINUE
      DO 230 IZ = 1, NX
C                                                       Returns zero answer.
        IF (NX .GT. 1) THEN
        F(IZ) = 0
        GOTO 230
        EndIf
C
C      F(IZ) = FAB * ( FAC * W0(IZ)
C     >            + 2.* (FAC -1.) * (-W1(IZ) + W2(IZ) ) )
C      F(Iz) = fab * ( Fac * (W0(Iz) - Y(Iz) / 2.)
C     >            + 2.* (fac -1.) * (-W1(Iz) + W2(Iz) + Y(Iz) / 12.))
C
  230 CONTINUE
C
      RETURN
C                        ****************************
      END

      subroutine iniread(Nini,qini,anfl,aknl,ordi,blam,Qm4,Qm5)

      implicit double precision (a-h,o-z)

      character*80 lin80

      parameter(NF0=5,Nshp=8)

      Common 
     >  /InpFin/ A(-NF0:NF0), B(-1:Nshp, -NF0:NF0), Ifun, Igfn(-NF0:NF0)
     >  /ShpNfl/ Mshp, Nfval, Nfsea 

      read(Nini,'(a80)')lin80
      read(Nini,'(a80)')lin80
      read(Nini,*)aknl,ordi,anfl,blam,qini,afun,Qm4,Qm5
      
      Ifun = Nint(afun)

      read(Nini,'(a80)')lin80
      read(lin80,*)lin80(1:62), Mshp, Nfval, Nfsea
      read(Nini,'(a80)')lin80
      do ifl=Nfval, Nfsea,-1
         read(Nini,*)Igfn(ifl),B(-1,ifl),A(ifl)
     >        , (B(iex,ifl),iex=1,Mshp)
      enddo
      return
C                      *******************************
      end
      SUBROUTINE INTEGR (NX, M, F,   G, IR)
C                                                   -=-=- integr
C
C     Computes (x/y) ** M * F(y)dy/y integrated over the range [x, 1];
C              Result is G.
C     Integand function F must be defined on an array of size NX which
C              covers the range [0, 1] of the variables x and y.
C     The output function G is returned on the same array.
C     The use of integration variable z defined by
C
C                  y = f(x)
C
C     where f(x) can be any monotonic function.
C
C     IR is an error return code: IR = 1   --  NX out of range
C                                    = 2   --  M  out of range
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER MSG*80
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
C
      DIMENSION   F(NX), G(NX)
C
      DATA IWRN1, IWRN2 / 0, 0 /
C
      IRR = 0
C
      IF (NX .LT. 1 .OR. XA(NX-1,1) .EQ. 0D0) THEN
        MSG = 'NX out of range in INTEGR call'
        CALL WARNI (IWRN1, NWRT, MSG, 'NX', NX, 0, MXX, 0)
        IRR = 1
      EndIf
C
      IF (M .LT. M1 .OR. M .GT. M2) THEN
        MSG ='Exponent M out of range in INTEGR'
        CALL WARNI (IWRN2, NWRT, MSG, 'M', M, M1, M2, 1)
        IRR = 2
      EndIf
C                                                                         NX
      G(NX) = 0D0
C                                                                       NX - 1
      TEM = H(1, NX-1, -M) * F(NX-2) + H(2, NX-1, -M) * F(NX-1)
     >    + H(3, NX-1, -M) * F(NX)
      IF (M .EQ. 0) THEN
         G(NX-1) = TEM
      Else
         G(NX-1) = TEM * XA(NX-1, M)
      EndIf
C                                                                     NX-2 : 2
      DO 10 I = NX-2, 2, -1
         TEM = TEM + H(1,I,-M)*F(I-1) + H(2,I,-M)*F(I)
     >             + H(3,I,-M)*F(I+1) + H(4,I,-M)*F(I+2)
         IF (M .EQ. 0) THEN
            G(I) = TEM
         Else
            G(I) = TEM * XA(I, M)
         EndIf
   10 CONTINUE
C                                                                            1
      TEM = TEM + H(2,1,-M)*F(1) + H(3,1,-M)*F(2) + H(4,1,-M)*F(3)
      IF (M .EQ. 0) THEN
         G(1) = TEM
      Else
         G(1) = TEM * XA(1, M)
      EndIf

      IR = IRR
C
      RETURN
C                        ****************************
      END

      SUBROUTINE KERNEL
     >(XX, FF1, FG1, GF1, GG1, PNSP, PNSM, FF2, FG2, GF2, GG2, NFL, IRT)

C     New version with 'regularized' kernel functions which are smooth, hence
C     are more suitable for interpolation.
C
C     Subroutine to calculate the values of the 1st and 2nd order evolution
C                kernel function at a given value of X.

C     Formulas used are from Furmanski & Petronzio, Phys.Lett. 97B, p437.
C     Notations are different than conventional.  FG <--> GF, NFL * wrong place
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      PARAMETER (D0 = 0.0, D1 = 1.0)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

C
      DATA IWRN /0/
C
      IRT = 0
C
      X = XX
      IF (X .LE. 0. .OR. X .GE. 1.) THEN
        CALL WARNR(IWRN, NWRT, 'X out of range in KERNEL', 'X', X,
     >             D0, D1, 1)
        IRT = 1
        RETURN
      EndIf
C                               1st order kernels with regularization factors,
C                                                   and excluding subtractions
      FF1 = PFF1_pn(x) * (1.- X)
      FG1 = PGF1_pn(x) * X
      GF1 = PFG1_pn(x)
      GG1 = PGG1_pn(x) * X * (1.-X)
C                                                  Begin 2nd order calculations
C                                                   with the non-singlet pieces
cpn 2001
      if (iknl.ge.2) then

cpn                                                         Non-singlet kernels
        PNSP = Fnsp_pn(x)
        PNSM = Fnsm_pn(x)
C                                                 Now the singlet kernel matrix
        FF2 = PFF2_pn(x) 
        FG2 = PGF2_pn(x)
        GF2 = PFG2_pn(x)
        GG2 = PGG2_pn(x)
         

      elseif(iknl.eq.-2) then

cpn                          The signs of PQQB in the expressions for PNSP and
cpn                          PNSM differ in the unpolarized and polarized case 

cpn                    The accuracy of the NLO evolution in the
cpn                     polarized case has not been checked yet, so that
cpn                     the following lines prevent the unknowing
cpn                     running of this part of the code
        print *,'----------------------------------------------------'
        print *,
     >    'The accuracy of the NLO evolution in the polarized case ', 
     >    'was not tested, so I exit. Comment out these lines in ',
     >    'the subroutine KERNEL if you would like to proceed'
        print *,'----------------------------------------------------'
        STOP
 

        PNSP = Fnsm_pn(x)
        PNSM = dFnsm_pn(x)
C                                          Now the singlet kernel matrix
CPN                                   NLO FF kernel without 1/(1-x)_{+} term

        FF2 = dPFF2_pn(x)
cpn                            FG(here)=GF (common notations) and vice versa
        FG2 = dPGF2_pn(x)
        GF2 = dPFG2_pn(x)
CPN                                 NLO GG kernel without 1/(1-x)_{+} piece
        GG2 = dPGG2_pn(x)
      endif  !iknl
C
C                                       Take out singular factors at both ends
      XLG = (LOG(1./(1.-X)) + 1.)
      XG2 = XLG ** 2

      PNSP=PNSP*(1.-x)
      PNSM=PNSM*(1.-x)

      FF2 = FF2 * X * (1.- X)
      FG2 = FG2 * X / XG2
      GF2 = GF2 * X / XG2
      GG2 = GG2 * X * (1.- X)

      RETURN
C                        ****************************
      END!kernel



      SUBROUTINE NSEVL (RHS, IKNL,NX,NT,JT,DT,TIN,NEFF,U0,UN)
C
C                               IKNL determines to which order in Alpha is the
C                                               calculation to be carried out.
C
C       Given the non-singlet parton distribution function U0 at some initial
C       QIN (Tt= 0) in the x interval (0, 1) covered by the array Iz = 1, NX,
C       this routine calculates the evoluted function U at Nt values of Tt at
C       intervals of Dt by numerically integrating the A-P equation using the
C       non-singlet kernel.
C
C       Un(Ix, Tt) = Y(x,es) at the sites Ix=0,..,Nx  (x = Ix * Dz);
C       Un(0, Tt) is obtained by quadratic extrapolation from Ix = 1, 2, 3
C                 for each Tt rather then by evolution because of possible
C                                               singular behavior at x = 0.
C
C       Data is stored at Tt = Is * Dt, Is = 0, ... , Nt.
C                       The function at Is = 0 is the input distribution.
C
C       Numerical iteration is performed with finer grain Ddt = Dt/Jt.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      Common / PdCntrl/ LPrt, LDbg
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
C
      DIMENSION  U0(NX), UN(0:NX, 0:NT)
      DIMENSION  Y0(MXX), Y1(MXX), YP(MXX), F0(MXX), F1(MXX), FP(MXX)

      external rhs
C
         If (Ldbg .Eq. 1) Then
            Write (Nwrt, '(A)') ' Non-Singlet:'
            N5 = Nx / 5 + 1
         Endif

      DDT = DT / JT
C
      IF (NX .GT. MXX) THEN
      WRITE (NOUT,*) 'Nx =', NX, ' greater than Max # of pts in NSEVL.'
      STOP 'Program stopped in NSEVL'
      EndIf
C
C                                       Compute an effective first order lamda
C                                       to be used in checking of moment evl.
      TMD = TIN + DT * NT / 2.
      AMU = EXP(TMD)
      TEM = 6./ (33.- 2.* NEFF) / ALPI(AMU)
      TLAM = TMD - TEM
C
C                       Fill first rows of output array, for initial value of Q
      DO 9 IX = 1, NX
      UN(IX, 0)  = U0(IX)
    9 CONTINUE
      UN(0, 0) = 3D0*U0(1) - 3D0*U0(2) - U0(1)
C                                                                   Initiation
      TT = TIN
      DO 10 IZ = 1, NX
      Y0(IZ)   = U0(IZ)
   10 CONTINUE
C                                                       loop in the Tt variable
      DO 20 IS = 1, NT
C                                                       fine- grained iteration
         DO 202 JS = 1, JT
C                                        Irnd is the counter of the Q-iteration
            IRND = (IS-1) * JT + JS
C                                               Use Runge-Katta the first round
            IF (IRND .EQ. 1) THEN
C
                CALL RHS (TT, Neff, Y0, F0)
C
                DO 250 IZ = 1, NX
                   Y0(IZ) = Y0(IZ) + DDT * F0(IZ)
  250           CONTINUE
C
                TT = TT + DDT
C
                CALL RHS (TT, NEFF, Y0, F1)
C
                DO 251 IZ = 1, NX
                   Y1(IZ) = U0(IZ) + DDT * (F0(IZ) + F1(IZ)) / 2D0
  251           CONTINUE
C                            What follows is a combination of the 2-step method
C                                  plus the Adams Predictor-Corrector Algorithm
            Else
C
                CALL RHS (TT, NEFF, Y1, F1)
C                                                                     Predictor
                DO 252 IZ = 1, NX
                   YP(IZ) = Y1(IZ) + DDT * (3D0 * F1(IZ) - F0(IZ)) / 2D0
  252           CONTINUE
C
C                       Increment of Tt at this place is part of the formalism
                TT = TT + DDT
C
                CALL RHS (TT, NEFF, YP, FP)
C                                                                     Corrector
                DO 253 IZ = 1, NX
                   Y1(IZ) = Y1(IZ) + DDT * (FP(IZ) + F1(IZ)) / 2D0
                   F0(IZ) = F1(IZ)
  253           CONTINUE
            EndIf
C
  202    CONTINUE
C                                                            Fill output array
         DO 260 IZ = 1, NX
            UN (IZ, IS) = Y1(IZ)
  260    CONTINUE
C
C               The value of the function at x=0 is obtained by extrapolation
         UN(0, IS) = 3D0*Y1(1) - 3D0*Y1(2) + Y1(3)
C                                                    Print out for Debugging
      If (LDbg .Eq. 1) Then
         Write (Nwrt, '(A, 5(1pE12.3))') '   :', (Un(Iz,Is), Iz=1,Nx,N5)
      Endif
C
   20 CONTINUE
C
      RETURN
C                        ****************************
      END

      SUBROUTINE NSRHSM (TT, NEFF, FI, FO)
C
C       Subroutine to compute the Right-Side of the Altarelli-Parisi Equation
C                This copy applies to the "NS-minus" piece -- (Qrk - Qrk-bar)
C
C       See comments in NSRHSP for details

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)

      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / KRNL01 / AFF2 (MXX),AFG2 (MXX), AGF2 (MXX), AGG2 (MXX),
     >                  ANSP (MXX), ANSM (MXX), ZGG2, ZFF2, ZQQB

      COMMON / KRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      DIMENSION G1(MXX), FI(NX), FO(NX)
      DIMENSION W0(MXX), W1(MXX), WH(MXX)
C
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = ALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
C
      CALL INTEGR (NX, 0, FI, W0, IR1)
      CALL INTEGR (NX, 1, FI, W1, IR2)
      CALL HINTEG (NX,    FI, WH)
C
      DO 230 IZ = 1, NX
      FO(IZ) = 2.* FI(IZ) + 4./3.* ( 2.* WH(IZ) - W0(IZ) - W1(IZ))
      FO(IZ) = CPL * FO(IZ)
  230 CONTINUE
C                                                  2nd order calculation
cpn 2001      IF (abs(IKNL) .EQ. 2) THEN
      IF (abs(IKNL) .GE. 2) THEN
      DZ = 1./ (NX - 1)
      DO 21 IX = 1, NX-1
        X = XV(IX)
C                                         Number of points in the I-th integral
        NP = NX - IX + 1
        IS = NP
C                                  Evaluate the integrand for the I-th integral
        DO 31 KZ = 2, NP
          IY = IX + KZ - 1
          IT = NX - IY + 1
C                                      XY = X / Y, already calculated in XARRAY
          XY = ZZ (IS, IT)
          G1(KZ) = PNS (IS,IT) * (FI(IY) - XY * FI(IX))
   31   CONTINUE
C
        TEM1 = SMPNOL (NP, DZ, G1, ERR)
C                                                 2nd order contribution
cpn 2001
          if (iknl.ge.2) then 
            TMP2 = (TEM1 - FI(IX) * ANSM(IX)) * CPL2
          elseif (iknl.eq.-2) then
            TMP2 = (TEM1 + FI(IX) * (ZQQB-ANSM(IX))) * CPL2
          endif                 !iknl

C
C                 TWONE = TMP2 / FO(IX)
C                 IF (TWONE .GE. 0.2) THEN
C                 PRINT '(A, 4(1PE15.4))',
C             >  'Second order F big in NSRHS: X, Q, 2nd, 1st = ',X,Q,TMP2,FO(I
C                 EndIf
C                                                  1st + 2nd order terms
        FO(IX) = FO(IX) + TMP2
C
   21 CONTINUE
C
      EndIf

      RETURN
C                        ****************************
      END

      SUBROUTINE NSRHSP (TT, NEFF, FI, FO)
C
C       Subroutine to compute the Right-Side of the Altarelli-Parisi Equation
C       This copy applies to the "NS-plus" piece involving (Qrk + Qrk-bar)
C
cpn 2001
C       IKNL = 1,2,12 1st & 2nd order evolution of the unpolarized case 
C                  in the MS-bar scheme      
C     -1, -2 :  1st and 2nd order of polarized case in the MS-bar scheme
C
C       Nx is the number of mesh-points, Tt is the Log Q variable.
C       Y is the input distribution function defined on the mesh points
C       F is the RHS value (which is also = dY/dt) defined on the same mesh
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)

      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / KRNL01 / AFF2 (MXX),AFG2 (MXX), AGF2 (MXX), AGG2 (MXX),
     >                  ANSP (MXX), ANSM (MXX), ZGG2, ZFF2, ZQQB

      COMMON / KRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      DIMENSION G1(MXX), FI(NX), FO(NX)
      DIMENSION W0(MXX), W1(MXX), WH(MXX)
C
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = ALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
C
      CALL INTEGR (NX, 0, FI, W0, IR1)
      CALL INTEGR (NX, 1, FI, W1, IR2)
      CALL HINTEG (NX,    FI, WH)
C
      DO 230 IZ = 1, NX
      FO(IZ) = 2.* FI(IZ) + 4./3.* ( 2.* WH(IZ) - W0(IZ) - W1(IZ))
      FO(IZ) = CPL * FO(IZ)
  230 CONTINUE
C                                                  2nd order calculation
cpn 2001      IF (abs(IKNL) .EQ. 2) THEN
      IF (abs(IKNL) .GE. 2) THEN
      DZ = 1./ (NX - 1)
      DO 21 IX = 1, NX-1
        X = XV(IX)
C                                        Number of points in the I-th integral
        NP = NX - IX + 1
C                                 Evaluate the integrand for the I-th integral
        DO 31 KZ = 2, NP
          IY = IX + KZ - 1
C                                     XY = X / Y, already calculated in XARRAY
          XY = ZZ (NX-IX+1, NX-IY+1)
          G1(KZ) = PNS (IX,IY) * (FI(IY) - XY * FI(IX))
   31   CONTINUE
C
        TEM1 = SMPNOL (NP, DZ, G1, ERR)
C                                       Assemble 2nd order contribution
        if (iknl.ge.2) then 
          TMP2 = (TEM1 + FI(IX) * (-ANSP(IX) + ZQQB)) * CPL2
        elseif (iknl.eq.-2) then
          TMP2 = (TEM1 - FI(IX) * ANSP(IX)) * CPL2
        endif 
C                                                  Sum 1st + 2nd order terms
        FO(IX) = FO(IX) + TMP2
C
   21 CONTINUE
C
      EndIf

      RETURN
C                        ****************************
      END

      FUNCTION PARDIS (IPRTN, XX, QQ)
C                                                   -=-=- pardis
C
C            (These general comments are enclosed in the
C            lead subprogram in order to survive forsplit.)

C ====================================================================
C GroupName: Pdffns
C Description: User Callable Parton distribution functions
C ListOfFiles: pardis prdis1 xf0 xf1
C ====================================================================
C #Header: /Net/u52/wkt/h22/2evl/RCS/Pdffns.f,v 7.6 99/08/26 23:25:34 wkt Exp $
C #Log:	Pdffns.f,v $
c Revision 7.6  99/08/26  23:25:34  wkt
c pardis.f (and its mirror prdis1) further "refined":
c use x**0.3 as interpolation variable; logistics tightened and more
c fully commented; extraneous calls and calculations eliminated.
c 
c Revision 7.5  99/08/23  23:57:34  wkt
c Fix small bug so that x=1 does not bomb.
c 
c Revision 7.4  99/08/22  22:33:52  wkt
c Small bug fix.
c 
c Revision 7.3  99/08/22  19:05:26  wkt
c See notes below:
c 
C --------------------------------------------------------
C ParDis refined to incorporate all the following features:
C (1) Smoothness: f(x,Q) is continuous and has continuous derivatives,
C         -- by applying 4-point (cubic) interpolation, 
C            and maintaining x in the interior interval;
C (2) Speed: optimized, so that it is even faster than the original,
C            by using Jon Pumplin's "in-line" code for interior intervals
C            rather than subroutine calls as before;
C            (Note; the original code reproduces Jon's results if the const.   
C            M is set = 3, rather than the default 2. It is just slower.)
C (3) End-interval handling for both interpolation and extrapolation:
C            Always use 4-point interpolation to ensure maximum reliability;
C            here I do use the more robust PolInt which allow the variable
C            (either x or Q) to be anywhere w.r.t. the grid points.
C            The cost in time is minimal, since end-intervals are rarely used.
C (4) Extrapolation for the interval 0 < x < Xmin:
C            Use x^2 *f(x,Q) for interpolation; so that x=0 can be included
C            as a regulating end point. This gives much more stable and 
C            "reasonable" results, compared to using only pts with x C= Xmin;
C (5) Same 4-pt extrapolation is applied to the regions Q < Qini, Q C Qmax;
C            This yields "reasonable" results, since Q-dep. is smooth
C            in t = Ln Ln Q, the natural evolution and interpolation variable
C            and it is not oscillatory.
C Also, the logic of the code has been much cleaned up.  The program should 
C be self-explanary.
C
C This code has been checked to reproduce the interior point results of Jon's
C program exactly.  Extrapolations in all directions in x and Q have also been
C tried and examined. Both the warning messages and the results appear 
C reasonable.
C ---------------------------------------------------------
c Revision 7.2  99/08/20  23:15:33  wkt
c X1ARAY common block synchronized.
c 
c Revision 7.1  99/08/19  17:36:58  wkt
c 1. Interpolation routine Pardis.f and Prdis1.f replaced by Jon Pumplin's
c improved (4-point interpolation) version with continuous derivatives.
c Speed is comparable to before, but crashes when arguments are outside
c the strict ranges of evolution.

c 2. Zf0(Ix) and Zf1(Ix) functions added in modules Xf0 and Xf1 to allow
c recall of the Z-array used in the evolution, just like for the X-array.
c 
c Revision 7.0  99/07/25  18:45:44  wkt
c Pavel's improved version: No real change. (Some cosmatic)
c 
c Revision 6.2  97/12/21  21:46:57  wkt
c KF-2 --> Nfmx  in  alternative module
c
c Revision 6.1  97/11/15  17:19:55  wkt
c v.6.0 + HLL revisions
c
C ====================================================================
C      FUNCTION PARDIS(IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find 
c  the parton distribution at an arbitray point in x and q.
c
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Character Msg*80
      LOGICAL LSTX
 
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      PARAMETER (Smll = 1D-9)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX

      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX), KF, Nelmt
 
      dimension fvec(4), fij(4)
      dimension xvpow(0:mxx)
      Data Iwrn1, Iwrn2, Iwrn3, OneP / 3*0, 1.00001 /

      parameter(xpow =  1.d0/3.d0)      !**** new choice of interpolation variable
c     parameter(xpow =  0.3d0)          !**** old choice of interpolation variable

      data nxsave / -999 /

      parameter(nqvec = 4)

      save nxsave, xvpow
      save xlast, qlast
      save jq, jx, JLx, JLq, ss, sy2, sy3, s23, ty2, ty3
      save const1 , const2, const3, const4, const5, const6
      save tt, t13, t12, t23, t34 , t24, tmp1, tmp2, tdet

c store the powers used for interpolation on first call or if nx has changed...
      if(nx .ne. nxsave) then
         nxsave = nx

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ

c Enforce thresholds immediately to save calculation time.
      if(q .lt. amass(iprtn)) then
         pardis = 0.d0
         return
      endif

c skip the initialization in x if same as in the previous call.
      if(x .eq. xlast) goto 100
      xlast = x

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1   
C
      If     (JLx .LE. -1) Then
!        Print '(A,1pE12.4)', 'Severe error: x <= 0 in ParDis! x = ', x 
!        Stop 
	ParDis=0
	return
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
         Msg = '0 < X < Xmin in ParDis; extrapolation used!'
         CALL WARNR (IWRN1, NWRT, Msg, 'X', X, Xmin, 1D0, 1)
      Elseif (JLx .LE. Nx-2) Then

C                For interior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
!        Print '(A,1pE12.4)', 'Severe error: x > 1 in ParDis! x = ', x 
!        Stop 
	ParDis=0
	return
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2 
      sy3 = ss - svec3 

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12 
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

100   continue

c skip the initialization in q if same as in the previous call.
      if(q .eq. qlast) goto 110
      qlast = q

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (Q .GE. QV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
         If (JLq .LT. 0) Then
          Msg = 'Q < Q0 in ParDis; extrapolation used!'
          CALL WARNR (IWRN2, NWRT, Msg, 'Q', Q, Qini, 1D0, 1)
         EndIf
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3
        If (JLq .GE. Nt) Then
         Msg = 'Q > Qmax in ParDis; extrapolation used!'
         CALL WARNR (IWRN3, NWRT, Msg, 'Q', Q, Qmax, 1D0, 1)
        Endif
      Endif
C                                   This is the interpolation variable in Q
      tt = log(log(q/al))

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf

110   continue

c get the pdf function values at the lattice points...

      jtmp = ((IPRTN + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1) 

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * Xa(1,2)
         fij(3) = Upd(J1+2) * Xa(2,2)
         fij(4) = Upd(J1+3) * Xa(3,2)
C 
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4 (XVpow(0), Fij(1), 4, ss, Fx, Dfx) 
         
         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 
C                                              Pdf is undefinged for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:
C        Jon's algorithm does not handle this case, give it to Polint also.

        Call Polint4 (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)

        Fvec(it) = Fx

       Else 
C                       for all interior points, use Jon's in-line function 
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)          
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         Fvec(it) = (const5*(Upd(J1)   - sf2*const1 + sf3*const2) 
     &             + const6*(Upd(J1+3) + sf2*const3 - sf3*const4) 
     &             + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...
C              Again, use Polint for end intervals, Jon for the interior

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4 (TV(0), Fvec(1), 4, tt, ff, Dfq)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4 (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
C                         Interior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 
     &      +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PARDIS = ff

      Return
C      				********************
      End

      SUBROUTINE PARPDF (IACT, NAME, VALUE, IRET)
C                                                   -=-=- parpdf

C            (These general comments are enclosed in the
C            lead subprogram in order to survive forsplit.)

C ====================================================================
C GroupName: Parpdf
C Description: routines to input and output evolution parameters
C ListOfFiles: parpdf evlpar evlget evlgt1 evlin evlset
C ====================================================================
C #Header: /Net/u52/wkt/h22/2evl/RCS/Parpdf.f,v 7.2 99/08/20 23:15:25 wkt Exp $
C #Log:	Parpdf.f,v $
c Revision 7.2  99/08/20  23:15:25  wkt
c X1ARAY common block synchronized.
c 
c Revision 7.1  99/08/20  09:53:37  wkt
c New variable NuIni (Unit # for .ini file) added to facilitate the use
c of user-defined .ini file to perform evolution (Iset = Isetin2 = 902).
c 
c Revision 7.0  99/07/25  18:44:26  wkt
c Pavel's improved version: No real change
c 
c Revision 6.2  97/11/16  00:28:18  wkt
c evlget and evlset modified to accommodate IpdMod & Iptn0
c Strictly speaking, these belong to the new pdfpac, not
c in this evlpac. But the former does not have its own
c input/output routine. Thus, the common /PdfSwh/ is shared
c in violation of our usual rules.
c
c Revision 6.1  97/11/15  17:19:48  wkt
c v.6.0 + HLL revisions
c
C ====================================================================

C      SUBROUTINE PARPDF (IACT, NAME, VALUE, IRET)
C
C               Actions: 0      type list of variables on unit VALUE.
C                        1      set variable with name NAME to VALUE, if
C                               it exists, Else set IRET to 0.
C                        2      find value of variable.  If it does not exist,
C                               set IRET to 0.
C                        3      request values of all parameters from terminal.
C                        4      type list of all values on unit VALUE.
c
C               IRET =   0      variable not found.
C                        1      successful search.
C                        2      variable found, but bad value.
C                        3      bad value for IACT.
C                        4      no variable search (i.e., IACT is 0, 3, or 4).
C
C               If necessary, VALUE is converted to integer by NINT
C
C              Use ILEVEL and START1 to start search for variable names close
C              to previous name to ensure effiency when reading in many values.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      CHARACTER NAME*(*), Uname*10
C
      LOGICAL START1
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA ILEVEL, LRET / 1, 1 /
C
      JRET = IRET
      CALL UPC (NAME, Ln, Uname)
      IF (IACT .EQ. 0 .OR. IACT .EQ. 4)
     >    IVALUE = NINT (VALUE)
      START1 = (IACT .NE. 1) .AND. (IACT .NE. 2)
      IF (START1)  ILEVEL = 1
C
      GOTO (1, 2), ILEVEL
C
    1 START1 = .TRUE.
      ILEVEL = 0
        print *, "inside call ParPDF"
      CALL PARQCD (IACT, Uname(1:Ln), VALUE, JRET)
              IF (JRET .EQ. 1)  GOTO 11
              IF (JRET .EQ. 2)  GOTO 12
              IF (JRET .EQ. 3)  GOTO 13
              IF (JRET .GT. 4)  GOTO 15
              ILEVEL =  ILEVEL + 1
    2 CALL EVLPAR (IACT, Uname(1:Ln), VALUE, JRET)
              IF (JRET .EQ. 1)  GOTO 11
              IF (JRET .EQ. 2)  GOTO 12
              IF (JRET .EQ. 3)  GOTO 13
              IF (JRET .GT. 4)  GOTO 15
              ILEVEL =  ILEVEL + 1
C
      IF (.NOT. START1) GOTO 1
C
      IF (JRET .EQ. 0)  GOTO 10
C
C                       Arrive here if IACT = 0, 3 and all is OK (IRET=4):
    9 CONTINUE
C                       WRITE (IVALUE, 100)
      GOTO 14
C                       Exits:
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
   13 CONTINUE
   14 CONTINUE
   15 CONTINUE
C                                   LRET is used for debugging purpose only
      IF (JRET .NE. 4) LRET = JRET
      IF (LRET.EQ.0 .OR. LRET.EQ.2 .OR. LRET.EQ.3) THEN
        PRINT *, 'Error in PARPDF: IRET, IACT, NAME, VALUE =',
     >  LRET, IACT, NAME, VALUE
        PAUSE
      EndIf
      IRET= JRET
      RETURN
C
  100 FORMAT (/)
C                        ****************************
      END


      SUBROUTINE PDFSET0 (I1, IHDRN, ALAM, TPMS, QINI, QMAX, XMIN,
C                                                   -=-=- pdfset
     >                   NU, HEADER, I2, I3, IRET, IRR, UserIn)
C
C                      ------------------------------
C       Arguments:
C
C               >  signifies an input  parameter,
C               <  signifies an output parameter,
C               *  indicates that the parameter is Input      if  I1=0 ;
C                                                  Output     if  I1=1 ;
C                                                  Inactive   if  I1=2 .
C
C       >       I1    =   0  if evolution calculation is to be carried out,
C                       1  if data table  is to be read from an existing file;
C                            (NU must be set to the Unit # of the file)
C                        -1  same as 1, except data is read to alternate
C                            common blocks for "second set" of PDF's
C                            (NU must be set to the Unit # of the file.)
C                         2  if data table in memory is to be written to an
C                            external file with no calculation (NU must be
C                            specified)
C
C       *       Ihdrn =   hadron ID code:  proton is 1,  anti-proton is -1;
C
C       *       Alam  =   QCD 'Lamda parameter' in the 'CWZ' scheme
C                         (See Collins & Tung 'Calculating heavy ... ',
C                          Fermilab-Pub-86/39-T for definition and
C                          conversion graph to '1st order' and 'MS-bar'
C                          Lamda values.)
C
C       *       Tpms  =   Mass of the Top-quark in GeV;
C
C       *       Qini  =   Q-initial (GeV), i.e. starting value of Q-evolution;
C
C       *       Qmax  =   Maxmum value of Q (GeV) to be used;
C
C       *       Xmin  =   Smallest value of X desired;
C
C       >       NU    =   Unit # to (from) which data is to be written (read);
C
C       *       Header=   Informational header for the file
C
C       >       I2    (if I1 = 0 is chosen):  code to choose initial distribu-
C                     tion functions at Q = Qini  (Ignored if I1 = 1, 2, -1)
C
C                     =  0  :     Use initial distribution function
C                                 UserIn (Iparton, X) to be supplied by user;
C                                 UserIn must be declared EXTERNAL in the
C                                 calling program.
C
C                  1 - 99:     Use "canned" parton distribution sets as input.
C                                 The list of parametrizations is updated
C                             periodically; it is given in the current version
C                                 of the PDF function.
C
C       >       I3    (if I1 = 0, 2 is chosen):  switch to write to file
C
C                     =  0  :     results of calculation not written to a file;
C
C                        1  :     results written to the file (with unit # NU)
C
C                      2  :     results written to the alternate common blocks
C
C       <       Iret        :     return error code
C
C                     =  0  :     no error
C
C                     =  1  :     Error opening file to read (check Irr)
C
C                     =  2  :     Error condition in reading file (check Irr)
C
C                     =  10 :     Xmin > 1. is not allowed
C
C                     =  11 :     Xmin (< 1E-7) is too small for use of the
C                                 current formalism; Default Nx, Ke kept
C
C                     =  20 :     I2 value is illegal
C
C                     =  22 :     All parameters are unchanged from previous
C                                 run; calculation by-passed
C
C                     =  30 :     I3 < 0 .OR. > 2 is illegal
C
C                     =  31 :     Error condition in opening file for write
C
C                     =  32 :     Write error from EVLWT  (check Irr value)
C
C                     =  33 :     Read error from EVLRD1  (check Irr value)
C
C                     =  40 :     I1 < 0 .or. > 1 is illegal
C
C                     =  98 :     Evolution calculation by-passed because
C                                 all relevant parameters remained the same
C                                 as before; data table in memory should
C                                 remain valid.
C
C       <       Irr   =  error code returned by Fortran when error condition
C                                 is encountered in reading from or writing
C                                 to files (check Fortran manual).
C
C                       ----------------------------
C
C       Additional parameters which can be 'dialed' by the user but which do
C                             appear in the above argument list because they
C                             are less likely to be changed:
C
C       Physics parameters:
C
C               NFL     is the number of quark flavors
C                       legal values: 0, 1, ... 6 ; default = 6;
C
C               M1, M2, M3, M4, M5  are masses of quarks u, d, s, c, b respec.
C
C              IKNL selects the evolution kernel on the RHS of the AP equation
C                  = 1 : 1st order unpolarized case
C                    2 : 2nd order unpolarized case in the MS-bar scheme
C                   12 : 2nd order unpolarized case in the DIS scheme
C               -1, -2 :      ......       polarized .. in the MS-bar scheme
C       GRID points in the variables X and T (=Log Q):
C
C               NX  (< MXX),  default adjusted according to choise of Xmin;
C               NT  (< MXQ),
C               (cf. program codes Parameter ... for values of Mxx and Mxq)
C
C              The (integer) parameter JT (with default set at 1) subdivides
C              the step-size of the variable Q in the solution of the evolu-
C              tion equation.  Using JT > 1 increases the numerical precision
C              of the results without increasing the size of the data block,
C              but, of course, at the expense of multiplying the execution
C              time.
C
C              A similar parameter (JX) will be added in a future version to
C              control the X-mesh points used in actual calulation,
C              independent of the declared array size.  This will allow us
C              to reduce the size of the data block without sacrificing
C              much numerical accuracy -- a desirable change for memory-
C              constrained applications.
C
C       All parameters listed in this section can be changed ('set') by the
C       following fortran call to a general purpose input/output routine:
C
C               Call ParPdf (1, 'NAME', VALUE, Iret)
C
C      prior to invoking the PDFSET subroutine.
C       Here NAME is the (character) name of the variable as listed above,
C       VALUE is the (real) value of the variable to be set, and Iret is an
C       error return code (e.g. Iret=0 means 'variable NAME cannot be found').
C
C       Likewise, the current value of any variable can be recalled by the
C       following call:
C
C               ParPdf (2, 'NAME', VALUE, Iret)
C
C      For other use of PARPDF and further details consult the program listing
C      in the section EVLPAR inside PDFIT.FOR.
C
C     For interactive use of these programs:
C     The entire list of adjustable parameters can be brought to the screen
C     and systematically changed in an interactive way by the Fortran call:
C
C               Call ChgEvl (InputSet)
C
C     where InputSet is an optional "setup subroutine" for any input initial
C     parton distribution function (at Q0) (denoted by UserIn in the above
C     argument list for PdfSet)   which needs to be called before the
C actual function evaluations (setting up common block of coefficients...etc.)
C     If needed, InputSet must be declared External in the calling program.
C     If not needed, it can be ignored.
C     The user will be prompted for possible inputs and other options.
C
C                      ----------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      CHARACTER HEADER*78, FORMT*11
C
      LOGICAL LWT, LOP

      Common / IOUnit / Nin, Nout, Nwrt
C
      EXTERNAL FPIN, UserIn
C
      IRET = 0
      IRR  = 0
      LWT = .FALSE.
C                                                 ----------------------------
C                                                       Read only section:
      IF (I1 .EQ. 1) THEN
C        OPEN (NU, FORM='FORMATTED', STATUS='OLD',IOSTAT=IRR,ERR=10)
         INQUIRE (NU, OPENED=LOP, FORM=FORMT)
         IF (.NOT. LOP)
     >      STOP 'Data table file not openned for read in PDFSET'
C                                                                     Read file
         CALL EVLRD (NU, HEADER, IRR)
         IF (IRR .NE. 0) THEN
            IRET = 2
            RETURN
         EndIf
C               Recall physical parameters used in generating data in the file
C
         CALL PARPDF (2, 'ALAM', ALAM, J1)
         CALL PARPDF (2, 'IHDN', HDRN, J2)
         CALL PARPDF (2, 'M6  ', TPMS, J3)
         CALL PARPDF (2, 'QINI', QINI, J4)
         CALL PARPDF (2, 'QMAX', QMAX, J5)
         CALL PARPDF (2, 'XMIN', XMIN, J6)
         IHDRN = NINT(HDRN)
C                                      Add Error testing codes (Js) if desired
   10    IF (IRR .NE. 0) IRET = 1
C
      ElseIF (I1 .EQ. -1) THEN
C
         INQUIRE (NU, OPENED=LOP, FORM=FORMT)
         IF (.NOT. LOP)
     >      STOP 'Data table file not openned for read in PDFSET'
C                                                                     Read file
         CALL EVLRD1 (NU, HEADER, IRR)
         IF (IRR .NE. 0) THEN
            IRET = 33
            RETURN
         EndIf
C               Recall physical parameters used in generating data in the file
C
         CALL PARPDF (2, 'ALAM', ALAM, J1)
         CALL PARPDF (2, 'IHDN', HDRN, J2)
         CALL PARPDF (2, 'M6  ', TPMS, J3)
         CALL PARPDF (2, 'QINI', QINI, J4)
         CALL PARPDF (2, 'QMAX', QMAX, J5)
         CALL PARPDF (2, 'XMIN', XMIN, J6)
         IHDRN = NINT(HDRN)
C                                      Add Error testing codes (Js) if desired
   12    IF (IRR .NE. 0) IRET = 1
C                                               ----------------------------
C                                               Main section to do evolution:
C                                               Set parameters for evolution
      ElseIF (I1 .EQ. 0) THEN
C
         HDRN = DBLE(IHDRN)
         CALL PARPDF (1, 'ALAM', ALAM, J1)
         CALL PARPDF (1, 'IHDN', HDRN, J2)
         CALL PARPDF (1, 'M6  ', TPMS, J3)
         CALL PARPDF (1, 'QINI', QINI, J4)
         CALL PARPDF (1, 'QMAX', QMAX, J5)
         CALL PARPDF (1, 'XMIN', XMIN, J5)
C                                            Add Error testing code if desired
C
         IF     (XMIN .GT. 1D0) THEN
            IRET = 10
            RETURN
         EndIf
C
C       Call evolution routine according to the choice of initial distribution
C
         IF     (I2 .LT. 0) THEN
            IRET = 20
            RETURN
         ElseIF (I2 .EQ. 0) THEN
            AI2 = 0.
            CALL PARPDF (1, 'IPD0', AI2, L1)
            IRR = 99
            Aout = Nout
            Print *, 'Start evolution with parameters:'
            CALL PARPDF (4, 'ALL', Aout, L1)
            CALL EVOLVE (UserIn, IRR)
C                                                       Add error checking here
         ElseIF (I2 .LE. 99) THEN
            AI2 = I2
            CALL PARPDF (1, 'IPD0', AI2, L1)
            CALL EVOLVE (FPIN, IRR)
            IF (IRR .EQ. 98) IRET = 98
         Else
            IRET = 20
            PRINT *, 'I2 out of range in PDFSET; I2 =', I2
            RETURN
         EndIf
         IF (IRR .EQ. 98) IRET = 22
C                                                 Write results to file if I3=1
         IF (I3 .LT. 0 .OR. I3 .GT. 2) THEN
            IRET = 30
         ElseIF (I3 .GT. 0) THEN
            LWT  = .TRUE.
         EndIf
C                                               ----------------------------
      ElseIF (I1 .EQ. 2) THEN
         LWT = .TRUE.
      Else
         IRET = 40
      EndIf
C
      IF (LWT) THEN

         INQUIRE (NU, OPENED=LOP, FORM=FORMT)
         IF (.NOT. LOP)
     >      STOP 'Data table file not openned for write in PDFSET'
C
        CALL EVLWT (NU, HEADER, IRR)
        IF (IRR .NE. 0) THEN
           PRINT *, 'Error Writing Data to File in PDFSET'
           IRET = 32
        EndIf
        REWIND (NU)
C
        IF (I3 .EQ. 2) THEN

         INQUIRE (NU, OPENED=LOP, FORM=FORMT)
         IF (.NOT. LOP)
     >      STOP 'Data table file not openned for read in PDFSET'
          CALL EVLRD1 (NU, HEADER, IRR)
          IF (IRR .NE. 0) THEN
            PRINT *,
     >      'Error Reading Data to the Alternate Block in PDFSET'
            IRET = 33
          EndIf
        EndIf

   30       IF (IRR .NE. 0) IRET = 31
      EndIf
C
      RETURN
C                        ****************************
      END

C                                                          =-=-= Pdffns
      FUNCTION PFF1 (XX)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LA, LB, LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      PARAMETER (MX = 3)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / KRNL00 / DZ, XL(MX), NNX
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / KRN1ST / FF1(0:MXX), FG1(0:MXX), GF1(0:MXX), GG1(0:MXX),
     >                  FF2(0:MXX), FG2(0:MXX), GF2(0:MXX), GG2(0:MXX),
     >                  PNSP(0:MXX), PNSM(0:MXX)

      SAVE
      DATA LA, LB / 2 * .FALSE. /
C
      LB = .TRUE.
      ENTRY TFF1(ZZ)
      LA = .TRUE.
      GOTO 2

      ENTRY QFF1 (XX)
      LB = .TRUE.
      ENTRY RFF1 (XX)
C
    2 IF (LA .AND. .NOT.LB) THEN
        Z = ZZ
        X = XFRMZ (Z)
      Else
        X = XX
      EndIf
C                            If Xmin < X < 1. interpolate in Z (equally spaced)
C                            If    0 < X < Xmin extrapolate in X (finite range)
      IF (X .GE. D1) THEN
        PFF1 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (FF1,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, FF1(1), MX, X, TEM, ERR)
      EndIf
C
      IF (LA) THEN
         IF (LB) THEN
            PFF1 = TEM / (1.-X)
            LB   =.FALSE.
         Else
            TFF1 = TEM / (1.-X) * DXDZ(Z)
         EndIf
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QFF1 = TEM
            LB   =.FALSE.
         Else
            RFF1 = TEM * X / (1.-X)
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY FNSP (XX)
      X = XX
      IF (X .GE. D1) THEN
        FNSP = 0.
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (PNSP,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, PNSP(1), MX, X, TEM, ERR)
      EndIf

      FNSP = TEM / (1.- X)
      RETURN
C                        ----------------------------
      ENTRY FNSM (XX)
      X = XX
      IF (X .GE. D1) THEN
        FNSM = 0.
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (PNSM,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, PNSM(1), MX, X, TEM, ERR)
      EndIf

      FNSM = TEM / (1.- X)
      RETURN
C                        ----------------------------
      ENTRY PFG1 (XX)
      LA = .TRUE.
      ENTRY QFG1 (XX)
      LB = .TRUE.
      ENTRY RFG1 (XX)
      X = XX
      IF (X .GE. D1) THEN
         PFG1 = 0
         RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (FG1,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, FG1(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PFG1 = TEM
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QFG1 = TEM * (1.- X)
            LB   =.FALSE.
         Else
            RFG1 = TEM * X
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY PGF1 (XX)
      LA = .TRUE.
      ENTRY QGF1 (XX)
      LB = .TRUE.
      ENTRY RGF1 (XX)
      X = XX
      IF (X .GE. D1) THEN
        PGF1= 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (GF1,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, GF1(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PGF1 = TEM / X
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QGF1 = TEM * (1.- X) / X
            LB   =.FALSE.
         Else
            RGF1 = TEM
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY PGG1 (XX)
      LA = .TRUE.
      ENTRY QGG1 (XX)
      LB = .TRUE.
      ENTRY RGG1 (XX)
      X = XX
      IF (X .GE. D1) THEN
        PGG1= 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (GG1,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, GG1(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PGG1 = TEM / X / (1.-X)
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QGG1 = TEM / X
            LB   =.FALSE.
         Else
            RGG1 = TEM / (1.-X)
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY PFF2 (XX)
      LA = .TRUE.
      ENTRY QFF2 (XX)
      LB = .TRUE.
      ENTRY RFF2 (XX)
      X = XX
C      XLG = LOG(X/(1.-X))**2
      IF (X .GE. D1) THEN
        PFF2 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (FF2,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, FF2(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PFF2 = TEM / X / (1.-X)
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QFF2 = TEM / X
            LB   =.FALSE.
         Else
            RFF2 = TEM / (1.-X)
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY PFG2 (XX)
      LA = .TRUE.
      ENTRY QFG2 (XX)
      LB = .TRUE.
      ENTRY RFG2 (XX)
      X = XX
      XM1 = 1.- X
      XLG = (LOG(1./XM1) + 1.)**2
      IF (X .GE. D1) THEN
        PFG2 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (FG2,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, FG2(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PFG2 = TEM / X * XLG
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QFG2 = TEM / X * XLG * XM1
            LB   =.FALSE.
         Else
            RFG2 = TEM * XLG
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY PGF2 (XX)
      LA = .TRUE.
      ENTRY QGF2 (XX)
      LB = .TRUE.
      ENTRY RGF2 (XX)
      X = XX
      XM1 = 1.- X
      XLG = (LOG(1./XM1) + 1.) ** 2
      IF (X .GE. D1) THEN
        PGF2 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (GF2,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, GF2(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PGF2 = TEM / X * XLG
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QGF2 = TEM / X * XLG * XM1
            LB   =.FALSE.
         Else
            RGF2 = TEM * XLG
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY PGG2 (XX)
      LA = .TRUE.
      ENTRY QGG2 (XX)
      LB = .TRUE.
      ENTRY RGG2 (XX)
      X = XX
C      XLG = LOG(X/(1.-X)) ** 2
      IF (X .GE. D1) THEN
        PGG2 = 0
        RETURN
      ElseIF (X .GE. XMIN) THEN
        Z = ZFRMX (X)
        TEM = FINTRP (GG2,  -DZ, DZ, NX,  Z, 1, ERR, IRT)
      Else
        CALL POLINT (XL, GG2(1), MX, X, TEM, ERR)
      EndIf
      IF (LA) THEN
         PGG2 = TEM / X / (1.-X)
         LA   =.FALSE.
      Else
         IF (LB) THEN
            QGG2 = TEM / X
            LB   =.FALSE.
         Else
            RGG2 = TEM / (1.-X)
         EndIf
      EndIf
      RETURN
C                        ----------------------------
      ENTRY VFF1 (XX)
      X = XX
      TEM = (1.+ X**2) / (1.- X)
      VFF1 = TEM *4./3.
      RETURN
C                        ****************************
      END
C                                                          =-=-= Evlaux
      function PFF1_pn(x)
cpn                                                       Returns P_{FF}^{(0)} 
      implicit NONE
      double precision PFF1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PFF1: x=',x,' is out of range'

      PFF1_pn=CF*(1.+x**2)/(1.-x) 


      Return
      End!PFF1_pn

cpn-------------------------------------------------------------------------


      function PFF2_pn(x)
cpn                    Returns P_{FF}^{(2)}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision PFF2_pn, x
      double precision  FFCF2, FFCFG,FFCFT,FFP,FFM,x2,xln1m,xln1m2,xln,
     >     xln2,spen2,spenc2,tem,xi,xm1i,xp1i

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PFF2: x=',x,' is out of range'

      XI = 1./ X
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      FFP = (1.+ X2) * XM1I
      FFM = (1.+ X2) * XP1I
 
      FFCF2 = - 1. + X + (1.- 3.*X) * XLN / 2. - (1.+ X) * XLN2 / 2.
     >     - FFP * (3.* XLN / 2. + 2.* XLN * XLN1M)
     >     + FFM * 2.* SPEN2
      FFCFG = 14./3.* (1.-X)
     >     + FFP * (11./6.* XLN + XLN2 / 2. + 67./18. - PI2 / 6.)
     >     - FFM * SPEN2
      FFCFT = - 16./3. + 40./3.* X + (10.* X + 16./3.* X2 + 2.) * XLN
     >     - 112./9.* X2 + 40./9./X - 2.* (1.+ X) * XLN2
     >     - FFP * (10./9. + 2./3. * XLN)

      tem = CFTR*NFlv * FFCFT + CF2 * FFCF2 + CFCG   * FFCFG

      PFF2_pn=tem

      Return
      End!PFF2_pn

cpn---------------------------------------------------------------------


      function PFG1_pn(x)
cpn                               Returns P_{FG}^{(0)}(conventional notations)
      implicit NONE
      double precision PFG1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PFG1: x=',x,' is out of range'

      PFG1_pn=2*Tr*Nflv *(1. - 2.* X + 2.* X**2)

      Return
      End!PGF1_pn

cpn-------------------------------------------------------------------------


      function PFG2_pn(x)
cpn                    Returns P_{FG}^{(2)} (conventional notations}

      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision PFG2_pn, x
      double precision  FGCFT,FGCGT,FGP,FGM,x2,xln1m,xln1m2,xln,xln2,
     >    spen2,spenc2,tem

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGF2: x=',x,' is out of range'

      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      FGP = 1. - 2.* X + 2.* X2
      FGM = 1. + 2.* X + 2.* X2
 
      FGCFT = 4.- 9.*X + (-1.+ 4.*X)*XLN + (-1.+ 2.*X)*XLN2 + 4.*XLN1M
     >     + FGP * (-4.*XLN*XLN1M + 4.*XLN + 2.*XLN2 - 4.*XLN1M
     >     + 2.*XLN1M2 - 2./3.* PI2 + 10.)
      FGCGT = 182./9.+ 14./9.*X + 40./9./X + (136./3.*X - 38./3.)*XLN
     >     - 4.*XLN1M - (2.+ 8.*X)*XLN2
     >     + FGP * (-XLN2 + 44./3.*XLN - 2.*XLN1M2 + 4.*XLN1M
     >     + PI2/3. - 218./9.)
     >     + FGM * 2. * SPEN2

      tem = CFTR*NFlv * FGCFT               + CGTR*nflv * FGCGT
      PFG2_pn=tem

      Return
      End!PFG2_pn

cpn----------------------------------------------------------------------------
      function PGF1_pn(x)
cpn                               Returns P_{GF}^{(0)}(conventional notations)
      implicit NONE
      double precision PGF1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg


      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGF1: x=',x,' is out of range'

      PGF1_pn=CF*(2.- 2.* X + X**2) / X

      Return
      End!PGF1_pn

cpn-------------------------------------------------------------------------

      function PGF2_pn(x)
cpn                    Returns P_{GF}^{(2)} (conventional notations}
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision PGF2_pn, x
      double precision  GFCF2, GFCFG,GFCFT,GFP,GFM,x2,xln1m,xln1m2,xln,
     >    xln2,spen2,spenc2,tem

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGF2: x=',x,' is out of range'


      X2 = X ** 2
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      GFP = (2.- 2.* X + X2) / X
      GFM = - (2.+ 2.* X + X2) / X
 
      GFCF2 = - 5./2.- 7./2.* X + (2.+ 7./2.* X) * XLN + (X/2.-1.)
     >     *XLN2- 2.* X * XLN1M- GFP * (3.* XLN1M + XLN1M ** 2)
      GFCFG = 28./9. + 65./18.* X + 44./9. * X2 - (12.+ 5.*X + 8./3.
     >     *X2)* XLN + (4.+ X) * XLN2 + 2.* X * XLN1M+ GFP * (-2.*XLN
     >     *XLN1M + XLN2/2. + 11./3.*XLN1M + XLN1M2- PI2/6. + 0.5)+
     >     GFM * SPEN2
      GFCFT = -4./3.* X - GFP * (20./9.+ 4./3.*XLN1M)


      tem= CFTR*NFlv * GFCFT + CF2 * GFCF2 + CFCG   * GFCFG
      PGF2_pn=tem

      Return
      End!PGF2_pn

cpn----------------------------------------------------------------------------
      function PGG1_pn(x)
cpn                                                       Returns P_{GG}^{(0)} 
      implicit NONE
      double precision PGG1_pn, x

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGG1: x=',x,' is out of range'

      PGG1_pn=2.*CG *(1./(1.-X) + 1./X - 2. + X - X**2)

      Return
      End!PGG1_pn

cpn---------------------------------------------------------------------------

      function PGG2_pn(x)
cpn                    Returns P_{GG}^{(2)}(x) 
      implicit NONE
      real pi,pi2
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)
      double precision PGG2_pn, x
      double precision  GGCFT,GGCGT,GGCG2,GGP,GGM,x2,xln1m,xln1m2,xln,
     >    xln2,spen2,spenc2,tem,xm1i,xp1i,xi

      external Spenc2 
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 
      integer nflv

      if (x.lt.0d0.or.x.gt.1d0) 
     >  print *,'Warning in the function PGG2: x=',x,' is out of range'

      XI = 1./ X
      X2 = X ** 2
      XM1I = 1./ (1.- X)
      XP1I = 1./ (1.+ X)
      XLN = LOG (X)
      XLN2 = XLN ** 2
      XLN1M = LOG (1.- X)
      xln1m2=xln1m**2
      SPEN2 = SPENC2 (X)

      GGP = XM1I + XI - 2. + X - X2
      GGM = XP1I - XI - 2. - X - X2

      GGCFT = -16.+ 8.*X + 20./3.*X2 + 4./3./X + (-6.-10.*X)*XLN
     >     - 2.* (1.+ X) * XLN2
      GGCGT = 2.- 2.*X + 26./9.*X2 - 26./9./X - 4./3.*(1.+X)*XLN
     >     - GGP * 20./9.
      GGCG2 = 27./2.*(1.-X) + 67./9.*(X2-XI) + 4.*(1.+X)*XLN2
     >     + (-25.+ 11.*X - 44.*X2)/3.*XLN
     >     + GGP * (67./9.- 4.*XLN*XLN1M + XLN2 - PI2/3.)
     >     + GGM * 2.* SPEN2
 
      tem = CFTR*NFlv * GGCFT + CG2 * GGCG2 + CGTR*NFlv * GGCGT

      PGG2_pn=tem

      Return
      End!PGG2_pn

cpn-------------------------------------------------------------------------

      FUNCTION PMSDIS (ISET, IHDRN, IPRTN, XX, QQ, IAT)
C                                                   -=-=- pmsdis

C            (These general comments are enclosed in the
C            lead subprogram in order to survive forsplit.)

C ====================================================================
C GroupName: Pdfmis
C Description: Misc.:DIS scheme <--> MS-bar scheme conversion
C non-essential x <--> z conversion: x2zconv
C Calculate moments: amomen
C ListOfFiles: pmsdis c2qx sqrk amomen x2zconv
C ====================================================================
C #Header: /Net/u52/wkt/h22/2evl/RCS/Pdfmis.f,v 7.0 99/07/25 18:46:45 wkt Exp $
C #Log:	Pdfmis.f,v $
c Revision 7.0  99/07/25  18:46:45  wkt
c Pavel's improved version: no change
c 
c Revision 6.2  97/12/11  10:29:13  wkt
c Setpdf --> Setevl
c
c Revision 6.1  97/11/15  17:20:02  wkt
c v.6.0 + HLL revisions
c
C ====================================================================
C      FUNCTION PMSDIS (ISET, IHDRN, IPRTN, XX, QQ, IAT)

C     Front-end function to DIS2MS which transforms between DIS and MS-bar
C     scheme distribution functions:
C     Also set up the common block MSTDIS for use by the integrand functions

C     On input,    IACT = 1 :  performs the DIS --> MS-BAR conversion
C                        -1 :       "       DIS <-- MS-BAR     "
C                         0 :  no conversion
C                      other:  error condition

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / MSTDS1 / X, Q, NF, JSET, JHDRN, JPRTN, NCNT
C
      EXTERNAL C2QX, QQY, QGY, GQY
      SAVE
C
      DATA TINY, D1P, DUM, IST / 1D-20, 1.001, -999., -5 /
      DATA AERR, RERR, SML, D1M / 0.d0, 1.D-4, 1.D-6, 0.999 /

      IACT = IAT
      IF (IACT .EQ. 0) THEN
         PMSDIS = PDF(ISET, IHDRN, IPRTN, XX, QQ, IR)
         PRINT *, 'IACT=0 in PMSDIS call, no convertion is done.'
         RETURN
      ElseIF (ABS(IACT) .GT. 1) THEN
         PRINT *, 'Illegal value of Iact in PMSDIS: IACT =', IACT
         PMSDIS = DUM
         RETURN
      EndIf
      IF (ISET .NE. IST .OR. IHDRN .NE. JHDRN) THEN
        IST = ISET
        JSET  = ISET
        JHDRN = IHDRN
        IF(IST.LE.0 .OR. IST.GT.31 .OR. (IST.GT.11.AND.IST.LT.21)) THEN
          PRINT *, 'ISET = ', ISET, ' out of range in PMSDIS.'
          PMSDIS = DUM
          IR  = 11
          RETURN
        EndIf
      EndIf

      JPRTN = IPRTN
      Q     = QQ
      NF    = NFL(Q)
      X     = XX

      IF (X .LE. TINY .OR. X .GT. D1P) THEN
        PRINT '(A,1PE13.4,A)', ' X =', X, ' out of range in PMSDIS.'
        PMSDIS = DUM
        IR  = 13
        RETURN
      ElseIF (X .GT. D1M) THEN
        PMSDIS = 0.
        RETURN
      EndIf

      AL = ALPI (Q)
C                                                             Begin computation
      TEM0 = PDF (JSET, JHDRN, JPRTN, X, Q, IR0)
      IF (TEM0 .LT. -SML)
     >   PRINT '(A, 3(1PE13.4))', 'PDF < 0 in PMSDIS', TEM0, X, Q

      IF (TEM0 .LT. SML) THEN
         PMSDIS = TEM0
         RETURN
      EndIf

      AERR = ABS(TEM0 / 10.)
      IF (JPRTN .EQ. 0) THEN
         TMGQ1=-ADZINT (C2QX,D0,X,D0,RERR,ER2,IR, 2,1) * SQRK(X)
         TMGQ2=-ADZINT (GQY, X, D1, AERR, RERR, ER3, IR, 1,2)
         TMGG =-ADZINT (QGY, X,D1, AERR,RERR,ER1,IR, 1,2) * 2.*NF

         TEM = TEM0 + AL * (TMGQ2 - TMGQ1 + TMGG) * IACT
      Else
         TMQQ1= ADZINT (C2QX, D0,X,D0,RERR,ER1,IR, 2,1) * TEM0
         TMQQ2= ADZINT (QQY, X, D1, AERR, RERR, ER2, IR, 1,2)
         TMQG = ADZINT (QGY, X, D1, AERR, RERR, ER3, IR, 1,2)
         TEM  = TEM0 + AL * (TMQQ2 - TMQQ1 + TMQG) * IACT
      EndIf

      PMSDIS = TEM
C
      RETURN
C                        ****************************
      END

      FUNCTION PrDis1 (IPRTN, XX, QQ)
C                                                   -=-=- prdis1

C                                     Copy of PARDIS for the alternate PDF set
C                                  Only the names of common blocks are changed
c
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Character Msg*80
      LOGICAL LSTX
 
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      PARAMETER (Smll = 1D-9)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT

      COMMON / X1ARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / Q1RAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / Q1RAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / E1LPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / P1VLDT / UPD(MXPQX), KF, Nelmt
C
C     Strictly speaking, this mixes the the common blocks of the two versions
C     Usuall, the array sized are the same, this does not matter.
C     Will cause trouble only under very rare occasions (alternate block is
C     used, and the block sizes are different.
C     I am going to bother fixing this now.        ---  wkt  8/22/99

      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
 
      dimension fvec(4), fij(4)
      dimension xvpow(0:mxx)
      Data Iwrn1, Iwrn2, Iwrn3, OneP / 3*0, 1.00001 /

      data xpow / 0.3d0 /	!**** choice of interpolation variable
      data ientry / 0 /

c store the powers used for interpolation on first call...
      if(ientry .eq. 0) then
         ientry = 1

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1   
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PrDis1! x = ', x 
        Stop 
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
         Msg = '0 < X < Xmin in PrDis1; extrapolation used!'
         CALL WARNR (IWRN1, NWRT, Msg, 'X', X, Xmin, 1D0, 1)
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PrDis1! x = ', x 
        Stop 
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2 
      sy3 = ss - svec3 

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12 
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (Q .GE. QV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
         If (JLq .LT. 0) Then
          Msg = 'Q < Q0 in PrDis1; extrapolation used!'
          CALL WARNR (IWRN2, NWRT, Msg, 'Q', Q, Qini, 1D0, 1)
         EndIf
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3
        If (JLq .GE. Nt) Then
         Msg = 'Q > Qmax in PrDis1; extrapolation used!'
         CALL WARNR (IWRN3, NWRT, Msg, 'Q', Q, Qmax, 1D0, 1)
        Endif
      Endif
C                                   This is the interpolation variable in Q
      tt = log(log(q/al))

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf

      Data nqvec / 4 /

c get the pdf function values at the lattice points...

      jtmp = ((IPRTN + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1) 

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * Xa(1,2)
         fij(3) = Upd(J1+2) * Xa(2,2)
         fij(4) = Upd(J1+3) * Xa(3,2)
C 
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint (XVpow(0), Fij(1), 4, ss, Fx, Dfx) 
         
         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 
C                                              Pdf is undefinged for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:
C        Jon's algorithm does not handle this case, give it to Polint also.

        Call Polint (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)

        Fvec(it) = Fx

       Else 
C                       for all interior points, use Jon's in-line function 
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)          
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1) 
     &               + const6*(Upd(J1+3)-g4) 
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...
C              Again, use Polint for end intervals, Jon for the interrior

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint (TV(0), Fvec(1), 4, tt, ff, Dfq)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 
     &	  +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PrDis1 = ff

      Return
C					********************
      End

      SUBROUTINE QARRAY (NINI)
C                                                   -=-=- qarray

C            (These general comments are enclosed in the
C            lead subprogram in order to survive forsplit.)

C ====================================================================
C GroupName: Evlutl
C Description: routines to manipulate grids and integration for evolution
C ListOfFiles: qarray xarray integr hinteg xfrmz zfrmx xfz zfxl dxdz
C ====================================================================
C #Header: /Net/u52/wkt/h22/2evl/RCS/Evlutl.f,v 7.1 99/08/19 17:45:24 wkt Exp $
C #Log:	Evlutl.f,v $
c Revision 7.1  99/08/19  17:45:24  wkt
c Slight modification in Xarray to return the exact Xv(1)=Xmin.
c 
c Revision 7.0  99/07/25  18:40:33  wkt
c Pavel's improved version: no change
c 
c Revision 6.4  97/12/11  11:14:19  wkt
c Some changes for consistency
c
c Revision 6.3  97/12/11  10:28:36  wkt
c cleanup unused local variables
c
c Revision 6.2  97/12/03  23:06:37  wkt
c xfrmz restored to the original log-linear version from the x**delta version
c Liang claims the latter has problems
c
c Revision 6.1  97/11/15  17:19:28  wkt
c v.6.0 + HLL revisions
c
C ====================================================================
C      SUBROUTINE QARRAY (NINI)
C
C       Given Qini, Qmax, and Nt, this routine go through the various
C       flavor thresholds; determines the step-sizes in-between each pair
C       of thresholds, computes a new NT, if necessary; and returns:
C        the variable NINI, NFMX as arguments;
C        the (NT+1)-dim Q- and T=lnlnQ/Lamda arrays in the common block QARAY1;
C        the Neff-dependent variable (see below) in QARAY2; and
C        KF in the common block EvlPac
C
C       Note in particular that the value of NT may increase by 1 or 2 in a
C       call to this routine because of the logistics of setting up the steps
C       for the flavor thresholds.
C
C      For given Neff, the variables TLN, DTN, NTL and NTN has the following
C      meaning:

C     Flavor thresholds
C     Nini        Nf                                          Nf+1  ....NFmx
C     |.. ........|<- DTN ->|<- DTN ->| .......................|..........|
C                  <------------  NTL(Nf) steps  ------------->
C     0         NTN(Nf)                                    NTN(Nf+1)     NT
C               TLN(Nf)       (TLN=lnln(Q/Lamda))          TNL(Nf+1)
C    TV(0).... TV(NTN(Nf)).....TV(NTN(Nf)+2).............TV(NTN(Nf+1))..TV(NT)
C    QV(0).... QV(NTN(Nf)).....QV(NTN(Nf)+2).............QV(NTN(Nf+1))..QV(NT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / QARAY2 / TLN(MXF), DTN(MXF), NTL(MXF), NTN(MXF)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      NCNT = 0

      IF (NT .GE. mxq) NT = mxq - 1

C                                        Set up overall t=log(Q) parameters

      S = LOG(QINI/AL)
      TINI = LOG(S)
      S = LOG(QMAX/AL)
      TMAX = LOG(S)
C                              Nt is the provisional # of mesh-points in t;
C                                    Dt0 is the approximate increment in T,
C                           Actual values of (Nt, Dt) are determined later.
    1 DT0 = (TMAX - TINI) / float(NT)
C                            Determine effective Nflv at Qini and at Qmax
      NINI = NFL(QINI)
      NFMX = NFL(QMAX)
      Call ParQcd (2, 'ORDER', Ord, Ir)
      Call ParQcd (2, 'ALAM', Al0, Ir)
      Call ParQcd (2, 'NFL', Afl0, Ir)
C                                      Set the total number of Quark flavors
C                                      in the QCD package to NfMx
      AFL = NfMx
      Call ParQcd (1, 'NFL', AFL, Ir)
C                                      and restore the QCD coupling
      Iordr = Nint (Ord)
      Ifl0  = Nint (Afl0)
      Call SetLam (Ifl0, Al0, Iordr)
C
C      Q2 evolution is carried out in stages separated by flavor thresholds
C
      NG = NFMX - NINI + 1
C                                                  -----------------------
C                                      Determine the threshold points and
C                                     Set up detailed  t- mesh structure
      QIN  = QINI
      QOUT = QINI
      S = LOG(QIN/AL)
      TIN  = LOG(S)
      TLN(1) = TIN
      NTL(1)  = 0
      QV(0) = QINI
      TV(0) = Tin
C
      DO 20 NEFF = NINI, NFMX
C
        ICNT = NEFF - NINI + 1
        IF (NEFF .LT. NFMX) THEN
          THRN = AMHATF (NEFF + 1)
          QOUN = MIN (QMAX, THRN)
        Else
          QOUN = QMAX
        EndIf
C
        IF (QOUN-QOUT .LE. 0.0001) THEN
          DT   = 0
          NITR = 0
        Else
          QOUT = QOUN
          S = LOG(QOUT/AL)
          TOUT = LOG(S)
          TEM = TOUT - TIN
C                               Nitr = Number of iterations in this stage
          NITR = INT (TEM / DT0) + 1
          DT  = TEM / NITR
        EndIf
C                                       Save book-keeping data on array
        DTN (ICNT) = DT
        NTN (ICNT) = NITR
        TLN (ICNT) = TIN
        NTL (ICNT+1) = NTL(ICNT) + NITR
C
C                      QV is the physical Q-value for this lattice point.
        IF (NITR .NE. 0) THEN
        DO 205 I = 1, NITR
           TV (NTL(ICNT)+I) = TIN + DT * I
           S = EXP (TV(NTL(ICNT)+I))
           QV (NTL(ICNT)+I) = AL * EXP (S)
  205   CONTINUE
        EndIf
C                                        Initialize for next iteration
        QIN = QOUT
        TIN = TOUT
C
   20 CONTINUE
C            Redefine Nt, as the actual number of t-points; may have been
C                             increased by 1 or 2 by the above algorithm
C                       If NT > MXQ, reduce NT until it is within bounds
      NCNT = NCNT + 1
      NTP = NTL (NG + 1)
      ND  = NTP - NT

      IF (NTP .GE. MXQ) THEN
         NT = MXQ - ND - NCNT
         GOTO 1
      EndIf
      NT = NTP
C
      RETURN
C                        ****************************
      END


      function  QFF1_pn (XX)
      implicit none
      double precision xx,qff1_pn,pff1_pn
      qff1_pn=pff1_pn(xx)*(1.-xx)
      RETURN
      End !qff1_pn

      function  QFF2_pn (XX)
      implicit none
      double precision xx,qff2_pn,pff2_pn
      qff2_pn=pff2_pn(xx)*(1.-xx)
      RETURN
      End!qff2_pn

      function  QFG1_pn (XX)
      implicit none
      double precision xx,qfg1_pn,pfg1_pn
      qfg1_pn=pfg1_pn(xx)*(1.-xx)
      RETURN
      End!qfg1_pn

      function  QFG2_pn (XX)
      implicit none
      double precision xx,qfg2_pn,pfg2_pn
      qfg2_pn=pfg2_pn(xx)*(1.-xx)
      RETURN
      End!qfg2_pn

      function  QGF1_pn (XX)
      implicit none
      double precision xx,qgf1_pn,pgf1_pn
      qgf1_pn=pgf1_pn(xx)*(1.-xx)
      RETURN
      End!qgf1_pn
   
      function  QGF2_pn (XX)
      implicit none
      double precision xx,qgf2_pn,pgf2_pn
      qgf2_pn=pgf2_pn(xx)*(1.-xx)
      RETURN
      End!qgf2_pn
   
      function QGG1_pn (XX)
      implicit none
      double precision xx,qgg1_pn,pgg1_pn
      qgg1_pn=pgg1_pn(xx)*(1.-xx)
      RETURN
      End!qgg1_pn

      function QGG2_pn (XX)
      implicit none
      double precision xx,qgg2_pn,pgg2_pn
      qgg2_pn=pgg2_pn(xx)*(1.-xx)
      RETURN
      End!qgg2_pn

      function  RFF1_pn (XX)
      implicit none
      double precision xx,rff1_pn,pff1_pn
      rff1_pn=pff1_pn(xx)*xx
      RETURN
      End!rff1_pn

cpn-------------------------------------------------------------------------


      function  RFF2_pn (XX)
      implicit none
      double precision xx,rff2_pn,pff2_pn
      rff2_pn=pff2_pn(xx)*xx
      RETURN
      End!rff2_pn

cpn-------------------------------------------------------------------------


      function  RFG1_pn (XX)
      implicit none
      double precision xx,rfg1_pn,pfg1_pn
      rfg1_pn=pfg1_pn(xx)*xx
      RETURN
      End!rfg1_pn

cpn-------------------------------------------------------------------------
 

      function  RFG2_pn (XX)
      implicit none
      double precision xx,rfg2_pn,pfg2_pn
      rfg2_pn=pfg2_pn(xx)*xx
      RETURN
      End!rfg2_pn

cpn-------------------------------------------------------------------------
 

      function  RGF1_pn (XX)
      implicit none
      double precision xx,rgf1_pn,pgf1_pn
      rgf1_pn=pgf1_pn(xx)*xx
      RETURN
      End!rgf1_pn

cpn-------------------------------------------------------------------------

      function  RGF2_pn (XX)
      implicit none
      double precision xx,rgf2_pn,pgf2_pn
      rgf2_pn=pgf2_pn(xx)*xx
      RETURN
      End!rgf2_pn

cpn-------------------------------------------------------------------------

      function RGG1_pn (XX)
      implicit none
      double precision xx,rgg1_pn,pgg1_pn
      rgg1_pn=pgg1_pn(xx)*xx
      RETURN
      End!rgg1_pn

cpn-------------------------------------------------------------------------

      function RGG2_pn (XX)
      implicit none
      double precision xx,rgg2_pn,pgg2_pn
      rgg2_pn=pgg2_pn(xx)*xx
      RETURN
      End!rgg2_pn

cpn-------------------------------------------------------------------------

      Subroutine SetEvl
C ====================================================================
C GroupName: Setevl
C Description: This is the lead Group with setup functions; DatEvl in setevl
C ListOfFiles: setevl pdfset
C ====================================================================
C $Header: /Net/u52/wkt/h22/2evl/RCS/Setevl.f,v 7.2 99/08/20 11:48:15 wkt Exp $
C $Log:	Setevl.f,v $
c Revision 7.2  99/08/20  11:48:15  wkt
c *** empty log message ***
c 
c Revision 7.1  99/08/12  09:43:33  wkt
c Default values of Nx, Xmin, .. in DatEvl changed to more commonly
c used settings.  No change in content of program.
c 
c Revision 7.0  99/07/25  19:07:35  wkt
c Pavel's improved version: No change in this module
c 
c Revision 6.8  98/08/25  00:20:55  wkt
c (i) data Xcr --> 1.5 to go along with the current zfrmx;
c (ii) parameter notation clarified in more comments in DatEvl.
c
c Revision 6.7  97/12/11  11:15:00  wkt
c Some changes for consistency
c
c Revision 6.6  97/12/11  10:29:37  wkt
c
c Revision 6.5  97/11/23  00:48:02  wkt
c Header section removed to become Evlpac.h
c
c Revision 6.3  97/11/16  10:58:54  wkt
c Added evlini, iniread, fitini to the Runevl group
c
c Revision 6.2  97/11/16  00:31:48  wkt
c one common / PdfMod / and its data statement added in block data module:
c
c Revision 6.1  97/11/15  17:20:15  wkt
c v.6.0 + HLL revisions
c
C ====================================================================
c Revision 5.96  96/10/28  00:21:13  wkt
c logistics
c
c Revision 5.93  96/10/27  15:25:46  wkt
c Parameter MxAdF (meant for v.6) value brought in line with the rest of pac.
c
c Revision 5.9  96/10/20  13:42:08  wkt
c 1. Consolidated with Liang's Cteq4 version;
c 2. Cleaned up numerous nuisances (variables not used, ..) which cause
c    compilation warnings
c 3. ZfrmX and XfrmZ and ancillaries replaced with new version based on
c    fractional power rather than log + linear transformation;
c 4. New Cteq, Mrs, Grv switches put in; but need their programs and
c    tables to run.
c
c Revision 5.6  96/06/02  23:36:31  wkt
c minor
c
c Revision 5.5  96/06/02  23:21:16  wkt
c Write and Read to file converted to 'Formatted'
c Commons EVLPAC and PEVLDT modified.
c
c Revision 5.1  95/07/20  15:40:23  wkt
c Evlini and InRead removed;
c ((They belong to FitPac; cause linking problem
c x2zconv retained.
c 14
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C                     Force link to DatEvl + other setup work if necessary
      External DatEvl

      Dummy = 0.

      Return
C                        ****************************
      END

      BLOCK DATA DatEvl
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL LSTX
C
      PARAMETER (Z = 1D -10, ZZ = 1D -20)

C              MxF is the maximum number of quark flavors
C              MxAdF = 1 gluon + 1 singlet quark =   2  (This version)
C                    ( + Any additional fake flavors
C         to save space for storing info on evolution for the next version)

C                  |                                      |
C                  |                                      |   _
C             SnQr |  Total singlet quark       (NfSn)    |   ^
C               t  |           .                ( MxF)    |
C               b  |           .                          |
C                  |           d                          |   M
C                  |           u         MxQ Q values     |   X
C                  |---------  G ------------->           |   P
C                / |          ubar          Q             |   N
C               /  |          dbar
C           x  /   |           .
C             /    |           .                (-MxF)    |   _
C     MxX x values

      PARAMETER (MxX = 204, MxQ = 25, MxF = 6, MxAdF = 2)
      PARAMETER (MxPN = MxF * 2 + MxAdF)
      PARAMETER (MxQX= MxQ * MxX,   MxPQX = MxQX * MxPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C                                These common blocks belong to evlpac proper
      Common / PdCntrl/ LPrt, LDbg
      COMMON / XXARAY / XCR, XMIN, XV(0:MxX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MxQ),TV(0:MxQ), NT,JT,NG
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX), KF, Nelmt

      Data LPrt, LDbg / 1, 0 /
      DATA QINI, QMAX, XMIN, XCR / 1.0, 1D4, .999999D-4, 1.5 /
      DATA KF, IKNL, IPD0, IHDN / 10, 2, 1, 1 /
      DATA NX, NT, JT, LSTX / 51, 10, 1, .FALSE. /
C
C                        ****************************
      END

      SUBROUTINE SNEVL(IKNL, NX, NT, JT, DT, TIN, NEFF, UI, GI, US, GS)
C
C       This is the singlet counter-part of the NSEVL subroutine. Refer to
C                       comments at the beginning of that program section.
C
C     Input parton distributions are Gi (for gluon) and Ui (for singlet quark)

C                               at Tt = 0; output distributions are Gs and Us
C                                       at Tt = IS*dt with IS = 1, 2, ... , Nt.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXQX= MXQ * MXX)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      Common / PdCntrl/ LPrt, LDbg
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
C
      DIMENSION UI(NX), US(0:NX, 0:NT)
      DIMENSION GI(NX), GS(0:NX, 0:NT)
      DIMENSION Y0(MXX), Y1(MXX), YP(MXX), F0(MXX), F1(MXX), FP(MXX)
      DIMENSION Z0(MXX), Z1(MXX), ZP(MXX), G0(MXX), G1(MXX), GP(MXX)

      DATA D0 / 0.0 /

         If (Ldbg .Eq. 1) Then
            Write (Nwrt, '(A)') 'Singlet:'
            N5 = Nx / 5 + 1
         Endif
C                                       Faster evolution of the singlet sector
C                                       requires finer iteration to achieve the
C                                       same accuracy  as in the non-singlet.
      JTT = 2 * JT
      DDT = DT / JTT
C
      IF (NX .GT. MXX) THEN
      WRITE (NOUT,*) 'Nx =', NX, ' greater than Max # of pts in SNEVL.'
      STOP 'Program stopped in SNEVL'
      EndIf
C                                       Compute an effective first order lamda
C                                       to be used in checking of moment evl.
      TMD = TIN + DT * NT / 2.
      AMU = EXP(TMD)
      TEM = 6./ (33.- 2.* NEFF) / ALPI(AMU)
      TLAM = TMD - TEM
C                                      Initialization (see previous subroutine)
      DO 9 IX = 1, NX
      US (IX, 0) = UI(IX)
      GS (IX, 0) = GI(IX)
    9 CONTINUE
      US ( 0, 0) = (UI(1) - UI(2))* 3D0 + UI(3)
      GS ( 0, 0) = (GI(1) - GI(2))* 3D0 + GI(3)
C
      TT = TIN
      DO 10 IZ = 1, NX
      Y0(IZ) = UI(IZ)
      Z0(IZ) = GI(IZ)
   10 CONTINUE
C                                                       loop in the Tt variable
      DO 20 IS = 1, NT
C                                                       fine- grained iteration
         DO 202 JS = 1, JTT
C                                          Irnd is the counter for Q-iterations
            IRND = (IS-1) * JTT + JS
C                                               Use Runge-Katta the first round
            IF (IRND .EQ. 1) THEN
C
                CALL SNRHS (TT, NEFF, Y0,Z0,  F0,G0)
C
                DO 250 IZ = 1, NX
                   Y0(IZ) = Y0(IZ) + DDT * F0(IZ)
                   Z0(IZ) = Z0(IZ) + DDT * G0(IZ)
  250           CONTINUE
C
                TT = TT + DDT
C
                CALL SNRHS (TT, NEFF, Y0, Z0,  F1, G1)
C
                DO 251 IZ = 1, NX
                   Y1(IZ) = UI(IZ) + DDT * (F0(IZ) + F1(IZ)) / 2D0
                   Z1(IZ) = GI(IZ) + DDT * (G0(IZ) + G1(IZ)) / 2D0
  251           CONTINUE
C                            What follows is a combination of the 2-step method
C                                   and the Adams Predictor-Corrector Algorithm
            Else
C
                CALL SNRHS (TT, NEFF, Y1, Z1,  F1, G1)
C                                                                     Predictor
                DO 252 IZ = 1, NX
                   YP(IZ) = Y1(IZ) + DDT * (3D0 * F1(IZ) - F0(IZ)) / 2D0
                   ZP(IZ) = Z1(IZ) + DDT * (3D0 * G1(IZ) - G0(IZ)) / 2D0
  252           CONTINUE
C                        Increment of Tt at this place is part of the formalism
                TT = TT + DDT
C
                CALL SNRHS (TT, NEFF, YP, ZP,  FP, GP)
C                                                                     Corrector
                DO 253 IZ = 1, NX
                   Y1(IZ) = Y1(IZ) + DDT * (FP(IZ) + F1(IZ)) / 2D0
                   Z1(IZ) = Z1(IZ) + DDT * (GP(IZ) + G1(IZ)) / 2D0
                   F0(IZ) = F1(IZ)
                   G0(IZ) = G1(IZ)
  253           CONTINUE
            EndIf
C
  202    CONTINUE
C                       Fill output array and restore factor of X, if necessary
C                                    For spin-averaged case, enforce positivity
cpn                    As of Dec. 11, 2001, the positivity is no more enforced
         DO 260 IX = 1, NX
cpn 2001           IF (IKNL .GT. 0) THEN
c           if (iknl .eq. 1 .or. iknl .eq. 2) then
c            US (IX, IS) = MAX(Y1(IX), D0)
c            GS (IX, IS) = MAX(Z1(IX), D0)
c           Else
            US (IX, IS) = Y1(IX)
            GS (IX, IS) = Z1(IX)
c           EndIf
  260    CONTINUE

C
C               The value of the function at x=0 is obtained by extrapolation
C
         US(0, IS) = 3D0*Y1(1) - 3D0*Y1(2) + Y1(3)
         GS(0, IS) = 3D0*Z1(1) - 3D0*Z1(2) + Z1(3)
C
C                                                    Print out for Debugging
      If (LDbg .Eq. 1) Then
         Write (Nwrt, '(A, 5(1pE12.3))') ' SQ:',(Us(Iz,Is), Iz=1,Nx,N5)
         Write (Nwrt, '(A, 5(1pE12.3))') '  G:',(Gs(Iz,Is), Iz=1,Nx,N5)
      Endif

   20 CONTINUE
C
      RETURN
C                        ****************************
      END


      SUBROUTINE SNRHS (TT, NEFF, FI, GI,  FO, GO)
C
C       Subroutine to compute the Right-Side of the Altarelli-Parisi Equation
C                                                       for the Singlet sector:
C       See comments in NSRHSP for notes on IKNL
C
C       FI, Z are the input distributions for quark and gluon respectively;
C       FO, G are the output dY/dt, dZ/dt.
C       Nx is the number of mesh-points, Tt is the Log Q variable.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
      PARAMETER (PI = 3.141592653589793, PI2 = PI ** 2)

      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / KRNL01 / AFF2 (MXX),AFG2 (MXX), AGF2 (MXX), AGG2 (MXX),
     >                  ANSP (MXX), ANSM (MXX), ZGG2, ZFF2, ZQQB
      COMMON / KRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      DIMENSION GI(NX), GO(NX), G1(MXX), G2(MXX), G3(MXX), G4(MXX)
      DIMENSION FI(NX), FO(NX), W0(MXX), W1(MXX), WH(MXX), WM(MXX)
      DIMENSION R0(MXX), R1(MXX), R2(MXX), RH(MXX), RM(MXX)

      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg

cpn 2001 <- Transformation coefficients from the MS-bar to DIS scheme
cpn         (used if iknl==12)
      COMMON / KRNDIS / DFFG(MXX, MXX), DGGF(MXX, MXX), 
     >                  ADFF2(MXX), ADFG2(MXX), ADGF2(MXX), ADGG2(MXX),
     >                  ZDFF2, ZDGF2
      DIMENSION D1(MXX), D2(MXX), D3(MXX), D4(MXX)
cpn 2001 ->
  

C
      S = EXP(TT)
      Q = AL * EXP (S)
      CPL = ALPI(Q)
      CPL2= CPL ** 2 / 2. * S
      CPL = CPL * S
C
      CALL INTEGR (NX,-1, FI, WM, IR1)
      CALL INTEGR (NX, 0, FI, W0, IR2)
      CALL INTEGR (NX, 1, FI, W1, IR3)
      CALL INTEGR (NX,-1, GI, RM, IR4)
      CALL INTEGR (NX, 0, GI, R0, IR5)
      CALL INTEGR (NX, 1, GI, R1, IR6)
      CALL INTEGR (NX, 2, GI, R2, IR7)

      CALL HINTEG (NX,    FI, WH)
      CALL HINTEG (NX,    GI, RH)


      IF (IKNL .GT. 0) THEN     !LO unpolarized case      
        DO 230 IZ = 1, NX
          FO(IZ) = ( 3d0/2D0*Cf * FI(IZ)
     >         + Cf * ( 2D0 * WH(IZ) - W0(IZ) - W1(IZ) ))
     >         + 2.*Tr*NEFF * ( R0(IZ) - 2D0 * R1(IZ) + 2D0 * R2(IZ) )
          FO(IZ) = FO(IZ) * CPL
C
          GO(IZ) = Cf * ( 2D0 * WM(IZ) - 2D0 * W0(IZ)  + W1(IZ) )
     >         + (11D0*Cg/3d0 - 4d0/3d0 * Tr*Neff) / 2D0 * GI(IZ) + 2*Cg
     >         * (RH(IZ) + RM(IZ) - 2D0 * R0(IZ) + R1(IZ) - R2(IZ))
          GO(IZ) = GO(IZ) * CPL
 230    CONTINUE
      Else                      !LO polarized case  
        DO 240 IZ = 1, NX
          FO(IZ) = 2.*Tr*NEFF * (-R0(IZ) + 2.* R1(IZ) )
     >         + 1.5d0*Cf* FI(IZ) + Cf* ( 2.* WH(IZ) - W0(IZ) - W1(IZ) )
          FO(IZ) = FO(IZ) * CPL
C
          GO(IZ) = Cf* ( 2.* W0(IZ) - W1(IZ) )
     >         + (11D0*Cg/3d0 - 4d0/3d0 * Tr*Neff) / 2D0* GI(IZ) +
     >         2.*Cg*(RH(IZ) + R0(IZ) - 2.* R1(IZ))
          GO(IZ) = GO(IZ) * CPL
 240    CONTINUE
      EndIf

C
cpn 2001      IF (abs(IKNL) .EQ. 2) THEN
      IF (abs(IKNL) .GE. 2) THEN
        DZ = 1./(NX - 1)
        DO 21 I = 1, NX-1
          X = XV(I)
C                                        Number of points in the I-th integral
          NP = NX - I + 1
          IS = NP
C                                 Evaluate the integrand for the I-th integral
cpn                               Do it separately in the unpolarized and
c                                 polarized case
cpn 2001
          if (iknl.ge.2) then
cpn
            g2(1)=0d0
            g3(1)=0d0
          
            DO KZ = 2, NP
              IY = I + KZ - 1
              IT = NX - IY + 1
C                            Quantities associated with kernel functions
              XY = ZZ (IS, IT)
C                                                       2nd order terms
              G1(KZ) = FFG(I, IY) * (FI(IY) - XY**2 *FI(I))

              G4(KZ) = GGF(I, IY) * (GI(IY) - XY**2 *GI(I))

              G2(KZ) = FFG(IS,IT) * (GI(IY) - xy*GI(I)) !FG
              G3(KZ) = GGF(IS,IT) * (FI(IY) - XY*FI(I)) !GF (usual notations)

cpn 2001 <- Transformation terms from the MS-bar to DIS scheme
           if (iknl.eq.12) then
             
             D1(KZ) =  DFFG(I, IY) * (FI(IY) -XY**2*FI(I))
             D2(KZ) =  DFFG(IS,IT) * (GI(IY) - xy*GI(I))      !FG 
             D3(KZ) =  DGGF(I, IY) * (FI(IY) - XY**2 *FI(I))  !GF (usual 
                                                             ! notations) 
             D4(KZ) =  DGGF(IS,IT) * (GI(IY) - XY*GI(I))    

cpn Alternative calculation of the integrands: slow, but more accurate
c             D1(kz) = DelFF2_pn(xv(i)/xv(iy))*dxtz(iy) * 
c     >         (FI(IY) - XY**2 *FI(I))
c             D2(KZ) =  DelFG2_pn(xv(i)/xv(iy))*dxtz(iy)* 
c     >         (GI(IY) - xy*GI(I))                         !FG 
c             D3(KZ) = DelGF2_pn(xv(i)/xv(iy))*dxtz(iy) *
c     >         (FI(IY) - XY**2 *FI(I))                     !GF (usual 
c                                                           ! notations) 
c             D4(KZ) = DelGG2_pn(xv(i)/xv(iy))*dxtz(iy) * 
c     >         (GI(IY) - XY*GI(I))    
  
           endif                !iknl.eq.12
cpn 2001 ->

            enddo               !kz
C
            TEM1 = SMPNOL (NP, DZ, G1, ERR)
            TEM2 = SMPSNA (NP, DZ, G2, ERR)
            TEM3 = SMPSNA (NP, DZ, G3, ERR)
            TEM4 = SMPNOL (NP, DZ, G4, ERR)

C        PRINT '(A, F12.3/ A, 4(1PE15.3))', ' Q =', Q, ' Int',
C     >          TEM1, TEM4, TEM2, TEM3

            TEM1 = TEM1 - FI(I) * (AFF2(I) + ZFF2)
            TEM4 = TEM4 - GI(I) * (AGG2(I) + ZGG2)
cpn
            tem2 = tem2 + GI(I)*AFG2(I)
            tem3=  tem3 + FI(I)*AGF2(I)

cpn 2001 <- In DIS scheme, add transformation terms
            if (iknl.eq.12) then
                                ! Integrals from x to 1
c$$$              tem1 = tem1 + SMPNOL (NP, DZ, D1, ERR)
c$$$              tem2 = tem2 + SMPSNA (NP, DZ, D2, ERR)
c$$$              tem3 = tem3 + SMPNOL (NP, DZ, D3, ERR)
c$$$              tem4 = tem4 + SMPSNA (NP, DZ, D4, ERR)
c$$$                                ! End-point terms
c$$$              TEM1 = TEM1 - FI(I) *(ADFF2(I) - ZDFF2)
c$$$              tem2 = tem2 + GI(I) *ADFG2(I)          
c$$$              TEM3 = TEM3 - FI(I) *(ADGF2(I) - ZDGF2)
c$$$              tem4 = tem4 + GI(I) *ADGG2(I)

cpn 2001 Same code written in a longer way for debugging purposes
              tem11 = SMPNOL (NP, DZ, D1, ERR)
              tem22 = SMPSNA (NP, DZ, D2, ERR)
              tem33 = SMPNOL (NP, DZ, D3, ERR)
              tem44 = SMPSNA (NP, DZ, D4, ERR)
                                ! End-point terms
              TEM11 = TEM11 - FI(I) *(ADFF2(I) - ZDFF2)
              tem22 = tem22 + GI(I) *ADFG2(I)          
              TEM33 = TEM33 - FI(I) *(ADGF2(I) - ZDGF2)
              tem44 = tem44 + GI(I) *ADGG2(I)

              tem1 = tem1 + tem11
              tem2 = tem2 + tem22
              tem3 = tem3 + tem33
              tem4 = tem4 + tem44

            endif !iknl.eq.12
cpn 2001 ->

          elseif(iknl.eq.-2) then
            g1(1)=0d0
            g2(1)=0d0
            g3(1)=0d0
            g4(1)=0d0

            DO KZ = 2, NP
              IY = I + KZ - 1
              IT = NX - IY + 1
C                            Quantities associated with kernel functions
              XY = ZZ (IS, IT)
C                                                       2nd order terms
              G1(KZ) = FFG(I, IY) * (FI(IY) - XY*FI(I)) !FF
              G2(KZ) = FFG(IS,IT) * (GI(IY) - xy*GI(I)) !FG
              G3(KZ) = GGF(IS,IT) * (FI(IY) - XY*FI(I)) !GF (usual notations)
              G4(KZ) = GGF(I, IY) * (GI(IY) - XY*GI(I)) !GG
            ENDDO               !kz
C
            TEM1 = SMPSNA (NP, DZ, G1, ERR)
            TEM2 = SMPSNA (NP, DZ, G2, ERR)
            TEM3 = SMPSNA (NP, DZ, G3, ERR)
            TEM4 = SMPSNA (NP, DZ, G4, ERR)
  
cpn                      Add 1/(1-x)_{+} terms to the FF and GG convolutions
            tem1=tem1+((67./9.- PI**2 / 3.)*CfCg
     >           -20d0/9d0*CfTr*Neff)*WH(I)
            tem4=tem4 +(- 20./9.*CgTr*Neff+(67./9.-Pi**2/3.)*Cg2)*RH(I)
  
cpn                                                    Add delta(1-x)  terms 
            TEM1 = TEM1 + FI(I) * (AFF2(I) + zff2)
            tem2 = tem2 + GI(I) *  AFG2(I)
            tem3=  tem3 + FI(I) *  AGF2(I)
            TEM4 = TEM4 + GI(I) * (AGG2(I) + zgg2)            

          endif                 !iknl

          TMF = TEM1 + TEM2
          TMG = TEM3 + TEM4

          FO(I) = FO(I) + TMF * CPL2
          GO(I) = GO(I) + TMG * CPL2
C
 21     CONTINUE
      EndIf

      RETURN
C                         *************************
      END

      FUNCTION SPENCE (X)
C                                                   -=-=- spence
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      EXTERNAL SPNINT, SPN2IN
C
      COMMON / SPENCC / XX
C
      DATA U1, AERR, RERR / 1.D0, 1.E-7, 5.E-3 /
C
      XX = X
      TEM = GausInt(SPNINT, XX, U1, AERR, RERR, ERR, IRT)
      SPENCE = TEM
C
      RETURN
C                        ----------------------------
      ENTRY SPENC2 (X)
C
      XX = X
      TEM = GausInt(SPN2IN, XX, U1, AERR, RERR, ERR, IRT)
      SPENC2 = TEM + LOG (XX) ** 2 / 2.
C
      RETURN
C                        ****************************
      END

      FUNCTION SPNINT (ZZ)
C                                                   -=-=- spnint
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / SPENCC / X
C
      Z = ZZ
      TEM = LOG (1.- Z) / Z
      SPNINT = TEM
C
      RETURN
C                        ----------------------------
      ENTRY SPN2IN (ZZ)
C
      Z = ZZ
      TEM = LOG (1.+ X - Z) / Z
      SPN2IN = TEM
C
      RETURN
C                        ****************************
      END

      FUNCTION SQRK (Y)
C                                                   -=-=- sqrk
C                                                    Returns singlet quark PDF
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / MSTDS1 / X, Q, NF, JSET, JHDRN, JPRTN, NCNT

      TEM = 0
      DO 201 IFL = 1, NF
         TEM = TEM + PDF(JSET, JHDRN, IFL, Y, Q, IR)
     >             + PDF(JSET, JHDRN,-IFL, Y, Q, IR)
  201 CONTINUE
      SQRK = TEM
      RETURN
C                        ****************************
      END
C
       subroutine StpQCD
C  Sets up various group factors
       implicit NONE
       common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
       double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg

       Cf=4d0/3d0
       Cg=3d0
       Tr=0.5
       CfTr=Cf*Tr
       CgTr=Cg*Tr
       Cg2=Cg**2
       CfCg=Cf*Cg
       Cf2=Cf*Cf

       End! StpQCD
cpn---------------------------------------------------------------------------

      SUBROUTINE STUPKL (NFL)
C                                 Set up the common block containing the arrays
C                                        for the first and second order kernels
C
C                                 Also calculates the integrals and constants
C                                   needed to complete the [p(x)]sub+ integrals

C                         Real calculation in done in the routine KERNEL which
C                           follows strictly the Furmanski-Petronzio notation.
C                     This routine converts their convention to my convention.
C                                              FG (mine) = GF (their)
C                                              GF (mine) = FG (their)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MX = 3)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / KRN1ST / FF1(0:MXX), FG1(0:MXX), GF1(0:MXX), GG1(0:MXX),
     >     FF2(0:MXX), FG2(0:MXX), GF2(0:MXX), GG2(0:MXX),
     >     PNSP(0:MXX), PNSM(0:MXX)
      COMMON / KRN2ND / FFG(MXX, MXX), GGF(MXX, MXX), PNS(MXX, MXX)
      COMMON / KRNL00 / DZ, XL(MX), NNX
      COMMON / KRNL01 / AFF2 (MXX),AFG2 (MXX), AGF2 (MXX), AGG2 (MXX),
     >     ANSP (MXX), ANSM (MXX), ZGG2, ZFF2, ZQQB
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      
      EXTERNAL PFF1_pn, RGG1_pn, RFF2_pn, RGG2_Pn
      external pfg1_pn,pfg2_pn,pgf1_pn,pgf2_pn
      EXTERNAL FNSP_pn, FNSM_pn
      external dPFF2_pn,dPGF2_pn,dPFG2_pn,dPGG2_pn  
cpn3 
c      external dFNSM_pn
cpn
      dimension aff1(mxx),agg1(mxx)
      PARAMETER (PI = 3.141592653589793d0, PI2 = PI**2)
C
      common/qcdpar2/Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      double precision Cf,Cg,Tr,Cf2,Cg2,CfTr,CgTr,CfCg
      common/nflav/nflv 

cpn 2001 <- Transformation coefficients from the MS-bar to DIS scheme
cpn         (used if iknl==12)
      COMMON / KRNDIS / DFFG(MXX, MXX), DGGF(MXX, MXX), 
     >                  ADFF2(MXX), ADFG2(MXX), ADGF2(MXX), ADGG2(MXX),
     >                  ZDFF2, ZDGF2
      
      dimension  DFF2(0:MXX), DFG2(0:MXX), DGF2(0:MXX), DGG2(0:MXX)
      external RDFF2_pn, RDGF2_pn, DelFG2_pn, DelGG2_pn
cpn 2001 ->

cpn2 
      double precision wdff2_pn, wdfg2_pn, wdgf2_pn, wdgg2_pn
      external wdff2_pn, wdfg2_pn, wdgf2_pn, wdgg2_pn

      data iwrn /0/

      SAVE
      DATA AERR, RERR / 0.0d0, 1.d-4 /

      nflv=nfl 
C                         I = 1, NX  corresponds to  Z = 0, 1  or  X = Xmin, 1
C                           The point I = 0 corresponding to  Z = -DZ.
C                           The point X = 0. is singular, it cannot be reached.
C
C                          First past these quantities to FUNCTION PFF1 ... etc
C                   via / KRNL00 / for interpolation and extrapolation purposes
      NNX = NX
      DZ = 1./ (NX - 1)
      DO 5 I0 = 1, MX
        XL(I0) = XV(I0)
    5 CONTINUE

C                                          Calculation & switching FG <--> GF
      DO 10 I = 1, NX-1
        XZ = XV(I)
        CALL KERNEL (XZ, FF1(I), GF1(I), FG1(I), GG1(I), PNSP(I),
     >       PNSM(I), FF2(I), GF2(I), FG2(I), GG2(I), NFL, IRT)

cpn 2001 <- For the DIS scheme
        if (iknl.eq.12) 
     >     CALL DELDIS(XZ, DFF2(I), DFG2(I), DGF2(I), DGG2(I), NFL, IRT)
cpn 2001 ->


 10   CONTINUE
C
      FF1(0) = FF1(1) * 3. - FF1(2) * 3. + FF1(3)
      FG1(0) = FG1(1) * 3. - FG1(2) * 3. + FG1(3)
      GF1(0) = GF1(1) * 3. - GF1(2) * 3. + GF1(3)
      GG1(0) = GG1(1) * 3. - GG1(2) * 3. + GG1(3)
      PNSP(0) = PNSP(1) * 3. - PNSP(2) * 3. + PNSP(3)
      PNSM(0) = PNSM(1) * 3. - PNSM(2) * 3. + PNSM(3)
      FF2(0) = FF2(1) * 3. - FF2(2) * 3. + FF2(3)
      FG2(0) = FG2(1) * 3. - FG2(2) * 3. + FG2(3)
      GF2(0) = GF2(1) * 3. - GF2(2) * 3. + GF2(3)
      GG2(0) = GG2(1) * 3. - GG2(2) * 3. + GG2(3)
C
      FF1(NX) = FF1(NX-1) * 3. - FF1(NX-2) * 3. + FF1(NX-3)
      FG1(NX) = FG1(NX-1) * 3. - FG1(NX-2) * 3. + FG1(NX-3)
      GF1(NX) = GF1(NX-1) * 3. - GF1(NX-2) * 3. + GF1(NX-3)
      GG1(NX) = GG1(NX-1) * 3. - GG1(NX-2) * 3. + GG1(NX-3)
      PNSM(NX) = PNSM(NX-1) * 3. - PNSM(NX-2) * 3. + PNSM(NX-3)
      PNSP(NX) = PNSP(NX-1) * 3. - PNSP(NX-2) * 3. + PNSP(NX-3)
      FF2(NX) = FF2(NX-1) * 3. - FF2(NX-2) * 3. + FF2(NX-3)
      FG2(NX) = FG2(NX-1) * 3. - FG2(NX-2) * 3. + FG2(NX-3)
      GF2(NX) = GF2(NX-1) * 3. - GF2(NX-2) * 3. + GF2(NX-3)
      GG2(NX) = GG2(NX-1) * 3. - GG2(NX-2) * 3. + GG2(NX-3)


cpn 2001 <- For the DIS scheme
      if (iknl.eq.12) then
         DFF2(0) = DFF2(1) * 3. - DFF2(2) * 3. + DFF2(3)
         DFG2(0) = DFG2(1) * 3. - DFG2(2) * 3. + DFG2(3)
         DGF2(0) = DGF2(1) * 3. - DGF2(2) * 3. + DGF2(3)
         DGG2(0) = DGG2(1) * 3. - DGG2(2) * 3. + DGG2(3)
         DFF2(NX) = DFF2(NX-1) * 3. - DFF2(NX-2) * 3. + DFF2(NX-3)
         DFG2(NX) = DFG2(NX-1) * 3. - DFG2(NX-2) * 3. + DFG2(NX-3)
         DGF2(NX) = DGF2(NX-1) * 3. - DGF2(NX-2) * 3. + DGF2(NX-3)
         DGG2(NX) = DGG2(NX-1) * 3. - DGG2(NX-2) * 3. + DGG2(NX-3)
      endif                     !iknl.eq.12
cpn 2001 -> For the DIS scheme


cpn                            Calculate the integrals in the end-point terms
      if (iknl.ge.2) then       !Unpolarized case 
C                                           Int of PFF1, PNSP, and PNSM
C                                       Int of x * Kernel for FF2, GG1, and GG2
        RER = RERR *4.
        AFF1 (1) = AdzInt(PFF1_PN, D0, XV(1), AERR, RERR, ER1, IRT,2,2)
        
        DGG1     = NFL / 3.
        TMPG     = AdzInt(RGG1_pn, D0, XV(1), AERR, RERR, ER3, IRT,2,2)
        AGG1 (1) = TMPG + DGG1
        
        ANSM (1) = AdzInt(FNSM_pn, D0, XV(1), AERR, RER, ER2, IRT,2,2)
        ANSP (1) = AdzInt(FNSP_pn, D0, XV(1), AERR, RER, ER2, IRT,2,2)

        AER = AFF1(1) * RER
        AFF2 (1) = AdzInt(RFF2_PN, D0, XV(1),  AER, RER, ER2, IRT,2,2)

        AER = AGG1(1) * RER
        AGG2 (1) = AdzInt(RGG2_PN, D0, XV(1),  AER, RER, ER4, IRT,2,2)

cpn  >>     
        AFG2(Nx)=d0
        AGF2(Nx)=d0
            
        tem = AdzInt(PFG1_pn,XV(NX-1),D1, AERR, RERR, ER1, IRT,2,2)
        aer=abs(tem*rer) 
        AFG2 (NX-1) = AdzInt(PFG2_pn,XV(NX-1),D1, AER, RER, ER2, 
     >                       IRT,2,2)

        tem = AdzInt(PGF1_pn, XV(NX-1),D1, AERR, RERR, ER1, IRT,2,2)
        aer=abs(tem*rer) 
        AGF2(Nx-1)=AdzInt(PGF2_pn,XV(NX-1),D1, AER, RER, ER2, IRT,2,2)
cpn <<
    
        DO I2 = 2, NX-1

          TEM =AdzInt(PFF1_PN,XV(I2-1),XV(I2),AERR,RERR,ER1,IRT,2,2)
          AFF1(I2) = TEM + AFF1(I2-1)
          
          AER = ABS(TEM * RER)
          AFF2(I2)=AdzInt(RFF2_PN,XV(I2-1),XV(I2),AER,RER,ER2,IRT,2,2)
     >         +AFF2(I2-1)
          
          TEM      = AdzInt(RGG1_pn,XV(I2-1),XV(I2),AERR,RERR,ER3,IRT,2
     >         ,2)
          TMPG     = TMPG + TEM
          AGG1(I2) = TMPG + DGG1
        
          AER = ABS(TEM * RER)
          AGG2(I2)=AdzInt(RGG2_PN,XV(I2-1),XV(I2),AER,RER,ER4,IRT,2,2)
     >         +AGG2(I2-1)

          ANSP(I2)=AdzInt(FNSP_pn,XV(I2-1),XV(I2),AERR,RER,ER4,IRT,2,2)
     >         +ANSP(I2-1)
          ANSM(I2)=AdzInt(FNSM_pn,XV(I2-1),XV(I2),AERR,RER,ER4,IRT,2,2)
     >         +ANSM(I2-1)

cpn  >>    
          tem = AdzInt(PFG1_pn,XV(Nx-I2),XV(Nx-I2+1),AERR,RERR, 
     >         ER1,IRT,2,2)
          aer=abs(tem*rer) 
          AFG2(NX-I2) = AdzInt(PFG2_Pn,XV(NX-I2),XV(NX-I2+1),AER,RER,ER2
     >         ,IRT,2,2)+ AFG2(NX-I2+1)
        
          tem = AdzInt(PGF1_pn, XV(Nx-I2),XV(Nx-I2+1),AERR,RERR, 
     >         ER1,IRT,2,2)
          aer=abs(tem*rer) 
          AGF2(Nx-I2) = Adzint(PGF2_pn,XV(Nx-I2),XV(Nx-I2+1),AER,
     >         RER,ER2,IRT,2,2)+ AGF2(NX-I2+1) 
cpn <<

C                                               Compare with exact ans
C                   X = XV(I2)
C                   XLM1 = LOG(1./(1.-X))
C                   TEMF = (2.*XLM1 - X - X**2/2.) * 4./3.
C                   TEMG = (XLM1 - X**2 + X**3/3. - X**4/4.) * 6.

C                  PRINT '(/I6, 5(1PE14.3))', I2, X, AFF1(I2), TEMF, TMPG, TEMG
        enddo                   !i2

C                          Factor relevant to the full PNSP calculation
        ANSP(NX)=AdzInt(FNSP_pn,XV(NX-1),D1,AERR,RER,ERR, IRT,2,2)
     >       +ANSP(NX-1)
        ANSM(NX)=AdzInt(FNSM_pn,XV(NX-1),D1,AERR,RER,ERR, IRT,2,2)
     >       +ANSM(NX-1)
       
cpn                              Finally, calculate the integrals of
cpn                              x*P_FG(x) and  x*P_GF(x)  between 0 and 1
        zff2=-28./27.*Cf2+94./27.*CfCg -52./27.*CfTr*Nfl
        zgg2= 37./27.*CfTr*Nfl + 35./54.*CgTr*Nfl

cpn 2001 <-        
cpn 2001 In the DIS scheme, calculate the antiderivatives and end-point terms 
        if(iknl.eq.12) then
cpn                    In this part, I estimate the absolute error
cpn                     as the value of the integral, estimated
cpn                     using the trapesoid formula, times the
cpn                     relative error
          AER = 0.5*rdff2_pn(xv(1))*xv(1)*rer
          ADFF2 (1) = AdzInt(RDFF2_PN, D0, XV(1),AER, RER, ER2,IRT,2,2)
          AER = rdgf2_pn(xv(1))*xv(1)*rer
          ADGF2 (1) = AdzInt(RDGF2_PN, D0, XV(1),AER, RERR, ER,IRT,2,2)

          ADFG2(Nx)=d0
          ADGG2(Nx)=d0
          
cpn From the Mathematica study, I know that the absolute values of
cpn                     ADFG2 and ADGF2 do not exceed 2 and 20,
cpn                     respectively. Hence I just set AER=RERR
          
          AER = RERR
          ADFG2 (NX-1) = AdzInt(DelFG2_pn,XV(NX-1),D1, AER, RER, ER2, 
     >      IRT,2,2)

          ADGG2(Nx-1)=AdzInt(DelGG2_pn,XV(NX-1),D1, AER, RER, ER2, 
     >      IRT,2,2)
          
          do i2=2, nx-1
            AER = 0.5*(rdff2_pn(xv(i2-1)) + rdff2_pn(xv(i2)))*
     >        (xv(i2)-xv(i2-1))*rerr
            ADFF2(I2)=AdzInt(RDFF2_PN,XV(I2-1),XV(I2),AER,RER,ER2,IRT,2
     >        ,2)+ADFF2(I2-1)

            AER = 0.5*(rdgf2_pn(xv(i2-1)) + rdgf2_pn(xv(i2)))*
     >        (xv(i2)-xv(i2-1))*rerr
            ADGF2(I2)=AdzInt(RDGF2_PN,XV(I2-1),XV(I2),AER,RER,ER2,IRT,2
     >        ,2)+ADGF2(I2-1)

            AER = RERR
            ADFG2(NX-I2) = AdzInt(DelFG2_Pn,XV(NX-I2),XV(NX-I2+1),AER
     >        ,RER,ER2,IRT,2,2)+ ADFG2(NX-I2+1)
        
            ADGG2(Nx-I2) = Adzint(DelGG2_pn,XV(Nx-I2),XV(Nx-I2+1),AER,
     >         RER,ER2,IRT,2,2)+ ADGG2(NX-I2+1)             
          enddo                 !i2
          

cpn                    Alternatively, the antiderivatives can be
cpn                     calculated using their explicit expressions,
cpn                     which is more accurate than the numerical
cpn                     integration. I use the numerical integration
cpn                     because I want to calculate antiderivatives 
cpn                     with the same method for all values of iknl, 
cpn                     but feel free to use the more accurate
cpn                     calculation below. This choice does not make
cpn                     any difference as soon as the overall accuracy
cpn                     of the evolution package is about 1% 
c$$$          do i2 = 1,nx-1
c$$$            ADFF2(i2) = wdff2_pn(xv(i2))
c$$$            ADFG2(i2)= wdfg2_pn(xv(i2))
c$$$            ADGF2(i2)= wdgf2_pn(xv(i2))
c$$$            ADGG2(i2)= wdgg2_pn(xv(i2))
c$$$
c$$$          enddo

            
cpn                    Calculate the end-point terms
          ZDFF2 = (-11*CfCg)/36. - (4*Cf*Nfl*Tr)/9.
          ZDGF2 = - ZDFF2

        endif                   !iknl.eq.12
cpn 2001 ->
 
      elseif (iknl.eq.-2) then  !polarized case
C                                       Int of PNSP, PNSM
C                                       Int of FF2, FG2, GF2 and GG2 without
c                                       1./(1.-x) terms
        RER = RERR/20d0 

        ANSM (1) = AdzInt(FNSP_pn, D0, XV(1), AERR, RER, ER2,
     >       IRT,2,2)
        ANSP (1) = AdzInt(FNSM_pn, D0, XV(1), AERR, RER, ER2,
     >       IRT,2,2)

cpn3 >>
c        ANSM (Nx-1) = AdzInt(dFNSM_PN, XV(Nx-1),D1, AERR, RER, ER2, IRT
c     >       ,2,2)
c        aer=rer*ansm(Nx-1)        
c
c        if (er2.gt.aer) 
c     >       ansm(Nx-1)=AdzInt(dFNSM_PN, XV(Nx-1),D1, AER, RER, ER2, IRT
c     >       ,2,2)
cpn3 <<

cpn                             Calculate AFF2 using AERR as an absolute error.
cpn                             Then, refine the calculation by assuming
cpn                             AER=RER*AFF2
        AFF2(Nx)=d0
        AFG2(Nx)=d0
        AGF2(Nx)=d0
        AGG2(Nx)=d0
        
        AFF2 (Nx-1) = AdzInt(dPFF2_PN, XV(Nx-1),D1, AERR, RER, ER2, IRT
     >       ,2,2)
        aer=rer*aff2(Nx-1)        

        if (er2.gt.aer) 
     >       aff2(Nx-1)=AdzInt(dPFF2_PN, XV(Nx-1),D1, AER, RER, ER2, IRT
     >       ,2,2)


cpn                                     Repeat for the other splitting kernels
        AFG2 (Nx-1) = AdzInt(dPFG2_PN, XV(Nx-1),D1, AERR, RER, ER2, IRT
     >       ,2,2)
        aer=rer*afg2(Nx-1)        
        if (er2.gt.aer) afg2(Nx-1)=AdzInt(dPFG2_PN, XV(Nx-1),D1, AER,
     >       RER,ER2, IRT,2,2)
       
        AGF2 (Nx-1) = AdzInt(dPGF2_PN,XV(Nx-1),D1, AERR, RER, ER2, IRT
     >       ,2,2)
        aer=rer*agf2(Nx-1)        
        if (er2.gt.aer) agf2(Nx-1)=AdzInt(dPGF2_PN, XV(Nx-1),D1, AER,
     >       RER,ER2, IRT,2,2)
       

        AGG2 (Nx-1) = AdzInt(dPGG2_PN, XV(Nx-1),D1, AERR, RER, ER2, IRT
     >       ,2,2)
        aer=rer*agg2(Nx-1)        
        if (er2.gt.aer) agg2(Nx-1)=AdzInt(dPGG2_PN, XV(Nx-1),D1, AER,
     >       RER,ER2, IRT,2,2)
       

        DO I2 = 2, NX-1

          ANSP(I2)=AdzInt(FNSM_pn,XV(I2-1),XV(I2),AERR,RER,ER4
     >         ,IRT,2,2)+ANSP(I2-1)
         ANSM(I2)=AdzInt(FNSP_pn,XV(I2-1),XV(I2),AERR,RER,ER4
     >         ,IRT,2,2)+ANSM(I2-1)

cpn3>>          
c          tem = AdzInt(dFNSM_PN, XV(Nx-i2),xv(nx-i2+1), AERR,RER, ER2
c     >         ,IRT,2,2)
c          aer=rer*tem        
c          if (er2.gt.aer) 
c     >         tem=AdzInt(dFNSM_PN, XV(Nx-i2),xv(Nx-i2+1), AER, RER,
c     >         ER2,IRT,2,2)
c 
c         ANSM(Nx-i2) = tem +ANSM(Nx-i2+1)
cpn3<<

           tem = AdzInt(dPFF2_PN, XV(Nx-i2),xv(nx-i2+1), AERR,RER, ER2
     >         ,IRT,2,2)
          aer=rer*tem        
          if (er2.gt.aer) 
     >         tem=AdzInt(dPFF2_PN, XV(Nx-i2),xv(Nx-i2+1), AER, RER,
     >         ER2,IRT,2,2)
          AFF2(Nx-i2) = tem +AFF2(Nx-i2+1)
          
          tem = AdzInt(dPFG2_PN, XV(Nx-i2),xv(nx-i2+1), AERR,RER, ER2
     >         ,IRT,2,2)
          aer=rer*tem        
          if (er2.gt.aer) 
     >         tem=AdzInt(dPFG2_PN, XV(Nx-i2),xv(Nx-i2+1), AER, RER,
     >         ER2,IRT,2,2)
          AFG2(Nx-i2) = tem +AFG2(Nx-i2+1)
          
          tem = AdzInt(dPGF2_PN, XV(Nx-i2),xv(nx-i2+1), AERR,RER, ER2
     >         ,IRT,2,2)
          aer=rer*tem        
          if (er2.gt.aer) 
     >         tem=AdzInt(dPGF2_PN, XV(Nx-i2),xv(Nx-i2+1), AER, RER,
     >         ER2,IRT,2,2)
          AGF2(Nx-i2) = tem +AGF2(Nx-i2+1)
          
          tem = AdzInt(dPGG2_PN, XV(Nx-i2),xv(nx-i2+1), AERR,RER, ER2
     >         ,IRT,2,2)
          aer=rer*tem        
          if (er2.gt.aer) 
     >         tem=AdzInt(dPGG2_PN, XV(Nx-i2),xv(Nx-i2+1), AER, RER,
     >         ER2,IRT,2,2)
          AGG2(Nx-i2) = tem +AGG2(Nx-i2+1)
          
        enddo                   !i2


C                          Factor relevant to the full PNSP calculation
cpn                        We do not need them, really 
cpn         ANSP(NX)=AdzInt(FNSP_pn,XV(NX-1),D1,AERR,RER,ERR, IRT,2,2)
c     >        +ANSP(NX-1)
c         ANSM(NX)=AdzInt(FNSM_pn,XV(NX-1),D1,AERR,RER,ERR, IRT,2,2)
c     >        +ANSM(NX-1)

cpn                    End-points contributions copied from the papers
cpn                     of Vogelsang and Ellis, NP B, 475, 47 and CERN-TH/96-50
        zff2=Cf2*2.6525392+CfCg*3.133587545-CfTr*Nfl*2.359912089
        zgg2=(-4./3.)*CgTr*Nfl+(-1.)*CfTr*Nfl+6.27283737*Cg2   
                                !6.27283737=3*zeta(3)+8/3
      endif                     !iknl
cpn 
C                      The numeric value is from the integral
C                      2*P_qqb between 0 to 1, or (13-2*Pi^2+8*Zeta[3])/2
      ZQQB=1.43862321154902*(Cf2-0.5*CfCg)

cpn                    Check the accuracy of the calculation
c      open(unit=37,file='/home/nadolsky/wkt/dat/ansm.dat',status
c     >     ='unknown')
c      do ix=1,nx-1,2
c        write(37,*) xv(ix),ansm(ix)
c      enddo 
cpn                    
c      close(37)
c      stop
      

cpn <<


C                                                 ----------------------------
C          Set up kernel functions for direct use in 2nd-order calculations
C
C                1 . . . . . . . . .               1 . . . . . . . . .
C                . 1 . . . . . . . .               . 1 . . . . . . . .
C                . . 1             .               . . 1             .
C                .     1    FF2    .               .     1     GG2   .
C                .       1         .               .       1         .
C                .         .       .               .         .       .
C                .   FG2     1     .               .   GF2     1     .
C                .             1   .               .             1   .
C                . . . . . . . . 1 .               . . . . . . . . 1 .
C                . . . . . . . . . 1               . . . . . . . . . 1
C
C       The same applies to the array PNS containing PNSM \ PNSP

      DO 21 IX = 1, NX-1
        X = XV(IX)
C                                        Number of points in the I-th integral
        NP = NX - IX + 1
        IS = NP
C                                   Evaluate the FG & GF kernel along diagonal
        XG2 = (LOG(1./(1.-X)) + 1.) ** 2
        FFG (IS, IS) = FG2(NX) * DXTZ(I) * XG2
        GGF (IS, IS) = GF2(NX) * DXTZ(I) * XG2
        PNS (IS, IS) = PNSM(NX) * DXTZ(I)

cpn 2001 For the DIS scheme <-
        if (iknl.eq.12) then
          DFFG (IS, IS) = DFG2(NX) * DXTZ(I) * XG2
          DGGF (IS, IS) = DGG2(NX) * DXTZ(I) * XG2
        endif !iknl.eq.12
cpn 2001 ->

        DO 31 KZ = 2, NP
          IY = IX + KZ - 1
          IT = NX - IY + 1
          XY = X / XV(IY)
          XM1 = 1.- XY
          XG2 = (LOG(1./XM1) + 1.) ** 2
           
          Z  = ZZ (IX, IY)
C                                           Quantities needed for interpolation
          TZ = (Z + DZ) / DZ
          IZ = TZ
          IZ = MAX (IZ, 0)
          IZ = MIN (IZ, NX-1)
          DT = TZ - IZ
C                 DXTZ is the Jacobian due to the change of variable dx/x -> dZ
C                                   2nd order terms -- by interpolation
          TEM = (FF2(IZ) * (1.- DT) + FF2(IZ+1) * DT) / XM1 / XY
          FFG (IX, IY) = TEM * DXTZ(IY)
          
          TEM = (FG2(IZ) * (1.- DT) + FG2(IZ+1) * DT) * XG2 / XY
          FFG (IS, IT) = TEM * DXTZ(IY)
          
          TEM = (GF2(IZ) * (1.- DT) + GF2(IZ+1) * DT) * XG2 / XY
          GGF (IS, IT) = TEM * DXTZ(IY)

          TEM = (GG2(IZ) * (1.- DT) + GG2(IZ+1) * DT) / XM1 / XY
          GGF (IX, IY) = TEM * DXTZ(IY)
          
          TEM = (PNSP(IZ) * (1.- DT) + PNSP(IZ+1) * DT) / XM1
          PNS (IX, IY) = TEM * DXTZ(IY)
          
          TEM = (PNSM(IZ) * (1.- DT) + PNSM(IZ+1) * DT) / XM1
          PNS (IS, IT) = TEM * DXTZ(IY)
          
cpn 2001 <-
          if (iknl.eq.12) then
            TEM = (DFF2(IZ) * (1.- DT) + DFF2(IZ+1) * DT) / XM1/ XY
            DFFG (IX, IY) = TEM * DXTZ(IY)
            
            TEM = (DFG2(IZ) * (1.- DT) + DFG2(IZ+1) * DT) * XG2 / XY
            DFFG (IS, IT) = TEM * DXTZ(IY)

            TEM = (DGF2(IZ) * (1.- DT) + DGF2(IZ+1) * DT) / XM1/XY
            DGGF (IX, IY) = TEM * DXTZ(IY)

            TEM = (DGG2(IZ) * (1.- DT) + DGG2(IZ+1) * DT) * XG2 / XY
            DGGF (IS, IT) = TEM * DXTZ(IY)
          endif                 !iknl.eq.12
cpn 2001 ->

 31     CONTINUE
 21   CONTINUE
      RETURN
C                        ****************************
      END !stupkl

      function  TFF1_PN(ZZ)
      implicit NONE 
      double precision PFF1_pn, x,xfrmz,tff1_PN,zz,dxdz
      external xfrmz

      x=xfrmz(zz)
      tff1_PN=pff1_pn(x)*dxdz(zz)
      Return
      End!TFF1_pn

      function  VFF1_pn (XX)
      implicit none
      double precision x,xx,vff1_pn,tem
      X = XX
      TEM = (1.+ X**2) / (1.- X)
      VFF1_pn = TEM *4./3.
      RETURN
      END

cpn **************************************************************************
cpn                    Functions containing polarized kernels
cpn **************************************************************************
cpn---------------------------------------------------------------------------

      SUBROUTINE WTUPD (UPD, NTL, NDAT, IRET)
C                                                   -=-=- wtupd
C
C                       The I/O operation is made a stand-alone subprogram
C                       so that fast execution is achieved with block
C                       data transfer for the actual size of the array
C                       (rather than the declared size in the main program).
      DOUBLE PRECISION UPD
      DIMENSION UPD (NTL)
C
      WRITE (NDAT, '(1pE13.5, 5E13.5)', IOSTAT=IRET) UPD
C
      RETURN
C                       ----------------------------
C
      ENTRY RDUPD (UPD, NTL, NDAT, IRET)
C
      READ (NDAT, *, IOSTAT=IRET) UPD
C                             To check read-in result against saved file.
C      Print '(1pE13.5, 5E13.5)', UPD
C
      RETURN
C                        ****************************
      END

      Subroutine X2Zconv (Nx, Xx, Zz, Xmn, Ixz)
C                                                   -=-=- x2zconv

C				Convert x to z:	Ixz = 1 : x --> z ;
C						            2 : z --> x

C                        Xmn is the Xmin parameter in the conversion
C                        formula (along with Xcr which has 1.5 as default).
C                        Xmn = 1e-3 or 1e-4 are good choices;
C                        Z = (0,1) for x = (Xmn, 1); if x<xmn then Z<0
C                        (which is perfectly ok if you don't mind!)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Dimension Xx(Nx), Zz(Nx)

      Call SetEvl
C						Recall current value of xmin
      Call ParPdf(2, 'Xmin', Xm, Ir)
C						Set Xmin to input grid value
      Call ParPdf(1, 'Xmin', Xmn, Ir)

      Do 11 Ix = 1, Nx
         If     (Ixz .eq. 1) Then
            Zz(Ix) = zfrmx (Xx(Ix))
         Elseif (Ixz .eq. 2) Then
            Xx(Ix) = xfrmz (Zz(Ix))
         Else
            Stop 'Ixz must be 1 or 2 in X2Zconv!'
         Endif
   11 Continue

C						Restore Xmin to original value
      Call ParPdf(1, 'Xmin', Xm, Ir)

      Return
C                       ***********************8
      End
	

cpn *************************************************************************
cpn                    Added June, 2000
cpn *************************************************************************
C------------------------------------------------------------------------------

      SUBROUTINE XARRAY
C                                                   -=-=- xarray
C
C       Given NX, this routine fills the x-arrays of the common blocks
C       related to the x-variable for the use of various other routines.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
C
      PARAMETER (D0 = 0.0, D10=10.0)
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      Character Msg*80
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / VARBAB / GB(NDG, NDH, MXX), H(NDH, MXX, M1:M2)
      COMMON / HINTEC / GH(NDG, MXX)
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)
C
      DIMENSION G1(NDG,NDH), G2(NDG,NDH), A(NDG)
C
      DATA F12, F22, F32 / 1D0, 1D0, 1D0 /
      DATA (G1(I,NDH), G2(I,1), I=1,NDG) / 0.0,0.0,0.0,0.0,0.0,0.0 /
      DATA PUNY / 1D-30 /
C
C                                                ----------------------------
C                                                   Change of variable X <--> Z
C                        Range of Z is [0, 1] for X in [XM, 1] and I in [1, NX]
C                                                 Use evenly spaced points in Z
      XV(0) = 0D0
      DZ = 1D0 / (NX-1)
      DO 10 I = 1, NX - 1
         Z = DZ * (I-1)
         ZV(I) = Z
         X = XFRMZ (Z)
C                               DXDZ is the Jacobian for the change of variable
         DXTZ(I) = DXDZ(Z) / X
C                                                               Fill x - arrays
         XV (I)  = X
         XA(I, 1) = X
         XA(I, 0) = LOG (X)
         DO 20 L = L1, L2
          IF (L .NE. 0 .AND. L .NE. 1)  XA(I, L) = X ** L
   20    CONTINUE
   10 CONTINUE
C                   Over-write XV(1) by hand, since inversion gives error of ~ 10^-4
         XV(1) = Xmin
C                                                         Fill last point, I=Nx
         XV(NX) = 1D0
         ZV(Nx) = 1D0
         DXTZ(NX) = DXDZ(1.D0)
         DO 21 L = L1, L2
            XA (NX, L) = 1D0
   21    CONTINUE
         XA (NX, 0) = 0D0
C                                                       Fill ELY = Log (1. - x)

      DO 11 I = 1, NX-1
         ELY(I) = LOG(1D0 - XV(I))
   11 CONTINUE
C                     Log(1-x) is infinite at x=1, for the purpose of numerical
C                       Calculations, the following extrapolated value is used
C                       to avoid an artificial discontinuity.
C
       ELY(NX) = 3D0* ELY(NX-1) - 3D0* ELY(NX-2) + ELY(NX-3)
C                                                  ----------------------------
C                                    Matrix elements for 2nd order calculations
C
C                                     1 . . . . . . . . .
C                                     . 1 . . . . . . . .
C                                     . . 1             .
C                                     .     1   Z (X/Y) .
C                                     .       1         .
C                                     .         .       .
C                                     .  X/Y      1     .
C                                     .             1   .
C                                     . . . . . . . . 1 .
C                                     . . . . . . . . . 1
      DO 17 IX = 1, NX
      ZZ (IX, IX) = 1.
      DO 17 IY = IX+1, NX
         XY = XV(IX) / XV(IY)
         ZZ (IX, IY) = ZFRMX (XY)
         ZZ (NX-IX+1, NX-IY+1) = XY
   17 CONTINUE
C                                                ------------------------------
C                                     Start of x - loop to compute ceefficients
C                                   for integrals used in INTEG, AMOM, ... etc.
      DO 30 I = 1, NX-1
C                                                  "F" matrix {a(i)} --> {f(i)}
      IF (I .NE. NX-1) THEN
        F11 = 1D0/XV(I)
        F21 = 1D0/XV(I+1)
        F31 = 1D0/XV(I+2)
        F13 = XV(I)
        F23 = XV(I+1)
        F33 = XV(I+2)
C                                                        Determinant for matrix
C
        DET = F11*F22*F33 + F21*F32*F13 + F31*F12*F23
     >      - F31*F22*F13 - F21*F12*F33 - F11*F32*F23
        IF (ABS(DET) .LT. PUNY) THEN
           Msg='Determinant close to zero; will be arbitrarily set to:'
           CALL WARNR(IWRN, NWRT, Msg, 'DET', PUNY, D0, D0, 0)
           DET = PUNY
        EndIf
C                                            Compute "G" matrix -- inverse of F
C                                             G1 is only needed from I=2 and on
        G2(1,2) = (F22*F33 - F23*F32) / DET
        G2(1,3) = (F32*F13 - F33*F12) / DET
        G2(1,4) = (F12*F23 - F13*F22) / DET
C
        G2(2,2) = (F23*F31 - F21*F33) / DET
        G2(2,3) = (F33*F11 - F31*F13) / DET
        G2(2,4) = (F13*F21 - F11*F23) / DET
C
        G2(3,2) = (F21*F32 - F22*F31) / DET
        G2(3,3) = (F31*F12 - F32*F11) / DET
        G2(3,4) = (F11*F22 - F12*F21) / DET
C                                               Compute coefficients for HINTEG
        B2 = LOG (XV(I+2)/XV(I))
        B3 = XV(I) * (B2 - 1.) + XV(I+2)
        GH (1,I) = B2 * G2 (2,2) + B3 * G2 (3,2)
        GH (2,I) = B2 * G2 (2,3) + B3 * G2 (3,3)
        GH (3,I) = B2 * G2 (2,4) + B3 * G2 (3,4)
      EndIf
C                                         "G-bar" is the "average" of G1 and G2
        DO 51 J = 1, NDH
           DO 52 L = 1, NDG
C                                                                First interval
              IF     (I .EQ. 1) THEN
                 GB(L,J,I) = G2(L,J)
C                                                                 last interval
              ElseIF (I .EQ. NX-1) THEN
                 GB(L,J,I) = G1(L,J)
C                                                        intermidiate intervals
              Else
                 GB(L,J,I) = (G1(L,J) + G2(L,J)) / 2D0
              EndIf
   52      CONTINUE
   51   CONTINUE
C                                                            Compute "A" matrix
        DO 35 MM = M1, M2
           DO 40 K = 1, NDG
             KK = K + MM - 2
             IF (KK .EQ. 0) THEN
               A(K) = XA(I+1, 0) - XA(I, 0)
             Else
               A(K) = (XA(I+1, KK) - XA(I, KK)) / DBLE(KK)
             EndIf
   40      CONTINUE
C                                                            Compute "H" matrix
           DO 41 J = 1, NDH
             TEM = 0
             DO 43 L = 1, NDG
               TEM = TEM + A(L) * GB(L,J,I)
   43        CONTINUE
             H(J,I,MM) = TEM
   41      CONTINUE
   35   CONTINUE
C                                               ------------------------------
C                                                         Initialize G1 matrix
      DO 42 J = 1, NDG
        DO 44 L = 1, NDG
           G1(L,J) = G2(L,J+1)
   44 CONTINUE
   42 CONTINUE
C
   30 CONTINUE
C                                    End of x - loop to calculate coefficients
C                                               ------------------------------
      LSTX = .TRUE.
      RETURN
C                        ****************************
      END
      FUNCTION XF0 (IX)
C                                                   -=-=- xf0
C                      The following functions return (X, Q) lattice values
C                      and the values of PDF at these lattice points;
C                      The "0" set is for the standard set of evolved distr.;
C                      The "1" set is for the 'alternate' set.
C                                         Value of X at site IX
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX= MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / QARAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / PEVLDT / UPD(MXPQX), KF, Nelmt
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)

      DATA IW1, IW2 / 2 * 0 /
      save

      IF (IX .LE. 0) THEN
        CALL WARNI(IW1,NWRT,'IX out of range in XF0', 'IX',IX,1,NX,1)
        XF0 = 0.
      ElseIF (IX .LE. NX) THEN
        XF0 = XV(IX)
      Else
        CALL WARNI(IW1,NWRT,'IX out of range in XF0', 'IX',IX,1,NX,1)
        XF0 = 0.
      EndIf

      RETURN
C                        ****************************
      Entry Zf0 (Ix)

      IF (IX .LE. 0) THEN
        CALL WARNI(IW1,NWRT,'IX out of range in Zf0', 'IX',IX,1,NX,1)
        Zf0 = 0.
      ElseIF (IX .LE. NX) THEN
        Zf0 = ZV(IX)
      Else
        CALL WARNI(IW1,NWRT,'IX out of range in Zf0', 'IX',IX,1,NX,1)
        Zf0 = 0.
      EndIf

      RETURN
C                        ****************************
      ENTRY QF0 (IQ)
C                                   Value of QQ at site IQ
      IF (IQ .LT. 0) THEN
        CALL WARNI(IW2,NWRT,'IQ out of range in QF0', 'IQ',IQ,0,NT,1)
        QF0 = 0.
      ElseIF (IQ .LE. NT) THEN
        QF0 = QV(IQ)
      Else
        CALL WARNI(IW2,NWRT,'IQ out of range in QF0', 'IQ',IQ,0,NT,1)
        QF0 = 0.
      EndIf

      RETURN
C                        ****************************

      ENTRY PDF0 (IPRTN, IX, IQ)
C                                   Value of PDF at lattice point (IX, IQ)
      JFL = IPRTN + (KF - 2) / 2
      J0  = JFL * (NT+1) * (NX+1)

      IUPD = J0 + IQ * (NX + 1) + IX + 1
      PDF0 = UPD(IUPD)

      RETURN
C                        ****************************
      END

      FUNCTION XF1 (IX)
C                                                   -=-=- xf1
C                                          Value of X at site IX - other set
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
C
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (MXPN = MXF * 2 + 2)
      PARAMETER (MXQX = MXQ * MXX,   MXPQX = MXQX * MXPN)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / X1ARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / Q1RAY1 / QINI,QMAX, QV(0:MXQ),TV(0:MXQ), NT,JT,NG
      COMMON / E1LPAC / AL, IKNL, IPD0, IHDN, NfMx
      COMMON / P1VLDT / UPD(MXPQX), KF, Nelmt
      COMMON / XYARAY / ZZ(MXX, MXX), ZV(0:MXX)

      save
      DATA IW1, IW2 / 2 * 0 /

      IF (IX .LE. 0) THEN
        CALL WARNI(IW1,NWRT,'IX out of range in XF1', 'IX',IX,1,NX,1)
        XF1 = 0.
      ElseIF (IX .LE. NX) THEN
        XF1 = XV(IX)
      Else
        CALL WARNI(IW1,NWRT,'IX out of range in XF1', 'IX',IX,1,NX,1)
        XF1 = 0.
      EndIf

      RETURN
C                        ****************************
      Entry Zf1 (Ix)

      IF (IX .LE. 0) THEN
        CALL WARNI(IW1,NWRT,'IX out of range in Zf1', 'IX',IX,1,NX,1)
        Zf1 = 0.
      ElseIF (IX .LE. NX) THEN
        Zf1 = ZV(IX)
      Else
        CALL WARNI(IW1,NWRT,'IX out of range in Zf1', 'IX',IX,1,NX,1)
        Zf1 = 0.
      EndIf

      RETURN
C                        ****************************
      ENTRY QF1 (IQ)
C                                   Value of QQ for site IQ  - alternate set
      IF (IQ .LE. 0) THEN
        CALL WARNI(IW2,NWRT,'IQ out of range in QF1', 'IQ',IQ,1,NT,1)
        QF1 = 0.
      ElseIF (IQ .LE. NT) THEN
        QF1 = XV(IQ)
      Else
        CALL WARNI(IW2,NWRT,'IQ out of range in QF1', 'IQ',IQ,1,NT,1)
        QF1 = 0.
      EndIf

      RETURN
C                        ****************************

      ENTRY PDF1 (IPRTN, IX, IQ)
C                                   Value of PDF on the lattice - alternate set
      JFL = IPRTN + (KF - 2) / 2
      J0  = JFL * (NT+1) * (NX+1)

      IUPD = J0 + IQ * (NX + 1) + IX + 1
      PDF1 = UPD(IUPD)

      RETURN
C                        ****************************
      END
C                                                          =-=-= Runevl
      FUNCTION XFRMZ (Z)
C                                                   -=-=- xfrmz
C                                       Invert the equation Z = ZFRMX (X)
C                                       XM and XC are the two parameters in
C                                       the transformation formula
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (MXX = 204, Rer = 1.d-4)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      COMMON / INVERT / ZA
C
      EXTERNAL ZFXL
C
      DATA ZLOW, ZHIGH, Tem, IWRN2 / -10.0, 1.00002, 1.d-4, 0 /
      Save Tem
 
C                              Mathematically, the full range of Z is allowed.
C                                Physically, Z < 1. is the only constraint.
C                                Note however,  Z < 0  corresponds to X < Xmin.
C                                We allow the range ZLOW < Z < 1.
      EPS = TEM * RER
      ZA = Z
C
      IF (Z .LE. ZHIGH .AND. Z .GT. ZLOW) THEN
          XLA = LOG (XMIN) * 1.5
          XLB = 0.00001
          TEM = ZBRNT (ZFXL, XLA, XLB, EPS, IRT)
      Else
C
        CALL WARNR (IWRN2, NWRT, 'Z out of range in XFRMZ, X set=0.',
     >              'Z', Z, ZLOW, ZHIGH, 1)
        TEM = 0
      EndIf
C
      XFRMZ = EXP(TEM)
C
      RETURN
C                        ****************************
      END

      FUNCTION XFZ (Z)
C                                                   -=-=- xfz
C                           Interpolates XV(I) to return value of X for given Z
C                                      Same as XFRMZ (Z), but should be faster.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0)
      PARAMETER (MXX = 204, MXQ = 25, MXF = 6)
      PARAMETER (M1=-3, M2=3, NDG=3, NDH=NDG+1, L1=M1-1, L2=M2+NDG-2)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / VARIBX / XA(MXX, L1:L2), ELY(MXX), DXTZ(MXX)
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
      DATA KX / 0 /
C
      IF (.NOT.LSTX) CALL XARRAY
      IF (KX .NE. NX) THEN
          KX = NX
          DZ = 1./(NX-1)
      EndIf

      XFZ = FINTRP (XV(1), D0, DZ, NX-1, Z, 1, ER, IR)

      RETURN
C                                                  ---------- Jacobian -------
      ENTRY DXTZF(Z)
C                                Jacobian function for (dx/x) <--> DXTZF(z) dz
      IF (.NOT.LSTX) CALL XARRAY
      IF (KX .NE. NX) THEN
          KX = NX
          DZ = 1./(NX-1)
      EndIf

      DXTZF = FINTRP (DXTZ, D0, DZ, NX-1, Z, 1, ER, IR)

      RETURN
C                        ****************************
      END
      function xLi(n,x)
C                                                   -=-=- xli
CPN-------------------------------------------------------------------------
C     Polylogarithm function for real x,-1<x<1.For n=2 the precision is about
C     12 significant figures. For n>2 the precision is smaller, but still
C     better than 4 digits.
C     Dilogarithm part is the translation of Carl Schmidt's C code.

      implicit NONE

      integer NCUT, i,n,m3
      double precision xLi,Out,x,pi2by6,zeta3,c1,c2
      double precision r,xt,L,xln1m

      parameter (m3=8)
      dimension c1(2:m3),c2(2:m3)

      data NCUT/27/
      data c1/0.75,-0.5833333333333333,0.454861111111111,
     >        -0.3680555555555555,0.3073611111111111,
     >        -0.2630555555555555,0.2294880243764172/

      data c2/-0.5,0.5,-0.4583333333333333,0.416666666666666,
     >        -0.3805555555555555,0.35,-0.3241071428571428/

      data zeta3,pi2by6 /1.20205690315959,1.64493406684823/

      L=0.0
      i=0
      r=1.0

C     Check if x lies in the correct range
      if (abs(x).gt.r) then
        print *,'Li: x out of range (-1,1) , x=',x
        STOP
      endif

      if (n.lt.0) then
       print *,'The polilogarithm Li undefined for n=',n
       STOP
      elseif (n.eq.0) then
       Out=x/(1d0-x)
      elseif (n.eq.1) then
       Out=-dlog(1-x)
      elseif (n.eq.2) then
                                                !Calculate dilogarithm
                                                !separately for x<0.5 and x>0.5
      if (x.ge.(-0.5).and.x.le.0.5) then

         do while(i.le.NCUT)
       	  i=i+1
          r=r*x

          L=L+r/i/i
         enddo
         Out=L
       elseif (x.eq.0) then
         Out=0d0
       elseif(x.gt.0.5) then !n.eq.2,x>0.5
         xt = 1.0-x
         L = pi2by6 - dlog(x)*dlog(xt)

         do while(i.le.NCUT)
          i=i+1
          r=r*xt

          L=L-r/i/i
         enddo
         Out=L
       elseif (x.lt.(-0.5)) then
         xt=-x/(1d0-x)
         L=-0.5*dlog(1-x)**2

         do while (i.le.NCUT)
          i=i+1
          r=r*xt
          L=L-r/i/i
         enddo
         Out=L
       endif
      elseif (n.eq.3.and.x.ge.0.8) then !use the expansion of Li3 near x=1
       L=zeta3+pi2by6*dlog(x)
       xt=(1d0-x)
       xln1m=dlog(xt)

       do i=2,m3
        L=L+(c1(i)+c2(i)*xln1m)*xt**i
       enddo
       Out=L
      else !n>3 or x=3,x<0.8

         do while(i.le.NCUT)
          i=i+1
          r=r*x

          L=L+r/dble(i)**dble(n)
         enddo
         Out=L

      endif

      xLi=Out
C                      **********************
      End ! xLi

C                                                          =-=-= Evlutl
      FUNCTION ZFRMX (XX)
C                                                   -=-=- zfrmx
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LSTX
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MXX = 204)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / XXARAY / XCR, XMIN, XV(0:MXX), LSTX, NX
C
      DATA IWRN1, HUGE, TINY / 0, 1.E35, 1.E-35 /
C
      F(X) = (XCR-XMIN) * LOG (X/XMIN) + LOG (XCR/XMIN) * (X-XMIN)
      D(X) = (XCR-XMIN) / X          + LOG (XCR/XMIN)
C
C    7 IF (XMIN .GE. XCR) THEN
C         PRINT *,'XMIN > XCR not allowed in ZFRMX; XMIN, XCR =',XMIN,XCR
C         STOP
C      EndIf
C
      X = XX
C
      IF (X .GE. XMIN) THEN
         TEM = F(X) / F(D1)
      ElseIF (X .GE. D0) THEN
         X = MAX (X, TINY)
         TEM = F(X) / F(D1)
      Else
         CALL WARNR(IWRN1, NWRT, 'X out of range in ZFRMX, Z set to 99.'
     >             , 'X', X, TINY, HUGE, 1)
         TEM = 99.
         STOP
      EndIf
C
      ZFRMX = TEM
C
      RETURN
C                                                ----------------------------
      ENTRY DZDX (XX)
C
C   77 IF (XMIN .GE. XCR) THEN
C         PRINT *,'XMIN > XCR not allowed in DZDX; XMIN, XCR =',XMIN,XCR
C         STOP
C      EndIf
C
      X = XX
      IF (X .GE. XMIN) THEN
         TEM = D(X) / F(D1)
      ElseIF (X .GE. D0) THEN
         X = MAX (X, TINY)
         TEM = D(X) / F(D1)
C         CALL WARNR(IWRN2, NWRT, 'X < Xmin in DZDX', 'X',X,D0,D1,0)
      Else
         CALL WARNR(IWRN1, NWRT, 'X out of range in DZDX, Z set to 99.'
     >             , 'X', X, TINY, HUGE, 1)
         TEM = 99.
         STOP
      EndIf
C
      DZDX = TEM

      RETURN
C                        ****************************
      END

      FUNCTION ZFXL (XL)
C                                                   -=-=- zfxl
C                                                  Used in XFRMZ for inversion
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / INVERT / ZA
C
      X = EXP(XL)
      TT = ZFRMX (X) - ZA
      ZFXL = TT
C
      RETURN
C                        ****************************
      END

