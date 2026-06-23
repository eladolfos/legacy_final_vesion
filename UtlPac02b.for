      Subroutine AChSh (In, Ch, Sh, y)
C                                                   -=-=- achsh

C===========================================================================
C GroupName: Miscutl
C Description:   Misc. functions
C ListOfFiles: achsh bta euler gama gamma loggam pi trngle chept nextun
C===========================================================================

C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Miscutl.f,v 1.2 98/03/09 01:11:07 wkt Exp $
C #Log:	Miscutl.f,v $
c Revision 1.2  98/03/09  01:11:07  wkt
c NextUn.f slightly modified.
c 
c Revision 1.1  97/12/21  21:19:19  wkt
c Initial revision
c 

C                                   Inverts the Cosh(y) / Sinh(y) functions
C                                   and returns Sinh(y) / Cosh(y) along y.
C     In = 1 : input is Ch = Cosh(y)
C          2 : input is Sh = Sinh(y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Data Sml / 1D-10 /

      If     (In .Eq. 1) Then
         If (Ch .Ge. 1D0+Sml) Then
           Sh = SqRt (Ch **2 - 1D0)
         ElseIf (Ch .Ge. 1D0-Sml) Then
           Ch = 1D0
           Sh = 0D0
         Else
           Print '(/A, 1pE11.3)'
     >     , ' Cosh value out of range in Achsh: Ch =', Ch
           Stop
         Endif
      ElseIf (In .Eq. 2) Then
         Ch = SqRt (Sh **2 + 1D0)
      Else
         Print '(A, I5)', 'In must be 1,2  in AChSh; I =', In
         Stop
      EndIf

      Ey = Ch + Sh
      y  = Log (Ey)

      Return
C                           ***************
      End

      FUNCTION AdfInt 
C                                                   -=-=- adfint
     >(F, FF, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

C===========================================================================
C GroupName: Adfint
C Description: Adaptive Integration with function argument
C ListOfFiles: adfint adfcal adfspl infgnd sgfint toftal
C===========================================================================
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adfint.f,v 1.1 97/12/21 21:18:32 wkt Exp $
C #Log:	Adfint.f,v $
c Revision 1.1  97/12/21  21:18:32  wkt
c Initial revision
c 

C                             Integration with integrand function F
C                             of the form F(x,FF) where FF is an external
C                             function F needs for its calculation.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F, FF
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / AdfWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / AdfWrk /
      DATA SMLL, Sml / 1E-20, 1E-12 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in AdfInt!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in AdfInt call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      DDX = B - A
      If (DDX .Le. 0D0) Then
        AdfInt = 0D0
        Ier = 4
        If (DDX .Lt. 0D0) 
     >     Print '(/A/)', 'B < A in AdzInt; check limits!!'
        Return
      ElseIf (DDX .Le. Sml) Then
        AdfInt = F(A + DDX/2, FF) * DDX
        Ier = 5
        Return
      EndIf

      NUMINT = 3
      DX = DDX/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1), FF)
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2., FF)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I), FF)
             ELSE
               IB = I
               FB = F(B-DX/2., FF)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I), FF)
          ENDIF
          CALL AdfCal(F, FF, I)
   10     CONTINUE
       CALL Toftal
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET)
     $               CALL AdfSpl(F, FF, I, IER)
   40             CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      AdfInt = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE AdfSpl (F, FF, I, IER)
C                                                   -=-=- adfspl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F, FF
      PARAMETER (MAXINT = 1000)
      COMMON / AdfWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / AdfWrk /
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
     
      CALL AdfCal (F, FF, I)
      CALL AdfCal (F, FF, NUMINT)
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
 
      SUBROUTINE AdfCal (F, FF, I)
C                                                   -=-=- adfcal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F, FF
      PARAMETER (MAXINT = 1000)
      COMMON / AdfWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / AdfWrk /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4., FF)

        CALL SgfInt (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4., FF)
        CALL SgfInt (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W, FF)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                      Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE SgfInt (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sgfint

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
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
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
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
     
      FUNCTION InfGnd (X, FX)
C                                                   -=-=- infgnd
C                    Return number of distinct points used in AdfInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)

      COMMON / AdfWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      Dimension X(MAXINT), FX(MAXINT)

      SAVE / AdfWrk /

       Call Sort2 (NumInt, U, Fu)
       Call Sort2 (NumInt, V, Fv)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

       InfGnd = NUMINT + 1

      RETURN
C                        ****************************
      END

      SUBROUTINE Toftal
C                                                   -=-=- toftal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / AdfWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / AdfWrk /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adgint
      FUNCTION AdgInt 
C                                                   -=-=- adgint
     >(F, FF, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)

C===========================================================================
C GroupName: Adgint
C Description: second copy of adfint
C ListOfFiles: adgint adgcal adgspl inggnd sggint togtal
C===========================================================================
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adgint.f,v 1.1 97/12/21 21:18:57 wkt Exp $
C #Log:	Adgint.f,v $
c Revision 1.1  97/12/21  21:18:57  wkt
c Initial revision
c 

C                             Integration with integrand function F
C                             of the form F(x,FF) where FF is an external
C                             function F needs for its calculation.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F, FF
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / AdgWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / AdgWrk /
      DATA SMLL / 1E-20 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in AdgInt!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in AdgInt call', 
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
               FU(1) = F(U(1), FF)
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2., FF)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I), FF)
             ELSE
               IB = I
               FB = F(B-DX/2., FF)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I), FF)
          ENDIF
          CALL AdgCal(F, FF, I)
   10     CONTINUE
       CALL Togtal
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET)
     $               CALL AdgSpl(F, FF, I, IER)
   40             CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      AdgInt = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE AdgSpl (F, FF, I, IER)
C                                                   -=-=- adgspl

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F, FF
      PARAMETER (MAXINT = 1000)
      COMMON / AdgWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / AdgWrk /
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
     
      CALL AdgCal (F, FF, I)
      CALL AdgCal (F, FF, NUMINT)
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
 
      SUBROUTINE AdgCal (F, FF, I)
C                                                   -=-=- adgcal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F, FF
      PARAMETER (MAXINT = 1000)
      COMMON / AdgWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / AdgWrk /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4., FF)

        CALL SggInt (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4., FF)
        CALL SggInt (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W, FF)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE SggInt (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sggint

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
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
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
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
     
      SUBROUTINE Togtal
C                                                   -=-=- togtal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / AdgWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / AdgWrk /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Simpgaus
      FUNCTION IngGnd (X, FX)
C                                                   -=-=- inggnd
C                    Return number of distinct points used in AdgInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)

      COMMON / AdgWrk / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
      Dimension X(MAXINT), FX(MAXINT)

      SAVE / AdgWrk /

       Call Sort2 (NumInt, U, Fu)
       Call Sort2 (NumInt, V, Fv)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

       IngGnd = NUMINT + 1

      RETURN
C                        ****************************
      END

      FUNCTION ADAPTIVEINTEGRATE(F, A, B, RERR, AERR)
      INTEGER MAXRECURSIONS, NSUBINTERVALS, MAXERROR, K
      INTEGER FINDMAX
      PARAMETER (MAXRECURSIONS = 100)
      REAL*8 ERR(MAXRECURSIONS), LOWERBOUND(MAXRECURSIONS)
      REAL*8 UPPERBOUND(MAXRECURSIONS), ERROR, INTLIST(MAXRECURSIONS)
      REAL*8 AREA, ERRORSUM, ERRORMAX, ERRORBOUND
      REAL*8 LOWER1, LOWER2, UPPER1, UPPER2, ERROR1, ERROR2
      REAL*8 AREA1, AREA2, AREA12, ERROR12
      REAL*8 F, ADAPTIVEINTEGRATE
      REAL*8 A,B, RERR, AERR
      EXTERNAL F, FINDMAX

      ADAPTIVEINTEGRATE = GAUSSQUAD(F,A,B,ERROR)
      ERRORBOUND = MAX(AERR,
     >                 RERR*DABS(ADAPTIVEINTEGRATE))
      IF(ERROR < ERRORBOUND) RETURN

      AREA = ADAPTIVEINTEGRATE
      ERRORSUM = ERROR
      ERRORMAX = ERROR

      MAXERROR = 1

      DO NSUBINTERVALS = 2,MAXSUBINTERVALS
        LOWER1 = LOWERBOUND(MAXERROR)
        UPPER2 = UPPERBOUND(MAXERROR)

        UPPER1 = (LOWER1+UPPER2)*0.5
        LOWER2 = UPPER1

        AREA1 = GAUSSQUAD(F, LOWER1, UPPER1, ERROR1)
        AREA2 = GAUSSQUAD(F, LOWER2, UPPER2, ERROR2)

        AREA12 = AREA1+AREA2
        ERROR12 = ERROR1 + ERROR2
        ERRORSUM = ERRORSUM + ERROR12 - ERRORMAX
        AREA = AREA + AREA12 - INTLIST(MAXERROR)

        INTLIST(MAXERROR) = AREA1
        INTLIST(NSUBINTERVALS) = LOWER1

        IF(ERROR2 > ERROR1) THEN
            LOWERBOUND(NSUBINTERVALS) = LOWER1
            LOWERBOUND(MAXERROR) = LOWER2
            UPPERBOUND(NSUBINTERVALS) = UPPER1
            INTLIST(MAXERROR) = AREA2
            INTLIST(NSUBINTERVALS) = AREA1
            ERR(MAXERROR) = ERROR2
            ERR(NSUBINTERVALS) = ERROR1
        ELSE
            LOWERBOUND(NSUBINTERVALS) = LOWER2
            UPPERBOUND(MAXERROR) = UPPER1
            UPPERBOUND(NSUBINTERVALS) = UPPER2
            ERR(MAXERROR) = ERROR1
            ERR(NSUBINTERVALS) = ERROR2
        ENDIF

        MAXERROR = FINDMAX(ERR)

        ERRORBOUND = MAX(AERR, RERR*DABS(AREA))

        IF(ERRORSUM .le. ERRORBOUND) EXIT
            
      ENDDO

      IF(NSUBINTERVALS .eq. MAXSUBINTERVALS) PRINT*, "WARNING: Reached",
     > "maximum recursions without converging!"

      ADAPTIVEINTEGRATE = 0

      DO K=1,NSUBINTERVALS
        ADAPTIVEINTEGRATE = ADAPTIVEINTEGRATE + INTLIST(K)
      ENDDO

      END

      FUNCTION GAUSSQUAD(F, A, B, ERR)
      INTEGER J, JJ
      REAL*8 HALFLENGTH, CENTER, FCENTER, RESULTGAUSS, RESULTKRONROD
      REAL*8 ABSCISSA, F1, F2, FSUM, F, A, B, ERR
      REAL*8 RESULTMEANKRONROD, ABSDIFFINTEGRAL
      REAL*8 KRONRODWGTS(8), GAUSSWGTS(4), ABSC(7)
      EXTERNAL F
      DATA KRONRODWGTS /0.022935322010529,0.063092092629979,0.104790010322250
     >,0.140653259715525,0.169004726639267,0.190350578064785,
     >0.204432940075298,0.209482141084728/
      DATA GAUSSWGTS /0.129484966168870,0.279705391489277,0.381830050505119
     >,0.417959183673469/
      DATA ABSC /0.991455371120813,0.949107912342759,0.864864423359769,
     >0.741531185599394,0.586087235467691,0.405845151377397,
     >0.207784955007898/

      HALFLENGTH = (B-A)*0.5
      CENTER = (A+B)*0.5
      FCENTER = F(CENTER)

      RESULTGAUSS = GAUSSWGTS(4)*FCENTER
      RESULTKRONROD = KRONRODWGTS(8)*FCENTER

      DO J = 1,3
        JJ = 2*J
        ABSCISSA = HALFLENGTH*ABSC(JJ)

        F1 = F(CENTER-ABSCISSA)
        F2 = F(CENTER+ABSCISSA)
        FSUM = F1+F2

        RESULTGAUSS = RESULTGAUSS + GAUSSWGTS(J)*FSUM
        RESULTKRONROD = RESULTKRONROD + KRONRODWGTS(JJ)*FSUM
      ENDDO

      DO J = 0,3
        JJ = 2*J+1
        ABSCISSA = HALFLENGTH*ABSC(JJ)

        F1 = F(CENTER-ABSCISSA)
        F2 = F(CENTER+ABSCISSA)
        FSUM = F1+F2

        RESULTKRONROD = RESULTKRONROD + KRONRODWGTS(JJ)*FSUM
      ENDDO

      RESULTMEANKRONROD = RESULTKRONROD*0.5
      ABSDIFFINTEGRAL=KRONRODWGTS(8)*(DABS(FCENTER-RESULTMEANKRONROD))
      DO J = 1, 6
        ABSCISSA = HALFLENGTH*ABSC(J)
        ABSDIFFINTEGRAL =
     >    ABSDIFFINTEGRAL+DABS(F(CENTER-ABSCISSA)-RESULTMEANKRONROD)*
     >    KRONRODWGTS(J)
        ABSDIFFINTEGRAL =
     >    ABSDIFFINTEGRAL+DABS(F(CENTER+ABSCISSA)-RESULTMEANKRONROD)*
     >    KRONRODWGTS(J)
      ENDDO

      ABSDIFFINTEGRAL = ABSDIFFINTEGRAL*HALFLENGTH

      ERR =
     >    ABSDIFFINTEGRAL*MIN(1.0,(200*DABS(RESULTKRONROD-RESULTGAUSS)
     >   *HALFLENGTH/ABSDIFFINTEGRAL)**1.5)

      IF(ERR.NE.ERR) THEN
          GAUSSQUAD = 0
          RETURN
      ENDIF

      GAUSSQUAD = RESULTKRONROD*HALFLENGTH

      END

      FUNCTION ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adz2nt
 
C===========================================================================
C GroupName: Adz2nt
C Description: second copy of adzint
C ListOfFiles: adz2nt adz2pl adz2al int2sz sgl2nt tot2lz
C=========================================================================== 
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adz2nt.f,v 1.1 97/12/21 21:19:00 wkt Exp $
C #Log:	Adz2nt.f,v $
c Revision 1.1  97/12/21  21:19:00  wkt
c Initial revision
c 

C List of GLOBAL Symbols

C     FUNCTION   ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ2PL (F, I, IER)
C     SUBROUTINE ADZ2AL (F,I)
C     SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT2LZ
C     FUNCTION   INT2SZ (X, FX)
C
C     COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
C    > ICTA, ICTB, NUMINT, IB
C                   ------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA SMLL / 1E-20 /
      data Ikmx /10/
      data fac0 /2d0/ 

	double precision ERS_OLD
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZ2NT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
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
          CALL ADZ2AL(F,I)
   10     CONTINUE
       CALL TOT2LZ
C                                                   Adaptive procedure:
cLai revised to ensure the accuracy requested and improve the speed
       Do
 20       TARGET = ABS(AERR) + ABS(RERR * RES)
          ERS_OLD=ERS
          if(ERS <= TARGET .or. Ier>0) goto 30
          NUMOLD = NUMINT
          facNum = Sqrt(dble(NUMOLD)*fac0)
          Ik=0
          Do
             DO I = 1, NUMOLD
                IF (ERR(I)*facNum .GT. TARGET) CALL ADZ2PL(F,I,IER)
             Enddo              ! I
             if(NUMOLD .ne. NUMINT) then
                goto 20
             else
                Ik=Ik+1
                facNum=facNum*fac0
                if(Ik>Ikmx) goto 30
             endif
          Enddo !do loop if NUMOLD==NUMINT
       enddo !do loop if .not.(ERS<=Target .or. Ier>0)
 30    continue
      ADZ2NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZ2PL (F, I, IER)
C                                                   -=-=- adz2pl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
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
     
      CALL ADZ2AL (F, I)
      CALL ADZ2AL (F, NUMINT)
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
 
      SUBROUTINE ADZ2AL (F,I)
C                                                   -=-=- adz2al
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZ2RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL2NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL2NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
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

      SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sgl2nt

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
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
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
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
     
      SUBROUTINE TOT2LZ
C                                                   -=-=- tot2lz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adz3nt
      FUNCTION INT2SZ (X, FX)
C                                                   -=-=- int2sz
C                    Return number of distinct points used in AdzInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)

      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      Dimension X(MAXINT), FX(MAXINT)

      SAVE / ADZ2RK /

       Call Sort2 (NumInt, U, Fu)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue

       Call Sort2 (NumInt, V, Fv)
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

      INT2SZ = NUMINT + 1
      RETURN
C                        ****************************
      END
C
      FUNCTION adz3nt (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adz3nt
 
C===========================================================================
C GroupName: Adz3nt
C Description: second copy of adzint
C ListOfFiles: adz3nt adz3pl adz3al int3sz sgl3nt tot3lz
C=========================================================================== 
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adz3nt.f,v 1.2 98/03/09 01:24:08 wkt Exp $
C #Log:	Adz3nt.f,v $
c Revision 1.2  98/03/09  01:24:08  wkt
c typo fixed
c 

C----------------------------------------------------------
C List of GLOBAL Symbols

C     FUNCTION   adz3nt (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE adz3pl (F, I, IER)
C     SUBROUTINE adz3al (F,I)
C     SUBROUTINE sgl3nt (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE tot3lz
C     FUNCTION   int3sz (X, FX)
C
C     COMMON / ADZ3RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
C    > ICTA, ICTB, NUMINT, IB
C                   ------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZ3RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ3RK /
      DATA SMLL / 1E-20 /
      data Ikmx /10/
      data fac0 /2d0/ 

	double precision ERS_OLD
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in adz3nt!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in adz3nt call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in adz3nt call', 
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
          CALL adz3al(F,I)
   10     CONTINUE
       CALL tot3lz
C                                                   Adaptive procedure:
cLai revised to ensure the accuracy requested and improve the speed
       Do
 20       TARGET = ABS(AERR) + ABS(RERR * RES)
          ERS_OLD=ERS
          if(ERS <= TARGET .or. Ier>0) goto 30
          NUMOLD = NUMINT
          facNum = Sqrt(dble(NUMOLD)*fac0)
          Ik=0
          Do
             DO I = 1, NUMOLD
                IF (ERR(I)*facNum .GT. TARGET ) CALL ADZ3PL(F,I,IER)
             Enddo              ! I
             if(NUMOLD .ne. NUMINT) then
                goto 20
             else
                Ik=Ik+1
                facNum=facNum*fac0
                if(Ik>Ikmx) goto 30
             endif
          Enddo !do loop if NUMOLD==NUMINT
       enddo !do loop if .not.(ERS<=Target .or. Ier>0)
 30    continue
      adz3nt = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE adz3pl (F, I, IER)
C                                                   -=-=- adz3pl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C_yfu change the MAXINT
C      PARAMETER (MAXINT = 2000)
      COMMON / ADZ3RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ3RK /
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
     
      CALL adz3al (F, I)
      CALL adz3al (F, NUMINT)
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
 
      SUBROUTINE adz3al (F,I)
C                                                   -=-=- adz3al
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ3RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZ3RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL sgl3nt (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL sgl3nt (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
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

      SUBROUTINE sgl3nt (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sgl3nt

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
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
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
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
     
      SUBROUTINE tot3lz
C                                                   -=-=- tot3lz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ3RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ3RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adfint
      FUNCTION int3sz (X, FX)
C                                                   -=-=- int3sz
C                    Return number of distinct points used in AdzInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)

      COMMON / ADZ3RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      Dimension X(MAXINT), FX(MAXINT)

      SAVE / ADZ3RK /

       Call Sort2 (NumInt, U, Fu)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue

       Call Sort2 (NumInt, V, Fv)
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

      int3sz = NUMINT + 1
      RETURN
C                        ****************************
      END
C
      FUNCTION ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adzint

C===========================================================================
C GroupName: Adzint
C Description: adaptive integration
C ListOfFiles: adzint adzcal adzspl intusz sglint totalz
C===========================================================================
C                                  Authors: Wu-Ki Tung and John C. Collins
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adzint.f,v 1.1 97/12/21 21:19:04 wkt Exp $
C #Log:	Adzint.f,v $
c Revision 1.1  97/12/21  21:19:04  wkt
c Initial revision
c 

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ (X, FX)
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ------------------------

C     Adaptive integration routine which allows the integrand to be 
C     indeterminant at the lower and/or the upper ends of integration. 

C     Can self-adjust to any integrable singularity at the ends and compute 
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).
 
C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm 
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)
 
C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
CError code is IER: 0 :  o.k.
C                1 :  maximum calls to function reached before the 
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C                4 :  Error on Limits : B < A ; zero result returned.
C                5 :  Range of integration DX zero or close to roundoff
C                     returns DX * F(A+DX/2)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL,Sml / 1d-20,1d-12 /
      data Ikmx /10/
      data fac0 /2d0/ 

      double precision ERS_OLD

      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZINT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      DDX = B - A
      If (DDX .Le. 0D0) Then
        AdzInt = 0D0
        Ier = 4
        If (DDX .Lt. 0D0) 
     >     Print '(/A/)', 'B < A in AdzInt; check limits!!'
        Return
      ElseIf (DDX .Le. Sml) Then
        AdzInt = F(A + DDX/2) * DDX
        Ier = 5
        Return
      EndIf

      NUMINT = 3
      DX = DDX/ NUMINT
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
          CALL ADZCAL(F,I)
   10     CONTINUE
       CALL TOTALZ
C                                                   Adaptive procedure:
cLai revised to ensure the accuracy requested and improve the speed
       Do
 20       TARGET = ABS(AERR) + ABS(RERR * RES)
	  ERS_OLD=ERS
          if(ERS <= TARGET .or. Ier>0) goto 30
          NUMOLD = NUMINT
          facNum = Sqrt(dble(NUMOLD)*fac0)
          Ik=0
          Do
             DO I = 1, NUMOLD
                IF (ERR(I)*facNum .GT. TARGET) CALL ADZSPL(F,I,IER)
             Enddo              ! I
             if(NUMOLD .ne. NUMINT) then
                goto 20
             else
                Ik=Ik+1
                facNum=facNum*fac0
                if(Ik>Ikmx) goto 30
             endif
          Enddo !do loop if NUMOLD==NUMINT
       enddo !do loop if .not.(ERS<=Target .or. Ier>0)
 30    continue
      ADZINT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE ADZSPL (F, I, IER)
C                                                   -=-=- adzspl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
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
     
      CALL ADZCAL (F, I)
      CALL ADZCAL (F, NUMINT)
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
 
      SUBROUTINE ADZCAL (F,I)
C                                                   -=-=- adzcal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZWRK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGLINT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGLINT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
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

      SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sglint

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
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
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
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
     
      SUBROUTINE TOTALZ
C                                                   -=-=- totalz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adz2nt
      FUNCTION INTUSZ (X, FX)
C                                                   -=-=- intusz

C                    Return number of distinct points used in AdzInt
C                  Also returns the x- and F(x)-values (ie. the integrand)
C                       at these points (perhaps for plotting purposes).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      Dimension X(MAXINT), FX(MAXINT)

      SAVE / ADZWRK /

       Call Sort2 (NumInt, U, Fu)
       Do 10 Ix = 1, NUMINT
         X(Ix)  =  U(Ix)
         FX(Ix) = Fu(Ix)
   10  Continue

       Call Sort2 (NumInt, V, Fv)
       X (NumInt + 1) = V (NumInt)
       FX(NumInt + 1) = Fv(NumInt)

      INTUSZ = NUMINT + 1

      RETURN
C                        ****************************
      END

      FUNCTION ASK(QUERY)
C                                                   -=-=- ask
C
C===========================================================================
C GroupName: Askwrn
C Description:  a few Input prompt utilities and warning utilities
C ListOfFiles: ask inquir inqure warni warnr
C===========================================================================
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Askwrn.f,v 1.2 98/03/09 01:09:10 wkt Exp $
C #Log:	Askwrn.f,v $
c Revision 1.2  98/03/09  01:09:10  wkt
c minor revision to take away character string across lines.
c 
c Revision 1.1  97/12/21  21:19:06  wkt
c Initial revision
c 

      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      LOGICAL ASK
      CHARACTER QUEND*11,  QUERY*(*), CH*1
C
      PARAMETER (QUEND= ' (Y OR N)? ' )
C
      CALL RTB (QUERY, LEN)
1     WRITE(NOUT, 90) QUERY(1:LEN), QUEND
      READ(NIN, 91) CH
      CALL UPCASE (CH)
      ASK = .FALSE.
      IF ( (CH.EQ.'Y') .OR. (CH.EQ.' ') ) THEN
         ASK = .TRUE.
         RETURN
      ELSE IF (CH.EQ.'N') THEN
         RETURN
      ELSE
         WRITE(NOUT, 92)
      ENDIF
      GOTO 1
 90     FORMAT ('$', 3A)
 91     FORMAT (A1)
 92     FORMAT (' BAD ANSWER--TRY AGAIN ')
C               *************************
      END

      FUNCTION BTA (X, Y)
C                                                   -=-=- bta
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      BTA = GAMA (X) * GAMA (Y) / GAMA (X+Y)
C
      RETURN
      END
C                       ****************************
C
      SUBROUTINE CHEPT (I, Msg, Var)
C                                                   -=-=- chept
C
C                               THIS ROUTINE PROVIDES CHECK-POINTS IN FORTRAN
C                               PROGRAMS FOR DEBUGGING PURPOSES
C
      Integer I
      Character*(*) Msg
      Double Precision Var

      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      WRITE (NOUT, 900) I, Msg, Var
  900 FORMAT ( 1X, 'CHECK-POINT  ', I2, 2x, A, G13.3 /)
C
      RETURN
C               *************************
      END
C
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
C                                                   -=-=- covsrt
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
C                        ****************************
      END
C             -------------------------------------------------
      FUNCTION EULER()
C                                                   -=-=- euler
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (E=0.57721566)
      ENTRY EULERF()
      EULER = E
      RETURN
C                       ****************************
      END

      FUNCTION FINTR2 (FF,  X0, DX, NX, Y0, DY, NY,  XV, YV,  ERR, IR)
C                                                   -=-=- fintr2
 
C     Two variable interpolation  --  Double Precision Version
C
c     Given array FF defined on evenly spaced lattice (0:Nx; 0:Ny), this
c     function routine calculates an interpolated value for the function at an
C     interrior point XV, YV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;  similarly for y:
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2
 
C     It uses (MX-1)th ((MY-1)th) order polynomial fits to MX (MY) neighoring
C     points in the X (Y) direction.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX=3, MY=3, MDX=MX/2, MDY=MY/2, AMDX=MDX, AMDY=MDY)
      DIMENSION FF (0:NX, 0:NY), XX(MX), YY(MY), ZZ(MX,MY)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA XX, YY / 0.0, 1.0, 2.0, 0.0, 1.0, 2.0 /
      DATA  IW1, IW2, IW3, IW4, IW5, IW6 / 6 * 0 /
C                                                               Initiation
      IR = 0
      X = XV
      Y = YV
      ERR = 0.
      ANX = NX
      ANY = NY
      FINTR2 = 0.
C                                                               Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTR2.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN (NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTR2.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      IF (NY .LT. 1) THEN
         CALL WARNI(IW2, NWRT, 'NY < 1, error in FINTR2.',
     >              'NY', NY, 1, 256, 1)
         IR = 3
         RETURN
      ELSE
         MNY = MIN(NY+1, MY)
      ENDIF
C
      IF (DY .LE. 0) THEN
         CALL WARNR(IW4, NWRT, 'DY < 0, error in FINTR2.',
     >              'DY', DY, D0, D1, 1)
         IR = 4
         RETURN
      ENDIF
C                                        Set up interpolation point and 3-array
      XM = X0 + DX * NX
      IF (X .LT. X0 .OR. X .GT. XM) THEN
        CALL WARNR(IW5,NWRT,
     >    'X out of range in FINTR2, Extrapolation used.',
     >    'X', X, X0, XM, 1)
      ENDIF
C
      TX = (X - X0) / DX
      IF (TX .LE. AMDX) THEN
        IX = 0
      ELSEIF (TX .GE. ANX-AMDX) THEN
        IX = NX - MX + 1
      ELSE
        IX = TX - MDX + 1
      ENDIF
      DDX = TX - IX
C                                  Set up interpolation point and 3-array for Y
      YM = Y0 + DY * NY
      ytiny = 1.e-5
c ytiny prevents complaints whose source is in small round off errors
c and not anything real
      IF (Y .LT. (Y0-ytiny) .OR. Y .GT. (YM+ytiny)) THEN
        CALL WARNR(IW6,NWRT,
     >    'Y out of range in FINTR2, Extrapolation used.',
     >    'Y',Y,Y0,YM,1)
      ENDIF
C                                        Set up interpolation point and 3-array
      TY = (Y - Y0) / DY
      IF (TY .LE. AMDY) THEN
        IY = 0
      ELSEIF (TY .GE. ANY-AMDY) THEN
        IY = NY - MY + 1
      ELSE
        IY = TY - MDY + 1
      ENDIF
      DDY = TY - IY
C                                     POLIN2 is taken from "Numerical Recipe"
      DO 15 JX = 1, MNX
      DO 15 JY = 1, MNY
      ZZ (JX, JY) = FF (IX+JX-1, IY+JY-1)
   15 CONTINUE
 
      CALL POLIN2 (XX, YY, ZZ, MNX, MNY, DDX, DDY, TEM, ERR)
 
      FINTR2 = TEM
C
      RETURN
C
      END
C                       ****************************
 
C                                                          =-=-= Plot4
      FUNCTION FINTRP (FF,  X0, DX, NX,  XV,  Iint, ERR, IR)
C                                                   -=-=- fintrp
C   > FF : function array
C   > X0 : starting value of x-array
C   > Dx : spacing of x-array (assumed to be equally spaced)
C   > Nx : number of array points
C   > XV : value of x where the function is to be evaluated (interpolation point)
C   > Iint : Switch to choose between interpolation method:
C            1 : 3-point (quadratic) polynomial interpolation, using PolInt
C            2 : 3-point rational function interpolation, using RatInt
C   > Err  : estimated error
C   > Return error code (see the actual code below).

C===========================================================================
C GroupName: Fintrp
C Description:  several interpolation routines
C ListOfFiles: fintrp fntr2p fintr2
C===========================================================================
 
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Fintrp.f,v 1.2 98/03/09 01:10:48 wkt Exp $
C #Log:	Fintrp.f,v $

C Revision Sunday, April 4, 2004 at 12:06: Added switch Iint to choose between
C          PolInt or RatInt

c Revision 1.2  98/03/09  01:10:48  wkt
c minor revision to take away character string across lines.
c 
c Revision 1.1  97/12/21  21:19:13  wkt
c Initial revision
c 

C     Single variable interpolation  --  Double Precision Version
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MX = 3)
C
c     Given array FF defined on evenly spaced lattice (0:Nx), this function
c     routine calculates an interpolated value for the function at an
C     interrior point XV.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;
C
C            TX         0                Tx              Nx
C            IZ         0  1  2   I-1  I         Nx-2    Nx
C                       |--|--| ... |--|--|--| ... |--|--|
C             X        X0                X               XM
C            XX                        0  1  2
 
C     It uses 2nd order polynomial fit from three neighoring points.
 
      DIMENSION FF (0:NX), XX(MX)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DATA SML, XX / 1.D-5,  0., 1.0, 2.0 /
      DATA  IW1, IW3, IW5 / 3 * 0 /
C                                                                     Initiation
      IR = 0
      X = XV
      ERR = 0.
      ANX = NX
      FINTRP = 0.
C                                                                          Tests
      IF (NX .LT. 1) THEN
         CALL WARNI(IW1, NWRT, 'Nx < 1, error in FINTRP.',
     >              'NX', NX, 1, 256, 1)
         IR = 1
         RETURN
      ELSE
         MNX = MIN(NX+1, MX)
      ENDIF
C
      IF (DX .LE. 0) THEN
         CALL WARNR(IW3, NWRT, 'DX < 0, error in FINTRP.',
     >              'DX', DX, D0, D1, 1)
         IR = 2
         RETURN
      ENDIF
C
      XM = X0 + DX * NX
      IF (X .LT. X0-SML .OR. X .GT. XM+SML) THEN
        CALL WARNR(IW5,NWRT,
     >     'X out of range in FINTRP, Extrapolation used.',
     >     'X',X,X0,XM,1)
      IR = 3
      ENDIF
C                                        Set up interpolation point and 3-array
      TX = (X - X0) / DX
      IF (TX .LE. 1.) THEN
        IX = 0
      ELSEIF (TX .GE. ANX-1.) THEN
        IX = NX - 2
      ELSE
        IX = TX
      ENDIF
      DDX = TX - IX
C                                        POLINT/RATINT are taken from "Numerical Recipe"
      If (Iint.eq.1) then
        CALL POLINT (XX, FF(IX), MNX, DDX, TEM, ERR)
      Elseif (Iint.eq.2) then
        CALL RATINT (XX, FF(IX), MNX, DDX, TEM, ERR)
      Else
      Print *, 'Iint variable in Fintrp.f must be 1 or 2; Iint = ', Iint
        Stop 
      EndIf
 
      FINTRP = TEM
C
      RETURN
C
      END
C                       ****************************
 
      FUNCTION FNTR2P (FF, X0,DX,NX, P0,DP,NP, X, P)
C                                                   -=-=- fntr2p
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
c     Given function ff defined on the two-dimensional array, this function
c     routine calculates an interpolated value for the function at any
c     given point.
c     The lowest value for the variable x is x0, the mesh-size is dx and the
c     array-size is 0:NX;  the minimum value for the variable p is p0, the
c     increment it dp, and the array-size is 0:NP.
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DIMENSION FF(0:NX, 0:NP)
C
      IF (X.EQ.X0+DX*NX.AND.P.EQ.P0+DP*NP) THEN
                        FNTR2P = FF(NX,NP)
                        RETURN
      END IF
C
      TX = (X - X0) / DX
      IX = TX
      DDX = TX - IX
      IF (IX .LT. 0 .OR. IX .GT. NX) GOTO 90
C
      TP = (P - P0) / DP
      IP = TP
      DDP = TP - IP
      IF (IP .LT. 0 .OR. IP .GT. NP) GOTO 90
C
      FNTR2P = (1.-DDX) * (1.-DDP) * FF(IX  , IP  )
     >         +    DDX  * (1.-DDP) * FF(IX+1, IP  )
     >         +(1.-DDX) *     DDP  * FF(IX  , IP+1)
     >         +    DDX  *     DDP  * FF(IX+1, IP+1)
C
      RETURN
C
   90 WRITE (NOUT, 990) X, P
  990 FORMAT (' value(s) of x, and/or p out of range in interp'/
     > / '  x =', F10.3, ' p =', F10.3 / )
C
      STOP
C                       *******************************
      END

      FUNCTION GAMA (X)
C                                                   -=-=- gama
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Z = X
      CALL LOGGAM (Z, U)
      GAMA = EXP (U)
C
      RETURN
C                       ****************************
      END
C
      SUBROUTINE GAMMA(Z,G)
C                                                   -=-=- gamma
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      X = Z
      CALL LOGGAM(X,U)
      G = EXP(U)
C
      RETURN
C                       ****************************
      END
C
      FUNCTION GAMMLN(XX)
C                                                   -=-=- gammln
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
C                        ****************************
      END

      FUNCTION GAMMQ(A,X)
C                                                   -=-=- gammq
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IF(X.LT.0..OR.A.LE.0.) PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
      ENDIF
      RETURN
C                        ****************************
      END

      FUNCTION GausInt(F,XL,XR,AERR,RERR,ERR,IRT)
C                                                   -=-=- gausint

C                                           Adptive Gaussian integration
C     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
        DIMENSION XLIMS(100), R(93), W(93)
        INTEGER PTR(4),NORD(4),NIN,NOUT,NWRT

        external f

        COMMON/IOUNIT/NIN,NOUT,NWRT
        DATA PTR,NORD/4,10,22,46,  6,12,24,48/
        DATA R/.2386191860,.6612093865,.9324695142,
     1 .1252334085,.3678314990,.5873179543,.7699026742,.9041172563,
     1 .9815606342,.0640568929,.1911188675,.3150426797,.4337935076,
     1 .5454214714,.6480936519,.7401241916,.8200019860,.8864155270,
     1 .9382745520,.9747285560,.9951872200,.0323801710,.0970046992,
     1 .1612223561,.2247637903,.2873624873,.3487558863,.4086864820,
     1 .4669029048,.5231609747,.5772247261,.6288673968,.6778723796,
     1 .7240341309,.7671590325,.8070662040,.8435882616,.8765720203,
     1 .9058791367,.9313866907,.9529877032,.9705915925,.9841245837,
     1 .9935301723,.9987710073,.0162767488,.0488129851,.0812974955,
     1 .1136958501,.1459737146,.1780968824,.2100313105,.2417431561,
     1 .2731988126,.3043649444,.3352085229,.3656968614,.3957976498,
     1 .4254789884,.4547094222,.4834579739,.5116941772,.5393881083,
     1 .5665104186,.5930323648,.6189258401,.6441634037,.6687183100,
     1 .6925645366,.7156768123,.7380306437,.7596023411,.7803690438,
     1 .8003087441,.8194003107,.8376235112,.8549590334,.8713885059,
     1 .8868945174,.9014606353,.9150714231,.9277124567,.9393703398,
     1 .9500327178,.9596882914,.9683268285,.9759391746,.9825172636,
     1 .9880541263,.9925439003,.9959818430,.9983643759,.9996895039/
        DATA W/.4679139346,.3607615730,.1713244924,
     1 .2491470458,.2334925365,.2031674267,.1600783285,.1069393260,
     1 .0471753364,.1279381953,.1258374563,.1216704729,.1155056681,
     1 .1074442701,.0976186521,.0861901615,.0733464814,.0592985849,
     1 .0442774388,.0285313886,.0123412298,.0647376968,.0644661644,
     1 .0639242386,.0631141923,.0620394232,.0607044392,.0591148397,
     1 .0572772921,.0551995037,.0528901894,.0503590356,.0476166585,
     1 .0446745609,.0415450829,.0382413511,.0347772226,.0311672278,
     1 .0274265097,.0235707608,.0196161605,.0155793157,.0114772346,
     1 .0073275539,.0031533461,.0325506145,.0325161187,.0324471637,
     1 .0323438226,.0322062048,.0320344562,.0318287589,.0315893308,
     1 .0313164256,.0310103326,.0306713761,.0302999154,.0298963441,
     1 .0294610900,.0289946142,.0284974111,.0279700076,.0274129627,
     1 .0268268667,.0262123407,.0255700360,.0249006332,.0242048418,
     1 .0234833991,.0227370697,.0219666444,.0211729399,.0203567972,
     1 .0195190811,.0186606796,.0177825023,.0168854799,.0159705629,
     1 .0150387210,.0140909418,.0131282296,.0121516047,.0111621020,
     1 .0101607705,.0091486712,.0081268769,.0070964708,.0060585455,
     1 .0050142027,.0039645543,.0029107318,.0018539608,.0007967921/
        DATA TOLABS,TOLREL,NMAX/1.E-35,5.E-4,100/
C
C
        TOLABS=AERR
        TOLREL=RERR
     
        GausInt=0.
        NLIMS=2
        XLIMS(1)=XL
        XLIMS(2)=XR
C
10      AA=(XLIMS(NLIMS)-XLIMS(NLIMS-1))/2D0
        BB=(XLIMS(NLIMS)+XLIMS(NLIMS-1))/2D0
        TVAL=0.
        DO 15 I=1,3
15      TVAL=TVAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        TVAL=TVAL*AA
        DO 25 J=1,4
        VAL=0.
        DO 20 I=PTR(J),PTR(J)-1+NORD(J)
20      VAL=VAL+W(I)*(F(BB+AA*R(I))+F(BB-AA*R(I)))
        VAL=VAL*AA
        TOL=MAX(TOLABS,TOLREL*ABS(VAL))
        IF (ABS(TVAL-VAL).LT.TOL) THEN
                GausInt=GausInt+VAL
                NLIMS=NLIMS-2
                IF (NLIMS.NE.0) GO TO 10
                RETURN
                END IF
25      TVAL=VAL
        IF (NMAX.EQ.2) THEN
                GausInt=VAL
                RETURN
                END IF
        IF (NLIMS.GT.(NMAX-2)) THEN
                WRITE(NOUT,50) GausInt,NMAX,BB-AA,BB+AA
                RETURN
                END IF
        XLIMS(NLIMS+1)=BB
        XLIMS(NLIMS+2)=BB+AA
        XLIMS(NLIMS)=BB
        NLIMS=NLIMS+2
        GO TO 10
C
50      FORMAT (' GausInt FAILS, GausInt,NMAX,XL,XR=',G15.7,I5,2G15.7)
C                        ****************************
        END
 
C                                                          =-=-= Fintrp
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
C                                                   -=-=- gaussj
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
C                        ****************************
      END

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
C                                                   -=-=- gcf
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
C                        ****************************
      END

      SUBROUTINE GLFIT
C                                                   -=-=- glfit
     > (X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,FUNCS)

C===========================================================================
C GroupName: Nurcpe
C Description:  Modules adapted from Numerical Recipe
C ListOfFiles: glfit covsrt gammln gammq gser gcf gaussj linfit sort2 ratint polint polin2 romint zbrnt
C===========================================================================
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Nurcpe.f,v 1.2 98/03/09 01:11:45 wkt Exp $
C #Log:	Nurcpe.f,v $
c Revision 1.2  98/03/09  01:11:45  wkt
c Set function=0 for error returns in RomInt and Zbrnt.
c 
c Revision 1.1  97/12/21  21:19:22  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                             Fit data to linear combination of given functions
C                                            From LFIT of Numerical Recipes
      PARAMETER (MMAX=50)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *    COVAR(NCVM,NCVM),BETA(MMAX),AFUNC(MMAX)

      external funcs

      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
          PAUSE 'Improper set in LISTA in GLFIT'
        ENDIF
12    CONTINUE
      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA in GLFIT'
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE
      DO 18 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
      DO 24 I=1,NDATA
        CALL FUNCS(X(I),AFUNC,MA)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE
      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
C                        ****************************
      END

      SUBROUTINE GRDATD (NCUR, NX, NY, NUMPT, XPT, YPT)
C                                                   -=-=- grdatd
C                                       Loads Array with data points
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*1 A, SYMBOL(10)
      COMMON / GRAPH / A(130,130)
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
      DATA SYMBOL /'1','2','3','4','5','6','7','8','9','0'/
      DO 20 I = 1, NUMPT
         IX =  (XPT(I) - XMIN) / DX + 0.5
         IF (IX.LT.1) IX = 1
         IF (IX.GT.NX) IX = NX
         DO 10 J = 1, NCUR
            YY = YPT(I, J)
            IY = (YMAX - YY) / DY + 0.5
            IF (IY.LT.1) IY = 1
            IF (IY.GT.NY) IY = NY
     
            A(IX,IY) = SYMBOL(J)
10       CONTINUE
20    CONTINUE
      RETURN
C                       ****************************
      END
C
      SUBROUTINE GRNULL(NX,NY)
C                                                   -=-=- grnull
C                                               LOADS ARRAY WITH BLANK SPACES
      CHARACTER*1  A, B
      COMMON / GRAPH / A(130,130)
      DATA B/' '/
      DO 10 I = 1,NX
         DO 10 J = 1,NY
            A (I, J) = B
10    CONTINUE
      RETURN
C                       ****************************
      END
C
      SUBROUTINE GRPRNT(NX,NY)
C                                                   -=-=- grprnt

      CHARACTER*1 A, AA, BB(0:9), CC
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / GRAPH / A(130,130)
      DATA AA,CC /'-','|'/
      DATA BB /'0','1','2','3','4','5','6','7','8','9' /
C
      WRITE(NOUT,99) AA,(AA, I=1,NX),AA
      DO 30 I = 1,NY
         IPRIME = MOD(NY+1 - I,10)
         WRITE(NOUT,99) BB(IPRIME), (A(J,I), J = 1,NX), CC
30    CONTINUE
      WRITE(NOUT,99) AA,(BB(MOD(J,10)), J = 1,NX), AA
99    FORMAT(' ', 130A1)
      RETURN
C                       ----------------------------
C
      ENTRY GRFILE (NNX, NNY, NUNIT)
C
C      OPEN (NUNIT, STATUS='OLD')
      WRITE (NUNIT, 9) AA,(AA, I=1,NNX),AA
C
      DO 31 I = 1, NNY
         IPRIME = MOD( NNY+1 - I,10)
         WRITE (NUNIT, 9) BB(IPRIME), (A(J,I), J = 1,NNX), CC
   31 CONTINUE
C
      WRITE (NUNIT, 9) AA,(BB(MOD(J,10)), J = 1,NNX), AA
   9  FORMAT(' ', 130A1)
C
      RETURN
C                       ****************************
      END

      SUBROUTINE GSER(GAMSER,A,X,GLN)
C                                                   -=-=- gser
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
C                        ****************************
      END

      SUBROUTINE INPT2I (V1, V2, PROMPT, VN1, VX1, VN2, VX2, IR)
C                                                   -=-=- inpt2i
C
C                       Given 'Prompt', reads 2 Integer input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V2.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT INTEGER (T, V)
      CHARACTER*(*) PROMPT, REPLY*4
      LOGICAL TEST
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(2), TN(2), TX(2)
C
      TEST = .TRUE.
      IF (IR .EQ. 99) TEST = .FALSE.
      IR = 0

      T(1) = V1
      T(2) = V2
      TN(1) = VN1
      TN(2) = VN2
      TX(1) = VX1
      TX(2) = VX2
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/1X, A, 2I10)', PROMPT (1:LEN),
     >' or type / to keep the defaults :   ', (T(I), I=1,2)
C
      READ (NIN, *, ERR=99) (T(I), I=1,2)
C
      IF (.NOT.TEST) GOTO 8
      DO 10 I = 1, 2
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, I10, A, I10, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C               *************************
      END
C
      SUBROUTINE INPT2R (V1, V2, PROMPT, VN1, VX1, VN2, VX2, IR)
C                                                   -=-=- inpt2r
C
C                       Given 'Prompt', reads 2 REAL*8 input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V2.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(2), TN(2), TX(2)
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      T(1) = V1
      T(2) = V2
      TN(1) = VN1
      TN(2) = VN2
      TX(1) = VX1
      TX(2) = VX2
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/A, 2(1pE15.3))', PROMPT (1:LEN),
     >      ' or type / to keep the defaults :   ', (T(I), I=1,2)
C
      READ (NIN, *, ERR=99) (T(I), I=1,2)
C
      IF (ITST .EQ. 0) GOTO 8
      DO 10 I = 1, 2
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, 1pE10.3, A, E10.3, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C                       ****************************
      END
C
      SUBROUTINE INPT3I (V1, V2, V3, PROMPT,
C                                                   -=-=- inpt3i
     >                   VN1, VX1, VN2, VX2, VN3, VX3, IR)
C
C                       Given 'Prompt', reads 3 Integer input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V3.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT INTEGER (T, V)
      CHARACTER*(*) PROMPT, REPLY*4
      LOGICAL TEST
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(3), TN(3), TX(3)
C
      TEST = .TRUE.
      IF (IR .EQ. 99) TEST = .FALSE.
      IR = 0

      IR = 0
      T(1) = V1
      T(2) = V2
      T(3) = V3
      TN(1) = VN1
      TN(2) = VN2
      TN(3) = VN3
      TX(1) = VX1
      TX(2) = VX2
      TX(3) = VX3
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/1X, A, 3I10)', PROMPT (1:LEN),
     >' or type / to keep the defaults :   ', (T(I), I=1,3)
C
      READ (NIN, *, ERR=99) (T(I), I=1,3)
C
      IF (.NOT.TEST) GOTO 8
      DO 10 I = 1, 3
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, I10, A, I10, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
      V3 = T(3)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C               *************************
      END
C
      SUBROUTINE INPT3R (V1, V2, V3, PROMPT,
C                                                   -=-=- inpt3r
     >VN1, VX1, VN2, VX2, VN3, VX3, IR)
C
C                       Given 'Prompt', reads 3 REAL*8 input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V3.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(3), TN(3), TX(3)
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      T(1) = V1
      T(2) = V2
      T(3) = V3
      TN(1) = VN1
      TN(2) = VN2
      TN(3) = VN3
      TX(1) = VX1
      TX(2) = VX2
      TX(3) = VX3
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/ A, 3(1pE15.3))', PROMPT (1:LEN),
     >      ' or type / to keep the defaults : ', (T(I), I=1,3)
C
      READ (NIN, *, ERR=99) (T(I), I=1,3)
C
      DO 10 I = 1, 3
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, 1pE10.3, A, E10.3, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
      V3 = T(3)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C                       ****************************
      END
C
      SUBROUTINE INPT5R (V1, V2, V3, V4, V5, PROMPT,
C                                                   -=-=- inpt5r
     >VN1, VX1, VN2, VX2, VN3, VX3, VN4, VX4, VN5, VX5, IR)
C
C                       Given 'Prompt', reads 5 REAL*8 input numbers from
C                       terminal as the values for the REAL dummy variables
C                       V1 - V5.   The values are checked against
C                       the allowed ranges (VnI, VxI) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION T(5), TN(5), TX(5)
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      T(1) = V1
      T(2) = V2
      T(3) = V3
      T(4) = V4
      T(5) = V5
      TN(1) = VN1
      TN(2) = VN2
      TN(3) = VN3
      TN(4) = VN4
      TN(5) = VN5
      TX(1) = VX1
      TX(2) = VX2
      TX(3) = VX3
      TX(4) = VX4
      TX(5) = VX5
      CALL RTB (PROMPT, LEN)
C
    1 KPRT = 0
      PRINT '(/1X, A/ A / 5(1pE15.3))', PROMPT (1:LEN),
     >      ' or type / to keep the defaults : ', (T(I), I=1,5)
C
      READ (NIN, *, ERR=99) (T(I), I=1,5)
C
      IF (ITST .EQ. 0) GOTO 8
      DO 10 I = 1, 5
      IF (T(I) .LT. TN(I) .OR. T(I) .GT. TX(I)) THEN
         PRINT '(1X, A, I2, A, 1pE10.3, A, E10.3, a)',
     >  'Input #', I, ' outside the range [', TN(I), ', ', TX(I), ']'
         KPRT = 1
      ENDIF
   10 CONTINUE
C
      IF (KPRT .EQ. 1) THEN
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 99) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 V1 = T(1)
      V2 = T(2)
      V3 = T(3)
      V4 = T(4)
      V5 = T(5)
C
      RETURN
C
   99 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C                       ****************************
      END
 
C                                                          =-=-= Miscutl
      SUBROUTINE INPUTC (NAME, PROMPT, IR)
C                                                   -=-=- inputc
C
C===========================================================================
C GroupName: Inputs
C Description:  A variety of input prompt routines
C ListOfFiles: inputc inputi inputr inpt2i inpt2r inpt3i inpt3r inpt5r
C===========================================================================
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Inputs.f,v 1.1 97/12/21 21:19:16 wkt Exp $
C #Log:	Inputs.f,v $
c Revision 1.1  97/12/21  21:19:16  wkt
c Initial revision
c 

C                       Given 'Prompt', reads character input from
C                       terminal as the value for the char variable
C                       'Name'.  A return of '/' leave 'Name' unchanged
C
      CHARACTER*(*) NAME, PROMPT, TEM*80
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      IR = 0

      CALL RTB (PROMPT, LEN)
    1 PRINT '(/ 1X, A/ 2A)', PROMPT (1:LEN),
     >      ' or type / to keep the default :   ', NAME
C
      READ (NIN, '(A)', ERR=10) TEM
C
      IF (TEM(1:1) .NE. '/') NAME = TEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Error in INPUTC, Try again!'
      GOTO 1
C
C               *************************
      END

      SUBROUTINE INPUTI (IVALUE, PROMPT, IMIN, IMAX, IR)
C                                                   -=-=- inputi
C
C                       Given 'Prompt', reads INTEGER value input
C                       from terminal as the value for the dummy
C                       variable IVALUE.   The answer is checked against
C                       the allowed range (Imin, Imax) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
C                       If IR = 99 upon input, the limits are ignored.

      CHARACTER*(*) PROMPT, REPLY*1
      LOGICAL TEST
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      TEST = .TRUE.
      IF (IR .EQ. 99) TEST = .FALSE.
      IR = 0

      CALL RTB (PROMPT, LEN)

      ITEM = IVALUE
    1 PRINT '(/1X, A/A, I8)', PROMPT (1:LEN),
     >      ' or type / to keep the default :   ', ITEM
C
      READ (NIN, *, ERR=10) ITEM
C
      IF (.NOT.TEST) GOTO 8
      IF (ITEM .LT. IMIN .OR. ITEM .GT. IMAX) THEN
         PRINT *, 'Input outside the range [', IMIN, ', ', IMAX, ']'
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 10) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 IVALUE = ITEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C               *************************
      END
C
      SUBROUTINE INPUTR (VALUE, PROMPT, VMIN, VMAX, IR)
C                                                   -=-=- inputr
C
C                       Given 'Prompt', reads DOUBLE PRECISION input
C                       from terminal as the value for the dummy
C                       variable VALUE.   The answer is checked against
C                       the allowed range (Vmin, Vmax) first.
C                       An error code of Ir = 1 is returned if the input
C                       is outside the range.
C
C          An input value of IR = 99 will cause the limit test to be by-passed.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) PROMPT, REPLY*4
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      ITST = 1
      IF (IR .EQ. 99) ITST = 0
      IR = 0

      CALL RTB (PROMPT, LEN)
      TEM = VALUE
    1 PRINT '(/1X, A/A, 1PE15.4)', PROMPT (1:LEN),
     >      ' or type / to keep the default :   ', TEM
C
      READ (NIN, *, ERR=10) TEM
C
      IF (ITST .EQ. 0) GOTO 8
      IF (TEM .LT. VMIN .OR. TEM .GT. VMAX) THEN
         PRINT *, 'Input outside the range [', VMIN, ', ', VMAX, ']'
    2    PRINT *, 
     > 'Type ''G'' to try again, or type ''I'' to Ignore and go on.'
         READ (NIN, '(A)', ERR= 10) REPLY
         CALL UPCASE (REPLY)
         IF     (REPLY(1:1) .EQ. 'G') THEN
                GOTO 1
         ELSEIF (REPLY(1:1) .NE. 'I') THEN
                PRINT *, 'You must answer ''G'' or ''I''; try again!'
                GOTO 2
         ENDIF
         IR = 1
      ENDIF
C
    8 VALUE = TEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Data type Error, Try again!'
      GOTO 1
C
C                       ****************************
      END
C
      SUBROUTINE INQUIR (PROMPT, REPLY, LEN)
C                                                   -=-=- inquir
      CHARACTER*(*) PROMPT, REPLY
C
      COMMON /IOUNIT/NIN, NOUT, NWRT
C
      WRITE (NOUT, 100) PROMPT
      READ (NIN, 110) REPLY
      CALL RTB (REPLY, LEN)
      RETURN
C
 100  FORMAT ('$', A, ': ')
 110  FORMAT (A)
C               *************************
      END
C
      SUBROUTINE INQURE (NAME, PROMPT, IR)
C                                                   -=-=- inqure
C
C                       Given 'Prompt', reads character input from
C                       terminal as the value for the char variable
C                       'Name'.  A return of '/' leave 'Name' unchanged
C
      CHARACTER*(*) NAME, PROMPT, TEM*80
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
    1 WRITE (NOUT, '(/ ''$'', 2A)') PROMPT, 
     >   ' or type / to keep the default :   '
C
      READ (NIN, '(A)', ERR=10) TEM
C
      IF (TEM(1:1) .NE. '/') NAME = TEM
C
      RETURN
C
   10 WRITE (NOUT, *) 'Error in INQURE, Try again!'
      GOTO 1
C
C               *************************
      END

      Subroutine invmatrix(mx,n,np,imx)
      Integer n,np
      Double Precision mx(np,np),imx(np,np)
      Integer i,j,indx(n)
      Double Precision d,a(np,np)

      do i=1,n
         do j=1,n
            imx(i,j)=0D0
            a(i,j)=mx(i,j)
         enddo
         imx(i,i)=1D0
      enddo
      call ludcmp(a,n,np,indx,d)
      do j=1,n
         call lubksb(a,n,np,indx,imx(1,j))
      enddo
      return
      end

      Subroutine ludcmp(a,n,np,indx,d)
      Integer n,np,indx(n)
      Double Precision d,a(np,np),TINY
      Parameter (TINY=1.0D-20)
      Integer i,imax,j,k
      Double Precision aamax,dum,sum,vv(n)
      Double Precision d0,d1
      Parameter (d0=0D0, d1=1D0)

      d=d1
      Do i=1,n
         aamax=d0
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo
         if (aamax.eq.d0) pause 'singular matrix in ludcmp'
         vv(i)=d1/aamax
      enddo
      
      do j=1,n
         do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
          enddo
          aamax=d0
          do i=j,n
             sum=a(i,j)
             do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
             enddo
             a(i,j)=sum
             dum=vv(i)*abs(sum)
             if (dum.ge.aamax) then
                imax=i
                aamax=dum
             endif
          enddo
          if (j.ne.imax) then
             do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
             enddo
             d=-d
             vv(imax)=vv(j)
          endif
          indx(j)=imax
          if(a(j,j).eq.d0) a(j,j)=TINY
          if(j.ne.n) then
             dum=d1/a(j,j)
             do i=j+1,n
                a(i,j)=a(i,j)*dum
             enddo
          endif
       enddo
       return
       end

      Subroutine lubksb(a,n,np,indx,b)
      Integer n,np,indx(n)
      Double Precision a(np,np),b(n)
      Integer i,ii,j,ll
      Double Precision sum

      ii=0
      do i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0) then
            do j=ii,i-1
               sum=sum-a(i,j)*b(j)
            enddo
         elseif(sum.ne.0d0) then
            ii=i
         endif
         b(i)=sum
      enddo
      do i=n,1,-1
         sum=b(i)
         do j=i+1,n
            sum=sum-a(i,j)*b(j)
         enddo
         b(i)=sum/a(i,i)
      enddo
      return
      end

      SUBROUTINE LINFIT (X,Y,NDATA,SIG,MWT, A,B,SIGA,SIGB,CHI2,Q)
C                                                   -=-=- linfit

C                                           Routine FIT from Numerical Recipes
C     Linear fit to Y(i) = A * X(i) + B

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION X(NDATA), Y(NDATA), SIG(NDATA)

      SX=0.
      SY=0.
      ST2=0.
      B=0.
      IF(MWT.NE.0) THEN
        SS=0.
        DO 11 I=1,NDATA
          WT=1./(SIG(I)**2)
          SS=SS+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
      ELSE
        DO 12 I=1,NDATA
          SX=SX+X(I)
          SY=SY+Y(I)
12      CONTINUE
        SS=FLOAT(NDATA)
      ENDIF
      SXOSS=SX/SS
      IF(MWT.NE.0) THEN
        DO 13 I=1,NDATA
          T=(X(I)-SXOSS)/SIG(I)
          ST2=ST2+T*T
          B=B+T*Y(I)/SIG(I)
13      CONTINUE
      ELSE
        DO 14 I=1,NDATA
          T=X(I)-SXOSS
          ST2=ST2+T*T
          B=B+T*Y(I)
14      CONTINUE
      ENDIF
      B=B/ST2
      A=(SY-SX*B)/SS
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
      SIGB=SQRT(1./ST2)
      CHI2=0.
      IF(MWT.EQ.0) THEN
        DO 15 I=1,NDATA
          CHI2=CHI2+(Y(I)-A-B*X(I))**2
15      CONTINUE
        Q=1.
        SIGDAT=SQRT(CHI2/(NDATA-2))
        SIGA=SIGA*SIGDAT
        SIGB=SIGB*SIGDAT
      ELSE
        DO 16 I=1,NDATA
          CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2
16      CONTINUE
        AAA = 0.5*(NDATA-2)
        BBB = 0.5*CHI2
        Q=GAMMQ(AAA, BBB)
      ENDIF
      RETURN
C               ******************
      END

      SUBROUTINE LOGGAM(X, U)
C                                                   -=-=- loggam
C
C       SUBROUTINE LOGGAM--TRANSCRIBED FROM NYU FAP ROUTINE OF MAX
C       GOLDSTEIN WRITTEN FOR FORTRAN IV
C       CALL IS CALL LOGGAM(X, U) WHERE X IS REAL
C        AND U IS THE FUNCTION VALUE
C        MUTILATION TO REAL VALUES COURTESY OF PORTER JOHNSON JUNE 1982
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION H(7)
C
      DATA H/ 2.69488974, 1.517473649, 1.011523068, 0.525606469,
     1          0.2523809524, 0.033333333, 0.0833333333 /
C
      B2 = 0.D0
      J = 2
      X2 = X
3797    IF(X) 2794, 9999, 100
100   T = X**2
5793    B7 =  T
C
C       REAL PART OF LOG
C
      T1 = 5D-1 * LOG(B7)
      IF(X .GE. 2.D0) GO TO 3793
      B2 = B2 + T1
      X  = X  + 1.D0
      J  = 1
      GOTO 3797
C
3793    T3 = T1* (X - 5D-1) -X + 0.9189385332
      T4 = X
      T1 = B7
      DO 200 I = 1,7
      T = H(I)/T1
      T4 = T*T4 + X
200     T1 = T4**2
C
      T3 = T4 -X + T3
      GO TO (8795,4794) , J
C
8795    T3 = T3 - B2
C
4794    IF(X2) 4796, 4795, 4795
4795    U = T3
      X = X2
      RETURN
C
4796    U = T3 -E4
      X = X2
      RETURN
C
C               X IS NEGATIVE
C
2794    E4 = 0.D0
      IE6 = 0
5797    E4 = E4 + LOG(ABS(X))
      IE6 = IE6 + 1
      X = X + 1.D0
      IF( X .LT. 0.D0) GO TO 5797
      GO TO 3797
C
9999    WRITE(1,9990) X2
9990    FORMAT(' ATTEMPTED TO TAKE LOGGAMOF X = ',F12.5)
C
      STOP
C                       *********************************
      END
C
      FUNCTION NEXTUN()
C                                                   -=-=- nextun
C                                    Returns an unallocated FORTRAN i/o unit.
      LOGICAL EX
C
      DO 1  N = 10, 98
          INQUIRE (UNIT=N, OPENED=EX)
          IF (.NOT. EX) then
             nextun = n
             RETURN
          end if
   1  CONTINUE

      Stop 'NextUnit number not found!  Stopped in NextUn().'

C               *************************
      END
C
C                                                          =-=-= Charutl
      FUNCTION PI()
C                                                   -=-=- pi
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (P=3.1415927)
      ENTRY PIFN()
      PI = P
      RETURN
      END
C                       ****************************
C
      SUBROUTINE PLT4D (NCUR, NUMPT, XPT, YPT, NAMEX, NAMEY, NUU)
C                                                   -=-=- plt4d
C
C===========================================================================
C GroupName: Plot4
C Description:   Poor-man's Graphs
C ListOfFiles: plt4d pltchd pltind grdatd grnull grprnt sortd sorti
C===========================================================================
C                                          By Porter Johnson and Wu-Ki Tung
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Plot4.f,v 1.1 97/12/21 21:19:25 wkt Exp $
C #Log:	Plot4.f,v $
c Revision 1.1  97/12/21  21:19:25  wkt
c Initial revision
c 

C       03 10 87
     
C                       Main Switching routine for plotting
C                       Determines ranges of variables
C                       Sorts data points in descending order in y
C                       Performs re-scaling if  desired
C                       Switches to selected plotting device
C                       Returns to calling program
C
C                       Mxxpt is the maximum limit of Ncur * Numpt
     
C       Input parameters:
     
C       Ncur    :       Number of curves = Number of columns of Y-array
     
C       Numpt   :       Number of x-points
     
C       Xpt     :       Array of x - values   [ Xpt (1 : Numpt) ]
     
C       Ypt     :       2-dim array of y - values    [Numpt x Ncur]
     
C       NameX   :       character const. = Name of x - variable
     
C       NameY   :       character array  = Names of y - variables
     
C       Nuu     :       Unit number to which records of plots and
C                       tables are to be written (on demand).
C
C                       Nuu = 0  means that no output file exists
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXCUR = 20, MXXPT = 2500)
     
      CHARACTER*130 COMENT, HEADNG
      CHARACTER*2 AL*4
      CHARACTER*(*) NAMEX, NAMEY(NCUR), NMY*20
      CHARACTER*11 NAMX, NAMY(MXCUR), NAM
      LOGICAL ASK, LDPD, ZX, ZY
     
      DIMENSION SCALE(MXCUR)
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
      DIMENSION SX (MXXPT), SY (MXXPT)
     
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY
     
      DATA XSML, YSML / 1E-5, 1E-5 /
      DATA SML /  1E-5 /
C                                Initiations
      NU = NUU
      AL = 'None'
      IACT = 0
     
      DO 19 I = 1, MXCUR
19       SCALE(I)= 1.
C                                Process the labels of fit the heading format
      DO 50 J = 0, NCUR
         IF (J .NE. 0) THEN
            NMY = NAMEY(J)
         ELSE
            NMY = NAMEX
         ENDIF
     
         CALL TRMSTR (NMY, KNM)
         KPD = (11-KNM)/2
         IF (KPD .LE. 0) THEN
            NAM = NMY
         ELSE
            DO 55 I = 1, KPD+1
               NAM(I:I) = ' '
   55       CONTINUE
            NAM(KPD+1:) = NMY
         ENDIF
     
         IF (J .EQ. 0) THEN
            NAMX = NAM
         ELSE
            NAMY(J) = NAM
         ENDIF
   50 CONTINUE
C                    Make sure that there is enough room for data points
     
      IF ( (NCUR * NUMPT)  .GT. MXXPT) THEN
         PRINT *,
     >   'Input data points exceed maximum array size in PLOTT routine'
         STOP
      ENDIF
C                       SORTD does two main tasks:
C                       It computes the limits of x- and y- coordinates
C                       and gives the results in  RANGED Common.
C                       It sorts the Y- values into descending order,
C                       and outputs the results as ONE-DIM arrays.
     
    1 CALL SORTD (NCUR, NUMPT, XPT, YPT, SX, SY, IACT)
      PRINT *, 'The limits for the X,Y-axes and the curves are:'
     
      WRITE (NOUT, 931) XMAX, XMIN, YMAX, YMIN
  931 FORMAT (/' Xmax = ', T10, 1PE10.2, 10X, 'Xmin = ', T40, E10.2
     >        /' Ymax = ', T10, 1PE10.2, 10X, 'Ymin = ', T40, E10.2 )
      WRITE (NOUT, *) '------------------------------------------------'
C
      WRITE (NOUT, 933)
     >    ( J, SY(1 + (J-1)*NUMPT), J, SY(NUMPT*J), J = 1, NCUR )
  933 FORMAT(' Y',I2,'mx = ' , T10, 1PE10.2, 10X,
     >          'Y',I2,'mn = ' , T40, 1PE10.2 )
      WRITE (NOUT, *)
     
      IF (ASK('Do you wish to re-scale any of the curves')) THEN
         PRINT *, 'Enter the scaling factors, in order: ',
     >                   'Terminate with / <cr>'
         READ *, SCALE
         DO 71 JCUR = 1, NCUR
            DO 81 IX = 1, NUMPT
                YPT(IX, JCUR) = YPT(IX, JCUR) * SCALE(JCUR)
                SY (IX + (JCUR-1)*NUMPT) =
     >                  SY (IX + (JCUR-1)*NUMPT) * SCALE(JCUR)
81           CONTINUE
71       CONTINUE
         IACT = 1
         GOTO 1
      ENDIF
     
      YX = YMAX
      YN = YMIN
701   WRITE (NOUT, *) 'Choose the Y-range for the graph; defaults ',
     >  'are:  Ymax, Ymin  ='
      WRITE (NOUT, *) YMAX, YMIN
     
      READ (NIN, *, ERR=701) YX, YN
      IF (YN .GE. YX) THEN
         WRITE (NOUT, *) 'Ymin > Ymax not allowed, try again!'
         GOTO 701
      ENDIF
C      If (Yn .Lt. Ymin-Abs(Ymin*Sml) .or. Yx .Gt. Ymax+Abs(Ymax*Sml))
C     >   Then
C         Write (Nout, *) 'New range cannot exceed original range;'
C         Write (Nout, *) 'Defaults are:', Ymax, Ymin, '   Try again!'
C         Goto 701
C      Endif
     
      DO 72 JC = 1, NCUR
         DO 82 IX = 1, NUMPT
            YPT (IX, JC) = MAX (YPT(IX, JC), YN)
            YPT (IX, JC) = MIN (YPT(IX, JC), YX)
            ST = SY (IX + (JC - 1) * NUMPT)
            ST = MAX (ST, YN)
            ST = MIN (ST, YX)
            SY (IX + (JC - 1) * NUMPT) = ST
  82     CONTINUE
  72  CONTINUE
     
      YMIN = YN
      YMAX = YX
     
   7  CALL PLTCHD (NCUR, NUMPT, XPT, YPT, NAMEX, NAMEY, SCALE, AL, NU)
     
C                               Restore the scale factors and Log flag
      DO 29 I = 1, MXCUR
   29    SCALE(I)= 1.
     
      IF (ASK('Do you wish to print out the results in tabular form'))
     >   THEN
         NCR = MIN (NCUR, 6)
         WRITE (NOUT, 1003) NAMX, (NAMY(I), I = 1, NCR)
 1003    FORMAT(/ 3X, 6(1X, A11, 1X) /)
     
         DO 543 IPT = 1, NUMPT
            WRITE (NOUT, 1004) XPT(IPT), (YPT(IPT, ICUR), ICUR=1, NCR)
  543    CONTINUE
 1004    FORMAT ( 1PE13.4, 6E13.4 )
     
         IF (NCUR .GT. 6) THEN
            WRITE (NOUT, 1005) NAMX, (NAMY(I), I = 6, NCUR)
 1005       FORMAT(// 2X, 6(1X, A11, 1X) /)
     
            DO 544 IPT = 1, NUMPT
               WRITE (NOUT, 1006) XPT(IPT), (YPT(IPT,ICUR),ICUR=6,NCUR)
 1006          FORMAT( 1PE13.4, 6E13.4 )
  544       CONTINUE
         ENDIF
         IF (NU .NE. 0) THEN
            LDPD=ASK('Do you wish to dump this table to the extl file')
            IF (LDPD) THEN
  401          WRITE (NOUT, *) 'Input Heading for the table: '
               WRITE (NOUT, *)
               READ (NIN, '(A)', IOSTAT=IOS, ERR=501) HEADNG
  501          IF (IOS .GT. 0) THEN
                  PRINT *, 'Read Error, IOSTAT = ', IOS, ' Try again!'
                  GOTO 401
               ENDIF
               WRITE (NU, 1013) HEADNG, NAMX, (NAMY(I), I = 1, NCR)
 1013          FORMAT(/ 1X, A / 2X, 6(1X, A11, 1X) /)
     
               DO 556 IPT = 1, NUMPT
                  WRITE (NU, 1014) XPT(IPT), (YPT(IPT,ICUR), ICUR=1,NCR)
 1014             FORMAT( 1PE13.4, 6E13.4 )
  556          CONTINUE
     
               IF (NCUR .GT. 6) THEN
                  WRITE (NU, 1015) NAMX, (NAMY(I), I = 6, NCUR)
 1015             FORMAT(/ 2X, 6(1X, A11, 1X) /)
C
                  DO 557 IPT = 1, NUMPT
                     WRITE (NU,1016)
     >                  XPT(IPT),(YPT(IPT,ICUR),ICUR=6,NCUR)
 1016                FORMAT( 1PE13.4, 6E13.4 )
  557             CONTINUE
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C                                       Initiation for the next section
      AL = 'None'
      ZX = .TRUE.
      ZY = .TRUE.
     
      IF ( ASK('Do you wish to convert either axis to Log-scale') )
     >          THEN
   46    WRITE (NOUT, 941)
941      FORMAT(' Type   LY   for log-scale y-axis,' /
     >                 '   "    LX   for log-scale x-axis,' /
     >                 '   "    LL   for log-log plot.    ')
         READ  (NIN, '(A)')  AL
         CALL UPCASE (AL)
     
C           XSML = MAX(XMIN, XSML)
C           YSML = MAX(YSML, YMIN)
         WRITE  (NOUT, 1101) XMIN, YMIN
 1101    FORMAT (
     >     ' If x or y becomes negative, taking log is illegal;' /
     >     ' Alternatively, if Xmin or Ymin is very close to zero, ',
     >         'the logarithm'/
     >    ' may extend too far into the negative range than desirable.'
     >    / ' Current values of XMIN and YMIN are: ',2(1PE15.2))
         Ymin = 1D-20
         CALL INPT2R (XSML, YSML, 'Specify cutoffs for X & Y',
     >          SML*SML, XMAX, Ymin, YMAX, J9)
     
C     >   /' Specify desired cutoff for  X  and  Y ;    Defaults are: '/
C     >   2(1PE15.2) )
C         READ (NIN, *) XSML, YSML
     
         IF (AL .EQ. 'LX' .OR. AL .EQ. 'LL') THEN
            IF (XMIN .LT. 0) PRINT *, 'Negative argument in Log X.'
            DO 61 IX = 1, NUMPT
               XX = XPT(IX)
               IF (XX .LT. XSML) THEN
                  IF (ZX) THEN
                     PRINT *, 'X < XSMALL; X set = ', XSML
                     ZX = .FALSE.
                  ENDIF
                  XX = XSML
               ENDIF
               XPT (IX) = LOG10 (XX)
61          CONTINUE
            IACT = 1
         ENDIF
         IF (AL .EQ. 'LY' .OR. AL .EQ. 'LL') THEN
94          IF (YMIN .LE. 0)  THEN
               PRINT *,
     >          'Negative argument in Log Y, converting to Log AbsY'
            ENDIF
            DO 602 JC = 1, NCUR
               DO 62 IX = 1, NUMPT
                  YY = ABS( YPT(IX,JC) )
                  YS = ABS( SY(IX+(JC-1)*NUMPT) )
                  IF (YY .LT. YSML) THEN
                     IF (ZY) THEN
                        PRINT *, 'Y < YSMALL; set Y = ', YSML
                        ZY = .FALSE.
                     ENDIF
                     YY = YSML
                  ENDIF
                  IF (YS .LT. YSML) YS = YSML
                  YPT(IX,JC) = LOG10 (YY)
                  SY(IX+(JC-1)*NUMPT) = LOG10 (YS)
62             CONTINUE
602         CONTINUE
            IACT = 1
         ENDIF
     
         IF (IACT .EQ. 0) THEN
            PRINT *, 'Illegal Choice of LX, LY and LL, Try again!'
            GOTO 46
         ELSE
            GOTO 1
         ENDIF
      ENDIF
     
  705 IF (ASK('Add a comment to the above results for the record'))THEN
         WRITE (NOUT, *)
         READ (NIN, '(A)', ERR=707) COMENT
         WRITE (NU, *) COMENT
         WRITE (NOUT, *)
      ENDIF
     
      RETURN
     
707   WRITE (NOUT, *) 'Error reading input; Try again!'
      GOTO 705
C                       ****************************
      END
     
      SUBROUTINE PLTCHD
C                                                   -=-=- pltchd
     >(NCUR, NUMPT, XPT, YPT, NAMEX, NAMEY, SCALE, AL, NUU)
C
C                               Produces Plots of Ncur curves in ASC charaters
C
C               Input parameters:
C
C               Ncur, Numpt, Xpt, Ypt, NameX, NameY, Nuu:    Same as in PLOT2
C
C               Scale   :       Number array of Scale factor for the cureves
C
C               Al      :       Flag for Logrithmic X- and/or Y-axis
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER*(*) AL
      CHARACTER*(*) NAMEX, NAMEY(NCUR)
      CHARACTER HEADNG*130
      LOGICAL LDP, ASK
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY
C
      DIMENSION SCALE(10)
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
C
      DATA NX, NY / 76, 20 /
C
      NU = NUU
C
      WRITE (NOUT, 99)
99    FORMAT(/ ' The following functions are being plotted' / )
      DO 30 I = 1, NCUR
         IF (I .EQ. 10) THEN
            J = 0
         ELSE
            J = I
         ENDIF
         WRITE (NOUT, 98)  NAMEY(I), J
98       FORMAT(' The function ', A,' --- symbol is ', I1)
30    CONTINUE
C
40    WRITE(NOUT,999) NUMPT
999   FORMAT( 38X, ' Number of data points is = ', I4, //
     >  ' Put in Nx (x-mesh), Ny (y-mesh), ',
     >  'or  Type  /  to use current default values.' / )
      PRINT *, 'Default:  Nx =', NX,  '         Ny =', NY
      READ (NIN,*) NX, NY
      IF ( (NX.LT.10) .OR. (NX.GT.130) ) GOTO 101
C
      DX = (XMAX - XMIN) / FLOAT(NX-1)
      DY = (YMAX - YMIN) / FLOAT(NY-1)
      IF (DY.LE.0) GOTO 50
C
      CALL GRNULL (NX, NY)
      CALL GRDATD (NCUR, NX, NY, NUMPT, XPT, YPT)
      CALL GRPRNT (NX, NY)
C
      WRITE(NOUT,997) XMIN, NAMEX, XMAX,   YMIN, YMAX
997     FORMAT( 1PE12.3, '<= X (', A, ') <= ', 1PE12.3,
     >          T50, 1PE12.3, '<= Y <= ', 1PE12.3 )
C
      IF (NU .NE. 0) THEN
C
         LDP = ASK('Do you wish to dump this graph to the extl file')
         IF (LDP) THEN
  401       WRITE (NOUT, *) 'Input Heading for the Plot:'
            READ  (NIN, '(A)', ERR=401) HEADNG
C
            WRITE (NU, *) HEADNG
            CALL GRFILE (NX, NY, NU)
            WRITE (NU, 977) XMIN, NAMEX, XMAX,   YMIN, YMAX
977         FORMAT( 1PE12.3, '<= X (', A, ') <= ', 1PE12.3,
     >          T49, 1PE12.3, '<= Y <= ', 1PE12.3 )
            IF (AL .EQ. 'LX' .OR. AL .EQ. 'LL')
     >          WRITE (NU, *) 'X-axis in Ln scale'
            IF (AL .EQ. 'LY' .OR. AL .EQ. 'LL') THEN
               DO 97 I = 1, NCUR
                  NAMEY(I) = 'Ln '//NAMEY(I)
97             CONTINUE
            ENDIF
            WRITE (NU, '(/)')
            WRITE (NU, 198) (I, NAMEY(I), SCALE(I), I = 1, NCUR)
198         FORMAT ( 1X, I2, ' :  ', A, '*', F8.1 )
         ENDIF
      ENDIF
C
50    CONTINUE
C
      IF (ASK (
     >    'Do you wish to change the meshsize for the same plot'))
     >  GOTO 40
C
      RETURN
C
101   WRITE (NOUT, 901)
901   FORMAT ( ' Input mesh-size outside the acceptable range!'//
     >         ' 10 < Nx < 130 .    Try Again!'/)
      GOTO 40
C                       ****************************
      END
C
      SUBROUTINE PLTIND (MXCUR, MXPT, NCUR, NUMPT, YP, YPT)
C                                                   -=-=- pltind
C               Converts the two-dim array Yp (Mxpt, Mxcur) with actual non-
C               zero dimension (Numpt, Ncur) into an one-dim array Ypt with
C               no unwanted zeros in the middle.   Helps avoid confusion and
C               mis-handling of data in the plotting routines which follow.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION YP(MXPT, MXCUR), YPT(MXPT*MXCUR)
      I = 0
      DO 10 ICUR = 1, NCUR
         DO 20 IX   = 1, NUMPT
           I = I + 1
           YPT(I) = YP(IX, ICUR)
20       CONTINUE
10    CONTINUE
      RETURN
C                       ****************************
      END
C
      SUBROUTINE POLIN2(X1A,X2A,YA,M,N,X1,X2,Y,DY)
C                                                   -=-=- polin2
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=20,MMAX=20)
      DIMENSION X1A(M),X2A(N),YA(M,N),YNTMP(NMAX),YMTMP(MMAX)
 
      DO 12 J=1,M
        DO 11 K=1,N
          YNTMP(K)=YA(J,K)
11      CONTINUE
        CALL POLINT(X2A,YNTMP,N,X2,YMTMP(J),DY)
12    CONTINUE
      CALL POLINT(X1A,YMTMP,M,X1,Y,DY)
      RETURN
      END
C                        ****************************
 
      SUBROUTINE POLINT3 (XA,YA,N,X,Y,DY)
c fast version of polint, valid only for N=3
c Have explicitly unrolled the loops.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      PARAMETER (NMAX=3)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

	if(n .ne. 3) then
	   print *,'fatal polint3 call',n
	   stop
	endif

      NS=1
      DIF=ABS(X-XA(1))

        DIFT=ABS(X-XA(1))
        IF (DIFT.LT.DIF) THEN
          NS=1
          DIF=DIFT
        ENDIF
        C(1)=YA(1)
        D(1)=YA(1)

        DIFT=ABS(X-XA(2))
        IF (DIFT.LT.DIF) THEN
          NS=2
          DIF=DIFT
        ENDIF
        C(2)=YA(2)
        D(2)=YA(2)

        DIFT=ABS(X-XA(3))
        IF (DIFT.LT.DIF) THEN
          NS=3
          DIF=DIFT
        ENDIF
        C(3)=YA(3)
        D(3)=YA(3)


      Y=YA(NS)
      NS=NS-1


          HO=XA(1)-X
          HP=XA(2)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN


          HO=XA(2)-X
          HP=XA(3)-X
          W=C(3)-D(2)
          DEN=W/(HO-HP)
          D(2)=HP*DEN
          C(2)=HO*DEN


        IF (2*NS.LT.2)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY


          HO=XA(1)-X
          HP=XA(3)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN

        IF (2*NS.LT.1)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY

      RETURN
      END
      SUBROUTINE POLINT4 (XA,YA,N,X,Y,DY)
c fast version of polint, valid only for N=4
c Have explicitly unrolled the loops.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      PARAMETER (NMAX=4)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

	if(n .ne. 4) then
	   print *,'fatal polint4 call',n
	   stop
	endif

      NS=1
      DIF=ABS(X-XA(1))

        DIFT=ABS(X-XA(1))
        IF (DIFT.LT.DIF) THEN
          NS=1
          DIF=DIFT
        ENDIF
        C(1)=YA(1)
        D(1)=YA(1)

        DIFT=ABS(X-XA(2))
        IF (DIFT.LT.DIF) THEN
          NS=2
          DIF=DIFT
        ENDIF
        C(2)=YA(2)
        D(2)=YA(2)

        DIFT=ABS(X-XA(3))
        IF (DIFT.LT.DIF) THEN
          NS=3
          DIF=DIFT
        ENDIF
        C(3)=YA(3)
        D(3)=YA(3)

        DIFT=ABS(X-XA(4))
        IF (DIFT.LT.DIF) THEN
          NS=4
          DIF=DIFT
        ENDIF
        C(4)=YA(4)
        D(4)=YA(4)


      Y=YA(NS)
      NS=NS-1


          HO=XA(1)-X
          HP=XA(2)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN


          HO=XA(2)-X
          HP=XA(3)-X
          W=C(3)-D(2)
          DEN=W/(HO-HP)
          D(2)=HP*DEN
          C(2)=HO*DEN


          HO=XA(3)-X
          HP=XA(4)-X
          W=C(4)-D(3)
          DEN=W/(HO-HP)
          D(3)=HP*DEN
          C(3)=HO*DEN


        IF (2*NS.LT.3)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY



          HO=XA(1)-X
          HP=XA(3)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN


          HO=XA(2)-X
          HP=XA(4)-X
          W=C(3)-D(2)
          DEN=W/(HO-HP)
          D(2)=HP*DEN
          C(2)=HO*DEN



        IF (2*NS.LT.2)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY


          HO=XA(1)-X
          HP=XA(4)-X
          W=C(2)-D(1)
          DEN=W/(HO-HP)
          D(1)=HP*DEN
          C(1)=HO*DEN

        IF (2*NS.LT.4-3)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY


      RETURN

      END
      SUBROUTINE POLINT (XA,YA,N,X,Y,DY)
C                                                   -=-=- polint
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.) THEN
CCPY            Print*, "Error: DEN = 0 in PolInt", HO, HP, XA(I), X,
CCPY     >      XA(I+M+1), I, M
            PAUSE
          ENDIF
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
C                        ****************************
 
      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
C                                                   -=-=- ratint
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.0
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)
          IF(DD.EQ.0.)PAUSE
          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
C                        ****************************
      END

      SUBROUTINE RBL (NU1, NU2, ILIN)
C                                                   -=-=- rbl

C===========================================================================
C GroupName: Charutl
C Description:   A collection of routines to handle characters
C ListOfFiles: rbl rtb trmstr uc upc
C===========================================================================

C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Charutl.f,v 1.1 97/12/21 21:19:10 wkt Exp $
C #Log:	Charutl.f,v $
c Revision 1.1  97/12/21  21:19:10  wkt
c Initial revision
c 

C                 Remove Trailing blanks    = 0   and put one blank at the end
C                 Remove Blank Lines:  ILIN = 1   removes all blank lines;
C                                             2   removes multiple blank lines
C                 Truncate to 80 record len = 3
C                 Also removes trailing blanks from all lines.

C                 NU1 : unit # of input  file to be processed.
C                 NU2 : unit # of output file
      CHARACTER LINE*132

      IF (ILIN .EQ. 3) THEN
        JLIN = -1
      ELSE
        JLIN = ILIN
      ENDIF

      Rewind (NU1)
      Rewind (NU2)

      DO 3 I = 1, 10000

      DO 5 IBLK = 0, 500
        Read (NU1, '(A)', END=10) LINE
        CALL RTB (LINE, LEN)
        IF (JLIN .EQ. -1) LEN = MIN (LEN, 80)
        IF (LEN .NE. 0) GOTO 6
    5 CONTINUE

    6 IF (IBLK .GT. 0 .AND. JLIN .NE. 1) WRITE (NU2, *)

      IF (JLIN .EQ. 0) THEN
        WRITE (NU2, '(A)') LINE(1:LEN)//' '
      ELSE
        WRITE (NU2, '(A)') LINE(1:LEN)
      ENDIF

    3 CONTINUE

   10 CONTINUE

      RETURN
C                        ****************************
      END
      FUNCTION ROMINT(F, A, B, N, ESTER, IRET)
C                                                   -=-=- romint
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
C                                       Romberg Integration of function F from
C                                       A to B using 2 ** N points.
C
C                                       EstEr is the estimated error.
C                                       Iret  is a return code
C                                            0   ---         O.K.
C                                            1   ---   N   is too big
C                                            2   ---   Lower limit A is bigger
C                                                      then upper limit B
      PARAMETER (MXN = 10)
      external f
C
      DIMENSION R(2, MXN+1)
C
      IRET = 0
      IF (N .GT. MXN) THEN
         IRET = 1
         ROMINT = 0
         RETURN
      ENDIF
C
      IF     (A .GT. B) THEN
         IRET = 2
         ROMINT = 0
         RETURN
      ELSEIF (A .EQ. B) THEN
         ROMINT = 0
         RETURN
      ENDIF
C
      H = B - A
      R(1, 1) = (F(A) + F(B)) * H / 2.D0
C
      DO 10 I = 2, N+1
C
         RT = 0
         XT = A - H / 2.D0
C
         DO 15 K = 1, 2**(I-2)
            RT = RT + F(XT + H * K)
   15    CONTINUE
C
         R(2, 1) = (R(1, 1) + H * RT) / 2.D0
C
         DO 16 J = 2, I
            AJ = 4.D0 ** (J - 1)
            R(2, J) = (AJ * R(2, J-1) - R(1, J-1)) / (AJ - 1.D0)
   16    CONTINUE
C
         H = H / 2.D0
         DO 17 J = 1, I
            R(1, J) = R(2, J)
   17    CONTINUE
C
   10 CONTINUE
C
      ROMINT = R(1, N+1)
      ESTER  = R(1, N+1) - R(1, N)
C
      RETURN
C                        ****************************
      END
C
      SUBROUTINE RTB (ch, lench)
C                                                   -=-=- rtb
C                           Set LENCH = length of CH, not counting trailing
C                           blanks and nulls
      CHARACTER*(*) ch
      integer lenmax,lench
C                         slightly change excluding trailing Ctrl key
      lenmax=LEN(ch)
      do lench=lenmax,1,-1
        if(ch(lench:lench).gt.' ') return
      enddo
      lench=0
      return
C               *************************
      END

      Subroutine SetUTL
C                                                   -=-=- setutl

C===========================================================================
C GroupName: Setutl
C Description: setup module and block data
C ListOfFiles: setutl
C===========================================================================
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      External DatUtl

      Dummy = 0.

      Return
C                        ****************************
      END

      BLOCK DATA DATUTL

      COMMON / IOUNIT / NIN, NOUT, NWRT

      DATA NIN, NOUT, NWRT / 5, 6, 56 /

C                         *************************
      END

C                                                          =-=-= Adzint
      FUNCTION SIMP (NX, DX, F)
C                                                   -=-=- simp

C===========================================================================
C GroupName: Simpgaus
C Description:  many variants of Simpson Integration + an adaptive gausint
C ListOfFiles: simp smpnol smpnor smpsna smpsnf gausint
C===========================================================================
C
C #Header: /Net/d2a/wkt/1hep/1utl/RCS/Simpgaus.f,v 1.3 98/07/21 16:29:16 wkt Exp $
C #Log:	Simpgaus.f,v $
c Revision 1.3  98/07/21  16:29:16  wkt
c make gausint consistent
c 
c Revision 1.2  98/03/09  01:12:37  wkt
c Set function=0 for error returns in RomInt and Zbrnt.
c 
c Revision 1.1  97/12/21  21:19:32  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      DIMENSION F(NX)
C
      IF (NX .EQ. 1) Then
        simp = 0
        RETURN
      EndIf

      IF (NX .LT. 0 .OR. NX .GT. 10000) GOTO 99
      IF (NX .GT. 4) GOTO 50
      GOTO (20, 30, 40), NX-1
C
   20 SIMP = (F(1) + F(2)) * DX / 2D0
      RETURN
   30 SIMP = (F(1) + 4D0 * F(2) + F(3)) * DX / 3D0
      RETURN
   40 SIMP = (( F(1) + 4D0 * F(2) +     F(3)) / 3D0
     >       +(-F(2) + 8D0 * F(3) + 5D0 * F(4)) / 12D0 ) * DX
      RETURN
C
   50 SE = F(2)
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + F(I)
      SO = SO + F(IM1)
   60 CONTINUE
      MS = MOD (NX, 2)
C
      IF (MS .EQ. 1) THEN
         SIMP = (F(1) + 4D0 * SE + 2D0 * SO + F(NX)) * DX / 3D0
      ELSE
         SIMP = (F(1) + 4D0 * SE + 2D0 * SO + F(NM1)) * DX / 3D0
     >         + (-F(NM1-1) + 8D0 * F(NM1) + 5D0 * F(NX)) * DX / 12D0
      END IF
C
      RETURN
C
   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >              'out of range in D-SIMP INTEGRATION ROUTINE')
      STOP
C                        ****************************
      END
 
      FUNCTION SMPNOL (NX, DX, FN, ERR)
C                                                   -=-=- smpnol
C                                            DP Left-Open Simpson Integration
C     Inputs:  Nx, Dx, Fn(1:Nx)                             [F(1) is not used]
C     Output:  Err                                          (Error estimate)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)
 
      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN SMPNOL!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(2)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(2) * 2.
      ELSE
         IF (MS .EQ. 0) THEN
            TEM = DX * (23.* FN(2) - 16.* FN(3) + 5.* FN(4)) / 12.
            TMP = DX * (3.* FN(2) - FN(3)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + SMPSNA (NX-1, DX, FN(2), ER1)
            ERR = ABS(ER1) + ERR
         ELSE
            TEM = DX * (8.* FN(2) - 4.* FN(3) + 8.* FN(4)) / 3.
            TMP = DX * (3.* FN(2) + 2.* FN(3) + 3.* FN(4)) / 2.
            ERR = ABS(TEM - TMP)
            TEM = TEM + SMPSNA (NX-4, DX, FN(5), ER1)
            ERR = ABS(ER1) + ERR
         ENDIF
      ENDIF
 
      SMPNOL = TEM
      RETURN
C                        ****************************
      END
C
      FUNCTION SMPNOR (NX, DX, FN, ERR)
C                                                   -=-=- smpnor
C                                              DP Right-Open Simpson Integration
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FN(NX)
 
      MS = MOD(NX, 2)
      IF (NX .LE. 1 .OR. NX .GT. 1000) THEN
         PRINT *, 'NX =', NX, ' OUT OF RANGE IN SMPNOR!'
         STOP
      ELSEIF (NX .EQ. 2) THEN
         TEM = DX * FN(Nx-1)
      ELSEIF (NX .EQ. 3) THEN
         TEM = DX * FN(NX-1) * 2.
      ELSE
        IF (MS .EQ. 0) THEN
         TEM = DX * (23.* FN(NX-1) - 16.* FN(NX-2) + 5.* FN(NX-3)) / 12.
         TMP = DX * (3.* FN(NX-1) - FN(NX-2)) / 2.
         ERR = ABS(TEM - TMP)
         TEM = TEM + SMPSNA (NX-1, DX, FN(1), ER1)
         ERR = ER1 + ERR
        ELSE
         TEM = DX * (8.* FN(NX-1) - 4.* FN(NX-2) + 8.* FN(NX-3)) / 3.
         TMP = DX * (3.* FN(NX-1) + 2.* FN(NX-2) + 3.* FN(NX-3)) / 2.
         ERR = ABS(TEM - TMP)
         TEM = TEM + SMPSNA (NX-4, DX, FN(1), ER1)
         ERR = ER1 + ERR
        ENDIF
      ENDIF
 
      SMPNOR = TEM
      RETURN
C                        ****************************
      END
 
      FUNCTION SMPSN1(FN, A, B, NX, ERR)
c modified version of smpsnf (jcp 11/27/01)
C
C                       Does integral of FN(X)*dx from A TO B by SIMPSON'S METHOD
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / IOUNIT / NIN, NOUT, NWRT

      PARAMETER (MXPT = 10000)

      DIMENSION X(MXPT)
      external fn

      IF ((NX .LT. 0) .OR. (NX .GT. MXPT)) GOTO 99
C
      DX = (B - A) / FLOAT(NX-1)

c print warning if A > B; but routine SHOULD be ok either way.
c original routine smpsnf set integral to zero if this happened.
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSN1, DX =', DX
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
c fast processing for very small NX -- give ultra-conservative error estimates --
      GOTO (20, 30, 40), NX-1
   20 SMPSN1 = (FN(X(1)) + FN(X(2))) * DX / 2.D0
      ERR = ABS(SMPSN1)
      RETURN
   30 SMPSN1 = (FN(X(1)) + 4.D0 * FN(X(2)) + FN(X(3))) * DX / 3.D0
      ERR = ABS(SMPSN1)
      RETURN
   40 SMPSN1 = (( FN(X(1)) + 4.D0 * FN(X(2)) +     FN(X(3))) / 3.D0
     > + (-FN(X(2)) + 8.D0 * FN(X(3)) + 5.D0 * FN(X(4))) / 12.D0 ) * DX
      ERR = ABS(SMPSN1)
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSN1 = (FN(X(1)) + 4.D0*SE + 2.D0*SO + FN(X(NX))) * DX/3.D0
        TRPZ = (FN(X(1)) + 2.D0*(SE + SO) + FN(X(NX))) * DX/2.D0
      ELSE
        SMPSN1 =(FN(X(1)) + 4.D0*SE + 2.D0*SO + FN(X(NM1))) * DX/3.D0
     > +(-FN(X(NM1-1)) + 8.D0*FN(X(NM1)) + 5.D0*FN(X(NX))) * DX/12.D0
        TRPZ = (FN(X(1)) + 2.D0*(SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX/2.D0
      ENDIF

      ERR = SMPSN1 - TRPZ
c ======================================================================
c print the points...
c	do i = 1, nx
c	   xx = x(i)
c	   ff = fn(xx)
c	   write(nout,666) i, xx, ff
666	   format(1x,'smpsn1:',i5,1x,e12.5,1x,e12.5)
c	enddo
c ======================================================================

      RETURN

   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >  'out of range in SIMPSON INTEGRATION SMPSN1')
      STOP
C                        ****************************
      END

      FUNCTION SMPSN2(FN, A, B, NX, ERR)
c modified version of smpsnf (jcp 11/27/01)
C
C                       Does integral of FN(X)*dx from A TO B by SIMPSON'S METHOD
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / IOUNIT / NIN, NOUT, NWRT

      PARAMETER (MXPT = 10000)

      DIMENSION X(MXPT)
      external fn

      IF ((NX .LT. 0) .OR. (NX .GT. MXPT)) GOTO 99
C
      DX = (B - A) / FLOAT(NX-1)

c print warning if A > B; but routine SHOULD be ok either way.
c original routine smpsnf set integral to zero if this happened.
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSN2, DX =', DX
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
c fast processing for very small NX -- give ultra-conservative error estimates --
      GOTO (20, 30, 40), NX-1
   20 SMPSN2 = (FN(X(1)) + FN(X(2))) * DX / 2.D0
      ERR = ABS(SMPSN2)
      RETURN
   30 SMPSN2 = (FN(X(1)) + 4.D0 * FN(X(2)) + FN(X(3))) * DX / 3.D0
      ERR = ABS(SMPSN2)
      RETURN
   40 SMPSN2 = (( FN(X(1)) + 4.D0 * FN(X(2)) +     FN(X(3))) / 3.D0
     > + (-FN(X(2)) + 8.D0 * FN(X(3)) + 5.D0 * FN(X(4))) / 12.D0 ) * DX
      ERR = ABS(SMPSN2)
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSN2 = (FN(X(1)) + 4.D0*SE + 2.D0*SO + FN(X(NX))) * DX/3.D0
        TRPZ = (FN(X(1)) + 2.D0*(SE + SO) + FN(X(NX))) * DX/2.D0
      ELSE
        SMPSN2 =(FN(X(1)) + 4.D0*SE + 2.D0*SO + FN(X(NM1))) * DX/3.D0
     > +(-FN(X(NM1-1)) + 8.D0*FN(X(NM1)) + 5.D0*FN(X(NX))) * DX/12.D0
        TRPZ = (FN(X(1)) + 2.D0*(SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX/2.D0
      ENDIF

      ERR = SMPSN2 - TRPZ
c ======================================================================
c print the points...
c	do i = 1, nx
c	   xx = x(i)
c	   ff = fn(xx)
c	   write(nout,666) i, xx, ff
666	   format(1x,'smpsn2:',i5,1x,e12.5,1x,e12.5)
c	enddo
c ======================================================================

      RETURN

   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >  'out of range in SIMPSON INTEGRATION SMPSN2')
      STOP
C                        ****************************
      END

      FUNCTION SMPSNA (NX, DX, F, ERR)
C                                                   -=-=- smpsna
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (MAXX = 1000)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DIMENSION F(NX)
      DATA IW1, IW2, TINY / 2*0, 1.E-35 /
C 
      IF (DX .LE. 0.) THEN
        CALL WARNR(IW2,NWRT,'DX cannot be < 0. in SMPSNA', 'DX', DX,
     >               D0, D1, 0)
        SMPSNA = 0.
        RETURN
      ENDIF
 
      IF (NX .LE. 0 .OR. NX .GT. MAXX) THEN
        CALL WARNI(IW1, NWRT, 'NX out of range in SMPSNA', 'NX', NX,
     >               1, MAXX, 1)
        SIMP = 0.
      ELSEIF (NX .EQ. 1) THEN
        SIMP = 0.
      ELSEIF (NX .EQ. 2) THEN
        SIMP = (F(1) + F(2)) / 2.
        ERRD = (F(1) - F(2)) / 2.
      ELSE
        MS = MOD(NX, 2)
 
C For odd # of intervels, compute the diff between the Simpson 3/8-rule result
C for the last three bins and the regular Simpson rule result for the next to 
C last two bins.  The problem is thereby reduced to a even-bin one in all cases.
 
        IF (MS .EQ. 0) THEN
          ADD = (9.*F(NX) + 19.*F(NX-1) - 5.*F(NX-2) + F(NX-3)) / 24.
          NZ = NX - 1
        ELSE
          ADD = 0.
          NZ = NX
        ENDIF
 
        IF (NZ .EQ. 3) THEN
          SIMP = (F(1) + 4.* F(2) + F(3)) / 3.
          TRPZ = (F(1) + 2.* F(2) + F(3)) / 2.
        ELSE
          SE = F(2)
          SO = 0
          NM1 = NZ - 1
          DO 60 I = 4, NM1, 2
            IM1 = I - 1
            SE = SE + F(I)
            SO = SO + F(IM1)
   60     CONTINUE
          SIMP = (F(1) + 4.* SE + 2.* SO + F(NZ)) / 3.
          TRPZ = (F(1) + 2.* (SE + SO) + F(NZ)) / 2.
        ENDIF
 
        ERRD = TRPZ - SIMP 
        SIMP = SIMP + ADD
 
      ENDIF
C
      SMPSNA = SIMP * DX
 
      IF (ABS(SIMP) .GT. TINY) THEN
        ERR = ERRD / SIMP
      ELSE
        ERR = 0.
      ENDIF
C
      RETURN
C                        ****************************
      END
 
      FUNCTION SMPSNF (FN, A, B, NX, ERR, IER)
C                                                   -=-=- smpsnf
C
C                       Does integral of F(X)dx from A TO B by the SIMPSON METHO
C
C                       Double precision version of SMPSN
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
C
      PARAMETER (MXPT = 1000)
C
      DIMENSION X(MXPT)
      external fn

C
      IF (NX .EQ. 1) Then
         SmpsnF = 0
         RETURN
      EndIf

      IF (NX .LT. 0 .OR. NX .GT. MXPT) GOTO 99
C
      DX = (B - A) / (NX-1)
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSNF, DX =', DX
        SMPSNF = 0
        RETURN
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
      GOTO (20, 30, 40), NX-1
   20 SMPSNF = (FN(X(1)) + FN(X(2))) * DX / 2D0
      RETURN
   30 SMPSNF = (FN(X(1)) + 4D0 * FN(X(2)) + FN(X(3))) * DX / 3D0
      RETURN
   40 SMPSNF = (( FN(X(1)) + 4D0 * FN(X(2)) +     FN(X(3))) / 3D0
     > + (-FN(X(2)) + 8D0 * FN(X(3)) + 5D0 * FN(X(4))) / 12D0 ) * DX
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSNF = (FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NX))) * DX / 3D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO) + FN(X(NX))) * DX / 2D0
      ELSE
        SMPSNF =(FN(X(1)) + 4D0 * SE + 2D0 * SO + FN(X(NM1))) * DX / 3D0
     > +(-FN(X(NM1-1)) + 8D0 * FN(X(NM1)) + 5D0 * FN(X(NX))) * DX / 12D0
        TRPZ = (FN(X(1)) + 2D0 * (SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX / 2D0
      ENDIF
C
      ERR = SMPSNF - TRPZ
C
      RETURN
C
   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >              'out of range in SIMPSON INTEGRATION ROUTINE')
      STOP
C                        ****************************
      END

      SUBROUTINE SORT2(N,RA,RB)
C                                                   -=-=- sort2
C                                                  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION RA(N),RB(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
C                    ************************
      END

      SUBROUTINE SORTD (NCUR, NUMPT, XPT, YPT, SX, SY, IACT)
C                                                   -=-=- sortd
     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
     
      DIMENSION XPT(NUMPT), YPT(NUMPT, NCUR)
      DIMENSION SX (NUMPT, NCUR), SY (NUMPT, NCUR)
     
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / RANGED / XMIN, XMAX, DX, YMIN, YMAX, DY
     
      DATA BIGNUM / 1.0 E30 /
     
      IF (IACT .NE. 0 .AND. IACT .NE. 1)
     >   WRITE (NOUT, *) 'Illegel value of Iact in SORTD, Iact =', IACT
     
      IF (IACT .EQ. 0) THEN
         DO 10 IC = 1, NCUR
            SY(1, IC) = YPT(1, IC)
            SX(1, IC) = XPT(1)
            DO 20 IP = 2, NUMPT
               DO 30 IP1 = 1, IP - 1
                  IF (YPT(IP, IC) .GE. SY(IP1, IC)) THEN
                     YTEM = SY(IP1, IC)
                     XTEM = SX(IP1, IC)
                     SY(IP1, IC) = YPT(IP, IC)
                     SX(IP1, IC) = XPT(IP)
                     DO 40 IP2 = IP1 + 1, IP
                        XX = SX(IP2, IC)
                        YY = SY(IP2, IC)
                        SX(IP2, IC) = XTEM
                        SY(IP2, IC) = YTEM
                        XTEM = XX
                        YTEM = YY
   40                CONTINUE
                     GOTO 20
                  END IF
   30          CONTINUE
               SX(IP, IC) = XPT(IP)
               SY(IP, IC) = YPT(IP, IC)
   20       CONTINUE
   10    CONTINUE
      ENDIF
C                               Find range in x (that is, find xmin and xmax)
1     XMIN =  BIGNUM
      XMAX = -BIGNUM
      DO 31 I = 1,NUMPT
         XMIN = MIN(XMIN, XPT(I))
         XMAX = MAX(XMAX, XPT(I))
31    CONTINUE
c                               Find range in y (ymin and ymax)
      YMIN =  BIGNUM
      YMAX = -BIGNUM
      DO 21 J = 1, NCUR
         YMIN = MIN (YMIN, SY(NUMPT, J))
         YMAX = MAX (YMAX, SY(1, J))
21    CONTINUE
     
      RETURN
C                       ****************************
      END
C
      SUBROUTINE SORTI (NCUR, NUMPT, MPT, LPT, MST, LST)
C                                                   -=-=- sorti
C                                Integer version of SORT1 in the PLOTT package.
C                      Sorts the Arrays MPT and LPT in decreasing order of LPT.
C                                                           Output MST and LST
      DIMENSION MPT(NUMPT), LPT(NUMPT, NCUR)
      DIMENSION MST (NUMPT, NCUR), LST (NUMPT, NCUR)

      COMMON / IOUNIT / NIN, NOUT, NWRT

      DATA BIGNUM / 999999 /

C      IF (IACT .NE. 0 .AND. IACT .NE. 1)
C     >   WRITE (NOUT, *) 'Illegel value of Iact in SORT1, Iact =', IACT

      DO 10 IC = 1, NCUR
      LST(1, IC) = LPT(1, IC)
      MST(1, IC) = MPT(1)
                DO 20 IP = 2, NUMPT
                        DO 30 IP1 = 1, IP - 1
                        IF (LPT(IP, IC) .GE. LST(IP1, IC)) THEN
                                YTEM = LST(IP1, IC)
                                XTEM = MST(IP1, IC)
                                LST(IP1, IC) = LPT(IP, IC)
                                MST(IP1, IC) = MPT(IP)
                                DO 40 IP2 = IP1 + 1, IP
                                        XX = MST(IP2, IC)
                                        YY = LST(IP2, IC)
                                        MST(IP2, IC) = XTEM
                                        LST(IP2, IC) = YTEM
                                        XTEM = XX
                                        YTEM = YY
   40                           CONTINUE
                                GOTO 20
                        END IF
   30           CONTINUE
                        MST(IP, IC) = MPT(IP)
                        LST(IP, IC) = LPT(IP, IC)
   20   CONTINUE
   10 CONTINUE
C                             Find range in x (that is, find MMin and MMax)
1     MMIN =  BIGNUM
      MMAX = -BIGNUM
      DO 31 I = 1,NUMPT
                MMIN = MIN(MMIN, MPT(I))
                MMAX = MAX(MMAX, MPT(I))
31    CONTINUE
c                               Find range in y (LMin and LMax)
      LMIN =  BIGNUM
      LMAX = -BIGNUM
      DO 21 J = 1, NCUR
        LMIN = MIN (LMIN, LST(NUMPT, J))
        LMAX = MAX (LMAX, LST(1, J))
21    CONTINUE

      RETURN
C                       **********************
      END
C                                                          =-=-= Askwrn
       SUBROUTINE TRMSTR(STRING,ILEN)
C                                                   -=-=- trmstr

C      Removes leading spaces and returns true length of a character string

       CHARACTER STRING*(*),SPACE*1

       DATA SPACE/' '/

       ILEN=0

       IF (STRING.EQ.SPACE) RETURN
C                                          Remove leading spaces
1      IF (STRING(1:1).NE.SPACE) GOTO 2
          STRING=STRING(2:)
          GOTO 1
2      CONTINUE
C                                          Count up trailing spaces
       DO 3 I=LEN(STRING),1,-1
          IF (STRING(I:I).NE.SPACE) THEN
             ILEN=I
             RETURN
          END IF
3      CONTINUE
C               *************************
       END

      FUNCTION TRNGLE (X, Y, Z, Irt)
C                                                   -=-=- trngle
C                                       "Triangle Function" of the three sides. 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DATA RDOFF / 1E-11 /

      IRT = 0

      Entry TNGLE (X, Y, Z)

      AMX = MAX (X*X, Y*Y, Z*Z) 

      TMP = X*X + Y*Y + Z*Z - 2.* (X*Y + Y*Z + Z*X)
      ATMP= ABS(TMP) 
      
      IF     (ATMP .LT. AMX * RDOFF) THEN
        TMP = ATMP
        IRT = 1
      ELSEIF (TMP .LT. 0.) THEN
        PRINT '(A, 4(1PE12.3))', 'X,Y,Z, TMP =', X, Y, Z, TMP
        STOP 'Negative argument in TRNGLE function; check for errors!'
      ENDIF

      Tmp = SQRT (TMP)

      TNGLE  = Tmp
      TRNGLE = Tmp
 
      RETURN
C                  ****************************
      END

        SUBROUTINE UC(A)
C                                                   -=-=- uc
C               Converts A to all upper case.
C               System dependent: assumes ICHAR(uc letter) - ICHAR(lc letter)
C               is constant.
        CHARACTER A*(*), C*(1)
        INTEGER I

        ENTRY UPCASE(A)
        DO 1 I=1, LEN(A)
                C = A(I:I)
C
C                       Use ASCII ordering for detecting lc:
                IF ( LGE(C, 'a') .AND. LLE(C, 'z') )THEN
                        A(I:I) = CHAR(ICHAR(C)+ICHAR('A')-ICHAR('a'))
                        ENDIF
 1              CONTINUE
        RETURN
C               *************************
        END
C
        CHARACTER*(*) FUNCTION UCASE(A)
C               Converts A to all upper case.
        CHARACTER*(*) A
        UCASE = A
        CALL UC(UCASE)
        RETURN
C                           **********************      
        END
C
      SUBROUTINE UpC (A, La, UpA)
C                                                   -=-=- upc

C  5/29/94 WKT

C  This is a variation of the two old routines UC(A) and UpCase(A). Here
C  the converted value is return to the new variable UpA, rather than
C  the input dummy variable A, to avoid bombing the program with error
C  "access violation, reason mask=04" when the routine is given a CONSTANT
C  actual argument (rather than a character variable argument).

C  The inconvenience of this new version is that the calling program must
C  declare the length of UpA explicitly; and it better be .Ge. Len(A).

C  We trim any excess trailing characters in UpA and replace with spaces.
C  To be exact, the returned value should be considered as UpA(1:La).

      CHARACTER A*(*), UpA*(*), C*(1)
      INTEGER I, La, Ld

      La = Len(A)
      Lb = Len(UpA)

      If (Lb .Lt. La) Stop 'UpCase conversion length mismatch!'

      Ld = ICHAR('A')-ICHAR('a')

      DO 1 I = 1, Lb

        If (I .Le. La) Then
         c = A(I:I)
         IF ( LGE(C, 'a') .AND. LLE(C, 'z') ) THEN

           UpA (I:I) = CHAR(Ichar(c) + ld)
         Else
           UpA (I:I) = C

         ENDIF
        Else
         UpA (I:I) = ' '
        Endif

 1    CONTINUE
      
      RETURN
C               *************************
      END

C                                                          =-=-= Nurcpe
      SUBROUTINE WARNI (IWRN, NWRT, MSG, NMVAR, IVAB,
C                                                   -=-=- warni
     >                  IMIN, IMAX, IACT)

C     t++++++++++++++++++++++++++++++++     Routines to handle Warnings
C                                                  Integer version
      CHARACTER*(*) MSG, NMVAR

      Save Iw

      Data Nmax / 100 /

      IW = IWRN
      IV = IVAB
      
      IF  (IW .EQ. 0) THEN
         PRINT '(1X,A/1X, 2A,I10 /A,I4)', MSG, NMVAR, ' = ', IV,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2I10)', ' The limits are: ', IMIN, IMAX
         WRITE (NWRT,'(A/2I10)') ' The limits are: ', IMIN, IMAX
         ENDIF
      ENDIF

      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(1X,A/1X,2A, I10)') MSG, NMVAR, ' = ', IV
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     'Too many warnings, Message suppressed !!'
      Endif

      IWRN = IW + 1

      RETURN
C               *************************
      END

      SUBROUTINE WARNR (IWRN, NWRT, MSG, NMVAR, VARIAB,
C                                                   -=-=- warnr
     >                  VMIN, VMAX, IACT)
 
C      Subroutine to handle warning messages.  Writes the (warning) message
C      and prints out the name and value of an offending variable to SYS$OUT
C      the first time, and to output file unit # NWRT in subsequent times.
C      
C      The switch IACT decides whether the limits (VMIN, VMAX) are active or
C      not.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
 
      CHARACTER*(*) MSG, NMVAR

      Save Iw

CCPY      Data Nmax / 100 /
      Data Nmax / 5 /

      IW = IWRN
      VR = VARIAB
 
      IF  (IW .EQ. 0) THEN
         PRINT '(1X, A/1X,2A,1PD16.7/A,I4)', MSG, NMVAR, ' = ', VR,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2(1PE15.4))', ' The limits are: ', VMIN, VMAX
         WRITE (NWRT,'(A/2(1PE15.4))') ' The limits are: ', VMIN, VMAX
         ENDIF
      ENDIF
 
      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(I5, 2A/1X,2A,1PD16.7)') IW, '   ', MSG,
     >                  NMVAR, ' = ', VR
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     '!! Too many warnings, Message suppressed from now on !!'
      Endif

      IWRN = IW + 1
 
      RETURN
C                  ****************************
      END

C                                                          =-=-= Inputs
      FUNCTION ZBRNT(FUNC, X1, X2, TOL, IRT)
C                                                   -=-=- zbrnt
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)
      external func
 
      IRT = 0
      TOL = ABS(TOL)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        PRINT *, 'Root must be bracketed for ZBRNT. Set = 0'
        IRT = 1
        ZBRNT=0.
        RETURN
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
          ZBRNT=B
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
      PRINT *, 'ZBRNT exceeding maximum iterations.'
      IRT = 2
      ZBRNT=B
      RETURN
      END
C                        ****************************
