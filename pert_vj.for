CPY Oct 04, 2008, Add the PERT part for the processes Z0_RA_UUB and Z0_RA_DDB 
C Feb 10, 2005, B(2) for bb-->HB process is added, and eliminate 
C "omega" in HB process which is not needed if the input
C B quark mass is the MS-bar running mass at M_B scale. 
C Error: April 14, 2004 (correct ALEPASY for A0, etc)
c-----------------------------------------------------------------------
c %%%%%%%%%%%%%%%%%           VERSION 6.1           %%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
C
CsB   Nov. 17, 1997 q Q -> Z Z process is added.
C
CsB   Mar. 30, 1997
C     This file and the numerical results for the Y piece were tested against
C     C.-P.'s version. This file contains all information what is in C.-P.'s
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
C     - the AA and AG pieces use 'fake' (DY type) Y pieces,
C       and in order to compare to experiments the real NLO pieces have to
C       be implemented. (It is under way as of 3/30/97),
C     - for AA and AG the A1-3 pieces in ppBar, pp and pN collisions are zero,
C     - for H0 production all the perturbative pieces are the same
C       (for any hadronic initial state),
C     - for pp collision the A1 and A3 pieces are zero (for all bosons)
C       when y = 0,
C     - agreement of GG and GL options with C.-P.
C
C
CCPY  May 1996 -- a major bug was found for (P N) process
C                 inside PNX1X2INT(X1,X2...)
C
CsB   Feb.  7, 1996 -- L0 and A3 pieces now ckecked to cancel in ppBar,
C                      pp and pN cases for W+/-, Z and photon.
c
CsB   Jan. 27, 1996 -- Implemented full regular perturbative pieces A1-4.
C                      L0,A3 cancels with L0,A3 of asmptotic piece for low qT.
C                      A2 checks against Mirkes.
c
CsB   Dec. 12, 1995 -- Separation of qqBar and qG processes implemented
c                      by the switch I_Proc
c
CsB   Dec.  5, 1995 -- Installed asymmetric pieces in parton luminosities.
c
CCPY  Oct. 23, 1995 -- Added LTOPT option for LO result
C                      (refer to WZNLO.FOR)
CCPY  Oct. 14, 1995 -- Added virtual gluon "resonance"
C
CCPY This code combines my original PERT.FOR and PERT_PP.FOR.
C    It does calculations for either ppbar or pp or pN collisions.
C
CJI  Sep. 4, 2015 -- Created new file from pert.for for processes with
C                    jets
      function fgetpert_vj(vmas,vpt,rapin,rapin2,gees,q0in,fresum,pma,
     >                     pert,asy,ier_pert,iflavor)
      implicit none
      include 'common.for'
      integer ier_pert
      real*8 fgetpert_vj
      real*8 vmas,vpt,rapin,rapin2,gees(3),q0in,fresum,pma,pert,asy
      integer flavor, iflavor
      common/flavor/flavor

      flavor = iflavor

      if(ibeam.eq.-1) then
C For P Pbar scattering
        call spertppbar_vj(vmas,vpt,rapin,rapin2,gees,q0in,fresum,pma,
     >                     pert,asy,ier_pert)
      else
C For P-P, P-Nucleus, pion-Nucleus scatterings
        call spertpN_vj(vmas,vpt,rapin,rapin2,gees,q0in,fresum,pma,
     >                  pert,asy,ier_pert)
      endif
c returns the contribution of PERTURBATIVE-ASMYPTOTIC
      fgetpert_vj=pma
	call flush()
      return
      end

      subroutine spertppbar_vj(vmas,vpt,rapin,rapin2,gees,q0in,
     >                     fresum,pma,gal_pert,gal_asymp,ier_pert)
Cgal: OCT 1993, JULY 1994
CCPY OCTOBER 1992, JULY 1993
C THIS IS MODEIFIED TO ONLY WORK FOR P PBAR MACHINE
C THIS IS PERT_PPBAR.FOR

      IMPLICIT NONE
      INCLUDE 'common.for'

      INTEGER NF_EFF
      REAL*8 PERT
      REAL*8 ASYMP,AMU
      INTEGER IRUN,NRUN
csb      REAL*8 EW_ALFA
      LOGICAL DUMB
cgal:
      integer ier_pert
      real*8 vmas,vpt,rapin,gees(3),q0in,fresum,pma,gal_pert,gal_asymp
      REAL*8 rapin2
      Real*8 Coup_M, Omega
      Double Precision RUN_MASS_TOT
      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC

      Integer KinCorr
      Common / CorrectKinematics / KinCorr
  
      real*8 pyalem
      external pyalem

      Logical First_1
      Data First_1 /.True./
      
      real*8 Z_beta0, Z_beta1

CJI Jan 2015: Add in hj calculation
      Integer HasJet
      REAL*8 D1s, R, ptj,mTran,s,t,u
      REAL*8 H0Approx,H0,H1,HB1,tb,ub
      Common /Jet/ D1s, R, t, ptj, HasJet
      REAL*8 KK,dsigmadt,Kgg,tau,GHJ2,xLi,Kgq
      External H0Approx,xLi
      Real*8 Li1, Li2, Li3,z,H

      integer iflavor
      common/flavor/iflavor

      ier_pert=0
Cgal: no longer do we enforce lowest order QCD; technically wrong here, but
c     numerically insignificant
c      IF(NORDER.EQ.2.and.rdflag) THEN
c       WRITE(*,*) ' FORCE NORDER = 1 IN PERT.FOR FOR ALL PDF'
c       WRITE(NWRT,*) ' FORCE NORDER = 1 IN PERT.FOR FOR ALL PDF'
c       NORDER=1
c     ENDIF

      vmas = 60
      QT_v = 10
      y_v = 0
      y2_v = 0
      ECM = 7000

      y_v=rapin
      y2_v=rapin2
      ptj=vmas
      mTran=dsqrt(ptj**2+mh**2)
      Q_V=dsqrt(mh**2+2d0*ptj**2+2d0*ptj*mTran*cosh(y_v-y2_v))
      S_Q=Q_V

      S_LAM=ALAMBD(S_Q)
      NF_EFF=NFL(S_Q)
      CALL VJSETUP(NF_EFF)
c
      S_BETA1=(33.0-2.0*NF_EFF)/12.0

c      ECM=ECM*UNIT
      ECM=ECM
      ECM2=ECM**2

CCPY SET NO_SIGMA0=0 FOR INCLUDING CONST=\sigma_0 IN THE RATES
C    SET NO_SIGMA0=1 FOR NOT INCLUDING CONST=\sigma_0 IN THE RATES
CsB        NO_SIGMA0=0

CCPY SET LTOPT=0 FOR CALCULATING THE ASYMPTOTIC PART
C    SET LTOPT=1 FOR CALCULATING THE DELTA_SIGMA FROM QT=0 TO PT
C    SET LTOPT=-1 FOR CALCULATING THE LEADING ORDER RESULT

CCPY
       IF(TYPE_V.eq.'HJ') THEN
          KK = 1d0 + ALPI(S_Q)*11d0/4d0
          tau = 4d0 * mt**2/mh**2
          GHJ2 = asin(1d0/sqrt(tau))**2
          GHJ2 = 1d0+(1d0-tau)*GHJ2
          GHJ2 = 4d0*dsqrt(2d0)*(ALPI(S_Q)/4d0)**2*GMU*tau**2*GHJ2**2
        if(iflavor.eq.0) then
          Kgg = 1d0/4d0/(NC**2-1d0)**2
          dsigmadt=GHJ2*ALPI(S_Q)/4d0*Kgg
        else if(iflavor .eq.1) then
          Kgq = 1d0/4d0/(NC**2-1d0)/NC
          dsigmadt=GHJ2*ALPI(S_Q)/4d0*Kgq
        endif
          const = PI*dsigmadt/S_Q**4
      ELSE
        PRINT*,' WRONG TYPE_V'
        CALL QUIT
      ENDIF

      DUMB=.false.
      if(DUMB) then
      WRITE(nout,*) ' ECM,IBEAM,LEPASY,LTOPT,No_Sigma0'
      WRITE(nout,*) ECM,IBEAM,LEPASY,LTOPT,NO_SIGMA0
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

cgal: ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cgal: enters the values of the parameters AFTER the setup routine

	qt_v=vpt

cgal: ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cgal:      READ(NINR,*) QT_V,Y_V
      IF(MOD(IRUN,5).EQ.0) THEN
       WRITE(*,*)QT_V,Y_V,IRUN
      ENDIF

      PERT=0.D0
      ASYMP=0.D0

      CALL YMAXIMUM
      IF(DABS(Y_V).GT.YMAX) THEN
C       WRITE(NWRT,*) 'ABS(Y_V)  > YMAX'
C       WRITE(NWRT,*)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
CSM        WRITE(NOUT,*) 'ABS(Y_V)  > YMAX'
CSM        WRITE(NOUT,910)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
       ier_pert=1
       GOTO 900
      elseIF(abs(abs(ymax)-ABS(Y_V)).LT.YMAX/100) THEN
c if rapidity is close to YMAX, then we are in a region where the program
c   is not accurate; we flag this situation
       ier_pert=2
      ENDIF

      IF(DABS(Y2_V).GT.YMAX) THEN
C       WRITE(NWRT,*) 'ABS(Y_V)  > YMAX'
C       WRITE(NWRT,*)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
CSM        WRITE(NOUT,*) 'ABS(Y_V)  > YMAX'
CSM        WRITE(NOUT,910)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
       ier_pert=1
       GOTO 900
      elseIF(abs(abs(ymax)-ABS(Y2_V)).LT.YMAX/100) THEN
c if rapidity is close to YMAX, then we are in a region where the program
c   is not accurate; we flag this situation
       ier_pert=2
      ENDIF

cpn Feb 2004       ACC_RERR=1.D-3
       ACC_RERR=1.D-4
C       ACC_RERR=1.D-5
       
       FIT3=.FALSE.

CsB___Definition of momentum fractions x1 and x2
C      If (KinCorr.Eq.0) then
        X_A=(ptj*DEXP(Y_V)+mTran*DEXP(Y2_V))/ECM
        X_B=(ptj*DEXP(-Y_V)+mTran*DEXP(-Y2_V))/ECM
C      Else If (KinCorr.Eq.1) then
C        X_A=Sqrt(Q_V**2+QT_V**2)/ECM*DEXP(Y_V)
C        X_B=Sqrt(Q_V**2+QT_V**2)/ECM/DEXP(Y_V)
C      End If

CCPY      AMU=C4*Q_V
      AMU=C4*S_Q

        s = x_A*X_B*ECM2
        t = -x_A*ECM*ptj*Dexp(-y_v)
        u = -X_B*ECM*ptj*Dexp(y_v)
        tb = mh**2-t
        ub = mh**2-u

       IF(TYPE_V.eq.'HJ') THEN
           H0 = NC*(NC**2-1)*(s**4+t**4+u**4+mh**8)/(s*t*u)

           z=1d0-mh**2/s
           if(dabs(z).gt.1d0) then
               Li1 = -xLi(2,1d0/z)-PI2/6d0-0.5*dlog(-z)**2
           else
               Li1 = xLi(2,z)
           endif

           z=t/mh**2
           if(dabs(z).gt.1d0) then
               Li2 = -xLi(2,1d0/z)-PI2/6d0-0.5*dlog(-z)**2
           else
               Li2 = xLi(2,z)
           endif

           z=u/mh**2
           if(dabs(z).gt.1d0) then
               Li3 = -xLi(2,1d0/z)-PI2/6d0-0.5*dlog(-z)**2
           else
               Li3 = xLi(2,z)
           endif

           H1 = CA*(67D0/9D0-23d0*nf/54d0+PI2/2d0+Log(1d0/R**2)*
     >          Log(s/ptj**2)+Log(s/ptj**2)**2+2d0/3d0*S_BETA1*
     >          Log(s/(ptj**2*R**2))+Log(tb/mh**2)**2-
     >          Log(-tb/t)**2-2d0*Log(-t/s)*Log(-u/s)+
     >          Log(ub/mh**2)**2-Log(-ub/u)**2+2d0*Li1+
     >          2d0*Li2+2d0*Li3)
           HB1 = mh**2*NC*(NC**2-1)*(NC-NF_EFF)*(s*t*u+mh**2*
     >          (s*t+s*u+t*u))/(3d0*s*t*u)
           H = H0*(1d0+ALPI(S_Q)/2d0*(H1+HB1/H0))
       endif

      CALL ASYMPTO_PN(X_A,X_B,AMU,ASYMP)

      IF(DEBUG) THEN
       WRITE(NWRT,*) ' ACC_RERR =',ACC_RERR
       WRITE(NWRT,*) ' AMU,C4,Q_V,S_Q,C1,C2,C3'
       WRITE(NWRT,*) AMU,C4,Q_V,S_Q,C1,C2,C3
       IF(FIT3) THEN
        WRITE(NWRT,*) 'FIT3 = TRUE'
       ELSE
        WRITE(NWRT,*) 'FIT3 = FALSE'
       ENDIF
      ENDIF

C      WRITE(NWRT,*)'QT_V,Y_V,Q_V,PERT,ASYMP,YMAX,X_A,X_B,NORDER'
C      WRITE(NWRT,*)QT_V,Y_V,Q_V,PERT,ASYMP,YMAX,X_A,X_B,NORDER
C      WRITE(NWRT,*)' '

CsB      WRITE(NOUT,920)QT_V,Y_V,PERT,ASYMP,YMAX,NORDER
920   FORMAT(F7.2,2X,F7.2,2X,D16.4,2X,D16.4,2X,F7.2,2X,I3)
C      IF(ASK('MORE')) GOTO 10

900   CONTINUE

cgal: send back d(SIGMA)/d(Q^2)d(QT)d(Y)
      gal_asymp=asymp

      if(HASJET.EQ.1) then
          gal_asymp = gal_asymp*x_a*x_b*ECM2*H
      endif

      END

C --------------------------------------------------------------------------
      subroutine spertpN_vj(vmas,vpt,rapin,rapin2,gees,q0in,
     >                    fresum,pma,gal_pert,gal_asymp,ier_pert)
C --------------------------------------------------------------------------
c
c The subroutine version
c
C THIS IS PERT_PN.FOR
C ADD DRELL-YAN (PHOTON) PROCESS
Cgal: OCT 1993, JULY 1994
CCPY OCTOBER 1992, JULY 1993
C THIS IS MODEIFIED TO ONLY WORK FOR (P P) AND (P NUCLEUS) MACHINE
      IMPLICIT NONE
      INCLUDE 'common.for'

      INTEGER NF_EFF
      REAL*8 PERT
      REAL*8 ASYMP,AMU
      INTEGER IRUN,NRUN
csb      REAL*8 EW_ALFA
      LOGICAL DUMB
cgal:
      integer ier_pert
      real*8 vmas,vpt,rapin,gees(3),q0in,fresum,pma,gal_pert,gal_asymp
      REAL*8 rapin2
      Real*8 Coup_M, Omega
      Double Precision RUN_MASS_TOT
      REAL*8 XMTOP,XMBOT,XMC
      COMMON/XMASS/XMTOP,XMBOT,XMC

      Integer KinCorr
      Common / CorrectKinematics / KinCorr
  
      real*8 pyalem
      external pyalem

      Logical First_1
      Data First_1 /.True./
      
      real*8 Z_beta0, Z_beta1

CJI Jan 2015: Add in hj calculation
      Integer HasJet
      REAL*8 D1s, R, ptj,mTran,s,t,u
      REAL*8 H0Approx,H0,H1,HB1,tb,ub
      Common /Jet/ D1s, R, t, ptj, HasJet, u
      REAL*8 KK,dsigmadt,Kgg,tau,GHJ2,xLi,Kgq
      External H0Approx,xLi
      Real*8 Li1, Li2, Li3,z,H

      integer iflavor
      common/flavor/iflavor

      ier_pert=0
Cgal: no longer do we enforce lowest order QCD; technically wrong here, but
c     numerically insignificant
c      IF(NORDER.EQ.2.and.rdflag) THEN
c       WRITE(*,*) ' FORCE NORDER = 1 IN PERT.FOR FOR ALL PDF'
c       WRITE(NWRT,*) ' FORCE NORDER = 1 IN PERT.FOR FOR ALL PDF'
c       NORDER=1
c     ENDIF

      y_v=rapin
      y2_v=rapin2
      ptj=vmas
      mTran=dsqrt(ptj**2+mh**2+QT_V**2)
      Q_V=dsqrt(mh**2+2d0*ptj**2+2d0*ptj*mTran*cosh(y_v-y2_v))
      S_Q=Q_V
      t = -ptj**2-exp(-y_v+y2_v)*ptj*mTran
      u = -ptj**2-exp(y_v-y2_v)*ptj*mTran

      S_LAM=ALAMBD(S_Q)
      NF_EFF=NFL(S_Q)
      CALL VJSETUP(NF_EFF)
c
      S_BETA1=(33.0-2.0*NF_EFF)/12.0

c      ECM=ECM*UNIT
      ECM=ECM
      ECM2=ECM**2

CCPY SET NO_SIGMA0=0 FOR INCLUDING CONST=\sigma_0 IN THE RATES
C    SET NO_SIGMA0=1 FOR NOT INCLUDING CONST=\sigma_0 IN THE RATES
CsB        NO_SIGMA0=0

CCPY SET LTOPT=0 FOR CALCULATING THE ASYMPTOTIC PART
C    SET LTOPT=1 FOR CALCULATING THE DELTA_SIGMA FROM QT=0 TO PT
C    SET LTOPT=-1 FOR CALCULATING THE LEADING ORDER RESULT

CCPY
       IF(TYPE_V.eq.'HJ') THEN
          KK = 1d0! + ALPI(S_Q)*11d0/4d0
          tau = 4d0 * mt**2/mh**2
          GHJ2 = asin(1d0/sqrt(tau))**2
          GHJ2 = 1d0+(1d0-tau)*GHJ2
          GHJ2 = 4d0*dsqrt(2d0)*(ALPI(125d0)/4d0)**2*GMU*4.0/9.0!tau**2*GHJ2**2
        if(iflavor.eq.0) then
          Kgg = 1d0/4d0/(NC**2-1d0)**2
          dsigmadt=GHJ2*ALPI(125d0)/4d0*Kgg
        else if(iflavor .eq.1) then
          Kgq = 1d0/4d0/(NC**2-1d0)/NC
          dsigmadt=GHJ2*ALPI(125d0)/4d0*Kgq
        endif
          const = PI*dsigmadt/S_Q**4
      ELSE
        PRINT*,' WRONG TYPE_V'
        CALL QUIT
      ENDIF

      DUMB=.false.
      if(DUMB) then
      WRITE(nout,*) ' ECM,IBEAM,LEPASY,LTOPT,No_Sigma0'
      WRITE(nout,*) ECM,IBEAM,LEPASY,LTOPT,NO_SIGMA0
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

cgal: ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cgal: enters the values of the parameters AFTER the setup routine

	qt_v=vpt

cgal: ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cgal:      READ(NINR,*) QT_V,Y_V
      IF(MOD(IRUN,5).EQ.0) THEN
       WRITE(*,*)QT_V,Y_V,IRUN
      ENDIF

      PERT=0.D0
      ASYMP=0.D0

      CALL YMAXIMUM
      IF(DABS(Y_V).GT.YMAX) THEN
C       WRITE(NWRT,*) 'ABS(Y_V)  > YMAX'
C       WRITE(NWRT,*)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
CSM        WRITE(NOUT,*) 'ABS(Y_V)  > YMAX'
CSM        WRITE(NOUT,910)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
       ier_pert=1
       GOTO 900
      elseIF(abs(abs(ymax)-ABS(Y_V)).LT.YMAX/100) THEN
c if rapidity is close to YMAX, then we are in a region where the program
c   is not accurate; we flag this situation
       ier_pert=2
      ENDIF

      IF(DABS(Y2_V).GT.YMAX) THEN
C       WRITE(NWRT,*) 'ABS(Y_V)  > YMAX'
C       WRITE(NWRT,*)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
CSM        WRITE(NOUT,*) 'ABS(Y_V)  > YMAX'
CSM        WRITE(NOUT,910)QT_V,Y_V,RESUM,RESUM_1,I,RESUM_2
       ier_pert=1
       GOTO 900
      elseIF(abs(abs(ymax)-ABS(Y2_V)).LT.YMAX/100) THEN
c if rapidity is close to YMAX, then we are in a region where the program
c   is not accurate; we flag this situation
       ier_pert=2
      ENDIF

cpn Feb 2004       ACC_RERR=1.D-3
       ACC_RERR=1.D-4
C       ACC_RERR=1.D-5
       
       FIT3=.FALSE.

CsB___Definition of momentum fractions x1 and x2
C      If (KinCorr.Eq.0) then
        X_A=(ptj*DEXP(Y_V)+mTran*DEXP(Y2_V))/ECM
        X_B=(ptj*DEXP(-Y_V)+mTran*DEXP(-Y2_V))/ECM
C      Else If (KinCorr.Eq.1) then
C        X_A=Sqrt(Q_V**2+QT_V**2)/ECM*DEXP(Y_V)
C        X_B=Sqrt(Q_V**2+QT_V**2)/ECM/DEXP(Y_V)
C      End If

CCPY      AMU=C4*Q_V
      AMU=C4*S_Q

        s = x_A*X_B*ECM2
        t = -x_A*ECM*ptj*Dexp(-y_v)
        u = -X_B*ECM*ptj*Dexp(y_v)
        tb = mh**2-t
        ub = mh**2-u

       IF(TYPE_V.eq.'HJ') THEN
           if(iflavor.eq.0) then
           H0 = NC*(NC**2-1)*(s**4+t**4+u**4+mh**8)/(s*t*u)

           z=1d0-mh**2/s
           if(dabs(z).gt.1d0) then
               Li1 = -xLi(2,1d0/z)-PI2/6d0-0.5*dlog(-z)**2
           else
               Li1 = xLi(2,z)
           endif

           z=t/mh**2
           if(dabs(z).gt.1d0) then
               Li2 = -xLi(2,1d0/z)-PI2/6d0-0.5*dlog(-z)**2
           else
               Li2 = xLi(2,z)
           endif

           z=u/mh**2
           if(dabs(z).gt.1d0) then
               Li3 = -xLi(2,1d0/z)-PI2/6d0-0.5*dlog(-z)**2
           else
               Li3 = xLi(2,z)
           endif

           H1 = CA*(67D0/9D0-23d0*nf/54d0+PI2/2d0+Log(1d0/R**2)*
     >          Log(s/ptj**2)+Log(s/ptj**2)**2+2d0/3d0*S_BETA1*
     >          Log(s/(ptj**2*R**2))+Log(tb/mh**2)**2-
     >          Log(-tb/t)**2-2d0*Log(-t/s)*Log(-u/s)+
     >          Log(ub/mh**2)**2-Log(-ub/u)**2+2d0*Li1+
     >          2d0*Li2+2d0*Li3)
           HB1 = mh**2*NC*(NC**2-1)*(NC-NF_EFF)*(s*t*u+mh**2*
     >          (s*t+s*u+t*u))/(3d0*s*t*u)
           H = H0*(1d0+ALPI(S_Q)/2d0*(H1+HB1/H0))
           else if(iflavor.eq.1) then
           H0 = CA*CF/(-u)*(s**2+t**2)+CA*CF/(-t)*(s**2+u**2)
           endif
           CONST=CONST*H0*X_A*X_B
       endif

      CALL ASYMPTO_PN(X_A,X_B,AMU,ASYMP)

      IF(DEBUG) THEN
       WRITE(NWRT,*) ' ACC_RERR =',ACC_RERR
       WRITE(NWRT,*) ' AMU,C4,Q_V,S_Q,C1,C2,C3'
       WRITE(NWRT,*) AMU,C4,Q_V,S_Q,C1,C2,C3
       IF(FIT3) THEN
        WRITE(NWRT,*) 'FIT3 = TRUE'
       ELSE
        WRITE(NWRT,*) 'FIT3 = FALSE'
       ENDIF
      ENDIF

C      WRITE(NWRT,*)'QT_V,Y_V,Q_V,PERT,ASYMP,YMAX,X_A,X_B,NORDER'
C      WRITE(NWRT,*)QT_V,Y_V,Q_V,PERT,ASYMP,YMAX,X_A,X_B,NORDER
C      WRITE(NWRT,*)' '

CsB      WRITE(NOUT,920)QT_V,Y_V,PERT,ASYMP,YMAX,NORDER
920   FORMAT(F7.2,2X,F7.2,2X,D16.4,2X,D16.4,2X,F7.2,2X,I3)
C      IF(ASK('MORE')) GOTO 10

900   CONTINUE

cgal: send back d(SIGMA)/d(Q^2)d(QT)d(Y)
      gal_asymp=asymp

      if(HASJET.EQ.1) then
CJI     modify to pass back a flatter distribution
          gal_asymp = gal_asymp*ECM2*PTJ**2/(2d0*QT_V)
      endif
      end
