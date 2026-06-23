      subroutine alekhin(x,q,mode,uv,dv,usea,dsea,str,chm,btm,glu)
      implicit real*8(a-h,o-z)
      real*8 pdfs(0:9),dpdfs(0:9,15),delpdf(0:9)

c  the ALEKHIN02 VFN nominal sets: mode=1,2,3 for LO,NLO.NNLO
c  requires the following data files to be in the same directory
c     a02.pdfs_1_vfn       a02.pdfs_2_vfn      a02.pdfs_3_vfn
c     a02.dpdfs_1_vfn      a02.dpdfs_2_vfn     a02.dpdfs_3_vfn
c  returns x times parton distributions

      q2=q*q
      xb=x
      call a02(xb,q2,pdfs,dpdfs,npdf,npar,mode,1,0)
      uv=pdfs(1)
      dv=pdfs(2)
      usea=pdfs(4)
      dsea=pdfs(6)
      str=pdfs(5)
      chm=pdfs(7)
      btm=pdfs(8)
      glu=pdfs(3)

      return 
      end

      subroutine a02(xb,q2,PDFS,DPDFS,NPDF,NPAR,KORD,KSCHEM,KSET)
c
c     This is a code for the parton distributions with account 
c     of their experimental (stat+syst) and theoretical uncertainties. 
c     The Q**2 range is 1.4d0 < Q**2 < 2d8, the x range is 1d-7 < x < 1d0. 
c
c  Input parameters:
c        KORD=1 -- the LO PDFs
c        KORD=2 -- the NLO PDFs
c        KORD=3 -- the NNLO PDFs
c      
c        KSCHEM=0 -- the fixed-flavor-number (FFN) scheme 
c        KSCHEM=1 -- the variable-flavor-number (VFN) scheme
c
c        KSET=0 -- nominal PDFs
c        KSET=1 -- PDFs with mass of c-quark increased from 1.5 to 1.75 GeV
c        KSET=2 -- PDFs with the strange sea suppression factor increased from 
c                  0.42 to 0.52
c        KSET=3 -- PDFs with the choice B (slow evolution) for the NNLO kernel 
c                  (used with KORD=3 only)
c
c  Output parameters:
c     The array PDFS contains fitted values of the strong coupling constant 
c     and the parton distributions at given x and Q:
c        PDFS(0) -- \alpha_s
c        PDFS(1) -- valence u-quarks 
c        PDFS(2) -- valence d-quarks
c        PDFS(3) -- gluons 
c        PDFS(4) -- sea u-quarks 
c        PDFS(5) -- s-quarks 
c        PDFS(6) -- sea d-quarks 
c        PDFS(7) -- c-quarks
c        PDFS(8) -- b-quarks
c        PDFS(9) -- t-quarks
c     NPDF is the number of PDFs returned (NPDF=6 for the FFN PDFs and 9 for 
c     the VFN ones).
c     Output array DPDFS(0:ipdf,ipar) contains derivatives of \alpha_s and
c     the PDFs on the fitted parameters with the number of the parameters 
c     returned in NPAR. With the derivatives of \alpha_s included one can take 
c     into account the correlations of the fitted PDFs with \alpha_s as well.
c     All derivatives are transformed to the orthonormal 
c     basis of eigenvectors of the parameters error matrix. For this reason 
c     the variation of the PDFs in the derivatives directions can be performed 
c     independently. For example the dispersion of the i-th PDF can be stored 
c     in DELPDF using the code 
c
c-----------------
c          DELPDF=0.
c          do k=1,npar
c            DELPDF=DELPDF+dpdfs(i,k)**2
c          end do
c-----------------
c     and its random value can be stored in RPDF using the code 
c-----------------
c          RPDF=pdfs(i)          
c          do k=1,npar
c            s=0.
c            do k=1,96
c              s=s+(2*rndm(xxx)-1)/sqrt(32.)
c            end do
c            RPDF=RPDF+s*dpdfs(i,k)
c          end do
c-----------------
c          
c         Reference: Phys. Rev. D67, 14002 (2003) [hep-ph/0211096]
c      
c         Comments: alekhin@sirius.ihep.su                      
c                                                               
c     Initial version: Nov 2002      
c     Revision of May 2003: interpolation scheme was simplified without 
c                           loosing of the quality; value of \alpha_s(Q)
c                           is now stored in PDFS(0), the error in 
c                           \alpha_s(Q) and their correlation with the   
c                           PDFs are included in the derivatives stored in 
c                           DPDFS; the grid is expanded to the lower and
c                           higher values of Q.
c    Revision of June 2003: minor change in the calling sequence (xb and Q2
c                           are now left unchanged even if they are out of the
c                           range allowed by the parameterization).   
c    Revision of Sep 2003:  the grid spacing is changed; parabolic 
c                           interpolation has been implemented; the grid was 
c                           expanded down to Q^2=0.8 Gev^2.

      implicit real*8(a-h,o-z)
      parameter(nxb=99,nq=20,np=9,nvar=15)
      real*4 f(0:np,nxb,nq+1)
      real*8 pdfs(0:np),dpdfs(0:np,nvar)
      real*4 df(nvar,0:np,nxb,nq+1)

      character locdir*2
      character pdford*1
      dimension pdford(3)
      character pdfschem*3
      dimension pdfschem(0:1)
      character pdfset*3
      dimension pdfset(0:3)

      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.8d0,2d8/
      data KORDS,KSCHEMS,KSETS /-1,-1,-1/

c I/O channel to read the data
      data nport/1/
c put in your local address of the PDFs files in LOCDIR
      data locdir /'./'/
c                   1234567890123456789012345678901234567890123456

      data pdford/'1','2','3'/
      data pdfschem /'ffn','vfn'/
      data pdfset /'   ','_mc','_ss','_kr'/

      save kords,kschems,ksets,f,df,dels,delx,x1,delx1,xlog1,nxbb

      if (kschem.eq.0) then 
        npdf=6
      else 
        npdf=9
      end if
      npar=nvar

      if(kords.eq.kord.and.kschems.eq.kschem.and.ksets.eq.kset) goto 10

      kords=kord
      kschems=kschem
      ksets=kset      

      dels=(log(log(qsqmax/0.04))-log(log(qsqmin/0.04)))/(nq-1)

      nxbb=nxb/2
      x1=0.3
      xlog1=log(x1)
      delx=(log(x1)-log(xmin))/(nxbb-1)
      DELX1=(1-x1)**2/(nxbb+1)

      open(unit=nport,status='old',err=199
     ,    ,file=locdir//'a02.pdfs_'//pdford(kord)//'_'
     /     //pdfschem(kschem)//pdfset(kset))
      do n=1,nxb-1
        do m=1,nq
          read(nport,100) (f(i,n,m),i=0,npdf)
        end do
      end do
      close(unit=nport)
  100 format (13f11.5)

      open(unit=nport,status='old'
     ,    ,file=locdir//'a02.dpdfs_'//pdford(kord)//'_'
     /                              //pdfschem(kschem))
      do n=1,nxb-1
        do m=1,nq
          do i=0,npdf 
            read (nport,*) (df(k,i,n,m),k=1,npar)
          end do
        end do
      end do
      close(unit=nport)

      do i=1,npdf
        do m=1,nq
          f(i,nxb,m)=0d0
          do k=1,npar
            df(k,i,nxb,m)=0d0
          end do
        end do
      end do
      do m=1,nq
        f(0,nxb,m)=f(0,nxb-1,m)
        do k=1,npar
          df(k,0,nxb,m)=df(k,0,nxb-1,m)
         end do
      end do

  10  continue

      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,qsq
      if(xb.lt.xmin.or.xb.gt.xmax)       print 98,x
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ')

      x=max(xb,xmin)
      x=min(xb,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)

      if (x.gt.x1) then
        xd=(1-x1)**2-(1-x)**2
        n=int(xd/delx1)+nxbb
        a=xd/delx1-n+nxbb
      else
        xd=log(x)-xlog1
        n=nxbb+int(xd/DELX)-1
        a=xd/delx-n+nxbb
      end if

      ss=log(log(qsq/0.04))-log(log(qsqmin/0.04))
      m=int(ss/dels)+1
      b=ss/dels-m+1

      do i=0,npdf
        if (n.gt.1.and.m.gt.1.and.n.ne.49) then 
        pdfs(i)= f(i,n,m)*(1+a*b-a**2-b**2) + f(i,n+1,m+1)*a*b
     +    +      f(i,n+1,m)*a*(a-2*b+1)/2. + f(i,n,m+1)*b*(b-2*a+1)/2.
     +    +      f(i,n-1,m)*a*(a-1)/2. + f(i,n,m-1)*b*(b-1)/2.
        else
        pdfs(i)= (1-a)*(1-b)*f(i,n,m) + (1-a)*b*f(i,n,m+1)
     .    +       a*(1-b)*f(i,n+1,m) + a*b*f(i,n+1,m+1)
        end if
        do k=1,npar
        if (n.gt.1.and.m.gt.1.and.n.ne.49) then 
        dpdfs(i,k)= df(k,i,n,m)*(1+a*b-a**2-b**2) + df(k,i,n+1,m+1)*a*b
     +    +  df(k,i,n+1,m)*a*(a-2*b+1)/2. + df(k,i,n,m+1)*b*(b-2*a+1)/2.
     +    +  df(k,i,n-1,m)*a*(a-1)/2. + df(k,i,n,m-1)*b*(b-1)/2.
        else 
        dpdfs(i,k)=(1-a)*(1-b)*df(k,i,n,m)+(1-a)*b*df(k,i,n,m+1)
     .  +       a*(1-b)*df(k,i,n+1,m) + a*b*df(k,i,n+1,m+1)
        end if
        end do
      end do

      return

 199  print *,'The PDF set is inavailable (FILE:'
     ,     ,locdir//'a02.pdfs_'//pdford(kord)//'_'
     /     //pdfschem(kschem)//pdfset(kset),')'
      return

      end





      FUNCTION Ctq1Pd (Iset, Iparton, X, Q, Irt)
C                                                   -=-=- ctq1pd

C===========================================================================
C GroupName: Cteq
C Description: CTEQ parton distributions
C ListOfFiles: ctq1pd ctq1pdf ctq2df ctq2pd ctq2pds ctq3df ctq3pd ctq3pds ctq4pdf readtbl partonx setctq4
C===========================================================================
C #Header: /Net/d2a/wkt/1hep/2pdf/prz/RCS/Cteq.f,v 6.3 98/08/14 09:24:13 wkt Exp $
C #Log:	Cteq.f,v $
c Revision 6.3  98/08/14  09:24:13  wkt
c setctq4 improved
c 
c Revision 6.1  97/11/15  18:05:08  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c 
c Revision 5.0  94/11/08  22:40:36  wkt
c Start of a new 5.* version; Files reorganized according to logic.
c 
c Revision 1.2  93/02/26  10:42:43  wkt
c Version with heavy quark threshold factor and faster algorithm.
c 
c Revision 1.1  93/02/14  17:30:21  botts
c The new Faster version.
c Revision 1.0  93/02/08  18:35:25  wkt
c Initial revision

C========================================================================

C     CTEQ distribution function in a parametrized form.  

C     (No data tables are needed.)

C The returned function value is the PROBABILITY density for a given FLAVOR.

C  !! A companion function (next module), which this one depends on, 
C  !!        Ctq1Pdf (Iset, Iparton, X, Q, Irt)
C  !! gives the VALENCE and SEA MOMENTUM FRACTION distributions. 

C  \\  A parallel (independent) program CtqPds (not included in this file) 
C  ||  in Subroutine form is also available. 
C  ||  It returns ALL the parton flavors at once in an array form.
C  //  See details in that separate file if you are interested.

C Ref.: "CTEQ Parton Distributions and Flavor Dependence of the Sea Quarks"
C       by: J. Botts, J.G. Morfin, J.F. Owens, J. Qiu, W.K. Tung & H. Weerts
C       MSUHEP-92-27, Fermilab-Pub-92/371, FSU-HEP-92-1225, ISU-NP-92-17
C       Now published in Phys. Lett.B304, 159 (1993).

C     Since this is an initial distribution, and there may be updates, it is 
C     useful for the authors to maintain a record of the distribution list.
C     Please do not freely distribute this program package; instead, refer any 
C     interested colleague to direct their request for a copy to:
C     Botts@msupa.pa.msu.edu  or  Botts@msupa (bitnet)  or  MSUHEP::Botts

C    If you have any questions concerning these distributions, direct inquires 
C    to Jim Botts or Wu-Ki Tung (username Tung at same E-mail nodes as above).
c 
C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C     where Iset (= 1, 2, ..., 5) is the set label; 

C       Name convention for CTEQ distributions:  CTEQnSx  where
C           n : version number                      (currently n = 1)
C           S : factorization scheme label: = [M D L] for [MS-bar DIS LO]  
c               resp.
C           x : special characteristics, if any 
C                    (e.g. S for singular gluon, L for "LEP lambda value")

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)

C X, Q are the usual x, Q; Irt is a return error code (not implemented yet).

C --> Iset = 1, 2, 3, 4, 5 correspond to the following CTEQ global fits:
C     cteq1M, cteq1MS, cteq1ML, cteq1D, cteq1L  respectively.

C --> QCD parameters for parton distribution set Iset can be obtained inside
C         the user's program by:
C     Dum = Prctq1 
C    >        (Iset, Iord, Ischeme, MxFlv,
C    >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
C    >         Xmin, Qini, Qmax, ExpNor)
C     where all but the first argument are output parameters.
C     They should be self-explanary -- see details in next module.

C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-5 < x < 1 ; 2 GeV < Q < 1 TeV.  The  
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   an additional threshold factor is applied which simulate the Q-dependence
C   of the QCD evolution in that region.

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Ifl = Iparton
      JFL = ABS(Ifl)
C                                                             Valence
      IF (Ifl .LE. 0) THEN
        VL = 0
      ELSEIF (Ifl .LE. 2) THEN
        VL = Ctq1Pdf(Iset, Ifl, X, Q, Irt)
      ELSE
        VL = 0
      ENDIF
C                                                                         Sea
      SEA = Ctq1Pdf (Iset, -JFL, X, Q, Irt)
C                                              Full (probability) Distribution 
      Ctq1Pd = (VL + SEA) / X
       
      Return
C                         *************************
      END

      FUNCTION Ctq1Pdf (Iset, Iprtn, XX, QQ, Irt)
C                                                   -=-=- ctq1pdf

C            Returns xf(x,Q) -- the momentum fraction distribution !!
C            Returns valence and sea rather than combined flavor distr.

C            Iset : PDF set label

C            Iprtn  : Parton label:   2, 1 = d_ and u_ valence
C                                     0 = gluon
C                            -1, ... -6 = u, d, s, c, b, t sea quarks

C            XX  : Bjorken-x
C            QQ  : scale parameter "Q"
C            Irt : Return code

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (Nex = 5, MxFl = 6, Npn = 3, Nst = 30, Nexpt=20)
      Parameter (Nst4 = Nst*4)

      DIMENSION
     >   Iord(Nst), Isch(Nst), Nqrk(Nst),Alm(Nst)
     > , Vlm(4:6,Nst), Qms(4:6, Nst)
     > , Xmn(Nst), Qmn(Nst), Qmx(Nst), Nexp(Nexpt)
     > , Mex(Nst), Mpn(Nst), ExpN(Nexpt, Nst), ExpNor(Nexpt)

c                                                             CTEQ1M
      DATA 
     >  Isch(1), Iord(1), Nqrk(1), Alm(1) /  1,  2,  6,  .152 / 
     >  (Vlm(I,1), I=4,6) / .231,    .152,    .059  /
     >  (Qms(I,1), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(1), Qmn(1), Qmx(1) /  1.E-5,  2.00,  1.E3  /
     >  Mex(1), Mpn(1), Nexp(1) /  5, 3, 8  /
     >  (ExpN(I, 1), I=1,8)
     >  / 0.989, 1.00, 1.02, 0.978, 1.10, 0.972, 0.987, 0.846 /
c                                                             CTEQ1MS
      DATA 
     >  Isch(2), Iord(2), Nqrk(2), Alm(2) /  1,  2,  6, .152  / 
     >  (Vlm(I,2), I=4,6) / .231,    .152,    .059  /
     >  (Qms(I,2), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(2), Qmn(2), Qmx(2) /  1.E-5,  2.00,  1.E3  /
     >  Mex(2), Mpn(2), Nexp(2) /  5, 3, 8  /
     >  (ExpN(I, 2), I=1,8 )
     >  / 0.989, 1.00, 1.02, 0.984, 1.05, 0.891, 0.923, 0.824 /
c                                                             CTEQ1ML
      DATA 
     >  Isch(3), Iord(3), Nqrk(3), Alm(3) /  1,  2,  6, .220  / 
     >  (Vlm(I,3), I=4,6) / .322,    .220,     .088  /
     >  (Qms(I,3), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(3), Qmn(3), Qmx(3) /  1.E-5,  2.00,  1.E3  /
     >  Mex(3), Mpn(3), Nexp(3) /  5, 3, 8 /
     >  (ExpN(I, 3), I=1,8 )
     >  / 0.985, 1.00, 1.01, 0.977, 1.07, 1.31, 1.19, 1.09 /

c                                                             CTEQ1D
      DATA 
     >  Isch(4), Iord(4), Nqrk(4), Alm(4) /  2,  2,  6, .164  / 
     >  (Vlm(I,4), I=4,6) / .247,    .164,    .064  /
     >  (Qms(I,4), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(4), Qmn(4), Qmx(4) /  1.E-5,  2.00,  1.E3  /
     >  Mex(4), Mpn(4), Nexp(4) /  5, 3, 8 /
     >  (ExpN(I, 4), I=1,8 )
     >  / 0.983, 1.00, 1.01, 0.975, 0.964, 1.23, 1.00, 1.12 /
c                                                             CTEQ1L
      DATA 
     >  Isch(5), Iord(5), Nqrk(5), Alm(5) /  1,  1,  6, .125  / 
     >  (Vlm(I,5), I=4,6) / .168,    .125,     .063   /
     >  (Qms(I,5), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(5), Qmn(5), Qmx(5) /  1.E-5,  2.00,  1.E3  /
     >  Mex(5), Mpn(5), Nexp(5) /  5, 3, 8  /
     >  (ExpN(I, 5), I=1,8 )
     >  / 0.982, 1.01, 1.00, 0.972, 0.840, 0.959, 0.930, 0.861 /

      Data ist, lp, qsto, Aln2 / 0, -10, 1.2345, 0.6931 /

      X  = XX
      if(iset.eq.ist.and.iprtn.eq.lp.and.qsto.eq.qq) goto 100

      Irt = 0
      ip = abs(iprtn)
      If (Ip.GE.5.and.QQ.LE.Qms(Ip, Iset)) Then
         Ctq1Pdf = 0.0
         Return
      Endif

c if heavy parton, different logarithmic scale

      if(ip.ge.5) then
         Qi = qms(ip,iset)
      else
         Qi   = Qmn (Iset)
      endif

      Alam = Alm (Iset)
      sta = log(qq/alam)
      stb = log(qi/alam)

c      SBL = LOG(QQ/Alam) / LOG(Qi/Alam)
      sbl = sta/stb
      SB = LOG (SBL)
      SB2 = SB*SB
      SB3 = SB2*SB

      iflv = 3 - iprtn

      Goto (1, 2, 3, 4, 5), Iset

 1    Goto(11,12,13,14,15,16,17,18,19)Iflv    
c   ifl =     2
 11   A0=0.3636E+01*(1.0 + 0.3122E+00*SB+0.1396E+00*SB2+0.4251E+00*SB3)
      A1=0.6930E+00-.2574E-01*SB+0.1047E+00*SB2-.2794E-01*SB3
      A2=0.3195E+01+0.4045E+00*SB-.3737E+00*SB2-.1677E+00*SB3
      A3=0.1009E+00*(1.0 -.1784E+01*SB+0.6263E+00*SB2+0.7337E-01*SB3)
     $  -1.0
      A4=0.2910E+00-.2793E+00*SB+0.6155E-01*SB2+0.5150E-02*SB3
      A5=0.0000E+00+0.3185E+00*SB+0.1953E+00*SB2+0.4184E-01*SB3
      goto 100
c   ifl =     1
 12   A0=0.2851E+00*(1.0 + 0.3617E+00*SB-.4526E+00*SB2+0.5787E-01*SB3)
      A1=0.2690E+00+0.1104E-01*SB+0.1888E-01*SB2-.1031E-01*SB3
      A2=0.3766E+01+0.7850E+00*SB-.3053E+00*SB2+0.1822E+00*SB3
      A3=0.2865E+02*(1.0 -.9774E+00*SB+0.5958E+00*SB2-.1234E+00*SB3)
     $  -1.0
      A4=0.8230E+00-.3612E+00*SB+0.5520E-01*SB2+0.1571E-01*SB3
      A5=0.0000E+00+0.2145E-01*SB+0.2289E+00*SB2-.4947E-01*SB3
      goto 100
c   ifl =     0
 13   A0=0.2716E+01*(1.0 -.2092E+01*SB+0.1500E+01*SB2-.3703E+00*SB3)
      A1=-.3100E-01-.7963E+00*SB+0.1129E+01*SB2-.4191E+00*SB3
      A2=0.8015E+01+0.1168E+01*SB-.1625E+01*SB2-.1130E+01*SB3
      A3=0.4813E+02*(1.0 -.4951E+00*SB-.8715E+00*SB2+0.5893E+00*SB3)
     $  -1.0
      A4=0.2773E+01-.6329E+00*SB-.1048E+01*SB2+0.1418E+00*SB3
      A5=0.0000E+00+0.5048E+00*SB+0.2390E+01*SB2-.4159E+00*SB3
      goto 100
c   ifl =    -1
 14   A0=0.3085E+00*(1.0 + 0.9422E+00*SB-.2606E+01*SB2+0.1364E+01*SB3)
      A1=0.5000E-02-.6433E+00*SB+0.4980E+00*SB2-.1780E+00*SB3
      A2=0.7490E+01+0.9112E+00*SB-.2047E+01*SB2+0.1456E+01*SB3
      A3=0.1145E-01*(1.0 + 0.4610E+01*SB+0.1699E+01*SB2+0.1296E+00*SB3)
     $  -1.0
      A4=0.6030E+00-.8081E+00*SB+0.9410E+00*SB2-.4458E+00*SB3
      A5=0.0000E+00-.1736E+01*SB+0.2863E+01*SB2-.1268E+01*SB3
      goto 100
c   ifl =    -2
 15   A0=0.1324E+00*(1.0 -.1050E+01*SB+0.4844E+00*SB2-.1043E+00*SB3)
      A1=-.1580E+00+0.1672E+00*SB-.4100E+00*SB2+0.1793E+00*SB3
      A2=0.8559E+01-.7351E-01*SB+0.5898E+00*SB2-.2655E+00*SB3
      A3=0.2378E+02*(1.0 -.1108E+00*SB-.1646E-01*SB2+0.1129E-01*SB3)
     $  -1.0
      A4=0.1477E+01+0.3312E-01*SB-.2191E+00*SB2+0.9588E-01*SB3
      A5=0.0000E+00+0.1850E+01*SB-.1481E+01*SB2+0.6222E+00*SB3
      goto 100
c   ifl =    -3
 16   A0=0.3208E+00*(1.0 -.4755E+00*SB-.4003E+00*SB2+0.2300E+00*SB3)
      A1=-.3200E-01-.3357E+00*SB+0.3222E-01*SB2+0.5011E-01*SB3
      A2=0.1164E+02+0.1048E+01*SB-.1097E+01*SB2-.4431E+00*SB3
      A3=0.5065E+02*(1.0 + 0.2484E+00*SB-.9235E+00*SB2+0.1935E+00*SB3)
     $  -1.0
      A4=0.3300E+01-.6785E+00*SB+0.5337E+00*SB2-.4035E+00*SB3
      A5=0.0000E+00-.2496E+00*SB+0.3903E+00*SB2+0.1392E+00*SB3
      goto 100
c   ifl =    -4
 17   A0=0.7967E-06*(1.0 + 0.1587E+01*SB+0.1812E+02*SB2-.1333E+02*SB3)
     $ *sqrt(sta - stb)
      A1=0.1096E+01-.1236E+01*SB+0.1014E+02*SB2+0.1940E+01*SB3
      A2=0.4366E+00+0.1197E+02*SB-.5471E+00*SB2-.5427E+01*SB3
      A3=0.4650E+03*(1.0 + 0.1310E+02*SB-.1918E+02*SB2+0.6791E+01*SB3)
     $  -1.0
      A4=-.8486E+00+0.7457E+00*SB-.1083E+02*SB2-.1210E+01*SB3
      A5=0.3494E+01-.3511E+01*SB-.1766E+01*SB2+0.3442E+01*SB3
      goto 100
c   ifl =    -5
 18   A0=0.1713E-03*(1.0 + 0.2562E+02*SB-.2988E+02*SB2+0.4798E+01*SB3)
     $ *sqrt(sta - stb)
      A1=-.5276E-01+0.4105E+00*SB-.1079E+01*SB2+0.6278E+00*SB3
      A2=0.4515E+01+0.8369E+01*SB-.1192E+02*SB2+0.3403E+01*SB3
      A3=0.1756E+01*(1.0 + 0.1325E+02*SB-.2997E+02*SB2+0.1758E+02*SB3)
     $  -1.0
      A4=0.3557E-01+0.4159E+01*SB-.6947E+01*SB2+0.2982E+01*SB3
      A5=0.2551E+01+0.2168E+01*SB-.5119E+01*SB2+0.3739E+01*SB3
      goto 100
c   ifl =    -6
 19   A0=0.7510E-04*(1.0 + 0.2836E+02*SB-.3000E+02*SB2-.2979E+02*SB3)
     $ *sqrt(sta - stb)
      A1=-.1855E+00+0.4543E+00*SB-.1448E+01*SB2+0.2009E-01*SB3
      A2=0.6775E+01-.4210E+01*SB-.1221E+01*SB2+0.1199E+02*SB3
      A3=0.1070E+01*(1.0 + 0.8356E+01*SB-.2992E+02*SB2+0.2433E+02*SB3)
     $  -1.0
      A4=-.4601E-01+0.4248E+01*SB-.1736E+01*SB2+0.1187E+02*SB3
      A5=0.2771E+01+0.1382E+01*SB-.4797E+01*SB2+0.1273E+01*SB3
      goto 100


 2    Goto(21,22,23,24,25,26,27,28,29)Iflv    
c                                                             CTEQ1MS
c   ifl =     2
 21   A0=0.1828E+01*(1.0 -.8698E+00*SB+0.2906E+00*SB2-.2003E-01*SB3)
      A1=0.6060E+00+0.8595E-01*SB-.4934E-01*SB2+0.2221E-01*SB3
      A2=0.3454E+01-.3115E+00*SB+0.1321E+01*SB2-.3490E+00*SB3
      A3=0.2616E+00*(1.0 -.1670E+01*SB+0.2333E+01*SB2+0.7730E-01*SB3)
     $  -1.0
      A4=0.8920E+00-.8500E-02*SB+0.4960E+00*SB2-.4045E-01*SB3
      A5=0.0000E+00+0.1091E+01*SB-.1613E+00*SB2+0.3773E-01*SB3
      goto 100
c   ifl =     1
 22   A0=0.2885E+00*(1.0 + 0.3388E+00*SB-.4550E+00*SB2+0.6005E-01*SB3)
      A1=0.2730E+00+0.1198E-01*SB+0.1880E-01*SB2-.1077E-01*SB3
      A2=0.3736E+01+0.7687E+00*SB-.2731E+00*SB2+0.1638E+00*SB3
      A3=0.2741E+02*(1.0 -.9585E+00*SB+0.5925E+00*SB2-.1239E+00*SB3)
     $  -1.0
      A4=0.8040E+00-.3546E+00*SB+0.6123E-01*SB2+0.1086E-01*SB3
      A5=0.0000E+00+0.4277E-01*SB+0.2187E+00*SB2-.4646E-01*SB3
      goto 100
c   ifl =     0
 23   A0=0.8416E-01*(1.0 -.1996E+01*SB+0.1903E+01*SB2-.6722E+00*SB3)
      A1=-.4790E+00-.5459E+00*SB+0.1638E+01*SB2-.4342E+00*SB3
      A2=0.5071E+01+0.1470E+01*SB-.2401E+01*SB2+0.1273E+01*SB3
      A3=0.2847E+02*(1.0 + 0.1124E+00*SB-.1338E+01*SB2+0.7115E+00*SB3)
     $  -1.0
      A4=0.4990E+00-.7208E+00*SB+0.3333E-03*SB2-.2354E+00*SB3
      A5=0.0000E+00-.4480E+00*SB+0.3720E+01*SB2-.1838E+01*SB3
      goto 100
c   ifl =    -1
 24   A0=0.4378E+00*(1.0 -.1244E+01*SB+0.3278E+01*SB2-.2098E+01*SB3)
      A1=0.3500E-01-.1298E+01*SB+0.1229E+01*SB2-.3665E+00*SB3
      A2=0.6781E+01+0.4078E+01*SB-.9711E+00*SB2-.1536E+01*SB3
      A3=0.1527E-03*(1.0 + 0.1430E+02*SB+0.3000E+02*SB2+0.2771E+02*SB3)
     $  -1.0
      A4=0.3060E+00+0.1011E+01*SB-.2045E+01*SB2+0.9422E+00*SB3
      A5=0.0000E+00-.3205E+01*SB+0.2683E+01*SB2-.1746E+00*SB3
      goto 100
c   ifl =    -2
 25   A0=0.7413E-01*(1.0 + 0.1291E+01*SB-.2667E+01*SB2+0.1076E+01*SB3)
      A1=-.2730E+00-.1206E+00*SB+0.1828E+00*SB2-.1001E+00*SB3
      A2=0.7719E+01+0.1537E+01*SB-.6410E+00*SB2-.3920E-01*SB3
      A3=0.1799E+02*(1.0 -.1334E+01*SB+0.1916E+01*SB2-.8878E+00*SB3)
     $  -1.0
      A4=0.1167E+01-.9176E-01*SB+0.5132E+00*SB2-.3460E+00*SB3
      A5=0.0000E+00-.5023E+00*SB+0.1951E+01*SB2-.8427E+00*SB3
      goto 100
c   ifl =    -3
 26   A0=0.6551E+00*(1.0 -.5968E-01*SB+0.5621E-02*SB2-.2074E+00*SB3)
      A1=0.2800E-01-.1138E+01*SB+0.1178E+01*SB2-.4425E+00*SB3
      A2=0.7553E+01+0.3996E+01*SB-.4448E+01*SB2+0.1673E+01*SB3
      A3=0.9264E-01*(1.0 -.1760E+01*SB+0.1634E+01*SB2-.4067E+00*SB3)
     $  -1.0
      A4=0.1970E+00+0.5256E+00*SB-.9775E+00*SB2+0.4488E+00*SB3
      A5=0.0000E+00-.3668E+01*SB+0.4757E+01*SB2-.1717E+01*SB3
      goto 100
c   ifl =    -4
 27   A0=0.1486E-03*(1.0 + 0.2107E+01*SB-.1056E+02*SB2+0.1403E+02*SB3)
     $ * sqrt(sta - stb)
      A1=0.2115E+00-.1702E+01*SB+0.2571E+01*SB2-.1177E+01*SB3
      A2=0.3533E+01+0.1367E+01*SB-.3397E+01*SB2+0.6260E+01*SB3
      A3=0.1096E+02*(1.0 + 0.9213E+01*SB-.2020E+02*SB2+0.1084E+02*SB3)
     $  -1.0
      A4=0.7041E+00-.7236E+00*SB+0.2766E-01*SB2+0.7352E+00*SB3
      A5=0.3904E+01-.4398E+01*SB+0.7056E+01*SB2-.3722E+01*SB3
      goto 100
c   ifl =    -5
 28   A0=0.1201E-03*(1.0 + 0.5408E+01*SB-.1489E+02*SB2+0.1667E+02*SB3)
     $  * sqrt(sta - stb)
      A1=0.1420E-01-.1525E+01*SB+0.2408E+01*SB2-.1154E+01*SB3
      A2=0.4254E+01+0.2836E+01*SB-.6018E+00*SB2+0.4133E+00*SB3
      A3=0.5696E+01*(1.0 + 0.9451E+01*SB-.2029E+02*SB2+0.1033E+02*SB3)
     $  -1.0
      A4=0.4775E+00-.6695E+00*SB+0.2747E+00*SB2-.1051E+00*SB3
      A5=0.3330E+01-.5133E+01*SB+0.6921E+01*SB2-.3283E+01*SB3
      goto 100
c   ifl =    -6
 29   A0=0.7697E-04*(1.0 + 0.2801E+02*SB-.1901E+02*SB2-.2880E+02*SB3)
     $ *sqrt(sta - stb)
      A1=-.2249E+00+0.4432E+00*SB-.1454E+01*SB2+0.3509E-01*SB3
      A2=0.6642E+01-.2702E+01*SB+0.8229E+01*SB2+0.8243E+01*SB3
      A3=0.1146E+01*(1.0 + 0.8104E+01*SB-.2998E+02*SB2+0.2812E+02*SB3)
     $  -1.0
      A4=-.6421E-01+0.4246E+01*SB-.2908E+01*SB2+0.9686E-02*SB3
      A5=0.2606E+01+0.1261E+01*SB-.4933E+01*SB2+0.3476E+00*SB3
      goto 100


 3    Goto(31,32,33,34,35,36,37,38,39)Iflv    
c                                                             CTEQ1ML
c   ifl =     2
 31   A0=0.3777E+01*(1.0 + 0.6986E+00*SB-.20655E+01*SB2+.10334E+01*SB3)
      A1=0.7100E+00+.2880E-01*SB-.7930E-01*SB2+0.5600E-01*SB3
      A2=0.3259E+01+0.1508E+01*SB-.3932E+01*SB2+0.20613E+01*SB3
      A3=0.1304E+00*(1.0 -.2016E+00*SB-.30015E+01*SB2+0.19118E+01*SB3)
     $     -1.0
      A4=0.2890E+00-0.4311E+00*SB+0.7387E+00*SB2-.3697E+00*SB3
      A5=0.0000E+00+0.4320E+00*SB+0.2449E+00*SB2-0.6670E-01*SB3
      goto 100
c   ifl =     1
 32   A0=0.2780E+00*(1.0 + 0.4355E+00*SB-0.4584E+00*SB2+0.4390E-01*SB3)
      A1=0.2760E+00+0.1420E-01*SB+0.1480E-01*SB2-.9800E-02*SB3
      A2=0.3710E+01+0.8250E+00*SB-.3581E+00*SB2+0.1978E+00*SB3
      A3=0.2928E+02*(1.0 -.10154E+01*SB+0.6037E+00*SB2-.1175E+00*SB3)
     $     -1.0
      A4=0.8070E+00-.3575E+00*SB+0.4920E-01*SB2+0.1584E-01*SB3
      A5=0.0000E+00+0.1860E-01*SB+0.2080E+00*SB2-.450E-01*SB3
      goto 100
c   ifl =     0
 33   A0=0.2924E+01*(1.0 -.18916E+01*SB+0.1191E+01*SB2-.2492E+00*SB3)
      A1=0.0000E+00-.9167E+00*SB+0.11147E+01*SB2-.3329E+00*SB3
      A2=0.8529E+01+0.7080E+00*SB-.11345E+01*SB2-.10563E+01*SB3
      A3=0.1420E+03*(1.0 -.15346E+01*SB+0.7261E+00*SB2-.5730E-01*SB3)
     $     -1.0
      A4=0.3396E+01-.11541E+01*SB-.8834E+00*SB2+0.2430E+00*SB3
      A5=0.0000E+00+0.1645E+00*SB+0.19041E+01*SB2+0.1474E+00*SB3
      goto 100
c   ifl =    -1
 34   A0=0.3471E+00*(1.0- 0.1753E+00*SB-.9189E+00*SB2+0.6211E+00*SB3)
      A1=0.1900E-01-.4579E+00*SB+0.2112E+00*SB2-.6180E-01*SB3
      A2=0.7301E+01-.17308E+01*SB+.13666E+01*SB2-.6400E-02*SB3
      A3=0.1853E-04*(1.0 -.18260E+02*SB-.2872E+02*SB2-.23456E+02*SB3)
     $     -1.0
      A4=0.4400E+00-.4672E+00*SB+0.6532E+00*SB2-.3222E+00*SB3
      A5=0.0000E+00-.4679E+00*SB+0.10741E+01*SB2-.5663E+00*SB3
      goto 100
c   ifl =    -2
 35   A0=0.1702E+00*(1.0 -.1041E+01*SB+0.4064E+00*SB2-.5888E-01*SB3)
      A1=-.9300E-01-.4742E-01*SB-.1959E+00*SB2+0.1039E+00*SB3
      A2=0.9119E+01-.7331E-01*SB+0.3506E+00*SB2-.2081E+00*SB3
      A3=0.2981E+02*(1.0 -.1912E+00*SB-.8947E-02*SB2+0.8805E-02*SB3)
     $     -1.0
      A4=0.1668E+01-.6678E-02*SB-.2894E+00*SB2+0.1221E+00*SB3
      A5=0.0000E+00+0.1245E+01*SB-.7843E+00*SB2+0.3724E+00*SB3
      goto 100
c   ifl =    -3
 36   A0=0.3910E+00*(1.0 -.1103E+01*SB+0.5383E+00*SB2-.1083E+00*SB3)
      A1=-.1400E-01-.2471E+00*SB-.8042E-01*SB2+0.7193E-01*SB3
      A2=0.9812E+01-.4860E+01*SB+0.5958E+01*SB2-.2342E+01*SB3
      A3=0.3749E+00*(1.0 -.3569E+01*SB+0.5456E+01*SB2-.2344E+01*SB3)
     $     -1.0
      A4=0.4940E+00+0.2772E+00*SB-.2732E+00*SB2+0.6466E-01*SB3
      A5=0.0000E+00+0.3927E+00*SB-.3216E+00*SB2+0.2164E+00*SB3
      goto 100
c   ifl =    -4
 37   A0=0.3815E-02*(1.0 + 0.2039E+02*SB-.2834E+02*SB2+0.1070E+02*SB3)
     $ * sqrt(sta - stb)
      A1=-.2789E-01-.7345E-03*SB-.3251E+00*SB2+0.1946E+00*SB3
      A2=0.3223E+01-.4268E+00*SB+0.4387E+01*SB2-.2401E+01*SB3
      A3=0.3338E-01*(1.0 -.1163E+02*SB+0.2995E+02*SB2-.1471E+02*SB3)
     $     -1.0
      A4=0.3646E+00-.5767E+00*SB+0.6088E+00*SB2-.2514E+00*SB3
      A5=0.1200E+01+0.2178E+00*SB-.4230E+00*SB2+0.4739E+00*SB3
      goto 100
c   ifl =    -5
 38   A0=0.1666E-02*(1.0 + 0.9518E+01*SB-.4715E+01*SB2-.1060E+01*SB3)
     $ * sqrt(sta - stb)
      A1=-.1231E+00+0.1656E+00*SB-.5219E+00*SB2+0.2750E+00*SB3
      A2=0.3693E+01+0.4922E+01*SB-.1200E+02*SB2+0.7929E+01*SB3
      A3=0.1778E+00*(1.0 + 0.3036E+01*SB-.1184E+02*SB2+0.7940E+01*SB3)
     $     -1.0
      A4=0.5353E+00-.1401E+01*SB+0.1970E+01*SB2-.9405E+00*SB3
      A5=0.1590E+01+0.1025E+01*SB-.2318E+01*SB2+0.1380E+01*SB3
      goto 100
c   ifl =    -6
 39   A0=0.4319E-03*(1.0 + 0.1100E+02*SB-.9520E+00*SB2+0.1434E+02*SB3)
     $ * sqrt(sta - stb)
      A1=-.2512E+00+0.3554E+00*SB-.4120E+00*SB2+0.1328E+00*SB3
      A2=0.4764E+01-.3513E+00*SB+0.1199E+02*SB2-.8290E+01*SB3
      A3=0.8458E-01*(1.0 + 0.2618E+01*SB+0.4407E+01*SB2+0.2991E+02*SB3)
     $    -1.0
      A4=0.3991E+00-.1363E+01*SB+0.1526E+01*SB2-.3179E+01*SB3
      A5=0.1981E+01+0.1496E+01*SB-.1501E+01*SB2+0.3880E+01*SB3
      goto 100

 4    Goto(41,42,43,44,45,46,47,48,49)Iflv
c                                                             CTEQ1D
c   ifl =     2
 41   A0=0.1634E+01*(1.0 -.8336E+00*SB+0.1640E+00*SB2+0.1530E+00*SB3)
      A1=0.5790E+00+0.8587E-01*SB-.6087E-01*SB2+0.1361E-01*SB3
      A2=0.2839E+01+0.3720E+00*SB+0.5264E+00*SB2+0.3538E-01*SB3
      A3=0.1095E+00*(1.0 -.4830E+00*SB+0.3708E+01*SB2-.6165E+00*SB3)
     $  -1.0
      A4=0.8010E+00-.1432E+00*SB+0.1442E+01*SB2-.1286E+01*SB3
      A5=0.0000E+00+0.1035E+01*SB-.5910E-01*SB2-.1982E+00*SB3
      goto 100
c   ifl =     1
 42   A0=0.3535E+00*(1.0 + 0.4352E+00*SB-.2095E+00*SB2-.8455E-02*SB3)
      A1=0.2660E+00-.4096E-03*SB+0.1502E-01*SB2-.1163E-01*SB3
      A2=0.3514E+01+0.8219E+00*SB-.2330E+00*SB2+0.1055E+00*SB3
      A3=0.2200E+02*(1.0 -.9716E+00*SB+0.4552E+00*SB2-.8202E-01*SB3)
     $  -1.0
      A4=0.9000E+00-.3207E+00*SB-.4808E-01*SB2+0.3492E-01*SB3
      A5=0.0000E+00-.6273E-01*SB+0.1497E+00*SB2-.5683E-01*SB3
      goto 100
c   ifl =     0
 43   A0=0.2743E+01*(1.0 -.2027E+01*SB+0.1517E+01*SB2-.4145E+00*SB3)
      A1=0.7000E-02-.9431E+00*SB+0.1231E+01*SB2-.4834E+00*SB3
      A2=0.8200E+01+0.1827E+01*SB-.3453E+01*SB2+0.6763E+00*SB3
      A3=0.4975E+02*(1.0 -.2329E+00*SB-.1245E+01*SB2+0.7194E+00*SB3)
     $  -1.0
      A4=0.2387E+01-.4077E+00*SB-.5542E+00*SB2-.9677E-02*SB3
      A5=0.0000E+00+0.2702E+00*SB+0.2389E+01*SB2-.8274E+00*SB3
      goto 100
c   ifl =    -1
 44   A0=0.2015E+00*(1.0 -.2133E+00*SB-.6770E+00*SB2+0.5011E+00*SB3)
      A1=-.7700E-01-.7104E-01*SB-.3720E+00*SB2+0.2159E+00*SB3
      A2=0.8008E+01-.2049E+01*SB+0.1800E+01*SB2-.4660E+00*SB3
      A3=0.2923E-05*(1.0 + 0.2327E+02*SB+0.1500E+02*SB2+0.2633E+02*SB3)
     $  -1.0
      A4=0.9020E+00-.9191E+00*SB+0.1104E+01*SB2-.5863E+00*SB3
      A5=0.0000E+00+0.5840E+00*SB-.8720E+00*SB2+0.4234E+00*SB3
      goto 100
c   ifl =    -2
 45   A0=0.9117E-01*(1.0 -.4089E+00*SB-.4361E+00*SB2+0.2512E+00*SB3)
      A1=-.2370E+00+0.2492E+00*SB-.3267E+00*SB2+0.1055E+00*SB3
      A2=0.8447E+01+0.6009E+00*SB+0.1003E+01*SB2-.1287E+01*SB3
      A3=0.3106E+02*(1.0 -.3901E-01*SB+0.1443E+00*SB2-.3433E+00*SB3)
     $  -1.0
      A4=0.1629E+01+0.7855E-01*SB-.1573E+00*SB2-.8595E-01*SB3
      A5=0.0000E+00+0.1558E+01*SB-.6295E+00*SB2+0.1847E+00*SB3
      goto 100
c   ifl =    -3
 46   A0=0.3997E+00*(1.0 -.1046E+01*SB+0.6194E+00*SB2-.1342E+00*SB3)
      A1=0.2000E-02-.2544E+00*SB-.1958E+00*SB2+0.1458E+00*SB3
      A2=0.9613E+01-.3919E+01*SB+0.9573E+01*SB2-.5623E+01*SB3
      A3=0.3620E+00*(1.0 -.1858E+01*SB+0.8312E+01*SB2-.5900E+01*SB3)
     $  -1.0
      A4=0.3840E+00+0.3572E+00*SB-.1191E+01*SB2+0.7310E+00*SB3
      A5=0.0000E+00+0.3351E+00*SB-.7709E+00*SB2+0.4296E+00*SB3
      goto 100
c   ifl =    -4
 47   A0=0.2156E-03*(1.0 + 0.2879E+02*SB-.2310E+02*SB2+0.9812E+01*SB3)
     $ * sqrt(sta - stb)
      A1=0.9086E-01-.1250E+00*SB-.7373E-01*SB2-.2201E-01*SB3
      A2=0.3588E+01+0.4518E+01*SB-.8930E-01*SB2+0.9163E-02*SB3
      A3=0.5216E+01*(1.0 + 0.5912E+00*SB-.4111E+00*SB2+0.7330E+00*SB3)
     $  -1.0
      A4=0.3145E+00+0.1233E+01*SB-.7478E+00*SB2+0.4657E+00*SB3
      A5=0.2723E+01-.4110E+00*SB+0.4868E-01*SB2-.3075E+00*SB3
      goto 100
c   ifl =    -5
 48   A0=0.7476E-03*(1.0 + 0.1454E+02*SB-.2509E+02*SB2+0.1184E+02*SB3)
     $ * sqrt(sta - stb)
      A1=-.1955E-01-.1712E+00*SB-.1686E+00*SB2+0.2339E+00*SB3
      A2=0.4616E+01-.6859E+00*SB-.3959E+01*SB2+0.5530E+01*SB3
      A3=0.9881E+01*(1.0 -.1239E+02*SB+0.2721E+02*SB2-.1850E+02*SB3)
     $  -1.0
      A4=0.1200E+02-.1133E+02*SB+0.8138E+01*SB2+0.1199E+02*SB3
      A5=0.2226E+01-.5738E+00*SB+0.5239E+00*SB2+0.3825E+00*SB3
      goto 100
c   ifl =    -6
 49   A0=0.8392E-06*(1.0 + 0.1844E+02*SB-.1110E+02*SB2-.2504E+02*SB3)
     $ * sqrt(sta - stb)
      A1=0.2127E+00-.5602E+00*SB+0.4777E+01*SB2-.1014E+02*SB3
      A2=0.1229E+01+0.7495E+01*SB-.5024E+01*SB2-.1200E+02*SB3
      A3=0.2868E+02*(1.0 + 0.7634E+01*SB-.2916E+02*SB2+0.2953E+02*SB3)
     $  -1.0
      A4=0.5970E+00+0.1138E+01*SB-.1439E+01*SB2-.1966E+01*SB3
      A5=0.6429E+01-.6673E+00*SB+0.7008E+01*SB2-.1157E+02*SB3
      goto 100

 5    Goto(51,52,53,54,55,56,57,58,59)Iflv
c                                                             CTEQ1L
c   ifl =     2
 51   A0=  1.791*(1.0 -0.449*SB-0.445*SB2+  0.401*SB3)
      A1=  0.608+  0.069*SB+  0.005*SB2-0.037*SB3
      A2=  3.470-0.375*SB+  2.267*SB2-1.261*SB3
      A3=  0.315*(1.0 -2.628*SB+  6.481*SB2-3.834*SB3)-1.0
      A4=  1.007-0.732*SB+  1.490*SB2-0.966*SB3
      A5=  0.000+  0.741*SB+  0.563*SB2-0.525*SB3
      goto 100
c   ifl =     1
 52   A0=  0.513*(1.0 +   0.032*SB-0.120*SB2+  0.013*SB3)
      A1=  0.276+  0.052*SB+  0.000*SB2-0.006*SB3
      A2=  3.579+  0.763*SB-0.135*SB2+  0.083*SB3
      A3= 17.993*(1.0 -0.725*SB+  0.241*SB2-0.020*SB3)-1.0
      A4=  1.120-0.357*SB+  0.008*SB2+  0.028*SB3
      A5=  0.000+  0.311*SB+  0.029*SB2-0.010*SB3
      goto 100
c   ifl =     0
 53   A0=  2.710*(1.0 -1.773*SB+  0.970*SB2-0.149*SB3)
      A1= -0.010-1.636*SB+  2.087*SB2-0.637*SB3
      A2=  7.174+  2.102*SB-2.209*SB2-0.420*SB3
      A3= 29.904*(1.0 -0.756*SB-0.506*SB2+  0.605*SB3)-1.0
      A4=  2.572-0.437*SB-0.968*SB2+  0.243*SB3
      A5=  0.000-1.776*SB+  4.266*SB2-0.335*SB3
      goto 100
c   ifl =    -1
 54   A0=  0.278*(1.0 - 1.022*SB+  0.6457*SB2-0.1824*SB3)
      A1=  0.0862*SB-0.8657*SB2+  0.4185*SB3
      A2= 11.000-1.2809*SB+ 1.2516*SB2+0.061*SB3
      A3= 37.338*(1.0 - 0.9404*SB+  0.2517*SB2+0.1364*SB3)-1.0
      A4=  1.960-  0.3385*SB-0.3422*SB2+0.3653*SB3
      A5=  0.000+1.424*SB-2.7503*SB2+  1.2226*SB3
      goto 100
c   ifl =    -2
 55   A0=  0.154*(1.0 -0.659*SB+  0.005*SB2+  0.061*SB3)
      A1= -0.128+  0.279*SB-0.786*SB2+  0.363*SB3
      A2=  8.649+  0.071*SB+  0.351*SB2-0.051*SB3
      A3= 43.685*(1.0 -0.603*SB+  0.037*SB2+  0.134*SB3)-1.0
      A4=  2.238-0.338*SB-0.199*SB2+  0.157*SB3
      A5=  0.000+  1.681*SB-2.068*SB2+  0.975*SB3
      goto 100
c   ifl =    -3
 56   A0=  0.372*(1.0 -1.939*SB+  1.504*SB2-0.440*SB3)
      A1=  0.009+  0.610*SB-1.387*SB2+  0.579*SB3
      A2= 10.273-4.833*SB+  6.583*SB2-2.633*SB3
      A3=  0.160*(1.0 +  10.325*SB-2.027*SB2+  1.571*SB3)-1.0
      A4=  0.819-1.660*SB+  1.845*SB2-0.829*SB3
      A5=  0.000+  3.558*SB-3.940*SB2+  1.302*SB3
      goto 100
c   ifl =    -4
 57   A0=  (7.5242E-5)*(1.0+22.0905*SB+7.1209*SB2-8.303*SB3)*
     $     sqrt(sta - stb)
      A1=  0.125-0.3027*SB+0.1564*SB2-0.091*SB3
      A2=  2.0388+1.2161*SB+11.5296*SB2-8.0659*SB3
      A3=  14.849*(1.0 -2.556*SB+3.5268*SB2-1.6353*SB3)-1.0
      A4=  0.3061-0.0901*SB+0.953*SB2-0.4871*SB3
      A5=  2.7352+0.1811*SB-0.5167*SB2+0.0543*SB3
      goto 100
c   ifl =    -5
 58   A0=  (3.751E-4)*(1.0 + 21.5993*SB+3.1379*SB2-18.8328*SB3)*
     $     sqrt(sta - stb)
      A1= -0.0256-0.7717*SB+ 1.1499*SB2-0.5037*SB3
      A2=  4.9241+4.0107*SB-4.7012*SB2+0.1097*SB3
      A3=  2.842*(1.0 -2.2184*SB+  2.0293*SB2-0.6907*SB3)-1.0
      A4=  -0.1352+ 0.8753*SB-1.2626*SB2+  0.667*SB3
      A5=  1.5627-0.4917*SB+ 1.5927*SB2-0.351*SB3
      goto 100
c   ifl =    -6
 59   A0=(2.725E-4)*(1.0 +  18.8497*SB-26.5797*SB2-29.0774*SB3)*
     $     sqrt(sta - stb)
      A1= -0.2204-1.0048*SB+0.9415*SB2-0.4274*SB3
      A2=  11.034-9.8362*SB-11.1034*SB2-9.1977*SB3
      A3=  2.084*(1.0 -2.881*SB+1.2778*SB2-2.9328*SB3)-1.0
      A4= -0.0872+  0.200*SB-1.6187*SB2-1.6058*SB3
      A5=  0.8684+4.7047*SB-1.4614*SB2-5.2309*SB3
      goto 100

 100  Ctq1Pdf = A0*(x**A1)*((1.-x)**A2)*(1.+A3*(x**A4))
     $     *(log(1.+1./x))**A5

      if(Ctq1Pdf.lt.0.0) Ctq1Pdf = 0.0

      Ist = Iset

      Lp  = Iprtn
      Qsto = QQ

      Return
C                                  -----------------------
      ENTRY WLAMBD (ISET, IORDER)

      IORDER = IORD (ISET)
      WLAMBD = ALM  (ISET)

      RETURN
C                                  -----------------------
      Entry PrCtq1 
     >        (Iset, Iordr, Ischeme, MxFlv,
     >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
     >         Xmin, Qini, Qmax, ExpNor)

C                           Return QCD parameters and Fitting parameters
C                           associated with parton distribution set Iset.
C    Iord    : Order Of Fit
C    Ischeme : (0, 1, 2)  for  (LO, MS-bar-NLO, DIS-NLO) resp.
C    MxFlv   : Maximum number of flavors included
C    Alam_i  : i = 4,5,6  Effective lambda for i-flavors 

C    Amas_i  : i = 4,5,6  Mass parameter for flavor i
C    Xmin, Qini, Qmax : self explanary
C    ExpNor(I) : Normalization factor for the experimental data set used in
C                obtaining the best global fit for parton distributions Iset:
C     I = 1,     2,      3,     4,     5,     6,     7,     8
C      BCDMS   NMC90  NMC280  CCFR   E605    WA70   E706   UA6

      Iordr  = Iord (Iset)
      Ischeme= Isch (Iset)
      MxFlv  = Nqrk (Iset)

      Alam4  = Vlm(4,Iset)
      Alam5  = Vlm(5,Iset)
      Alam6  = Vlm(6,Iset)

      Amas4  = Qms(4,Iset)
      Amas5  = Qms(5,Iset)
      Amas6  = Qms(6,Iset)

      Xmin   = Xmn  (Iset)
      Qini   = Qmn  (Iset)
      Qmax   = Qmx  (Iset)

      Do 101 Iexp = 1, Nexp(Iset)
         ExpNor(Iexp) = ExpN(Iexp, Iset)
  101 Continue

      Ctq1Pdf=0D0
      Return
C                         *************************
      END
      FUNCTION Ctq2df (Iset, Iprtn, XX, QQ, Irt)
C                                                   -=-=- ctq2df

C            Returns xf(x,Q) -- the momentum fraction distribution !!
C            Returns valence and sea rather than combined flavor distr.

C            Iset : PDF set label

C            Iprtn  : Parton label:   2, 1 = d_ and u_ valence
C                                     0 = gluon
C                            -1, ... -6 = u, d, s, c, b, t sea quarks

C            XX  : Bjorken-x
C            QQ  : scale parameter "Q"
C      Irt : Return code
C      0 : no error
C      1 : parametrization is slightly negative; reset to 0.0.
C          (This condition happens rarely -- only for large x where the 
C           absolute value of the parton distribution is extremely small.) 


      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (Nex = 5, MxFl = 6, Npn = 3, Nst = 30, Nexpt=20)
      Parameter (Nst4 = Nst*4)

      DIMENSION
     >   Iord(Nst), Isch(Nst), Nqrk(Nst),Alm(Nst)
     > , Vlm(4:6,Nst), Qms(4:6, Nst)
     > , Xmn(Nst), Qmn(Nst), Qmx(Nst), Nexp(Nexpt)
     > , Mex(Nst), Mpn(Nst), ExpN(Nexpt, Nst), ExpNor(Nexpt)

c                                          Run le26 - CTEQ2M
c
      DATA 
     >  Isch(1), Iord(1), Nqrk(1), Alm(1) /  1,  2,  6, .213  / 
     >  (Vlm(I,1), I=4,6) / .213,    .139,     .053   /
     >  (Qms(I,1), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(1), Qmn(1), Qmx(1) /  1.E-5,  1.60,  1.E3  /
     >  Mex(1), Mpn(1), Nexp(1) /  5, 3, 10 /
     >  (ExpN(I,1), I=1,10 )
     >  / 0.990,1.012,1.022,0.980,1.062,0.870,0.843
     >   ,0.815,0.974,1.029 /

c                                          Run sa17 - CTEQ2MS
c
      DATA 
     >  Isch(2), Iord(2), Nqrk(2), Alm(2) /  1,  2,  6, .208  / 
     >  (Vlm(I,2), I=4,6) / .208,    .135,     .051   /
     >  (Qms(I,2), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(2), Qmn(2), Qmx(2) /  1.E-5,  1.60,  1.E3  /
     >  Mex(2), Mpn(2), Nexp(2) /  5, 3, 10 /
     >  (ExpN(I,2), I=1,10 )
     >  / 0.992,1.017,1.023,0.982,1.079,0.879,0.845
     >   ,0.814,0.984,1.036 /

c                                          Run fa06 - CTEQ2MF
c
      DATA 
     >  Isch(3), Iord(3), Nqrk(3), Alm(3) /  1,  2,  6, .208  / 
     >  (Vlm(I,3), I=4,6) / .208,    .135,     .051   /
     >  (Qms(I,3), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(3), Qmn(3), Qmx(3) /  1.E-5,  1.60,  1.E3  /
     >  Mex(3), Mpn(3), Nexp(3) /  5, 3, 10 /
     >  (ExpN(I,3), I=1,10 )
     >  / 0.989,1.014,1.021,0.977,1.099,0.812,0.789
     >   ,0.822,0.909,0.950 /

c                                          Run ll25 - CTEQ2ML
c
      DATA 
     >  Isch(4), Iord(4), Nqrk(4), Alm(4) /  1,  2,  6, .322  / 
     >  (Vlm(I,4), I=4,6) / .322,    .220,     .088   /
     >  (Qms(I,4), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(4), Qmn(4), Qmx(4) /  1.E-5,  1.60,  1.E3  /
     >  Mex(4), Mpn(4), Nexp(4) /  5, 3, 10 /
     >  (ExpN(I,4), I=1,10 )
     >  / 0.990,1.009,1.022,0.978,1.050,1.027,0.955
     >   ,0.848,0.985,1.053 /

c                                          Run lo24 - CTEQ2L
c
      DATA 
     >  Isch(5), Iord(5), Nqrk(5), Alm(5) /  0,  1,  6, .190  / 
     >  (Vlm(I,5), I=4,6) / .190,    .143,     .072   /
     >  (Qms(I,5), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(5), Qmn(5), Qmx(5) /  1.E-5,  1.60,  1.E3  /
     >  Mex(5), Mpn(5), Nexp(5) /  5, 3, 10 /
     >  (ExpN(I,5), I=1,10 )
     >  / 0.990,1.018,1.022,0.982,0.725,0.931,0.855
     >   ,0.795,0.965,1.007 /

c                                          Run da06 - CTEQ2D
c
      DATA 
     >  Isch(6), Iord(6), Nqrk(6), Alm(6) /  2,  2,  6, .235  / 
     >  (Vlm(I,6), I=4,6) / .235,    .155,     .060   /
     >  (Qms(I,6), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(6), Qmn(6), Qmx(6) /  1.E-5,  1.60,  1.E3  /
     >  Mex(6), Mpn(6), Nexp(6) /  5, 3, 10 /
     >  (ExpN(I,6), I=1,10 )
     >  / 0.989,1.019,1.021,0.978,0.961,0.958,0.904
     >   ,0.857,0.965,1.022 /

      Data Ist, Lp, Qsto / 0, -10, 1.2345 /

      save Ist, Lp, Qsto
      save SB, SB2, SB3

      X  = XX
      Irt = 0
      if(Iset.eq.Ist .and. Qsto.eq.QQ) then
C                                             if only change is in x:
        if (Iprtn.eq.Lp) goto 100
C                         if change in flv is within "light" partons:
        if (Iprtn.ge.-3 .and. Lp.ge.-3) goto 501
      endif

      Ip = abs(Iprtn)
C                                                  Set up Qi for SB
      If (Ip .GE. 4) then
         If (QQ .LE. Qms(Ip, Iset)) Then
           Ctq2df = 0.0
           Return
         Endif
         Qi = Qms(ip, Iset)
      Else
         Qi = Qmn(Iset)
      Endif
C                   Use "standard lambda" of parametrization program
      Alam = Alm (Iset)

      SBL = LOG(QQ/Alam) / LOG(Qi/Alam)
      SB = LOG (SBL)
      SB2 = SB*SB
      SB3 = SB2*SB

 501  Iflv = 3 - Iprtn

      Goto (1,2,3,4,5,6, 311) Iset

 1    Goto(11,12,13,14,15,16,17,18,19)Iflv    
c   Ifl =   2
  11  A0=Exp( 0.2143E+00+0.8417E+00*SB -0.2451E+01*SB2+0.9875E+00*SB3)
      A1= 0.5209E+00-0.2384E+00*SB +0.5086E+00*SB2-0.2123E+00*SB3 
      A2= 0.3178E+01+0.5258E+01*SB -0.8102E+01*SB2+0.3334E+01*SB3 
      A3=-0.8537E+00+0.5921E+01*SB -0.1007E+02*SB2+0.4146E+01*SB3 
      A4= 0.1821E+01+0.2822E-01*SB +0.1662E+00*SB2-0.1058E+00*SB3 
      A5= 0.0000E+00-0.1090E+01*SB +0.3136E+01*SB2-0.1301E+01*SB3 
      goto 100
c   Ifl =   1
  12  A0=Exp(-0.1314E+01-0.1342E-01*SB +0.1136E+00*SB2-0.1557E+00*SB3)
      A1= 0.2780E+00+0.2558E-01*SB +0.4467E-02*SB2-0.2472E-02*SB3 
      A2= 0.3672E+01+0.5324E+00*SB +0.3531E-01*SB2+0.7928E-03*SB3 
      A3= 0.2957E+02-0.2000E+02*SB +0.5929E+01*SB2+0.3390E+00*SB3 
      A4= 0.8069E+00-0.2877E+00*SB +0.3574E-01*SB2+0.5622E-02*SB3 
      A5= 0.0000E+00+0.2287E+00*SB -0.4052E-01*SB2+0.5589E-01*SB3 
      goto 100
c   Ifl =   0
  13  A0=Exp(-0.1059E+00-0.1461E+01*SB -0.2544E+00*SB2+0.4526E-01*SB3)
      A1=-0.2578E+00+0.1385E+00*SB -0.1383E+00*SB2+0.3811E-01*SB3 
      A2= 0.5195E+01+0.9648E+00*SB -0.2103E+00*SB2-0.6701E-01*SB3 
      A3= 0.5131E+01+0.2151E+01*SB -0.2880E+01*SB2+0.6608E+00*SB3 
      A4= 0.1118E+01+0.2636E+00*SB -0.5140E+00*SB2+0.1613E+00*SB3 
      A5= 0.0000E+00+0.2456E+01*SB -0.8741E+00*SB2+0.2136E+00*SB3 
      goto 100
c   Ifl =  -1
  14  A0=Exp(-0.2732E+00-0.3523E+01*SB +0.3657E+01*SB2-0.1415E+01*SB3)
      A1=-0.3807E+00+0.1211E+00*SB -0.1231E+00*SB2+0.3753E-01*SB3 
      A2= 0.9698E+01-0.2596E+01*SB +0.2412E+01*SB2-0.9257E+00*SB3 
      A3=-0.6165E+00+0.1120E+01*SB -0.1708E+01*SB2+0.6383E+00*SB3 
      A4= 0.7292E-01-0.1339E+00*SB +0.2104E+00*SB2-0.7987E-01*SB3 
      A5=-0.1370E+01+0.2452E+01*SB -0.1804E+01*SB2+0.6459E+00*SB3 
      goto 100
c   Ifl =  -2
  15  A0=Exp(-0.2319E+01-0.3182E+01*SB +0.3572E+01*SB2-0.1431E+01*SB3)
      A1=-0.2622E+00+0.3085E+00*SB -0.4394E+00*SB2+0.1496E+00*SB3 
      A2= 0.9481E+01-0.3627E+01*SB +0.5640E+01*SB2-0.2265E+01*SB3 
      A3= 0.5000E+02-0.1851E+02*SB +0.2640E+01*SB2-0.6001E+00*SB3 
      A4= 0.1566E+01-0.7375E+00*SB +0.8736E+00*SB2-0.3449E+00*SB3 
      A5=-0.7983E-01+0.3236E+01*SB -0.3373E+01*SB2+0.1236E+01*SB3 
      goto 100
c   Ifl =  -3
  16  A0=Exp(-0.1855E+01-0.5302E+01*SB +0.8433E+00*SB2-0.1236E+00*SB3)
      A1=-0.4000E-02-0.1345E+01*SB +0.1192E+01*SB2-0.3039E+00*SB3 
      A2= 0.6870E+01+0.1246E+01*SB -0.8968E+00*SB2-0.9791E-01*SB3 
      A3= 0.0000E+00+0.4616E+01*SB +0.1026E+02*SB2+0.2844E+02*SB3 
      A4= 0.1000E-02+0.4098E+00*SB -0.4250E+00*SB2+0.1100E+00*SB3 
      A5= 0.0000E+00-0.2151E+01*SB +0.2991E+01*SB2-0.7717E+00*SB3 
      goto 100
c   Ifl =  -4
  17  A0=SB** 0.7722E+00*Exp(-0.7241E+01-0.7885E-01*SB -0.1124E+01*SB2)
      A1=-0.3971E+00+0.9132E+00*SB -0.1175E+01*SB2+0.3573E+00*SB3 
      A2= 0.6367E+01-0.6565E+01*SB +0.8114E+01*SB2-0.2666E+01*SB3 
      A3= 0.2878E+02-0.2000E+02*SB +0.7000E+00*SB2+0.3000E+02*SB3 
      A4= 0.1010E+00-0.4592E+00*SB +0.5877E+00*SB2-0.1472E+00*SB3 
      A5= 0.1749E+00+0.3875E+01*SB -0.3768E+01*SB2+0.1316E+01*SB3 
      goto 100
c   Ifl =  -5
  18  A0=SB** 0.1299E+00*Exp(-0.4868E+01-0.4339E+01*SB +0.7080E+00*SB2)
      A1=-0.1705E+00-0.3381E+00*SB +0.5287E+00*SB2-0.2644E+00*SB3 
      A2= 0.5610E+01-0.1365E+01*SB +0.1835E+01*SB2-0.5655E+00*SB3 
      A3=-0.1001E+01+0.3044E+01*SB +0.2680E+01*SB2+0.1426E+02*SB3 
      A4= 0.3814E-02+0.3430E+00*SB -0.6926E+00*SB2+0.3486E+00*SB3 
      A5= 0.1156E+01+0.2016E+01*SB -0.1674E+01*SB2+0.5981E+00*SB3 
      goto 100
c   Ifl =  -6
  19  A0=SB** 0.9819E+00*Exp(-0.7859E+01+0.6819E+00*SB -0.3386E+01*SB2)
      A1=-0.1055E+00-0.1413E+01*SB +0.3451E+01*SB2-0.2466E+01*SB3 
      A2= 0.4055E+01+0.8107E+01*SB -0.1576E+02*SB2+0.8094E+01*SB3 
      A3= 0.3799E+01+0.9616E+01*SB -0.1984E+02*SB2+0.2641E+02*SB3 
      A4= 0.3619E+00-0.8627E+00*SB -0.9390E-01*SB2+0.9196E+00*SB3 
      A5= 0.3779E+01-0.6073E+01*SB +0.9999E+01*SB2-0.4304E+01*SB3 
      goto 100

 2    Goto(21,22,23,24,25,26,27,28,29)Iflv    
c   Ifl =   2
  21  A0=Exp( 0.2790E+00+0.7294E+00*SB -0.2202E+01*SB2+0.8599E+00*SB3)
      A1= 0.5380E+00-0.2261E+00*SB +0.4636E+00*SB2-0.1871E+00*SB3 
      A2= 0.3259E+01+0.2141E+01*SB -0.2947E+01*SB2+0.1245E+01*SB3 
      A3=-0.8390E+00+0.1448E+01*SB -0.2331E+01*SB2+0.8658E+00*SB3 
      A4= 0.1847E+01-0.3943E+01*SB +0.5998E+01*SB2-0.2191E+01*SB3 
      A5= 0.0000E+00-0.9719E+00*SB +0.2830E+01*SB2-0.1137E+01*SB3 
      goto 100
c   Ifl =   1
  22  A0=Exp(-0.1318E+01+0.2328E-01*SB +0.5179E-01*SB2-0.1305E+00*SB3)
      A1= 0.2760E+00+0.4429E-01*SB -0.2626E-01*SB2+0.7143E-02*SB3 
      A2= 0.3660E+01+0.5232E+00*SB +0.5491E-01*SB2-0.4115E-02*SB3 
      A3= 0.2910E+02-0.2000E+02*SB +0.6631E+01*SB2-0.3050E-01*SB3 
      A4= 0.8010E+00-0.2688E+00*SB +0.1051E-01*SB2+0.1195E-01*SB3 
      A5= 0.0000E+00+0.2887E+00*SB -0.1398E+00*SB2+0.8194E-01*SB3 
      goto 100
c   Ifl =   0
  23  A0=Exp(-0.1623E+01-0.7232E+00*SB +0.1889E+00*SB2+0.1140E+00*SB3)
      A1=-0.5000E+00+0.8611E-01*SB +0.2203E-01*SB2-0.1401E-01*SB3 
      A2= 0.3821E+01+0.8976E+00*SB +0.1400E+00*SB2-0.9163E-01*SB3 
      A3= 0.5809E+01-0.5060E+01*SB +0.3808E+00*SB2+0.2519E+00*SB3 
      A4= 0.4500E+00-0.5121E+00*SB +0.1979E+00*SB2-0.2705E-01*SB3 
      A5= 0.0000E+00+0.1210E+01*SB -0.2921E+00*SB2+0.1240E+00*SB3 
      goto 100
c   Ifl =  -1
  24  A0=Exp(-0.6986E-01-0.5954E+00*SB -0.1582E+01*SB2+0.5104E+00*SB3)
      A1=-0.8461E+00+0.2127E+00*SB +0.9425E-01*SB2-0.5264E-01*SB3 
      A2= 0.1200E+02+0.1659E+01*SB -0.5354E+01*SB2+0.1795E+01*SB3 
      A3= 0.2958E+02+0.3000E+02*SB +0.3000E+02*SB2-0.1965E+02*SB3 
      A4= 0.4000E+01-0.4865E+00*SB +0.9460E+00*SB2+0.3432E+00*SB3 
      A5=-0.3378E+01+0.1656E+01*SB +0.1123E+01*SB2-0.4667E+00*SB3 
      goto 100
c   Ifl =  -2
  25  A0=Exp(-0.1929E+01-0.2626E+01*SB +0.2926E+01*SB2-0.1297E+01*SB3)
      A1=-0.6627E+00+0.4561E+00*SB -0.3818E+00*SB2+0.1239E+00*SB3 
      A2= 0.9506E+01-0.2724E+01*SB +0.4283E+01*SB2-0.1804E+01*SB3 
      A3= 0.1897E+02+0.1642E+01*SB -0.8390E+01*SB2+0.3894E+01*SB3 
      A4= 0.1024E+01-0.1786E+00*SB +0.4535E+00*SB2-0.2075E+00*SB3 
      A5=-0.1746E+01+0.3572E+01*SB -0.2908E+01*SB2+0.1093E+01*SB3 
      goto 100
c   Ifl =  -3
  26  A0=Exp(-0.4913E+00-0.6866E+01*SB +0.1432E+01*SB2-0.1749E+00*SB3)
      A1=-0.1157E+00-0.1567E+01*SB +0.1439E+01*SB2-0.3724E+00*SB3 
      A2= 0.7730E+01+0.9748E+00*SB -0.1157E+01*SB2-0.8358E-02*SB3 
      A3=-0.6050E+00+0.1835E+01*SB +0.3788E+01*SB2+0.3000E+02*SB3 
      A4= 0.1620E-08+0.4590E+00*SB -0.4070E+00*SB2+0.8900E-01*SB3 
      A5=-0.7048E+00-0.2505E+01*SB +0.4000E+01*SB2-0.1161E+01*SB3 
      goto 100
c   Ifl =  -4
  27  A0=SB** 0.7393E+00*Exp(-0.6518E+01-0.3998E+00*SB -0.1111E+01*SB2)
      A1=-0.6482E+00+0.1125E+01*SB -0.1290E+01*SB2+0.3940E+00*SB3 
      A2= 0.8487E+01-0.9235E+01*SB +0.9353E+01*SB2-0.2913E+01*SB3 
      A3= 0.2265E+02-0.1999E+02*SB +0.4105E+01*SB2+0.2144E+02*SB3 
      A4= 0.8990E-01-0.4372E+00*SB +0.5941E+00*SB2-0.1469E+00*SB3 
      A5=-0.9690E+00+0.5068E+01*SB -0.4368E+01*SB2+0.1503E+01*SB3 
      goto 100
c   Ifl =  -5
  28  A0=SB** 0.9880E+00*Exp(-0.7180E+01-0.2494E+01*SB +0.3561E-01*SB2)
      A1=-0.4301E+00-0.2611E+00*SB +0.3914E+00*SB2-0.1638E+00*SB3 
      A2= 0.5137E+01+0.1506E+01*SB -0.9588E+00*SB2-0.1596E+00*SB3 
      A3= 0.1483E+02+0.2998E+02*SB +0.2357E+02*SB2-0.9353E+01*SB3 
      A4= 0.2426E+00+0.1371E+00*SB -0.3791E+00*SB2+0.1948E+00*SB3 
      A5= 0.1463E+01+0.1907E+00*SB +0.3557E+00*SB2+0.2097E-01*SB3 
      goto 100
c   Ifl =  -6
  29  A0=SB** 0.1005E+01*Exp(-0.5255E+01-0.9866E-01*SB -0.2737E+01*SB2)
      A1=-0.3140E+00-0.2055E+00*SB +0.5594E+00*SB2-0.2960E+00*SB3 
      A2= 0.9227E+01-0.4569E+01*SB -0.9724E+01*SB2+0.1026E+02*SB3 
      A3= 0.1131E+02-0.1972E+02*SB -0.1107E+02*SB2+0.2311E+02*SB3 
      A4= 0.1488E+01+0.1737E+01*SB +0.4323E+01*SB2-0.9925E+01*SB3 
      A5= 0.1895E+01-0.7350E+00*SB +0.3780E+01*SB2-0.1408E+01*SB3 
      goto 100

 3    Goto(31,32,33,34,35,36,37,38,39)Iflv    
c   Ifl =   2
  31  A0=Exp(-0.7913E+00-0.2789E+01*SB -0.7289E-01*SB2+0.1770E+00*SB3)
      A1= 0.4942E+00-0.7886E-01*SB +0.9057E-01*SB2-0.5259E-01*SB3 
      A2= 0.3727E+01+0.1089E+01*SB -0.1004E+01*SB2+0.4345E+00*SB3 
      A3= 0.1944E+01+0.7846E+01*SB +0.7984E+01*SB2+0.5548E+01*SB3 
      A4= 0.2940E-02+0.8428E-04*SB +0.1266E+00*SB2-0.3517E-01*SB3 
      A5=-0.1060E+00-0.1192E-01*SB +0.1130E+01*SB2-0.4527E+00*SB3 
      goto 100
c   Ifl =   1
  32  A0=Exp(-0.1344E+01+0.7859E-02*SB +0.4623E-01*SB2-0.1273E+00*SB3)
      A1= 0.2760E+00+0.4201E-01*SB -0.1795E-01*SB2+0.3212E-02*SB3 
      A2= 0.3660E+01+0.5247E+00*SB +0.4405E-01*SB2+0.1391E-02*SB3 
      A3= 0.2981E+02-0.2000E+02*SB +0.6566E+01*SB2+0.2479E-01*SB3 
      A4= 0.7950E+00-0.2732E+00*SB +0.2470E-01*SB2+0.6157E-02*SB3 
      A5= 0.0000E+00+0.2793E+00*SB -0.9197E-01*SB2+0.5953E-01*SB3 
      goto 100
c   Ifl =   0
  33  A0=Exp( 0.9746E+00-0.3252E+01*SB +0.1664E+01*SB2-0.6410E+00*SB3)
      A1=-0.5271E-02-0.3198E+00*SB +0.1279E+00*SB2-0.1256E-02*SB3 
      A2= 0.5740E+01-0.3139E+01*SB +0.3841E+01*SB2-0.1415E+01*SB3 
      A3= 0.7161E-01-0.4363E+01*SB +0.4925E+01*SB2-0.1614E+01*SB3 
      A4= 0.1860E+01+0.1342E+01*SB -0.2234E+01*SB2+0.1047E+01*SB3 
      A5= 0.7409E-01+0.2390E+01*SB -0.1457E+01*SB2+0.5853E+00*SB3 
      goto 100
c   Ifl =  -1
  34  A0=Exp(-0.8454E+00-0.3334E+01*SB +0.3591E+01*SB2-0.1485E+01*SB3)
      A1=-0.2826E-02-0.2810E+00*SB -0.3809E-01*SB2+0.6585E-01*SB3 
      A2= 0.9139E+01-0.2811E+01*SB +0.4730E+01*SB2-0.2157E+01*SB3 
      A3=-0.3120E+00+0.1217E+01*SB -0.1726E+01*SB2+0.6220E+00*SB3 
      A4= 0.1793E-01-0.4608E-01*SB +0.5294E-01*SB2-0.1709E-01*SB3 
      A5=-0.1471E+00+0.1104E+01*SB -0.1358E+01*SB2+0.7200E+00*SB3 
      goto 100
c   Ifl =  -2
  35  A0=Exp(-0.1398E+01-0.3536E+01*SB +0.3849E+01*SB2-0.1549E+01*SB3)
      A1=-0.1332E-01-0.2155E-01*SB -0.3404E+00*SB2+0.1569E+00*SB3 
      A2= 0.9981E+01-0.3499E+01*SB +0.5448E+01*SB2-0.2198E+01*SB3 
      A3= 0.3736E+02-0.2000E+02*SB +0.6675E+01*SB2-0.7276E+00*SB3 
      A4= 0.1705E+01-0.1013E+01*SB +0.1122E+01*SB2-0.4057E+00*SB3 
      A5=-0.1189E-01+0.2698E+01*SB -0.3429E+01*SB2+0.1389E+01*SB3 
      goto 100
c   Ifl =  -3
  36  A0=Exp(-0.2979E+01-0.6085E+01*SB +0.2428E+01*SB2-0.6482E+00*SB3)
      A1=-0.1372E+00-0.1281E+00*SB +0.1587E+00*SB2-0.9637E-01*SB3 
      A2= 0.7009E+01-0.1609E+01*SB +0.2765E+01*SB2-0.1177E+01*SB3 
      A3= 0.1308E+01+0.9583E+01*SB +0.2360E+02*SB2+0.2999E+02*SB3 
      A4= 0.2509E-01+0.2106E+00*SB -0.4405E+00*SB2+0.2075E+00*SB3 
      A5=-0.2069E-01+0.1971E+01*SB -0.1615E+01*SB2+0.6039E+00*SB3 
      goto 100
c   Ifl =  -4
  37  A0=SB** 0.8072E+00*Exp(-0.6920E+01-0.5031E+00*SB -0.9965E+00*SB2)
      A1=-0.2118E+00+0.7930E+00*SB -0.1101E+01*SB2+0.3302E+00*SB3 
      A2= 0.8039E+01-0.7170E+01*SB +0.8657E+01*SB2-0.2893E+01*SB3 
      A3= 0.2926E+02-0.1993E+02*SB +0.1841E+01*SB2+0.2996E+02*SB3 
      A4= 0.1339E+00-0.5531E+00*SB +0.6505E+00*SB2-0.1595E+00*SB3 
      A5= 0.7439E+00+0.3307E+01*SB -0.3284E+01*SB2+0.1152E+01*SB3 
      goto 100
c   Ifl =  -5
  38  A0=SB** 0.9925E+00*Exp(-0.2190E+01-0.3393E+01*SB -0.8631E+00*SB2)
      A1=-0.1261E+00-0.2368E+00*SB +0.4143E+00*SB2-0.1577E+00*SB3 
      A2= 0.4585E+01+0.5227E+01*SB -0.3248E+01*SB2-0.2599E+00*SB3 
      A3=-0.1094E+01+0.4927E+00*SB -0.9921E+00*SB2+0.3138E+01*SB3 
      A4= 0.1396E+00+0.2562E+00*SB +0.1844E+00*SB2-0.1599E+00*SB3 
      A5= 0.8621E+00+0.4715E+00*SB +0.2547E+01*SB2-0.8429E+00*SB3 
      goto 100
c   Ifl =  -6
  39  A0=SB** 0.1016E+01*Exp(-0.5397E+01-0.1979E+01*SB -0.2441E+00*SB2)
      A1=-0.1426E+00-0.2861E+00*SB +0.7434E+00*SB2-0.5214E+00*SB3 
      A2= 0.6363E+01+0.4028E+00*SB -0.8356E+01*SB2+0.6814E+01*SB3 
      A3=-0.2526E+00+0.2425E+01*SB -0.1407E+02*SB2+0.3000E+02*SB3 
      A4= 0.1125E+00-0.1089E+01*SB +0.9977E+01*SB2+0.1000E+02*SB3 
      A5= 0.2669E+01-0.6366E+00*SB +0.4355E+01*SB2-0.2919E+01*SB3 
      goto 100

 4    Goto(41,42,43,44,45,46,47,48,49)Iflv    
c   Ifl =   2
  41  A0=Exp( 0.3760E+00+0.5491E+00*SB -0.1845E+01*SB2+0.6803E+00*SB3)
      A1= 0.5650E+00-0.1953E+00*SB +0.3761E+00*SB2-0.1419E+00*SB3 
      A2= 0.3464E+01+0.3817E+01*SB -0.5384E+01*SB2+0.2057E+01*SB3 
      A3=-0.5850E+00+0.5566E+01*SB -0.9000E+01*SB2+0.3433E+01*SB3 
      A4= 0.2322E+01-0.1431E+00*SB +0.3901E+00*SB2-0.1678E+00*SB3 
      A5= 0.0000E+00-0.7370E+00*SB +0.2310E+01*SB2-0.8743E+00*SB3 
      goto 100
c   Ifl =   1
  42  A0=Exp(-0.1324E+01+0.1169E-01*SB +0.1969E-01*SB2-0.7583E-01*SB3)
      A1= 0.2890E+00+0.5832E-01*SB -0.2921E-01*SB2+0.4701E-02*SB3 
      A2= 0.3580E+01+0.5291E+00*SB -0.5662E-02*SB2+0.2746E-01*SB3 
      A3= 0.3021E+02-0.1999E+02*SB +0.6250E+01*SB2-0.3035E+00*SB3 
      A4= 0.7990E+00-0.2531E+00*SB +0.5556E-02*SB2+0.8272E-02*SB3 
      A5= 0.0000E+00+0.3674E+00*SB -0.1383E+00*SB2+0.4665E-01*SB3 
      goto 100
c   Ifl =   0
  43  A0=Exp(-0.1920E+00-0.7015E+00*SB -0.9113E+00*SB2+0.2352E+00*SB3)
      A1=-0.2120E+00+0.1133E-01*SB -0.1553E-01*SB2+0.2822E-02*SB3 
      A2= 0.4549E+01+0.1250E+01*SB -0.4647E+00*SB2+0.9617E-01*SB3 
      A3= 0.1197E+02-0.4156E+01*SB +0.1413E+00*SB2+0.1607E+00*SB3 
      A4= 0.1616E+01+0.1082E+00*SB -0.6651E+00*SB2+0.2356E+00*SB3 
      A5= 0.0000E+00+0.1824E+01*SB -0.2063E+00*SB2+0.1148E-01*SB3 
      goto 100
c   Ifl =  -1
  44  A0=Exp(-0.1388E+01-0.7408E+00*SB -0.6454E+00*SB2+0.2373E+00*SB3)
      A1=-0.2928E+00-0.1726E-01*SB +0.4033E-01*SB2-0.2514E-01*SB3 
      A2= 0.9975E+01-0.2048E+01*SB -0.6060E+00*SB2+0.5225E+00*SB3 
      A3= 0.2687E+02-0.4683E+01*SB -0.1999E+02*SB2+0.1188E+02*SB3 
      A4= 0.4000E+01-0.6773E+00*SB +0.4301E+00*SB2+0.4524E+00*SB3 
      A5=-0.7164E+00+0.7488E+00*SB +0.5766E+00*SB2-0.2609E+00*SB3 
      goto 100
c   Ifl =  -2
  45  A0=Exp(-0.2272E+01-0.2998E+01*SB +0.3282E+01*SB2-0.1203E+01*SB3)
      A1=-0.2062E+00+0.3320E+00*SB -0.5074E+00*SB2+0.1655E+00*SB3 
      A2= 0.9667E+01-0.3497E+01*SB +0.5271E+01*SB2-0.1984E+01*SB3 
      A3= 0.4996E+02-0.3241E+01*SB -0.1425E+02*SB2+0.3849E+01*SB3 
      A4= 0.1619E+01-0.5354E+00*SB +0.5753E+00*SB2-0.2238E+00*SB3 
      A5= 0.8755E-01+0.3195E+01*SB -0.3496E+01*SB2+0.1197E+01*SB3 
      goto 100
c   Ifl =  -3
  46  A0=Exp(-0.1864E+01-0.5258E+01*SB +0.1034E+01*SB2-0.1550E+00*SB3)
      A1= 0.1000E-02-0.1090E+01*SB +0.8345E+00*SB2-0.1887E+00*SB3 
      A2= 0.6898E+01-0.4951E+00*SB +0.4279E+00*SB2-0.2727E+00*SB3 
      A3= 0.0000E+00+0.4322E+01*SB +0.8181E+01*SB2+0.2309E+02*SB3 
      A4= 0.1000E-02+0.3550E+00*SB -0.3220E+00*SB2+0.7294E-01*SB3 
      A5= 0.0000E+00-0.1347E+01*SB +0.1896E+01*SB2-0.4491E+00*SB3 
      goto 100
c   Ifl =  -4
  47  A0=SB** 0.7528E+00*Exp(-0.7684E+01+0.6791E-01*SB -0.9094E+00*SB2)
      A1=-0.3732E+00+0.8408E+00*SB -0.1020E+01*SB2+0.3046E+00*SB3 
      A2= 0.4984E+01-0.5534E+01*SB +0.6418E+01*SB2-0.1856E+01*SB3 
      A3= 0.3761E+02-0.1999E+02*SB -0.3358E+01*SB2+0.2999E+02*SB3 
      A4= 0.1161E+00-0.4680E+00*SB +0.5567E+00*SB2-0.1633E+00*SB3 
      A5= 0.3028E+00+0.3339E+01*SB -0.3004E+01*SB2+0.9160E+00*SB3 
      goto 100
c   Ifl =  -5
  48  A0=SB** 0.1011E+01*Exp(-0.7217E+01-0.2288E+01*SB +0.3450E+00*SB2)
      A1=-0.1955E+00-0.3371E+00*SB +0.5111E+00*SB2-0.2210E+00*SB3 
      A2= 0.4302E+01-0.1214E+01*SB +0.3104E+01*SB2-0.1408E+01*SB3 
      A3= 0.1487E+02+0.1549E+02*SB +0.2875E+02*SB2-0.1922E+02*SB3 
      A4= 0.8935E-02+0.3571E+00*SB -0.6668E+00*SB2+0.3037E+00*SB3 
      A5= 0.1570E+01+0.7105E+00*SB -0.6070E+00*SB2+0.3796E+00*SB3 
      goto 100
c   Ifl =  -6
  49  A0=SB** 0.9986E+00*Exp(-0.5847E+01-0.2798E+00*SB -0.9882E+00*SB2)
      A1=-0.2154E+00-0.8282E-01*SB +0.3611E-01*SB2+0.2623E-01*SB3 
      A2= 0.3250E+01+0.9635E+01*SB -0.1274E+02*SB2+0.4453E+01*SB3 
      A3=-0.2594E+01+0.9097E+01*SB +0.1581E+02*SB2-0.9123E+01*SB3 
      A4= 0.1768E+01-0.2749E+01*SB +0.9999E+01*SB2+0.9995E+01*SB3 
      A5= 0.2521E+01-0.1802E-01*SB +0.4820E+00*SB2+0.2004E+00*SB3 
      goto 100

 5    Goto(51,52,53,54,55,56,57,58,59)Iflv    
c   Ifl =   2
  51  A0=Exp( 0.7248E-01+0.3941E+00*SB -0.1772E+01*SB2+0.7629E+00*SB3)
      A1= 0.4964E+00-0.1224E+00*SB +0.3646E+00*SB2-0.1685E+00*SB3 
      A2= 0.3000E+01+0.2780E+01*SB -0.4028E+01*SB2+0.1816E+01*SB3 
      A3=-0.1064E+01+0.3062E+01*SB -0.5927E+01*SB2+0.2785E+01*SB3 
      A4= 0.3193E+01+0.1499E+01*SB -0.2765E+01*SB2+0.1019E+01*SB3 
      A5= 0.1524E-01-0.4541E+00*SB +0.2281E+01*SB2-0.1033E+01*SB3 
      goto 100
c   Ifl =   1
  52  A0=Exp(-0.1794E+01-0.2055E+00*SB -0.3350E-01*SB2-0.5084E-01*SB3)
      A1= 0.1748E+00+0.4637E-01*SB -0.2048E-01*SB2+0.2596E-02*SB3 
      A2= 0.3321E+01+0.6253E+00*SB +0.2148E-01*SB2+0.1288E-01*SB3 
      A3= 0.4355E+02-0.2000E+02*SB +0.5486E+01*SB2+0.1536E+00*SB3 
      A4= 0.9586E+00-0.3217E+00*SB +0.4458E-01*SB2-0.1404E-03*SB3 
      A5=-0.6595E-02+0.3499E+00*SB -0.7048E-01*SB2+0.2619E-01*SB3 
      goto 100
c   Ifl =   0
  53  A0=Exp(-0.6194E+00-0.2643E+00*SB -0.1875E+01*SB2+0.6011E+00*SB3)
      A1=-0.2600E+00+0.8704E-01*SB -0.7375E-01*SB2+0.1876E-01*SB3 
      A2= 0.4620E+01+0.1578E+01*SB -0.8411E+00*SB2+0.1527E+00*SB3 
      A3= 0.1604E+02-0.1230E+02*SB +0.6939E+01*SB2-0.2012E+01*SB3 
      A4= 0.1255E+01+0.4769E+00*SB -0.9915E+00*SB2+0.3439E+00*SB3 
      A5= 0.1116E-02+0.2409E+01*SB -0.4442E+00*SB2+0.3431E-01*SB3 
      goto 100
c   Ifl =  -1
  54  A0=Exp(-0.1571E+01-0.1905E+00*SB -0.8672E+00*SB2+0.2070E+00*SB3)
      A1=-0.3266E+00+0.6428E-01*SB -0.8694E-01*SB2+0.1778E-01*SB3 
      A2= 0.8921E+01-0.5010E+00*SB -0.9658E+00*SB2+0.3893E+00*SB3 
      A3= 0.1329E+02+0.4652E+01*SB -0.2000E+02*SB2+0.1001E+02*SB3 
      A4= 0.3283E+01-0.3400E+00*SB -0.1957E+00*SB2+0.8063E+00*SB3 
      A5=-0.5701E+00+0.4042E+00*SB +0.5239E+00*SB2-0.1665E+00*SB3 
      goto 100
c   Ifl =  -2
  55  A0=Exp(-0.2281E+01-0.2768E+01*SB +0.3137E+01*SB2-0.1278E+01*SB3)
      A1=-0.2624E+00+0.4142E+00*SB -0.5936E+00*SB2+0.1937E+00*SB3 
      A2= 0.9438E+01-0.3179E+01*SB +0.5107E+01*SB2-0.2179E+01*SB3 
      A3= 0.5000E+02-0.1802E+02*SB -0.7515E+01*SB2+0.2991E+01*SB3 
      A4= 0.1809E+01-0.9121E+00*SB +0.8854E+00*SB2-0.3582E+00*SB3 
      A5= 0.4056E-01+0.3033E+01*SB -0.3431E+01*SB2+0.1253E+01*SB3 
      goto 100
c   Ifl =  -3
  56  A0=Exp(-0.2318E+01-0.4104E+01*SB -0.1502E+00*SB2+0.1693E+00*SB3)
      A1=-0.2251E-01-0.1101E+01*SB +0.1037E+01*SB2-0.3290E+00*SB3 
      A2= 0.6989E+01+0.1794E+01*SB -0.1811E+01*SB2+0.3061E+00*SB3 
      A3= 0.7972E+00+0.7806E+01*SB +0.1869E+02*SB2+0.2999E+02*SB3 
      A4= 0.4795E-01+0.1622E+00*SB -0.3977E+00*SB2+0.1920E+00*SB3 
      A5=-0.5275E-01-0.2616E+01*SB +0.3076E+01*SB2-0.7425E+00*SB3 
      goto 100
c   Ifl =  -4
  57  A0=SB** 0.8431E+00*Exp(-0.6539E+01-0.1875E+00*SB -0.1346E+01*SB2)
      A1=-0.4970E+00+0.9062E+00*SB -0.1169E+01*SB2+0.3703E+00*SB3 
      A2= 0.4939E+01-0.2995E+01*SB +0.4483E+01*SB2-0.1704E+01*SB3 
      A3= 0.3113E+02-0.1997E+02*SB +0.1540E+01*SB2+0.3000E+02*SB3 
      A4= 0.1349E+00-0.5418E+00*SB +0.6142E+00*SB2-0.1360E+00*SB3 
      A5=-0.8590E+00+0.3956E+01*SB -0.3612E+01*SB2+0.1401E+01*SB3 
      goto 100
c   Ifl =  -5
  58  A0=SB** 0.2639E-01*Exp(-0.2099E+01-0.2681E+01*SB +0.2925E+00*SB2)
      A1=-0.2243E+00-0.5343E-01*SB -0.1953E-01*SB2+0.1586E-01*SB3 
      A2= 0.4294E+01+0.1102E+01*SB -0.1822E+00*SB2-0.2481E+00*SB3 
      A3=-0.9998E+00+0.8275E-01*SB +0.5494E+00*SB2-0.1982E+00*SB3 
      A4= 0.5904E-04+0.9222E-01*SB -0.9293E-01*SB2+0.9159E-01*SB3 
      A5= 0.2657E+00+0.1770E+01*SB -0.7111E+00*SB2+0.2525E+00*SB3 
      goto 100
c   Ifl =  -6
  59  A0=SB** 0.1009E+01*Exp(-0.7032E+01+0.4562E+01*SB -0.9081E+01*SB2)
      A1=-0.1412E+00-0.5076E+00*SB +0.9513E+00*SB2-0.4326E+00*SB3 
      A2= 0.5385E+01+0.3023E+01*SB -0.1162E+02*SB2+0.7006E+01*SB3 
      A3= 0.4997E+01-0.1600E+02*SB +0.1342E+02*SB2+0.1197E+02*SB3 
      A4= 0.5825E+00+0.3994E+00*SB -0.1255E+01*SB2+0.6486E+00*SB3 
      A5= 0.3365E+01-0.4026E+01*SB +0.8385E+01*SB2-0.2260E+01*SB3 
      goto 100

 6    Goto(61,62,63,64,65,66,67,68,69)Iflv    
c   Ifl =   2
  61  A0=Exp( 0.1590E+00+0.5580E+00*SB -0.1838E+01*SB2+0.7018E+00*SB3)
      A1= 0.5110E+00-0.1625E+00*SB +0.3547E+00*SB2-0.1412E+00*SB3 
      A2= 0.3158E+01+0.3962E+01*SB -0.5866E+01*SB2+0.2375E+01*SB3 
      A3=-0.6000E+00+0.6144E+01*SB -0.1056E+02*SB2+0.4345E+01*SB3 
      A4= 0.2306E+01-0.4669E-01*SB +0.2711E+00*SB2-0.1640E+00*SB3 
      A5= 0.0000E+00-0.6638E+00*SB +0.2239E+01*SB2-0.8843E+00*SB3 
      goto 100
c   Ifl =   1
  62  A0=Exp(-0.1182E+01+0.1449E+00*SB +0.2753E-01*SB2-0.1009E+00*SB3)
      A1= 0.2540E+00+0.2686E-01*SB -0.1546E-01*SB2+0.5396E-02*SB3 
      A2= 0.3442E+01+0.5576E+00*SB +0.1937E-01*SB2+0.6696E-02*SB3 
      A3= 0.2545E+02-0.2000E+02*SB +0.7355E+01*SB2-0.7058E+00*SB3 
      A4= 0.9170E+00-0.3090E+00*SB +0.1705E-01*SB2+0.8534E-02*SB3 
      A5= 0.0000E+00+0.1449E+00*SB -0.7821E-01*SB2+0.6405E-01*SB3 
      goto 100
c   Ifl =   0
  63  A0=Exp(-0.3410E+00-0.9613E+00*SB -0.4969E+00*SB2+0.9360E-01*SB3)
      A1=-0.2400E+00+0.1473E+00*SB -0.1593E+00*SB2+0.4538E-01*SB3 
      A2= 0.4841E+01+0.9311E+00*SB +0.1601E-03*SB2-0.1331E+00*SB3 
      A3= 0.7427E+01-0.1397E+01*SB +0.1489E+00*SB2-0.2848E+00*SB3 
      A4= 0.9600E+00+0.3697E+00*SB -0.4246E+00*SB2+0.1032E+00*SB3 
      A5= 0.0000E+00+0.2484E+01*SB -0.9908E+00*SB2+0.2568E+00*SB3 
      goto 100
c   Ifl =  -1
  64  A0=Exp( 0.1176E+00-0.3418E+01*SB +0.3529E+01*SB2-0.1367E+01*SB3)
      A1=-0.3654E+00+0.1914E+00*SB -0.2192E+00*SB2+0.6933E-01*SB3 
      A2= 0.1099E+02-0.4281E+01*SB +0.3729E+01*SB2-0.1254E+01*SB3 
      A3=-0.7514E+00+0.7696E+00*SB -0.1134E+01*SB2+0.4245E+00*SB3 
      A4= 0.7690E-01-0.6558E-01*SB +0.8726E-01*SB2-0.3345E-01*SB3 
      A5=-0.1447E+01+0.2617E+01*SB -0.2094E+01*SB2+0.7536E+00*SB3 
      goto 100
c   Ifl =  -2
  65  A0=Exp(-0.2412E+01-0.2522E+01*SB +0.3126E+01*SB2-0.1305E+01*SB3)
      A1=-0.2353E+00+0.3118E+00*SB -0.4864E+00*SB2+0.1689E+00*SB3 
      A2= 0.9017E+01-0.2437E+01*SB +0.4659E+01*SB2-0.2044E+01*SB3 
      A3= 0.5000E+02-0.1158E+02*SB -0.9260E+01*SB2+0.2847E+01*SB3 
      A4= 0.1726E+01-0.6849E+00*SB +0.7864E+00*SB2-0.3300E+00*SB3 
      A5= 0.5080E-01+0.2858E+01*SB -0.3297E+01*SB2+0.1246E+01*SB3 
      goto 100
c   Ifl =  -3
  66  A0=Exp(-0.1966E+01-0.4405E+01*SB +0.2436E+00*SB2+0.4576E-01*SB3)
      A1=-0.4000E-02-0.1229E+01*SB +0.1118E+01*SB2-0.2988E+00*SB3 
      A2= 0.6902E+01+0.1266E+01*SB -0.1068E+01*SB2+0.3062E-01*SB3 
      A3= 0.0000E+00+0.3987E+01*SB +0.9389E+01*SB2+0.1881E+02*SB3 
      A4= 0.1000E-02+0.3528E+00*SB -0.4201E+00*SB2+0.1248E+00*SB3 
      A5= 0.0000E+00-0.2149E+01*SB +0.2925E+01*SB2-0.7609E+00*SB3 
      goto 100
c   Ifl =  -4
  67  A0=SB** 0.7561E+00*Exp(-0.6960E+01+0.5634E-01*SB -0.1170E+01*SB2)
      A1=-0.4232E+00+0.9269E+00*SB -0.1161E+01*SB2+0.3470E+00*SB3 
      A2= 0.6057E+01-0.5790E+01*SB +0.7352E+01*SB2-0.2435E+01*SB3 
      A3= 0.2941E+02-0.1999E+02*SB -0.8345E+00*SB2+0.3000E+02*SB3 
      A4= 0.1069E+00-0.4620E+00*SB +0.5614E+00*SB2-0.1336E+00*SB3 
      A5=-0.1865E+00+0.3953E+01*SB -0.3791E+01*SB2+0.1315E+01*SB3 
      goto 100
c   Ifl =  -5
  68  A0=SB** 0.5661E-02*Exp(-0.2123E+01-0.3026E+01*SB +0.1912E+00*SB2)
      A1=-0.2011E+00-0.1338E-01*SB -0.3974E-01*SB2+0.1948E-01*SB3 
      A2= 0.4906E+01+0.1740E+01*SB -0.1387E+01*SB2+0.1263E+00*SB3 
      A3=-0.1000E+01+0.5767E-01*SB +0.6377E+00*SB2+0.4736E-01*SB3 
      A4= 0.5927E-04+0.1039E+00*SB -0.9797E-01*SB2+0.6881E-01*SB3 
      A5= 0.4017E+00+0.1981E+01*SB -0.7758E+00*SB2+0.2916E+00*SB3 
      goto 100
c   Ifl =  -6
  69  A0=SB** 0.1008E+01*Exp(-0.7211E+01+0.3273E+01*SB -0.6979E+01*SB2)
      A1=-0.1026E+00-0.4948E+00*SB +0.1188E+01*SB2-0.8016E+00*SB3 
      A2= 0.5397E+01+0.2135E+01*SB -0.9531E+01*SB2+0.6115E+01*SB3 
      A3= 0.4966E+01-0.1111E+02*SB +0.4732E+01*SB2+0.1568E+02*SB3 
      A4= 0.5345E+00-0.1935E+00*SB +0.5816E+00*SB2-0.6794E+00*SB3 
      A5= 0.3569E+01-0.3477E+01*SB +0.8756E+01*SB2-0.4139E+01*SB3 
      goto 100

 311  stop 'This option is not currently supported.'

 100  Ctq2df = A0 *(x**A1) *((1.-x)**A2) *(1.+A3*(x**A4))
     $            *(log(1.+1./x))**A5

      if(Ctq2df.lt.0.0) then
        Ctq2df = 0.0
        Irt=1
      endif

      Ist = Iset

      Lp  = Iprtn
      Qsto = QQ

      Return
C                                  -----------------------
      ENTRY Wlamd2 (Iset, Iorder, Neff)

      Iorder = IORD (Iset)
      Wlamd2 = VLM  (Neff, Iset)

      RETURN
C                                  -----------------------
      Entry PrCtq2
     >        (Iset, Iordr, Ischeme, MxFlv,
     >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
     >         Xmin, Qini, Qmax, ExpNor)

C                           Return QCD parameters and Fitting parameters
C                           associated with parton distribution set Iset.
C    Iord    : Order Of Fit
C    Ischeme : (0, 1, 2)  for  (LO, MS-bar-NLO, DIS-NLO) resp.
C    MxFlv   : Maximum number of flavors included
C    Alam_i  : i = 4,5,6  Effective lambda for i-flavors

C    Amas_i  : i = 4,5,6  Mass parameter for flavor i
C    Xmin, Qini, Qmax : self explanary
C    ExpNor(I) : Normalization factor for the experimental data set used in
C                obtaining the best global fit for parton distributions Iset:
C     I = 1,     2,      3,     4,     5,     6,     7,     8,    9,    10
C      BCDMS   NMC90  NMC280  CCFR   E605    WA70   E706   UA6    H1   ZEUS

      Iordr  = Iord (Iset)
      Ischeme= Isch (Iset)
      MxFlv  = Nqrk (Iset)

      Alam4  = Vlm(4,Iset)
      Alam5  = Vlm(5,Iset)
      Alam6  = Vlm(6,Iset)

      Amas4  = Qms(4,Iset)
      Amas5  = Qms(5,Iset)
      Amas6  = Qms(6,Iset)

      Xmin   = Xmn  (Iset)
      Qini   = Qmn  (Iset)
      Qmax   = Qmx  (Iset)

      Do 201 Iexp = 1, Nexp(Iset)
 201     ExpNor(Iexp) = ExpN(Iexp, Iset)

      Ctq2df=0D0
      Return
C                         *************************
      END


C     Version 2 CTEQ distribution function in a parametrized form.

C By: J. Botts, J. Huston, H.L. Lai, J.G. Morfin, J.F. Owens, J. Qiu, W.K. Tung
C     & H. Weerts

C  To avoid the proliferation of parton distribution functions, we recommend 
C  that these distributions should replace Version 1 CTEQ distributions for 
C  general usage. The differences between the two sets of distributions, as 
C  briefly described below, do not significantly affect most applications for 
C  fixed-target and hadron-collider applications. 

C  Both CTEQ1 and CTEQ2 distributions fit existing DIS, Drell-Yan and Direct 
C  phton data with excellent chi-squares. They represent two distinct ways to  
C  resolve the inconsistency between CCFR and NMC measurements on F2 at small-x    
C  and the neutrino dimuon measurements of s(x), as discovered by the CTEQ1 
C  analysis: in the CTEQ1 analysis, the F2 measurements of CCFR and NMC are  
C  taken seriously, leaving out the dimuon information on s(x); whereas in the 
C  CTEQ2 analysis, thelatest CCFR NLO dimuon analysis of s(x) is used (within 
C  errors) as input butleaving out the small-x F2 data which conflict with this 
C  input. The small-xbehavior of the parton distributions are contrained in the 
C  latter case by thenewly released HERA data.  For details, see our 
C  forthcoming paper.

C     This file contains three versions of the same CTEQ2 parton distributions: 
C 
C Two "front-end" subprograms:    
C     FUNCTION Ctq2Pd (Iset, Iparton, X, Q, Irt) 
C         returns the PROBABILITY density for a GIVEN flavor;
C     SUBROUTINE Ctq2Pds(Iset, Pdf, XX, QQ, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C One lower-level subprogram:
C     FUNCTION Ctq2df (Iset, Iprtn, XX, QQ, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.
C Supplementary functions to return the relevant QCD parameters and other
C information concerning these distributions are also included (see below).      

C     Since this is an initial distribution of version 2, it is
C     useful for the authors to maintain a record of the distribution list in
C     case there are revisions or corrections.  
C     In the interest of maintaining the integrity of this package,
C     please do not freely distribute this program package; instead, refer any
C     interested colleagues to direct their request for a copy to:
C     Botts@hades.ifh.de or Lai@cteq11.pa.msu.edu

C If you have detailed questions concerning these CTEQ2 distributions, direct 
C inquires to Botts, Lai (see above) or Wu-Ki Tung (Tung@msupa.pa.msu.edu).

C     -------------------------------------------
C     Detailed instructions follow.

C     Name convention for CTEQ distributions:  CTEQnSx  where
C           n : version number                      (currently n = 2)
C           S : factorization scheme label: = [M D L] for [MS-bar DIS LO] 
c               resp.
C           x : special characteristics, if any
C            (e.g. S(F) for singular (flat) small-x, L for "LEP lambda value")

C    Explanation of functional arguments:

C    Iset is the set label; in this version, Iset = 1, 2, 3, 4, 5, 6 
C                           correspond to the following CTEQ global fits:

C          cteq2M  : best fit in the MS-bar scheme 
C          cteq2MS : singular small-x
C          cteq2MF : flat small-x
C          cteq2ML : large lambda (Lambda(5) = 220 MeV)

C          cteq2L  : best fit in Leading order QCD

C          cteq2D  : best fit in the DIS scheme

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)

C   X, Q are the usual x, Q; 
C   Irt is a return error code (see individual modules for explanation).
C       
C     ---------------------------------------------
C --> QCD parameters for parton distribution set Iset can be obtained inside
C         the user's program by:
C     Dum = Prctq2
C    >        (Iset, Iord, Ischeme, MxFlv,
C    >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
C    >         Xmin, Qini, Qmax, ExpNor)
C     where all but the first argument are output parameters.
C     They should be self-explanatory -- see details under ENTRY Prctq2.

C  Since the QCD Lambda value for the various sets are needed more often than
C  the other parameters in most applications, a special function
C     Wlamd2 (Iset, Iorder, Neff)                    is provided
C  which returns the lambda value for Neff = 4,5,6 effective flavors as well as
C  the order these values pertain to.

C     ----------------------------------------------
C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for fixed target
C     experiments and 0.0001 < x < 0.01 from first official data of HERA.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-5 < x < 1 ; 1.6 GeV < Q < 1 TeV.  The 
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   a threshold factor is applied to A0 which simulates the proper Q-dependence
C   of the QCD evolution in that region according to the renormalization
C   scheme defined in Collins-Tung, Nucl. Phys. B278, 934 (1986).

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.
C    ------------------------------------------------

      FUNCTION Ctq2Pd (Iset, Iparton, X, Q, Irt)
C                                                   -=-=- ctq2pd

C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C   --- the PROBABILITY density

C   (Iset, Iparton, X, Q): explained above;

C    Irt : return error code: see module Ctq2df for explanation.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Ifl = Iparton
      JFL = ABS(Ifl)
C                                                             Valence
      IF (Ifl.Eq.1 .or. Ifl.Eq.2) THEN
        VL = Ctq2df(Iset, Ifl, X, Q, Irt)
      ELSE
        VL = 0.
      ENDIF
C                                                             Sea
      SEA = Ctq2df (Iset, -JFL, X, Q, Irt)
C                                              Full (probability) Distribution
      Ctq2pd = (VL + SEA) / X
      
      Return
C                         *************************
      END
C
C
      SUBROUTINE Ctq2Pds(Iset, Pdf, X, Q, Irt)
C                                                   -=-=- ctq2pds

C   This function returns the CTEQ parton distributions xf^Iset_Iprtn/proton
C   --- the Momentum density in array form
c
C    (Iset, X, Q): explained in header comment lines;

C     Irt : return error code -- cumulated over flavors: 
C           see module Ctq2df for explanation on individual flavors.
C     Pdf (Iparton);  
C         Iparton = -6, -5, ...0, 1, 2 ... 6
C               has the same meaning as explained in the header comment lines.
    
      Implicit Double Precision (A-H, O-Z)
      Dimension Pdf (-6:6)

      Irt=0
      do 10 I=-6,2
         if(I.le.0) then
            Pdf(I) = Ctq2df(Iset,I,X,Q,Irt1)
            Pdf(-I)= Pdf(I)
         else
            Pdf(I) = Ctq2df(Iset,I,X,Q,Irt1) + Pdf(-I)
         endif
         Irt=Irt+Irt1
  10  Continue

      Return
C                         *************************
      End

      FUNCTION Ctq3df (Iset, Iprtn, XX, QQ, Irt)
C                                                   -=-=- ctq3df

C            Returns xf(x,Q) -- the momentum fraction distribution !!
C            Returns valence and sea rather than combined flavor distr.

C            Iset : PDF set label

C            Iprtn  : Parton label:   2, 1 = d_ and u_ valence
C                                     0 = gluon
C                            -1, ... -6 = u, d, s, c, b, t sea quarks

C            XX  : Bjorken-x
C            QQ  : scale parameter "Q"
C      Irt : Return code
C      0 : no error
C      1 : parametrization is slightly negative; reset to 0.0.
C          (This condition happens rarely -- only for large x where the 
C          absolute value of the parton distribution is extremely small.) 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      Parameter (Nst = 3)

      DIMENSION
     >   Iord(Nst), Isch(Nst), Nqrk(Nst),Alm(Nst)
     > , Vlm(4:6,Nst), Qms(4:6, Nst)
     > , Xmn(Nst), Qmn(Nst), Qmx(Nst)

c                                          --------- CTEQ3M
c
      DATA 
     >  Isch(1), Iord(1), Nqrk(1), Alm(1) /  1,  2,  6, .239  / 
     >  (Vlm(I,1), I=4,6) / .239,    .158,     .063   /
     >  (Qms(I,1), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(1), Qmn(1), Qmx(1) /  1.E-6,  1.60,  1.E4  /

c                                          --------- CTEQ3L
c
      DATA 
     >  Isch(2), Iord(2), Nqrk(2), Alm(2) /  1,  1,  6, .177  / 
     >  (Vlm(I,2), I=4,6) / .177,    .132,     .066   /
     >  (Qms(I,2), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(2), Qmn(2), Qmx(2) /  1.E-6,  1.60,  1.E4  /

c                                          --------- CTEQ3D
c
      DATA 
     >  Isch(3), Iord(3), Nqrk(3), Alm(3) /  1,  2,  6, .247  / 
     >  (Vlm(I,3), I=4,6) / .247,    .164,     .066   /
     >  (Qms(I,3), I=4,6) / 1.60,   5.00,  180.0 /
     >  Xmn(3), Qmn(3), Qmx(3) /  1.E-6,  1.60,  1.E4  /


      Data Ist, Lp, Qsto / 0, -10, 1.2345 /

      save Ist, Lp, Qsto
      save SB, SB2, SB3

      X  = XX
      Irt = 0
      if(Iset.eq.Ist .and. Qsto.eq.QQ) then
C                                             if only change is in x:
        if (Iprtn.eq.Lp) goto 100
C                         if change in flv is within "light" partons:
        if (Iprtn.ge.-3 .and. Lp.ge.-3) goto 501
      endif

      Ip = abs(Iprtn)
C                                                  Set up Qi for SB
      If (Ip .GE. 4) then
         If (QQ .LE. Qms(Ip, Iset)) Then
           Ctq3df = 0.0
           Return
         Endif
         Qi = Qms(ip, Iset)
      Else
         Qi = Qmn(Iset)
      Endif
C                   Use "standard lambda" of parametrization program
      Alam = Alm (Iset)

      SBL = LOG(QQ/Alam) / LOG(Qi/Alam)
      SB = LOG (SBL)
      SB2 = SB*SB
      SB3 = SB2*SB

 501  Iflv = 3 - Iprtn

      Goto (1,2,3, 311) Iset

 1    Goto(11,12,13,14,15,16,17,18,19)Iflv    
c   Ifl =   2
  11  A0=Exp(-0.7266E+00-0.1584E+01*SB +0.1259E+01*SB2-0.4305E-01*SB3)
      A1= 0.5285E+00-0.3721E+00*SB +0.5150E+00*SB2-0.1697E+00*SB3 
      A2= 0.4075E+01+0.8282E+00*SB -0.4496E+00*SB2+0.2107E+00*SB3 
      A3= 0.3279E+01+0.5066E+01*SB -0.9134E+01*SB2+0.2897E+01*SB3 
      A4= 0.4399E+00-0.5888E+00*SB +0.4802E+00*SB2-0.1664E+00*SB3 
      A5= 0.3678E+00-0.8929E+00*SB +0.1592E+01*SB2-0.5713E+00*SB3 
      goto 100
c   Ifl =   1
  12  A0=Exp( 0.2259E+00+0.1237E+00*SB +0.3035E+00*SB2-0.2935E+00*SB3)
      A1= 0.5085E+00+0.1651E-01*SB -0.3592E-01*SB2+0.2782E-01*SB3 
      A2= 0.3732E+01+0.4901E+00*SB +0.2218E+00*SB2-0.1116E+00*SB3 
      A3= 0.7011E+01-0.6620E+01*SB +0.2557E+01*SB2-0.1360E+00*SB3 
      A4= 0.8969E+00-0.2429E+00*SB +0.1811E+00*SB2-0.6888E-01*SB3 
      A5= 0.8636E-01+0.2558E+00*SB -0.3082E+00*SB2+0.2535E+00*SB3 
      goto 100
c   Ifl =   0
  13  A0=Exp(-0.2318E+00-0.9779E+00*SB -0.3783E+00*SB2+0.1037E-01*SB3)
      A1=-0.2916E+00+0.1754E+00*SB -0.1884E+00*SB2+0.6116E-01*SB3 
      A2= 0.5349E+01+0.7460E+00*SB +0.2319E+00*SB2-0.2622E+00*SB3 
      A3= 0.6920E+01-0.3454E+01*SB +0.2027E+01*SB2-0.7626E+00*SB3 
      A4= 0.1013E+01+0.1423E+00*SB -0.1798E+00*SB2+0.1872E-01*SB3 
      A5=-0.5465E-01+0.2303E+01*SB -0.9584E+00*SB2+0.3098E+00*SB3 
      goto 100
c   Ifl =  -1
  14  A0=Exp(-0.2906E+01-0.1069E+00*SB -0.1055E+01*SB2+0.2496E+00*SB3)
      A1=-0.2875E+00+0.6571E-01*SB -0.1987E-01*SB2-0.1800E-02*SB3 
      A2= 0.9854E+01-0.2715E+00*SB -0.7407E+00*SB2+0.2888E+00*SB3 
      A3= 0.1583E+02-0.7687E+01*SB +0.3428E+01*SB2-0.3327E+00*SB3 
      A4= 0.9763E+00+0.7599E-01*SB -0.2128E+00*SB2+0.6852E-01*SB3 
      A5=-0.8444E-02+0.9434E+00*SB +0.4152E+00*SB2-0.1481E+00*SB3 
      goto 100
c   Ifl =  -2
  15  A0=Exp(-0.2328E+01-0.3061E+01*SB +0.3620E+01*SB2-0.1602E+01*SB3)
      A1=-0.3358E+00+0.3198E+00*SB -0.4210E+00*SB2+0.1571E+00*SB3 
      A2= 0.8478E+01-0.3112E+01*SB +0.5243E+01*SB2-0.2255E+01*SB3 
      A3= 0.1971E+02+0.3389E+00*SB -0.5268E+01*SB2+0.2099E+01*SB3 
      A4= 0.1128E+01-0.4701E+00*SB +0.7779E+00*SB2-0.3506E+00*SB3 
      A5=-0.4708E+00+0.3341E+01*SB -0.3375E+01*SB2+0.1353E+01*SB3 
      goto 100
c   Ifl =  -3
  16  A0=Exp(-0.3780E+01+0.2499E+01*SB -0.4962E+01*SB2+0.1936E+01*SB3)
      A1=-0.2639E+00-0.1575E+00*SB +0.3584E+00*SB2-0.1646E+00*SB3 
      A2= 0.8082E+01+0.2794E+01*SB -0.5438E+01*SB2+0.2321E+01*SB3 
      A3= 0.1811E+02-0.2000E+02*SB +0.1951E+02*SB2-0.6904E+01*SB3 
      A4= 0.9822E+00+0.4972E+00*SB -0.8690E+00*SB2+0.3415E+00*SB3 
      A5= 0.1772E+00-0.6078E+00*SB +0.3341E+01*SB2-0.1473E+01*SB3 
      goto 100
c   Ifl =  -4
  17  A0=SB** 0.1122E+01*Exp(-0.4232E+01-0.1808E+01*SB +0.5348E+00*SB2)
      A1=-0.2824E+00+0.5846E+00*SB -0.7230E+00*SB2+0.2419E+00*SB3 
      A2= 0.5683E+01-0.2948E+01*SB +0.5916E+01*SB2-0.2560E+01*SB3 
      A3= 0.2051E+01+0.4795E+01*SB -0.4271E+01*SB2+0.4174E+00*SB3 
      A4= 0.1737E+00+0.1717E+01*SB -0.1978E+01*SB2+0.6643E+00*SB3 
      A5= 0.8689E+00+0.3500E+01*SB -0.3283E+01*SB2+0.1026E+01*SB3 
      goto 100
c   Ifl =  -5
  18  A0=SB** 0.9906E+00*Exp(-0.1496E+01-0.6576E+01*SB +0.1569E+01*SB2)
      A1=-0.2140E+00-0.6419E-01*SB -0.2741E-02*SB2+0.3185E-02*SB3 
      A2= 0.5781E+01+0.1049E+00*SB -0.3930E+00*SB2+0.5174E+00*SB3 
      A3=-0.9420E+00+0.5511E+00*SB +0.8817E+00*SB2+0.1903E+01*SB3 
      A4= 0.2418E-01+0.4232E-01*SB -0.1244E-01*SB2-0.2365E-01*SB3 
      A5= 0.7664E+00+0.1794E+01*SB -0.4917E+00*SB2-0.1284E+00*SB3 
      goto 100
c   Ifl =  -6
  19  A0=SB** 0.1000E+01*Exp(-0.8460E+01+0.1154E+01*SB +0.8838E+01*SB2)
      A1=-0.4316E-01-0.2976E+00*SB +0.3174E+00*SB2-0.1429E+01*SB3 
      A2= 0.4910E+01+0.2273E+01*SB +0.5631E+01*SB2-0.1994E+02*SB3 
      A3= 0.1190E+02-0.2000E+02*SB -0.2000E+02*SB2+0.1292E+02*SB3 
      A4= 0.5771E+00-0.2552E+00*SB +0.7510E+00*SB2+0.6923E+00*SB3 
      A5= 0.4402E+01-0.1627E+01*SB -0.2085E+01*SB2-0.6737E+01*SB3 
      goto 100

 2    Goto(21,22,23,24,25,26,27,28,29)Iflv    
c   Ifl =   2
  21  A0=Exp( 0.1141E+00+0.4764E+00*SB -0.1745E+01*SB2+0.7728E+00*SB3)
      A1= 0.4275E+00-0.1290E+00*SB +0.3609E+00*SB2-0.1689E+00*SB3 
      A2= 0.3000E+01+0.2946E+01*SB -0.4117E+01*SB2+0.1989E+01*SB3 
      A3=-0.1302E+01+0.2322E+01*SB -0.4258E+01*SB2+0.2109E+01*SB3 
      A4= 0.2586E+01-0.1920E+00*SB -0.3754E+00*SB2+0.2731E+00*SB3 
      A5=-0.2251E+00-0.5374E+00*SB +0.2245E+01*SB2-0.1034E+01*SB3 
      goto 100
c   Ifl =   1
  22  A0=Exp( 0.1907E+00+0.4205E-01*SB +0.2752E+00*SB2-0.3171E+00*SB3)
      A1= 0.4611E+00+0.2331E-01*SB -0.3403E-01*SB2+0.3174E-01*SB3 
      A2= 0.3504E+01+0.5739E+00*SB +0.2676E+00*SB2-0.1553E+00*SB3 
      A3= 0.7452E+01-0.6742E+01*SB +0.2849E+01*SB2-0.1964E+00*SB3 
      A4= 0.1116E+01-0.3435E+00*SB +0.2865E+00*SB2-0.1288E+00*SB3 
      A5= 0.6659E-01+0.2714E+00*SB -0.2688E+00*SB2+0.2763E+00*SB3 
      goto 100
c   Ifl =   0
  23  A0=Exp(-0.7631E+00-0.7241E+00*SB -0.1170E+01*SB2+0.5343E+00*SB3)
      A1=-0.3573E+00+0.3469E+00*SB -0.3396E+00*SB2+0.9188E-01*SB3 
      A2= 0.5604E+01+0.7458E+00*SB -0.5082E+00*SB2+0.1844E+00*SB3 
      A3= 0.1549E+02-0.1809E+02*SB +0.1162E+02*SB2-0.3483E+01*SB3 
      A4= 0.9881E+00+0.1364E+00*SB -0.4421E+00*SB2+0.2051E+00*SB3 
      A5=-0.9505E-01+0.3259E+01*SB -0.1547E+01*SB2+0.2918E+00*SB3 
      goto 100
c   Ifl =  -1
  24  A0=Exp(-0.2740E+01-0.7987E-01*SB -0.9015E+00*SB2-0.9872E-01*SB3)
      A1=-0.3909E+00+0.1244E+00*SB -0.4487E-01*SB2+0.1277E-01*SB3 
      A2= 0.9163E+01+0.2823E+00*SB -0.7720E+00*SB2-0.9360E-02*SB3 
      A3= 0.1080E+02-0.3915E+01*SB -0.1153E+01*SB2+0.2649E+01*SB3 
      A4= 0.9894E+00-0.1647E+00*SB -0.9426E-02*SB2+0.2945E-02*SB3 
      A5=-0.3395E+00+0.6998E+00*SB +0.7000E+00*SB2-0.6730E-01*SB3 
      goto 100
c   Ifl =  -2
  25  A0=Exp(-0.2449E+01-0.3513E+01*SB +0.4529E+01*SB2-0.2031E+01*SB3)
      A1=-0.4050E+00+0.3411E+00*SB -0.3669E+00*SB2+0.1109E+00*SB3 
      A2= 0.7470E+01-0.2982E+01*SB +0.5503E+01*SB2-0.2419E+01*SB3 
      A3= 0.1503E+02+0.1638E+01*SB -0.8772E+01*SB2+0.3852E+01*SB3 
      A4= 0.1137E+01-0.1006E+01*SB +0.1485E+01*SB2-0.6389E+00*SB3 
      A5=-0.5299E+00+0.3160E+01*SB -0.3104E+01*SB2+0.1219E+01*SB3 
      goto 100
c   Ifl =  -3
  26  A0=Exp(-0.3640E+01+0.1250E+01*SB -0.2914E+01*SB2+0.8390E+00*SB3)
      A1=-0.3595E+00-0.5259E-01*SB +0.3122E+00*SB2-0.1642E+00*SB3 
      A2= 0.7305E+01+0.9727E+00*SB -0.9788E+00*SB2-0.5193E-01*SB3 
      A3= 0.1198E+02-0.1799E+02*SB +0.2614E+02*SB2-0.1091E+02*SB3 
      A4= 0.9882E+00-0.6101E+00*SB +0.9737E+00*SB2-0.4935E+00*SB3 
      A5=-0.1186E+00-0.3231E+00*SB +0.3074E+01*SB2-0.1274E+01*SB3 
      goto 100
c   Ifl =  -4
  27  A0=SB** 0.1122E+01*Exp(-0.3718E+01-0.1335E+01*SB +0.1651E-01*SB2)
      A1=-0.4719E+00+0.7509E+00*SB -0.8420E+00*SB2+0.2901E+00*SB3 
      A2= 0.6194E+01-0.1641E+01*SB +0.4907E+01*SB2-0.2523E+01*SB3 
      A3= 0.4426E+01-0.4270E+01*SB +0.6581E+01*SB2-0.3474E+01*SB3 
      A4= 0.2683E+00+0.9876E+00*SB -0.7612E+00*SB2+0.1780E+00*SB3 
      A5=-0.4547E+00+0.4410E+01*SB -0.3712E+01*SB2+0.1245E+01*SB3 
      goto 100
c   Ifl =  -5
  28  A0=SB** 0.9838E+00*Exp(-0.2548E+01-0.7660E+01*SB +0.3702E+01*SB2)
      A1=-0.3122E+00-0.2120E+00*SB +0.5716E+00*SB2-0.3773E+00*SB3 
      A2= 0.6257E+01-0.8214E-01*SB -0.2537E+01*SB2+0.2981E+01*SB3 
      A3=-0.6723E+00+0.2131E+01*SB +0.9599E+01*SB2-0.7910E+01*SB3 
      A4= 0.9169E-01+0.4295E-01*SB -0.5017E+00*SB2+0.3811E+00*SB3 
      A5= 0.2402E+00+0.2656E+01*SB -0.1586E+01*SB2+0.2880E+00*SB3 
      goto 100
c   Ifl =  -6
  29  A0=SB** 0.1001E+01*Exp(-0.6934E+01+0.3050E+01*SB -0.6943E+00*SB2)
      A1=-0.1713E+00-0.5167E+00*SB +0.1241E+01*SB2-0.1703E+01*SB3 
      A2= 0.6169E+01+0.3023E+01*SB -0.1972E+02*SB2+0.1069E+02*SB3 
      A3= 0.4439E+01-0.1746E+02*SB +0.1225E+02*SB2+0.8350E+00*SB3 
      A4= 0.5458E+00-0.4586E+00*SB +0.9089E+00*SB2-0.4049E+00*SB3 
      A5= 0.3207E+01-0.3362E+01*SB +0.5877E+01*SB2-0.7659E+01*SB3 
      goto 100

 3    Goto(31,32,33,34,35,36,37,38,39)Iflv    
c   Ifl =   2
  31  A0=Exp( 0.3961E+00+0.4914E+00*SB -0.1728E+01*SB2+0.7257E+00*SB3)
      A1= 0.4162E+00-0.1419E+00*SB +0.3680E+00*SB2-0.1618E+00*SB3 
      A2= 0.3248E+01+0.3028E+01*SB -0.4307E+01*SB2+0.1920E+01*SB3 
      A3=-0.1100E+01+0.2184E+01*SB -0.3820E+01*SB2+0.1717E+01*SB3 
      A4= 0.2082E+01-0.2756E+00*SB +0.3043E+00*SB2-0.1260E+00*SB3 
      A5=-0.4822E+00-0.5706E+00*SB +0.2243E+01*SB2-0.9760E+00*SB3 
      goto 100
c   Ifl =   1
  32  A0=Exp( 0.2148E+00+0.5814E-01*SB +0.2734E+00*SB2-0.2902E+00*SB3)
      A1= 0.4810E+00+0.1657E-01*SB -0.3800E-01*SB2+0.3125E-01*SB3 
      A2= 0.3509E+01+0.3923E+00*SB +0.4010E+00*SB2-0.1932E+00*SB3 
      A3= 0.7055E+01-0.6552E+01*SB +0.3466E+01*SB2-0.5657E+00*SB3 
      A4= 0.1061E+01-0.3453E+00*SB +0.4089E+00*SB2-0.1817E+00*SB3 
      A5= 0.8687E-01+0.2548E+00*SB -0.2967E+00*SB2+0.2647E+00*SB3 
      goto 100
c   Ifl =   0
  33  A0=Exp(-0.4665E+00-0.7554E+00*SB -0.3323E+00*SB2-0.2734E-04*SB3)
      A1=-0.3359E+00+0.2395E+00*SB -0.2377E+00*SB2+0.7059E-01*SB3 
      A2= 0.5451E+01+0.6086E+00*SB +0.8606E-01*SB2-0.1425E+00*SB3 
      A3= 0.1026E+02-0.9352E+01*SB +0.4879E+01*SB2-0.1150E+01*SB3 
      A4= 0.9935E+00-0.5017E-01*SB -0.1707E-01*SB2-0.1464E-02*SB3 
      A5=-0.4160E-01+0.2305E+01*SB -0.1063E+01*SB2+0.3211E+00*SB3 
      goto 100
c   Ifl =  -1
  34  A0=Exp(-0.3323E+01+0.2296E+00*SB -0.1109E+01*SB2+0.2223E+00*SB3)
      A1=-0.3410E+00+0.8847E-01*SB -0.1111E-01*SB2-0.5927E-02*SB3 
      A2= 0.9753E+01-0.5182E+00*SB -0.4670E+00*SB2+0.1921E+00*SB3 
      A3= 0.1977E+02-0.1600E+02*SB +0.9481E+01*SB2-0.1864E+01*SB3 
      A4= 0.9818E+00+0.2839E-02*SB -0.1188E+00*SB2+0.3584E-01*SB3 
      A5=-0.7934E-01+0.1004E+01*SB +0.3704E+00*SB2-0.1220E+00*SB3 
      goto 100
c   Ifl =  -2
  35  A0=Exp(-0.2714E+01-0.2868E+01*SB +0.3700E+01*SB2-0.1671E+01*SB3)
      A1=-0.3893E+00+0.3341E+00*SB -0.3897E+00*SB2+0.1420E+00*SB3 
      A2= 0.8359E+01-0.3267E+01*SB +0.5327E+01*SB2-0.2245E+01*SB3 
      A3= 0.2359E+02-0.5669E+01*SB -0.4602E+01*SB2+0.3153E+01*SB3 
      A4= 0.1106E+01-0.4745E+00*SB +0.7739E+00*SB2-0.3417E+00*SB3 
      A5=-0.5557E+00+0.3433E+01*SB -0.3390E+01*SB2+0.1354E+01*SB3 
      goto 100
c   Ifl =  -3
  36  A0=Exp(-0.3985E+01+0.2855E+01*SB -0.5208E+01*SB2+0.1937E+01*SB3)
      A1=-0.3337E+00-0.1150E+00*SB +0.3691E+00*SB2-0.1709E+00*SB3 
      A2= 0.7968E+01+0.3641E+01*SB -0.6599E+01*SB2+0.2642E+01*SB3 
      A3= 0.1873E+02-0.1999E+02*SB +0.1734E+02*SB2-0.5813E+01*SB3 
      A4= 0.9731E+00+0.5082E+00*SB -0.8780E+00*SB2+0.3231E+00*SB3 
      A5=-0.5542E-01-0.4189E+00*SB +0.3309E+01*SB2-0.1439E+01*SB3 
      goto 100
c   Ifl =  -4
  37  A0=SB** 0.1105E+01*Exp(-0.3952E+01-0.1901E+01*SB +0.5137E+00*SB2)
      A1=-0.3543E+00+0.6055E+00*SB -0.6941E+00*SB2+0.2278E+00*SB3 
      A2= 0.5955E+01-0.2629E+01*SB +0.5337E+01*SB2-0.2300E+01*SB3 
      A3= 0.1933E+01+0.4882E+01*SB -0.3810E+01*SB2+0.2290E+00*SB3 
      A4= 0.1806E+00+0.1655E+01*SB -0.1893E+01*SB2+0.6395E+00*SB3 
      A5= 0.4790E+00+0.3612E+01*SB -0.3152E+01*SB2+0.9684E+00*SB3 
      goto 100
c   Ifl =  -5
  38  A0=SB** 0.9818E+00*Exp(-0.1825E+01-0.7464E+01*SB +0.2143E+01*SB2)
      A1=-0.2604E+00-0.1400E+00*SB +0.1702E+00*SB2-0.8476E-01*SB3 
      A2= 0.6005E+01+0.6275E+00*SB -0.2535E+01*SB2+0.2219E+01*SB3 
      A3=-0.9067E+00+0.1149E+01*SB +0.1974E+01*SB2+0.4716E+01*SB3 
      A4= 0.3915E-01+0.5945E-01*SB -0.9844E-01*SB2+0.2783E-01*SB3 
      A5= 0.5500E+00+0.1994E+01*SB -0.6727E+00*SB2-0.1510E+00*SB3 
      goto 100
c   Ifl =  -6
  39  A0=SB** 0.1002E+01*Exp(-0.8553E+01+0.3793E+00*SB +0.9998E+01*SB2)
      A1=-0.5870E-01-0.2792E+00*SB +0.6526E+00*SB2-0.1984E+01*SB3 
      A2= 0.4716E+01+0.4473E+00*SB +0.1128E+02*SB2-0.1937E+02*SB3 
      A3= 0.1289E+02-0.1742E+02*SB -0.1983E+02*SB2-0.9274E+00*SB3 
      A4= 0.5647E+00-0.2732E+00*SB +0.1074E+01*SB2+0.5981E+00*SB3 
      A5= 0.4390E+01-0.1262E+01*SB -0.9026E+00*SB2-0.9394E+01*SB3 
      goto 100

 311  stop 'This option is not currently supported.'

 100  Ctq3df = A0 *(x**A1) *((D1-x)**A2) *(D1+A3*(x**A4))
     $            *(log(D1+D1/x))**A5

      if(Ctq3df.lt.D0) then
        Ctq3df = D0
        Irt=1
      endif

      Ist = Iset

      Lp  = Iprtn
      Qsto = QQ

      Return
C                                  -----------------------
      ENTRY Wlamd3 (Iset, Iorder, Neff)

C     Returns the EFFECTIVE QCD lambda values for order=Iorder and
C     effective # of flavors = Neff for each of the PDF sets.

      Iorder = Iord (Iset)
      Wlamd3 = VLM  (Neff, Iset)

      RETURN

C                         *************************
      END

C     Version 3 CTEQ distribution function in a parametrized form.

C   By: H.L. Lai, J. Botts, J. Huston, J.G. Morfin, J.F. Owens, J. Qiu,
C       W.K. Tung & H. Weerts;  Preprint MSU-HEP/41024, CTEQ 404 

C   This file contains three versions of the same CTEQ3 parton distributions: 
C 
C Two "front-end" subprograms:    
C     FUNCTION Ctq3Pd (Iset, Iparton, X, Q, Irt) 
C         returns the PROBABILITY density for a GIVEN flavor;
C     SUBROUTINE Ctq3Pds(Iset, Pdf, XX, QQ, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C One lower-level subprogram:
C     FUNCTION Ctq3df (Iset, Iprtn, XX, QQ, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.

C      One supplementary function to return the QCD lambda parameter 
C      concerning these distributions is also included (see below). 

C     Although DOUBLE PRECISION is used, conversion to SINGLE PRECISION
C     is straightforward by removing the 
C     Implicit Double Precision statements. 

C     Since this is an initial distribution of version 3, it is
C     useful for the authors to maintain a record of the distribution
C     list in case there are revisions or corrections.
C     In the interest of maintaining the integrity of this package,
C     please do not freely distribute this program package; instead, refer
C     any interested colleagues to direct their request for a copy to:
C     Lai@cteq11.pa.msu.edu or Tung@msupa.pa.msu.edu.

C   If you have detailed questions concerning these CTEQ3 distributions, 
C   or if you find problems/bugs using this initial distribution, direct 
C   inquires to Hung-Liang Lai or Wu-Ki Tung.

C     -------------------------------------------
C     Detailed instructions follow.

C     Name convention for CTEQ distributions:  CTEQnSx  where
C        n : version number                      (currently n = 3)
C        S : factorization scheme label: = [M L D] for [MS-bar LO DIS] 
c               resp.
C        x : special characteristics, if any
C        (e.g. S(F) for singular (flat) small-x, L for "LEP lambda value")
C        (not applicable to CTEQ3 since only three standard sets are given.)

C    Explanation of functional arguments:

C    Iset is the set label; in this version, Iset = 1, 2, 3 
C                           correspond to the following CTEQ global fits:

C          cteq3M  : best fit in the MS-bar scheme 
C          cteq3L  : best fit in Leading order QCD
C          cteq3D  : best fit in the DIS scheme

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)
C  *** WARNING: We use the parton label 2 as D-quark, and 1 as U-quark which 
C               might be different with your labels.

C   X, Q are the usual x, Q; 
C   Irt is a return error code (see individual modules for explanation).
C       
C     ---------------------------------------------

C  Since the QCD Lambda value for the various sets are needed more often than
C  the other parameters in most applications, a special function
C     Wlamd3 (Iset, Iorder, Neff)                    is provided
C  which returns the lambda value for Neff = 4,5,6 effective flavors as well as
C  the order these values pertain to.

C     ----------------------------------------------
C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for fixed target
C     experiments and 0.0001 < x < 0.1 from HERA data.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-6 < x < 1 ; 1.6 GeV < Q < 10 TeV.  The 
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   a threshold factor is applied to A0 which simulates the proper Q-dependence
C   of the QCD evolution in that region according to the renormalization
C   scheme defined in Collins-Tung, Nucl. Phys. B278, 934 (1986).

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.
C    ------------------------------------------------

      FUNCTION Ctq3Pd (Iset, Iparton, X, Q, Irt)
C                                                   -=-=- ctq3pd

C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C   --- the PROBABILITY density

C   (Iset, Iparton, X, Q): explained above;

C    Irt : return error code: see module Ctq3df for explanation.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Ifl = Iparton
      JFL = ABS(Ifl)
C                                                             Valence
      IF (Ifl.Eq.1 .or. Ifl.Eq.2) THEN
        VL = Ctq3df(Iset, Ifl, X, Q, Irt)
      ELSE
        VL = 0.
      ENDIF
C                                                             Sea
      SEA = Ctq3df (Iset, -JFL, X, Q, Irt)
C                                              Full (probability) Distribution
      Ctq3pd = (VL + SEA) / X
      
      Return
C                         *************************
      END
C
C
      SUBROUTINE Ctq3Pds(Iset, Pdf, X, Q, Irt)
C                                                   -=-=- ctq3pds

C   This function returns the CTEQ parton distributions xf^Iset_Iprtn/proton
C   --- the Momentum density in array form
c
C    (Iset, X, Q): explained in header comment lines;

C     Irt : return error code -- cumulated over flavors: 
C           see module Ctq3df for explanation on individual flavors.
C     Pdf (Iparton);  
C         Iparton = -6, -5, ...0, 1, 2 ... 6
C               has the same meaning as explained in the header comment lines.
    
      Implicit Double Precision (A-H, O-Z)
      Dimension Pdf (-6:6)

      Irt=0
      do 10 I=-6,2
         if(I.le.0) then
            Pdf(I) = Ctq3df(Iset,I,X,Q,Irt1)
            Pdf(-I)= Pdf(I)
         else
            Pdf(I) = Ctq3df(Iset,I,X,Q,Irt1) + Pdf(-I)
         endif
         Irt=Irt+Irt1
  10  Continue

      Return
C                         *************************
      End

      Function Ctq4Pdf (Iparton, X, Q)
C                                                   -=-=- ctq4pdf

C                CTEQ Parton Distribution Functions: Version 4.6
C                             June 21, 1996
C                   Modified: 10/17/96, 1/7/97, 1/15/97
C                             2/17/97, 2/21/97
C                   Last Modified on April 2, 1997
C
C   Ref[1]: "IMPROVED PARTON DISTRIBUTIONS FROM GLOBAL ANALYSIS OF RECENT DEEP
C         INELASTIC SCATTERING AND INCLUSIVE JET DATA"
C   By: H.L. Lai, J. Huston, S. Kuhlmann, F. Olness, J. Owens, D. Soper
C       W.K. Tung, H. Weerts
C       Phys. Rev. D55, 1280 (1997)
C
C   Ref[2]: "CHARM PRODUCTION AND PARTON DISTRIBUTIONS"
C   By: H.L. Lai and W.K. Tung
C       Z. Phys. C74, 463 (1997)
C
C   This package contains 13 sets of CTEQ4 PDF's. Details are:
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5   Table_File
C ---------------------------------------------------------------------------
C Ref[1]
C   1    CTEQ4M   Standard MSbar scheme   0.116     298   202    cteq4m.tbl
C   2    CTEQ4D   Standard DIS scheme     0.116     298   202    cteq4d.tbl
C   3    CTEQ4L   Leading Order           0.132     236   181    cteq4l.tbl
C   4    CTEQ4A1  Alpha_s series          0.110     215   140    cteq4a1.tbl
C   5    CTEQ4A2  Alpha_s series          0.113     254   169    cteq4a2.tbl
C   6    CTEQ4A3            ( same as CTEQ4M )
C   7    CTEQ4A4  Alpha_s series          0.119     346   239    cteq4a4.tbl
C   8    CTEQ4A5  Alpha_s series          0.122     401   282    cteq4a5.tbl
C   9    CTEQ4HJ  High Jet                0.116     303   206    cteq4hj.tbl
C   10   CTEQ4LQ  Low Q0                  0.114     261   174    cteq4lq.tbl
C ---------------------------------------------------------------------------
C Ref[2]
C   11   CTEQ4HQ  Heavy Quark             0.116     298   202    cteq4hq.tbl
C   12   CTEQ4HQ1 Heavy Quark:Q0=1,Mc=1.3 0.116     298   202    cteq4hq1.tbl
C        (Improved version of CTEQ4HQ, recommended)
C   13   CTEQ4F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=385)   cteq4f3.tbl
C   14   CTEQ4F4  Nf=4 FixedFlavorNumber  0.111     292   XXX    cteq4f4.tbl
C ---------------------------------------------------------------------------
C   
C   The available applied range is 10^-5 < x < 1 and 1.6 < Q < 10,000 (GeV) 
C   except CTEQ4LQ(4HQ1) for which Q starts at a lower value of 0.7(1.0) GeV.  
C   Lam5 (Lam4, Lam3) represents Lambda value (in MeV) for 5 (4,3) flavors. 
C   The matching alpha_s between 4 and 5 flavors takes place at Q=5.0 GeV,  
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C   
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq4(Iset) 
C   where Iset is the desired PDF specified in the above table.
C   
C   The function Ctq4Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton] 
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C      whereas CTEQ4F3 has, by definition, only 3 flavors and gluon;
C              CTEQ4F4 has only 4 flavors and gluon.
C   
C   For detailed information on the parameters used, e.q. quark masses, 
C   QCD Lambda, ... etc.,  see info lines at the beginning of the 
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single 
C   precision.
C   
C   If you have detailed questions concerning these CTEQ4 distributions, 
C   or if you find problems/bugs using this package, direct inquires to 
C   Hung-Liang Lai(Lai_H@pa.msu.edu) or Wu-Ki Tung(Tung@pa.msu.edu).
C   
C===========================================================================
C      Function Ctq4Pdf (Iparton, X, Q)

      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq4Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq4Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq4Pdf: '
     >              , Iparton
         Endif
         Ctq4Pdf = 0D0
         Return
      Endif

      Ctq4Pdf = PartonX (Iparton, X, Q)
      if(Ctq4Pdf.lt.0.D0)  Ctq4Pdf = 0.D0

      Return

C                             ********************
      End

c ---------------------------------------------------------------------
	subroutine ctq5ini (iset)
c ---------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	parameter (isetmax=7)

      Dimension Iorder(isetmax), Nflv(isetmax), Alambda(isetmax)

      data 
     > (Iorder(i), i=1,isetmax) / 2, 2, 1, 2, 2, 2, 2 /
     > (Nflv(i),   i=1,isetmax) / 5, 5, 5, 5, 5, 3, 4 /
     > (Alambda(i), i=1,isetmax)
     > / 0.226, 0.226, 0.146, 0.226, 0.226, 0.395, 0.325 /

	if((iset .ge. 1) .and. (iset .le. 7)) then
	   print *, Nflv(iset), Alambda(iset), Iorder(iset)
	   Call SetLam(Nflv(iset), Alambda(iset), Iorder(iset))
	   return
	else
	   print *,'undefined iset in ctq5ini: ', iset
	   stop
	endif

	end
c --------------------------------------------------------------------------
	double precision function ctq5L(ifl,x,q)
c Parametrization of cteq5L parton distribution functions (J. Pumplin 9/99).
c ifl: 1=u,2=d,3=s,4=c,5=b;0=gluon;-1=ubar,-2=dbar,-3=sbar,-4=cbar,-5=bbar.
c --------------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	integer ifl

	ii = ifl
	if(ii .gt. 2) then
	   ii = -ii
	endif

	if(ii .eq. -1) then
	   sum = faux5L(-1,x,q)
	   ratio = faux5L(-2,x,q)
	   ctq5L = sum/(1.d0 + ratio)

	elseif(ii .eq. -2) then
	   sum = faux5L(-1,x,q)
	   ratio = faux5L(-2,x,q)
	   ctq5L = sum*ratio/(1.d0 + ratio)

	elseif(ii .ge. -5) then
	   ctq5L = faux5L(ii,x,q)

	else
	   ctq5L = 0.d0 

	endif

	return
	end

c ---------------------------------------------------------------------
      double precision function faux5L(ifl,x,q)
c auxiliary function for parametrization of CTEQ5L (J. Pumplin 9/99).
c ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ifl

      parameter (nex=8, nlf=2)
      dimension am(0:nex,0:nlf,-5:2)
      dimension alfvec(-5:2), qmavec(-5:2)
      dimension mexvec(-5:2), mlfvec(-5:2)
      dimension ut1vec(-5:2), ut2vec(-5:2)
      dimension af(0:nex)



      data mexvec( 2) / 8 /
      data mlfvec( 2) / 2 /
      data ut1vec( 2) /  0.4971265E+01 /
      data ut2vec( 2) / -0.1105128E+01 /
      data alfvec( 2) /  0.2987216E+00 /
      data qmavec( 2) /  0.0000000E+00 /
      data (am( 0,k, 2),k=0, 2)
     & /  0.5292616E+01, -0.2751910E+01, -0.2488990E+01 /
      data (am( 1,k, 2),k=0, 2)
     & /  0.9714424E+00,  0.1011827E-01, -0.1023660E-01 /
      data (am( 2,k, 2),k=0, 2)
     & / -0.1651006E+02,  0.7959721E+01,  0.8810563E+01 /
      data (am( 3,k, 2),k=0, 2)
     & / -0.1643394E+02,  0.5892854E+01,  0.9348874E+01 /
      data (am( 4,k, 2),k=0, 2)
     & /  0.3067422E+02,  0.4235796E+01, -0.5112136E+00 /
      data (am( 5,k, 2),k=0, 2)
     & /  0.2352526E+02, -0.5305168E+01, -0.1169174E+02 /
      data (am( 6,k, 2),k=0, 2)
     & / -0.1095451E+02,  0.3006577E+01,  0.5638136E+01 /
      data (am( 7,k, 2),k=0, 2)
     & / -0.1172251E+02, -0.2183624E+01,  0.4955794E+01 /
      data (am( 8,k, 2),k=0, 2)
     & /  0.1662533E-01,  0.7622870E-02, -0.4895887E-03 /

      data mexvec( 1) / 8 /
      data mlfvec( 1) / 2 /
      data ut1vec( 1) /  0.2612618E+01 /
      data ut2vec( 1) / -0.1258304E+06 /
      data alfvec( 1) /  0.3407552E+00 /
      data qmavec( 1) /  0.0000000E+00 /
      data (am( 0,k, 1),k=0, 2)
     & /  0.9905300E+00, -0.4502235E+00,  0.1624441E+00 /
      data (am( 1,k, 1),k=0, 2)
     & /  0.8867534E+00,  0.1630829E-01, -0.4049085E-01 /
      data (am( 2,k, 1),k=0, 2)
     & /  0.8547974E+00,  0.3336301E+00,  0.1371388E+00 /
      data (am( 3,k, 1),k=0, 2)
     & /  0.2941113E+00, -0.1527905E+01,  0.2331879E+00 /
      data (am( 4,k, 1),k=0, 2)
     & /  0.3384235E+02,  0.3715315E+01,  0.8276930E+00 /
      data (am( 5,k, 1),k=0, 2)
     & /  0.6230115E+01,  0.3134639E+01, -0.1729099E+01 /
      data (am( 6,k, 1),k=0, 2)
     & / -0.1186928E+01, -0.3282460E+00,  0.1052020E+00 /
      data (am( 7,k, 1),k=0, 2)
     & / -0.8545702E+01, -0.6247947E+01,  0.3692561E+01 /
      data (am( 8,k, 1),k=0, 2)
     & /  0.1724598E-01,  0.7120465E-02,  0.4003646E-04 /

      data mexvec( 0) / 8 /
      data mlfvec( 0) / 2 /
      data ut1vec( 0) / -0.4656819E+00 /
      data ut2vec( 0) / -0.2742390E+03 /
      data alfvec( 0) /  0.4491863E+00 /
      data qmavec( 0) /  0.0000000E+00 /
      data (am( 0,k, 0),k=0, 2)
     & /  0.1193572E+03, -0.3886845E+01, -0.1133965E+01 /
      data (am( 1,k, 0),k=0, 2)
     & / -0.9421449E+02,  0.3995885E+01,  0.1607363E+01 /
      data (am( 2,k, 0),k=0, 2)
     & /  0.4206383E+01,  0.2485954E+00,  0.2497468E+00 /
      data (am( 3,k, 0),k=0, 2)
     & /  0.1210557E+03, -0.3015765E+01, -0.1423651E+01 /
      data (am( 4,k, 0),k=0, 2)
     & / -0.1013897E+03, -0.7113478E+00,  0.2621865E+00 /
      data (am( 5,k, 0),k=0, 2)
     & / -0.1312404E+01, -0.9297691E+00, -0.1562531E+00 /
      data (am( 6,k, 0),k=0, 2)
     & /  0.1627137E+01,  0.4954111E+00, -0.6387009E+00 /
      data (am( 7,k, 0),k=0, 2)
     & /  0.1537698E+00, -0.2487878E+00,  0.8305947E+00 /
      data (am( 8,k, 0),k=0, 2)
     & /  0.2496448E-01,  0.2457823E-02,  0.8234276E-03 /

      data mexvec(-1) / 8 /
      data mlfvec(-1) / 2 /
      data ut1vec(-1) /  0.3862583E+01 /
      data ut2vec(-1) / -0.1265969E+01 /
      data alfvec(-1) /  0.2457668E+00 /
      data qmavec(-1) /  0.0000000E+00 /
      data (am( 0,k,-1),k=0, 2)
     & /  0.2647441E+02,  0.1059277E+02, -0.9176654E+00 /
      data (am( 1,k,-1),k=0, 2)
     & /  0.1990636E+01,  0.8558918E-01,  0.4248667E-01 /
      data (am( 2,k,-1),k=0, 2)
     & / -0.1476095E+02, -0.3276255E+02,  0.1558110E+01 /
      data (am( 3,k,-1),k=0, 2)
     & / -0.2966889E+01, -0.3649037E+02,  0.1195914E+01 /
      data (am( 4,k,-1),k=0, 2)
     & / -0.1000519E+03, -0.2464635E+01,  0.1964849E+00 /
      data (am( 5,k,-1),k=0, 2)
     & /  0.3718331E+02,  0.4700389E+02, -0.2772142E+01 /
      data (am( 6,k,-1),k=0, 2)
     & / -0.1872722E+02, -0.2291189E+02,  0.1089052E+01 /
      data (am( 7,k,-1),k=0, 2)
     & / -0.1628146E+02, -0.1823993E+02,  0.2537369E+01 /
      data (am( 8,k,-1),k=0, 2)
     & / -0.1156300E+01, -0.1280495E+00,  0.5153245E-01 /

      data mexvec(-2) / 7 /
      data mlfvec(-2) / 2 /
      data ut1vec(-2) /  0.1895615E+00 /
      data ut2vec(-2) / -0.3069097E+01 /
      data alfvec(-2) /  0.5293999E+00 /
      data qmavec(-2) /  0.0000000E+00 /
      data (am( 0,k,-2),k=0, 2)
     & / -0.6556775E+00,  0.2490190E+00,  0.3966485E-01 /
      data (am( 1,k,-2),k=0, 2)
     & /  0.1305102E+01, -0.1188925E+00, -0.4600870E-02 /
      data (am( 2,k,-2),k=0, 2)
     & / -0.2371436E+01,  0.3566814E+00, -0.2834683E+00 /
      data (am( 3,k,-2),k=0, 2)
     & / -0.6152826E+01,  0.8339877E+00, -0.7233230E+00 /
      data (am( 4,k,-2),k=0, 2)
     & / -0.8346558E+01,  0.2892168E+01,  0.2137099E+00 /
      data (am( 5,k,-2),k=0, 2)
     & /  0.1279530E+02,  0.1021114E+00,  0.5787439E+00 /
      data (am( 6,k,-2),k=0, 2)
     & /  0.5858816E+00, -0.1940375E+01, -0.4029269E+00 /
      data (am( 7,k,-2),k=0, 2)
     & / -0.2795725E+02, -0.5263392E+00,  0.1290229E+01 /

      data mexvec(-3) / 7 /
      data mlfvec(-3) / 2 /
      data ut1vec(-3) /  0.3753257E+01 /
      data ut2vec(-3) / -0.1113085E+01 /
      data alfvec(-3) /  0.3713141E+00 /
      data qmavec(-3) /  0.0000000E+00 /
      data (am( 0,k,-3),k=0, 2)
     & /  0.1580931E+01, -0.2273826E+01, -0.1822245E+01 /
      data (am( 1,k,-3),k=0, 2)
     & /  0.2702644E+01,  0.6763243E+00,  0.7231586E-02 /
      data (am( 2,k,-3),k=0, 2)
     & / -0.1857924E+02,  0.3907500E+01,  0.5850109E+01 /
      data (am( 3,k,-3),k=0, 2)
     & / -0.3044793E+02,  0.2639332E+01,  0.5566644E+01 /
      data (am( 4,k,-3),k=0, 2)
     & / -0.4258011E+01, -0.5429244E+01,  0.4418946E+00 /
      data (am( 5,k,-3),k=0, 2)
     & /  0.3465259E+02, -0.5532604E+01, -0.4904153E+01 /
      data (am( 6,k,-3),k=0, 2)
     & / -0.1658858E+02,  0.2923275E+01,  0.2266286E+01 /
      data (am( 7,k,-3),k=0, 2)
     & / -0.1149263E+02,  0.2877475E+01, -0.7999105E+00 /

      data mexvec(-4) / 7 /
      data mlfvec(-4) / 2 /
      data ut1vec(-4) /  0.4400772E+01 /
      data ut2vec(-4) / -0.1356116E+01 /
      data alfvec(-4) /  0.3712017E-01 /
      data qmavec(-4) /  0.1300000E+01 /
      data (am( 0,k,-4),k=0, 2)
     & / -0.8293661E+00, -0.3982375E+01, -0.6494283E-01 /
      data (am( 1,k,-4),k=0, 2)
     & /  0.2754618E+01,  0.8338636E+00, -0.6885160E-01 /
      data (am( 2,k,-4),k=0, 2)
     & / -0.1657987E+02,  0.1439143E+02, -0.6887240E+00 /
      data (am( 3,k,-4),k=0, 2)
     & / -0.2800703E+02,  0.1535966E+02, -0.7377693E+00 /
      data (am( 4,k,-4),k=0, 2)
     & / -0.6460216E+01, -0.4783019E+01,  0.4913297E+00 /
      data (am( 5,k,-4),k=0, 2)
     & /  0.3141830E+02, -0.3178031E+02,  0.7136013E+01 /
      data (am( 6,k,-4),k=0, 2)
     & / -0.1802509E+02,  0.1862163E+02, -0.4632843E+01 /
      data (am( 7,k,-4),k=0, 2)
     & / -0.1240412E+02,  0.2565386E+02, -0.1066570E+02 /

      data mexvec(-5) / 6 /
      data mlfvec(-5) / 2 /
      data ut1vec(-5) /  0.5562568E+01 /
      data ut2vec(-5) / -0.1801317E+01 /
      data alfvec(-5) /  0.4952010E-02 /
      data qmavec(-5) /  0.4500000E+01 /
      data (am( 0,k,-5),k=0, 2)
     & / -0.6031237E+01,  0.1992727E+01, -0.1076331E+01 /
      data (am( 1,k,-5),k=0, 2)
     & /  0.2933912E+01,  0.5839674E+00,  0.7509435E-01 /
      data (am( 2,k,-5),k=0, 2)
     & / -0.8284919E+01,  0.1488593E+01, -0.8251678E+00 /
      data (am( 3,k,-5),k=0, 2)
     & / -0.1925986E+02,  0.2805753E+01, -0.3015446E+01 /
      data (am( 4,k,-5),k=0, 2)
     & / -0.9480483E+01, -0.9767837E+00, -0.1165544E+01 /
      data (am( 5,k,-5),k=0, 2)
     & /  0.2193195E+02, -0.1788518E+02,  0.9460908E+01 /
      data (am( 6,k,-5),k=0, 2)
     & / -0.1327377E+02,  0.1201754E+02, -0.6277844E+01 /



      if(q .le. qmavec(ifl)) then
         faux5L = 0.d0
         return
      endif

      if(x .ge. 1.d0) then
         faux5L = 0.d0
         return
      endif

      tmp = log(q/alfvec(ifl))
      if(tmp .le. 0.d0) then
         faux5L = 0.d0
         return
      endif

      sb = log(tmp)
      sb1 = sb - 1.2d0
      sb2 = sb1*sb1

      do i = 0, nex
         af(i) = 0.d0
         sbx = 1.d0
         do k = 0, mlfvec(ifl)
            af(i) = af(i) + sbx*am(i,k,ifl)
            sbx = sb1*sbx
         enddo
      enddo

      y = -log(x)
      u = log(x/0.00001d0)

      part1 = af(1)*y**(1.d0+0.01d0*af(4))*(1.d0+ af(8)*u)
      part2 = af(0)*(1.d0 - x) + af(3)*x 
      part3 = x*(1.d0-x)*(af(5)+af(6)*(1.d0-x)+af(7)*x*(1.d0-x))
      part4 = ut1vec(ifl)*log(1.d0-x) + 
     &	      AF(2)*log(1.d0+exp(ut2vec(ifl))-x)

      faux5L = exp(log(x) + part1 + part2 + part3 + part4)

c include threshold factor...
      faux5L = faux5L * (1.d0 - qmavec(ifl)/q)

      return
      end
c --------------------------------------------------------------------------
	double precision function ctq5MI(ifl,x,q)
c Parametrization of cteq5MI parton distribution functions (J. Pumplin 9/99).
c ifl: 1=u,2=d,3=s,4=c,5=b;0=gluon;-1=ubar,-2=dbar,-3=sbar,-4=cbar,-5=bbar.
c --------------------------------------------------------------------------
	implicit double precision (a-h,o-z)
	integer ifl

	ii = ifl
	if(ii .gt. 2) then
	   ii = -ii
	endif

	if(ii .eq. -1) then
	   sum = faux5MI(-1,x,q)
	   ratio = faux5MI(-2,x,q)
	   ctq5MI = sum/(1.d0 + ratio)

	elseif(ii .eq. -2) then
	   sum = faux5MI(-1,x,q)
	   ratio = faux5MI(-2,x,q)
	   ctq5MI = sum*ratio/(1.d0 + ratio)

	elseif(ii .ge. -5) then
	   ctq5MI = faux5MI(ii,x,q)

	else
	   ctq5MI = 0.d0 

	endif

	return
	end

c ---------------------------------------------------------------------
      double precision function faux5MI(ifl,x,q)
c auxiliary function for parametrization of CTEQ5MI (J. Pumplin 9/99).
c ---------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ifl

      parameter (nex=8, nlf=2)
      dimension am(0:nex,0:nlf,-5:2)
      dimension alfvec(-5:2), qmavec(-5:2)
      dimension mexvec(-5:2), mlfvec(-5:2)
      dimension ut1vec(-5:2), ut2vec(-5:2)
      dimension af(0:nex)




      data mexvec( 2) / 8 /
      data mlfvec( 2) / 2 /
      data ut1vec( 2) /  0.5141718E+01 /
      data ut2vec( 2) / -0.1346944E+01 /
      data alfvec( 2) /  0.5260555E+00 /
      data qmavec( 2) /  0.0000000E+00 /
      data (am( 0,k, 2),k=0, 2)
     & /  0.4289071E+01, -0.2536870E+01, -0.1259948E+01 /
      data (am( 1,k, 2),k=0, 2)
     & /  0.9839410E+00,  0.4168426E-01, -0.5018952E-01 /
      data (am( 2,k, 2),k=0, 2)
     & / -0.1651961E+02,  0.9246261E+01,  0.5996400E+01 /
      data (am( 3,k, 2),k=0, 2)
     & / -0.2077936E+02,  0.9786469E+01,  0.7656465E+01 /
      data (am( 4,k, 2),k=0, 2)
     & /  0.3054926E+02,  0.1889536E+01,  0.1380541E+01 /
      data (am( 5,k, 2),k=0, 2)
     & /  0.3084695E+02, -0.1212303E+02, -0.1053551E+02 /
      data (am( 6,k, 2),k=0, 2)
     & / -0.1426778E+02,  0.6239537E+01,  0.5254819E+01 /
      data (am( 7,k, 2),k=0, 2)
     & / -0.1909811E+02,  0.3695678E+01,  0.5495729E+01 /
      data (am( 8,k, 2),k=0, 2)
     & /  0.1889751E-01,  0.5027193E-02,  0.6624896E-03 /

      data mexvec( 1) / 8 /
      data mlfvec( 1) / 2 /
      data ut1vec( 1) /  0.4138426E+01 /
      data ut2vec( 1) / -0.3221374E+01 /
      data alfvec( 1) /  0.4960962E+00 /
      data qmavec( 1) /  0.0000000E+00 /
      data (am( 0,k, 1),k=0, 2)
     & /  0.1332497E+01, -0.3703718E+00,  0.1288638E+00 /
      data (am( 1,k, 1),k=0, 2)
     & /  0.7544687E+00,  0.3255075E-01, -0.4706680E-01 /
      data (am( 2,k, 1),k=0, 2)
     & / -0.7638814E+00,  0.5008313E+00, -0.9237374E-01 /
      data (am( 3,k, 1),k=0, 2)
     & / -0.3689889E+00, -0.1055098E+01, -0.4645065E+00 /
      data (am( 4,k, 1),k=0, 2)
     & /  0.3991610E+02,  0.1979881E+01,  0.1775814E+01 /
      data (am( 5,k, 1),k=0, 2)
     & /  0.6201080E+01,  0.2046288E+01,  0.3804571E+00 /
      data (am( 6,k, 1),k=0, 2)
     & / -0.8027900E+00, -0.7011688E+00, -0.8049612E+00 /
      data (am( 7,k, 1),k=0, 2)
     & / -0.8631305E+01, -0.3981200E+01,  0.6970153E+00 /
      data (am( 8,k, 1),k=0, 2)
     & /  0.2371230E-01,  0.5372683E-02,  0.1118701E-02 /

      data mexvec( 0) / 8 /
      data mlfvec( 0) / 2 /
      data ut1vec( 0) / -0.1026789E+01 /
      data ut2vec( 0) / -0.9051707E+01 /
      data alfvec( 0) /  0.9462977E+00 /
      data qmavec( 0) /  0.0000000E+00 /
      data (am( 0,k, 0),k=0, 2)
     & /  0.1191990E+03, -0.8548739E+00, -0.1963040E+01 /
      data (am( 1,k, 0),k=0, 2)
     & / -0.9449972E+02,  0.1074771E+01,  0.2056055E+01 /
      data (am( 2,k, 0),k=0, 2)
     & /  0.3701064E+01, -0.1167947E-02,  0.1933573E+00 /
      data (am( 3,k, 0),k=0, 2)
     & /  0.1171345E+03, -0.1064540E+01, -0.1875312E+01 /
      data (am( 4,k, 0),k=0, 2)
     & / -0.1014453E+03, -0.5707427E+00,  0.4511242E-01 /
      data (am( 5,k, 0),k=0, 2)
     & /  0.6365168E+01,  0.1275354E+01, -0.4964081E+00 /
      data (am( 6,k, 0),k=0, 2)
     & / -0.3370693E+01, -0.1122020E+01,  0.5947751E-01 /
      data (am( 7,k, 0),k=0, 2)
     & / -0.5327270E+01, -0.9293556E+00,  0.6629940E+00 /
      data (am( 8,k, 0),k=0, 2)
     & /  0.2437513E-01,  0.1600939E-02,  0.6855336E-03 /

      data mexvec(-1) / 8 /
      data mlfvec(-1) / 2 /
      data ut1vec(-1) /  0.5243571E+01 /
      data ut2vec(-1) / -0.2870513E+01 /
      data alfvec(-1) /  0.6701448E+00 /
      data qmavec(-1) /  0.0000000E+00 /
      data (am( 0,k,-1),k=0, 2)
     & /  0.2428863E+02,  0.1907035E+01, -0.4606457E+00 /
      data (am( 1,k,-1),k=0, 2)
     & /  0.2006810E+01, -0.1265915E+00,  0.7153556E-02 /
      data (am( 2,k,-1),k=0, 2)
     & / -0.1884546E+02, -0.2339471E+01,  0.5740679E+01 /
      data (am( 3,k,-1),k=0, 2)
     & / -0.2527892E+02, -0.2044124E+01,  0.1280470E+02 /
      data (am( 4,k,-1),k=0, 2)
     & / -0.1013824E+03, -0.1594199E+01,  0.2216401E+00 /
      data (am( 5,k,-1),k=0, 2)
     & /  0.8070930E+02,  0.1792072E+01, -0.2164364E+02 /
      data (am( 6,k,-1),k=0, 2)
     & / -0.4641050E+02,  0.1977338E+00,  0.1273014E+02 /
      data (am( 7,k,-1),k=0, 2)
     & / -0.3910568E+02,  0.1719632E+01,  0.1086525E+02 /
      data (am( 8,k,-1),k=0, 2)
     & / -0.1185496E+01, -0.1905847E+00, -0.8744118E-03 /

      data mexvec(-2) / 7 /
      data mlfvec(-2) / 2 /
      data ut1vec(-2) /  0.4782210E+01 /
      data ut2vec(-2) / -0.1976856E+02 /
      data alfvec(-2) /  0.7558374E+00 /
      data qmavec(-2) /  0.0000000E+00 /
      data (am( 0,k,-2),k=0, 2)
     & / -0.6216935E+00,  0.2369963E+00, -0.7909949E-02 /
      data (am( 1,k,-2),k=0, 2)
     & /  0.1245440E+01, -0.1031510E+00,  0.4916523E-02 /
      data (am( 2,k,-2),k=0, 2)
     & / -0.7060824E+01, -0.3875283E-01,  0.1784981E+00 /
      data (am( 3,k,-2),k=0, 2)
     & / -0.7430595E+01,  0.1964572E+00, -0.1284999E+00 /
      data (am( 4,k,-2),k=0, 2)
     & / -0.6897810E+01,  0.2620543E+01,  0.8012553E-02 /
      data (am( 5,k,-2),k=0, 2)
     & /  0.1507713E+02,  0.2340307E-01,  0.2482535E+01 /
      data (am( 6,k,-2),k=0, 2)
     & / -0.1815341E+01, -0.1538698E+01, -0.2014208E+01 /
      data (am( 7,k,-2),k=0, 2)
     & / -0.2571932E+02,  0.2903941E+00, -0.2848206E+01 /

      data mexvec(-3) / 7 /
      data mlfvec(-3) / 2 /
      data ut1vec(-3) /  0.4518239E+01 /
      data ut2vec(-3) / -0.2690590E+01 /
      data alfvec(-3) /  0.6124079E+00 /
      data qmavec(-3) /  0.0000000E+00 /
      data (am( 0,k,-3),k=0, 2)
     & / -0.2734458E+01, -0.7245673E+00, -0.6351374E+00 /
      data (am( 1,k,-3),k=0, 2)
     & /  0.2927174E+01,  0.4822709E+00, -0.1088787E-01 /
      data (am( 2,k,-3),k=0, 2)
     & / -0.1771017E+02, -0.1416635E+01,  0.8467622E+01 /
      data (am( 3,k,-3),k=0, 2)
     & / -0.4972782E+02, -0.3348547E+01,  0.1767061E+02 /
      data (am( 4,k,-3),k=0, 2)
     & / -0.7102770E+01, -0.3205337E+01,  0.4101704E+00 /
      data (am( 5,k,-3),k=0, 2)
     & /  0.7169698E+02, -0.2205985E+01, -0.2463931E+02 /
      data (am( 6,k,-3),k=0, 2)
     & / -0.4090347E+02,  0.2103486E+01,  0.1416507E+02 /
      data (am( 7,k,-3),k=0, 2)
     & / -0.2952639E+02,  0.5376136E+01,  0.7825585E+01 /

      data mexvec(-4) / 7 /
      data mlfvec(-4) / 2 /
      data ut1vec(-4) /  0.2783230E+01 /
      data ut2vec(-4) / -0.1746328E+01 /
      data alfvec(-4) /  0.1115653E+01 /
      data qmavec(-4) /  0.1300000E+01 /
      data (am( 0,k,-4),k=0, 2)
     & / -0.1743872E+01, -0.1128921E+01, -0.2841969E+00 /
      data (am( 1,k,-4),k=0, 2)
     & /  0.3345755E+01,  0.3187765E+00,  0.1378124E+00 /
      data (am( 2,k,-4),k=0, 2)
     & / -0.2037615E+02,  0.4121687E+01,  0.2236520E+00 /
      data (am( 3,k,-4),k=0, 2)
     & / -0.4703104E+02,  0.5353087E+01, -0.1455347E+01 /
      data (am( 4,k,-4),k=0, 2)
     & / -0.1060230E+02, -0.1551122E+01, -0.1078863E+01 /
      data (am( 5,k,-4),k=0, 2)
     & /  0.5088892E+02, -0.8197304E+01,  0.8083451E+01 /
      data (am( 6,k,-4),k=0, 2)
     & / -0.2819070E+02,  0.4554086E+01, -0.5890995E+01 /
      data (am( 7,k,-4),k=0, 2)
     & / -0.1098238E+02,  0.2590096E+01, -0.8062879E+01 /

      data mexvec(-5) / 6 /
      data mlfvec(-5) / 2 /
      data ut1vec(-5) /  0.1619654E+02 /
      data ut2vec(-5) / -0.3367346E+01 /
      data alfvec(-5) /  0.5109891E-02 /
      data qmavec(-5) /  0.4500000E+01 /
      data (am( 0,k,-5),k=0, 2)
     & / -0.6800138E+01,  0.2493627E+01, -0.1075724E+01 /
      data (am( 1,k,-5),k=0, 2)
     & /  0.3036555E+01,  0.3324733E+00,  0.2008298E+00 /
      data (am( 2,k,-5),k=0, 2)
     & / -0.5203879E+01, -0.8493476E+01, -0.4523208E+01 /
      data (am( 3,k,-5),k=0, 2)
     & / -0.1524239E+01, -0.3411912E+01, -0.1771867E+02 /
      data (am( 4,k,-5),k=0, 2)
     & / -0.1099444E+02,  0.1320930E+01, -0.2353831E+01 /
      data (am( 5,k,-5),k=0, 2)
     & /  0.1699299E+02, -0.3565802E+02,  0.3566872E+02 /
      data (am( 6,k,-5),k=0, 2)
     & / -0.1465793E+02,  0.2703365E+02, -0.2176372E+02 /



      if(q .le. qmavec(ifl)) then
         faux5MI = 0.d0
         return
      endif

      if(x .ge. 1.d0) then
         faux5MI = 0.d0
         return
      endif

      tmp = log(q/alfvec(ifl))
      if(tmp .le. 0.d0) then
         faux5MI = 0.d0
         return
      endif

      sb = log(tmp)
      sb1 = sb - 1.2d0
      sb2 = sb1*sb1

      do i = 0, nex
         af(i) = 0.d0
         sbx = 1.d0
         do k = 0, mlfvec(ifl)
            af(i) = af(i) + sbx*am(i,k,ifl)
            sbx = sb1*sbx
         enddo
      enddo

      y = -log(x)
      u = log(x/0.00001d0)

      part1 = af(1)*y**(1.d0+0.01d0*af(4))*(1.d0+ af(8)*u)
      part2 = af(0)*(1.d0 - x) + af(3)*x 
      part3 = x*(1.d0-x)*(af(5)+af(6)*(1.d0-x)+af(7)*x*(1.d0-x))
      part4 = ut1vec(ifl)*log(1.d0-x) + 
     &	      AF(2)*log(1.d0+exp(ut2vec(ifl))-x)

      faux5MI = exp(log(x) + part1 + part2 + part3 + part4)

c include threshold factor...
      faux5MI = faux5MI * (1.d0 - qmavec(ifl)/q)

      return
      end
C ==============================================================================
C          CTEQ5 Parton Distribution Functions in Parametrized Form
C                             
C               Preliminary version: May 25, 1999
C
C   Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C         CTEQ5 PPARTON DISTRIBUTIONS"
C   hep-ph/9903282
C
C   These parametrizations were obtained by Jon Pumplin.
C
C   Since this is an initial distribution of the parametrized version, it
C   is useful for the authors to know any possible problems.  Please send
C   any questions or comments to:   Tung@pa.msu.edu or Pumplin@pa.msu.edu
C
C   The calling sequences for the parametrization are the same as those 
C   that were used for the CTEQ2, CTEQ3 and CTEQ4 parton distributions.
C
C
C     Name convention for CTEQ distributions:  CTEQnSx  where
C        n : version number                      (currently n = 5)
C        S : factorization scheme label: = [M L D] for [MS-bar LO DIS] 
c               resp.
C        x : special characteristics, if any
C
C  7 sets of CTEQ5 PDF's have been distributed in tablular form. 
C ---------------------------------------------------------------------------
C  Iset   PDF        Description       Alpha_s(Mz)  Lam4  Lam5
C ---------------------------------------------------------------------------
C   1    CTEQ5M   Standard MSbar scheme   0.118     326   226
C   2    CTEQ5D   Standard DIS scheme     0.118     326   226
C   3    CTEQ5L   Leading Order           0.127     192   146
C   4    CTEQ5HJ  Large-x gluon enhanced  0.118     326   226
C   5    CTEQ5HQ  Heavy Quark             0.118     326   226
C   6    CTEQ5F3  Nf=3 FixedFlavorNumber  0.106     (Lam3=395)
C   7    CTEQ5F4  Nf=4 FixedFlavorNumber  0.112     309   XXX
C ---------------------------------------------------------------------------
C                    ****************************** 
C!!!!     For now, only the most widely used sets CTEQ5M and CTEQ5L    !!!!
C      !!!!     are available here in parametrized form.          !!!!
C                    ******************************
C
C   The following user-callable routines are provided:
C 
C     FUNCTION Ctq5Pd (Iset, Iprtn, X, Q, Irt) 
C         returns the PROBABILITY density for a GIVEN flavor;
C
C     FUNCTION Ctq5df (Iset, Iprtn, X, Q, Irt)
C         returns the MOMENTUM density of a GIVEN valence or sea distribution.
C
C     SUBROUTINE Ctq5Pds(Iset, Pdf, X, Q, Irt)
C         returns an array of MOMENTUM densities for ALL flavors;
C
C     FUNCTION WLAMD5 (ISET, IORDER, NEFF) 
C         returns the QCD lambda parameter appropriate to the distributions.
C
C   The arguments of these routines are as follows: 
C
C ******  Iset is the set number:  1 for CTEQ5M or 3 or CTEQ5L  ******
C
C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)
C  *** WARNING: We use the parton label 2 as D-quark and 1 as U-quark, 
C               which might be different from your labels.
C
C   X, Q are the usual x, Q; 
C
C   Irt is an error code: 0 if there was no error; 1 or more if (x,q) was 
C   outside the range of validity of the parametrization.
C       
C
C  Since the QCD Lambda value for the various sets are needed more often than
C  the other parameters in most applications, a special function
C     WLAMD5 (ISET, IORDER, NEFF)                    is provided
C  which returns the lambda value for Neff = 4,5,6 effective flavors as well as
C  the order these values pertain to.
C
C
C  Range of validity:
C  
C     The range of (x, Q) covered by this parametrization of the QCD evolved
C     parton distributions is 1E-6 < x < 1 ; 1.1 GeV < Q < 10 TeV.  Of course,
C     the PDF's are constrained by data only in a subset of that region; and 
C     the assumed DGLAP evolution is unlikely to be valid for all of it either.
C
C     The range of (x, Q) used in the CTEQ5 round of global analysis is 
C     approximately  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for 
C     fixed target experiments and 0.0001 < x < 0.1 from HERA data.  
C
C
C   DOUBLE PRECISION is used throughout in these routines, but conversion to 
C   SINGLE PRECISION is possible by removing the Implicit Double Precision statements. 
C
C **************************************************************************

C ********************************************************
      FUNCTION CTQ5PD(ISET, IPARTON, X, Q, IRT)
C ********************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

c if called at a point (x,q) that is outside the region that was 
c actually parametrized, return a value of 0, and set the error code IRT=1.  
c The user can remove the following IF statement to receive instead an 
c extrapolated value, which may be wildly unphysical.
      if((x .lt. 1.e-6). or. (x .gt. 1.) 
     &	 .or. (q .lt. .99) .or. (q .gt. 10000.)) then
         ctq5pd = 0.d0
         irt = 1
         return
      endif

      irt = 0
      if(iset .eq. 3) then
         ctq5pd = ctq5L(iparton,x,q)
      elseif(iset .eq. 1) then
         ctq5pd = ctq5Mi(iparton,x,q)
      else
         print *,'iset=',iset,' has not been parametrized.' 
	   print '(/A)', 'Use the interpolation-table version instead.'
         stop
      endif

      return
      end

C ********************************************************
      FUNCTION CTQ5DF(ISET, IFL, X, Q, IRT)
C ********************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CTQ5DF = X * CTQ5PD(ISET, IPARTON, X, Q, IRT)
        
      RETURN
      END

C ********************************************************
      SUBROUTINE CTQ5PDS(ISET, PDF, X, Q, IRT)
C ********************************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION PDF (-6:6)

      IRT = 0

      DO IFL= -6,2
         PDF(IFL) = CTQ5PD(ISET,IFL,X,Q,IRT1)
         IRT = IRT + IRT1

         IF (IFL .LE. -3) THEN
            PDF(-IFL) = PDF(IFL)
         ENDIF

      ENDDO

      RETURN
      END

      Function Ctq5Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
	Print *, 'X out of range in Ctq5Pdf: ', X
	Stop
      Endif
      If (Q .lt. Alambda) Then
	Print *, 'Q out of range in Ctq5Pdf: ', Q
	Stop
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
	     Warn = .false.
	     Print *, 'Warning: Iparton out of range in Ctq5Pdf: '
     >              , Iparton
         Endif
         Ctq5Pdf = 0D0
         Return
      Endif

      Ctq5Pdf = PartonX (Iparton, X, Q)
      if(Ctq5Pdf.lt.0.D0)  Ctq5Pdf = 0.D0

      Return

C                             ********************
      End

C============================================================================
C                CTEQ Parton Distribution Functions: version 6 & 6.5
C                             April 10, 2002, v6.01
C                             February 23, 2003, v6.1
C                             August 6, 2003, v6.11
C                             December 12, 2004, v6.12
C                             December 4, 2006, v6.5 (CTEQ6.5M series added)
C                             March 23, 2007, v6.51 (CTEQ6.5S/C series added)
C                             April 24, 2007, v6.52 (minor improvement)
C
C   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
C       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       JHEP 0207:012(2002), hep-ph/0201195
C
C   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
C       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
C       JHEP 0310:046(2003), hep-ph/0303013
C
C   Ref[3]: "Neutrino dimuon Production and Strangeness Asymmetry of the Nucleon"
C       By: F. Olness, J. Pumplin, S. Stump, J. Huston, P. Nadolsky, H.L. Lai, S. Kretzer, J.F. Owens, W.K. Tung
C       Eur. Phys. J. C40:145(2005), hep-ph/0312323
C
C   Ref[4]: "CTEQ6 Parton Distributions with Heavy Quark Mass Effects"
C       By: S. Kretzer, H.L. Lai, F. Olness, W.K. Tung
C       Phys. Rev. D69:114005(2004), hep-ph/0307022
C
C   Ref[5]: "Heavy Quark Mass Effects in Deep Inelastic Scattering and Global QCD Analysis"
C       By : W.K. Tung, H.L. Lai, A. Belyaev, J. Pumplin, D. Stump, C.-P. Yuan
C       JHEP 0702:053(2007), hep-ph/0611254
C
C   Ref[6]: "The Strange Parton Distribution of Nucleon: Global Analysis and Applications"
C       By : H.L. Lai, P. Nadolsky, J. Pumplin, D. Stump, W.K. Tung, C.-P. Yuan
C       hep-ph/0702268, to appear in JHEP
C
C   Ref[7]: "The Charm Content of the Nucleon"
C       By : J. Pumplin, H.L. Lai, W.K. Tung
C       hep-ph/0701220, to appear in Phys. Rev. D
C
C   This package contains
C   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
C   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
C   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].
C   (4) 5 special sets for strangeness study from Ref[3].
C   (5) 1 special set for heavy quark study from Ref[4].
C   (6) CTEQ6.5M and its 40 up/down eigenvector sets from Ref[5].
C   (7) 8 sets of PDFs resulting from the strangeness study, Ref[6].
C   (8) 7 sets of PDFs resulting from the charm study, Ref[7].
C

C  Details about calling convention are:
C ---------------------------------------------------------------------------
C  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
C ===========================================================================
C Standard, "best-fit", sets:          Ref.[1]
C --------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
C   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
C 200    CTEQ6.1M Ref.[2]: updated CTEQ6M (see below, under "uncertainty" section)
C --------------------------
C Special sets for strangeness study:  Ref.[3]
C --------------------------
C  11    CTEQ6A   Class A                 0.118     326   226    cteq6sa.pds
C  12    CTEQ6B   Class B                 0.118     326   226    cteq6sb.pds
C  13    CTEQ6C   Class C                 0.118     326   226    cteq6sc.pds
C  14    CTEQ6B+  Large [S-]              0.118     326   226    cteq6sb+.pds
C  15    CTEQ6B-  Negative [S-]           0.118     326   226    cteq6sb-.pds
C --------------------------
C Special set for Heavy Quark study:   Ref.[4]
C --------------------------
C  21    CTEQ6HQ                          0.118     326   226    cteq6hq.pds
C --------------------------
C Released sets for strangeness study:  Ref.[6]
C -------------------------- s=sbr
C  30    CTEQ6.5S0   Best-fit             0.118     326   226    ctq65.s+0.pds
C  31    CTEQ6.5S1   Low s+               0.118     326   226    ctq65.s+1.pds
C  32    CTEQ6.5S2   High s+              0.118     326   226    ctq65.s+2.pds
C  33    CTEQ6.5S3   Alt Low s+           0.118     326   226    ctq65.s+3.pds
C  34    CTEQ6.5S4   Alt High s+          0.118     326   226    ctq65.s+4.pds
C -------------------------- s!=sbr
C          strangeness asymmetry <x>_s-
C  35    CTEQ6.5S-0  Best-fit    0.0014    0.118     326   226    ctq65.s-0.pds
C  36    CTEQ6.5S-1  Low        -0.0010    0.118     326   226    ctq65.s-1.pds
C  37    CTEQ6.5S-2  High        0.0050    0.118     326   226    ctq65.s-2.pds
C --------------------------
C Released sets for charm study:  Ref.[7]
C --------------------------
C  40    CTEQ6.5C0   no intrinsic charm   0.118     326   226    ctq65.c0.pds
C  41    CTEQ6.5C1   BHPS model for IC    0.118     326   226    ctq65.c1.pds
C  42    CTEQ6.5C2   BHPS model for IC    0.118     326   226    ctq65.c2.pds
C  43    CTEQ6.5C3   Meson cloud model    0.118     326   226    ctq65.c3.pds
C  44    CTEQ6.5C4   Meson cloud model    0.118     326   226    ctq65.c4.pds
C  45    CTEQ6.5C5   Sea-like model       0.118     326   226    ctq65.c5.pds
C  46    CTEQ6.5C6   Sea-like model       0.118     326   226    ctq65.c6.pds
C     Momentum Fraction carried by c,cbar at Q0=1.3 GeV:
C    Iset:charm  ,cbar     | Iset:charm  ,cbar     | Iset:charm  ,cbar
C    41: 0.002857,0.002857 | 43: 0.003755,0.004817 | 45: 0.005714,0.005714
C    42: 0.010000,0.010000 | 44: 0.007259,0.009312 | 46: 0.012285,0.012285
C ============================================================================
C For uncertainty calculations using eigenvectors of the Hessian:
C ---------------------------------------------------------------
C     central + 40 up/down sets along 20 eigenvector directions
C                             -----------------------------
C                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
C                             -----------------------
C  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
C             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
C                              -----------------------
C  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
C             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects, Ref[5]:  central fit: CTEQ6.5M (=CTEQ65.00)
C                              -----------------------
C  3xx  CTEQ65.xx  +/- sets               0.118     326   226    ctq65.xx.pds
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 300      is CTEQ65.00 (=CTEQ6.5M),
C             301/302 are CTEQ65.01/02, +/- sets of 1st eigenvector, ... etc.
C ===========================================================================
C   ** ALL fits are obtained by using the same coupling strength
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in the references.
C
C   The table grids are generated for 10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV).
C   (For CTEQ6.5S/C series: 10^-7 < x < 1 and 1.3 < Q < 10^5 (GeV))
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   Pumplin@pa.msu.edu or Tung@pa.msu.edu.
C
C===========================================================================

      Function Ctq6Pdf (Iparton, X, Q)
      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0D0 .or. X .gt. 1D0) Then
!        Print *, 'X out of range in Ctq6Pdf: ', X
!        Stop
	Ctq6Pdf=0d0
      Endif
      If (Q .lt. Alambda) Then
!        Print *, 'Q out of range in Ctq6Pdf: ', Q
!        Stop
	Ctq6Pdf=0d0
      Endif
      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
             Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if(Ctq6Pdf.lt.0.D0) then 
         Ctq6Pdf = 0.D0
      Endif

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax0=7)
      Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
      Logical fmtpds
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.','cteq6s'
     >  ,'ctq65.' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
      Data Isetmin2,Isetmax2 /200,240/
      Data Isetmin3,Isetmax3 /300,340/
      Data IsetminS,IsetmaxS /11,15/
      Data IsetmnSp07,IsetmxSp07 /30,34/
      Data IsetmnSm07,IsetmxSm07 /35,37/
      Data IsetmnC07,IsetmxC07 /40,46/
      Data IsetHQ /21/
      Common /Setchange/ Isetch
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         fmtpds=.false.

         If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
C                                                     Iset = 1,2,3 for 6m, 6d, 6l
            Tablefile=Flnm(Iset)//'.tbl'
         Elseif (Iset.eq.4) Then
C                                                               4  (2nd LO fit)
            Tablefile=Flnm(Iset)//'1.tbl'
         Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
C                                                               101 - 140    
            write(nn,'(I3)') Iset
            Tablefile=Flnm(1)//nn//'.tbl'
         Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
C                                                               200 - 240   
            write(nn,'(I3)') Iset
            Tablefile=Flnm(5)//nn(2:3)//'.tbl'
         Elseif (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) Then
C                                                               11 - 15   
            fmtpds=.true.
            If(Iset.eq.11) then
               Tablefile=Flnm(6)//'a.pds'
            Elseif(Iset.eq.12) then
               Tablefile=Flnm(6)//'b.pds'
            Elseif(Iset.eq.13) then
               Tablefile=Flnm(6)//'c.pds'
            Elseif(Iset.eq.14) then
               Tablefile=Flnm(6)//'b+.pds'
            Elseif(Iset.eq.15) then
               Tablefile=Flnm(6)//'b-.pds'
            Endif
         Elseif (Iset.eq.IsetHQ) Then
C                                                               21   
            fmtpds=.true.
            TableFile='cteq6hq.pds'
         Elseif (Iset.ge.IsetmnSp07 .and. Iset.le.IsetmxSp07) Then
C                                                    (Cteq6.5S)  30 - 34   
            fmtpds=.true.
            write(nn,'(I2)') Iset
            Tablefile=Flnm(7)//'s+'//nn(2:2)//'.pds'
         Elseif (Iset.ge.IsetmnSm07 .and. Iset.le.IsetmxSm07) Then
C                                                    (Cteq6.5S)  35 - 37   
            fmtpds=.true.
            Is = Iset - 5
            write(nn,'(I2)') Is
            Tablefile=Flnm(7)//'s-'//nn(2:2)//'.pds'
         Elseif (Iset.ge.IsetmnC07 .and. Iset.le.IsetmxC07) Then
C                                                    (Cteq6.5C)  40 - 46
            fmtpds=.true.
            write(nn,'(I2)') Iset
            Tablefile=Flnm(7)//'c'//nn(2:2)//'.pds'
         Elseif (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) Then
C                                                    (Cteq6.5)  300 - 340   
            fmtpds=.true.
            write(nn,'(I3)') Iset
            Tablefile=Flnm(7)//nn(2:3)//'.pds'
         Else
            Print *, 'Invalid Iset number in SetCtq6 :', Iset
            Stop
         Endif
         IU= NextUn()
         Open(IU, File=Tablefile, Status='OLD', Err=100)
 21      Call Readpds (IU,fmtpds)
         Close (IU)
         Isetold=Iset
         Isetch=1
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq6!!'
      Stop
C                             ********************
      End

      Subroutine Readpds (Nu,fmtpds)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      Logical fmtpds
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
      If(fmtpds) then
C                                               This is the .pds (WKT) format
        Read  (Nu, *) N0, N0, N0, NfMx, MxVal, N0
        If(MxVal.gt.MaxVal) MxVal=3 !old .pds format (read in KF, not MxVal)
 
        Read  (Nu, '(A)') Line
        Read  (Nu, *) NX,  NT, N0, NG, N0

        Read  (Nu, '(A)') (Line,I=1,NG+2)
        Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
        XV(0)=0D0
      Else
C                                               This is the old .tbl (HLL) format
         MxVal=2
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log(Log (TV(Iq) /Al))
 11      Continue
      Endif

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End
C

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
      Implicit Double Precision (A-H,O-Z)

      PARAMETER (MXX = 105, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
 
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > /Setchange/ Isetch

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

      If((XX.eq.X).and.(QQ.eq.Q)) goto 99
c store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

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
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
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
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

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


c get the pdf function values at the lattice points...

 99   If (Iprtn .Gt. MxVal) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4F (XVpow(0), Fij(1), ss, Fx)

         If (x .GT. 0D0)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

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

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End

      SUBROUTINE POLINT4F (XA,YA,X,Y)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan. 
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN
      
      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
      END

c comment NextUn out, since it exists in the utility pack.
c      Function NextUn()
c                                 Returns an unallocated FORTRAN i/o unit.
c      Logical EX
c
c      Do 10 N = 10, 300
c         INQUIRE (UNIT=N, OPENED=EX)
c         If (.NOT. EX) then
c            NextUn = N
c            Return
c         Endif
c 10   Continue
c      Stop ' There is no available I/O unit. '
c               *************************
c      End

      FUNCTION DFLM (LPRTN, XX, QQ, NST)
C                                                   -=-=- dflm

C===========================================================================
C GroupName: MiscPd
C Description: Other parton distributions: dflm, morfin-tung, GRV ...
C ListOfFiles: dflm pdxmt pdzxmt pdzmt pdzzmt pdgrv grv94lo grv94ho grv94di fv fw fws
C===========================================================================
C #Header: /Net/d2a/wkt/1hep/2pdf/prz/RCS/MiscPd.f,v 6.2 98/08/14 09:35:47 wkt Exp $
C #Log:	MiscPd.f,v $
c Revision 6.2  98/08/14  09:35:47  wkt
c GRV comment lines standardized
c 
c Revision 6.1  97/11/15  18:07:05  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c 
C========================================================================
C      FUNCTION DFLM (LPRTN, XX, QQ, NST)

      DOUBLE PRECISION DFLM, XX, QQ
C                                     Note that the DFLM routines are in single
C                                   precision  and they call a single-precision
C                                     interpolating routine from KERNLIB
      TEM = 0.
      X  = XX

C     IF (X .GT. 0.95) THEN
C       DFLM = 0.
C       RETURN
C     EndIf

C     Q2 = QQ ** 2

C     IF     (LPRTN .EQ. 0) THEN
C       IF (NST .EQ. 1) CALL FXAVER(X, Q2, 'GLUON', TEM)
C       IF (NST .EQ. 2) CALL FXNLLA(X, Q2, 'GLUON', TEM)
C     ElseIF (LPRTN .EQ. 1) THEN
C       IF (NST .EQ. 1) THEN
C         CALL FXAVER(X, Q2, 'UBAR ', TEM1)
C         CALL FXAVER(X, Q2, 'UPVAL', TEM)
C       ElseIF (NST .EQ. 2) THEN
C         CALL FXNLLA(X, Q2, 'UBAR ', TEM1)
C         CALL FXNLLA(X, Q2, 'UPVAL', TEM)
C       EndIf
C       TEM = TEM + TEM1
C     ElseIF (LPRTN .EQ. 2) THEN
C       IF (NST .EQ. 1) THEN
C         CALL FXAVER(X, Q2, 'UBAR ', TEM1)
C         CALL FXAVER(X, Q2, 'DOVAL', TEM)
C       ElseIF (NST .EQ. 2) THEN
C         CALL FXNLLA(X, Q2, 'UBAR ', TEM1)
C         CALL FXNLLA(X, Q2, 'DOVAL', TEM)
C       EndIf
C       TEM = TEM + TEM1
C     ElseIF (LPRTN .EQ. -1 .OR. LPRTN .EQ. -2) THEN
C       IF (NST .EQ. 1) CALL FXAVER(X, Q2, 'UBAR ', TEM)
C       IF (NST .EQ. 2) CALL FXNLLA(X, Q2, 'UBAR ', TEM)
C     ElseIF (ABS(LPRTN) .EQ. 3) THEN
C       IF (NST .EQ. 1) CALL FXAVER(X, Q2, 'SBAR ', TEM)
C       IF (NST .EQ. 2) CALL FXNLLA(X, Q2, 'SBAR ', TEM)
C     ElseIF (ABS(LPRTN) .EQ. 4) THEN
C       IF (NST .EQ. 1) CALL FXAVER(X, Q2, 'CBAR ', TEM)
C       IF (NST .EQ. 2) CALL FXNLLA(X, Q2, 'CBAR ', TEM)
C     ElseIF (ABS(LPRTN) .EQ. 5) THEN
C       IF (NST .EQ. 1) CALL FXAVER(X, Q2, 'BBAR ', TEM)
C       IF (NST .EQ. 2) CALL FXNLLA(X, Q2, 'BBAR ', TEM)
C     ElseIF (ABS(LPRTN) .EQ. 6) THEN
C       IF (NST .EQ. 1) CALL FXAVER(X, Q2, 'TBAR ', TEM)
C       IF (NST .EQ. 2) CALL FXNLLA(X, Q2, 'TBAR ', TEM)
C     Else
C       PRINT *, 'LPRTN = ', LPRTN, ' out of range in DFLM!'
C     EndIf

      DFLM = TEM / X
      RETURN
C                        ****************************
      END
       FUNCTION FV (X, N, AK, BK, A, B, C, D)
C                                                   -=-=- fv
       IMPLICIT DOUBLE PRECISION (A - Z)
       DX = DSQRT (X)
       FV = N * X**AK * (1.+ A*X**BK + X * (B + C*DX)) * (1.- X)**D
       RETURN
       END
C *
       FUNCTION FW (X, S, AL, BE, AK, BK, A, B, C, D, E, ES)
C                                                   -=-=- fw
       IMPLICIT DOUBLE PRECISION (A - Z)
       LX = DLOG (1./X)
       FW = (X**AK * (A + X * (B + X*C)) * LX**BK + S**AL
     1      * DEXP (-E + DSQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
C *
       FUNCTION FWS (X, S, AL, BE, AK, AG, B, D, E, ES)
C                                                   -=-=- fws
       IMPLICIT DOUBLE PRECISION (A - Z)
       DX = DSQRT (X)
       LX = DLOG (1./X)
       FWS = S**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1       * DEXP (-E + DSQRT (ES * S**BE * LX))
       RETURN
       END


       SUBROUTINE GRV94DI (X, Q2, UV, DV, DEL, UDB, SB, GL)
C                                                   -=-=- grv94di
       IMPLICIT DOUBLE PRECISION (A - Z)
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       S  = DLOG (DLOG(Q2/LAM2) / DLOG(MU2/LAM2))
       DS = DSQRT (S)
       S2 = S * S
       S3 = S2 * S
C *...UV :
       NU  =  2.484 + 0.116 * S + 0.093 * S2
       AKU =  0.563 - 0.025 * S
       BKU =  0.054 + 0.154 * S
       AU  = -0.326 - 0.058 * S - 0.135 * S2
       BU  = -3.322 + 8.259 * S - 3.119 * S2 + 0.291 * S3
       CU  =  11.52 - 12.99 * S + 3.161 * S2
       DU  =  2.808 + 1.400 * S - 0.557 * S2 + 0.119 * S3
       UV  = FV (X, NU, AKU, BKU, AU, BU, CU, DU)
C *...DV :
       ND  =  0.156 - 0.017 * S
       AKD =  0.299 - 0.022 * S
       BKD =  0.259 - 0.015 * S
       AD  =  3.445 + 1.278 * S + 0.326 * S2
       BD  = -6.934 + 37.45 * S - 18.95 * S2 + 1.463 * S3
       CD  =  55.45 - 69.92 * S + 20.78 * S2
       DD  =  3.577 + 1.441 * S - 0.683 * S2 + 0.179 * S3
       DV  = FV (X, ND, AKD, BKD, AD, BD, CD, DD)
C *...DEL :
       NE  =  0.099 + 0.019 * S + 0.002 * S2 
       AKE =  0.419 - 0.013 * S
       BKE =  1.064 - 0.038 * S
       AE  = -44.00 + 98.70 * S - 14.79 * S2
       BE  =  28.59 - 40.94 * S - 13.66 * S2 + 2.523 * S3
       CE  =  84.57 - 108.8 * S + 31.52 * S2
       DE  =  7.469 + 2.480 * S - 0.866 * S2
       DEL = FV (X, NE, AKE, BKE, AE, BE, CE, DE)
C *...UDB :
       ALX =  1.215
       BEX =  0.466
       AKX =  0.326 + 0.150 * S
       BKX =  0.956 + 0.405 * S
       AGX =  0.272
       BGX =  3.794 - 2.359 * DS
       CX  =  2.014
       DX  =  7.941 + 0.534 * DS - 0.940 * S + 0.410 * S2
       EX  =  3.049 + 1.597 * S
       ESX =  4.396 - 4.594 * DS + 3.268 * S
       UDB = FW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
C *...SB :
       ALS =  0.175
       BES =  0.344
       AKS =  1.415 - 0.641 * DS
       AS  =  0.580 - 9.763 * DS + 6.795 * S  - 0.558 * S2
       BS  =  5.617 + 5.709 * DS - 3.972 * S
       DST =  13.78 - 9.581 * S  + 5.370 * S2 - 0.996 * S3
       EST =  4.546 + 0.372 * S2
       ESS =  5.053 - 1.070 * S  + 0.805 * S2
       SB  = FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
C *...GL :
       ALG =  1.258
       BEG =  1.846
       AKG =  2.423
       BKG =  2.427 + 1.311 * S  - 0.153 * S2
       AG  =  25.09 - 7.935 * S
       BG  = -14.84 - 124.3 * DS + 72.18 * S
       CG  =  590.3 - 173.8 * S
       DG  =  5.196 + 1.857 * S
       EG  = -1.648 + 3.988 * S  - 0.432 * S2
       ESG =  3.232 - 0.542 * S
       GL  = FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
C *
C *...FUNCTIONAL FORMS OF THE PARAMETRIZATIONS :
C *
       SUBROUTINE GRV94HO (X, Q2, UV, DV, DEL, UDB, SB, GL)
C                                                   -=-=- grv94ho
       IMPLICIT DOUBLE PRECISION (A - Z)
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       S  = DLOG (DLOG(Q2/LAM2) / DLOG(MU2/LAM2))
       DS = DSQRT (S)
       S2 = S * S
       S3 = S2 * S
C *...UV :
       NU  =  1.304 + 0.863 * S
       AKU =  0.558 - 0.020 * S
       BKU =          0.183 * S
       AU  = -0.113 + 0.283 * S - 0.321 * S2
       BU  =  6.843 - 5.089 * S + 2.647 * S2 - 0.527 * S3
       CU  =  7.771 - 10.09 * S + 2.630 * S2
       DU  =  3.315 + 1.145 * S - 0.583 * S2 + 0.154 * S3
       UV  = FV (X, NU, AKU, BKU, AU, BU, CU, DU)
C *...DV :
       ND  =  0.102 - 0.017 * S + 0.005 * S2
       AKD =  0.270 - 0.019 * S
       BKD =  0.260
       AD  =  2.393 + 6.228 * S - 0.881 * S2
       BD  =  46.06 + 4.673 * S - 14.98 * S2 + 1.331 * S3
       CD  =  17.83 - 53.47 * S + 21.24 * S2
       DD  =  4.081 + 0.976 * S - 0.485 * S2 + 0.152 * S3
       DV  = FV (X, ND, AKD, BKD, AD, BD, CD, DD)
C *...DEL :
       NE  =  0.070 + 0.042 * S - 0.011 * S2 + 0.004 * S3
       AKE =  0.409 - 0.007 * S
       BKE =  0.782 + 0.082 * S
       AE  = -29.65 + 26.49 * S + 5.429 * S2
       BE  =  90.20 - 74.97 * S + 4.526 * S2
       CE  =  0.0
       DE  =  8.122 + 2.120 * S - 1.088 * S2 + 0.231 * S3
       DEL = FV (X, NE, AKE, BKE, AE, BE, CE, DE)
C *...UDB :
       ALX =  0.877
       BEX =  0.561
       AKX =  0.275
       BKX =  0.0
       AGX =  0.997
       BGX =  3.210 - 1.866 * S
       CX  =  7.300
       DX  =  9.010 + 0.896 * DS + 0.222 * S2
       EX  =  3.077 + 1.446 * S
       ESX =  3.173 - 2.445 * DS + 2.207 * S
       UDB = FW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
C *...SB :
       ALS =  0.756
       BES =  0.216
       AKS =  1.690 + 0.650 * DS - 0.922 * S
       AS  = -4.329 + 1.131 * S
       BS  =  9.568 - 1.744 * S
       DST =  9.377 + 1.088 * DS - 1.320 * S + 0.130 * S2
       EST =  3.031 + 1.639 * S
       ESS =  5.837 + 0.815 * S
       SB  = FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
C *...GL :
       ALG =  1.014
       BEG =  1.738
       AKG =  1.724 + 0.157 * S
       BKG =  0.800 + 1.016 * S
       AG  =  7.517 - 2.547 * S
       BG  =  34.09 - 52.21 * DS + 17.47 * S
       CG  =  4.039 + 1.491 * S
       DG  =  3.404 + 0.830 * S
       EG  = -1.112 + 3.438 * S  - 0.302 * S2
       ESG =  3.256 - 0.436 * S
       GL  = FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
C *
C *...NLO PARAMETRIZATION (DIS) :
C *
       SUBROUTINE GRV94LO (X, Q2, UV, DV, DEL, UDB, SB, GL)
C                                                   -=-=- grv94lo
       IMPLICIT DOUBLE PRECISION (A - Z)
       MU2  = 0.23
       LAM2 = 0.2322 * 0.2322
       S  = DLOG (DLOG(Q2/LAM2) / DLOG(MU2/LAM2))
       DS = DSQRT (S)
       S2 = S * S
       S3 = S2 * S
C *...UV :
       NU  =  2.284 + 0.802 * S + 0.055 * S2
       AKU =  0.590 - 0.024 * S
       BKU =  0.131 + 0.063 * S
       AU  = -0.449 - 0.138 * S - 0.076 * S2
       BU  =  0.213 + 2.669 * S - 0.728 * S2
       CU  =  8.854 - 9.135 * S + 1.979 * S2
       DU  =  2.997 + 0.753 * S - 0.076 * S2
       UV  = FV (X, NU, AKU, BKU, AU, BU, CU, DU)
C *...DV :
       ND  =  0.371 + 0.083 * S + 0.039 * S2
       AKD =  0.376
       BKD =  0.486 + 0.062 * S
       AD  = -0.509 + 3.310 * S - 1.248 * S2
       BD  =  12.41 - 10.52 * S + 2.267 * S2
       CD  =  6.373 - 6.208 * S + 1.418 * S2
       DD  =  3.691 + 0.799 * S - 0.071 * S2
       DV  = FV (X, ND, AKD, BKD, AD, BD, CD, DD)
C *...DEL :
       NE  =  0.082 + 0.014 * S + 0.008 * S2
       AKE =  0.409 - 0.005 * S
       BKE =  0.799 + 0.071 * S
       AE  = -38.07 + 36.13 * S - 0.656 * S2
       BE  =  90.31 - 74.15 * S + 7.645 * S2
       CE  =  0.0
       DE  =  7.486 + 1.217 * S - 0.159 * S2
       DEL = FV (X, NE, AKE, BKE, AE, BE, CE, DE)
C *...UDB :
       ALX =  1.451
       BEX =  0.271
       AKX =  0.410 - 0.232 * S
       BKX =  0.534 - 0.457 * S
       AGX =  0.890 - 0.140 * S
       BGX = -0.981
       CX  =  0.320 + 0.683 * S
       DX  =  4.752 + 1.164 * S + 0.286 * S2
       EX  =  4.119 + 1.713 * S
       ESX =  0.682 + 2.978 * S
       UDB = FW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
C *...SB :
       ALS =  0.914
       BES =  0.577
       AKS =  1.798 - 0.596 * S
       AS  = -5.548 + 3.669 * DS - 0.616 * S
       BS  =  18.92 - 16.73 * DS + 5.168 * S
       DST =  6.379 - 0.350 * S  + 0.142 * S2
       EST =  3.981 + 1.638 * S
       ESS =  6.402
       SB  = FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
C *...GL :
       ALG =  0.524
       BEG =  1.088
       AKG =  1.742 - 0.930 * S
       BKG =        - 0.399 * S2
       AG  =  7.486 - 2.185 * S
       BG  =  16.69 - 22.74 * S  + 5.779 * S2
       CG  = -25.59 + 29.71 * S  - 7.296 * S2
       DG  =  2.792 + 2.215 * S  + 0.422 * S2 - 0.104 * S3
       EG  =  0.807 + 2.005 * S 
       ESG =  3.841 + 0.316 * S
       GL  = FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       RETURN
       END
C *
C *...NLO PARAMETRIZATION (MS(BAR)) :
C *
      subroutine cjeppe1(nx,my,xx,yy,ff,cc)
      implicit real*8(a-h,o-z)
      parameter(nnx=49,mmy=37)
      dimension xx(nx),yy(my),ff(nnx,mmy),ff1(nnx,mmy),ff2(nnx,mmy),
     xff12(nnx,mmy),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),
     xcl(16),cc(nx,my,4,4),iwt(16,16)

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
     x		  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
     x		  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
     x		  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
     x		  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
     x		  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
     x		  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
     x		  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
     x		  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
     x		  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
     x		  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
     x		  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
     x		  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


      do 42 m=1,my
      dx=xx(2)-xx(1)
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
      do 41 n=2,nx-1
      ff1(n,m)=cpolderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),
     xff(n+1,m))
   41 continue
   42 continue

      do 44 n=1,nx
      dy=yy(2)-yy(1)
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy
      dy=yy(my)-yy(my-1)
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
      do 43 m=2,my-1
      ff2(n,m)=cpolderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),
     xff(n,m+1))
   43 continue
   44 continue

      do 46 m=1,my
      dx=xx(2)-xx(1)
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
      dx=xx(nx)-xx(nx-1)
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
      do 45 n=2,nx-1
      ff12(n,m)=cpolderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),
     xff2(n+1,m))
   45 continue
   46 continue

      do 53 n=1,nx-1
      do 52 m=1,my-1
      d1=xx(n+1)-xx(n)
      d2=yy(m+1)-yy(m)
      d1d2=d1*d2

      yy0(1)=ff(n,m)
      yy0(2)=ff(n+1,m)
      yy0(3)=ff(n+1,m+1)
      yy0(4)=ff(n,m+1)

      yy1(1)=ff1(n,m)
      yy1(2)=ff1(n+1,m)
      yy1(3)=ff1(n+1,m+1)
      yy1(4)=ff1(n,m+1)

      yy2(1)=ff2(n,m)
      yy2(2)=ff2(n+1,m)
      yy2(3)=ff2(n+1,m+1)
      yy2(4)=ff2(n,m+1)

      yy12(1)=ff12(n,m)
      yy12(2)=ff12(n+1,m)
      yy12(3)=ff12(n+1,m+1)
      yy12(4)=ff12(n,m+1)

      do 47 k=1,4
      z(k)=yy0(k)
      z(k+4)=yy1(k)*d1
      z(k+8)=yy2(k)*d2
      z(k+12)=yy12(k)*d1d2
   47 continue

      do 49 l=1,16
      xxd=0.
      do 48 k=1,16
      xxd=xxd+iwt(k,l)*z(k)
   48 continue
      cl(l)=xxd
   49 continue
      l=0
      do 51 k=1,4
      do 50 j=1,4
      l=l+1
      cc(n,m,k,j)=cl(l)
   50 continue
   51 continue
   52 continue
   53 continue
      return
      end

      subroutine cjeppe2(x,y,nx,my,xx,yy,cc,z)
      implicit real*8(a-h,o-z)
      dimension xx(nx),yy(my),cc(nx,my,4,4)      

      n=llocx(xx,nx,x)
      m=llocx(yy,my,y)

      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))

      z=0.
      do 1 l=4,1,-1
      z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u
     .       +cc(n,m,l,2))*u+cc(n,m,l,1)
    1 continue
      return
      end

      integer function llocx(xx,nx,x)
      implicit real*8(a-h,o-z)
      dimension xx(nx)
      if(x.le.xx(1)) then
      llocx=1
      return
      endif
      if(x.ge.xx(nx)) then 
      llocx=nx-1  
      return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) go to 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
      jl=jm
      else
      ju=jm
      endif
      go to 1
    2 llocx=jl
      return
      end


      real*8 function  cpolderiv(x1,x2,x3,y1,y2,y3)
      implicit real*8(a-h,o-z)
      cpolderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*
     .(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end
C                                                                 Sep.12, 99
       SUBROUTINE MRSEB
C                                                   -=-=- mrseb
     $     (X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,xGLU) 
C***************************************************************C      
C                                                               C
C                                                               C
C     NEW MRS VERSIONS: S0', D0', D-'  (NOVEMBER  1992)         C
C         MODE 1 IS THE 1990 KMRS(B0) SET;                      C
C         MODES 2-4 ARE NEW SETS FITTED TO THE RECENT NMC       C
C     AND CCFR PRELIMINARY STRUCTURE FUNCTION DATA.             C
C     THE THREE NEW SETS HAVE LAMBDA(MSbar,NF=4) = 215 MeV      C
C         MODES 5-7 ARE NEW SETS FITTED TO THE RECENT NMC       C
C         AND CCFR ***FINAL***  STRUCTURE FUNCTION DATA.        C
C         DIFFERENCES BETWEEN THESE AND THE APRIL 1992          C
C         VERSIONS ARE SMALL! THE THREE NEW SETS HAVE           C
C         LAMBDA(MSbar,NF=4) = 230 MeV                          C
C                                                               C
C     THE REFERENCE IS: A.D. Martin, R.G. Roberts and           C
C     W.J. Stirling, RAL  preprint RAL-92-078 (1992)            C
C                                                               C
C        MODE 1: KMRS(B0) (Lambda(4) = 190 MeV)                 C
C        MODE 2: MRS(S0) (updated B0, Lambda(4) = 215 MeV)      C
C        MODE 3: MRS(D0) (... but with ubar not= dbar)          C
C        MODE 4: MRS(D-) (updated B-, ubar not= dbar)           C
C        MODE 5: MRS(S0') (updated B0, Lambda(4) = 230 MeV)     C
C        MODE 6: MRS(D0') (... but with ubar not= dbar)         C
C        MODE 7: MRS(D-') (updated B-, ubar not= dbar)          C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C    THE FIRST NUMBER IN THE "1" GRID IS 0.01727                C
C    THE FIRST NUMBER IN THE "2" GRID IS 0.01356                C
C    THE FIRST NUMBER IN THE "3" GRID IS 0.00527                C
C    THE FIRST NUMBER IN THE "4" GRID IS 0.00474                C
C    THE FIRST NUMBER IN THE "5" GRID IS 0.01617                C
C    THE FIRST NUMBER IN THE "6" GRID IS 0.00820                C
C    THE FIRST NUMBER IN THE "7" GRID IS 0.00678                C
C                                                               C
C    NOTE THE EXTRA ARGUMENT IN THIS SUBROUTINE MRSEB,          C
C    TO ALLOW FOR THE POSSIBILITY OF A *** DIFFERENT ***        C
C    UBAR AND DBAR SEA!                                         C
C                                                               C
C                         -*-                                   C
C                                                               C
C    (NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION        C
C    IS RETURNED I.E. G(X) = GLU/X ETC. IF IN DOUBT, CHECK THE  C
C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV)          C
C                                                               C
C                         -*-                                   C
C                                                               C
C     (THE RANGE OF APPLICABILITY IS CURRENTLY:                 C
C     10**-5 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C
C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C
C     - PROBLEMS, COMMENTS ETC TO WJS@UK.AC.DUR.HEP             C
C                                                               C
C                                                               C
C***************************************************************C

      IMPLICIT REAL*8(A-H,O-Z)

      IF(MODE.EQ.1) then
           CALL STRC78(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      ElseIF(MODE.EQ.2) then
           CALL STRC79(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,xGLU)
      ElseIF(MODE.EQ.3) then
           CALL STRC80(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,xGLU)
      ElseIF(MODE.EQ.4) then
           CALL STRC81(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,xGLU)
      ElseIF(MODE.EQ.5) then
           CALL STRC82(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      ElseIF(MODE.EQ.6) then
           CALL STRC83(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      ElseIF(MODE.EQ.7) then
           CALL STRC84(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      ElseIF(MODE.EQ.8) then
           CALL STRC33(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      ElseIF(MODE.EQ.9) then
           CALL STRC30(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      ElseIF(MODE.EQ.10) then
           CALL STRC37(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
      Elseif(mode.ge.11 .and. mode.le.14) then
        Nst=mode-10
        q2=SCALE**2
        call mrsX(x,q2,upv,dnv,usea,dsea,str,chm,bot,xglu,Nst,96) 
      Elseif(mode.ge.21 .and. mode.le.25) then
        Nst=mode-20
        q2=SCALE**2
        call mrsX(x,q2,upv,dnv,usea,dsea,str,chm,bot,xglu,Nst,98) 
      Elseif(mode.ge.101 .and. mode.le.112) then
        Nst=mode-100
        q2=SCALE**2
        call mrsX(x,q2,upv,dnv,usea,dsea,str,chm,bot,xglu,Nst,99) 
      Elseif(mode.ge.201 .and. mode.le.204) then
        Nst=mode-200
        call mrst2001(x,scale,Nst,upv,dnv,usea,dsea,str,chm,bot,xglu)
      Elseif(mode.ge.211 .and. mode.le.214) then
        Nst=mode-210
        call mrst23(x,scale,Nst,upv,dnv,usea,dsea,str,chm,bot,xglu)
      Elseif(mode.eq.901) then
        Nst=mode-899
        call alekhin(x,scale,Nst,upv,dnv,usea,dsea,str,chm,bot,xglu)
      Else
        Print *,'Call to MRSEB error: Mode variable out of range!',mode
      Endif 

c  The following two lines are added due to threshold effects on mass.
      If (BOT .lt. 0.) BOT=0.D0
      If (CHM .lt. 0.) CHM=0.D0

      RETURN
      END
C

      SUBROUTINE STRC78(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc78
 
C     THIS IS THE NORMAL KMRS(B0) FIT  (190 MeV)    USEA=DSEA=2*STR
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0,5,5/
      DATA INIT/0/
 
      xsave=x  ! don't let x be altered if it's out of range!!
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(78,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      XGLU=G(3)
      USEA=G(4)
      CHM=G(5)
      STR=G(6)
      BOT=G(7)
      DSEA=G(8)

      x=xsave  !restore x
 
      RETURN
      END
 
C
      SUBROUTINE STRC79(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc79
 
C     THIS IS THE NEW B0-TYPE FIT "S0"  WITH UBAR=DBAR (215 MeV)
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0,5,5/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(79,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      XGLU=G(3)
      USEA=G(4)
      CHM=G(5)

      STR=G(6)
      BOT=G(7)
      DSEA=G(8)
 
      x=xsave
 
      RETURN
      END
C
      SUBROUTINE STRC80(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc80
 
C     THIS IS THE NEW B0-TYPE FIT "D0"  WITH UBAR < DBAR (215 MeV)
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0,5,5/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(80,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      XGLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C
      SUBROUTINE STRC81(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc81
 
C     THIS IS THE NEW B--TYPE FIT "D-"  WITH UBAR < DBAR (215 MeV)
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,4,5,0,0,5,5/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(81,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      XGLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C                                    ======================  End of MRS program
      SUBROUTINE STRC82(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc82
 
C     THIS IS THE NEW  FIT "S0'" WITH UBAR = DBAR (230 MeV)
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,10,0,0,10,10/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
c io no 97 is being used because 82 is used by the package
      READ(97,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      XGLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C
      SUBROUTINE STRC83(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc83
 
C     THIS IS THE NEW B0-TYPE FIT "D0'"  WITH UBAR < DBAR (230 MeV)
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,10,0,0,10,10/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
c io no 98 is being used because 83 is being used by the package
      READ(98,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      XGLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C
      SUBROUTINE STRC84(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,XGLU)
C                                                    strc84
 
C     THIS IS THE NEW B--TYPE FIT "D-'"  WITH UBAR < DBAR (230 MeV)
 
      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,7,0,0,7,7/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
c io no. is 99 because 84 is being used by the program
      READ(99,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      XGLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END
C                       ****************************

      SUBROUTINE STRC33(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C                                                     -=-=- mrsaaj
C                                                    strc33

C************************ MRS A : Mode=8 ***********************C
C								C
C								C
C     This is the package for the new MRS(A) parton             C
C     distributions. The minimum Q^2  value is 5 GeV^2          C
C     although a version which evolves backwards down to        C
C     0.625 GeV^2 is also available. The x range is, as before  C
C     10^-5 < x < 1. MSbar factorization is used.               C
C     The package reads a grid, which is in a separate file,    C  
C     (for033.dat, ftn33, ...). Note that x times the parton    C
C     distribution is returned, Q is the scale in GeV,          C
C     and Lambda(MSbar,nf=4) = 230 MeV. MODE is not used here.  C
C								C
C     NOTE! An analytic (Duke-Owens style) approximation        C
C     is also now available. It is more compact, but takes      C
C     much longer to output the pdf values. Contact WJS if      C
C     interested.                                               C
C								C
C         The reference  to MRS(A) is :                         C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         University of Durham preprint  DTP/94/34  (1994)      C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE  GRID IS 0.00383              C
C								C
C***************************************************************C
C*********************** MRS A' : Mode=9 ***********************C
C*********************** MRS J : Mode=10 ***********************C
C ( 5/15/96 Liang                                               C
C   I ignore J' since it is not fitted to DIS/others data.      C
C   I don't consider it is useful. )                            C                         
C								C
C     This is a package for the new MRS(J,J') parton            C
C     distributions. The minimum Q^2  value is 5 GeV^2 and the  C
C     x range is, as before, 10^-5 < x < 1. MSbar factorization C
C     is used. The package reads 4 grids, which are in separate C
C     files  A'=for030.dat/ftn30, J=for037.dat/ftn37            C  
C            J'=for038.dat/ftn38                                C  
C     Note that x times the parton distribution is returned,    C
C     Q is the scale in GeV,                                    C
C     and Lambda(MSbar,nf=4) = 231 MeV for A'                   C
C                            = 344 MeV for J                    C
C                            = 507 MeV for J'                   C
C								C
C	MODE=1 for MRS(A')                                      C
C	MODE=2 for MRS(J)                                       C
C	MODE=3 for MRS(J')                                      C
C								C
C         The reference is :  E.W.N. Glover,                    C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         Durham preprint DTP/96/22 (1996)                      C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE 30 GRID IS 0.00341            C
C         THE FIRST NUMBER IN THE 37 GRID IS 0.00356            C
C         THE FIRST NUMBER IN THE 38 GRID IS 0.00150            C
C								C
C***************************************************************C
C     THIS IS THE NEW  "A" FIT -- May 1994 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(33,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END

      SUBROUTINE STRC30(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C                                                    strc30

C     THIS IS THE NEW  "Aprime" FIT -- Feb 1995 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(30,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
      END

      SUBROUTINE STRC37(X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
C                                                    strc37

C     THIS IS THE NEW  "J" FIT -- Mar 1996 -- standard Q^2 range

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(nx=47)
      parameter(ntenth=21)
      DIMENSION F(8,NX,20),G(8),XX(NX),N0(8)
      DATA XX/1.d-5,2.d-5,4.d-5,6.d-5,8.d-5,
     .        1.D-4,2.D-4,4.D-4,6.D-4,8.D-4,
     .        1.D-3,2.D-3,4.D-3,6.D-3,8.D-3,
     .        1.D-2,2.D-2,4.D-2,6.D-2,8.D-2,
     .     .1D0,.125D0,.15D0,.175D0,.2D0,.225D0,.25D0,.275D0,
     .     .3D0,.325D0,.35D0,.375D0,.4D0,.425D0,.45D0,.475D0,
     .     .5D0,.525D0,.55D0,.575D0,.6D0,.65D0,.7D0,.75D0,
     .     .8D0,.9D0,1.D0/
      DATA XMIN,XMAX,QSQMIN,QSQMAX/1.D-5,1.D0,5.D0,1310720.D0/
      DATA N0/2,5,5,9,0,0,9,9/
      DATA INIT/0/

      xsave=x
 
      IF(INIT.NE.0) GOTO 10
      INIT=1
      DO 20 N=1,nx-1
      DO 20 M=1,19
      READ(37,50)F(1,N,M),F(2,N,M),F(3,N,M),F(4,N,M),F(5,N,M),F(7,N,M),
     .          F(6,N,M),F(8,N,M)
C 1=UV 2=DV 3=GLUE 4=UBAR 5=CBAR 7=BBAR 6=SBAR 8=DBAR
         DO 25 I=1,8
  25     F(I,N,M)=F(I,N,M)/(1.D0-XX(N))**N0(I)
  20  CONTINUE
      DO 31 J=1,NTENTH-1
      XX(J)=DLOG10(XX(J))+1.1D0
      DO 31 I=1,8
      IF(I.EQ.7) GO TO 31
      DO 30 K=1,19
  30  F(I,J,K)=DLOG(F(I,J,K))*F(I,ntenth,K)/DLOG(F(I,ntenth,K))
  31  CONTINUE
  50  FORMAT(8F10.5)
      DO 40 I=1,8
      DO 40 M=1,19
  40  F(I,nx,M)=0.D0
  10  CONTINUE
      IF(X.LT.XMIN) X=XMIN
      IF(X.GT.XMAX) X=XMAX
      QSQ=SCALE**2
      IF(QSQ.LT.QSQMIN) QSQ=QSQMIN
      IF(QSQ.GT.QSQMAX) QSQ=QSQMAX
      XXX=X
      IF(X.LT.1.D-1) XXX=DLOG10(X)+1.1D0
      N=0
  70  N=N+1
      IF(XXX.GT.XX(N+1)) GOTO 70
      A=(XXX-XX(N))/(XX(N+1)-XX(N))
      RM=DLOG(QSQ/QSQMIN)/DLOG(2.D0)
      B=RM-DINT(RM)
      M=1+IDINT(RM)
      DO 60 I=1,8
      G(I)= (1.D0-A)*(1.D0-B)*F(I,N,M)+(1.D0-A)*B*F(I,N,M+1)
     .    + A*(1.D0-B)*F(I,N+1,M)  + A*B*F(I,N+1,M+1)
      IF(N.GE.ntenth) GOTO 65
      IF(I.EQ.7) GOTO 65
          FAC=(1.D0-B)*F(I,ntenth,M)+B*F(I,ntenth,M+1)
          G(I)=FAC**(G(I)/FAC)
  65  CONTINUE
      G(I)=G(I)*(1.D0-X)**N0(I)
  60  CONTINUE
      UPV=G(1)
      DNV=G(2)
      USEA=G(4)
      DSEA=G(8)
      STR=G(6)
      CHM=G(5)
      GLU=G(3)
      BOT=G(7)
 
      x=xsave
 
      RETURN
C				*******************
      END
C
      subroutine mrsX
     > (x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu,Nset,Nyr)
C***************************************************************C
C								C                 Nyr = 96   MRSR
C     This is a package for the new MRS(R1,R2,R3,R4) parton     C
C     distributions. There are several important changes from   C
C     earlier MRS packages:                                     C
C       -- the q**2 range is enlarged to 1.25d0 < q**2 < 1d7,   C
C          the x range is still 1d-5 < x < 1d0                  C
C       -- the interpolation routine has been slightly modified C
C       -- the call is now to mrs96() rather than to MRSEB()    C 
C     Note that the grid files which the program reads in       C
C     (mrsr1.dat,...) are now larger and more obviously named.  C
C								C
C     As before, x times the parton distribution is returned,   C
C     q is the scale in GeV, MSbar factorization is assumed,    C
C     and Lambda(MSbar,nf=4) = 241 MeV for R1 (mode=11)         C
C                            = 344 MeV for R2 (mode=12)         C
C                            = 241 MeV for R3 (mode=13)         C
C                            = 344 MeV for R4 (mode=14)         C
C								C
C         The reference is:                                     C
C         A.D. Martin, R.G. Roberts and W.J. Stirling,          C
C         University of Durham preprint DTP/96/44 (1996)        C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C             >>>>>>>>  CROSS CHECK  <<<<<<<<                   C
C                                                               C
C         THE FIRST NUMBER IN THE R1 GRID IS 0.00150            C
C         THE FIRST NUMBER IN THE R2 GRID IS 0.00125            C
C         THE FIRST NUMBER IN THE R3 GRID IS 0.00181            C
C         THE FIRST NUMBER IN THE R4 GRID IS 0.00085            C
C								C
C***************************************************************C
C****************************************************************C
C								 C                  Nyr = 98     MRST
C     This is a package for the new MRS 1998 parton              C
C     distributions. The format is similar to the previous       C
C     (1996) MRS-R series.                                       C
C								 C
C     As before, x times the parton distribution is returned,    C
C     q is the scale in GeV, MSbar factorization is assumed,     C
C     and Lambda(MSbar,nf=4) is given below for each set.        C
C								 C
C     TEMPORARY NAMING SCHEME:                                   C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C  1     mrs981  central gluon, a_s    300      0.1175   0.00561  C
C  2     mrs982  higher gluon          300      0.1175   0.00510  C
C  3     mrs983  lower gluon           300      0.1175   0.00408  C
C  4     mrs984  lower a_s             229      0.1125   0.00586  C
C  5     mrs985  higher a_s            383      0.1225   0.00410  C
C						                 C
C						                 C
C      The corresponding grid files are called mrs981.dat etc.    C
C							  	 C
C      The reference is:                                         C
C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
C      Univ. Durham preprint DTP/98/??, hep-ph/??????? (1998)    C
C                                                                C
C      Comments to : W.J.Stirling@durham.ac.uk                   C
C                                                                C
C								 C
C****************************************************************C
C****************************************************************C
C								 C                  Nyr = 99
C     This is a package for the new **corrected** MRST parton    C
C     distributions. The format is similar to the previous       C
C     (1998) MRST series.                                        C
C								 C
C     NOTE: 7 new sets are added here, corresponding to shifting C
C     the small x HERA data up and down by 2.5%, and by varying  C
C     the charm and strange distributions, and by forcing a      C
C     larger d/u ratio at large x.                               C
C								 C
C     As before, x times the parton distribution is returned,    C
C     q is the scale in GeV, MSbar factorization is assumed,     C
C     and Lambda(MSbar,nf=4) is given below for each set.        C
C								 C
C     NAMING SCHEME:                                             C
C						                 C
C  mode  set    comment             L(4)/MeV  a_s(M_Z)  grid#1   C
C  ----  ---    -------             --------  -------   ------   C
C								 C
C  1     COR01  central gluon, a_s    300      0.1175   0.00537  C
C  2     COR02  higher gluon          300      0.1175   0.00497  C
C  3     COR03  lower gluon           300      0.1175   0.00398  C
C  4     COR04  lower a_s             229      0.1125   0.00585  C
C  5     COR05  higher a_s            383      0.1225   0.00384  C
C  6     COR06  quarks up             303.3    0.1178   0.00497  C
C  7     COR07  quarks down           290.3    0.1171   0.00593  C
C  8     COR08  strange up            300      0.1175   0.00524  C
C  9     COR09  strange down          300      0.1175   0.00524  C
C  10    C0R10  charm up              300      0.1175   0.00525  C
C  11    COR11  charm down            300      0.1175   0.00524  C
C  12    COR12  larger d/u            300      0.1175   0.00515  C
C						                 C
C      The corresponding grid files are called cor01.dat etc.    C
C							  	 C
C      The reference is:                                         C
C      A.D. Martin, R.G. Roberts, W.J. Stirling, R.S Thorne      C
C      Univ. Durham preprint DTP/99/64, hep-ph/9907231 (1999)    C
C                                                                C
C      Comments to : W.J.Stirling@durham.ac.uk                   C
C                                                                C
C								 C
C****************************************************************C

      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,ntenth=23,np=8)
      character*20 mrsfile(10,20)
      real*8 f(np,nx,nq+1),qq(nq),xx(nx),g(np),n0(np)
      integer  Nver(10)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data n0/3,4,5,9,9,9,9,9/
      data init1, init2 / 0, 0 /
      data 
     > (Nver(I), I=1,3) / 96, 98, 99 /
     > (mrsfile(1,I), I=1,4) 
     >          /'mrsr1.dat','mrsr2.dat','mrsr3.dat','mrsr4.dat'/
     > (mrsfile(2,I), I=1,5)
     >          / 'mrs981.dat','mrs982.dat','mrs983.dat'
     >           ,'mrs984.dat','mrs985.dat' /
     > (mrsfile(3,I), I=1,12)
     >         / 'cor01.dat','cor02.dat','cor03.dat','cor04.dat'
     >          ,'cor05.dat','cor06.dat','cor07.dat','cor08.dat'
     >          ,'cor09.dat','cor10.dat','cor11.dat','cor12.dat' /

      save
      xsave=x
      q2save=qsq

      if(init1.eq.Nset .and. init2.eq.Nyr) goto 10

        Do I = 1,3
          If (Nyr .Eq. Nver(I)) Goto 3
        EndDo
        Print *, 
     > 'Nyr argument in mrsX should be (96,98,99); your Nyr =', Nyr
        Stop

  3     Nv = I
        IU= NextUn()
        open(unit=IU,file=mrsfile(Nv,Nset),status='old')

        do 20 n=1,nx-1
        do 20 m=1,nq
        read(IU,50)f(1,n,m),f(2,n,m),f(3,n,m),f(4,n,m),
     .		  f(5,n,m),f(7,n,m),f(6,n,m),f(8,n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
      do 25 i=1,np
  25	 f(i,n,m)=f(i,n,m)/(1d0-xx(n))**n0(i)
  20  continue
      do 31 j=1,ntenth-1
      xx(j)=dlog10(xx(j)/xx(ntenth))+xx(ntenth)
      do 31 i=1,8
      if(i.eq.5.or.i.eq.7) goto 31
      do 30 k=1,nq
  30  f(i,j,k)=dlog10(f(i,j,k)/f(i,ntenth,k))+f(i,ntenth,k)
  31  continue
  50  format(8f10.5)
      do 40 i=1,np
      do 40 m=1,nq
  40  f(i,nx,m)=0d0

      init1=Nset
      init2=Nyr
      close(IU)

  10  continue
      if(x.lt.xmin) x=xmin
      if(x.gt.xmax) x=xmax
      if(qsq.lt.qsqmin)	qsq=qsqmin
      if(qsq.gt.qsqmax)	qsq=qsqmax
      xxx=x
      if(x.lt.xx(ntenth)) xxx=dlog10(x/xx(ntenth))+xx(ntenth)
      n=0
  70  n=n+1
      if(xxx.gt.xx(n+1)) goto 70
      a=(xxx-xx(n))/(xx(n+1)-xx(n))
      m=0
  80  m=m+1
      if(qsq.gt.qq(m+1)) goto 80
      b=(qsq-qq(m))/(qq(m+1)-qq(m))
      do 60 i=1,np
      g(i)= (1d0-a)*(1d0-b)*f(i,n,m)   + (1d0-a)*b*f(i,n,m+1)
     .	  +       a*(1d0-b)*f(i,n+1,m) +       a*b*f(i,n+1,m+1)
      if(n.ge.ntenth) goto 65
      if(i.eq.5.or.i.eq.7) goto 65
        fac=(1d0-b)*f(i,ntenth,m)+b*f(i,ntenth,m+1)
 	  g(i)=fac*10d0**(g(i)-fac)
  65  continue
      g(i)=g(i)*(1d0-x)**n0(i)
  60  continue
      upv=g(1)
      dnv=g(2)
      usea=g(4)
      dsea=g(8)
      str=g(6)
      chm=g(5)
      glu=g(3) 
      bot=g(7)
        x=xsave
        qsq=q2save
      return
      end
      
      subroutine mrst2001(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      data Icall / 0 /

      If (Icall .le. 10) then 
        Print *, 'Call to dummy MRST2001 routine.'
        Print *, 'It occupies too much space in the library.'
        Print *, 'MRST2001.o must be linked manually when needed!'
      EndIf

      return
      end
      subroutine mrst23(x,q,mode,upv,dnv,usea,dsea,str,chm,bot,glu)
C=====================================================================
C WKT 04.09.28
C The two recent packages from MRST have been combined into a simplified
C package that (i) allows linking to both 2002 and 2003 (which the
C original versions did not, because of duplicating subroutine names;
C (ii) streamlined the function calls by merging the four copies of the
C same main subroutine (2 copis each of mrst1,2) into a single one.

C  Mode=1 gives the 2002 NLO set with alpha_s(M_Z,NLO) = 0.1197      C
C  Mode=2 gives the 2002 NNLO set with alpha_s(M_Z,NNLO) = 0.1154    C
C  Mode=3 gives the 2003 NLO set with Lambda(4) = 278 MeV            C
C  Mode=4 gives the 2003 NNLO set with lambda(4) = 231.8 MeV         C
C======================================================================
C      Original MRST headers:
C***************************************************************C
C								C
C  This is a package for the new MRST 2002 updated NLO and      C
C  NNLO parton distributions.                                   C
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0211080                                  C
C                                                               C
C  There are 2 pdf sets corresponding to mode = 1, 2            C
C                                                               C
C  Mode=1 gives the NLO set with alpha_s(M_Z,NLO) = 0.1197      C
C  This set reads a grid whose first number is 0.00949          C
C                                                               C
C  Mode=2 gives the NNLO set with alpha_s(M_Z,NNLO) = 0.1154    C
C  This set reads a grid whose first number is 0.00685          C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
C***************************************************************C
C                                                               C
C  This is a package for the MRST 2003 'conservative' NLO and   C
C  NNLO parton distributions.                                   C
C  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
C  R.S. Thorne, hep-ph/0307262                                  C
C                                                               C
C  There are 2 pdf sets corresponding to mode = 1, 2            C
C                                                               C
C  Mode=1 gives the NLO set with Lambda(4) = 278 MeV            C
C  This set reads a grid called mrst2003cnlo.dat                C
C  whose first number is 0.01057                                C
C                                                               C
C  Mode=2 gives the NNLO set with lambda(4) = 231.8 MeV         C
C  This set reads a grid called mrst2003cnnlo.dat               C
C  whose first number is 0.00563                                C
C                                                               C
C  Both fits are variants of the MRST2002 fits, but here the    C
C  range of DIS structure function data fitted is restricted to C
C  x > 0.005, W^2 > 15 GeV^2 and Q^2 > 10 (7) GeV^2 for the NLO C
C  (NNLO) fits.                                                 C
C                                                               C
C         Comments to : W.J.Stirling@durham.ac.uk               C
C                                                               C
C***************************************************************C
c-------------------------------------------------------------------
c  HLL 08/31/06
c  Remove subrountines/functions 
c     jeppe1, jeppe2, locx, polderiv 
c  which duplicate with LHAPDF
c-------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      q2=q*q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(x.lt.xmin.or.x.gt.xmax)       print 98,x

      DATA IU / 53 /
C      IU = NextUn()
      if    (mode.eq.1) then
        open(unit=IU,file='mrst2002nlo.dat',status='old')
      elseif(mode.eq.2) then
        open(unit=IU,file='mrst2002nnlo.dat',status='old')
      elseif(mode.eq.3) then
        open(unit=IU,file='mrst2003cnlo.dat',status='old')
      elseif(mode.eq.4) then
        open(unit=IU,file='mrst2003cnnlo.dat',status='old')
      endif
      call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,IU)

  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
      return
      end

      subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm,bot,glu,NU)
      implicit real*8(a-h,o-z)
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26)
      real*8 f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq),
     .f6(nx,nq),f7(nx,nq),f8(nx,nq),fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4),
     .cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4),
     .ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,
     .	      1d-4,2d-4,4d-4,6d-4,8d-4,
     .	      1d-3,2d-3,4d-3,6d-3,8d-3,
     .	      1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,
     .	   .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,
     .	   .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,
     .	   .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,
     .	   .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,
     .        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,
     .        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,
     .        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,
     .        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      data init/0/
      save
      xsave=x
      q2save=qsq
      if(init.ne.0) goto 10
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(NU,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),
     .		  f5(n,m),f7(n,m),f6(n,m),f8(n,m)
c notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea
  20  continue
      close(NU)
      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
  40  continue
      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      call cjeppe1(nx,nq,xxl,qql,f1,cc1)
      call cjeppe1(nx,nq,xxl,qql,f2,cc2)
      call cjeppe1(nx,nq,xxl,qql,f3,cc3)
      call cjeppe1(nx,nq,xxl,qql,f4,cc4)
      call cjeppe1(nx,nq,xxl,qql,f6,cc6)
      call cjeppe1(nx,nq,xxl,qql,f8,cc8)

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      !call cjeppe1(nx,nqc,xxl,qqlc,fc,ccc)

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      !call cjeppe1(nx,nqb,xxl,qqlb,fb,ccb)


      init=1
   10 continue

      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call cjeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
      call cjeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
      call cjeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
      call cjeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
      call cjeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
      call cjeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)

      chm=0.d0
      if(qsq.gt.emc2) then
      call cjeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then
      call cjeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
      endif

      x=xsave
      qsq=q2save
      return
   50 format(8f10.5)
      end

      FUNCTION PartonX (IPRTN, X, Q)
C
C   Given the parton distribution function in the array Upd in
C   COMMON / CtqPar1 / , this routine fetches u(fl, x, q) at any value of
C   x and q using Mth-order polynomial interpolation for x and Ln(Q/Lambda).
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      PARAMETER (M= 2, M1 = M + 1)
C
      Logical First
      Common 
     > / CtqPar1 / Al, XV(0:MXX), QL(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
C
      Dimension Fq(M1), Df(M1)

      Data First /.true./
      save First
C                                                 Work with Log (Q)
      QG  = LOG (Q/AL)

C                           Find lower end of interval containing X
      JL = -1
      JU = Nx+1
 11   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (X .GT. XV(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif

      Jx = JL - (M-1)/2
      If (X .lt. Xmin .and. First ) Then
         First = .false.
         Print '(A, 2(1pE12.4))', 
     >     ' WARNING: X < Xmin, extrapolation used; X, Xmin =', X, Xmin
         If (Jx .LT. 0) Jx = 0
      Elseif (Jx .GT. Nx-M) Then
         Jx = Nx - M
      Endif
C                                    Find the interval where Q lies
      JL = -1
      JU = NT+1
 12   If (JU-JL .GT. 1) Then
         JM = (JU+JL) / 2
         If (QG .GT. QL(JM)) Then
            JL = JM
         Else
            JU = JM
         Endif
         Goto 12
      Endif

      Jq = JL - (M-1)/2
      If (Jq .LT. 0) Then
         Jq = 0
         If (Q .lt. Qini)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q < Qini, extrapolation used; Q, Qini =', Q, Qini
      Elseif (Jq .GT. Nt-M) Then
         Jq = Nt - M
         If (Q .gt. Qmax)  Print '(A, 2(1pE12.4))', 
     >     ' WARNING: Q > Qmax, extrapolation used; Q, Qmax =', Q, Qmax
      Endif

      If (Iprtn .GE. 3) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
C                             Find the off-set in the linear array Upd
      JFL = Ip + NfMx
      J0  = (JFL * (NT+1) + Jq) * (NX+1) + Jx
C
C                                           Now interpolate in x for M1 Q's
      Do 21 Iq = 1, M1
         J1 = J0 + (Nx+1)*(Iq-1) + 1
         Call Polint (XV(Jx), Upd(J1), M1, X, Fq(Iq), Df(Iq))
 21   Continue
C                                          Finish off by interpolating in Q
      Call Polint (QL(Jq), Fq(1), M1, QG, Ftmp, Ddf)

      PartonX = Ftmp
C
      RETURN
C                        ****************************
      END
      FUNCTION Pdf (Jset, Ihadron, Iparton, X, Q, Ir)

C========================================================================
C GroupName: Pdfs
C Description: user callable parton distribution functions
C ListOfFiles: pdf pdfh pdfp pdftst
C========================================================================
C $Header: /Net/u52/wkt/h22/2prz/RCS/Pdfs.f,v 7.1 99/08/20 22:07:38 wkt Exp $
C $Log:	Pdfs.f,v $
c Revision 7.1  99/08/20  22:07:38  wkt
c common PdfSwh fixed
c 
c Revision 7.0  99/08/20  10:56:38  wkt
c Renamed Version 7.0 since common block /PdfSwh/ has changed;
c and in order to synchronize with EvlPac7. Iset = 902 (used .ini) added.
c 
c Revision 6.5  1999/08/19  21:23:24  wkt
c Add Iset = 902 option to evolve from .ini file in non-interactive mode
c and with user setable grids.
c
c Revision 6.4  1999/08/13  11:19:14  wkt
c pdfp.f modified to synchronize with new setpdf concerning the added
c option of Iset=Isetin0=900 -- cf. comments in Setpdf.f
c
c Revision 6.3  97/12/17  23:03:41  wkt
c Iset 1410 extended to 1414
c 
c Revision 6.2  97/11/16  00:18:14  wkt
c Revised pdf, pdfh, pdfp : added IpdMod,Iptn0 switches in / PdfSwh /
c See comments in pdf.f module.
c 
c Revision 6.1  97/11/15  18:07:17  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c 
C========================================================================

C These callable parton distributions are, in order of generality:

C      FUNCTION Pdf (Jset, Ihadron, Iparton, X, Q, Ir)
C         It is now only a shell to make it compatible to existing 
C         programs calling the PDF function in this format.

C      FUNCTION PdfH (Ihadron, Iparton, X, Q)
C         This is mainly for Hadron targets other than the proton.

C      FUNCTION PdfP (LPRTN, XD, QD)
C         Parton distribution in proton.

C ========================================================================

C      FUNCTION Pdf (Jset, Ihadron, Iparton, X, Q, Ir)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Common / PdfSwh / Iset, IpdMod, Iptn0, NuIni

C      Iset = Jset		This is now commented out to avoid confusion in some
C						Unusual situations.
C		It is very important to make sure that
C       	Setpdf(Iset) 
C		is called before a specific set of pdf is to be used.

C                           The (IpdMod, Iptn0) switches turn on/off 
C                           specified parton flavors. Cf. pdfp module.

      Pdf = PdfH (Ihadron, Iparton, X, Q)

      Return
C                  *****************************
      End

      FUNCTION PdfH (Ihadron, Iparton, X, Q)
C     -------------------------------------------------
C     Parton Distribution Functions inside Hadrons.

C     Revised 4/4/94 by HLL & WKT: 
C        PdfH retains its name and argument list for compatibility with all
C        existing programs which Call this function; PdfH switches between
C        target hadrons; 
C        PdfP is for proton target; it switches between different Iset's.
C
C     Ihadron = -1,  0,  1,  2,  3,  4,  5,  6 : 
C              pbar, n,  p,  D, Cu,  C, ~D, Fe 

C               5 is "isoscalar-corrected iron" hence = D

C         In all cases, adjust the Iparton label
C         to convert to the corresponding proton distribution which is
C         given in Function PdfP;

C     ---------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     Character Msg*80

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / EVLPAC / AL, IKNL, IPD0, IHDN, NfMx

      Common / PdfSwh / Iset, IpdMod, Iptn0, NuIni

      DATA
     1   HUGE, DUM,  D0m,   D1p 
     1 / 1D10, 0.D0, -1D-6, 1.000001D0 /
     1   IW1,  IW2    / 2*0  /
C                                                -------  Check x-range
      IF (X.LE.D0m .or. X.GT.D1p) Then
        Call WARNR(IW1,NWRT,'X out of range in PdfH.', 'X', X,
     >               D0, D1, 1)
        TEM = DUM
      EndIf
C                   --- Conversion of  Ihadron  to proton distributions,
C                                     if necessary ---- 
      If (Ihadron.gt.6 .or. Ihadron.lt.-1) then
         Call WARNI(IW, NWRT,
     >    'Only Ihardon=-1,0,1,2,...,6 (pbar,n,p,D,*,) are active',
     >    'Ihadron', Ihadron, -1,6,1)
       PdfH=0.D0
       Return
      Endif
C                                 
      Jp=Abs(Iparton)
      Neff = NFL(Q)
c nfl(q) returns the number of `light' flavors at scale Q - effective
      If ( Jp .gt. NEFF) then
c 		if Jp > Neff, then set PdfH=0 and return
c         Call WARNI(IW, NWRT,
c     >    'Iparton out of range',
c     >    'Iparton', Iparton, -Neff,Neff,1)
	 PdfH = 0D0
         Return
      Endif

      If (Jp.eq.1 .or. Jp.eq.2) Then
C                                   Use Isospin symmetry n<->p  == u<->d
         Ipartner=3-Jp
         If (Iparton.lt.0) Ipartner=-Ipartner

         If (Ihadron.eq.1) then
            Tem= PdfP(Iparton, X, Q)
         Elseif (Ihadron.eq.-1) then
            Tem= PdfP(-Iparton, X, Q)
         Elseif (Ihadron.eq.0) then
            Tem= PdfP(Ipartner, X, Q)
         Elseif (Ihadron.eq.2 .or.Ihadron.eq.4 .or.Ihadron.eq.5) then
            Tem=( PdfP(Iparton, X, Q)
     >           +PdfP(Ipartner, X, Q) )/2.D0
         Elseif (Ihadron.eq.3) then
            Tem=( 29.D0* PdfP(Iparton, X, Q)
     >           +35.D0* PdfP(Ipartner, X, Q) )/64.D0
         Elseif (Ihadron.eq.6) then
            Tem=( 26.D0* PdfP(Iparton, X, Q)
     >           +30.D0* PdfP(Ipartner, X, Q) )/56.D0
         Endif
      Else
         Tem= PdfP(Iparton, X, Q)
      Endif
C                                      --- Make sure PdfH >= 0 --------
C                                       (unless Iknl<0 - polarized pdf)
      IF (TEM .LT. D0 .and. Iknl .ge. 1) Then
        IF (TEM .LT. D0m .AND. X .LE. 0.55D0) Then
        Call WARNR(IW2,NWRT,'PdfH < 0; Set --> 0', 'PDF',TEM,D0,HUGE,1)
CCPY        WRITE (NWRT, '(A, 2I5, 2(1PE15.3))')
CCPY     >      ' ISET, Iparton, X, Q = ', ISET, Iparton, X, Q
        EndIf
        TEM = D0
      EndIf
C                        -------- Return function value and error code
      PdfH = TEM

      RETURN
C                        ****************************
      END
      FUNCTION PdfP (LPRTN, XD, QD)
C
C $Header: /home/wkt/common/flib/2prz/RCS/pdfp.f,v 8.1 2003/06/29 21:37:19 wkt Exp $
C $Log: pdfp.f,v $
C Revision 8.1  2003/06/29 21:37:19  wkt
C Option Iset = 903 to use .pds input has been added.
C
C Revision 8.0  2003/06/29 19:44:13  wkt
C Start another round of archiving.
C
C  **** If change the Iset settings, have to change in SetPdf also.
C
C    This routine gives the parton ( IPARTN ) distribution function inside
C    the proton in a chosen evolved or parametrized form ( Iset )

C  2. For Iset = Isetev0 , IsetevX, or Isetini:
C                steer to results from QCD evolution program;

C    It gives the probability distribution, not the momemtum-weighted one.
C

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      Parameter
     >(Isetev0=10, Isetev1=11,
     > Isetin0=900, Isetin1=901, Isetin2=902, Isetin3=903,
     > Isettbl=911)

      CHARACTER MSG*75
C
      COMMON
     >/ IOUNIT / NIN, NOUT, NWRT
     >/ PdfSwh / Iset, IpdMod, Iptn0, NuIni
C                           IpdMod: 0 normal mode
C                                   1 only flavor +/- Iptn0 is turned on
C                                  -1      flavor +/- Iptn0 is turned off
C
     >/ XQold / Xold, Qold  ! LHAPDF
      Dimension f(-6:6)     ! LHAPDF
      DATA IW1 / 0 /
      Data Xold, Qold /-1d0, -1d0 /
CJI 2021 Add in lhapdf
      integer initlhapdf, iset_save
      common/lhapdf/initlhapdf, iset_save

      IER = 0
C
      Iprtn = LPRTN

      If((IpdMod.eq.1).and.(Iprtn.Ne.Iptn0.and.Iprtn.Ne.-Iptn0)) Then
         PdfP = 0.
         Return
      Elseif (IpdMod.eq.-1.and.(Iprtn.Eq.Iptn0.or.Iprtn.Eq.-Iptn0)) Then
         PdfP = 0
         Return
      EndIf
      x = XD
      Q = QD
C                              ====  Begin of overall Iset IfBlock =====
      if(initlhapdf==1) then
        if(iprtn == 1) then
            iprtn = 2
        elseif(iprtn == 2) then
            iprtn = 1
        elseif(iprtn == -1) then
            iprtn = -2
        elseif(iprtn == -2) then
            iprtn = -1
        endif
        call lhapdf_xfxq(1, iset_save, iprtn, x, q, tmp);
        tmp = tmp/x
      elseIf     (Iset .eq.  0) Then
C                                          -----   Experimental set ------

         Tmp = PdfTst (Iprtn, X, Q)
      Elseif (Iset.eq.Isetev0 .or.
     >  (Iset.Ge.Isetin0 .and. Iset.Le.Isetin3)) then
C                                                     evolved: 10, 900,901,902
                            Tmp = ParDis (Iprtn, X, Q)
      Elseif (Iset .eq. Isetev1) then
C                                                             11 evolved
                            Tmp = PrDis1 (Iprtn, X, Q)
      Elseif (Iset .EQ. 801) then
c                                                             801 EHLQ1
                            Tmp = PSHK (Iprtn, X, Q, 1)
      Elseif (Iset .EQ. 802) then
c                                                             802 Duke-Owens 1
                            Tmp = PSHK (Iprtn, X, Q, 3)
      ElseIF (Iset .EQ. 803) then
c                                                             803 DFLM
                            Tmp = DFLM (Iprtn, X, Q, 2)
      ElseIF (Iset .EQ. 804) then
c                                                             804 MT90-S1
                            Tmp = Pdxmt (1, Iprtn, X, Q, ir)
      ElseIf (Iset.gt.2000 .and. Iset.lt.3000) then
c                                                             2001- MRS
                            Nset = Iset - 2000
                            Tmp = PDMRS (Iprtn, X, Q, Nset)
      Elseif (Iset.gt.3000 .and. Iset.lt.4000) then
c                                                             3001- GRV
c  3001:LO; 3002:NLO(MSbar); 3003:NLO(DIS)
                            Nset = Iset - 3000
                            Tmp = PDGRV (Iprtn, X, Q, Nset)
      Elseif (Iset.gt.1100 .and. Iset.le.1106) then
c                                                             1101- CTEQ1
                            Nset = Iset-1100
                            tmp = Ctq1Pd(Nset, Iprtn, X, Q, ir)
      Elseif (Iset.gt.1200 .and. Iset.le.1206)  then
c                                                             1201- CTEQ2
                            Nset = Iset - 1200
                            Tmp = Ctq2Pd(Nset, Iprtn, X, Q, ir)
      ElseIf (Iset.gt.1300 .and. Iset.le.1303) then
c                                                             1301- CTEQ3
                            Nset = Iset - 1300
                            Tmp = Ctq3Pd(Nset, Iprtn, X, Q, ir)
      Elseif (Iset.gt.1400 .and. Iset.le.1414) then
c                                                             1401- CTEQ4
                            Tmp = Ctq4Pdf (Iprtn, X, Q)
      Elseif (Iset.eq. Isettbl) then
                            Tmp = Ctq4Pdf (Iprtn, X, Q)

      Elseif (Iset.gt.1500 .and. Iset.le.1509) then
c                                                       1501-07 CTEQ5 tbl
                            Tmp = Ctq5Pdf (Iprtn, X, Q)
      Elseif (Iset.Eq.1511) then
c                                                  CTEQ5M parametrized
                            Tmp = Ctq5Pd(1, Iprtn, X, Q, ir)
      Elseif (Iset.Eq.1512) then
c                                                  CTEQ5L parametrized
                            Tmp = Ctq5Pd(3, Iprtn, X, Q, ir)
      Elseif ((Iset>=1600.and.Iset<=1699).or.
     >        (Iset>=1900.and.Iset<=1999)) then
c                                                  1601-1699 CTEQ6 tbl/pds
c                                                  1900-1999 CTEQ6.5 pds
                            Tmp = Ctq6Pdf (Iprtn, X, Q)
      Elseif (Iset.ge.11900 .and. Iset.le.11952) then !ZL
        Tmp=CT14Pdf(Iprtn, X, Q)

      Elseif (Iset.ge.20000 .and. Iset.le.99999) then
c                                                       LHAPDF
         If (X.ne.Xold .or. Q.ne.Qold) then
            Call evolvePDF(X,Q,f)
            Xold=X
            Qold=Q
         Endif
         If(Iprtn==1 .or. Iprtn==2) then
                Tmp=f(3-Iprtn)/X
         Elseif(Iprtn==-1 .or. Iprtn==-2) then
                Tmp=f(-3-Iprtn)/X
         Else
                Tmp=f(Iprtn)/X
         Endif

      Else
         MSG=
     >'Iset chosen is currently inactive. PdfP set equal to zero.'
         CALL WARNI (IW1, NWRT, MSG, 'Iset', Iset, 1, 9999, 0)
                            Tmp = 0.D0
      EndIf

      PdfP = Tmp

   10 RETURN
C                        ****************************
      END
      FUNCTION PdfTst (IPRTN, XX, QQ)
C                                        Place-holder for a test pdf

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Print *, 
     >'You just called a dummy test function PdfTst in PrtPdf!'

      PdfTst = 0
C
      RETURN
C                        ****************************
      END

       Function PDGRV(Iprtn, X, Q, Nset)
C                                                   -=-=- pdgrv
c                                                             61- GRV
c  61:LO; 62:NLO(MSbar); 63:NLO(DIS)

       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
          Q2= Q * Q
          if(Nset.eq.1) then
             Call GRV94LO (X, Q2, UV, DV, DEL, UDB, SB, GL)
          elseif(Nset.eq.2) then
             Call GRV94HO (X, Q2, UV, DV, DEL, UDB, SB, GL)
          elseif(Nset.eq.3) then
             Call GRV94DI (X, Q2, UV, DV, DEL, UDB, SB, GL)
          Endif

          If (Iprtn.eq.2) then
             PDGRV = DV + (UDB + DEL) / 2D0
          ElseIf (Iprtn.eq.1) then
             PDGRV = UV + (UDB - DEL) / 2D0
          ElseIf (Iprtn.eq.0) then
             PDGRV = GL
          ElseIf (Iprtn.eq.-1) then
             PDGRV = (UDB - DEL) / 2D0
          ElseIf (Iprtn.eq.-2) then
             PDGRV = (UDB + DEL) / 2D0
          ElseIf (Iprtn.eq.-3 .or. Iprtn.eq.3) then
             PDGRV = SB
          Else
             PDGRV = 0D0
          Endif

          PDGRV = PDGRV / X

          Return
C     ===============================================
          End

      FUNCTION PDMRS (LPRTN, X, Q, NST)
C                                                   -=-=- pdmrs
C
C===========================================================================
C GroupName: Mrs
C Description: MRS parton distributions
C ListOfFiles: pdmrs mrseb strcxxx mrsr mrs98 mrs98x 
C===========================================================================
C #Header: /Net/d2a/wkt/1hep/2pdf/prz/RCS/Mrs.f,v 1.1 98/12/22 23:04:02 wkt Exp $
C #Log:	Mrs.f,v $
c Revision 1.1  98/12/22  23:04:02  wkt
c Initial revision
c 
c Revision 6.1  97/11/15  18:07:13  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c 
C========================================================================

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      DATA IWRN / 0 /
C
      MST = NST
      XX = X
      QQ = Q
C
      CALL MRSEB (XX, QQ, MST, UV, DV, USEA, DSEA, STR, CHM, BTM, GLN)
C
      IF     (LPRTN .EQ. 0) THEN
                TEM = GLN
      ElseIF (LPRTN .EQ. 1) THEN
                TEM = UV + USEA
      ElseIF (LPRTN .EQ. 2) THEN
                TEM = DV + DSEA
      ElseIF (LPRTN .EQ. -1) Then
                TEM = USEA
      ElseIF (LPRTN .EQ. -2) THEN
                TEM = DSEA
      ElseIF (LPRTN .EQ. 3 .OR. LPRTN .EQ. -3) THEN
                TEM = STR
      ElseIF (LPRTN .EQ. 4 .OR. LPRTN .EQ. -4) THEN
                TEM = CHM
      ElseIF (LPRTN .EQ. 5 .OR. LPRTN .EQ. -5) THEN
                TEM = BTM
      Else
      CALL WARNI(IWRN, NWRT, 'MRS only has 5 flavors', 'Iparton',
     >           LPRTN, -5, 5, 1)
                TEM = 0
      EndIf
                PDMRS = TEM / X
      RETURN
C                        ****************************
      END
C                                                                 MAY 30 90
      FUNCTION PDXMT (ISET, IPARTON, X, Q, IRT)
C                                                   -=-=- pdxmt

C             For ISET = 1, 2 .. , returns sets of Parton Distributions
C	      (in the proton) with parton label Iparton (6, 5, ...,0, ...-6)
C	      for (t, b, c, s, d, u, g, u-bar, ... t-bar), and kinematic
C	      variables (X, Q).   IRT is a return error code.
C
C     Iset =  1, 2, 3, 4 corresponds to the S1, B1, B2, and E1 fits of Morfin-
C	      Tung (Fermilab-Pub-90/24, IIT-90/11) to NLO in the DIS scheme.
C
C	      5 (Set S1M) corresponds to the same set as 1 (S1) but expressed
C	      in the MS-bar scheme.
C	      
C             All the above sets assume a SU(3)-symmetric sea.

C	      6 is the CTEQ Singular distribution
C
C	      7 is the CTEQ Small Lambda distribution
C
C	      8 is the CTEQ Normal no experimental normalization
c               constraint fit
c
C	      9 is the CTEQ Normal with experimental normalization
c               constraint fit
C
C  The "lambda" parameter (4-flavors) for each parton distribution set can be
C  obtained by making the following FUNCTION call:
C             Alam = Vlambd (Iset, Iorder)
C  where Iset is the (input) set #, Iorder is the (output) order of the fit (1
C  for set 9, 2 for all the others), and Alam is the value of the effective QCD
C  lambda for 4 flavors.

C             Details about the 1 - 5 distributions are
C	      given in the above-mentioned preprint.
C	      
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION  THRSLD(0:6)

      DATA (THRSLD(I), I=0,6) / 4*0.0, 1.5, 5.0, 180.0 /

      IFL = IPARTON
      JFL = ABS(IFL)
C                                                   Return 0 if below threshold
      IF (Q .LE. THRSLD(JFL)) THEN
        PDXMT = 0.0
        RETURN
      EndIf
C                                                                       Valence
      IF (IFL .LE. 0) THEN
         VL = 0
      ElseIF (IFL .LE. 2) THEN
         VL = PDZXMT(ISET, IFL, X, Q, IRT)
      Else
         VL = 0
      EndIf
C                                                                         Sea
      SEA = PDZXMT (ISET, -JFL, X, Q, IRT)
 
      PDXMT = VL + SEA
       
      RETURN
C                         *************************
      END
      
      FUNCTION PDZMT (IU, LP, XX, QQ, IRT)
C                                                   -=-=- pdzmt

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      IFL = LP
C                                                                       Valence
      IF (IFL .LE. 0) THEN
        VL = 0
      ElseIF (IFL .LE. 2) THEN
        VL = PDZZMT(IU, IFL, XX, QQ, IRT)
      Else
        VL = 0
      EndIf
C                                                                         Sea
      KP = -ABS (IFL)
      SEA = PDZZMT (IU, KP, XX, QQ, IRT)

      PDZMT = VL + SEA

      RETURN
C                         *************************
      END
C
      FUNCTION PDZXMT (IST, LP, XX, QQ, IRT)
C                                                   -=-=- pdzxmt

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NEX = 3, MXFL = 6, NPN = 3, NST = 30)
      parameter (nexnu = 5)

      DIMENSION
     1 AC(0:NEXNU, 0:NPN, -MXFL:2, NST), A(0:NEXNU), 
     1 ALM(NST), Q02(NST), MEX(NST), MPN(NST), MQRK(NST), Iord(NST)

      DATA MEX, MPN, MQRK / 10*NEX,20*NEXNU, NST*2, NST*6 /
C                                         Set S-1:    PDF parameters from /L352
      DATA IORD(1), ALM(1), Q02(1) / 2, 0.212, 4.00 /
     > (((AC(IEX,IPN,IFL,1), IEX=0,3), IPN=0,2), IFL=-6,2)
     > /  1.34,  0.15,
     >    5.30, -1.96, -0.57,  0.16,  0.43,  1.08, -0.08, -0.02,
     >    0.06, -0.03,  1.62,  0.11,  3.68, -1.94, -0.33,  0.14,
     >    0.53,  0.87, -0.10, -0.01,  0.03,  0.02,  1.88, -0.33,
     >    7.52, -1.34, -2.78,  0.10, -1.13,  2.92,  0.13, -0.04,
     >    0.04, -0.49, -0.99, -0.33,  8.53, -1.55, -1.54,  0.03,
     >   -1.08,  2.02,  0.10, -0.03,  0.39, -0.39, -0.99, -0.33,
     >    8.53, -1.55, -1.54,  0.03, -1.08,  2.02,  0.10, -0.03,
     >    0.39, -0.39, -0.99, -0.33,  8.53, -1.55, -1.54,  0.03,
     >   -1.08,  2.02,  0.10, -0.03,  0.39, -0.39, -3.98, -0.15,
     >    7.46,  0.35,  0.72, -0.06,  0.96,  0.89, -0.63,  0.00,
     >   -0.30, -0.04, -6.28, -0.18,  6.56,  0.65,  2.62,  0.02,
     >    1.40,  1.13, -1.18, -0.03, -0.38, -0.16,-13.08, -0.40,
     >   15.35, -0.43,  8.54,  0.31,-11.83,  3.18, -2.70, -0.12,
     >    4.16, -0.82 /
C                                          Set B1:    PDF parameters from /L212
      DATA IORD(2), ALM(2), Q02(2) / 2, 0.194, 4.00 /
     > (((AC(IEX,IPN,IFL,2), IEX=0,3), IPN=0,2), IFL=-6,2)
     > /  1.30,  0.19,
     >    5.24, -1.81, -0.57,  0.15,  0.44,  1.06, -0.09, -0.02,
     >    0.05, -0.02,  1.59,  0.14,  3.65, -1.81, -0.34,  0.13,
     >    0.53,  0.86, -0.10, -0.01,  0.03,  0.02,  1.48, -0.14,
     >    6.75, -0.50, -2.49, -0.11, -0.54,  2.13,  0.04,  0.03,
     >   -0.15, -0.24, -1.08, -0.13,  8.40, -0.88, -1.33, -0.21,
     >   -0.51,  1.18, -0.03,  0.06,  0.07, -0.05, -1.08, -0.13,
     >    8.39, -0.88, -1.33, -0.21, -0.50,  1.18, -0.03,  0.06,
     >    0.07, -0.05, -1.08, -0.13,  8.39, -0.88, -1.33, -0.21,
     >   -0.50,  1.18, -0.03,  0.06,  0.07, -0.05, -4.22, -0.02,
     >    7.29,  0.90,  0.88, -0.17,  1.08,  0.50, -0.69,  0.03,
     >   -0.39,  0.08, -6.42, -0.09,  6.47,  1.03,  2.67, -0.03,
     >    1.39,  1.00, -1.21, -0.02, -0.42, -0.14,-12.92, -0.36,
     >   15.74, -0.30,  8.33,  0.32,-12.73,  3.35, -2.68, -0.13,
     >    4.51, -0.91 /
C                                            Set B2:  PDF parameters from /L261
      DATA IORD(3), ALM(3), Q02(3) / 2, 0.191, 4.00 /
     > (((AC(IEX,IPN,IFL,3), IEX=0,3), IPN=0,2), IFL=-6,2)
     > /  1.38,  0.18,
     >    5.40, -1.91, -0.59,  0.16,  0.42,  1.11, -0.08, -0.02,
     >    0.06, -0.03,  1.64,  0.09,  3.74, -2.02, -0.33,  0.14,
     >    0.54,  0.88, -0.10, -0.01,  0.03,  0.02,  1.52, -0.72,
     >    7.75, -2.18, -2.71,  0.45, -1.56,  3.75,  0.15, -0.15,
     >    0.16, -0.76, -0.85, -0.82,  9.19, -2.76, -1.43,  0.35,
     >   -0.92,  2.56, -0.03, -0.09,  0.12, -0.40, -0.85, -0.82,
     >    9.19, -2.76, -1.43,  0.35, -0.92,  2.56, -0.03, -0.10,
     >    0.12, -0.40, -0.85, -0.82,  9.19, -2.76, -1.43,  0.35,
     >   -0.92,  2.56, -0.03, -0.10,  0.12, -0.40, -3.74, -0.58,
     >    9.63, -1.09,  0.21,  0.24, -1.13,  2.10, -0.50, -0.07,
     >    0.25, -0.33, -6.07, -0.52,  8.33, -0.52,  2.33,  0.22,
     >    0.28,  1.91, -1.15, -0.07, -0.28, -0.31,-12.08, -0.73,
     >   21.14, -1.92,  7.31,  0.54,-19.17,  4.59, -2.35, -0.18,
     >    6.64, -1.25 /
C                                            Set E1:  PDF parameters from /L152
      DATA IORD(4), ALM(4), Q02(4) / 2, 0.155, 4.00 /
     > (((AC(IEX,IPN,IFL,4), IEX=0,3), IPN=0,2), IFL=-6,2)
     > /  1.43,  0.16,
     >    6.17, -1.94, -0.65,  0.16,  0.43,  1.12, -0.08, -0.02,
     >    0.06, -0.02,  1.69,  0.11,  3.69, -1.99, -0.33,  0.14,
     >    0.54,  0.90, -0.11, -0.01,  0.03,  0.02,  2.11, -0.33,
     >    7.93, -1.51, -3.01,  0.10, -1.40,  3.14,  0.18, -0.04,
     >    0.09, -0.55, -0.84, -0.32,  8.96, -1.70, -1.65,  0.02,
     >   -1.24,  2.15,  0.12, -0.03,  0.45, -0.43, -0.84, -0.32,
     >    8.96, -1.70, -1.65,  0.02, -1.24,  2.15,  0.12, -0.03,
     >    0.45, -0.43, -0.84, -0.32,  8.96, -1.70, -1.65,  0.02,
     >   -1.24,  2.15,  0.12, -0.03,  0.45, -0.43, -3.87, -0.15,
     >    7.83,  0.21,  0.85, -0.07,  1.00,  0.93, -0.73,  0.00,
     >   -0.36, -0.03, -6.09, -0.17,  6.75,  0.54,  2.81,  0.01,
     >    1.74,  1.15, -1.34, -0.03, -0.56, -0.16,-12.56, -0.38,
     >   14.62, -0.41,  8.69,  0.30,-11.27,  3.19, -2.93, -0.12,
     >    4.29, -0.87 /
C                               Set S1M:  PDF parameters from /L352 -- MS-Bar
      DATA IORD(5), ALM(5), Q02(5) / 2, 0.212, 4.00 /
     > (((AC(IEX,IPN,IFL,5), IEX=0,3), IPN=0,2), IFL=-6,2)
     > /  1.75,  0.11,
     >    6.20, -2.35, -1.02,  0.26, -0.41,  1.68,  0.05, -0.06,
     >    0.29, -0.24,  2.03,  0.06,  4.43, -2.35, -0.78,  0.24,
     >   -0.18,  1.52,  0.03, -0.04,  0.22, -0.19,  1.09, -0.24,
     >    5.97, -0.64, -2.41,  0.08, -0.90,  2.71, -0.12,  0.02,
     >   -0.35, -0.20, -0.14, -0.49, 10.24, -2.57, -1.98,  0.02,
     >   -1.43,  2.32,  0.23, -0.02,  0.44, -0.47, -0.14, -0.49,
     >   10.24, -2.57, -1.98,  0.02, -1.44,  2.32,  0.23, -0.02,
     >    0.45, -0.47, -0.15, -0.49, 10.23, -2.57, -1.98,  0.02,
     >   -1.44,  2.32,  0.23, -0.02,  0.45, -0.47, -2.36, -0.49,
     >    9.00, -1.74, -1.42,  0.44, -0.46,  3.93,  0.21, -0.22,
     >    0.29, -1.34, -2.19, -1.07, 11.30, -4.85, -3.86,  1.56,
     >   -7.20, 10.51,  1.57, -0.73,  3.85, -4.36,-24.77,  7.52,
     >  -99.51, 36.02,-23.00,  0.48,-16.45, 16.51, 34.44, -6.26,
     >   97.19,-40.40 /

C                                            Set B0:  PDF parameters from /P154
      DATA IORD(9), ALM(9), Q02(9) / 1, 0.144, 4.00 /
     > (((AC(IEX,IPN,IFL,9), IEX=0,3), IPN=0,2), IFL=-6,2)
     > /  1.38,  0.16,
     >    5.40, -1.97, -0.62,  0.19,  0.59,  1.24, -0.10, -0.02,
     >    0.03, -0.05,  1.67,  0.08,  3.75, -2.09, -0.33,  0.17,
     >    0.70,  0.98, -0.13, -0.01,  0.00,  0.02,  1.52, -0.25,
     >    7.01, -0.79, -3.17, -0.01, -0.90,  2.90,  0.25,  0.00,
     >   -0.08, -0.54, -0.81, -0.07,  9.19, -0.89, -1.13, -0.46,
     >    0.35,  0.33, -0.26,  0.16, -0.49,  0.40, -0.81, -0.07,
     >    9.19, -0.89, -1.13, -0.46,  0.35,  0.33, -0.26,  0.16,
     >   -0.49,  0.40, -0.81, -0.07,  9.19, -0.89, -1.13, -0.46,
     >    0.35,  0.33, -0.26,  0.16, -0.49,  0.40, -3.62, -0.06,
     >    8.30,  0.16,  0.03, -0.21, -0.60,  1.26, -0.48,  0.05,
     >    0.25, -0.15, -6.16, -0.11,  6.49,  0.71,  2.37, -0.05,
     >    1.28,  1.37, -1.24, -0.02, -0.41, -0.26,-12.68, -0.35,
     >   14.87, -0.17,  8.36,  0.28,-12.56,  3.39, -2.89, -0.12,
     >    4.75, -0.96 /


      IRT = 0
      IFL = -4-LP
      
      X  = XX
      Q0 = SQRT (Q02(IST))
      ALAM = ALM(IST)
      SQ = LOG ( LOG(QQ/ALAM) / LOG(Q0/ALAM) )

      if (ist.ge.10) then
         ipnno = 3
      Else
         ipnno = 2
      EndIf

      DO 20 IEX = 0, MEX(IST)
          A(IEX) = AC(IEX, 0, IFL, IST)
          DO 21 IPN = 1, ipnno
 21          A(IEX) = A(IEX) + AC(IEX, IPN, IFL, IST) * SQ **IPN
 20   continue

c if ist >= 10, then a new fit with the Modified Traditional parametrization
c is used, otherwise a MT param. is correct.

      if(ist.ge.10) then
         pdf = exp(a(0))*x**(a(1))*(1.-x)**(a(2))
     $        *(1. + (exp(a(3)) - 1.) * (x**(a(4))) )
     $        *(log(1.+1./x))**(a(5))
      Else
         pdf = exp(a(0))*x**(a(1))*(1.-x)**(a(2))
     $        *(log(1.+1./x)**(a(3)))
      EndIf


      PDZXMT = PDF / X        
c parton distributions, Not momentum distributions

      RETURN
C                                                -----------------------------
      ENTRY VLAMBD (ISET, IORDER)

      IORDER = IORD (ISET)
      VLAMBD = ALM  (ISET)

      RETURN
C                         *************************
      END

      FUNCTION PDZZMT (IU, LP, XX, QQ, IRT)
C                                                   -=-=- pdzzmt
C                                                                 Jan 15 90
C                                                                 DEC 25 89
C                Parton Distribution for flavor IFL inside the proton at (X,Q)
C             by computing from coefficient parameters in files Xnnn.cfn
C
C            IU :    Fortran Unit # of a file containing the parameter values
C                                                   and table of coefficients

C            IRT:    Return codes for error and warning conditions
C                 0  O.K.
C                 1  Q value < QMIN,  Q is set equal to QMIN;
C                 2  Q value > QMAX,
C                 3  IFL out of range, Zero returned
C
C                        ----------------------------

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      CHARACTER CH6*6, LIN8*80

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NEX = 4, MXFL = 6, NPN = 3, NST = 10)

      DIMENSION
     1 AC(0:NEX, 0:NPN, -MXFL:2, NST), A(0:NEX), FX(0:NEX),
     1 IUN(NST), ALM(NST), Q02(NST), MEX(NST), MPN(NST), MQRK(NST)

      DATA NU / 0 /

      IRT = 0
      IFL = LP
C                                                       Initiation, section
C                                            If coefficients not read, do so;
C                             If already read, determine value of IST to use
      IF (NU .GE. 1) THEN
        DO 7 IST = 1, NU
          IF (IU .EQ. IUN(IST)) THEN
            GOTO 8
          EndIf
    7  CONTINUE
      EndIf

        NU = NU + 1
        IUN(NU) = IU
        IST = NU
C                                            Read in coefficients from Xnnn.cfn
        READ (IU, '(A)') LIN8                                             Title
           PRINT '(2A)', LIN8
        READ (IU, '(A)') LIN8                                             Label

        READ (IU, *, ERR=101)
     1     ALM(NU), NFL, NORD, NT, XMIN, Q02(NU), QMAX,
     1     MEX(NU), MPN(NU), MQRK(NU)
        QMIN = SQRT(Q02(NU))
           PRINT '(/A, I1, A, F6.3, A,I2, A/)',
     >  ' Order-', NORD, ' Lamda =', ALM(NU), ' for ', NFL, ' flavors.'
           PRINT '(A, 3(1PE11.3) / 2A)',
     >' Nominal kinematic range: Xmin, Qmin, Qmax =', XMIN,QMIN,QMAX,
     >' Parametrization provides natural smooth extrapolation',
     1' beyond this range.'

        READ (IU, '(A)') (LIN8, I=1,7)                                   Labels

        DO 11 IEX = 0, MEX(NU)
          READ (IU, '(A)') LIN8                                           Label
          READ (IU, '(A6, 9F7.2)', ERR=101)
     1    (CH6, (AC(IEX,IPN,IF,NU),IF=2,-MQRK(NU),-1),IPN=0,MPN(NU))
   11   CONTINUE
C                                                  --------------------------
    8 CONTINUE

      X  = XX
      Q0 = SQRT (Q02(IST))
      ALAM = ALM(IST)
      SQ = LOG ( LOG(QQ/ALAM) / LOG(Q0/ALAM) )

      FX(0) = EXP(D1)
      FX(1) = X
      FX(2) = 1.- X
      FX(3) = LOG (1.+ 1./X)

      PDF = 1.
      DO 20 IEX = 0, MEX(IST)
        A(IEX) = AC(IEX, 0, IFL, IST)
        DO 21 IPN = 1, MPN(IST)
          A(IEX) = A(IEX) + AC(IEX, IPN, IFL, IST) * SQ **IPN
21      CONTINUE
        PDF = PDF * FX(IEX) **(A(IEX))
   20 CONTINUE

      PDZZMT = PDF / X

      RETURN
  101 STOP 'Error reading coefficients from data table in PDZZMT!!'
C                         *************************
      END
C----------------------------------------------------------------------------

      FUNCTION PSHK (LPARTN, X, Q, NSET)
C                                                   -=-=- pshk
C
C===========================================================================
C GroupName: DoEhlq
C Description: Peskin's package of EHLQ and Duke-Owens pdf's
C ListOfFiles: pshk xfd xfu xfg xfsea xfs xfc xfb xft tcheby
C===========================================================================
C #Header: /Net/cteq06/users/wkt/1hep/2pdf/pdf/v6/RCS/DoEhlq.f,v 6.1 97/11/15 18:06:41 wkt Exp $
C #Log:	DoEhlq.f,v $
c Revision 6.1  97/11/15  18:06:41  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c 
C========================================================================

C                                  Original Author of this section is M. Peskin
C                                  New data tables for 1986 version of EHLQ
C                                  courtesy of E. Eichten
C
C       Revised by J. Qiu in March 1989:
C
C       1, Qmin^2 = 4 GeV^2 for Duke_Owens distributions,
C                           not 5 GeV^2 as in old version;
C       2, ADOpting Ed Berger's higher precision parameters
C                           for Duke_Owens distributions.
C------------------------------------------------------------------------------
C       DIRECTIONS FOR PARTON:
C                                          MAY 28, 1984
C
C    THIS FILE CONTAINS PARMETRIZATIONS OF THE PROTON STRUCTURE
C     FUNCTIONS RECENTLY OBTAINED BY
C          EICHTEN, HINCHLIFFE, LANE, AND QUIGG
C               AND BY
C           DUKE AND OWENS
C
C
C    THESE STRUCTURE FUNCTIONS ARE FIT TO DEEP INELASTIC SCATTERING
C     DATA AT CURRENT VALUES OF Q**2 AND EXTRAPOLATED TO HIGH Q**2
C     BY INTEGRATING THE ALTERELLI-PARISI EQUATIONS.  ONLY THE
C     LEADING-ORDER TERMS IN THE EVOLUTION KERNEL ARE INCLUDED; THIS
C     CORRESPONDS TO THE LEADING-LOG APPROXIMATION.
C
C   FOR A DETAILED DISCUSSION OF THE PARAMETRIZATIONS
C      AND THEIR DERIVATION,          SEE:
C          E. EICHTEN, ET. AL., FERMILAB PREPRINT  (1984)
C             AND
C         D. W. DUKE AND J. F. OWENS, PHYS REV D30, 49 (1984)
C
C  THE EHLQ CLAIM THAT THEIR PARAMETRIZATION IS VALID OVER THE RANGE
C       SQRT(Q**2) =  SQRT(5) GEV   TO   SQRT(Q**2) = 10**4 GEV
C         X =  10**(-4)  TO 1
C   DUKE AND OWENS ARE NOT SO PRECISE ABOUT LIMITS OF VALIDITY
C       THEY CLAIM THAT THEIR PARAMETRIZATION IS VALID AT THE FEW
C         PERCENT LEVEL FOR SQRT(Q**2) UP TO ABOUT 1 TEV AND
C             OVER THE BULK OF THE X RANGE FROM 0 TO 1
C
C   EACH GROUP PROVIDES TWO SETS OF STRUCTURE FUNCTIONS, CORRESPONDING
C     TO TWO DIFFERENT VALUES OF LAMBDA.
C        (NOTE: THIS IS LAMBDA (LEADING ORDER), NOT LAMBDA (MSBAR))
C     THE TWO SETS ARE GENERATED WITH THE FOLLOWING VALUES OF LAMBDA:
C         EHLQ:       NSET = 1:  LAMBDA = 200 MEV,  QMIN^2 = 5 GEV^2
C                     NSET = 2:  LAMBDA = 290 MEV,  QMIN^2 = 5 GEV^2
C         D AND O:    NSET = 3:  LAMBDA = 200 MEV,  QMIN^2 = 4 GEV^2
C                     NSET = 4:  LAMBDA = 400 MEV,  QMIN^2 = 4 GEV^2
C
C
C    ACCESS THIS SET OF STRUCTURE FUNCTIONS IN THE FOLLOWING WAY:
C
C       THIS FILE CONTAINS A FUNCTION SUBROUTINE CORRESPONDING TO
C         EACH INDEPENDENT STRUCTURE FUNCTION  XF(X):
C
C              XF(X) FOR:              CALL THE FUNCTION:
C
C      U QUARKS                        XFU(X,Q)
C      D QUARKS                        XFD(X,Q)
C      LIGHT SEA QUARKS:
C             U, D, U BAR, OR D BAR    XFSEA(X,Q)
C      S OR S BAR                      XFS(X,Q)
C      C OR C BAR                      XFC(X,Q)
C      B OR B BAR                      XFB(X,Q)
C      T OR T BAR                      XFT(X,Q)
C      GLUONS                          XFG(X,Q)
C
C              (AND, FOR EP ENTHUSIASTS:
C      PHOTONS IN THE ELECTRON         XFP(X,ECM)
C         WHICH SIMPLY COMPUTES THE WEISZACKER-WILLIAMS PHOTON SPECTRUM)
C
C
C     (THE MOMENTUM SUM RULE READS:
C
C           1 = INT DX (XFU + XFD + 2 XFSEA
C                           + 2 (XFS + XFC + XFB + XFT) + XFG). )
C
C        DUKE AND OWENS ASSUME THAT XFS = XFSEA  AND DO NOT PROVIDE
C           VALUES OF XFB AND XFT
C
C    THESE FUNCTIONS TAKE THE ARGUMENTS:
C           X = LONGITUDINAL FRACTION
C           Q = SQRT(Q**2)(IN GEV)
C     THE FUNCTIONS ARE WRITTEN WITH DOUBLE PRECISION ARITHMETIC
C
C !!! MAKE SURE THAT Q IS IN THE RANGE  SQRT(5) < Q < 10**4 !!!
C !!!  MAKE SURE THAT X IS IN THE RANGE  10**-4  < X < 1   !!!
C
C       HEAVY QUARK STRUCTURE FUNCTIONS RETURN THE VALUE 0
C         UNTIL  Q = QMIN (OF ORDER 2 * M QUARK)
C               THE HEAVY QUARK MASSES ASSUMED ARE
C                    MB = 5.5 GEV  ;  MT = 30 GEV
C
C !!!  TO CHOOSE A SET OF STRUCTURE FUNCTIONS, SET  NSET = 1, 2, 3 OR 4;
C             (THE DEFAULT IS NSET = 2;    LAMBDA= 290 MEV)
C        THE FUNCTIONS ALL INCLUDE A STATEMENT:
C
C     COMMON/NSET/ NSET
C
C          AFTER THAT, JUST RELAX; WE DO IT ALL !!!!
C
C               GOOD LUCK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / NSET / ISET
C
      ISET = NSET
C
      IF     (LPARTN .EQ. 0) THEN
                PSHK = XFG(X, Q)
      ElseIF (LPARTN .EQ. 1) THEN
                PSHK = XFU(X, Q)
      ElseIF (LPARTN .EQ. 2) THEN
                PSHK = XFD(X, Q)
      ElseIF (LPARTN .EQ. -1 .OR. LPARTN .EQ. -2) THEN
                PSHK = XFSEA(X, Q)
      ElseIF (LPARTN .EQ. 3 .OR. LPARTN .EQ. -3) THEN
                PSHK = XFS(X, Q)
      ElseIF (LPARTN .EQ. 4 .OR. LPARTN .EQ. -4) THEN
                PSHK = XFC(X, Q)
      ElseIF (LPARTN .EQ. 5 .OR. LPARTN .EQ. -5) THEN
                PSHK = XFB(X, Q)
      ElseIF (LPARTN .EQ. 6 .OR. LPARTN .EQ. -6) THEN
                PSHK = XFT(X, Q)
      Else
      CALL WARNI(IWRNT,NWRT,'IPARTN out of range in PSHK','IPARTN',
     > LPARTN, -6, 6, 1)
                PSHK = 0
      EndIf
                PSHK = PSHK / X
      RETURN
C                        ****************************
      END
      Subroutine ReadTbl (Nu)
      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      PARAMETER (MXX = 105, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      MxVal=2
      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
C                                               This is the old .tbl (HLL) format
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log (TV(Iq) /Al)
 11      Continue

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End
      Subroutine SetCtq4 (Iset)
C                                                   -=-=- setctq4
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=14)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq4m.tbl', 'cteq4d.tbl', 'cteq4l.tbl'
     > , 'cteq4a1.tbl', 'cteq4a2.tbl', 'cteq4m.tbl', 'cteq4a4.tbl'
     > , 'cteq4a5.tbl', 'cteq4hj.tbl', 'cteq4lq.tbl'
     > , 'cteq4hq.tbl', 'cteq4hq1.tbl', 'cteq4f3.tbl', 'cteq4f4.tbl' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
         IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
	    goto 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22       continue

         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
            Print *, 'Invalid Iset number in SetCtq4 :', Iset
            Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
         Endif
         Call ReadTbl (IU)
         Close (IU)
         Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq4!!'
      Stop
C                             ********************
      End
C                                                          =-=-= MiscPd
      Subroutine SetCtq5 (Iset)
      Implicit Double Precision (A-H,O-Z)
      Parameter (Isetmax=7)
      Character Flnm(Isetmax)*12, Tablefile*40
      Data (Flnm(I), I=1,Isetmax)
     > / 'cteq5m.tbl', 'cteq5d.tbl', 'cteq5l.tbl', 'cteq5hj.tbl'
     > , 'cteq5hq.tbl', 'cteq5f3.tbl', 'cteq5f4.tbl' /
      Data Tablefile / 'test.tbl' /
      Data Isetold, Isetmin, Isettest / -987, 1, 911 /
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
	 IU= NextUn()
         If (Iset .eq. Isettest) then
            Print* ,'Opening ', Tablefile
 21         Open(IU, File=Tablefile, Status='OLD', Err=101)
            goto 22
 101        Print*, Tablefile, ' cannot be opened '
            Print*, 'Please input the .tbl file:'
            Read (*,'(A)') Tablefile
            Goto 21
 22         continue
         ElseIf (Iset.lt.Isetmin .or. Iset.gt.Isetmax) Then
	    Print *, 'Invalid Iset number in SetCtq5 :', Iset
	    Stop
         Else
            Tablefile=Flnm(Iset)
            Open(IU, File=Tablefile, Status='OLD', Err=100)
	 Endif
         Call ReadTbl (IU)
         Close (IU)
	 Isetold=Iset
      Endif
      Return

 100  Print *, ' Data file ', Tablefile, ' cannot be opened '
     >//'in SetCtq5!!'
      Stop
C                             ********************
      End

      Subroutine SetPdf (Jset)
C                                                   -=-=- setpdf
C
C $Header: /home/wkt/common/flib/2prz/RCS/setpdf.f,v 8.2 2004/01/24 01:44:19 wkt Exp $
C $Log: setpdf.f,v $
C Revision 8.2  2004/01/24 01:44:19  wkt
C  Same as MrsEB; add Alekhin/WJS as a special case
C
C Revision 8.1  2003/06/29 21:43:17  wkt
C Option Iset = 903 to use .pds files as input added.
C
C Revision 8.0  2003/06/29 19:43:22  wkt
C Start another round of archiving
C
C===========================================================================
C GroupName: Setpdf
C Description: setup package for parametrizationd of Parton distributions
C ListOfFiles: setpdf pdfp
C===========================================================================

C Revision 1.4  2002/02/04 03:30:31  wkt
C Update to include CTEQ6 PDF's
C
C Revision 1.2  2001/11/21 04:21:59  wkt
C Added Mrst2001
C
c Revision 7.1  99/08/20  14:38:01  wkt
c Bug fix for Iset=900 option
c
c Revision 7.0  99/08/20  10:56:41  wkt
c Renamed Version 7.0 since common block /PdfSwh/ has changed;
c and in order to synchronize with EvlPac7. Iset = 902 (used .ini) added.
c
c Revision 6.13  99/08/20  10:30:56  wkt
c Use the added variable (cf. EvlPac7.2) NuIni for setting Iset = Isetin2
c (902) up.  Calling program must open the .ini file and assign Unit#=NuIni.
c
c Revision 6.12  99/08/19  21:25:04  wkt
c Improved error handling if input .ini file cannot be openned.
c
c Revision 6.11 1999/08/19  21:23:24  wkt
c Add Iset = 902 option to evolve from .ini file in non-interactive mode
c and with user setable grids.
c
c Revision 6.10  99/08/12  22:50:23  wkt
c Iset = 900 (Isetini0) option added. This will generate pdf's from a
c pre-defined file "pdf.ini" and default Nx, Xmin, Nq, Qmax without prompt.
c
c Revision 6.9  1998/12/22  23:04:28  wkt
c MRS98 added; modules re-ordered
c
c Revision 6.8  97/12/19  21:35:56  wkt
c 1410 extended to 1414
c
c Revision 6.7  97/12/19  12:23:00  wkt
c *** empty log message ***
c
c Revision 6.1  97/11/15  18:05:08  wkt
c OverAll pdf functions + parametrized pdf's (with evolution package excluded)
c SetPdf written by HHL, modified for v.6. by wkt

C========================================================================
C                                      (This file includes Block Data DatPdf)
C      Subroutine SetPdf (Jset)
C
C  **** If change the Iset order, have to change in PdfP also.
C
C                                 Force link to DatPdf + other setup work

C  =======================================================================
C                Explanation for Iset of each PDF's
C  =======================================================================
C  Setup routine for the pdfpac ...
C    for the evolution package:
C                             force link to datpdf to initialize parameters
C    for parton distributions using .ini  or  .tbl input files:
C                call the appropriate setup routines to read in input files
C    for canned CTEQ, MRS and other distributions:
C                                   setup the appropriate QCD lambda values
C  -----------------------------------------------------------------------
C  Iset assignments in this routine are coordinated with those in PdfP
C       Any changes made here must be reflected there at the same time.

C  Current Iset assignments:

C  Iset      Description       FortranFun needed     Data/Input File(s) needed
C ------------------------------------------------------------------------
C     0       Test                  --                    --
C  -----------------------------------------------------------------------
C             Evolved from
C  Isetev0,1    input fn       Fini.f (DummyArg)        (Block Data DatPdf)
C  (10, 11)

C             Evolved from
C  Isetin0   input para tbl .ini                         pdf.ini
C  (900)    (use default grid pts)
C  Isetin1   input para tbl .ini                         xxx.ini
C  (901)    (prompt for grid pts)
C  Isetin2   input para tbl .ini                         xxx.ini
C  (902)    (open xxx.ini, and set grid pts by ParPdf calls)
C  Isetin3   input pdf tbl .pds                         xxx.pds
C  (903)    (prompt for xxx.pds; read in pdf table)

C             Read in table
C  Isettbl  fr evolved results    (Tblxxx.f)            xxx.tbl
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
C  1501 -09    CTEQ5 M,D,L,
C  1511/13     CTEQ5M/L                (Ctq5Pd.f)
C  1600 -40    CTEQ6M + 40 up/down sets                   cteq6m[100-140].tbl
C  1641 -44    CTEQ6 M, D,L,L1         (Ctq6Pdf.f)        cteq6[m,d,l,l1].tbl
C  1645 -49    CTEQ6Sn                 (Ctq6Pdf.f)        cteq6[a,b,c,b+,b-].pds
C  1650 -90    CTEQ6.1M + 40 up/down sets                 ctq61.[200-240].tbl
C  1691        CTEQ6HQ                 (Ctq6Pdf.f)        cteq6hq.pds
C  1900 -40    CTEQ6.5M + 40 up/down sets (Ctq6Pdf.f)     ctq65.xx.pds  (xx = 00:40)
C  For other CTEQ5/6 distributions, use 900-902 with the corresponding.ini file
C                                    or 903 with the corresponding .pds file
C
C  2001 - 4    KMRS, MRSS0,D0,D-      Strc78 - 81.f       For078 - 81.Dat
C  2005 - 7      MRSS0',D0',D-'       Strc82 - 84.f       For097 - 99.Dat
C  2008 - 10     MRSA, A', J          Strc33,30,37.f      For033,30,37.Dat
C  2011 - 14     MRSR1 - 4            mrsr.f            mrsr1 - 4.Dat
C  2021 - 25     MRS98 1 - 5          Mrs98.f; Mrs98x.f     Ftxx.dat
C  2101 - 12     MRS99 1 - 12           mrs99.f           Cornn.dat
C  2201 - 4      MRST01 1 -4            mrst2001.f        alfnnn.dat j121.dat
C  2211 - 4      MRST02 1,2 & 2003c 1,2     mrst2001.f    mrst2002nlo.dat, ... etc.
C                             (all MRS's need PDMRS.f and MRSEB.f)
C
c  2901        Alekhin NLO (2002) as implemented by J. Stirling
C
C  3001 - 3    GRV LO,NLO,DIS ('94)  PDGRV.f + GRV94LO,HO,DI.f

C  - - - - - - - - - - - - - - Some old PDF's - - - - - - - - - - - - - -
C    801         EHLQ1                  PSHK.f
C    802         Duke-Owens 1           PSHK.f
C    803         DFLM                   DFLM.f
C    804         MT90-S1                Pdxmt.f
C

C  -----------------------------------------------------------------------
C  =======================================================================
C   More explanation for CTEQ and MRS PDF's in their respective modules
C  =======================================================================

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Parameter
     >(Isetev0=10, Isetev1=11,
     > Isetin0=900, Isetin1=901, Isetin2=902, Isetin3=903,
     > Isettbl=911)

      Character FileIni*40, Header*80
      Character*20 PARM(20)  ! LHAPDF
      Dimension VALUE(20)    ! LHAPDF

      Common
     > / PdfSwh / Iset, IpdMod, Iptn0, NuIni
     > / QCDtable /  Alambda, Nfl, Iorder
     > / XQold / Xold, Qold ! LHAPDF

      Data
     > Xmin, Qmax, Nx, Nq / .9999999D-5, 1.0000001d4, 95, 21 /
     > Len, FileIni / 3, 'pdf' /

cZL
      CHARACTER*40 PDF_FILE
      COMMON/FILINP/ PDF_FILE
      character*2 nn

      integer initlhapdf, iset_save
      COMMON/LHAPDF/initlhapdf, iset_save

      Iset = Jset
      Zm = VbnMas(4)

      If (Iset.eq.Isetev0  .or. Iset.eq.Isetev1
     > .or.Iset.eq.Isetin0 .or. Iset.eq.Isetin1
     > .or.Iset.eq.Isettbl)  Call SetEvl        ! initialize the defaults


      if(initlhapdf.eq.1) then
        print*, "here getorderas"
        Call lhapdf_getorderas(1, iset_save, Nord)
        Ord=Nord+1
        Call ParQcd(1,'ORDR',Ord,Ir)
        call lhapdf_alphasq(1, iset_save, Zm, alfs)
        Call Alfset(Zm,Alfs)

       elseIf (Iset.Ge.Isetin0 .and. Iset.Le.Isetin2) then

        If (Iset .Eq. Isetin0) Then
C                                    use pdf.ini and defaults for  Xmin, Qmax, Nx, Nq
          IU= NextUn()
          OPEN (IU, FILE='pdf.ini', STATUS='OLD', ERR=900)
          Goto 99
 900      Stop 'Cannot open the file /pdf.ini/; program stopped.'

        ElseIf (Iset .Eq. Isetin1) Then
c                                                             .ini input
         Print '(2A)', 'Please input the .ini file:   ', FileIni
 902     Read '(A)', FileIni
         Call TrmStr(FileIni, Len)
         IU= NextUn()
         OPEN (IU, FILE=FileIni(1:Len)//'.ini', STATUS='OLD', ERR=901)
         Goto 903
 901     Print *, 'Cannot open the file ', FileIni(1:Len)//'.ini'
         Print *, 'Input valid ini filename again; or ^C to stop.'
         Goto 902

 903     Print '(A, 2(1pE12.3), 2I6)',
     >   'Input Xmin, Qmax, Nx, Nq; defaults:', Xmin, Qmax, Nx, Nq
         Read (*,*) Xmin, Qmax, Nx, Nq

        ElseIf (Iset .Eq. Isetin2) Then
             Call ParPdf (2, 'Xmin', Xmin, Ir)
             Call ParPdf (2, 'Qmax', Qmax, Ir)
             Call ParPdf (2, 'Nx', Anx, Ir)
             Call ParPdf (2, 'Nt', Anq, Ir)
         Call ParPdf (2, 'NuIni', aIu, Ir)
             Nx = Nint(Anx)
             Nq = Nint(Anq)
         IU = Nint(aIu)
C                                         It is assumed that Unit #IU has been openned
C                                         for read, and it is a valid .ini file
        EndIf

 99     Call Evlini(IU, Xmin, Qmax, Nx, Nq)
        close(IU)

	  Iset = 10				! Reset Iset=10, so Pardis will be used in PDF calls.

      Elseif (Iset .eq. Isetin3) then
c         Print '(2A)', 'Please input the .pds file:   ', FileIni
c 802     Read '(A)', FileIni
c         Call TrmStr(FileIni, Len)
c         IU= NextUn()
c         OPEN (IU, FILE=FileIni(1:Len)//'.pds', STATUS='OLD', ERR=801)
c         Goto 803
c 801     Print *, 'Cannot open the file ', FileIni(1:Len)//'.pds'
c         Print *, 'Input valid .pds filename again; or ^C to stop.'
c         Goto 802
cZL
         Call TrmStr(PDF_FILE, Len)
         IU= NextUn()
         OPEN (UNIT=IU, FILE=PDF_FILE(1:Len), STATUS='OLD', ERR=801)
 803	   Call PdfSet0 (1, IHDRN, ALAM, TPMS, QINI, QMAX, XMIN,
     >             IU, HEADER, I2, I3, IRET, IRR, datpdf)
	   Print *, 'PDF values loaded; header =', Header
cZL	   Open (IU, FILE=FileIni(1:Len)//'.pds', STATUS='OLD', ERR=801)
	   Open (UNIT=IU, FILE=PDF_FILE(1:Len), STATUS='OLD', ERR=801)
         Read (IU, '(A)') (header, I = 1,2)
         Read (IU, *) order, fl, Alambda
	   Iorder = Nint(order)
	   Nfl = Nint(fl)
         Call SetLam (Nfl,Alambda,Iorder)

	   Iset = 10				! Reset Iset=10, so Pardis will be used in PDF calls.
	
      Elseif (Iset .eq. Isettbl) then
c                                                             .tbl input
         Call SetCtq4 (Iset)
         Call SetLam (Nfl,Alambda,Iorder)
      Elseif (Iset .EQ. 801) then
c                                                             801 EHLQ1
         Vlam = .2D0
         Call SetLam (4, Vlam, 2)
      Elseif (Iset .EQ. 802) then
c                                                             802 Duke-Owens 1
         Vlam = .2D0
         Call SetLam (4, Vlam, 2)
      Else IF (Iset .EQ. 803) then
c                                                             803 DFLM
         Vlam = .3D0
         Call SetLam (4, Vlam, 2)
      Else IF (Iset .EQ. 804) then
c                                                             804 MT90-S1
         Vlam = Vlambd (5, IorPdf)
         Call SetLam (4, Vlam, IorPdf)
      Else If (Iset.gt.2000 .and. Iset.lt.3000) then
c                                                             2001- MRS
         If (Iset.eq.2041) then
                Vlam = .19D0
         ElseIf (Iset.eq.2010 .or. Iset.eq.2012 .or. Iset.eq.2014) then
                Vlam = .344D0
         ElseIf (Iset.eq.2011 .or. Iset.eq.2013) then
                Vlam = .241D0
C                                                                                                     !MRS98
         ElseIf (Iset.eq.2021 .or. Iset.eq.2022 .or. Iset.eq.2023) then
                Vlam = .300D0
         ElseIf (Iset.eq.2024) then
                Vlam = .229D0
         ElseIf (Iset.eq.2025) then
                Vlam = .383D0
C                                                     !MRS99
         ElseIf (Iset.Ge.2101 .and. Iset.Le.2103) then
                Vlam = .300D0
         ElseIf (Iset.Ge.2108 .and. Iset.Le.2112) then
                Vlam = .300D0
         ElseIf (Iset.eq.2104) then
                Vlam = .229D0
         ElseIf (Iset.eq.2105) then
                Vlam = .383D0
         ElseIf (Iset.eq.2106) then
                Vlam = .3033D0
         ElseIf (Iset.eq.2107) then
                Vlam = .2903D0
C                                         !MRST2001
         ElseIf (Iset.eq.2201) then
                Vlam = 0.323D0
         ElseIf (Iset.eq.2202) then
                Vlam = 0.290D0
         ElseIf (Iset.eq.2203) then
                Vlam = 0.362D0
         ElseIf (Iset.eq.2204) then
                Vlam = 0.353D0
C                                         !MRST2002
         ElseIf (Iset.eq.2211) then
                Call AlfSet(Zm, 0.1197D0)
                Vlam = ALAMF(4)
         ElseIf (Iset.eq.2212) then
                Call AlfSet(Zm, 0.1154D0)
                Vlam = ALAMF(4)
C                                         !MRST2003
         ElseIf (Iset.eq.2213) then
                Vlam = 0.278D0
         ElseIf (Iset.eq.2214) then
                Vlam = 0.231D0
         ElseIf (Iset.eq.2901) then						! tmp: Alekhin / WJS
                Vlam = 0.323D0
         Else
                Vlam = .23D0
         Endif
         Call SetLam (4, Vlam, 2)
      Elseif (Iset.gt.3000 .and. Iset.lt.4000) then
c                                                             3001- GRV
         If (Iset.eq.3001) then
             Vlam= .232D0
             Call SetLam (3, Vlam, 1)
         Else
             Vlam= .248D0
             Call SetLam (3, Vlam, 2)
         Endif
      Elseif (Iset.gt.1100 .and. Iset.le.1106) then
c                                                             1101- CTEQ1
         Vlam = Wlambd (Iset-1100, IorPdf)
         Call SetLam (5, Vlam, IorPdf)
      Elseif (Iset.gt.1200 .and. Iset.le.1206)  then
c                                                             1201- CTEQ2
         Vlam = Wlamd2 (Iset-1200, IorPdf, 5)
         Call SetLam (5, Vlam, IorPdf)
      ElseIf (Iset.gt.1300 .and. Iset.le.1303) then
c                                                             1301- CTEQ3
         Vlam = Wlamd3 (Iset-1300, IorPdf, 5)
         Call SetLam (5, Vlam, IorPdf)
      Elseif (Iset.gt.1400 .and. Iset.le.1414) then
c                                                             1401- CTEQ4
c                             Call SetCtq4 to initialize
         Call SetCtq4 (Iset-1400)
         Call SetLam (Nfl,Alambda,Iorder)
      Elseif (Iset.gt.1500 .and. Iset.le.1509) then
c                                                   1501-07 CTEQ5 tables
         Call SetCtq5 (Iset-1500)
         Call SetLam (Nfl,Alambda,Iorder)
      Elseif (Iset .Eq. 1511) then
c                                                  CTEQ5M parametrized
         Call SetLam (5 , 0.226D0, 2)
      Elseif (Iset .Eq. 1512) then
c                                                  CTEQ5L parametrized
         Call SetLam (5 , 0.146D0, 1)
      Elseif (Iset.ge.1600 .and. Iset.le.1699) then
c                                                   1601 - 1650 CTEQ6 tables
         If (Iset.ge.1600 .and. Iset.le.1640) then
            Call SetCtq6 (Iset-1500)    ! for Iset=1600-1640, SetCtq6 expects 100-140
C                                         CTEQ6M is 1600
         ElseIf (Iset.ge.1641 .and. Iset.le.1644) then
            Call SetCtq6 (Iset-1640)    ! for Iset=1641-1644, SetCtq6 expects 1-4 (M,D,L,L1)
         ElseIf (Iset.ge.1645 .and. Iset.le.1649) then
            Call SetCtq6 (Iset-1634)    ! for Iset=1644-1649, SetCtq6 expects 11-15 (a,b,c,b+,b-))
         ElseIf (Iset.ge.1650 .and. Iset.le.1690) then
            Call SetCtq6 (Iset-1450)    ! for Iset=1650-1690, SetCtq6 expects 200-240
         ElseIf (Iset.eq.1691) then
C                                         CTEQ6.1M is 1650
            Call SetCtq6 (Iset-1670)                           ! CTEQ6HQ, SetCtq6 expects 21
         EndIf

         Call SetLam (Nfl,Alambda,Iorder)
      ElseIf (Iset>=1900 .and. Iset<=1940) then
         Call SetCtq6 (Iset-1600) ! for 6.5M series, SetCtq6 expects 300-340
         Call SetLam (Nfl,Alambda,Iorder)
      ElseIf (Iset>=1950 .and. Iset<=1957) then
         Call SetCtq6 (Iset-1920) ! for 6.5S series, SetCtq6 expects 30-37
         Call SetLam (Nfl,Alambda,Iorder)
      ElseIf (Iset>=1960 .and. Iset<=1966) then
         Call SetCtq6 (Iset-1920) ! for 6.5C series, SetCtq6 expects 40-46
         Call SetLam (Nfl,Alambda,Iorder)
C==============
      ElseIf (Iset>=11900 .and. Iset<=11952) then !ZL
CZL         Call InitialPDF (Iset-1900) ! Use CT14PDF.f for CT10 PDFs
         if(Iset<=11909) then
           Write(nn,'(I1)') Iset-11900
           PDF_FILE='pdf0'//nn(1:1)//'.pds'
         else
           Write(nn,'(I2)') Iset-11900
           PDF_FILE='pdf'//nn(1:2)//'.pds'
         endif 
         Call SetCT14(PDF_FILE(1:9)) ! Use CT14PDF.f for CT10 PDFs
CCPY March 2014
         Call CT14GetPars(xmin,Qini,Qmax,Nloops,Nfl)
        print *,'xmin,Qini,Qmax,Nloops,Nfl=',xmin,Qini,Qmax,Nloops,Nfl
CCPY MARCH 2014: FOR ALPHA_S(MZ) = 0.118 
C  Ordr, Nfl, lambda = 2.   5.  0.2262   
C         Call SetLam (5,0.2262d0,2) !0.118
C         Call SetLam (5,0.2018d0,2) !0.116
C         Call SetLam (5,0.2526d0,2) !0.120
C TO MAKE IT ASLO WORK AT THE NNLO, I MODIFIED THE FUNCTION ALPI AS:
C ALPI=CT14Alphas(AMU)/PI 
C==============      
      Elseif (Iset.ge.20000 .and. Iset.le.99999) then  ! LHAPDF
        Call SetLHAPARM('SILENT')
        Call SetLHAPARM('NOSTAT')
        PARM(1)='DEFAULT'
        VALUE(1)=Iset
C        Call PDFSET(PARM,VALUE)
        Xold=-1d0   ! reset Xold and Qold to avoid mismatched at pdfp
        Qold=-1d0
        Call GetOrderAs(Nord)
        Ord=Nord+1
        Call ParQcd(1,'ORDR',Ord,Ir)
        Alfs=AlphasPDF(Zm)
        Call Alfset(Zm,Alfs)
      Elseif (Iset.Ne.0 .and. Iset.Ne.Isetev0 .and.
     >          Iset.Ne.Isetev1) then
         Print *,
     >   'Illegal argument in Setpdf(Iset) call; Iset =', Iset
         Stop
      Endif

      Return
 801	print *, "pdf.pds name is not right"
	stop
C                        ****************************
      END
      SUBROUTINE TCHEBY(T,X)
C                                                   -=-=- tcheby
C      THIS SUBROUTINE GENERATES THE FIRST 6 TCHEBYSHEV POLYNOMIALS,
C        EVALUATED AT THE POINT X
C      IT ASSIGNS TX(N) = T SUB (N-1) (X)
C          USING THE CONVENTIONS OF ABRAMOWITZ AND STEGUN:
C                 T SUB N (COS A) = COS(NA)
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DIMENSION T(6)
      XX = X**2
      T(1) =  1.0D0
      T(2) =  X
      T(3) =  2.0D0*XX - 1.0D0
      T(4) = (4.0D0*XX - 3.0D0)*X
      T(5) = 8.0D0*(XX- 1.0D0)*XX + 1.0D0
      T(6) = ((16.0D0*XX - 20.0D0)*XX + 5.0D0)* X
      RETURN
      END
C                       ****************************
C
C                                                          =-=-= Mrs
      FUNCTION XFB(X,Q)
C                                                   -=-=- xfb
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2),TMIN(2),TMINN(2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.2D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (LAMBDA) 101,101,104
  101 CONTINUE
C  101 PRINT 102
  102 FORMAT('   Duke-Owens has no B_ or T_quark; '/
     * '         I HAVE SET NSET = 2 (EHLQ SET 2) for XFB')
      NSET = 2
      LAMBDA = 0.29D0
  104 TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C         IN THIS CASE, A LARGER TMIN MUST BE SPECIFIED:
C          E, L, AND H USE A DIFFERENT VALUE FOR X > 0.1 AND X < 0.1
      DATA TMIN/ 8.1905D0, 7.4474D0/
      DATA TMINN/ 8.06604D0, 7.4474D0/
C
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      TMN = TMIN(NSET)
      IF (T - TMN)  60,60,61
   60 XFB = 0.0D0
      RETURN
C
   61 XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMN))/(TMAX - TMN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.9010E-02,-.1401E-01,0.7150E-02,-.4130E-02,0.1260E-02,-.1040E-02,
     *0.6280E-02,-.9320E-02,0.4780E-02,-.2890E-02,0.9100E-03,-.8200E-03,
     *-.2930E-02,0.4090E-02,-.1890E-02,0.7600E-03,-.2300E-03,0.1400E-03,
     *0.3900E-03,-.1200E-02,0.4400E-03,-.2500E-03,0.2000E-04,-.2000E-04,
     *0.2600E-03,0.1400E-03,-.8000E-04,0.1000E-03,0.1000E-04,0.1000E-04,
     *-.2600E-03,0.3200E-03,0.1000E-04,-.1000E-04,0.1000E-04,-.1000E-04,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.8980E-02,-.1459E-01,0.7510E-02,-.4410E-02,0.1310E-02,-.1070E-02,
     *0.5970E-02,-.9440E-02,0.4800E-02,-.3020E-02,0.9100E-03,-.8500E-03,
     *-.3050E-02,0.4440E-02,-.2100E-02,0.8500E-03,-.2400E-03,0.1400E-03,
     *0.5300E-03,-.1300E-02,0.5600E-03,-.2700E-03,0.3000E-04,-.2000E-04,
     *0.2000E-03,0.1400E-03,-.1100E-03,0.1000E-03,0.0000E+00,0.0000E+00,
     *-.2600E-03,0.3200E-03,0.0000E+00,-.3000E-04,0.1000E-04,-.1000E-04/
C
      XFB = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFB = XFB + C(I,J,NSET)* TX(I) *TT(J)
      XFB = (1.0D0 - X)**7 *  XFB
      RETURN
C
   12 TMN = TMINN(NSET)
      IF (T - TMN)  62,62,63
   62 XFB = 0.0D0
      RETURN
C
   63 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMN))/(TMAX - TMN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     *0.8029    ,-1.075    ,0.3792    ,-.7843E-01,0.1007E-01,-.1090E-02,
     *0.7903    ,-1.099    ,0.4153    ,-.9301E-01,0.1317E-01,-.1410E-02,
     *-.1704E-01,-.1130E-01,0.2882E-01,-.1341E-01,0.3040E-02,-.3600E-03,
     *-.7200E-03,0.7230E-02,-.5160E-02,0.1080E-02,-.5000E-04,-.4000E-04,
     *0.3050E-02,-.4610E-02,0.1660E-02,-.1300E-03,-.1000E-04,0.1000E-04,
     *-.4360E-02,0.5230E-02,-.1610E-02,0.2000E-03,-.2000E-04,0.0000E+00,
C
C   D FOR NSET = 1 IS ABOVE; D FOR NSET = 2 IS BELOW
C
     *0.8672    ,-1.174    ,0.4265    ,-.9252E-01,0.1244E-01,-.1460E-02,
     *0.8500    ,-1.194    ,0.4630    ,-.1083    ,0.1614E-01,-.1830E-02,
     *-.2241E-01,-.5630E-02,0.2815E-01,-.1425E-01,0.3520E-02,-.4300E-03,
     *-.7300E-03,0.8030E-02,-.5780E-02,0.1380E-02,-.1300E-03,-.4000E-04,
     *0.3460E-02,-.5380E-02,0.1960E-02,-.2100E-03,0.1000E-04,0.1000E-04,
     *-.4850E-02,0.5950E-02,-.1890E-02,0.2600E-03,-.3000E-04,0.0000E+00/
C
      XFB = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFB = XFB + D(I,J,NSET)* TY(I) *TT(J)
      XFB = (1.0D0 - X)**7 *  XFB
      RETURN
      END
C                       ****************************
C
      FUNCTION XFC(X,Q)
C                                                   -=-=- xfc
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2)
      DIMENSION AAR(3,2), AR(3,2), BR(3,2),
     *                 ALPR(3,2), BETR(3,2), GAMR(3,2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.20D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (NSET.EQ.3)  LAMBDA = 0.20D0
      IF (NSET.EQ.4)  LAMBDA = 0.40D0
      IF (LAMBDA) 110,110,102
  102 IF (NSET-3) 104, 105, 105
C
C                   DUKE AND OWENS PARAMETRIZATION
C
  105 NCC = NSET - 2
      DATA AAR/ 0.0D0, 0.13479D0, -0.074693D0,
     *                    0.0D0, 0.067368D0, -0.030574D0/
      DATA AR/ -0.0355D0, -0.22237D0, -0.057685D0,
     *                 -0.11989D0, -0.23293D0, -0.023273D0/
      DATA BR/ 6.3494D0, 3.2649D0, -0.90945D0,
     *                 3.5087D0, 3.6554D0, -0.45313D0/
      DATA ALPR/ 0.0D0,-3.0331D0, 1.5042D0,
     *                       0.0D0, -0.47369D0, 0.35793D0/
      DATA BETR/ 0.0D0, 17.431D0, -11.255D0,
     *                       0.0D0,  9.5041D0, -5.4303D0/
      DATA GAMR/ 0.0D0,-17.861D0, 15.571D0,
     *                       0.0D0, -16.563D0,  15.524D0/
      S = LOG( LOG(Q/LAMBDA) / LOG(2.0D0/LAMBDA))
      SS = S**2
      AA = AAR(1,NCC) + AAR(2,NCC)*S + AAR(3,NCC)* SS
      A = AR(1,NCC) + AR(2,NCC)*S + AR(3,NCC) * SS
      B = BR(1,NCC) + BR(2,NCC)*S + BR(3,NCC) * SS
      ALP = ALPR(1,NCC) + ALPR(2,NCC)*S + ALPR(3,NCC)*SS
      BET = BETR(1,NCC) + BETR(2,NCC)*S + BETR(3,NCC)*SS
      GAM = GAMR(1,NCC) + GAMR(2,NCC) * S + GAMR(3,NCC) * SS
      PREF = 1.0D0 + ALP*X + BET * X**2 + GAM* X**3
      XFC = AA * X**A * (1.0D0 - X)**B *PREF
      RETURN
  110 CONTINUE
C  110 PRINT 111
  111 FORMAT('     TROUBLE:  NSET HAS NOT BEEN SET TO 1, 2, 3, OR 4'/
     * '                 I ASSUME NSET = 2')
      NSET = 2
      LAMBDA = 0.29D0
C
C        EHLQ PARAMETRIZATION
C
  104 TMIN = 2.0D0 *LOG(QMIN/LAMBDA)
      TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.9270E-02,-.1817E-01,0.9590E-02,-.6390E-02,0.1690E-02,-.1540E-02,
     *0.5710E-02,-.1188E-01,0.6090E-02,-.4650E-02,0.1240E-02,-.1310E-02,
     *-.3960E-02,0.7100E-02,-.3590E-02,0.1840E-02,-.3900E-03,0.3400E-03,
     *0.1120E-02,-.1960E-02,0.1120E-02,-.4800E-03,0.1000E-03,-.4000E-04,
     *0.4000E-04,-.3000E-04,-.1800E-03,0.9000E-04,-.5000E-04,-.2000E-04,
     *-.4200E-03,0.7300E-03,-.1600E-03,0.5000E-04,0.5000E-04,0.5000E-04,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.9980E-02,-.1945E-01,0.1055E-01,-.6870E-02,0.1860E-02,-.1560E-02,
     *0.5700E-02,-.1203E-01,0.6250E-02,-.4860E-02,0.1310E-02,-.1370E-02,
     *-.4490E-02,0.7990E-02,-.4170E-02,0.2050E-02,-.4400E-03,0.3300E-03,
     *0.1470E-02,-.2480E-02,0.1460E-02,-.5700E-03,0.1200E-03,-.1000E-04,
     *-.9000E-04,0.1500E-03,-.3200E-03,0.1200E-03,-.6000E-04,-.4000E-04,
     *-.4200E-03,0.7600E-03,-.1400E-03,0.4000E-04,0.7000E-04,0.5000E-04/
C
      XFC = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFC = XFC + C(I,J,NSET)* TX(I) *TT(J)
      XFC = (1.0D0 - X)**7 *  XFC
      IF (NSET.EQ.1) XFC = XFC * (1.0D0 -X)
      RETURN
C
   12 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     *0.8098    ,-1.042    ,0.3398    ,-.6824E-01,0.8760E-02,-.9000E-03,
     *0.8961    ,-1.217    ,0.4339    ,-.9287E-01,0.1304E-01,-.1290E-02,
     *0.3058E-01,-.1040    ,0.7604E-01,-.2415E-01,0.4600E-02,-.5000E-03,
     *-.2451E-01,0.4432E-01,-.1651E-01,0.1430E-02,0.1200E-03,-.1000E-03,
     *0.1122E-01,-.1457E-01,0.2680E-02,0.5800E-03,-.1200E-03,0.3000E-04,
     *-.7730E-02,0.7330E-02,-.7600E-03,-.2400E-03,0.1000E-04,0.0000E+00,
C
C  D FOR NSET = 1 IS ABOVE, D FOR NSET = 2 IS BELOW
C
     *0.8698    ,-1.131    ,0.3836    ,-.8111E-01,0.1048E-01,-.1300E-02,
     *0.9626    ,-1.321    ,0.4854    ,-.1091    ,0.1583E-01,-.1700E-02,
     *0.3057E-01,-.1088    ,0.8022E-01,-.2676E-01,0.5590E-02,-.5600E-03,
     *-.2845E-01,0.5164E-01,-.1918E-01,0.2210E-02,-.4000E-04,-.1500E-03,
     *0.1311E-01,-.1751E-01,0.3310E-02,0.5100E-03,-.1200E-03,0.5000E-04,
     *-.8590E-02,0.8380E-02,-.9200E-03,-.2600E-03,0.1000E-04,-.1000E-04/
C
      XFC = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFC = XFC + D(I,J,NSET)* TY(I) *TT(J)
      XFC = (1.0D0 - X)**7 *  XFC
      IF (NSET.EQ.1) XFC = XFC * (1.0D0 -X)
      RETURN
      END

      FUNCTION XFD(X,Q)
C                                                   -=-=- xfd
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2)
      DIMENSION ETA3R(3,2), ETA4R(3,2), GAMDR(3,2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.20D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (NSET.EQ.3)  LAMBDA = 0.20D0
      IF (NSET.EQ.4)  LAMBDA = 0.40D0
      IF (LAMBDA) 110,110,112
  112 IF (NSET-3) 104, 105, 105
C
C                   DUKE AND OWENS PARAMETRIZATION
C
  105 NCC = NSET -2
      DATA ETA3R/ 0.763D0, -0.23696D0, 0.025836D0,
     *                         0.7608D0,-0.2317D0, 0.023232D0/
      DATA ETA4R/ 4.00D0, 0.62664D0, -0.019163D0,
     *                         3.83D0, 0.62746D0, -0.019155D0/
      DATA GAMDR/  0.0D0, -0.42068D0, 0.032809D0,
     *                            0.0D0, -0.41843D0, 0.035972D0/
      S = LOG( LOG(Q/LAMBDA) / LOG(2.0D0/LAMBDA))
      SS = S**2
      ETA3 = ETA3R(1,NCC) + ETA3R(2,NCC)*S + ETA3R(3,NCC) * SS
      ETA4 = ETA4R(1,NCC) + ETA4R(2,NCC) * S + ETA4R(3,NCC) * SS
      GAMD = GAMDR(1,NCC) + GAMDR(2,NCC) * S + GAMDR(3,NCC) * SS
      BTA = GAMA(ETA3)*GAMA(ETA4+1.0D0)/GAMA(ETA3+ETA4 + 1.0D0)
      XND = 1.0D0/(BTA * (1.0D0 + GAMD*ETA3/(ETA3 + ETA4 + 1.0D0)))
      XFDV = XND*(X**ETA3) * (1.0D0 - X)**ETA4 * (1.0D0+ GAMD * X)
      XFD =   XFDV + XFSEA(X,Q)
      RETURN
  110 CONTINUE
C  110 PRINT 111
  111 FORMAT('     TROUBLE:  NSET HAS NOT BEEN SET TO 1, 2, 3, OR 4'/
     * '                 I ASSUME NSET = 2')
      NSET = 2
      LAMBDA = 0.29D0
C
C        EHLQ PARAMETRIZATION
C
  104 TMIN = 2.0D0 *LOG(QMIN/LAMBDA)
      TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.38130   ,-.8090E-01,-.16336   ,-.2185E-01,-.8430E-02,-.6200E-03,
     *-.29475   ,-.14348   ,0.16650   ,0.6638E-01,0.1473E-01,0.4080E-02,
     *0.12518   ,0.10422   ,-.4722E-01,-.3683E-01,-.1038E-01,-.2860E-02,
     *-.5478E-01,-.5678E-01,0.8900E-02,0.1484E-01,0.5340E-02,0.1520E-02,
     *0.2220E-01,0.2567E-01,-.3000E-04,-.4970E-02,-.2160E-02,-.6500E-03,
     *-.9530E-02,-.1204E-01,-.1510E-02,0.1510E-02,0.8300E-03,0.2700E-03,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.3578    ,-.8622E-01,-.1480    ,-.1840E-01,-.7820E-02,-.4500E-03,
     *-.2925    ,-.1304    ,0.1696    ,0.6243E-01,0.1353E-01,0.3750E-02,
     *0.1318    ,0.1041    ,-.5486E-01,-.3872E-01,-.1038E-01,-.2850E-02,
     *-.6162E-01,-.6143E-01,0.1303E-01,0.1740E-01,0.5940E-02,0.1670E-02,
     *0.2643E-01,0.2957E-01,-.1490E-02,-.6450E-02,-.2630E-02,-.7700E-03,
     *-.1218E-01,-.1497E-01,-.1260E-02,0.2240E-02,0.1120E-02,0.3500E-03/
C
      XFDVAL = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFDVAL = XFDVAL + C(I,J,NSET)* TX(I) *TT(J)
      XFDVAL = (1.0D0 - X)**4 *  XFDVAL
      XFD = XFDVAL + XFSEA(X,Q)
      RETURN
C
   12 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     *0.12613   ,0.13542   ,0.3958E-01,0.8240E-02,0.1660E-02,0.4500E-03,
     *0.3890E-02,-.1159E-01,-.1625E-01,-.9610E-02,-.3710E-02,-.1260E-02,
     *-.1910E-02,-.5600E-03,0.1590E-02,0.1590E-02,0.8400E-03,0.3900E-03,
     *0.6400E-03,0.4900E-03,-.1500E-03,-.2900E-03,-.1800E-03,-.1000E-03,
     *-.2000E-03,-.1900E-03,0.0000E+00,0.6000E-04,0.4000E-04,0.3000E-04,
     *0.7000E-04,0.8000E-04,0.2000E-04,-.1000E-04,-.1000E-04,-.1000E-04,
C
C   D FOR NSET = 1 IS ABOVE; D FOR NSET = 2 IS BELOW
C
     *0.1263    ,0.1334    ,0.3732E-01,0.7070E-02,0.1260E-02,0.3400E-03,
     *0.3660E-02,-.1357E-01,-.1795E-01,-.1031E-01,-.3880E-02,-.1280E-02,
     *-.2100E-02,-.3600E-03,0.2050E-02,0.1920E-02,0.9800E-03,0.4400E-03,
     *0.7700E-03,0.5400E-03,-.2400E-03,-.3900E-03,-.2400E-03,-.1300E-03,
     *-.2600E-03,-.2300E-03,0.2000E-04,0.9000E-04,0.6000E-04,0.4000E-04,
     *0.9000E-04,0.1000E-03,0.2000E-04,-.2000E-04,-.2000E-04,-.1000E-04/
C
      XFDVAL = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFDVAL = XFDVAL + D(I,J,NSET)* TY(I) *TT(J)
      XFDVAL = (1.0D0 - X)**4 *  XFDVAL
      XFD = XFDVAL + XFSEA(X,Q)
      RETURN
      END
C                       ****************************
C
      FUNCTION XFG(X,Q)
C                                                   -=-=- xfg
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2)
      DIMENSION AAR(3,2), AR(3,2), BR(3,2),
     *                 ALPR(3,2), BETR(3,2), GAMR(3,2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.20D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (NSET.EQ.3)  LAMBDA = 0.20D0
      IF (NSET.EQ.4)  LAMBDA = 0.40D0
      IF (LAMBDA) 110,110,102
  102 IF (NSET-3) 104, 105, 105
C
C                   DUKE AND OWENS PARAMETRIZATION
C
  105 NCC = NSET - 2
      DATA AAR/ 1.564D0, -1.7112D0, 0.63751D0,
     *                    0.8789D0, -0.97093D0, 0.43388D0/
      DATA AR/ 0.0D0, -0.94892D0, 0.32505D0,
     *                 0.0D0, -1.1612D0,  0.4759D0/
      DATA BR/ 6.00D0, 1.4345D0, -0.10485D1,
     *                 4.0D0, 1.2271D0, -0.25369D0/
      DATA ALPR/ 9.0D0,-7.1858D0,  0.25494D0,
     *                       9.0D0,-5.6354D0, -0.81747D0/
      DATA BETR/ 0.0D0, -16.457D0, 10.947D0,
     *                       0.0D0, -7.5438D0,  5.5034D0/
      DATA GAMR/ 0.0D0, 15.261D0, -10.085D0,
     *                       0.0D0,-0.59649D0,  0.12611D0/
      S = LOG( LOG(Q/LAMBDA) / LOG(2.0D0/LAMBDA))
      SS = S**2
      AA = AAR(1,NCC) + AAR(2,NCC)*S + AAR(3,NCC)* SS
      A = AR(1,NCC) + AR(2,NCC)*S + AR(3,NCC) * SS
      B = BR(1,NCC) + BR(2,NCC)*S + BR(3,NCC) * SS
      ALP = ALPR(1,NCC) + ALPR(2,NCC)*S + ALPR(3,NCC)*SS
      BET = BETR(1,NCC) + BETR(2,NCC)*S + BETR(3,NCC)*SS
      GAM = GAMR(1,NCC) + GAMR(2,NCC) * S + GAMR(3,NCC) * SS
      PREF = 1.0D0 + ALP*X + BET * X**2 + GAM* X**3
      XFG = AA * X**A * (1.0D0 - X)**B *PREF
      RETURN
  110 CONTINUE
C  110 PRINT 111
  111 FORMAT('     TROUBLE:  NSET HAS NOT BEEN SET TO 1, 2, 3, OR 4'/
     * '                 I ASSUME NSET = 2')
      NSET = 2
      LAMBDA = 0.29D0
C
C        EHLQ PARAMETRIZATION
C
  104 TMIN = 2.0D0 *LOG(QMIN/LAMBDA)
      TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.94819   ,-.95779   ,0.10085   ,-.10510   ,0.3456E-01,-.3054E-01,
     *-.96265   ,0.53790   ,0.33684   ,-.9525E-01,0.1488E-01,-.2051E-01,
     *0.43004   ,-.8306E-01,-.33719   ,0.4902E-01,-.9160E-02,0.1041E-01,
     *-.19249   ,-.1790E-01,0.21830   ,0.7490E-02,0.4140E-02,-.1860E-02,
     *0.8183E-01,0.1926E-01,-.10718   ,-.1944E-01,-.2770E-02,-.5200E-03,
     *-.3884E-01,-.1234E-01,0.5410E-01,0.1879E-01,0.3350E-02,0.1040E-02,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     * 2.367    ,0.4453    ,0.3660    ,0.9467E-01,0.1341    ,0.1661E-01,
     *-3.170    ,-1.795    ,0.3313E-01,-.2874    ,-.9827E-01,-.7119E-01,
     * 1.823    , 1.457    ,-.2465    ,0.3739E-01,0.6090E-02,0.1814E-01,
     *-1.033    ,-.9827    ,0.2136    ,0.1169    ,0.5001E-01,0.1684E-01,
     *0.5133    ,0.5259    ,-.1173    ,-.1139    ,-.4988E-01,-.2021E-01,
     *-.2881    ,-.3145    ,0.5667E-01,0.9161E-01,0.4568E-01,0.1951E-01/
C
      XFG = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFG = XFG + C(I,J,NSET)* TX(I) *TT(J)
      XFG = (1.0D0 - X)**5 *  XFG
      IF (NSET.EQ.2)  XFG = XFG * (1.0D0 - X)
      RETURN
C
   12 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     * 29.47734 ,-39.02468 , 14.6357  ,-3.33516  ,0.50538   ,-.5915E-01,
     * 25.5896  ,-39.54527 , 16.6142  ,-4.29861  ,0.69036   ,-.8243E-01,
     *-1.66291  , 1.17624  , 1.11844  ,-.70986   ,0.19481   ,-.2404E-01,
     *-.21679   ,0.81705   ,-.71688   ,0.18507   ,-.1924E-01,-.3250E-02,
     *0.20880   ,-.43547   ,0.22391   ,-.2446E-01,-.3620E-02,0.1910E-02,
     *-.9097E-01,0.16009   ,-.5681E-01,-.2500E-02,0.2580E-02,-.4700E-03,
C
C   D FOR NSET = 1 IS ABOVE; D FOR NSET = 2 IS BELOW
C
     * 30.36    ,-40.62    , 15.78    ,-3.699    ,0.6020    ,-.7031E-01,
     * 27.00    ,-41.67    , 17.70    ,-4.804    ,0.7862    ,-.1060    ,
     *-1.909    , 1.357    , 1.127    ,-.7181    ,0.2232    ,-.2481E-01,
     *-.2488    ,0.9781    ,-.8127    ,0.2094    ,-.2997E-01,-.4710E-02,
     *0.2506    ,-.5427    ,0.2672    ,-.3103E-01,-.1800E-02,0.2870E-02,
     *-.1128    ,0.2087    ,-.6972E-01,-.2480E-02,0.2630E-02,-.8400E-03/
C
      XFG = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFG = XFG + D(I,J,NSET)* TY(I) *TT(J)
      XFG = (1.0D0 - X)**5 *  XFG
      IF (NSET.EQ.2)  XFG = XFG *(1.0D0 -X)
      RETURN
      END
C                       ****************************
C
      FUNCTION XFSEA(X,Q)
C                                                   -=-=- xfsea
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2)
      DIMENSION AAR(3,2), AR(3,2), BR(3,2),
     *                 ALPR(3,2), BETR(3,2), GAMR(3,2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.20D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (NSET.EQ.3)  LAMBDA = 0.20D0
      IF (NSET.EQ.4)  LAMBDA = 0.40D0
      IF (LAMBDA) 110,110,102
  102 IF (NSET-3) 104, 105, 105
C
C                   DUKE AND OWENS PARAMETRIZATION
C
  105 NCC = NSET - 2
      DATA AAR/ 1.265D0, -1.1323D0, 0.29268D0,
     *                    1.6714D0, -1.9168D0, 0.58175D0/
      DATA AR/ 0.0D0, -0.37162D0, -0.028977D0,
     *                 0.0D0, -0.27307D0, -0.16392D0/
      DATA BR/ 8.05D0, 1.5877D0, -0.15291D0,
     *                 9.145D0, 0.53045D0, -0.76271D0/
      DATA ALPR/ 0.0D0, 6.3059D0, -0.27342D0,
     *                       0.0D0, 15.665D0, -2.8341D0/
      DATA BETR/ 0.0D0, -10.543D0, -3.1674D0,
     *                       0.0D0, -100.63D0,  44.658D0/
      DATA GAMR/ 0.0D0, 14.698D0, 9.798D0,
     *                       0.0D0, 223.24D0, -116.76D0/
      S = LOG( LOG(Q/LAMBDA) / LOG(2.0D0/LAMBDA))
      SS = S**2
      AA = AAR(1,NCC) + AAR(2,NCC)*S + AAR(3,NCC)* SS
      A = AR(1,NCC) + AR(2,NCC)*S + AR(3,NCC) * SS
      B = BR(1,NCC) + BR(2,NCC)*S + BR(3,NCC) * SS
      ALP = ALPR(1,NCC) + ALPR(2,NCC)*S + ALPR(3,NCC)*SS
      BET = BETR(1,NCC) + BETR(2,NCC)*S + BETR(3,NCC)*SS
      GAM = GAMR(1,NCC) + GAMR(2,NCC) * S + GAMR(3,NCC) * SS
      PREF = 1.0D0 + ALP*X + BET * X**2 + GAM* X**3
      XFSEA = AA * X**A * (1.0D0 - X)**B *PREF/ (6.0D0)
      RETURN
  110 CONTINUE
C  110 PRINT 111
  111 FORMAT('     TROUBLE:  NSET HAS NOT BEEN SET TO 1, 2, 3, OR 4'/
     * '                 I ASSUME NSET = 2')
      NSET = 2
      LAMBDA = 0.29D0
C
C        EHLQ PARAMETRIZATION
C
  104 TMIN = 2.0D0 *LOG(QMIN/LAMBDA)
      TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.6870E-01,-.6861E-01,0.2973E-01,-.5400E-02,0.3780E-02,-.9700E-03,
     *-.1802E-01,0.1400E-03,0.6490E-02,-.8540E-02,0.1220E-02,-.1750E-02,
     *-.4650E-02,0.1480E-02,-.5930E-02,0.6000E-03,-.1030E-02,-.8000E-04,
     *0.6440E-02,0.2570E-02,0.2830E-02,0.1150E-02,0.7100E-03,0.3300E-03,
     *-.3930E-02,-.2540E-02,-.1160E-02,-.7700E-03,-.3600E-03,-.1900E-03,
     *0.2340E-02,0.1930E-02,0.5300E-03,0.3700E-03,0.1600E-03,0.9000E-04,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.1008    ,-.7100E-01,0.1973E-01,-.5710E-02,0.2930E-02,-.9900E-03,
     *-.5271E-01,-.1823E-01,0.1792E-01,-.6580E-02,0.1750E-02,-.1550E-02,
     *0.1220E-01,0.1763E-01,-.8690E-02,-.8800E-03,-.1160E-02,-.2100E-03,
     *-.1190E-02,-.7180E-02,0.2360E-02,0.1890E-02,0.7700E-03,0.4100E-03,
     *-.9100E-03,0.2040E-02,-.3100E-03,-.1050E-02,-.4000E-03,-.2400E-03,
     *0.1190E-02,-.1700E-03,-.2000E-03,0.4200E-03,0.1700E-03,0.1000E-03/
C
      XFSEA = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFSEA = XFSEA + C(I,J,NSET)* TX(I) *TT(J)
      XFSEA = (1.0D0 - X)**7 *  XFSEA
      IF (NSET.EQ.1)  XFSEA = XFSEA * (1.0D0 - X)
      RETURN
C
   12 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     * 1.01386  ,-1.10585  ,0.33739   ,-.7444E-01,0.8850E-02,-.8700E-03,
     *0.92334   ,-1.28541  ,0.44755   ,-.9786E-01,0.1419E-01,-.1120E-02,
     *0.4888E-01,-.12708   ,0.8606E-01,-.2608E-01,0.4780E-02,-.6000E-03,
     *-.2691E-01,0.4887E-01,-.1771E-01,0.1620E-02,0.2500E-03,-.6000E-04,
     *0.7040E-02,-.1113E-01,0.1590E-02,0.7000E-03,-.2000E-03,0.0000E+00,
     *-.1710E-02,0.2290E-02,0.3800E-03,-.3500E-03,0.4000E-04,0.1000E-04,
C
C  D FOR NSET = 1 IS ABOVE, D FOR NSET = 2 IS BELOW
C
     * 1.081    ,-1.189    ,0.3868    ,-.8617E-01,0.1115E-01,-.1180E-02,
     *0.9917    ,-1.396    ,0.4998    ,-.1159    ,0.1674E-01,-.1720E-02,
     *0.5099E-01,-.1338    ,0.9173E-01,-.2885E-01,0.5890E-02,-.6500E-03,
     *-.3178E-01,0.5703E-01,-.2070E-01,0.2440E-02,0.1100E-03,-.9000E-04,
     *0.8970E-02,-.1392E-01,0.2050E-02,0.6500E-03,-.2300E-03,0.2000E-04,
     *-.2340E-02,0.3010E-02,0.5000E-03,-.3900E-03,0.6000E-04,0.1000E-04/
C
      XFSEA = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFSEA = XFSEA + D(I,J,NSET)* TY(I) *TT(J)
      XFSEA = (1.0D0 - X)**7 *  XFSEA
      IF (NSET.EQ.1)  XFSEA = XFSEA * (1.0D0 - X)
      RETURN
      END
C                       ****************************
C
      FUNCTION XFS(X,Q)
C                                                   -=-=- xfs
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2)
      DIMENSION AAR(3,2), AR(3,2), BR(3,2),
     *                 ALPR(3,2), BETR(3,2), GAMR(3,2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.20D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (NSET.EQ.3)  LAMBDA = 0.20D0
      IF (NSET.EQ.4)  LAMBDA = 0.40D0
      IF (LAMBDA) 110,110,102
  102 IF (NSET-3) 104, 105, 105
C
C                   DUKE AND OWENS PARAMETRIZATION
C
  105 NCC = NSET - 2
      DATA AAR/ 1.265D0, -1.1323D0, 0.29268D0,
     *                    1.6714D0, -1.9168D0, 0.58175D0/
      DATA AR/ 0.0D0, -0.37162D0, -0.028977D0,
     *                 0.0D0, -0.27307D0, -0.16392D0/
      DATA BR/ 8.05D0, 1.5877D0, -0.15291D0,
     *                 9.145D0, 0.53045D0, -0.76271D0/
      DATA ALPR/ 0.0D0, 6.3059D0, -0.27342D0,
     *                       0.0D0, 15.665D0, -2.8341D0/
      DATA BETR/ 0.0D0, -10.543D0, -3.1674D0,
     *                       0.0D0, -100.63D0,  44.658D0/
      DATA GAMR/ 0.0D0, 14.698D0, 9.798D0,
     *                       0.0D0, 223.24D0, -116.76D0/
      S = LOG( LOG(Q/LAMBDA) / LOG(2.0D0/LAMBDA))
      SS = S**2
      AA = AAR(1,NCC) + AAR(2,NCC)*S + AAR(3,NCC)* SS
      A = AR(1,NCC) + AR(2,NCC)*S + AR(3,NCC) * SS
      B = BR(1,NCC) + BR(2,NCC)*S + BR(3,NCC) * SS
      ALP = ALPR(1,NCC) + ALPR(2,NCC)*S + ALPR(3,NCC)*SS
      BET = BETR(1,NCC) + BETR(2,NCC)*S + BETR(3,NCC)*SS
      GAM = GAMR(1,NCC) + GAMR(2,NCC) * S + GAMR(3,NCC) * SS
      PREF = 1.0D0 + ALP*X + BET * X**2 + GAM* X**3
      XFS = AA * X**A * (1.0D0 - X)**B *PREF/ (6.0D0)
      RETURN
  110 CONTINUE
C  110 PRINT 111
  111 FORMAT('     TROUBLE:  NSET HAS NOT BEEN SET TO 1, 2, 3, OR 4'/
     * '                 I ASSUME NSET = 2')
      NSET = 2
      LAMBDA = 0.29D0
C
C        EHLQ PARAMETRIZATION
C
  104 TMIN = 2.0D0 *LOG(QMIN/LAMBDA)
      TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.4968E-01,-.4173E-01,0.2102E-01,-.3270E-02,0.3240E-02,-.6700E-03,
     *-.6150E-02,-.1294E-01,0.6740E-02,-.6890E-02,0.9000E-03,-.1510E-02,
     *-.8580E-02,0.5050E-02,-.4900E-02,-.1600E-03,-.9400E-03,-.1500E-03,
     *0.7840E-02,0.1510E-02,0.2220E-02,0.1400E-02,0.7000E-03,0.3500E-03,
     *-.4410E-02,-.2220E-02,-.8900E-03,-.8500E-03,-.3600E-03,-.2000E-03,
     *0.2520E-02,0.1840E-02,0.4100E-03,0.3900E-03,0.1600E-03,0.9000E-04,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.6478E-01,-.4537E-01,0.1643E-01,-.3490E-02,0.2710E-02,-.6700E-03,
     *-.2223E-01,-.2126E-01,0.1247E-01,-.6290E-02,0.1120E-02,-.1440E-02,
     *-.1340E-02,0.1362E-01,-.6130E-02,-.7900E-03,-.9000E-03,-.2000E-03,
     *0.5080E-02,-.3610E-02,0.1700E-02,0.1830E-02,0.6800E-03,0.4000E-03,
     *-.3580E-02,0.6000E-04,-.2600E-03,-.1050E-02,-.3800E-03,-.2300E-03,
     *0.2420E-02,0.9300E-03,-.1000E-03,0.4500E-03,0.1700E-03,0.1100E-03/
C
      XFS = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFS = XFS + C(I,J,NSET)* TX(I) *TT(J)
      XFS = (1.0D0 - X)**7 *  XFS
      IF (NSET.EQ.1) XFS = XFS * (1.0D0 - X)
      RETURN
C
   12 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     *0.92351   ,-1.08483  ,0.34642   ,-.7210E-01,0.9140E-02,-.9100E-03,
     *0.93146   ,-1.27376  ,0.45122   ,-.9775E-01,0.1380E-01,-.1310E-02,
     *0.4739E-01,-.12960   ,0.8482E-01,-.2642E-01,0.4760E-02,-.5700E-03,
     *-.2653E-01,0.4953E-01,-.1735E-01,0.1750E-02,0.2800E-03,-.6000E-04,
     *0.6940E-02,-.1132E-01,0.1480E-02,0.6500E-03,-.2100E-03,0.0000E+00,
     *-.1680E-02,0.2340E-02,0.4200E-03,-.3400E-03,0.5000E-04,0.1000E-04,
C
C   D FOR NSET = 1 IS ABOVE; D FOR NSET = 2 IS BELOW
C
     *0.9868    ,-1.171    ,0.3940    ,-.8459E-01,0.1124E-01,-.1250E-02,
     * 1.001    ,-1.383    ,0.5044    ,-.1152    ,0.1658E-01,-.1830E-02,
     *0.4928E-01,-.1368    ,0.9021E-01,-.2935E-01,0.5800E-02,-.6600E-03,
     *-.3133E-01,0.5785E-01,-.2023E-01,0.2630E-02,0.1600E-03,-.8000E-04,
     *0.8840E-02,-.1416E-01,0.1900E-02,0.5800E-03,-.2500E-03,0.1000E-04,
     *-.2300E-02,0.3080E-02,0.5500E-03,-.3700E-03,0.7000E-04,0.1000E-04/
C
      XFS = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFS = XFS + D(I,J,NSET)* TY(I) *TT(J)
      XFS = (1.0D0 - X)**7 *  XFS
      IF (NSET.EQ.1) XFS = XFS * (1.0D0 - X)
      RETURN
      END
C
      FUNCTION XFT(X,Q)
C                                                   -=-=- xft
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2),TMIN(2),TMINN(2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.2D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (LAMBDA) 101,101,104
  101 CONTINUE
C  101 PRINT 102
  102 FORMAT('   Duke_Owens has no T_quark content; '/
     * '         I HAVE SET NSET = 2 (EHLQ SET 2) FOR XFT.')
      NSET = 2
      LAMBDA = 0.29D0
  104 TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C         IN THIS CASE, A LARGER TMIN MUST BE SPECIFIED:
C          E, L, AND H USE A DIFFERENT VALUE FOR X > 0.1 AND X < 0.1
      DATA TMIN/ 11.5528D0, 11.4283D0/
      DATA TMINN/ 10.8097D0, 10.8097D0/
C
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      TMN = TMIN(NSET)
      IF (T - TMN)  60,60,61
   60 XFT = 0.0D0
      RETURN
C
   61 XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMN))/(TMAX - TMN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.4410E-02,-.7480E-02,0.3770E-02,-.2580E-02,0.7300E-03,-.7100E-03,
     *0.3840E-02,-.6050E-02,0.3030E-02,-.2030E-02,0.5800E-03,-.5900E-03,
     *-.8800E-03,0.1660E-02,-.7500E-03,0.4700E-03,-.1000E-03,0.1000E-03,
     *-.8000E-04,-.1500E-03,0.1200E-03,-.9000E-04,0.3000E-04,0.0000E+00,
     *0.1300E-03,-.2200E-03,-.2000E-04,-.2000E-04,-.2000E-04,-.2000E-04,
     *-.7000E-04,0.1900E-03,-.4000E-04,0.2000E-04,0.0000E+00,0.0000E+00,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.4260E-02,-.7530E-02,0.3830E-02,-.2680E-02,0.7600E-03,-.7300E-03,
     *0.3640E-02,-.6050E-02,0.3030E-02,-.2090E-02,0.5900E-03,-.6000E-03,
     *-.9200E-03,0.1710E-02,-.8200E-03,0.5000E-03,-.1200E-03,0.1000E-03,
     *-.5000E-04,-.1600E-03,0.1300E-03,-.9000E-04,0.3000E-04,0.0000E+00,
     *0.1300E-03,-.2100E-03,-.1000E-04,-.2000E-04,-.2000E-04,-.1000E-04,
     *-.8000E-04,0.1800E-03,-.5000E-04,0.2000E-04,0.0000E+00,0.0000E+00/
C
      XFT = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFT = XFT + C(I,J,NSET)* TX(I) *TT(J)
      XFT = (1.0D0 - X)**7 *  XFT
      RETURN
C
   12 TMN = TMINN(NSET)
      IF (T - TMN)  62,62,63
   62 XFT = 0.0D0
      RETURN
C
   63 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMN))/(TMAX - TMN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     *0.6623    ,-.9248    ,0.3519    ,-.7930E-01,0.1110E-01,-.1180E-02,
     *0.6380    ,-.9062    ,0.3582    ,-.8479E-01,0.1265E-01,-.1390E-02,
     *-.2581E-01,0.2125E-01,0.4190E-02,-.4980E-02,0.1490E-02,-.2100E-03,
     *0.7100E-03,0.5300E-03,-.1270E-02,0.3900E-03,-.5000E-04,-.1000E-04,
     *0.3850E-02,-.5060E-02,0.1860E-02,-.3500E-03,0.4000E-04,0.0000E+00,
     *-.3530E-02,0.4460E-02,-.1500E-02,0.2700E-03,-.3000E-04,0.0000E+00,
C
C  D FOR NSET = 1 IS ABOVE, D FOR NSET = 2 IS BELOW
C
     *0.7146    ,-1.007    ,0.3932    ,-.9246E-01,0.1366E-01,-.1540E-02,
     *0.6856    ,-.9828    ,0.3977    ,-.9795E-01,0.1540E-01,-.1790E-02,
     *-.3053E-01,0.2758E-01,0.2150E-02,-.4880E-02,0.1640E-02,-.2500E-03,
     *0.9200E-03,0.4200E-03,-.1340E-02,0.4600E-03,-.8000E-04,-.1000E-04,
     *0.4230E-02,-.5660E-02,0.2140E-02,-.4300E-03,0.6000E-04,0.0000E+00,
     *-.3890E-02,0.5000E-02,-.1740E-02,0.3300E-03,-.4000E-04,0.0000E+00/
C
      XFT = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFT = XFT + D(I,J,NSET)* TY(I) *TT(J)
      XFT = (1.0D0 - X)**7 *  XFT
      RETURN
      END
C                       ****************************
C
      FUNCTION XFU(X,Q)
C                                                   -=-=- xfu
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION LAMBDA
      DIMENSION TX(6),TT(6),TY(6),C(6,6,2),D(6,6,2)
      DIMENSION ETA1R(3,2), ETA2R(3,2), GAMUDR(3,2), ETA3R(3,2),
     *         ETA4R(3,2), GAMDR(3,2)
      COMMON/NSET/ NSET
      IF(NSET.LE.2) THEN
         QMIN = DSQRT(5.0D0)
      Else
         QMIN = DSQRT(4.0D0)
      EndIf
      QMAX = 1.0D4
C  CHECK WHETHER VALUES OF X AND Q ARE IN THE RIGHT RANGE
    1 FORMAT('    TROUBLE: Q IS TOO SMALL (LESS THAN SQRT(5) GEV)')
    2 FORMAT('    TROUBLE: Q IS TOO LARGE (GREATER THAN 10**4 GEV)')
    3 FORMAT('    TROUBLE: X IS TOO SMALL (LESS THAN 10**(-4))')
    4 FORMAT('    TROUBLE: X IS TOO LARGE (GREATER THAN 1)')
C  DEFINE VARIOUS PARAMETERS
      LAMBDA = 0.0D0
      IF (NSET.EQ.1)  LAMBDA = 0.20D0
      IF (NSET.EQ.2)  LAMBDA = 0.29D0
      IF (NSET.EQ.3)  LAMBDA = 0.20D0
      IF (NSET.EQ.4)  LAMBDA = 0.40D0
      IF (LAMBDA) 110,110,102
  102 IF (NSET-3) 104, 105, 105
C
C                   DUKE AND OWENS PARAMETRIZATION
C
  105 NCC = NSET - 2
      DATA ETA1R/ 0.419D0, 0.004383D0, -0.007412D0,
     *                         0.3743D0, 0.013946D0,  0.00031695D0/
      DATA ETA2R/ 3.46D0, 0.72432D0, -0.065998D0,
     *                         3.329D0, 0.75343D0, -0.076125D0/
      DATA GAMUDR/ 4.40D0, -4.8644D0, 1.3274D0,
     *                              6.032D0, -6.2153D0, 1.5561D0/
      DATA ETA3R/ 0.763D0, -0.23696D0, 0.025836D0,
     *                         0.7608D0,-0.2317D0, 0.023232D0/
      DATA ETA4R/ 4.00D0, 0.62664D0, -0.019163D0,
     *                         3.83D0, 0.62746D0, -0.019155D0/
      DATA GAMDR/  0.0D0, -0.42068D0, 0.032809D0,
     *                            0.0D0, -0.41843D0, 0.035972D0/
      S = LOG( LOG(Q/LAMBDA) / LOG(2.0D0/LAMBDA))
      SS = S**2
      ETA1 = ETA1R(1,NCC) + ETA1R(2,NCC)*S + ETA1R(3,NCC) * SS
      ETA2 = ETA2R(1,NCC) + ETA2R(2,NCC) * S + ETA2R(3,NCC) * SS
      GAMUD = GAMUDR(1,NCC) + GAMUDR(2,NCC) * S + GAMUDR(3,NCC) * SS
      ETA3 = ETA3R(1,NCC) + ETA3R(2,NCC) * S + ETA3R(3,NCC) * SS
      ETA4 = ETA4R(1,NCC) + ETA4R(2,NCC) * S + ETA4R(3,NCC) * SS
      GAMD = GAMDR(1,NCC) + GAMDR(2,NCC) * S + GAMDR(3,NCC) * SS
      BTA = GAMA(ETA1)*GAMA(ETA2+ 1.0D0)/GAMA(ETA1+ETA2 + 1.0D0)
      XNUD = 3.0D0/(BTA * (1.0D0 + GAMUD*ETA1/(ETA1 + ETA2 + 1.0D0)))
      XFUD = XNUD*(X**ETA1) * (1.0D0 - X)**ETA2 * (1.0D0 + GAMUD * X)
      BTA = GAMA(ETA3)*GAMA(ETA4+ 1.0D0)/GAMA(ETA3+ETA4 + 1.0D0)
      XND = 1.0D0/(BTA * (1.0D0 + GAMD*ETA3/(ETA3 + ETA4 + 1.0D0)))
      XFDV = XND*(X**ETA3) * (1.0D0 - X)**ETA4 * (1.0D0 + GAMD * X)
      XFUV = XFUD -XFDV
      XFU =   XFUV + XFSEA(X,Q)
      RETURN
  110 CONTINUE
c  110 PRINT 111
  111 FORMAT('     TROUBLE:  NSET HAS NOT BEEN SET TO 1, 2, 3, OR 4'/
     * '                 I ASSUME NSET = 2')
      NSET = 2
      LAMBDA = 0.29D0
C
C        EHLQ PARAMETRIZATION
C
  104 TMIN = 2.0D0 *LOG(QMIN/LAMBDA)
      TMAX = 2.0D0 *LOG(QMAX/LAMBDA)
      T  =  2.0D0 * LOG(Q/LAMBDA)
C
C  EHL & Q GIVE DIFFERENT PARAMETRIZATIONS FOR X >  0.1 AND X < 0.1
      IF(X.LT.0.1D0) GO TO 12
C
      XP = (2.0D0 * X - 1.1D0)/0.9D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TX,XP)
      CALL TCHEBY(TT,TP)
C
      DATA C/
     *0.76772   ,-.20874   ,-.33026   ,-.2517E-01,-.1570E-01,-.1000E-03,
     *-.53259   ,-.26612   ,0.32007   ,0.11918   ,0.2434E-01,0.7620E-02,
     *0.21618   ,0.18812   ,-.8375E-01,-.6515E-01,-.1743E-01,-.5040E-02,
     *-.9211E-01,-.9952E-01,0.1373E-01,0.2506E-01,0.8770E-02,0.2550E-02,
     *0.3670E-01,0.4409E-01,0.9600E-03,-.7960E-02,-.3420E-02,-.1050E-02,
     *-.1549E-01,-.2026E-01,-.3060E-02,0.2220E-02,0.1240E-02,0.4100E-03,
C
C  C FOR NSET = 1 IS ABOVE, C FOR NSET = 2 IS BELOW
C
     *0.7237    ,-.2189    ,-.2995    ,-.1909E-01,-.1477E-01,0.2500E-03,
     *-.5314    ,-.2425    ,0.3283    ,0.1119    ,0.2223E-01,0.7070E-02,
     *0.2289    ,0.1890    ,-.9859E-01,-.6900E-01,-.1747E-01,-.5080E-02,
     *-.1041    ,-.1084    ,0.2108E-01,0.2975E-01,0.9830E-02,0.2830E-02,
     *0.4394E-01,0.5116E-01,-.1410E-02,-.1055E-01,-.4230E-02,-.1270E-02,
     *-.1991E-01,-.2539E-01,-.2780E-02,0.3430E-02,0.1720E-02,0.5500E-03/
C
      XFUVAL = 0.0D0
      DO 15 I = 1,6
      DO 15 J = 1,6
   15 XFUVAL = XFUVAL + C(I,J,NSET)* TX(I) *TT(J)
      XFUVAL = (1.0D0 - X)**3 *  XFUVAL
      XFU = XFUVAL + XFSEA(X,Q)
      RETURN
C
   12 Y = (2.0D0 * LOG(X) + 11.51293D0)/6.90776D0
      TP = (2.0D0 * T - (TMAX + TMIN))/(TMAX - TMIN)
      CALL TCHEBY(TY,Y)
      CALL TCHEBY(TT,TP)
C
      DATA D/
     *0.23946   ,0.29055   ,0.9778E-01,0.2149E-01,0.3440E-02,0.5000E-03,
     *0.1751E-01,-.6090E-02,-.2687E-01,-.1916E-01,-.7970E-02,-.2750E-02,
     *-.5760E-02,-.5040E-02,0.1080E-02,0.2490E-02,0.1530E-02,0.7500E-03,
     *0.1740E-02,0.1960E-02,0.3000E-03,-.3400E-03,-.2900E-03,-.1800E-03,
     *-.5300E-03,-.6400E-03,-.1700E-03,0.4000E-04,0.6000E-04,0.4000E-04,
     *0.1700E-03,0.2200E-03,0.8000E-04,0.1000E-04,-.1000E-04,-.1000E-04,
C
C   D FOR NSET = 1 IS ABOVE; D FOR NSET = 2 IS BELOW
C
     *0.2410    ,0.2884    ,0.9369E-01,0.1900E-01,0.2530E-02,0.2400E-03,
     *0.1765E-01,-.9220E-02,-.3037E-01,-.2085E-01,-.8440E-02,-.2810E-02,
     *-.6450E-02,-.5260E-02,0.1720E-02,0.3110E-02,0.1830E-02,0.8700E-03,
     *0.2120E-02,0.2320E-02,0.2600E-03,-.4900E-03,-.3900E-03,-.2300E-03,
     *-.6900E-03,-.8200E-03,-.2000E-03,0.7000E-04,0.9000E-04,0.6000E-04,
     *0.2400E-03,0.3100E-03,0.1100E-03,0.0000E+00,-.2000E-04,-.2000E-04/
C
      XFUVAL = 0.0D0
      DO 16 I = 1,6
      DO 16 J = 1,6
   16 XFUVAL = XFUVAL + D(I,J,NSET)* TY(I) *TT(J)
      XFUVAL = (1.0D0 - X)**3 *  XFUVAL
      XFU = XFUVAL + XFSEA(X,Q)
      RETURN
      END
