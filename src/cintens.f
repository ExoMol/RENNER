C***********************************************
      SUBROUTINE CINTNS(EIG1,EIG2,
     1                   HELP,DIPMAT,
     1                   IEGVEC,
     1                   NRAUS,NFCTS,
     1                   JMAXP1,IEIGMX,MXFCCO,
     2                   MAXPSY,NRECUN,
     3                   JX21,N2XS,
     4                   INTFLG,INTFST,
     5                   NFIL6,NFIL53,
     6                   NFIL51,NFIL52,
     7                   NFIL41,NFIL42,NFIL43)
C***********************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8  XTHR,XCLG
      REAL*8  CLG
      REAL*8  THREEJ
C
      INTEGER IEIGMX,MXFCCO,N2XS
      INTEGER IEGVEC(2,2,2,2)
      INTEGER NFCTS(JMAXP1,4)
      INTEGER JMAXP1
      INTEGER INTFLG,INTF1,INTF2
      LOGICAL INTFST
      INTEGER NRAUS
      INTEGER JX2D,JX21,JX22
      INTEGER MAXTAU,ITAU1,ITAU2
      INTEGER MAXPSY,IPSYM1,IPSYM2
      INTEGER NRECUN,NUNIT
      INTEGER V21,NS1,ISURF1,ICABU1
      INTEGER V22,NS2,ISURF2,ICABU2
      CHARACTER*1 CABU1,CABU2,TRANS
      INTEGER START1,START2 
      INTEGER ANZ1,ANZ2
      INTEGER NR1,NR2  
      INTEGER LGTH1,LGTH2
      INTEGER INDSY1,INDSY2
      INTEGER INDSH1,INDSH2
      REAL*8  XLAB1,XLAB2
      REAL*8  E1,E2
      INTEGER N1,N2,IN1,IN2
      INTEGER K1,K2,IK1,IK2
      INTEGER NMIN1,NMAX1
      INTEGER NMIN2,NMAX2
      INTEGER S1,S2,L1,L2,H
      INTEGER NFEIG1,NFEIG2
      REAL*8 LINSTR,TR1,TR2
      REAL*8 EIG1(IEIGMX),EIG2(IEIGMX),HELP(MXFCCO)
      REAL*8 DIPMAT(MXFCCO,MXFCCO)
      INTEGER NFIL6,NFIL53
      INTEGER NFIL51,NFIL52
      INTEGER NFIL41,NFIL42,NFIL43
C
      include 'intens.h'
C
C
      IF (NRAUS.EQ.0) THEN
        NUNIT=NFIL53
        OPEN(UNIT=NFIL53,STATUS='SCRATCH',ACCESS='DIRECT',
     1       RECL=6*NRECUN,ERR=1009)
      ENDIF
C
C
      DO 10 JX2D=0,2,2
        IF (INTFST.AND.(JX2D.EQ.2)) GOTO 10
        JX22=JX21-JX2D
        NMIN1=IABS(JX21-N2XS)/2
        NMAX1=    (JX21+N2XS)/2
        NMIN2=IABS(JX22-N2XS)/2
        NMAX2=    (JX22+N2XS)/2
        INTF1=INTFLG
        IF (INTF1.EQ.0) THEN
          NFEIG1=NFIL51
        ELSE
          NFEIG1=NFIL52
        ENDIF
        IF (JX21.EQ.JX22) THEN
          INTF2=INTFLG
        ELSE
          INTF2=1-INTFLG
        ENDIF
        IF (INTF2.EQ.0) THEN
          NFEIG2=NFIL51 
        ELSE
          NFEIG2=NFIL52 
        ENDIF
        IF (JX21.EQ.JX22) THEN
          MAXTAU=0
        ELSE
          MAXTAU=1
        ENDIF
        DO 20 IPSYM1=1,MAXPSY
          IPSYM2=3-IPSYM1
          IF (MAXPSY.EQ.1) IPSYM2=1
          DO 30 ITAU1=0,MAXTAU
            ITAU2=1-ITAU1
            INDSY1=2*(IPSYM1-1)+(ITAU1+1)
            INDSY2=2*(IPSYM2-1)+(ITAU2+1)
            START1=IEGVEC(INTF1+1,IPSYM1,ITAU1+1,1)
            ANZ1  =IEGVEC(INTF1+1,IPSYM1,ITAU1+1,2)
            START2=IEGVEC(INTF2+1,IPSYM2,ITAU2+1,1)
            ANZ2  =IEGVEC(INTF2+1,IPSYM2,ITAU2+1,2)
            CALL MATPOS(NMIN1,NMAX1,NMAX1+1,NFCTS,
     1                    JMAXP1,INDSY1,LGTH1)
            CALL MATPOS(NMIN2,NMAX2,NMAX2+1,NFCTS,
     1                    JMAXP1,INDSY2,LGTH2)
            NR1=0
40          NR1=NR1+1
            NR2=0
            IF (NR1.GT.ANZ1) GOTO 30
            READ (NFEIG1,REC=START1+NR1)
     1             E1,XLAB1,(EIG1(II),II=1,LGTH1)
            CALL UNPACK(JX21,N1,K1,V21,NS1,ISURF1,
     1                    CABU1,INDSY1,XLAB1)
            IF (CABU1.EQ.'A') THEN
              ICABU1=1
            ELSE
              ICABU1=2
            ENDIF
50          NR2=NR2+1
            IF (NR2.GT.ANZ2) GOTO 40
              READ (NFEIG2,REC=START2+NR2)
     1             E2,XLAB2,(EIG2(II),II=1,LGTH2)
              CALL UNPACK(JX22,N2,K2,V22,NS2,ISURF2,
     1                    CABU2,INDSY2,XLAB2)
              IF (CABU2.EQ.'A') THEN
                 ICABU2=1
              ELSE
                 ICABU2=2
              ENDIF
              EDIFF=E2-E1
C
              IF (MIN(E1,E2).GT.3000.0) THEN
                IF ((E1.GT.3000.0).AND.(E2.GT.E1)) GOTO 40
                GOTO 50
              ENDIF
C
              IF ((EDIFF.GT.FRQHI1).AND.
     1            (EDIFF.GT.FRQHI2).AND.
     1            (EDIFF.GT.FRQHI3)) GOTO 40
              EDIFF=ABS(EDIFF)
              IF (((EDIFF.GT.FRQLO1).AND.(EDIFF.LT.FRQHI1)).OR.
     1            ((EDIFF.GT.FRQLO2).AND.(EDIFF.LT.FRQHI2)).OR.
     1            ((EDIFF.GT.FRQLO3).AND.(EDIFF.LT.FRQHI3))) 
     1            GOTO 200
100           GOTO 50
200           IF ((.NOT.(((INT1V2.EQ.V21   ).OR.(INT1V2.EQ.-1)).AND.
     1                   ((INT1NS.EQ.NS1   ).OR.(INT1NS.EQ.-1)).AND.
     1                   ((INT1ST.EQ.ICABU1).OR.(INT1ST.EQ.-1)).AND.
     1                   ((INT1SU.EQ.ISURF1).OR.(INT1SU.EQ.-1)).AND.
     1                   ((INT2V2.EQ.V22   ).OR.(INT2V2.EQ.-1)).AND.
     1                   ((INT2NS.EQ.NS2   ).OR.(INT2NS.EQ.-1)).AND.
     1                   ((INT2ST.EQ.ICABU2).OR.(INT2ST.EQ.-1)).AND.
     1                   ((INT2SU.EQ.ISURF2).OR.(INT2SU.EQ.-1))))
     1            .AND.
     1            (.NOT.(((INT2V2.EQ.V21   ).OR.(INT2V2.EQ.-1)).AND.
     1                   ((INT2NS.EQ.NS1   ).OR.(INT2NS.EQ.-1)).AND.
     1                   ((INT2ST.EQ.ICABU1).OR.(INT2ST.EQ.-1)).AND.
     1                   ((INT2SU.EQ.ISURF1).OR.(INT2SU.EQ.-1)).AND.
     1                   ((INT1V2.EQ.V22   ).OR.(INT1V2.EQ.-1)).AND.
     1                   ((INT1NS.EQ.NS2   ).OR.(INT1NS.EQ.-1)).AND.
     1                   ((INT1ST.EQ.ICABU2).OR.(INT1ST.EQ.-1)).AND.
     1                   ((INT1SU.EQ.ISURF2).OR.(INT1SU.EQ.-1))))
     1            ) GOTO 50
              LINSTR=0.0D+00
              DO 60 IN1=NMIN1,NMAX1
              DO 70 IN2=NMIN2,NMAX2
                TR2=0.0D+00
                DO 80 KDEL=-1,1,1
                DO 90 IK1=0,IN1 
                  IK2=IK1+KDEL
                  IF ((IK2.LT.0).OR.(IK2.GT.IN2)) GOTO 90
                  CALL MATPOS(NMIN1,IN1,IK1,NFCTS,JMAXP1,
     1                        INDSY1,S1)
                  CALL MATPOS(NMIN2,IN2,IK2,NFCTS,JMAXP1,
     1                        INDSY2,S2)
                  S1=S1+1
                  S2=S2+1
                  IF ((IK1.EQ.0).AND.(MOD(IN1,2).EQ.1)) THEN
                    IF (MOD(INDSY1,2).EQ.1) THEN 
                      INDSH1=INDSY1+1
                    ELSE
                      INDSH1=INDSY1-1
                    ENDIF
                  ELSE
                    INDSH1=INDSY1
                  ENDIF
                  L1=NFCTS(IK1+1,INDSH1)
                  IF (L1.EQ.0) GOTO 90 
                  IF ((IK2.EQ.0).AND.(MOD(IN2,2).EQ.1)) THEN 
                    IF (MOD(INDSY2,2).EQ.1) THEN
                      INDSH2=INDSY2+1
                    ELSE
                      INDSH2=INDSY2-1
                    ENDIF
                  ELSE
                    INDSH2=INDSY2
                  ENDIF
                  L2=NFCTS(IK2+1,INDSH2)
                  IF (L2.EQ.0) GOTO 90
                  IF (KDEL.EQ.   0) NFDP=NFIL42
                  IF (ITAU1.EQ.0) THEN
                    IF (KDEL.EQ.(-1)) NFDP=NFIL41
                    IF (KDEL.EQ.   1) NFDP=NFIL43
                  ELSE
                    IF (KDEL.EQ.(-1)) NFDP=NFIL43
                    IF (KDEL.EQ.   1) NFDP=NFIL41
                  ENDIF
                  IF (ITAU1.EQ.0) THEN
                    CALL RCNRDP(IPSYM1,
     1                          MOD(IN1,2),MOD(IN2,2),
     2                          IK1+1,IK2+1,KDEL,
     3                          NRCNR)
                  ELSE
                    CALL RCNRDP(IPSYM2,
     1                          MOD(IN2,2),MOD(IN1,2),
     2                          IK2+1,IK1+1,-KDEL,
     3                          NRCNR)
                  ENDIF
                  IF (ITAU1.EQ.0) THEN
                    READ (NFDP,REC=NRCNR)
     1                   ((DIPMAT(II,JJ),II=1,MXFCCO),JJ=1,MXFCCO)
                  ELSE
                    READ (NFDP,REC=NRCNR)
     1                   ((DIPMAT(JJ,II),II=1,MXFCCO),JJ=1,MXFCCO)
                  ENDIF
                  CALL DGEMM('T','N',L2,1,L1,1.0D+00,
     1                       DIPMAT,MXFCCO,
     1                       EIG1(S1),L1,0.0D+00,HELP,MXFCCO)
                  TR1=0.0D+00
                  CALL DGEMM('T','N',1,1,L2,1.0D+00,
     1                       HELP,MXFCCO,
     1                       EIG2(S2),L2,0.0D+00,TR1,1)
                  IF (KDEL.NE.0) TR1=TR1*(-KDEL)
                  IF (((IK1.EQ.0).AND.(IK2.NE.0)).OR.
     1                ((IK2.EQ.0).AND.(IK1.NE.0))) 
     2                   TR1=TR1*SQRT(2.0D+00)
                  IF (MOD(IK1,2).EQ.1) TR1=-TR1
C
                  XTHR= 
     1            THREEJ(
     1                DFLOAT(IN2),1.0D+00,DFLOAT(IN1),
     1                DFLOAT(IK2),DFLOAT(-KDEL),DFLOAT(-IK1)
     1            )
C
C                 XCLG=CLG(IN2,IN1-IN2,IK2,-KDEL)
C                 XCLG=XCLG/SQRT(2.0D+00*DFLOAT(IN1)+1.0D+00)
C                 IF (MOD(IABS(IN2-1+IK1),2).EQ.1) XCLG=-XCLG
C 
                  TR1=TR1*XTHR
C
                  TR2=TR2+TR1
90              CONTINUE
80              CONTINUE
                IF (MOD(IABS(IN1+IN2),2).EQ.1) TR2=-TR2
                TR2=TR2*SQRT((2.0D+00*DFLOAT(IN1)+1.0D+00)
     1                 *     (2.0D+00*DFLOAT(IN2)+1.0D+00))
     1             *SIXJ(DFLOAT(JX22)/2.0D+00,DFLOAT(IN2),
     1                   DFLOAT(N2XS)/2.00D+00,DFLOAT(IN1),
     1                   DFLOAT(JX21)/2.0D+00,1.0D+00)
                LINSTR=LINSTR+TR2
70            CONTINUE
60            CONTINUE
              LINSTR=(DFLOAT(JX21)+1.0D+00)
     1              *(DFLOAT(JX22)+1.0D+00)
     2              *LINSTR*LINSTR
              IF (LINSTR.GT.XINLIM) THEN
                NRAUS=NRAUS+1
                IF (E1.LT.E2) THEN
                   WRITE(NFIL53,REC=NRAUS) 
     1                  XLAB1,XLAB2,E1,E2-E1,LINSTR,
     1                  0.0D+00
                ELSE
                   WRITE(NFIL53,REC=NRAUS)
     1                  XLAB2,XLAB1,E2,E1-E2,LINSTR,
     1                  0.0D+00
                ENDIF
              ENDIF
              GOTO 50 
30        CONTINUE
20      CONTINUE
10    CONTINUE
C
      RETURN
1009  WRITE (NFIL6,2009) NUNIT
2009  FORMAT(1H0,'RENNER.CNT.ERR SCRATCH FILE COULD NOT BE O',
     1           'PENED, UNIT NUMBER IS ',I2)
      STOP
      END
C
C*********************************************
      SUBROUTINE CALINT(NRAUS,PARFUN,ISO,IACT)
C*********************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER PRTINT
      INTEGER IACT
C
C     IACT=1: Absorption
C     IACT=2: Emission
C
      include 'rensys.h'
      include 'intens.h'
C
      INTEGER NRAUS
      INTEGER ISO
      REAL*8 PARFUN
C
      REAL*8 TEMPCO
      REAL*8 TEMPFC
C
      REAL*8 XLAB1,XLAB2,
     1       E1,DE,
     2       LINSTR,INTENS
C
      INTEGER JX21,N1,K1,V21,NS1,ISURF1,INDSY1
      CHARACTER*1 CABU1
C
      INTEGER GNS
C
      TEMPCO=PLANCK*VELLGT/BOLTZ/TEMPRA
      IF (IACT.EQ.1) THEN
         TEMPFC=8.0D+00*PI*PI*PI*1.00D-36
     1         *AVOGNO/3.00D+00/PLANCK/VELLGT
      ELSE
         TEMPFC=64.0D+00*PI*PI*PI*PI*1.00D-36
     1         *AVOGNO*VELLGT/3.00D+00
      ENDIF
C
      DO 10 II=1,NRAUS
        READ (NFIL53,REC=II)
     1       XLAB1,XLAB2,E1,DE,LINSTR,INTENS
        CALL UNPACK(JX21,N1,K1,V21,NS1,ISURF1,
     1              CABU1,INDSY1,XLAB1)
        IF ((INDSY1.EQ.1).OR.(INDSY1.EQ.4)) THEN
          GNS=IGNS(1,ISO)
        ELSE
          GNS=IGNS(2,ISO)
        ENDIF
C
        IF (IACT.EQ.1) THEN 
          INTENS=TEMPFC*DE*DFLOAT(GNS)
     1          *DEXP(-E1*TEMPCO)
     1          *(1.0D+00-DEXP(-DE*TEMPCO))
     1          *LINSTR/PARFUN
        ELSE
          INTENS=TEMPFC*DE*DE*DE*DE*DFLOAT(GNS)
     1          *DEXP(-(E1+DE)*TEMPCO)
     1          *LINSTR/PARFUN
        ENDIF
C      
        WRITE(NFIL53,REC=II)
     1       XLAB1,XLAB2,E1,DE,LINSTR,INTENS
10    CONTINUE
C
      RETURN
      END
C
