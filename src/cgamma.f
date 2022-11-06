      SUBROUTINE CGAMMA ( KVAL , COSGM , SINGM , DGADRH , XMUZZ ,
     1                    ECOEF1 , ECOEF2 , F2COEF )
      IMPLICIT REAL*8 (A-H,O-Z)
C             :
C DATE        : 24.02.1994
C AUTHOR      : PER JENSEN (WITH HELP FROM R. BEARDSWORTH)
C UPDATES     :
C LANGUAGE    : FORTRAN
C             :
C
C     SUBROUTINE CALCULATES THE VALUES OF VA(RHO), VB(RHO),
C     AND MUZZ(RHO)= F2(RHO)/IPP(RHO).
C
C
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8 B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8 CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5       NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(128) , PARMAX ,
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      REAL*8 I0,I0P,I0PP,IRRVAL,IRRP,IRRPP
      REAL*8 CE,SE,C2R,S2R,CRME,SRME,ECOEF1,ECOEF2,
     1      F1VAL,F2VAL,GVAL,V0VAL,F2COEF,CRHOE
C
      include 'value.h'
      include 'molcul.h'
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
      include 'lzcomp.h'
      COMMON / MACHINE / MEPPOW , LSTPOW
C
      SQ2 = SQRT(2.0D+00)
      EPS=0.5D0*RHO+ECOEF1*ATAN(ECOEF2*TAN(0.5D0*RHO))
      CR=COS(RHO)
      SR=SIN(RHO)
      CE=COS(EPS)
      SE=SIN(EPS)
      C2R=CR*CR-SR*SR
      S2R=2.0E0*CR*SR
      CRME=CR*CE+SR*SE
      SRME=SR*CE-CR*SE
C
      CANG=1.0D+00-CR
      VA =  (((((((GA8*CANG+GA7)*CANG+GA6)*CANG+GA5)*CANG
     1      +GA4)*CANG+GA3)*CANG+GA2)*CANG+GA1)*CANG
      VB =  (((((((HA8*CANG+HA7)*CANG+HA6)*CANG+HA5)*CANG
     1     +HA4)*CANG+HA3)*CANG+HA2)*CANG+HA1)*CANG
C
C       CALCULATE MUZZ
C
      XMUZZ = (U1*CRME*CRME+U3*CE*CE+2.0E0*U13*CE*CRME)
     1      *F2COEF/(SR*SR)
C
      DELTA=VB-VA
      W12=XMUZZ*KVAL*XLAMB
C
      IF (ABS(W12) .EQ. 0.0D+00) THEN
               HELP=1.0D+00
               GOTO 100
      ENDIF
C
      AUX1=DELTA**2+4.0D+00*W12**2
C
      IF (ABS(DELTA) .GT. 0.0D+00) THEN
C
      IF (-LOG(ABS(VB-VA)) .GT. LSTPOW*LOG(2.0D+00)) THEN
                       HELP=0.0D+00
                       GOTO 100
      ENDIF
C
      IF (-LOG(ABS(2.0D+00*W12/(VB-VA))) .GT.
     1                            MEPPOW*LOG(2.0D+00)) THEN
                       HELP=1.0D+00
                       GOTO 100
      ENDIF
      ELSE
             IF (ABS(W12) .EQ. 0.0D+00) THEN
                       HELP=1.0D+00
             ELSE
                       HELP=0.0D+00
             ENDIF
             GOTO 100
      ENDIF
C

      HELP=DELTA/SQRT(AUX1)
C
100   COSGM=SQRT(1.0D+00+HELP)/SQ2
      SINGM=SQRT(1.0D+00-HELP)/SQ2
C
      DVADR=(((((((8.0D+00*GA8*CANG+7.0D+00*GA7)*CANG+
     1    6.0D+00*GA6)*CANG+5.0D+00*GA5)*CANG+4.0D+00*GA4)*CANG
     2    +3.0D+00*GA3)*CANG+2.0D+00*GA2)*CANG+GA1)*SR
      DVBDR=(((((((8.0D+00*HA8*CANG+7.0D+00*HA7)*CANG+
     1    6.0D+00*HA6)*CANG+5.0D+00*HA5)*CANG+4.0D+00*HA4)*CANG
     2    +3.0D+00*HA3)*CANG+2.0D+00*HA2)*CANG+HA1)*SR
C
      DDEDR=DVBDR-DVADR
C
      DEDR=(U1+U13*CR)/(U1+U3+2.0D+00*U13*CR)
C
      DWWDR=-2.0D+00*KVAL*XLAMB*F2COEF*((U1*CRME*CRME+U3*CE*CE
     1                         +2.0E0*U13*CE*CRME)*CR/SR
     1 + (U1*CRME*SRME*(1.0D+00-DEDR)
     1 + U3*CE*SE*DEDR
     1 + U13*(SE*CRME*DEDR+CE*SRME*(1.0D+00-DEDR))))/(SR*SR)
C
      IF (ABS(HELP-1.0D+00) .LT. 1.0D-10) THEN
            DCGDR=0.0D+00
            DSGDR=0.0D+00
            GOTO 150
      ENDIF
C
      DHEDR=(DDEDR*(1.0D+00-DELTA**2/AUX1)
     1     -4.0D+00*DWWDR*W12*DELTA/AUX1)/SQRT(AUX1)
C
      DCGDR= DHEDR/2.0D+00/SQ2/SQRT(1.0D+00+HELP)
      DSGDR=-DHEDR/2.0D+00/SQ2/SQRT(1.0D+00-HELP)
C
150   DGADRH=-(COSGM*DSGDR-SINGM*DCGDR)
      SINGM=-SINGM
C
      RETURN
      END
