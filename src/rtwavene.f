      SUBROUTINE WAVFUN ( EV2K , U , F1 , F2 , I0 , GRHO , RTI0 ,
     1                   VPV , VMV ,
     1                   PHIB , PHPR , XJACO , YRHS , XTX , XTY ,
     2                   CS , AS   , IA1   , IA2  ,
     3                   COSGM , SINGM ,
     3                   LENIAR , NFILBN )
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 AMASS(3,10)
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
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
C
      INTEGER LENREC,LENIAR,I,L,IERR,NUMNOD,NTRIES,NP,II,NPST,
     1       IFAIL,IPHST,IDEST,IEND,NRECS,KREC,KP1,NM1,
     2       NOFCT,IREMAI,IOFFS
C
      REAL*8 CONVIV,FACODD,FACEVE,EOLD,PCOUT,PCIN,OLDPHI,
     1      SUM,ALPHA,SUMM,SUMP,TRIALE,IRRVAL,IRRP,
     2      PHIDER,FCTVAL,HELP1,HELP3,NEWPHI
C
      REAL*8 EETAB(50)
C
      REAL*8 EV2K(V2MXP1,JMAXP1) , F1(NSTINT) , F2(NSTINT) ,
     1      U(NSTINT) , VPV(NSTINT) , VMV(NSTINT) ,
     1      I0(NSTINT) , GRHO(NSTINT) , RTI0(NSTINT) ,
     2      PHIB(NSTINT) , PHPR(NSTINT) , XJACO(NSERIN,NSEPP2) ,
C
     3      YRHS(NSERIN) , XTX(NSEPP2,NSEPP2) , XTY(NSEPP2) ,
C
     4      CS(NSEPP2) , AS(NSEQP1) ,
C
     4      DER(4) , EREST(4),
     1      COSGM(NSTINT*JMAXP1),SINGM(NSTINT*JMAXP1)
C
      INTEGER IA1(LENIAR),IA2(LENIAR)
C
C BEFORE COMMENCING THE NUMEROV-COOLEY NUMERICAL INTEGRATION WE
C STORE THE FOLLOWING FUNCTION IN THE F1 ARRAY:
C
C                       0                 _2
C     W  = F  (P ) + 2 I   (P ) V  (P ) / H
C      I    1   I       PP   I   0   I
C
C      IN THE RENNER-TELLER VERSION OF MORBID, THIS FUNCTION IS
C      EXTENDED WITH A TERM ORIGINATING IN THE NON-ZERO
C      ELECTRONIC ANGULAR MOMENTUM.
C
C                  -1                           2
C V  IS GIVEN IN CM  , AND I   IS GIVEN IN AMU A .
C  0                        PP
C
C THE CONVERSION FACTOR CONVIV CONVERTS UNITS SO THAT WI IS A
C DIMENSIONLESS NUMBER.
C
      include 'isotop.h'
      include 'value.h'
      include 'lzcomp.h'
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               Z1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
C
      IEND=2*NSTINT
      CONVIV=(8.0D+00*PI*PI*VELLGT)/(PLANCK*AVOGNO*1.0D+16)
      IF (NSURF .EQ. 1) THEN
	    VXMIN=VMINS1
      ELSE
	    VXMIN=VMINS2
      ENDIF
C
      REWIND NFIL46
      REWIND NFIL47
      READ (NFIL46) (COSGM(L), L=1,JMAXP1*NSTINT)
      READ (NFIL47) (SINGM(L), L=1,JMAXP1*NSTINT)
C
      DO 200 I=1,NSTINT
      RHO=DFLOAT(I)*HSTEP
      CANG=1.0D+00-COS(RHO)
      RTI0(I)=SQRT(I0(I))
      VA =  (((((((GA8*CANG+GA7)*CANG+GA6)*CANG+GA5)*CANG
     1      +GA4)*CANG+GA3)*CANG+GA2)*CANG+GA1)*CANG
      VB =  (((((((HA8*CANG+HA7)*CANG+HA6)*CANG+HA5)*CANG
     1     +HA4)*CANG+HA3)*CANG+HA2)*CANG+HA1)*CANG
      I0(I)=CONVIV*I0(I)
      VMV(I)=I0(I)*VA
200   VPV(I)=I0(I)*VB
C
C START LOOP OVER K TO CALCULATE WAVEFUNCTIONS
C
      DO 1002 KP1=1,JMAXP1
      INDB=(KP1-1)*NSTINT
C
C CALCULATE EFFECTIVE POTENTIAL
C

      KQUA=KP1-1
      XKSPLS=KQUA*KQUA+XLAMB*XLAMB
      XKL=KQUA*XLAMB
      IF (NSURF .EQ. 1) THEN
      DO 100 I=1,NSTINT
      INDR=INDB+I
100   U(I) = F1(I) + (COSGM(INDR)**2*VMV(I)+SINGM(INDR)**2*VPV(I)
     1     - I0(I)*VXMIN)
     1   + F2(I)*(XKSPLS+4.0D+00*COSGM(INDR)*SINGM(INDR)*XKL) 
      ELSE
      DO 102 I=1,NSTINT
      INDR=INDB+I
102   U(I) = F1(I) + (SINGM(INDR)**2*VMV(I)+COSGM(INDR)**2*VPV(I)
     1     - I0(I)*VXMIN)
     1   + F2(I)*(XKSPLS-4.0D+00*COSGM(INDR)*SINGM(INDR)*XKL) 
      ENDIF
C
C WE NOW START THE LOOP OVER V2 TO GENERATE THE K=0 BENDING
C WAVEFUNCTIONS
C
      FACODD=2.0D+00*HSTEP/3.0D+00
      FACEVE=2.0D+00*FACODD
      DO 1000 V2P1=1,V2MXP1
      DO 212 III=1,50
         EETAB(III)=-1.0
212   CONTINUE
      V2=V2P1-1
      EEL=0.0
      EEU=100000.0
      IF (EETAB(V2P1).GT.0) THEN
         EEL=EETAB(V2P1)-5.0
         EEU=EETAB(V2P1)+5.0
         IF (EEL.LT.0) EEL=0.0 
      ENDIF
      EDEL=1000.0
      NTRIES=0
210   EEM=(EEL+EEU)/2.0
211   EGUESS=EEM
      EOLD=EGUESS
      CALL NUMCOO ( F1 , F2 , U , I0 , GRHO , PHIB ,
     1             PHPR , XJACO , YRHS , XTX , XTY , CS ,
     2             AS , IA1 , IA2 ,
     3             LENIAR , PCOUT, IERR)
      IF (IERR .NE. 0) THEN
           IF (IERR .EQ. 1) THEN
                 EEM=EEM+100.0    
                 EEU=EEU+200.0
                 GOTO 211
           ELSE
                 STOP
           ENDIF
      ENDIF
C
C COUNT THE NUMBER OF NODES IN THE WAVEFUNCTION FOUND
C
      NUMNOD=0
      NM2=NSTINT-2
      OLDPHI=PHIB(1)
      NEWPHI=PHIB(3)
      DO 230 I=2,NM2
      IF ((OLDPHI*PHIB(I) .LT. 0.0D+00) .OR. (PHIB(I) .EQ. 0.0D+00
     1 .AND. NEWPHI*OLDPHI .LT. 0.0D+00)) NUMNOD=NUMNOD+1
      NEWPHI=PHIB(I+2)
230   OLDPHI=PHIB(I)
      IF ((IERR.EQ.0).AND.(NUMNOD+1.LT.51)) EETAB(NUMNOD+1)=EGUESS
      IF (NUMNOD .EQ. V2) GOTO 240
C     IF (IPRINT .NE. 0) WRITE (NFIL6,5000) KQUA,V2,NUMNOD,EGUESS,EOLD
      IF (ABS(EEU-EEL).LT.0.01) THEN
         NTRIES=0
         IF (V2P1.EQ.1) THEN
            EEL=0.0
         ELSE
            EEL=EV2K(V2P1-1,KQUA+1)
         ENDIF
         EEU=EEL+EDEL
         EDEL=EDEL+500.0
      ELSE
        IF (NUMNOD.LT.V2) THEN
           EEL=EEM
        ELSE
           EEU=EEM
        ENDIF
      ENDIF
238   EOLD=EGUESS
      NTRIES=NTRIES+1
      IF (NTRIES .GE. NSTNIN) GOTO 9000
      GOTO 210
240   EV2K(V2P1,KQUA+1)=EGUESS
C
C MULTIPLY THE WAVEFUNCTION WITH SQRT(IRR) (EQ. (6.4) OF JENSEN)
C
      DO 300 I=1,NSTINT
300   PHIB(I)=RTI0(I)*PHIB(I)
C
C NORMALIZE THE WAVEFUNCTION
C
      SUM=0.0E0
      DO 400 I=1,NSTINT,2
400   SUM=SUM+FACODD*PHIB(I)*PHIB(I)+FACEVE*PHIB(I+1)*PHIB(I+1)
      SUM=SQRT(SUM)
      DO 500 I=1,NSTINT
500   PHIB(I)=PHIB(I)/SUM
C
C WE DIFFERENTIATE THE FIRST NSERIN POINTS OF THE WAVEFUNCTION
C ANALYTICALLY
C
      RHO=0.0D+00
      IF (NSURF .EQ. 1) THEN
           ALPHA=ABS(KQUA-XLAMB)+0.5D+00
      ELSE
           ALPHA=ABS(KQUA+XLAMB)+0.5D+00
      ENDIF
      DO 700 NP=1,NSERIN
      RHO=RHO+HSTEP
      CR=COS(RHO)
      SR=SIN(RHO)
      IRRVAL=(U1*U3-U13*U13*CR*CR)/(U1+U3+2.0D+00*U13*CR)
      HELP1=U1+U13*CR
      HELP3=U3+U13*CR
      IRRP=2.0D+00*U13*SR*HELP1*HELP3/(HELP1+HELP3)**2
      SUMM=AS(1)
      SUMP=0.0D+00
      DO 600 II=2,NSEQP1
      I=II-1
      SUMM=SUMM+AS(II)*RHO**(I+I)
600   SUMP=SUMP+DFLOAT(I+I)*AS(II)*RHO**(I+I-1)
      PHIDER=(SUMM*ALPHA/RHO+SUMP)*EXP(LOG(RHO)*ALPHA)
      FCTVAL=SUMM*EXP(LOG(RHO)*ALPHA)
700   PHPR(NP)=RTI0(NP)*(PHIDER+0.5D+00*FCTVAL*IRRP/IRRVAL)
     1                                  /(SUM*PCOUT)
C
C THE REST OF THE WAY WE USE THE NAGLIB ROUTINE D04AAF
C
      NPST=NSERIN+1
      HBASE=1.0D+00
      DO 800 NP=NPST,NSTINT
      CALL D04AAF( NP , 1 , DER , EREST , PHIB , IFAIL , NSTINT)
      IF (IFAIL .NE. 0 .OR. ABS(DER(1)) .GT. THRSH4) GOTO 9010
      PHPR(NP)=DER(1)/HSTEP
800   CONTINUE
C
      KREC=KQUA*V2MXP1+V2P1
C
      WRITE (NFILBN,REC=KREC) (PHIB(I),I=1,NSTINT),
     1                       (PHPR(I),I=1,NSTINT)
1000  CONTINUE
1002  CONTINUE
C
      IF (IPRINT.EQ.0) RETURN
      WRITE (NFIL6,5010)
C
      DO 1200 KP1Z=1,JMAXP1
      WRITE (*,*) 'KP1=',KP1Z
      DO 1200 V2P1=1,V2MXP1
      V2=V2P1-1
      WRITE (NFIL6,5020) V2,EV2K(V2P1,KP1Z)
1200  CONTINUE
      WRITE (NFIL6,5030)
C
      RETURN
C
9000  WRITE (NFIL6,9001) NSTNIN
9001  FORMAT(1H0,' RENTEL.NMI.ERR  MAXIMUM NUMBER OF ATTEMPTS (',
     1      I3,') ATTAINED IN THE NUMEROV-COOLEY INTEGRATION')
      STOP
9010  CONTINUE
      IF (IFAIL .EQ. 0) GOTO 9020
      WRITE (NFIL6,9011) IFAIL
9011  FORMAT(1H0,' RENTEL.NMI.ERR   DIFFERENTIAL OF WAVE FUNCTIONS',
     1          ' YIELDS IFAIL =',I2/
     2      1H ,'                  (SUBROUTINE WAVFUN)'/)
      STOP
9020  WRITE (NFIL6,9021) DER(1),THRSH4,NP
9021  FORMAT(1H0,' RENTEL.NMI.ERR   ABSOLUTE ERROR ESTIMATE FOR THE ',
     1                          'FIRST DERIVATIVE OF THE WAVEFUNCTION'/
     2      1H ,'                  (',D11.3,') IS GREATER THAN THE ',
     3                             'LIMIT (',D11.3,')'/
     4      1H ,'                  IN STEP ',I4,' (SUBROUTINE WAVFUN)'/)
      STOP
C
5000  FORMAT(1H0,' RENTEL.NMI.WRN  FOR KQUA = ',I2,', V2 = ',I3,
     1      ' DOES NOT EQUAL NUMBER OF NODES = ',I3/
     2      1H ,' EGUESS = ',F10.3,', ETRIAL = ',F10.3)
5010  FORMAT('1',5X,12('*'),' ENERGIES FROM THE K=0 NUMEROV',
     1      '-COOLEY INTEGRATION (CM-1) ',12('*')//
     2      ' ',17X,30('-')/' ',17X,'I',9X,'I',18X,'I'/
     3      ' ',17X,'I',4X,'V2',3X,'I',8X,'E ',8X,'I'/
     4     ' ',17X,'I',9X,'I',18X,'I'/' ',17X,30('-'))
5020  FORMAT(' ',17X,'I',3X,I3,3X,'I',3X,F12.5,3X,'I')
5030  FORMAT(' ',17X,30('-')//)
C
      END
C
C
      SUBROUTINE NUMCOO ( F1 , F2 , U , I0 , GRHO , PHIB ,
     1                   PHPR , XJACO , YRHS , XTX , XTY , CS ,
     2                   AS , IA1 , IA2 ,
     3                   LENIAR , PCOUT, IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 AMASS(3,10)
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
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER LENIAR,KK,I,IC,IERR,NITER,NP1,NM1
      REAL*8 Y,HH,SUM,SUMOUT,SUMIN,PCOUT,PCIN,DE,YC,YCM1,YCP1
C
      REAL*8 F1(NSTINT) , F2(NSTINT) , U(NSTINT) ,
C
     1      I0(NSTINT) , GRHO(NSTINT) , PHIB(NSTINT) , PHPR(NSTINT) ,
C
     2      XJACO(NSERIN,NSEPP2) , YRHS(NSERIN) , XTX(NSEPP2,NSEPP2) ,
C
     3      XTY(NSEPP2) , CS(NSEPP2) , AS(NSEQP1)
C
      INTEGER IA1(LENIAR),IA2(LENIAR)
C
      include 'isotop.h'
      include 'value.h'
      include 'lzcomp.h'
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               Z1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
C
C     SUBROUTINE SOLVES EQN 6.5 , FOR THE GIVEN VALUE OF KQUA ,
C     BY THE NUMEROV - COOLEY TECHNIQUE.
C
C                  (                       0                  )
C      2           (                    2 I                   )
C     D            (                2      PP                 )
C     ---PHI (P) = ( F (P) + F (P) K  + ----- ( V   (P) - E ) ) PHI (P)
C       2   B      (  1       2         _2       EFF          )    B
C     DP           (                    H                     )
C
C     FOR THIS THE FOLLOWING DEFINTIONS (EQUATIONS 6.14 TO 6.17) ARE
C     NEEDED :-
C
C     PH  = PHI  (P )
C       I      B   I
C
C             0         _2
C     I  = 2 I   (P ) / H
C      I      PP   I
C
C                             2      0                   _2
C     U  = F  (P ) + F  (P ) K  + 2 I   (P ) V    (P ) / H
C      I    1   I     2   I          PP   I   EFF   I
C
C                 2
C     Y  = ( 1 - H  ( U  - I  E ) ) PH
C      I         --    I    I         I
C                12
C
C     STARTING FROM AN INITIAL GUESS FOR THE ENERGY , GIVEN IN WAVFUN
C     AN ITERATIVE SCHEME IS FOLLOWED. THIS IS PERFORMED BY USING THE
C     RECURSION RELATION (EQN 6.18) :-
C
C                           2
C     Y    + Y    - 2 Y  = H  ( U  - I  E ) PH
C      I+1    I-1      I         I    I       I
C
C     H IS HERE THE HSTEP LENGTH USED TO GENERATE THE F1 , F2 , I0'S
C
C     TO CALCULATE VALUES OF Y(I) AND PHIB(I) BOTH FROM RHO = 0
C     AND FROM RHO = RHOMAX. THESE TWO PARTS ARE CALLED THE OUTWARD
C     AND INWARD INTEGRATIONS RESPECTIVELY.
C     THE OUTWARD INTEGRATION IS CARRIED OUT UNTIL A MAXIMUM IS REACHED
C     IN THE WAVEFUNCTION. THIS CROSSING POINT (IC , YC , PC) IS THEN
C     THE STOPPING POINT USED FOR THE INWARD INTEGRATION.
C     BOTH SETS OF WAVEFUNCTIONS ARE SCALED SO THAT PC(OUT) = PC(IN) = 1
C     AN ERROR FOR THE ENERGY CAN BE CALCULATED FROM EQN 6.24 :-
C
C                                                         NSTINT
C               OUT           IN       2                    \--      2
C     D(E) = (-Y    + 2 Y  - Y    ) / H + U  - I  E ) Y  /   \  I  PH
C               C-1      C    C+1          C    C      C     /   I   I
C                                                           /--
C                                                           I=1
C
C     AND THIS IS ADDED TO THE EGUESS AT EACH ITERATION , WHEN D(E)
C     IS LESS THAN THRSH5 THE ITERATION HAS CONVERGED.
C
      Y(I)=(1.0D+00-HH*(U(I)-I0(I)*EGUESS)/12.0D+00)*PHIB(I)
C
      HH=HSTEP*HSTEP
C
      PHIB(NSTINT)=0.0D+00
      PHIB(NSTINT-1)=PNM1
C
      NITER=0
C
100   CALL INTOUT ( F1 , F2 , U , I0 , GRHO , PHIB ,
     1             PHPR , XJACO , YRHS , XTX , XTY , CS ,
     2             AS , IA1 , IA2 ,
     3             LENIAR , SUMOUT , IC ,
     4             PCOUT, IERR)
      IF (IERR.NE.0) RETURN
C
      CALL INTIN ( F1 , F2 , U , I0 , GRHO , PHIB ,
     3            SUMIN , IC , PCIN)
C
      SUM=SUMIN/(PCIN*PCIN)+SUMOUT/(PCOUT*PCOUT)
C
      YC=Y(IC)/PCOUT
C
      YCM1=Y(IC-1)/PCOUT
C
      YCP1=Y(IC+1)/PCIN
C
      DE=((-YCM1+YC+YC-YCP1)/HH+(U(IC)-I0(IC)*EGUESS))*YC/SUM
C
      EGUESS=EGUESS+DE
C
      NITER=NITER+1
C
      IF (ABS(DE).LE.THRSH5) GO TO 150
C
      IF (NITER.EQ.NSTNIN) GO TO 900
C
      GO TO 100
C
C
150   DO 200 I=1,IC
C
      PHIB(I)=PHIB(I)/PCOUT
C
200   CONTINUE
C
      NP1=IC+1
      NM1=NSTINT-1
C
      DO 300 I=NP1,NM1
C
      PHIB(I)=PHIB(I)/PCIN
C
300   CONTINUE
C
C
      RETURN
C
900   WRITE (NFIL6,1000) NSTNIN , DE , THRSH5
      IERR=2
      RETURN
C
1000  FORMAT(1H ,' RENTEL.NMI.ERR   ITERATION LIMIT REACHED IN',
     1          ' NUMERICAL INTEGRATION (NITER=',I2,')',/,
     2                      16X,'ENERGY DIFFERENCE =',D13.4,/,
     3                      16X,'CONVERGENCE THRESHOLD =',D13.4,/)
C
      END
C
C
      SUBROUTINE INTOUT ( F1 , F2 , U , I0 , GRHO , PHIB ,
     1                   PHPR , XJACO , YRHS , XTX , XTY , CS ,
     2                   AS , IA1 , IA2 ,
     3                   LENIAR , SUMOUT , IC ,
     4                   PCOUT, IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 AMASS(3,10)
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
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER LENIAR,KK,I,II,IC,NP1,NM1,IST,IEND,IERR
      REAL*8 HH,CONST,SUM,REDFAC,SUMOUT,PCOUT,YI,YIM1,YIP1
C
      REAL*8 F1(NSTINT) , F2(NSTINT) , U(NSTINT) ,
C
     1      I0(NSTINT) , GRHO(NSTINT) , PHIB(NSTINT) , PHPR(NSTINT) ,
C
     2      XJACO(NSERIN,NSEPP2) , YRHS(NSERIN) , XTX(NSEPP2,NSEPP2) ,
C
     3      XTY(NSEPP2) , CS(NSEPP2) , AS(NSEQP1)
C
      INTEGER IA1(LENIAR),IA2(LENIAR)
C
      include 'isotop.h'
      include 'value.h'
      include 'lzcomp.h'
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               Z1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
C
C     SUBROUTINE PERFORMS THE OUTWARD INTEGRATION OF NUMCOO ,
C     UNTIL THE FIRST MAXIMUM IN THE WAVE FUNCTION IS FOUND.
C     THE SUM I0(I)*PHIB(I)**2 IS SAVED FOR THE OUTWARD INTEGRATION
C     AND WILL LATER BR DIVIDED BY PC(OUT)**2.
C     THE FIRST NSERIN POINTS OF THE FUNCTION ARE CALCULATED BY
C     LSQCAP AS DESCRIBED THEREIN
C
      HH=HSTEP*HSTEP
C
      CALL LSQCAP ( F1 , F2 , U , I0 , GRHO , PHIB ,
     1             PHPR , XJACO , YRHS , XTX , XTY , CS ,
     2             AS , IA1 , IA2 , LENIAR , IERR)
      IF (IERR.NE.0) RETURN
C
      SUM=0.0D+00
C
      NM1=NSERIN-1
C
      DO 100 I=1,NM1
C
      SUM=SUM+I0(I)*PHIB(I)*PHIB(I)
C
100   CONTINUE
C
      IST=NSERIN
      IEND=NSTINT-2
C
      YIM1=PHIB(NM1)*(1.0D+00-HH*(U(NM1)-I0(NM1)*EGUESS)/12.0D+00)
C
      YI=PHIB(NSERIN)*(1.0D+00-HH*(U(NSERIN)-I0(NSERIN)*EGUESS)
     1                        /12.0D+00)
C
C
      DO 200 I=IST,IEND
C
      CONST=HH*(U(I)-I0(I)*EGUESS)
C
      PHIB(I)=YI/(1.0D+00-CONST/12.0D+00)
C
      YIP1=CONST*PHIB(I)+YI+YI-YIM1
C
      IF (PHIB(I).LT.PHIB(I-1)) GO TO 300
C
      YIM1=YI
      YI=YIP1
C
      SUM=SUM+I0(I)*PHIB(I)*PHIB(I)
C
      IF (SUM .LT. THRSH3) GOTO 200
C
      REDFAC=SQRT(THRSH3)
C
      DO 150 II=1,I
C
150   PHIB(II)=PHIB(II)/REDFAC
C
      DO 160 II=1,NSEQP1
C
160   AS(II)=AS(II)/REDFAC
C
      YI=YI/REDFAC
C
      YIM1=YIM1/REDFAC
C
      SUM=SUM/THRSH3
C
200   CONTINUE
C
      GO TO 900
C
300   SUMOUT=SUM
C
      IC=I-1
C
      PCOUT=PHIB(IC)
C
      RETURN
C
900   CONTINUE
      IF (IPRINT .NE. 0) WRITE (NFIL6,1000)
      IERR=1
C
      RETURN
C
1000  FORMAT(1H ,' RENTEL.NMI.ERR   NO TURNING POINT FOUND IN ',
     1          'OUTWARD INTEGRATION',/)
C
      END
C
C
C
      SUBROUTINE INTIN ( F1 , F2 , U , I0 , GRHO , PHIB ,
     3                  SUMIN , IC , PCIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 AMASS(3,10)
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
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER KK,I,IC,II,NP1,NM1,IST,IEND
      REAL*8 HH,CONST,SUM,SUMIN,PCIN,YI,YIM1,YIP1,CONS
C
      REAL*8 F1(NSTINT) , F2(NSTINT) , U(NSTINT) ,
C
     1      I0(NSTINT) , GRHO(NSTINT) , PHIB(NSTINT)
      include 'isotop.h'
      include 'value.h'
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               Z1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
C
C     SUBROUTINE PERFORMS THE INWARD INTEGRATION FROM RHO=RHOMAX
C     STOPPING AT THE MAXIMUM IN THE WAVEFUNCTION DETERMINED IN
C     INTOUT. THE RECURSION RELATIONS GIVEN IN NUMCOO ARE USED
C     IN CONJUNCTION WITH THE TWO STARTING VALUES :-
C
C     PHI(RHOMAX) = 0
C
C     PHI(RHOMAX-H) = PNM1 ( O 10**-32 )
C
C
      HH=HSTEP*HSTEP
C
      SUM=0.0D+00
C
      IST=IC+1
      IEND=NSTINT-1
C
      NM1=IEND
C
      YIP1=0.0D+00
C
      YI=PHIB(NM1)*(1.0D+00-HH*(U(NM1)-I0(NM1)*EGUESS)/12.0D+00)
C
      DO 200 II=IST,IEND
      I=IST+IEND-II
C
      CONST=HH*(U(I)-I0(I)*EGUESS)
C
      PHIB(I)=YI/(1.0D+00-CONST/12.0D+00)
C
      YIM1=CONST*PHIB(I)+YI+YI-YIP1
C
      YIP1=YI
      YI=YIM1
C
      SUM=SUM+I0(I)*PHIB(I)*PHIB(I)
C
200   CONTINUE
C
      SUMIN=SUM
C
      PCIN=YI/(1.0D+00-HH*(U(IC)-I0(IC)*EGUESS)/12.0D+00)
C
      RETURN
C
      END
C
C
      SUBROUTINE LSQCAP ( F1 , F2 , U , I0 , GRHO , PHIB ,
     1                   PHPR , XJACO , YRHS , XTX , XTY , CS ,
     2                   AS , IA1 , IA2 , LENIAR , IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M1,M2,M3,M,U1,U3,U13,V
C
      REAL*8 AMASS(3,10)
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
      INTEGER IQUANT(9,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
      INTEGER LENIAR,KK,I,NP1,NM1,IST,IEND,IERR,II,JJ,
     1       NPEND,NI,IRET,N,J
      REAL*8 HH,SUM,SUMOUT,PCOUT,YI,YIM1,YIP1,DET,TEST,
     1      RHOA,SUMRHS,ALPHA,COEFF
C
      REAL*8 F1(NSTINT) , F2(NSTINT) , U(NSTINT) ,
C
     1      I0(NSTINT) , GRHO(NSTINT) , PHIB(NSTINT) , PHPR(NSTINT) ,
C
     2      XJACO(NSERIN,NSEPP2) , YRHS(NSERIN) , XTX(NSEPP2,NSEPP2) ,
C
     3      XTY(NSEPP2) , CS(NSEPP2) , AS(NSEQP1)
C
      INTEGER IA1(LENIAR),IA2(LENIAR)
C
      include 'isotop.h'
      include 'value.h'
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               Z1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      include 'rentel.h'
      include 'integ.h'
      include 'dimen.h'
      include 'lsfit.h'
      include 'rensys.h'
      include 'bcoeff.h'
      include 'crcoef.h'
      include 'morse.h'
      include 'modim.h'
C
C     SUBROUTINE FITS THE FIRST NSERIN VALUES OF
C
C     T  = U  - I  E
C      I    I    I  GUESS
C
C     TO THE SERIES EXPANSION IN EQN. 6.30 :-
C
C                     P
C                   \--
C            C       \       2I
C     T(P) =  -2  +  /  C   P
C            ---    /--  2I
C             2     I=0
C            P
C
C     THIS IS DONE WITH A LEAST SQUARES METHOD , WITH THE EQUATIONS
C     BEING SOLVED BY THE ROUTINE SOLVD.
C     FROM THE VALUES OF C2I , THE VALUE ALPHA GIVEN BY (EQN. 6.34) :-
C
C                                    1/2
C     ALPHA = 1/2 ( 1 + ( 1 + 4 C   )    )
C                                -2
C
C     AND THE A2I'S GIVEN BY (EQN. 6.35) :-
C                                               P
C                                             \--
C    ((ALPHA+N+1)*(ALPHA+N+1) - C   ) A    =   \   C   A
C    (                           -2 )  N+2     /    2I  N-2I
C                                             /--
C                                             I=0
C
C     ARE CALCULATED.
C
C     FROM THE A2I'S VALUES OF PHIB(RHO) MAY BE CALCULATED FROM
C     EQN 6.32 :-
C                        Q
C                ALPHA \--       2I
C     PHIB(P) = P       \   A   P
C                       /    2I
C                      /--
C                      I=0
C
C
C
      IERR=0
      RHO=0.0D+00
C
      DO 50 N=1,NSERIN
C
      RHO=RHO+HSTEP
C
      DO 40 I=1,NSEPP2
C
40    XJACO(N,I)=RHO**(I+I-4)
C
      YRHS(N)=U(N)-I0(N)*EGUESS
C
50    CONTINUE
C
      DO 340 I=1,NSEPP2
      SUMRHS=0.0D+00
      DO 320 J=1,NSEPP2
      SUM=0.0D+00
      DO 310 N=1,NSERIN
310   SUM=SUM+XJACO(N,I)*XJACO(N,J)
320   XTX(I,J)=SUM
      DO 330 N=1,NSERIN
330   SUMRHS=SUMRHS+XJACO(N,I)*YRHS(N)
340   XTY(I)=SUMRHS
      CALL SOLVD (XTX,IA1,IA2,NSEPP2,NSEPP2,LENIAR,DET,TEST)
      IF (TEST .NE. 0.0D+00) GOTO 900
      DO 360 I=1,NSEPP2
      SUM=0.0D+00
      DO 350 J=1,NSEPP2
350   SUM=SUM+XTX(I,J)*XTY(J)
360   CS(I)=SUM
C
      IF (CS(1).LT.-0.25D+00) CS(1)=-0.25D+00
C
      ALPHA=(1.0D+00+SQRT(1.0D+00+4.0D+00*CS(1)))/2.0D+00
C
      AS(1)=1.0D+00
C
      IF (NSEQP1.LT.2) GO TO 120
C
      DO 100 JJ=2,NSEQP1
C
      N=JJ+JJ-4
C
      COEFF=(ALPHA+DFLOAT(N+1))*(ALPHA+DFLOAT(N+2))-CS(1)
C
      NPEND=MIN0(NSEPP2,JJ)
C
      SUM=0.0D+00
      DO 90 II=2,NPEND
      I=II-2
C
      NI=JJ-I-1
C
      SUM=SUM+CS(II)*AS(NI)
C
90    CONTINUE
C
      AS(JJ)=SUM/COEFF
C
100   CONTINUE
C
120   RHO=0.0D+00
C
      DO 200 N=1,NSERIN
C
      RHO=RHO+HSTEP
C
      RHOA=RHO**ALPHA
C
      SUM=0.0D+00
      DO 150 II=1,NSEQP1
      I=II-1
C
      SUM=SUM+AS(II)*RHO**(I+I)
C
150   CONTINUE
C
      PHIB(N)=RHOA*SUM
C
200   CONTINUE
C
      RETURN
C
900   IERR=5
      WRITE (NFIL6,5000) TEST
5000  FORMAT(1H0,' RENTEL.SSO.ERR   SOLVD HAS FAILED IN SUBROUTINE',
     1          ' LSQCAP, TEST = ',D11.3)
      RETURN
C
      END
