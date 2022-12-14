      SUBROUTINE GENSTR (                           IW1    ,
     1                   IW3    , OVER   , TMAT   , ESTR   ,
     1                   SMAT   , VMAT   , WMAT   , NOPERS )
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
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
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 , NFIL6 ,
     5       NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
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
      INTEGER        IOSC1,IOSC3,I,J,IR,JL
      REAL*8 G,GC1,GC3,GC1113,GC3133,GC1133,GC3113
      INTEGER IW1(LENIW),IW3(LENIW),V2L,V2R,ISTA1,ISTA3,IJMP1,
     1       IJMP3,IR1,IR3,IR11,IR33,IR13,IR111,IR333,IR113,IR133,
     2       IR1111,IR3333,IR1113,IR1333,IR1133,
     3       IV1,IV3,IV11,IV33,IV13,IV111,IV333,IV113,IV133,
     4       IV1111,IV3333,IV1113,IV1333,IV1133,
     5       IG0,IG1,IG3,IG11,IG33,IG13,
     6       IC11,IC33,IC13,IC31,IC111,IC333,IC133,IC311,
     7       IC113,IC313,IC1111,IC3333,IC1333,IC3111,IC1113,
     8       IC3133,IC1133,IC3113
      INTEGER INDEXL,INDEXR,IOFFS,JOFFS,INDVVP,INDPVV,
     1       LABELA,LABELB,LABELC,LVP1,LVP3,LPV1,LPV3,
     2       LVP11,LVP33,LPV11,LPV33,LVP13,LPV13,
     3       LVP111,LVP333,LPV111,LPV333,LVP113,LPV113,
     4       LVP133,LPV133,LABVVP,LABPVV
      INTEGER DIMA,KSTYP1,KSTYP2,LSTYP1,LSTYP2,IOFFSA,JOFFSA,
     1       IOFFSB,JOFFSB
      REAL*8                                     OVER(2,MBASP1,MBASP1),
     1       ESTR(NOBAS),TMAT(NOBAS,NOBAS),
     1       SMAT(NOBAS,NOBAS,NOPERS),WMAT(NOBAS,NOBAS),
     1       VMAT(NOBAS,NOBAS)
      REAL*8  CONCOR,CONFC1,CONG11,CONG33,CONG13,CORHOE
      COMMON /ISOTOP/ IQUANT,AMASS
C
      COMMON /EQUFCS/W11,W33,W13,W111,W333,W113,W133,
     1               W1111,W3333,W1113,W1333,W1133
C
C
      COMMON /VALUES/ RHO,EPS,EPSP,EPSPP,EPSPPP,
     1               CR,SR,CSE,SNE,CRE,SRE,CORO
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     1               F11,F33,F13,F111,F333,F113,F133,
     1               F1111,F3333,F1113,F1333,F1133,
     2               F1,F1A1,F2A1,F3A1,F4A1,F3,F1A3,F2A3,F3A3,F4A3,
     3               F1A11,F2A11,F3A11,F1A33,F2A33,F3A33,
     4               F1A13,F2A13,F3A13,
     5               F1A111,F2A111,F1A333,F2A333,
     6               F1A113,F2A113,F1A133,F2A133,
     7               FA1111,FA3333,FA1113,
     8               FA1333,FA1133, R12RF1, R32RF1, R12RF2, R32RF2,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
      COMMON /RENTEL/G1,G3,G11,G33,G13,G111,G333,G113,G133,
     1               G1111,G3333,G1113,G1333,G1133,
     2               GA1,GA2,GA3,GA4,GA5,GA6,GA7,GA8,
     2               G1A1,G2A1,G3A1,G4A1,G1A3,G2A3,G3A3,G4A3,
     3               G1A11,G2A11,G3A11,G1A33,G2A33,G3A33,
     4               G1A13,G2A13,G3A13,
     5               G1A111,G2A111,G1A333,G2A333,
     6               G1A113,G2A113,G1A133,G2A133,
     7               GA1111,GA3333,GA1113,
     8               GA1333,GA1133,
     1               HA1,HA2,HA3,HA4,HA5,HA6,HA7,HA8,
     2               H1A1,H2A1,H3A1,H4A1,H1A3,H2A3,H3A3,H4A3,
     3               H1A11,H2A11,H3A11,H1A33,H2A33,H3A33,
     4               H1A13,H2A13,H3A13,
     5               H1A111,H2A111,H1A333,H2A333,
     6               H1A113,H2A113,H1A133,H2A133,
     7               HA1111,HA3333,HA1113,
     8               HA1333,HA1133
C
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMINS1 , VMINS2 , ZTRIAL(2) , V0TYPE ,                EN01300
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
C
      COMMON /DIMEN/ NSURF,  MBASIS ,  V2MAX  , V2MXP1 ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ , NFINTA
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      include 'rensys.h'
C
      COMMON/BCOEFF/
     1      B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      COMMON/CRCOEF/
     1      CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      KSTYP1=KSTYPA(1)
      KSTYP2=KSTYPA(2)
      LSTYP1=LSTYPA(1)
      LSTYP2=LSTYPA(2)
      DIMA=V2MXP1*LSTYP1
      IOSC1=1
      IOSC3=2
      IF ( SYMM ) IOSC3=1
C
C 1. GENERATE MATRIX OF Y1
C
      CALL YONLY(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  2. GENERATE MATRIX OF Y3
C
      CALL YONLY(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  3. GENERATE MATRIX OF Y1*Y1
C
      CALL YY(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  4. GENERATE MATRIX OF Y3*Y3
C
      CALL YY(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C   5. GENERATE MATRIX OF Y1*Y3
C
      CALL Y1Y2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  6. GENERATE MATRIX OF Y1*Y1*Y1
C
      CALL YYY(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  7. GENERATE MATRIX OF Y3*Y3*Y3
C
      CALL YYY(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN(VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  8. GENERATE MATRIX OF Y1*Y1*Y3
C
      CALL Y1Y1Y2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  9. GENERATE MATRIX OF Y1*Y3*Y3
C
      CALL Y1Y1Y2(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  10. GENERATE MATRIX OF Y1*Y1*Y1*Y1
C
      CALL YYYY(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  11. GENERATE MATRIX OF Y3*Y3*Y3*Y3
C
      CALL YYYY(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  12. GENERATE MATRIX OF Y1*Y1*Y1*Y3
C
      CALL Y13Y2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  13. GENERATE MATRIX OF Y1*Y3*Y3*Y3
C
      CALL Y13Y2(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  14. GENERATE MATRIX OF Y1*Y1*Y3*Y3
C
      CALL Y12Y22(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  15. GENERATE MATRIX OF Y1*P1+P1*Y1
C
      CALL YP(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  16. GENERATE MATRIX OF Y3*P3+P3*Y3
C
      CALL YP(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  17. GENERATE MATRIX OF P1*Y3
C
      CALL P1Y2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  18. GENERATE MATRIX OF P3*Y1
C
      CALL P1Y2(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  19. GENERATE MATRIX OF Y1*Y1*P1+P1*Y1*Y1
C
      CALL Y2P(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  20. GENERATE MATRIX OF Y3*Y3*P3+P3*Y3*Y3
C
      CALL Y2P(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  21. GENERATE MATRIX OF P1*Y3*Y3
C
      CALL P1Y2Y2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  22. GENERATE MATRIX OF P3*Y1*Y1
C
      CALL P1Y2Y2(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  23. GENERATE MATRIX OF (Y1*P1+P1*Y1)*Y3
C
      CALL YPY(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  24. GENERATE MATRIX OF (Y3*P3+P3*Y3)*Y1
C
      CALL YPY(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  25. GENERATE MATRIX OF Y1*Y1*Y1*P1+P1*Y1*Y1*Y1
C
      CALL Y3P(IOSC1,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  26. GENERATE MATRIX OF Y3*Y3*Y3*P3+P3*Y3*Y3*Y3
C
      CALL Y3P(IOSC3,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  27. GENERATE MATRIX OF P1*Y3*Y3*Y3
C
      CALL PY3(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  28. GENERATE MATRIX OF P3*Y1*Y1*Y1
C
      CALL PY3(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  29. GENERATE MATRIX OF (Y1*Y1*P1+P1*Y1*Y1)*Y3
C
      CALL Y2PY(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  30. GENERATE MATRIX OF (Y3*Y3*P3+P3*Y3*Y3)*Y1
C
      CALL Y2PY(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  31. GENERATE MATRIX OF (Y1*P1+P1*Y1)*Y3*Y3
C
      CALL YPY2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  32. GENERATE MATRIX OF (Y3*P3+P3*Y3)*Y1
C
      CALL YPY2(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  33. GENERATE MATRIX OF P1*P3
C
      CALL P1P2(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  34. GENERATE MATRIX OF (Y1*P1+P1*Y1)*P3
C
      CALL PYYPP(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  35. GENERATE MATRIX OF (Y3*P3+P3*Y3)*P1
C
      CALL PYYPP(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  36. GENERATE MATRIX OF (Y1*Y1*P1+P1*Y1*Y1)*P3
C
      CALL YYPP(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  37. GENERATE MATRIX OF (Y3*Y3*P3+P3*Y3*Y3)*P1
C
      CALL YYPP(IOSC3,IOSC1,IW3,IW1,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C  38. GENERATE MATRIX OF (Y1*P1+P1*Y1)*(Y3*P3+P3*Y3)
C
      CALL YPPY(IOSC1,IOSC3,IW1,IW3,OVER,VMAT)
      KNDEX=KNDEX+1
      CALL SIMTRN (VMAT,TMAT,SMAT(1,1,KNDEX),WMAT,NOBAS)
C
C
C
      RETURN
      END
C
C
      SUBROUTINE SIMTRN ( A , U , B , W , N )
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,I,J,K
      REAL*8 A(N,N),U(N,N),B(N,N),W(N,N),POM
C
      DO 20 I=1,N
      DO 20 J=1,N
      POM=0.0E0
      DO 10 K=1,N
10    POM=POM+U(K,I)*A(K,J)
20    W(J,I)=POM
C
      DO 40 I=1,N
      DO 40 J=1,N
      POM=0.0E0
      DO 30 K=1,N
30    POM=POM+W(K,I)*U(K,J)
40    B(I,J)=POM
C
      RETURN
      END
