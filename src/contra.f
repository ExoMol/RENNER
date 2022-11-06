C**************************************************************
      SUBROUTINE RECPOS(JMAXP1,JP1,KP1,INDSYM,DELK,NMINP1,RECN)
C**************************************************************
C     Diese Prozedur berechnet die Position eines
C     zusammengesetzten K-Blocks
C     in der Datei NFIL34.
      INTEGER JMAXP1,JP1,KP1,NSURF,INDSYM,DELK,NMINP1,RECN
      RECN= ( JP1*(JP1-1)/2 - NMINP1*(NMINP1-1)/2 + (KP1-1) ) * 20
     1      +(5*(INDSYM-1)+DELK+2)+1
C     WRITE (*,100) 'JMAXP1=',JMAXP1,'JP1=',JP1,'KP1=',KP1,
C    1              'INDSYM=',INDSYM,
C    2              'DELK=',DELK,'RECN=',RECN
100   FORMAT(1X,7(A,I3))
      RETURN
      END

C*********************************************************
      SUBROUTINE CONTRA(BIG,EIGVAL,IFAIL,LREC34,
     1                  LREC35,DIMTOT,NABLK,NBBLK,
     1                  WRK,IWORK,ASSGN,MAXPSY,NFCTS,
     2                  IFLAG,ITAUP1,KP1,
     3                  COENMX,ISOMAX,ISO,
     5                  PA,PB,ESAVE,PABSV,
     5                  EIGDAT,
     6                  IBOOK,MXJVP1,NMINP1,
     7                  POSFCT,MXFCCO,NRFCCO,
     8                  SYMM)
C*********************************************************
C
      implicit real*8 (A-H,O-Z)
C
      include 'rensys.h'
      include 'dimen.h'
C
      INTEGER PRTINT,V2MAX,V2MXP1
C

      INTEGER        NSURF,  MBASIS ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ , NFINTA
C
      INTEGER DIMTOT
      LOGICAL SYMM
      REAL*8 COENMX(ISOMAX)
      REAL*8 BIG(LREC35,LREC35)
      REAL*8 EIGVAL(DIMTOT)
      INTEGER IFAIL(DIMTOT),INFO
      REAL*8 ASSGN(JMAXP1,2*MAXPSY,LREC35,5)
      INTEGER NFCTS(JMAXP1,4)
      INTEGER POSFCT(JMAXP1,4)
      REAL*8 WRK(DIMTOT*8)
      INTEGER IWORK(DIMTOT*5)
      INTEGER DELK
      INTEGER INDSYM,KP1
      CHARACTER*1 JOBZ,RANGE,UPLO
      CHARACTER*130 FLINE
      REAL*8 VL,VU
      INTEGER IL,IU
      REAL*8 ABSTOL
      INTEGER NF,II,JJ
      INTEGER M
      INTEGER IRECE,IR,JL
      INTEGER IVEC,MAXC
      INTEGER NDIMST,V2,NS
      INTEGER IFLAG,ITAUP1,STRSYM
      INTEGER NAGTB,NBGTA
      REAL*8 QA,QB
      INTEGER JLINE,KVAL,KAPART,KBPART
      REAL*8 PA(DIMTOT),PB(DIMTOT)
      REAL*8 ESAVE(DIMTOT,2),PABSV(DIMTOT,2)
      INTEGER EIGDAT(DIMTOT,2,3)
      INTEGER IBOOK(MXJVP1,2,2,12)
C
      CHARACTER   STSYM
      CHARACTER*2 TOTSYM
C
      INDSYM=2*(IFLAG-1)+ITAUP1
C
C     Diagonalisierung der Matrix BIG
C     
      JOBZ='V'
      RANGE='A'
      UPLO='U'
      VL=0.0D0
      VU=0.0D0
      IL=1
      IU=1
      ABSTOL=0.0D0
      M=1
      CALL DSYEVX(JOBZ,RANGE,UPLO,
     1            DIMTOT ,BIG,LREC35,
     2            VL,VU,IL,IU,
     3            ABSTOL,
     4            M,EIGVAL,BIG,LREC35,
     5            WRK,DIMTOT*8,IWORK,
     6            IFAIL,INFO)
      IF (INFO.NE.0) THEN
        WRITE (NFIL6,6100)
        STOP
      ENDIF
C
C     Analyse der Eigenvektoren
C
C     Bedeutung der Feldelemente von ASSGN:
C
C     <1> 1 oder 2, je nachdem, ob a- oder b-Charakter ueberwiegend
C     die anderen charakterisieren diejenige Basisfunktion
C     mit dem groessten Beitrag:
C     <2> NS
C     <3> V2
C     <4> 1 oder 2, je nachdem, ob A- oder B-Streckschwingungsfunktion
C     <5> E
C
      DO 40 IVEC=1,DIMTOT
         ASSGN(KP1,INDSYM,IVEC,5)=EIGVAL(IVEC)
         MAXC=IDAMAX(DIMTOT,BIG(1,IVEC),1)
         IF (MAXC.LE.NABLK) THEN
           ASSGN(KP1,INDSYM,IVEC,1)=1
           IF (MOD(KP1,2).EQ.1) THEN
             NDIMST=IBOOK(KP1,ITAUP1,IFLAG,7)
           ELSE
             NDIMST=IBOOK(KP1,ITAUP1,IFLAG,8)
           ENDIF
         ELSE
           ASSGN(KP1,INDSYM,IVEC,1)=2
           MAXC=MAXC-NABLK
           IF (MOD(KP1,2).EQ.1) THEN
             NDIMST=IBOOK(KP1,ITAUP1,IFLAG,9)
           ELSE
             NDIMST=IBOOK(KP1,ITAUP1,IFLAG,10)
           ENDIF
         ENDIF
         X=DFLOAT(MAXC-1)/NDIMST
         V2=INT(X)
         NS=MAXC-NDIMST*V2
         ASSGN(KP1,INDSYM,IVEC,3)=V2
         ASSGN(KP1,INDSYM,IVEC,2)=NS
         STRSYM=1
         ITAU=ITAUP1-1
         IF (IFLAG.EQ.1) THEN
            IF (MOD(KP1,2).EQ.0) STRSYM=2
         ELSE
            IF (MOD(KP1,2).EQ.1) STRSYM=2
         ENDIF
         IF (.NOT.SYMM) STRSYM=1
         ASSGN(KP1,INDSYM,IVEC,4)=STRSYM
C
 40   CONTINUE
C
         NF=0
         IF (COENMX(ISO).EQ.0.0) THEN
            NF=DIMTOT
         ELSE
           DO 52 II=1,DIMTOT
             IF (EIGVAL(II).LE.COENMX(ISO)) NF=II
 52        CONTINUE
         ENDIF
         NFCTS(KP1,INDSYM)=NF
C
C     Abspeichern der Eigenvektoren
C 
      IF (NF.GT.MXFCCO) MXFCCO=NF
      POSFCT(KP1,INDSYM)=NRFCCO
      DO 54 II=1,NF
        WRITE (NFIL35,REC=NRFCCO+II-1)
     1        (BIG(JJ,II),JJ=1,LREC35)
54    CONTINUE 
      NRFCCO=NRFCCO+NF     
C
C     Verschiedene Informationen ausdrucken
C
       NAGTB=0
       DO 70 II=1,DIMTOT
          QA=0.0D+00
          QB=0.0D+00
C
          DO 50 JJ=1,NABLK
             QA=QA+BIG(JJ,II)**2
 50       CONTINUE
C
          DO 60 JJ=1,NBBLK
             QB=QB+BIG(JJ+NABLK,II)**2
 60       CONTINUE
C
          IF (QA .GE. QB) NAGTB=NAGTB+1
C
          PA(II)=QA
          PB(II)=QB
 70    CONTINUE
C
       NBGTA=DIMTOT-NAGTB
C
       JJ=0
       DO 110 II=1,NAGTB
100    JJ=JJ+1
       IF (PA(JJ) .GE. PB(JJ)) THEN
           ESAVE(II,1)=EIGVAL(JJ)
	   PABSV(II,1)=PA(JJ)
	   EIGDAT(II,1,1)=NINT(ASSGN(KP1,INDSYM,JJ,1))
	   EIGDAT(II,1,2)=NINT(ASSGN(KP1,INDSYM,JJ,2))
	   EIGDAT(II,1,3)=NINT(ASSGN(KP1,INDSYM,JJ,3))
       ELSE
	   GOTO 100
       ENDIF
110    CONTINUE
C
       JJ=0
       DO 130 II=1,NBGTA
120    JJ=JJ+1
       IF (PA(JJ) .LT. PB(JJ)) THEN
	   ESAVE(II,2)=EIGVAL(JJ)
	   PABSV(II,2)=PA(JJ)
	   EIGDAT(II,2,1)=NINT(ASSGN(KP1,INDSYM,JJ,1))
	   EIGDAT(II,2,2)=NINT(ASSGN(KP1,INDSYM,JJ,2))
	   EIGDAT(II,2,3)=NINT(ASSGN(KP1,INDSYM,JJ,3))
       ELSE
	   GOTO 120
       ENDIF
130    CONTINUE
C
      IF (IPRINT.NE.0) THEN
        KVAL=KP1-1
C
        IF (MOD(MOD(KVAL,2)+(IFLAG-1),2).EQ.0) THEN
          STSYM='A'
        ELSE
          STSYM='B'
        ENDIF
        IF (INDSYM.EQ.1) TOTSYM='A1'
        IF (INDSYM.EQ.3) TOTSYM='B2'
        IF (INDSYM.EQ.2) TOTSYM='B1'
        IF (INDSYM.EQ.4) TOTSYM='A2'
C
        WRITE (NFIL6,1500) TOTSYM,KVAL,NAGTB,NBGTA,STSYM
        DO 132 JLINE=1,MAX(NAGTB,NBGTA)
          WRITE (FLINE(1:12),1600) JLINE
          IF (JLINE .GT. NAGTB) THEN
	      WRITE (FLINE(13:56),1700)
          ELSE
              KAPART=NINT(1.0D+02*PABSV(JLINE,1))
              KBPART=NINT(1.0D+02*(1.0D+00-PABSV(JLINE,1)))
	      WRITE (FLINE(13:56),1800) ESAVE(JLINE,1),KAPART,KBPART,
     1	            EIGDAT(JLINE,1,3),EIGDAT(JLINE,1,2),EIGDAT(JLINE,1,1)
          ENDIF
C
          IF (JLINE .GT. NBGTA) THEN
	      WRITE (FLINE(57:100),1700)
          ELSE
              KAPART=NINT(1.0D+02*PABSV(JLINE,2))
              KBPART=NINT(1.0D+02*(1.0D+00-PABSV(JLINE,2)))
	      WRITE (FLINE(57:100),1800) ESAVE(JLINE,2),KAPART,KBPART,
     1  	    EIGDAT(JLINE,2,3),EIGDAT(JLINE,2,2),EIGDAT(JLINE,2,1)
          ENDIF
          WRITE (NFIL6,1900) FLINE(1:100)
132     CONTINUE
        WRITE (NFIL6,2000)
      ENDIF
140   CONTINUE
C
      RETURN
1500  FORMAT('1',5X,12('*'),' ', A,
     1       ' CONTRACTED FUNCTIONS FOR K = ',I3,1X,12('*')//
     3 ' ',17X,'NUMBER OF A-TYPE FUNCTIONS                 ',I3/
     4 ' ',17X,'NUMBER OF B-TYPE FUNCTIONS                 ',I3/
     4 ' ',17X,'STRETCHING SYMMETRY                          ',A //
     1 ' ',10X,100('-')/
     2 ' ',10X,'I',10X,'I',43X,'I',43X,'I'/
     2 ' ',10X,'I',2X,'FCT. #',2X,'I',
     3  2(4X,'ENERGY/CM-1 ',2X,'%A ',2X,'%B ',3X,
     4       'V2 ',2X,'NS ',2X,'SUR',1X,'I')/
     1 ' ',10X,'I',10X,'I',43X,'I',43X,'I'/' ',10X,100('-')/
     1 ' ',10X,'I',43X,'I',43X,'I')
1600  FORMAT('I',3X,I4,3X,'I')
1700  FORMAT(43X,'I')
1800  FORMAT(2X,F14.5,2X,I3,2X,I3,2X,I3,2X,I3,2X,I3,2X,'I')
1900  FORMAT(' ',10X,A100)
2000  FORMAT(' ',10X,100('-')///)
6100  FORMAT(1H0,'  RENTEL.CTK.ERR  MATRIX DIAGONALIZATION  ',
     1             'FAILED IN CONTRA'/)
6200  FORMAT(1H0,'  RENTEL.CTK.WRN  FOR K = ',I3,' ONLY ',I3,
     1             ' A-TYPE BENDING FUNCTIONS WERE FOUND.',/,
     2       1H ,15X,'   V2CA = ',I3,' WILL BE LOWERED ACCORDINGLY')
6300  FORMAT(1H0,'  RENTEL.CTK.WRN  FOR K = ',I3,' ONLY ',I3,
     1             ' B-TYPE BENDING FUNCTIONS WERE FOUND.',/,
     2       1H ,15X,'   V2CB = ',I3,' WILL BE LOWERED ACCORDINGLY')
       END
C
C*********************************************************
      SUBROUTINE MATPOS(NMIN,N,K,NFCTS,JMAXP1,INDSYM,MPOS)
C*********************************************************
C     Diese Prozedur dient der korrekten Positionierung
C     eines Blockes in der Endmatrix.
      INTEGER NMIN,N,K,INDSYM
      INTEGER NFCTS(JMAXP1,4)
      INTEGER MPOS
      INTEGER NZ,I
      INTEGER NUP,N0
C
      MPOS=0
      NP1=N+1
      NMINP1=NMIN+1
      NZ=N-NMIN
      IF (K.EQ.0) THEN
        NUP=N-1
      ELSE
        NUP=N
      ENDIF
      DO 30 N0=NMIN,NUP,1
        IF (MOD(N0,2).EQ.1) THEN
          IF (MOD(INDSYM,2).EQ.0) THEN 
            INDSYH=INDSYM-1
          ELSE
            INDSYH=INDSYM+1
          ENDIF
        ELSE
          INDSYH=INDSYM
        ENDIF
        MPOS=MPOS+NFCTS(1,INDSYH)
30    CONTINUE
      DO 10 I=2,NP1-1
        IF (I.GT.NMINP1) NZ=NZ-1
        MPOS=MPOS+NZ*NFCTS(I,INDSYM)
10    CONTINUE
      DO 20 I=2,K
        MPOS=MPOS+NFCTS(I,INDSYM)
20    CONTINUE
      RETURN
      END
C
C******************************************************
      SUBROUTINE MAKMAT(BTRN,
     1                  HMAT,KMAT,WMAT,KBLK,
     2                  IRODIM,LREC34,LREC35,NMIN,NMAX,
     3                  IABLKE,IBBLKE,IABLKO,IBBLKO,
     4                  ITAUP1,INDSYM,NFCTS,
     5                  ASMAPP,ASMSPP,AASMPP,
     5                  ASMAMM,ASMSMM,AASMMM,
     5                  ASMAMP,ASMSMP,AASMMP,
     6                  LSPOS1,LSPOA1,
     6                  LSPOS2,LSPOA2,
     8                  ISPOSY,ISPOAS,
     9                  NOBAL,NOBAR,
     1                  NOPERS,NINTOP,
     2                  FACTS,JX2VAL,N2XS,ESS,XJVAL,
     3                  IBOOK,MXJVP1,IPSYM,
     4                  KSTY11,KSTY12,
     5                  NSPOSY,NSPOAS,SYMM,
     6                  POSFCT,MXFCCO,NRFCCO)
C******************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'rensys.h'
      include 'dimen.h'
C
      INTEGER PRTINT,V2MAX,V2MXP1
C
C
      INTEGER POSFCT(JMAXP1,4)
      REAL*8 BTRN(LREC35,MXFCCO)
      REAL*8 HMAT(IRODIM,IRODIM)
      REAL*8 KMAT(LREC35,LREC35),WMAT(LREC35,MXFCCO)
      REAL*8 KBLK(LREC35,LREC35)
      INTEGER NFCTS(JMAXP1,4)
      LOGICAL LEER
      INTEGER DIMR,DIMC
C
      INTEGER LSPOS1,LSPOA1,
     6        LSPOS2,LSPOA2,
     8        ISPOSY,ISPOAS,
     9        NOBAL,NOBAR,
     1        NOPERS,NINTOP
      REAL*8  ASMAPP(NOBAL,NOBAL,NOPERS),
     1        ASMSPP(LSPOS1),AASMPP(LSPOA1)
      REAL*8  ASMAMM(NOBAR,NOBAR,NOPERS),
     1        ASMSMM(LSPOS2),AASMMM(LSPOA2)
      REAL*8  ASMAMP(NOBAL,NOBAR,NINTOP),
     1        ASMSMP(ISPOSY),AASMMP(ISPOAS)
C
      LOGICAL SYMM
      INTEGER IBOOK(MXJVP1,2,2,12)
C
      INTEGER V2MP11,V2MP12
C
      NMINP1=NMIN+1
      NMAXP1=NMAX+1
      DO 10 N1P1=NMINP1,NMAXP1
      DO 20 NDEL=0,1,1
      N2P1=N1P1+NDEL
      IF ((N2P1.LT.NMINP1).OR.(N2P1.GT.NMAXP1)) GOTO 20
      NVAL=N1P1-1
      XNVAL=DFLOAT(NVAL)
      KNDSPP=0
      KNDAPP=0
      KNDSMM=0
      KNDAMM=0
      KNDSMP=0
      KNDAMP=0
      KNDSPM=0
      KNDAPM=0
      DO 30 K1P1=1,N1P1
      DO 40 KDEL=0,2,1
      K2P1=K1P1+KDEL
      IF ((K2P1.LT.1).OR.(K2P1.GT.N2P1)) GOTO 40
      KP1=K1P1
C
C     Die Dimensionen des Blockes bestimmen
C
      IF (MOD(K1P1,2).EQ.1) THEN
        IABLKR=IABLKE
        IBBLKR=IBBLKE
      ELSE
        IABLKR=IABLKO
        IBBLKR=IBBLKO
      ENDIF
      IF (MOD(K2P1,2).EQ.1) THEN
        IABLKC=IABLKE
        IBBLKC=IBBLKE
      ELSE
        IABLKC=IABLKO
        IBBLKC=IBBLKO
      ENDIF
      IF (K1P1.EQ.1) THEN
        IF (MOD(N1P1+ITAUP1,2).EQ.0) THEN
          IBBLKR=0
        ELSE
          IABLKR=0
        ENDIF
      ENDIF
      IF (K2P1.EQ.1) THEN
        IF (MOD(N2P1+ITAUP1,2).EQ.0) THEN
          IBBLKC=0
        ELSE
          IABLKC=0
        ENDIF
      ENDIF
      DIMR=IABLKR+IBBLKR
      DIMC=IABLKC+IBBLKC
C
C     Block loeschen
C
      DO 50 I=1,LREC35
      DO 50 J=1,LREC35
        KMAT(I,J)=0
50    CONTINUE
      LEER=.TRUE.
C
C     additive SO-Kopplung
C
      IF ((NDEL.EQ.0).AND.(KDEL.EQ.0)) THEN
      LEER=.FALSE.
C
C   DELTA-N = 0 MATRIX ELEMENTS
C
      NABLK=IABLKR
      NBBLK=IBBLKR
C
C      DO ALL MATRIX ELEMENTS COUPLING + WITH +
C
      FACTJ=FACTS*SQRT(DFLOAT((2*NVAL+1)**2))
     1           *((-1.0D+00)**((JX2VAL+N2XS)/2))
     1           *SIXJ( XNVAL , ESS   , XJVAL   ,
     1                  ESS   , XNVAL , 1.0D+00 )
C
      NOFSRO=0
      NOFSCO=0
      V2MP11=IBOOK(NVAL+1,ITAUP1,IPSYM,11)
      V2MP12=IBOOK(NVAL+1,ITAUP1,IPSYM,11)
      ITAURO=ITAUP1-1
      ITAUCO=ITAUP1-1
      NP1ROW=NVAL+1
      NP1COL=NVAL+1
      IF (MOD(KP1,2) .EQ. 1) THEN
      NDIMRO=IBOOK(NVAL+1,ITAUP1,IPSYM,7)
      NDIMCO=IBOOK(NVAL+1,ITAUP1,IPSYM,7)
      IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
      ELSE
          ISTOFF=KSTY11
          JSTOFF=KSTY11
      ENDIF
      ELSE
      NDIMRO=IBOOK(NVAL+1,ITAUP1,IPSYM,8)
      NDIMCO=IBOOK(NVAL+1,ITAUP1,IPSYM,8)
      IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
      ELSE
          ISTOFF=KSTY11
          JSTOFF=KSTY11
      ENDIF
      ENDIF
C
      IF (NABLK .GT. 0) THEN
      CALL SOME  ( ASMSPP , AASMPP , ASMAPP ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             LSPOS1 , LSPOA1 , ISTOFF , JSTOFF ,
     1             NOBAL  , NOBAL  , NOPERS , NSPOSY ,
     1             NSPOAS , FACTJ  , 0 , KP1 , KNDSPP , KNDAPP )
      ELSE 
           KNDSPP=KNDSPP+V2MP11*V2MP12
           IF (.NOT.SYMM) KNDAPP=KNDAPP+V2MP11*V2MP12
      ENDIF
C
C      DO ALL MATRIX ELEMENTS COUPLING - WITH -
C
      NOFSRO=NABLK
      NOFSCO=NABLK
      V2MP11=IBOOK(NVAL+1,ITAUP1,IPSYM,12)
      V2MP12=IBOOK(NVAL+1,ITAUP1,IPSYM,12)
      ITAURO=2-ITAUP1
      ITAUCO=2-ITAUP1
      NP1ROW=NVAL+1
      NP1COL=NVAL+1
      IF (MOD(KP1,2) .EQ. 1) THEN
      NDIMRO=IBOOK(NVAL+1,ITAUP1,IPSYM,9)
      NDIMCO=IBOOK(NVAL+1,ITAUP1,IPSYM,9)
      IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
      ELSE
          ISTOFF=KSTY12
          JSTOFF=KSTY12
      ENDIF
      ELSE
      NDIMRO=IBOOK(NVAL+1,ITAUP1,IPSYM,10)
      NDIMCO=IBOOK(NVAL+1,ITAUP1,IPSYM,10)
      IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
      ELSE
          ISTOFF=KSTY12
          JSTOFF=KSTY12
      ENDIF
      ENDIF
      IF (NBBLK .GT. 0) THEN
      CALL SOME  ( ASMSMM , AASMMM , ASMAMM ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             LSPOS2 , LSPOA2 , ISTOFF , JSTOFF ,
     1             NOBAR  , NOBAR  , NOPERS , NSPOSY ,
     1             NSPOAS , FACTJ  , 0 , KP1 , KNDSMM , KNDAMM )
      ELSE
           KNDSMM=KNDSMM+V2MP11*V2MP12
           IF (.NOT.SYMM) KNDAMM=KNDAMM+V2MP11*V2MP12
      ENDIF
C
C      DO ALL MATRIX ELEMENTS COUPLING + WITH -
C
      NOFSRO=0
      NOFSCO=NABLK
      V2MP11=IBOOK(NVAL+1,ITAUP1,IPSYM,11)
      V2MP12=IBOOK(NVAL+1,ITAUP1,IPSYM,12)
      ITAURO=ITAUP1-1
      ITAUCO=2-ITAUP1
      NP1ROW=NVAL+1
      NP1COL=NVAL+1
C
      IF (MOD(KP1,2) .EQ. 1) THEN
      NDIMRO=IBOOK(NVAL+1,ITAUP1,IPSYM,7)
      NDIMCO=IBOOK(NVAL+1,ITAUP1,IPSYM,9)
      IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
      ELSE
          ISTOFF=KSTY11
          JSTOFF=KSTY12
      ENDIF
      ELSE
      NDIMRO=IBOOK(NVAL+1,ITAUP1,IPSYM,8)
      NDIMCO=IBOOK(NVAL+1,ITAUP1,IPSYM,10)
      IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
      ELSE
          ISTOFF=KSTY11
          JSTOFF=KSTY12
      ENDIF
      ENDIF
      IF (NABLK*NBBLK .GT. 0) THEN
      CALL SOME  ( ASMSMP , AASMMP , ASMAMP ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             ISPOSY , ISPOAS , ISTOFF , JSTOFF ,
     1             NOBAL  , NOBAR  , NINTOP , NSPOSY ,
     1             NSPOAS , FACTJ  , 1 , KP1 , KNDSMP , KNDAMP )
      DO 200 II=1,NBBLK
      DO 200 JJ=1,NABLK
        KMAT(NABLK+II,JJ)=KMAT(JJ,NABLK+II)
200   CONTINUE     
      ELSE 
           KNDSMP=KNDSMP+V2MP11*V2MP12
           IF (.NOT.SYMM) KNDAMP=KNDAMP+V2MP11*V2MP12
      ENDIF
      ENDIF
C
C     reine SO-Kopplung
C
      IF ((NDEL.EQ.1).AND.(KDEL.EQ.0)) THEN
        LEER=.FALSE.
C
C   DELTA-N = 1 MATRIX ELEMENTS
C
      XNVP1=DFLOAT(NVAL+1)
C
C      DO ALL MATRIX ELEMENTS COUPLING + WITH +
C
      FACTJ=FACTS*SQRT(DFLOAT((2*NVAL+1)*(2*NVAL+3)))
     1           *((-1.0D+00)**(1+(JX2VAL+N2XS)/2))
     1           *SIXJ( XNVAL , ESS   , XJVAL   ,
     1                  ESS   , XNVP1 , 1.0D+00 )
C
      NOFSRO=0
      NOFSCO=0
      V2MP11=IBOOK(N1P1,ITAUP1,IPSYM,11)
      V2MP12=IBOOK(N2P1,ITAUP1,IPSYM,11)
      ITAURO=ITAUP1-1
      ITAUCO=ITAUP1-1
      NP1ROW=N1P1
      NP1COL=N2P1
      IF (MOD(KP1,2) .EQ. 1) THEN
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,7)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,7)
        IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
            ISTOFF=0
            JSTOFF=0
        ELSE
            ISTOFF=KSTY11
            JSTOFF=KSTY11
        ENDIF
      ELSE
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,8)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,8)
        IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
            ISTOFF=0
            JSTOFF=0
        ELSE
            ISTOFF=KSTY11
            JSTOFF=KSTY11
        ENDIF
      ENDIF
C
      IF ((IABLKC*IABLKR) .GT. 0) THEN
      CALL SOME  ( ASMSPP , AASMPP , ASMAPP ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             LSPOS1 , LSPOA1 , ISTOFF , JSTOFF ,
     1             NOBAL  , NOBAL  , NOPERS , NSPOSY ,
     1             NSPOAS , FACTJ  , 0 , KP1 , KNDSPP , KNDAPP )
      ELSE
        KNDSPP=KNDSPP+V2MP11*V2MP12
        IF (.NOT.SYMM) KNDAPP=KNDAPP+V2MP11*V2MP12
      ENDIF
C
C      DO ALL MATRIX ELEMENTS COUPLING - WITH -
C
      NOFSRO=IABLKR
      NOFSCO=IABLKC
      V2MP11=IBOOK(N1P1,ITAUP1,IPSYM,12)
      V2MP12=IBOOK(N2P1,ITAUP1,IPSYM,12)
      ITAURO=2-ITAUP1
      ITAUCO=2-ITAUP1
      NP1ROW=N1P1
      NP1COL=N2P1
      IF (MOD(KP1,2) .EQ. 1) THEN
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,9)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,9)
        IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
            ISTOFF=0
            JSTOFF=0
        ELSE
            ISTOFF=KSTY12
            JSTOFF=KSTY12
        ENDIF
      ELSE
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,10)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,10)
        IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
        ELSE
          ISTOFF=KSTY12
          JSTOFF=KSTY12
        ENDIF
      ENDIF
      IF ((IBBLKR*IBBLKC) .GT. 0) THEN
      CALL SOME  ( ASMSMM , AASMMM , ASMAMM ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             LSPOS2 , LSPOA2 , ISTOFF , JSTOFF ,
     1             NOBAR  , NOBAR  , NOPERS , NSPOSY ,
     1             NSPOAS , FACTJ  , 0 , KP1 , KNDSMM , KNDAMM )
      ELSE
        KNDSMM=KNDSMM+V2MP11*V2MP12
        IF (.NOT.SYMM) KNDAMM=KNDAMM+V2MP11*V2MP12    
      ENDIF
C
C      DO ALL MATRIX ELEMENTS COUPLING + WITH -
C
      NOFSRO=0
      NOFSCO=IABLKC
      V2MP11=IBOOK(N1P1,ITAUP1,IPSYM,11)
      V2MP12=IBOOK(N2P1,ITAUP1,IPSYM,12)
      ITAURO=ITAUP1-1
      ITAUCO=2-ITAUP1
      NP1ROW=N1P1
      NP1COL=N2P1
C
      IF (MOD(KP1,2) .EQ. 1) THEN
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,7)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,9)
        IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
            ISTOFF=0
            JSTOFF=0
        ELSE
            ISTOFF=KSTY11
            JSTOFF=KSTY12
        ENDIF
      ELSE
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,8)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,10)
        IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
        ELSE
          ISTOFF=KSTY11
          JSTOFF=KSTY12
        ENDIF
      ENDIF
      IF ((IABLKR*IBBLKC) .GT. 0) THEN
      CALL SOME  ( ASMSMP , AASMMP , ASMAMP ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             ISPOSY , ISPOAS , ISTOFF , JSTOFF ,
     1             NOBAL  , NOBAR  , NINTOP , NSPOSY ,
     1             NSPOAS , FACTJ  , 1 , KP1 , KNDSPM , KNDAPM )
      ELSE
        KNDSPM=KNDSPM+V2MP11*V2MP12
        IF (.NOT.SYMM) KNDAPM=KNDAPM+V2MP11*V2MP12
      ENDIF
C
C      DO ALL MATRIX ELEMENTS COUPLING - WITH +
C
      NOFSRO=IABLKR
      NOFSCO=0
      V2MP11=IBOOK(N1P1,ITAUP1,IPSYM,12)
      V2MP12=IBOOK(N2P1,ITAUP1,IPSYM,11)
      ITAURO=2-ITAUP1
      ITAUCO=ITAUP1-1
      NP1ROW=N1P1
      NP1COL=N2P1
C
      IF (MOD(KP1,2) .EQ. 1) THEN
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,9)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,7)
        IF (IPSYM .EQ. 1 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
        ELSE
          ISTOFF=KSTY12
          JSTOFF=KSTY11
        ENDIF
      ELSE
        NDIMRO=IBOOK(N1P1,ITAUP1,IPSYM,10)
        NDIMCO=IBOOK(N2P1,ITAUP1,IPSYM,8)
        IF (IPSYM .EQ. 2 .OR. (.NOT.SYMM)) THEN
          ISTOFF=0
          JSTOFF=0
        ELSE
          ISTOFF=KSTY12
          JSTOFF=KSTY11
        ENDIF
      ENDIF
      IF ((IBBLKR*IABLKC) .GT. 0) THEN
      CALL SOME  ( ASMSMP , AASMMP , ASMAMP ,
     1             KMAT , NOFSRO , NOFSCO ,
     1             NDIMRO , NDIMCO , V2MP11 , V2MP12 ,
     1             LREC35 ,
     1             ITAURO , ITAUCO , NP1ROW , NP1COL ,
     1             ISPOSY , ISPOAS , ISTOFF , JSTOFF ,
     1             NOBAL  , NOBAR  , NINTOP , NSPOSY ,
     1             NSPOAS , FACTJ  , 2 , KP1 , KNDSMP , KNDAMP )
      ELSE
        KNDSMP=KNDSMP+V2MP11*V2MP12
        IF (.NOT.SYMM) KNDAMP=KNDAMP+V2MP11*V2MP12
      ENDIF
      ENDIF
C
C     Transformieren und Speichern
C
C     IF (.NOT.LEER) THEN
C
C     Transformieren
C
      CALL BSSTRN(N1P1,K1P1,N2P1,K2P1,
     1            INDSYM,INDSYM,
     2            POSFCT,NFCTS,
     3            JMAXP1,
     4            NROW,NCOL,
     4            DIMR,DIMC,
     5            KMAT,BTRN,WMAT,
     6            LREC35,MXFCCO,
     7            NFIL35)
C     ENDIF
C
C     Evtl. Block hinzuaddieren
C
      IF ((NDEL.EQ.0).AND.((NROW*NCOL).GT.0)) THEN
        LEER=.FALSE.
        CALL RECPOS(JMAXP1,N1P1,K1P1,INDSYM,KDEL,NMINP1,IREC)
        READ (NFIL34,REC=IREC)
     1       ((KBLK(II,JJ),II=1,MXFCCO),JJ=1,MXFCCO)
        DO 61 JJ=1,MXFCCO
        DO 61 II=1,MXFCCO
          KMAT(II,JJ) = KMAT(II,JJ) + KBLK(II,JJ)
61      CONTINUE
      ENDIF
C
C     Ausgabe des Blockes
C
C       IF (NDEL.EQ.0) THEN
C         WRITE (*,*) 'INDSYM=',INDSYM,
C    1    ' N1P1=',N1P1,' K1P1=',K1P1,' KDEL=',KDEL
C         DO 99 II=1,NCOL
C           WRITE (*,*) (KMAT(JJ,II),JJ=1,NROW)
C99       CONTINUE
C       ENDIF
C
C     Abspeichern des Blockes
C
      IF (.NOT.LEER) THEN
        IF ((NCOL*NROW).GT.0) THEN
          CALL MATPOS(NMIN,N1P1-1,K1P1-1,NFCTS,JMAXP1,INDSYM,NPOSR)
          CALL MATPOS(NMIN,N2P1-1,K2P1-1,NFCTS,JMAXP1,INDSYM,NPOSC)
          DO 60 J=1,NCOL
          DO 60 I=1,NROW
            HMAT(NPOSR+I,NPOSC+J)=KMAT(I,J)
60        CONTINUE
        ENDIF
      ENDIF
40    CONTINUE
30    CONTINUE
20    CONTINUE
10    CONTINUE
      RETURN
      END
C
      SUBROUTINE BSSTRN(N1P1,K1P1,N2P1,K2P1,
     1                  INDSY1,INDSY2,
     2                  POSFCT,NFCTS,
     3                  JMAXP1,
     4                  NROW,NCOL,
     4                  DIMR,DIMC,
     5                  KMAT,BTRN,WMAT,
     6                  LREC35,MXFCCO,
     7                  NFIL35)
C     
      INTEGER POSFCT(JMAXP1,4)
      INTEGER NFCTS (JMAXP1,4)
C
      INTEGER DIMR,DIMC
C     
      REAL*8  KMAT(LREC35,LREC35)
      REAL*8  BTRN(LREC35,MXFCCO),WMAT(LREC35,MXFCCO)
C
C
C     linke Transformationsmatrix laden
C
        IF ((MOD(N1P1,2).EQ.0).AND.(K1P1.EQ.1)) THEN
          IF (MOD(INDSY1,2).EQ.0) THEN
            INDSYL=INDSY1-1
          ELSE
            INDSYL=INDSY1+1
          ENDIF
        ELSE
          INDSYL=INDSY1
        ENDIF
        IRECL=POSFCT(K1P1,INDSYL)
        NROW=NFCTS(K1P1,INDSYL)
        DO 55 II=1,NROW
           READ(NFIL35,REC=IRECL+II-1)
     1         (BTRN(JJ,II),JJ=1,LREC35)
55      CONTINUE
C
C      erste Multiplikation
C
        IF (NROW.GT.0) THEN
          CALL DGEMM('T','N',DIMC,NROW,DIMR,1.0D+00,KMAT,LREC35,BTRN,
     1          LREC35,0.0D+00,WMAT,LREC35)
        ENDIF
C
C       zweite Transformationsmatrix laden
C
        IF ((MOD(N2P1,2).EQ.0).AND.(K2P1.EQ.1)) THEN
          IF (MOD(INDSY2,2).EQ.0) THEN
            INDSYR=INDSY2-1
          ELSE
            INDSYR=INDSY2+1
          ENDIF
        ELSE
          INDSYR=INDSY2
        ENDIF
        IRECR=POSFCT(K2P1,INDSYR)
        NCOL=NFCTS(K2P1,INDSYR)
        DO 56 II=1,NCOL
           READ(NFIL35,REC=IRECR+II-1)
     1         (BTRN(JJ,II),JJ=1,LREC35)
56      CONTINUE
C
C     zweite Mutliplikation
C
        IF (NCOL.GT.0) THEN
          CALL DGEMM('T','N',NROW,NCOL,DIMC,1.0D+00,WMAT,LREC35,BTRN,
     1          LREC35,0.0D+00,KMAT,LREC35)
        ENDIF
C
      RETURN
      END
