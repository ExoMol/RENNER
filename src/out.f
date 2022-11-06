C***********************************************
      SUBROUTINE OUTINT(NRAUS,NFIL53,NFIL6,SYMM)
C***********************************************
      INTEGER NRAUS,NFIL53,NFIL6
      LOGICAL SYMM
C
      INTEGER II,III
      REAL*8 XLAB1,XLAB2,E1,DE,LINSTR,INTENS
      INTEGER JX21,N1,K1,V21,NS1,ISURF1,
     1        INDSY1
      INTEGER JX22,N2,K2,V22,NS2,ISURF2,
     1        INDSY2
      CHARACTER*1 CABU1,CABU2
      CHARACTER*2 TOTSY1,TOTSY2
      INTEGER ZZ
C
      include 'intens.h'
C
      IF (IOPT.EQ.1) THEN
        WRITE (NFIL6,2000)
        ZZ=1
        DO 10 II=1,NRAUS
          READ (NFIL53,REC=II)
     1         XLAB1,XLAB2,E1,DE,LINSTR,INTENS
          CALL UNPACK(JX21,N1,K1,V21,NS1,ISURF1,
     1                CABU1,INDSY1,XLAB1)
          IF (SYMM) THEN
            IF (INDSY1.EQ.1) TOTSY1='A1'
            IF (INDSY1.EQ.2) TOTSY1='B1'
            IF (INDSY1.EQ.3) TOTSY1='B2'
            IF (INDSY1.EQ.4) TOTSY1='A2'
          ELSE
            IF (INDSY1.EQ.1) TOTSY1='A'''
            IF (INDSY1.EQ.2) TOTSY1='A"'
          ENDIF
          CALL UNPACK(JX22,N2,K2,V22,NS2,ISURF2,
     1                CABU2,INDSY2,XLAB2)
          IF (SYMM) THEN
            IF (INDSY2.EQ.1) TOTSY2='A1'
            IF (INDSY2.EQ.2) TOTSY2='B1'
            IF (INDSY2.EQ.3) TOTSY2='B2'
            IF (INDSY2.EQ.4) TOTSY2='A2'
          ELSE
            IF (INDSY2.EQ.1) TOTSY2='A'''
            IF (INDSY2.EQ.2) TOTSY2='A"' 
          ENDIF
          IF (INTENS.GT.TRTLIM) THEN
            IF (MOD(JX21,2).EQ.0) THEN
              WRITE (NFIL6,2001)
     1              ZZ,
     2              JX22/2,TOTSY2,N2,
     1              K2,V22,NS2,ISURF2,CABU2,
     2              JX21/2,TOTSY1,N1,
     1              K1,V21,NS1,ISURF1,CABU1,
     3              E1,DE,LINSTR,INTENS
            ELSE
              WRITE (NFIL6,2002)
     1              ZZ,
     2              JX22,TOTSY2,N2,
     2              K2,V22,NS2,ISURF2,CABU2,
     2              JX21,TOTSY1,N1,
     2              K1,V21,NS1,ISURF1,CABU1,
     3              E1,DE,LINSTR,INTENS
            ENDIF
            ZZ=ZZ+1
          ENDIF
10      CONTINUE
        WRITE (NFIL6,2003)
      ENDIF
C
      IF (IOPT.EQ.2) THEN
        III=1
        WRITE (NFIL6,1000)
        DO 20 II=1,NRAUS
          READ (NFIL53,REC=II)
     1         XLAB1,XLAB2,E1,DE,LINSTR,INTENS
          IF (INTENS.GT.TRTLIM) THEN
            WRITE (NFIL6,1001) DE-1.0D-6,0.0 
            WRITE (NFIL6,1001) DE,INTENS
            WRITE (NFIL6,1001) DE+1.0D-6,0.0
            III=III+1
          ENDIF
20      CONTINUE
        WRITE (NFIL6,1002)
      ENDIF
C
1000  FORMAT('0',5X,49('-'),/,
     1       '0',5X,'I',47X,'I',/,       
     1       '0',5X,'I    ','  ',
     1              'ENERGY [cm-1]   ','  ',
     1              'INTENSITY [cm*mole-1]  I',/,
     1       '0',5X,'I',47X,'I',/,
     1       '0',5X,49('-'))
1001  FORMAT(5X,'  ',D16.8,'  ',D16.8)
1002  FORMAT('0',5X,49('-'))
C
2000  FORMAT('0',5X,121('-'),/,
     1       '0',5X,'I',119X,'I',/,
     1       '0',5X,'I','          UPPER                    '
     2                 ,'     LOWER                    '
     3                 ,54X, 'I',/,
     1       '0',5X,'I',119X,'I',/,
     2       '0',5X,'I  ','   No ','  J','   ','SYM',
     3              '  N KA V2 NS',1X,1X,1X,1X,'I',1X,
     4              '   ','  J','   ','SYM',
     5              '  N KA V2 NS',1X,1X,1X,1X,'I',1X,
     6              ' E(LOWER)   ','I',1X,
     6              ' TRANSITION ','I',1X,
     7              ' LINESTRGTH ','I',1X,
     8              ' INTENSITY  ','I',/,
     1       '0',5X,'I',33X,'I',29X,'I',13X,'I',13X,
     1              'I',13X,'I',13X,'I',/,
     1       '0',5X,121('-'))
2001  FORMAT('0',5X,'I  ',I5,' ',
     1              I3,'   ',A2,' ',
     1              4I3,1X,I1,A1,1X,'I',1X,
     1              '   ',I3,'   ',A2,' ',
     2              4I3,1X,I1,A1,1X,'I',1X,
     3              F11.4,1X,'I',1X,
     3              F11.4,1X,'I',1X,
     4              E11.4,1X,'I',1X,
     5              E11.4,1X,'I')
2002  FORMAT('0',5X,'I  ',I5,' ',
     1              I3,'/2 ',A2,' ',
     1              4I3,1X,I1,A1,1X,'I',1X,
     1              '   ',I3,'/2 ',A2,' ',
     2              4I3,1X,I1,A1,1X,'I',1X,
     3              F11.4,1X,'I',1X,
     3              F11.4,1X,'I',1X,
     4              E11.4,1X,'I',1X,
     5              E11.4,1X,'I')
2003  FORMAT('0',5X,121('-'))
C
      RETURN
      END
