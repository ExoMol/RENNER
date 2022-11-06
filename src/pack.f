c*********************************************************
      SUBROUTINE PACK(J2X,N,K,V2,NS,ISURF,CABU,INDSYM,RES)
c*********************************************************
C 
C     Labels in eine reelle Zahl packen 
C
      INTEGER J2X,N,K,V2,NS,ISURF,INDSYM
      CHARACTER*1 CABU
      INTEGER ICABU
      REAL*8  RES 
C
      IF (CABU.EQ.'A') THEN
         ICABU=0
      ELSE
         ICABU=1
      ENDIF
C
      RES=0
C
      RES=J2X
      RES=RES*100+N
      RES=RES*100+K
      RES=RES*100+V2
      RES=RES*100+NS
      RES=RES*10+ISURF
      RES=RES*10+ICABU
      RES=RES*10+INDSYM 
C
      RETURN
      END
C***********************************************************
      SUBROUTINE UNPACK(J2X,N,K,V2,NS,ISURF,CABU,INDSYM,RES)
C***********************************************************
C
C     verschluesselte Labels auspacken
C
      INTEGER J2X,N,K,V2,NS,ISURF,ICABU,INDSYM 
      REAL*8  RES
      REAL*8  RESH
      CHARACTER*1 CABU
C 
      RESH=RES
C
      RES=RES/10
      INDSYM=NINT((RES-DINT(RES))*10)
      RES=DINT(RES)
C
      RES=RES/10
      ICABU=NINT((RES-DINT(RES))*10)
      RES=DINT(RES)
C
      RES=RES/10
      ISURF=NINT((RES-DINT(RES))*10)
      RES=DINT(RES)
C
      RES=RES/100
      NS=NINT((RES-DINT(RES))*100)
      RES=DINT(RES)
C     
      RES=RES/100
      V2=NINT((RES-DINT(RES))*100)   
      RES=DINT(RES)
C     
      RES=RES/100
      K =NINT((RES-DINT(RES))*100)
      RES=DINT(RES)
C     
      RES=RES/100
      N =NINT((RES-DINT(RES))*100)
      RES=DINT(RES)
C     
      RES=RES/100
      J2X=NINT((RES-DINT(RES))*100)
      RES=DINT(RES)
C
      IF (ICABU.EQ.0) THEN
         CABU='A'
      ELSE
         CABU='B'
      ENDIF
C
      RES=RESH
C
      RETURN
      END
