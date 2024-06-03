      SUBROUTINE ORD2(A,NR,IMAX)                                        
C     ORDERS ELEMENTS IN ARRAY A INCREASING IN SIZE                     
C       REORDERS CORRESPONDING INDICES (NR)                             
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL CONT                                                      
      DIMENSION A(*),NR(*)                                              
      DO 50 I=1,IMAX                                                    
   50 NR(I)=I                                                           
  100 I=1                                                               
      CONT=.FALSE.                                                      
  110 I=I+1                                                             
      IF(A(I).LT.A(I-1))THEN                                            
C       INTERCHANGE I AND (I-1) ELEMENT IN ARRAYS A AND NR              
        HLP=A(I)                                                        
        A(I)=A(I-1)                                                     
        A(I-1)=HLP                                                      
        NHLP=NR(I)                                                      
        NR(I)=NR(I-1)                                                   
        NR(I-1)=NHLP                                                    
        CONT=.TRUE.                                                     
      ENDIF                                                             
      IF(I.LT.IMAX)GO TO 110                                            
      IF(CONT) GO TO 100                                                
      RETURN                                                            
      END                                                               
