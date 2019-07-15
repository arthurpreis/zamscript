        PROGRAM teste
        IMPLICIT DOUBLE PRECISION(E-Z)
        INTEGER I, N, NC, NFOUT
        REAL XI(201), DXI(201)

        WRITE(*,*) 'FORTRAN TESTE'
        N = 201
        NFOUT=0
        NC=185
        EPS1=0.035D0
        EPS1P1=EPS1+1.0D0
        EPS2=-0.15D0
        EPS2P1=EPS2+1.0D0
        N1=N-1
        N2=N-2
        XI(N1)=DLOG(1.D-10)
        TEMP1=(EPS1P1**(NC-1)-1.0D0)/EPS1
        TEMP2=EPS1P1**(NC-2)*EPS2P1*(EPS2P1**(N1-NC)-1.0D0)/EPS2
        WRITE(*,*) TEMP1, TEMP2
        DXI(1)=XI(N1)/(TEMP1+TEMP2)
        WRITE(*,*) 1, XI(1), DXI(1)
        DO 1 I=2,NC-1
            TEMP=EPS1P1**(I-1)
            DXI(I)=DXI(1)*TEMP
            XI(I)=DXI(1)*(TEMP-1.0D0)/EPS1
            WRITE(*,*) I, XI(I), DXI(I)
1       CONTINUE
        XI(NC)=XI(NC-1)+DXI(NC-1)
        DO 2 I=NC,N-3
            J=I-NC+1
            DXI(I)=DXI(NC-1)*EPS2P1**J
            TEMP=EPS2P1*(EPS2P1**J-1.0D0)/EPS2
            XI(I+1)=XI(NC)+DXI(NC-1)*TEMP
            WRITE(*,*) I, XI(I), DXI(I)
    2   CONTINUE
        DXI(N-2)=EPS2P1**(N-NC-1)
        DO 3 I=2,N1
            IF (XI(I).LT.DLOG(1.D0-QFIT).AND.NFOUT.EQ.0) THEN
                NFOUT=I
                GOTO 4
            ENDIF
3       CONTINUE
4       WRITE (*,*) ' NFOUT=',NFOUT
        END PROGRAM
