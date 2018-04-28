tic
% CLEANING WORKSPACE
clear

% OPENING FILES
F10 = fopen('NZOUT3D.txt','w+');
F17 = fopen('DIAGS3D.txt','w+');

% CONSTANTS AND VARIABLES

NX=34;
NY=34;
NZ=34;
NX1=NX-1;
NY1=NY-1;
NZ1=NZ-1;
NTEST=4;
NSTOP=256;
DELX=1./17;
DELY=1./17;
DELZ=1./17;
THINC=180;
PHINC=180;
ETHINC=1;
EPHINC=0;
AMP=1000;
BETA=64;
EPS0=8.854E-12;
XMU0=1.2566306E-6;
T=0;

% DEFINING ARRAYS
IDONE=zeros(NX,NY,NZ);
IDTWO=zeros(NX,NY,NZ);
IDTHRE=zeros(NX,NY,NZ);
EXS=zeros(NX,NY,NZ);
EYS=zeros(NX,NY,NZ);
EZS=zeros(NX,NY,NZ);
HXS=zeros(NX,NY,NZ);
HYS=zeros(NX,NY,NZ);
HZS=zeros(NX,NY,NZ);
STORE=zeros(1,4);
ESCTC=zeros(1,9);
EINCC=zeros(1,9);
EDEVCN=zeros(1,9);
ECRLX=zeros(1,9);
ECRLY=zeros(1,9);
ECRLZ=zeros(1,9);
EPS=zeros(1,9);
SIGMA=zeros(1,9);
EYSX1(1:4,1:NY1,1:NZ1)=zeros(4,33,33);
EYSX2(1:4,1:NY1,1:NZ1)=zeros(4,33,33);
EZSX1(1:4,1:NY1,1:NZ1)=zeros(4,33,33);
EZSX2(1:4,1:NY1,1:NZ1)=zeros(4,33,33);
EXSY1(1:NX1,1:4,1:NZ1)=zeros(33,4,33);
EXSY2(1:NX1,1:4,1:NZ1)=zeros(33,4,33);
EZSY1(1:NX1,1:4,1:NZ1)=zeros(33,4,33);
EZSY2(1:NX1,1:4,1:NZ1)=zeros(33,4,33);
EXSZ1(1:NX1,1:NY1,1:4)=zeros(33,33,4);
EXSZ2(1:NX1,1:NY1,1:4)=zeros(33,33,4);
EYSZ1(1:NX1,1:NY1,1:4)=zeros(33,33,4);
EYSZ2(1:NX1,1:NY1,1:4)=zeros(33,33,4);

for L=1:9,
        ESCTC(L)=0;
        EINCC(L)=0;
        EDEVCN(L)=0;
        ECRLX(L)=0;
        ECRLY(L)=0;
        ECRLZ(L)=0;
end

% RA BUILD SPHERE WITH CENTER(SC,SC,SC) AND RADIUS RA

MTYPE=2;
RA=8.2;
SC=17.5;

for I=1:NX,
  for J=1:NY,
    for K=1:NZ,
            R=sqrt((I-SC)^2+(J-SC)^2+(K-SC)^2);
            if (R <= RA), 
              IDONE(I,J,K)=MTYPE;
              IDONE(I,J,K+1)=MTYPE;
              IDONE(I,J+1,K+1)=MTYPE;
              IDONE(I,J+1,K)=MTYPE;
              IDTWO(I,J,K)=MTYPE;
              IDTWO(I+1,J,K)=MTYPE;
              IDTWO(I+1,J,K+1)=MTYPE;
              IDTWO(I,J,K+1)=MTYPE;
              IDTHRE(I,J,K)=MTYPE;
              IDTHRE(I+1,J,K)=MTYPE;
              IDTHRE(I+1,J+1,K)=MTYPE;
              IDTHRE(I,J+1,K)=MTYPE; 
            end
    end      
  end
end

% IF THE USER HAS SPECIFIED THE PROPER MATERIAL TYPES TO THE PROPER IDXXX 
% ARRAYS.

for I=1:NZ,
  for J=1:NY,
    for K=1:NX,
            if ((IDONE(I,J,K) >= 10) || (IDTWO(I,J,K) >= 10) || ...
                    (IDTHRE(I,J,K) >= 10))
                fprintf (F17,'ERROR OCCURED. ILLEGAL VALUE FOR DIELECTRIC \n');
                fprintf (F17,'TYPE AT LOCATION:');
                fprintf (F17,'%d %d %d','I','J','K','\n');
            end
    end
  end
end

% DEFINE PI AND C

C=1.0/sqrt(XMU0*EPS0);
PI=4.0*atan(1.0);

% CALCULATE DT - THE MAXIMUM TIME STEP ALLOWED BY THE COURANT STABILITY 
% CONDITION

DTXI=C/DELX;
DTYI=C/DELY;
DTZI=C/DELZ;
DT=1./sqrt(DTXI^2+DTYI^2+DTZI^2);
 
% PARAMETER ALPHA IS THE DECAY RATE DETERMINED BY BETA

ALPHA=(1./(BETA*DT/4))^2;

BETADT = BETA*DT;
PERIOD = 2*BETADT;

% SET OFFSET FOR COMPUTING INCIDENT FIELDS

OFF=1.0;

% VARIABLES FOR SMOOTH COSINE INCIDENT FUNCTION

COSTH=cos(PI*THINC/180);
SINTH=sin(PI*THINC/180);
COSPH=cos(PI*PHINC/180);
SINPH=sin(PI*PHINC/180);

% FIND AMPLITUDE OF INCIDENT FIELD COMPONENTS

AMPX=AMP*(ETHINC*COSTH*COSPH-EPHINC*SINPH);
AMPY=AMP*(ETHINC*COSTH*SINPH+EPHINC*COSPH);
AMPZ=AMP*(-ETHINC*SINTH);
 
% FIND RELATIVE SPATIAL DELAY FOR X, Y, Z CELL DISPLACEMENT

XDISP=-COSPH*SINTH;
YDISP=-SINPH*SINTH;
ZDISP=-COSTH;

% FOR PARTICULAR COMPONENTS OF FIELDS AS DETERMINED BY XXX

for I=1:9,
        EPS(I)=EPS0;
        SIGMA(I)=0;
end
     
% DEFINE EPS AND SIGMA FOR EACH MATERIAL HERE

EPS(2)=4*EPS0;
SIGMA(2)=0.005;

% GENERATE MULTIPLICATIVE CONSTANTS FOR FIELD UPDATE EQUATIONS 

% FREE SPACE

 DTEDX=DT/(EPS0*DELX);
 DTEDY=DT/(EPS0*DELY);
 DTEDZ=DT/(EPS0*DELZ);
 DTMDX=DT/(XMU0*DELX);
 DTMDY=DT/(XMU0*DELY);
 DTMDZ=DT/(XMU0*DELZ);

% LOSSY DIELECTRICS

for I=2:9,
        ESCTC(I)=EPS(I)/(EPS(I)+SIGMA(I)*DT);
        EINCC(I)=SIGMA(I)*DT/(EPS(I)+SIGMA(I)*DT);
        EDEVCN(I)=DT*(EPS(I)-EPS0)/(EPS(I)+SIGMA(I)*DT);
        ECRLX(I)=DT/((EPS(I)+SIGMA(I)*DT)*DELX);
        ECRLY(I)=DT/((EPS(I)+SIGMA(I)*DT)*DELY);
        ECRLZ(I)=DT/((EPS(I)+SIGMA(I)*DT)*DELZ);
end
   
% FIND MAXIMUM SPATIAL DELAY TO MAKE SURE PULSE PROPAGES INTO SPACE PROPERLY

DELAY=0;
if (XDISP <= 0) 
    DELAY=DELAY-XDISP*NX1*DELX;
end
if (YDISP <= 0) 
    DELAY=DELAY-YDISP*NY1*DELY;
end
if (ZDISP <= 0) 
    DELAY=DELAY-ZDISP*NZ1*DELZ;
end

% COMPUTE OUTER RADIATION BOUNDARY CONDITION (ORBC) CONSTANTS 

CXD=(C*DT-DELX)/(C*DT+DELX);
CYD=(C*DT-DELY)/(C*DT+DELY);
CZD=(C*DT-DELZ)/(C*DT+DELZ);
CXU=CXD;

% COMPUTE 2ND ORDER ORBC CONSTANTS

CXX=2*DELX/(C*DT+DELX);
CYY=2*DELY/(C*DT+DELY);
CZZ=2*DELZ/(C*DT+DELZ);

CXFYD=DELX*C*DT*C*DT/(2.*DELY*DELY*(C*DT+DELX));
CXFZD=DELX*C*DT*C*DT/(2.*DELZ*DELZ*(C*DT+DELX));
CYFZD=DELY*C*DT*C*DT/(2.*DELZ*DELZ*(C*DT+DELY));
CYFXD=DELY*C*DT*C*DT/(2.*DELX*DELX*(C*DT+DELY));
CZFXD=DELZ*C*DT*C*DT/(2.*DELX*DELX*(C*DT+DELZ));
CZFYD=DELZ*C*DT*C*DT/(2.*DELY*DELY*(C*DT+DELZ));

% WRITE SETUP DATA TO FILE
fprintf (F17,'THE PROBLEM SPACE %d %d %d CELLS X,Y,Z DIRECTIONS\n\n',NX,NY,NZ);
fprintf (F17,'CELL SIZE DELX %10.6f, DELY %10.6f ',DELX,DELY);
fprintf (F17,'DELZ %10.6f meter\n\n',DELZ);
fprintf (F17,'TIMESTEP %12.6E SECOND, ',DT);
fprintf (F17,'MAXIMUM NUMBER OF TIMESTEPS%d\n\n',NSTOP);
fprintf (F17,'INCIDENT GAUSSIAN PULSE AMPLITUDE %6.0f V/M, ',AMP);
fprintf (F17,'DECAY FACTOR ALPHA=%12.3E, ',ALPHA);
fprintf (F17,'WIDTH BETA=%6.0f\n\n',BETA);
fprintf (F17,'ıINCIDENT PLANE WAVE POLARIZATION\n');
fprintf (F17,'RELATIVE ELECTRIC FIELD THETA COMPONENT%4.1f\n',ETHINC);
fprintf (F17,'RELATIVE ELECTRIC FIELD PHI COMPONENT %4.1f\n',EPHINC);
fprintf (F17,'PLANE WAVE INCIDENT FROM THETA%d\nDERECE %d\n\n',THINC,PHINC');
fprintf (F17,'EX AMP %d V/M\nEY AMP%d V/M\nEZ AMP%d V/M\n\n',AMPX,AMPY,AMPZ);
fprintf (F17,'DELAY = %d\n\n',DELAY);

for N=1:NSTOP,

% SUBROUTINE EXSFLD

% THIS SUBROUTINE UPDATES THE EX SCATTERED FIELD 
 
    for K=2:NZ1,
     for J=2:NY1,
        for I=1:NX1,
            
% DETERMINE MATERIAL TYPE

             if (IDONE(I,J,K) == 0)
                 
% FREE SPACE

                 EXS(I,J,K)=EXS(I,J,K)+(HZS(I,J,K)-HZS(I,J-1,K))*DTEDY-(...
                 HYS(I,J,K)-HYS(I,J,K-1))*DTEDZ;
             elseif (IDONE(I,J,K) == 1)
                 
% PERFECT CONDUCTOR              
                 
                DIST=((I-1)*DELX+0.5*DELX*OFF)*XDISP+((J-1)*DELY)*YDISP+...
                    ((K-1)*DELZ)...
                *ZDISP + DELAY;
                source=0;
                TAU=T-DIST/C;
                if (TAU < 0) 
                elseif(TAU > PERIOD)
                else
                    source=exp(-ALPHA*((TAU-BETADT)^2));
                end

                EXI=AMPX*source;
                EXS(I,J,K)=-EXI;
                 
             else
                 
% LOSSY DIELECTRIC

                DIST=((I-1)*DELX+0.5*DELX*OFF)*XDISP+((J-1)*DELY)*...
                    YDISP+((K-1)*DELZ)*ZDISP + DELAY;
                source=0;
                dsrce=0;
                TAU=T-DIST/C;
                
                if (TAU < 0) 
                    elseif(TAU > PERIOD)
                else
                    source=exp(-ALPHA*((TAU-BETADT)^2));
                    dsrce=exp(-ALPHA*((TAU-BETADT)^2))*(-2*ALPHA*...
                        (TAU-BETADT));
                end
                
                % X component of coming wave
                exi=AMPX*source;
                  
                %
                dexi=AMPX*dsrce;

                EXS(I,J,K)=EXS(I,J,K)*ESCTC(IDONE(I,J,K))-...
                     EINCC(IDONE(I,J,K))*exi-...
                     EDEVCN(IDONE(I,J,K))*dexi+(HZS(I,J,K)-...
                     HZS(I,J-1,K))*ECRLY(IDONE(I,J,K))-(HYS(I,J,K)-...
                     HYS(I,J,K-1))*ECRLZ(IDONE(I,J,K));
             end
         end
     end
    end
    
% SUBROUTINE EYSFLD
    
    for K=2:NZ1,
     for J=1:NY1,
         for I=2:NX1,
             
% DETERMINE MATERIAL TYPE
             
             if (IDTWO(I,J,K) == 0)
             
% FREE SPACE
                 
             EYS(I,J,K)=EYS(I,J,K)+(HXS(I,J,K)-HXS(I,J,K-1))*DTEDZ...
                     -(HZS(I,J,K)-HZS(I-1,J,K))*DTEDX;
             elseif (IDTWO(I,J,K) == 1)
         
% PERFECT CONDUCTOR          

                DIST=((I-1)*DELX)*XDISP+((J-1)*DELY+0.5*DELY*OFF)*YDISP+...
                    ((K-1)*DELZ)*...
                ZDISP + DELAY;
                source=0;
                TAU=T-DIST/C;
                if (TAU < 0) 
                elseif(TAU > PERIOD)
                else
                    source=exp(-ALPHA*((TAU-BETADT)^2));
                end

                EYI=AMPX*source;
                EYS(I,J,K)=-EYI;
             else
                 
% LOSSY DIELECTRIC
                
                source=0;
                dsrce=0;
                DIST=((I-1)*DELX)*XDISP+((J-1)*DELY+0.5*DELY*OFF)*YDISP...
                    +((K-1)*DELZ)*...
                ZDISP + DELAY;
                TAU=T-DIST/C;

                if (TAU < 0) 
                    elseif(TAU > PERIOD)
                else
                    source=exp(-ALPHA*((TAU-BETADT)^2));
                    dsrce=exp(-ALPHA*((TAU-BETADT)^2))*(-2*ALPHA*...
                        (TAU-BETADT));
                end

                eyi=AMPY*source;
                deyi=AMPY*dsrce;
                
                EYS(I,J,K)=EYS(I,J,K)*ESCTC(IDTWO(I,J,K))-...
                     EINCC(IDTWO(I,J,K))*eyi-...
                     EDEVCN(IDTWO(I,J,K))*deyi+(HXS(I,J,K)-...
                     HXS(I,J,K-1))*ECRLZ(IDTWO(I,J,K))-(HZS(I,J,K)-...
                     HZS(I-1,J,K))*ECRLX(IDTWO(I,J,K));
             end
         end
     end
    end
    
% SUBROUTINE EZSFLD
      
    for K=1:NZ1,
        for J=2:NY1,
         for I=2:NX1,
             
% DETERMINE MATERIAL TYPE
           
            if (IDTHRE(I,J,K) == 0)
                
% FREE SPACE
                
                    EZS(I,J,K)=EZS(I,J,K)+(HYS(I,J,K)-HYS(I-1,J,K))*...
                        DTEDX-(HXS(I,J,K)-HXS(I,J-1,K))*DTEDY;
            elseif (IDTHRE(I,J,K) == 1)
                
% PERFECT CONDUCTOR
                    DIST=((I-1)*DELX)*XDISP+((J-1)*DELY)*YDISP+((K-1)*...
                        DELZ+0.5*DELZ*OFF)*ZDISP + DELAY;
                    source=0;
                    TAU=T-DIST/C;
                    if (TAU < 0) 
                    elseif(TAU > PERIOD)
                    else
                        source=exp(-ALPHA*((TAU-BETADT)^2));
                    end

                    EZI=AMPX*source;
                    EZS(I,J,K)=-EZI;
            else
                
% LOSSY DIELECTRIC
                    DIST=((I-1)*DELX)*XDISP+((J-1)*DELY)*YDISP+((K-1)*...
                        DELZ+0.5*DELZ*OFF)*ZDISP + DELAY;
                    dsrce=0;
                    source=0;
                    TAU=T-DIST/C;

                    if (TAU < 0) 
                    elseif(TAU > PERIOD)
                    else
                        source=exp(-ALPHA*((TAU-BETADT)^2));
                        dsrce=exp(-ALPHA*((TAU-BETADT)^2))*(-2*ALPHA*...
                            (TAU-BETADT));
                    end

                    ezi=AMPZ*source;
                    
                    dezi=AMPZ*dsrce;
                    
                    EZS(I,J,K)=EZS(I,J,K)*ESCTC(IDTHRE(I,J,K))-...
                    EINCC(IDTHRE(I,J,K))*ezi-...
                    EDEVCN(IDTHRE(I,J,K))*dezi+(HYS(I,J,K)-...
                    HYS(I-1,J,K))*ECRLX(IDTHRE(I,J,K))-(HXS(I,J,K)-...
                    HXS(I,J-1,K))*ECRLY(IDTHRE(I,J,K));
            end
         end
        end
     end

% SUBROUTINE RADEYX
% DO EDGES WITH FIRST ORDER ORBC
     
    for K=2:NZ1,
            J=1;
            EYS(1,J,K)=EYSX1(2,J,K)+CXD*(EYS(2,J,K)-EYSX1(1,J,K));
            EYS(NX,J,K)=EYSX1(3,J,K)+CXU*(EYS(NX1,J,K)-EYSX1(4,J,K));
            J=NY1;
            EYS(1,J,K)=EYSX1(2,J,K)+CXD*(EYS(2,J,K)-EYSX1(1,J,K));
            EYS(NX,J,K)=EYSX1(3,J,K)+CXU*(EYS(NX1,J,K)-EYSX1(4,J,K));
    end
    for J=2:NY1-1,
            K=2;
            EYS(1,J,K)=EYSX1(2,J,K)+CXD*(EYS(2,J,K)-EYSX1(1,J,K));
            EYS(NX,J,K)=EYSX1(3,J,K)+CXU*(EYS(NX1,J,K)-EYSX1(4,J,K));
            K=NZ1;
            EYS(1,J,K)=EYSX1(2,J,K)+CXD*(EYS(2,J,K)-EYSX1(1,J,K));
            EYS(NX,J,K)=EYSX1(3,J,K)+CXU*(EYS(NX1,J,K)-EYSX1(4,J,K));
    end

% DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES
    
    for K=3:NZ1-1,
            for J=2:NY1-1,
             EYS(1,J,K)=-EYSX2(2,J,K)+CXD*(EYS(2,J,K)+EYSX2(1,J,K))...
                  +CXX*(EYSX1(1,J,K)+EYSX1(2,J,K))+CXFYD*(EYSX1(1,J+1,K)...
                  -2.*EYSX1(1,J,K)+EYSX1(1,J-1,K)+EYSX1(2,J+1,K)-2.*...
                  EYSX1(2,J,K)+EYSX1(2,J-1,K))+CXFZD*(EYSX1(1,J,K+1)-2.*...
                  EYSX1(1,J,K)+EYSX1(1,J,K-1)+EYSX1(2,J,K+1)-2.*...
                  EYSX1(2,J,K)+EYSX1(2,J,K-1));
             EYS(NX,J,K)=-EYSX2(3,J,K)+CXD*(EYS(NX1 ,J,K)+EYSX2( 4,J,K))...
                  +CXX*(EYSX1(4,J,K)+EYSX1(3,J,K))+CXFYD*(EYSX1(4,J+1,K)...
                  -2.*EYSX1(4,J,K)+EYSX1(4,J-1,K)+EYSX1(3,J+1,K)-2.*...
                  EYSX1(3,J,K)+EYSX1(3,J-1,K))+CXFZD*(EYSX1(4,J,K+1)-2.*...
                  EYSX1(4,J,K)+EYSX1(4,J,K-1)+EYSX1(3,J,K+1)-2.*...
                  EYSX1(3,J,K)+EYSX1(3,J,K-1));
            end
    end
    
% SAVING THE VALUES
    for K=2:NZ1,
            for J=1:NY1,
              EYSX2(1,J,K)=EYSX1(1,J,K);
              EYSX2(2,J,K)=EYSX1(2,J,K);
              EYSX2(3,J,K)=EYSX1(3,J,K);
              EYSX2(4,J,K)=EYSX1(4,J,K);
              EYSX1(1,J,K)=EYS(1,J,K);
              EYSX1(2,J,K)=EYS(2,J,K);
              EYSX1(3,J,K)=EYS(NX1,J,K);
              EYSX1(4,J,K)=EYS(NX,J,K);
            end
    end

% SUBROUTINE RADEZX        
% DO EDGES WITH FIRST ORDER ORBC
    
    for K=1:NZ1,
        J=2;
        EZS(1,J,K)=EZSX1(2,J,K)+CXD*(EZS(2,J,K)-EZSX1(1,J,K));
        EZS(NX,J,K)=EZSX1(3,J,K)+CXU*(EZS(NX1,J,K)-EZSX1(4,J,K));
        J=NY1;
        EZS(1,J,K)=EZSX1(2,J,K)+CXD*(EZS(2,J,K)-EZSX1(1,J,K));
        EZS(NX,J,K)=EZSX1(3,J,K)+CXU*(EZS(NX1,J,K)-EZSX1(4,J,K));
    end
    for J=3:NY1-1,
        K=1;
        EZS(1,J,K)=EZSX1(2,J,K)+CXD*(EZS(2,J,K)-EZSX1(1,J,K));
        EZS(NX,J,K)=EZSX1(3,J,K)+CXU*(EZS(NX1,J,K)-EZSX1(4,J,K));
        K=NZ1;
        EZS(1,J,K)=EZSX1(2,J,K)+CXD*(EZS(2,J,K)-EZSX1(1,J,K));
        EZS(NX,J,K)=EZSX1(3,J,K)+CXU*(EZS(NX1,J,K)-EZSX1(4,J,K));
    end
    
% DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES 
    
    for K=2:NZ1-1,
        for J=3:NY1-1,
            EZS(1,J,K)=-EZSX2(2,J,K)+CXD*(EZS(2,J,K)+EZSX2(1,J,K))...
                +CXX*(EZSX1(1,J,K)+EZSX1(2,J,K))...
                +CXFYD*(EZSX1(1,J+1,K)-2.*EZSX1(1,J,K)+EZSX1(1,J-1,K)...
                +EZSX1(2,J+1,K)-2.*EZSX1(2,J,K)+EZSX1(2,J-1,K))...
                +CXFZD*(EZSX1(1,J,K+1)-2.*EZSX1(1,J,K)+EZSX1(1,J,K-1)...
                +EZSX1(2,J,K+1)-2.*EZSX1(2,J,K)+EZSX1(2,J,K-1));
            EZS(NX,J,K)=-EZSX2(3,J,K)+CXD*(EZS(NX1,J,K)+EZSX2(4,J,K))...
                +CXX*(EZSX1(4,J,K)+EZSX1(3,J,K))...
                +CXFYD*(EZSX1(4,J+1,K)-2.*EZSX1(4,J,K)+EZSX1(4,J-1,K)...
                +EZSX1(3,J+1,K)-2.*EZSX1(3,J,K)+EZSX1(3,J-1,K))...
                +CXFZD*(EZSX1(4,J,K+1)-2.*EZSX1(4,J,K)+EZSX1(4,J,K-1)...
                +EZSX1(3,J,K+1)-2.*EZSX1(3,J,K)+EZSX1(3,J,K-1));
        end
    end
    
% SAVING VALUES
    
    for K=1:NZ1,
        for J=2:NY1,
          EZSX2(1,J,K)=EZSX1(1,J,K);
          EZSX2(2,J,K)=EZSX1(2,J,K);
          EZSX2(3,J,K)=EZSX1(3,J,K);
          EZSX2(4,J,K)=EZSX1(4,J,K);
          EZSX1(1,J,K)=EZS(1,J,K);
          EZSX1(2,J,K)=EZS(2,J,K);
          EZSX1(3,J,K)=EZS(NX1,J,K);
          EZSX1(4,J,K)=EZS(NX,J,K);
        end
    end
    
% SUBROUTINE RADEZY
% DO EDGES WITH FIRST ORDER ORBC

    for K=1:NZ1,
        I=2;
        EZS(I,1,K)=EZSY1(I,2,K)+CYD*(EZS(I,2,K)-EZSY1(I,1,K));
        EZS(I,NY,K)=EZSY1(I,3,K)+CYD*(EZS(I,NY1,K)-EZSY1(I,4,K));
        I=NX1;
        EZS(I,1,K)=EZSY1(I,2,K)+CYD*(EZS(I,2,K)-EZSY1(I,1,K));
        EZS(I,NY,K)=EZSY1(I,3,K)+CYD*(EZS(I,NY1,K)-EZSY1(I,4,K));
    end
    for I=3:NX1-1,
        K=1;
        EZS(I,1,K)=EZSY1(I,2,K)+CYD*(EZS(I,2,K)-EZSY1(I,1,K));
        EZS(I,NY,K)=EZSY1(I,3,K)+CYD*(EZS(I,NY1,K)-EZSY1(I,4,K));
        K=NZ1;
        EZS(I,1,K)=EZSY1(I,2,K)+CYD*(EZS(I,2,K)-EZSY1(I,1,K));
        EZS(I,NY,K)=EZSY1(I,3,K)+CYD*(EZS(I,NY1,K)-EZSY1(I,4,K));
    end
    
% DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

    for K=2:NZ1-1,
        for I=3:NX1-1,
            EZS(I,1,K)=-EZSY2(I,2,K)+CYD*(EZS(I,2,K)+EZSY2(I,1,K))...
                +CYY*(EZSY1(I,1,K)+EZSY1(I,2,K))...
                +CYFXD*(EZSY1(I+1,1,K)-2.*EZSY1(I,1,K)+EZSY1(I-1,1,K)...
                +EZSY1(I+1,2,K)-2.*EZSY1(I,2,K)+EZSY1(I-1,2,K))...
                +CYFZD*(EZSY1(I,1,K+1)-2.*EZSY1(I,1,K)+EZSY1(I,1,K-1)...
                +EZSY1(I,2,K+1)-2.*EZSY1(I,2,K)+EZSY1(I,2,K-1));
          EZS(I,NY,K)=-EZSY2(I,3,K)+CYD*(EZS(I,NY1,K)+EZSY2(I,4,K))...
              +CYY*(EZSY1(I,4,K)+EZSY1(I,3,K))...
              +CYFXD*(EZSY1(I+1,4,K)-2.*EZSY1(I,4,K)+EZSY1(I-1,4,K)...
              +EZSY1(I+1,3,K)-2.*EZSY1(I,3,K)+EZSY1(I-1,3,K))...
              +CYFZD*(EZSY1(I,4,K+1)-2.*EZSY1(I,4,K)+EZSY1(I,4,K-1)...
              +EZSY1(I,3,K+1)-2.*EZSY1(I,3,K)+EZSY1(I,3,K-1));
        end
    end
    
% SAVING VALUES
    
    for K=1:NZ1,
        for I=2:NZ1,
            EZSY2(I,1,K)=EZSY1(I,1,K);
            EZSY2(I,2,K)=EZSY1(I,2,K);
            EZSY2(I,3,K)=EZSY1(I,3,K);
            EZSY2(I,4,K)=EZSY1(I,4,K);
            EZSY1(I,1,K)=EZS(I,1,K);
            EZSY1(I,2,K)=EZS(I,2,K);
            EZSY1(I,3,K)=EZS(I,NY1,K);
            EZSY1(I,4,K)=EZS(I,NY,K);
        end
    end
    
% SUBROUTINE RADEXY
% DO EDGES WITH FIRST ORDER ORBC

    for K=2:NZ1,
        I=1;
        EXS(I,1,K)=EXSY1(I,2,K)+CYD*(EXS(I,2,K)-EXSY1(I,1,K));
        EXS(I,NY,K)=EXSY1(I,3,K)+CYD*(EXS(I,NY1,K)-EXSY1(I,4,K));
        I=NX1;
        EXS(I,1,K)=EXSY1(I,2,K)+CYD*(EXS(I,2,K)-EXSY1(I,1,K));
        EXS(I,NY,K)=EXSY1(I,3,K)+CYD*(EXS(I,NY1,K)-EXSY1(I,4,K));
    end
    for I=2:NX-1,
        K=2;
        EXS(I,1,K)=EXSY1(I,2,K)+CYD*(EXS(I,2,K)-EXSY1(I,1,K));
        EXS(I,NY,K)=EXSY1(I,3,K)+CYD*(EXS(I,NY1,K)-EXSY1(I,4,K));
        K=NZ1;
        EXS(I,1,K)=EXSY1(I,2,K)+CYD*(EXS(I,2,K)-EXSY1(I,1,K));
        EXS(I,NY,K)=EXSY1(I,3,K)+CYD*(EXS(I,NY1,K)-EXSY1(I,4,K));
    end
    
% DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES
    
    for K=3:NZ1-1,
        for I=2:NX1-1,
            EXS(I,1,K)=-EXSY2(I,2,K)+CYD*(EXS(I,2,K)+EXSY2(I,1,K))...
                +CYY*(EXSY1(I,1,K)+EXSY1(I,2,K))...
                +CYFXD*(EXSY1(I+1,1,K)-2.*EXSY1(I,1,K)+EXSY1(I-1,1,K)...
                +EXSY1(I+1,2,K)-2.*EXSY1(I,2,K)+EXSY1(I-1,2,K))...
                +CYFZD*(EXSY1(I,1,K+1)-2.*EXSY1(I,1,K)+EXSY1(I,1,K-1)...
                +EXSY1(I,2,K+1)-2.*EXSY1(I,2,K)+EXSY1(I,2,K-1));
          EXS(I,NY,K)=-EXSY2(I,3,K)+CYD*(EXS(I,NY1,K)+EXSY2(I,4,K))...
              +CYY*(EXSY1(I,4,K)+EXSY1(I,3,K))...
              +CYFXD*(EXSY1(I+1,4,K)-2.*EXSY1(I,4,K)+EXSY1(I-1,4,K)...
              +EXSY1(I+1,3,K)-2.*EXSY1(I,3,K)+EXSY1(I-1,3,K))...
              +CYFZD*(EXSY1(I,4,K+1)-2.*EXSY1(I,4,K)+EXSY1(I,4,K-1)...
              +EXSY1(I,3,K+1)-2.*EXSY1(I,3,K)+EXSY1(I,3,K-1));
        end
    end
    
% SAVING VALUES

    for K=2:NZ1,
        for I=1:NX1,
            EXSY2(I,1,K)=EXSY1(I,1,K);
            EXSY2(I,2,K)=EXSY1(I,2,K);
            EXSY2(I,3,K)=EXSY1(I,3,K);
            EXSY2(I,4,K)=EXSY1(I,4,K);
            EXSY1(I,1,K)=EXS(I,1,K);
            EXSY1(I,2,K)=EXS(I,2,K);
            EXSY1(I,3,K)=EXS(I,NY1,K);
            EXSY1(I,4,K)=EXS(I,NY,K);
        end
    end

% SUBROUTINE RADEXZ
% DO EDGES WITH FIRST ORDER ORBC

    for J=2:NY1,
        I=1;
        EXS(I,J,1)=EXSZ1(I,J,2)+CZD*(EXS(I,J,2)-EXSZ1(I,J,1));
        EXS(I,J,NZ)=EXSZ1(I,J,3)+CZD*(EXS(I,J,NZ1)-EXSZ1(I,J,4));
        I=NX1;
        EXS(I,J,1)=EXSZ1(I,J,2)+CZD*(EXS(I,J,2)-EXSZ1(I,J,1));
        EXS(I,J,NZ)=EXSZ1(I,J,3)+CZD*(EXS(I,J,NZ1)-EXSZ1(I,J,4));
    end
    for I=2:NX1-1,
        J=2;
        EXS(I,J,1)=EXSZ1(I,J,2)+CZD*(EXS(I,J,2)-EXSZ1(I,J,1));
        EXS(I,J,NZ)=EXSZ1(I,J,3)+CZD*(EXS(I,J,NZ1)-EXSZ1(I,J,4));
        J=NY1;
        EXS(I,J,1)=EXSZ1(I,J,2)+CZD*(EXS(I,J,2)-EXSZ1(I,J,1));
        EXS(I,J,NZ)=EXSZ1(I,J,3)+CZD*(EXS(I,J,NZ1)-EXSZ1(I,J,4));
    end
    
% DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES

    for J=3:NY1-1,
        for I=2:NX1-1,
            EXS(I,J,1)=-EXSZ2(I,J,2)+CZD*(EXS(I,J,2)+EXSZ2(I,J,1))...
                +CZZ*(EXSZ1(I,J,1)+EXSZ1(I,J,2))...
                +CZFXD*(EXSZ1(I+1,J,1)-2.*EXSZ1(I,J,1)+EXSZ1(I-1,J,1)...
                +EXSZ1(I+1,J,2)-2.*EXSZ1(I,J,2)+EXSZ1(I-1,J,2))...
                +CZFYD*(EXSZ1(I,J+1,1)-2.*EXSZ1(I,J,1)+EXSZ1(I,J-1,1)...
                +EXSZ1(I,J+1,2)-2.*EXSZ1(I,J,2)+EXSZ1(I,J-1,2));
         EXS(I,J,NZ)=-EXSZ2(I,J,3)+CZD*(EXS(I,J,NZ1)+EXSZ2(I,J,4))...
             +CZZ*(EXSZ1(I,J,4)+EXSZ1(I,J,3))...
             +CZFXD*(EXSZ1(I+1,J,4)-2.*EXSZ1(I,J,4)+EXSZ1(I-1,J,4)...
             +EXSZ1(I+1,J,3)-2.*EXSZ1(I,J,3)+EXSZ1(I-1,J,3))...
             +CZFYD*(EXSZ1(I,J+1,4)-2.*EXSZ1(I,J,4)+EXSZ1(I,J-1,4)...
             +EXSZ1(I,J+1,3)-2.*EXSZ1(I,J,3)+EXSZ1(I,J-1,3));
        end
    end
    
% SAVING VALUES
    
    for J=2:NY1,
        for I=1:NX1,
            EXSZ2(I,J,1)=EXSZ1(I,J,1);
            EXSZ2(I,J,2)=EXSZ1(I,J,2);
            EXSZ2(I,J,3)=EXSZ1(I,J,3);
            EXSZ2(I,J,4)=EXSZ1(I,J,4);
            EXSZ1(I,J,1)=EXS(I,J,1);
            EXSZ1(I,J,2)=EXS(I,J,2);
            EXSZ1(I,J,3)=EXS(I,J,NZ1);
            EXSZ1(I,J,4)=EXS(I,J,NZ);
        end
    end

% SUBROUTINE RADEYZ
% DO EDGES WITH FIRST ORDER ORBC   
    
    for J=1:NY1,
        I=2;
        EYS(I,J,1)=EYSZ1(I,J,2)+CZD*(EYS(I,J,2)-EYSZ1(I,J,1));
        EYS(I,J,NZ)=EYSZ1(I,J,3)+CZD*(EYS(I,J,NZ1)-EYSZ1(I,J,4));
        I=NX1;
        EYS(I,J,1)=EYSZ1(I,J,2)+CZD*(EYS(I,J,2)-EYSZ1(I,J,1));
        EYS(I,J,NZ)=EYSZ1(I,J,3)+CZD*(EYS(I,J,NZ1)-EYSZ1(I,J,4));
    end
    for I=3:NX1-1,
        J=1;
        EYS(I,J,1)=EYSZ1(I,J,2)+CZD*(EYS(I,J,2)-EYSZ1(I,J,1));
        EYS(I,J,NZ)=EYSZ1(I,J,3)+CZD*(EYS(I,J,NZ1)-EYSZ1(I,J,4));
        J=NY1;
        EYS(I,J,1)=EYSZ1(I,J,2)+CZD*(EYS(I,J,2)-EYSZ1(I,J,1));
        EYS(I,J,NZ)=EYSZ1(I,J,3)+CZD*(EYS(I,J,NZ1)-EYSZ1(I,J,4));
    end
    
% DO 2ND ORDER ORBC ON REMAINING PORTIONS OF FACES 
    
    for J=2:NY1-1,
        for I=3:NX1-1,
            EYS(I,J,1)=-EYSZ2(I,J,2)+CZD*(EYS(I,J,2)+EYSZ2(I,J,1))...
                +CZZ*(EYSZ1(I,J,1)+EYSZ1(I,J,2))...
                +CZFXD*(EYSZ1(I+1,J,1)-2.*EYSZ1(I,J,1)+EYSZ1(I-1,J,1)...
                +EYSZ1(I+1,J,2)-2.*EYSZ1(I,J,2)+EYSZ1(I-1,J,2))...
                +CZFYD*(EYSZ1(I,J+1,1)-2.*EYSZ1(I,J,1)+EYSZ1(I,J-1,1)...
                +EYSZ1(I,J+1,2)-2.*EYSZ1(I,J,2)+EYSZ1(I,J-1,2));
          EYS(I,J,NZ)=-EYSZ2(I,J,3)+CZD*(EYS(I,J,NZ1)+EYSZ2(I,J,4))...
              +CZZ*(EYSZ1(I,J,4)+EYSZ1(I,J,3))...
              +CZFXD*(EYSZ1(I+1,J,4)-2.*EYSZ1(I,J,4)+EYSZ1(I-1,J,4)...
              +EYSZ1(I+1,J,3)-2.*EYSZ1(I,J,3)+EYSZ1(I-1,J,3))...
              +CZFYD*(EYSZ1(I,J+1,4)-2.*EYSZ1(I,J,4)+EYSZ1(I,J-1,4)...
              +EYSZ1(I,J+1,3)-2.*EYSZ1(I,J,3)+EYSZ1(I,J-1,3));
        end
    end
    
% SAVING VALUES
    
    for J=1:NY1,
        for I=2:NX1,
            EYSZ2(I,J,1)=EYSZ1(I,J,1);
            EYSZ2(I,J,2)=EYSZ1(I,J,2);
            EYSZ2(I,J,3)=EYSZ1(I,J,3);
            EYSZ2(I,J,4)=EYSZ1(I,J,4);
            EYSZ1(I,J,1)=EYS(I,J,1);
            EYSZ1(I,J,2)=EYS(I,J,2);
            EYSZ1(I,J,3)=EYS(I,J,NZ1);
            EYSZ1(I,J,4)=EYS(I,J,NZ);
        end
    end
    
    T=T+DT/2;
    
% SUBROUTINE HXSFLD
% UPDATES THE HX SCATTERED FIELD COMPONENTS

    for K=1:NZ1,
        for J=1:NY1,
            for I=2:NX1,
                HXS(I,J,K)=HXS(I,J,K)-(EZS(I,J+1,K)-EZS(I,J,K))*DTMDY...
                    +(EYS(I,J,K+1)-EYS(I,J,K))*DTMDZ;
            end
        end
    end

% SUBROUTINE HYSFLD    
% UPDATES THE HY SCATTERED FIELD COMPONENTS
    
    for K=1:NZ1,
        for J=2:NY1,
            for I=1:NX1,
                HYS(I,J,K)=HYS(I,J,K)-(EXS(I,J,K+1)-EXS(I,J,K))*DTMDZ...
                    +(EZS(I+1,J,K)-EZS(I,J,K))*DTMDX;
            end
        end
    end
    
% SUBROUTINE HZSFLD
% UPDATES THE HZ SCATTERED FIELD COMPONENTS


    for K=2:NZ1,
        for J=1:NY1,
            for I=1:NX1,
                HZS(I,J,K)=HZS(I,J,K)-(EYS(I+1,J,K)-EYS(I,J,K))*DTMDX...
                    +(EXS(I,J+1,K)-EXS(I,J,K))*DTMDY;
            end
        end
    end
    
    T=T+DT/2;
    
% SAVING VALUES

    if (N ~= 1) 
        
% DEFINE NTEST TEST POINT CELL LOCATION(S) HERE
        
        for NPT=1:NTEST,
            I=IOBS(NPT);
            J=JOBS(NPT);
            K=KOBS(NPT);
            
% DISTRIBUTED FIELDS

            if (NTYPE(NPT) == 1) 
                STORE(NPT)=EXS(I,J,K);
            end
            if (NTYPE(NPT) == 2)
                STORE(NPT)=EYS(I,J,K);
            end
            if (NTYPE(NPT) == 3)
                STORE(NPT)=EZS(I,J,K);
            end
            if (NTYPE(NPT) == 4) 
                STORE(NPT)=HXS(I,J,K);
            end
            if (NTYPE(NPT) == 5)
                STORE(NPT)=HYS(I,J,K);
            end
            if (NTYPE(NPT) == 6)
                STORE(NPT)=HZS(I,J,K);
            end
            
% DIRECTIONS
% X DIRECTION

            if (NTYPE(NPT) == 7) 
                STORE(NPT)=0;
                for KK=K:K+1,
                    for JJ=J:J+1,
                        STORE(NPT)=STORE(NPT)+(-HYS(I,JJ,KK)+...
                            HYS(I,JJ,KK-1))*DELY+(HZS(I,JJ,KK)-...
                            HZS(I,JJ-1,KK))*DELZ;
                    end
                end
            end
            
% Y DIRECTION
            
            if (NTYPE(NPT) == 8)
                STORE(NPT)=0;
                for KK=K:K+1,
                    for II=I:I+1,
                        STORE(NPT)=STORE(NPT)+(-HZS(II,J,KK)+...
                            HZS(II-1,J,KK))*DELZ+(HXS(II,J,KK)-...
                            HXS(II,J,KK-1))*DELX;
                    end
                end
            end
            
% Z DIRECTION
            
            if (NTYPE(NPT) == 9)
                STORE(NPT)=0;
                for JJ=J:J+1,
                    for II=I:I+1,
                        STORE(NPT)=STORE(NPT)+(-HXS(II,JJ,K)+...
                            HXS(II,JJ-1,K))*DELX+(HYS(II,JJ,K)-...
                            HYS(II-1,JJ,K))*DELY;
                    end
                end
            end
        end
    else
        
% DECISIONS OF VALUES WHICH WILL BE CALCULATED
%           1 = EXS (X component of distributed electric field)
%           2 = EYS (Y component of distributed electric field)
%           3 = EZS (Z component of distributed electric field)
%           4 = HXS (X component of distributed magnetic field)
%           5 = HYS (Y component of distributed magnetic field)
%           6 = HZS (Z componont of distributed magnetic field)
%           7 = IX (x-component of current through rectangular loop of H)
%           8 = IY (y-component of current through rectangular loop of H)
%           9 = IZ (z-component of current through rectangular loop of H)

        NTYPE(1)=1;
        NTYPE(2)=1;
        NTYPE(3)=5;
        NTYPE(4)=5;
        
% DEFINING COORDINATES OF TEST CELLS
        IOBS(1)=17;
        JOBS(1)=18;
        KOBS(1)=25;
        IOBS(2)=17;
        JOBS(2)=18;
        KOBS(2)=18;
        IOBS(3)=17;
        JOBS(3)=18;
        KOBS(3)=24;
        IOBS(4)=17;
        JOBS(4)=18;
        KOBS(4)=17;
        
        NPTS=NTEST;
        
% OUTPUT 

        fprintf (F10,'%d %d %d %d %d %d %d',DELX,DELY,DELZ,DT,NSTOP,NPTS);
        fprintf (F10,'\n');
        fprintf (F17,'Calculations are done and locations are saved to file %d\n',NPTS);
        
% NTYPE DATA TYPE CHECK

        for NPT=1:NPTS,     
fprintf (F17,'%d %d %d %d %d %d\n',NPT,NTYPE(NPT),IOBS(NPT),JOBS(NPT),KOBS(NPT));
            if (IOBS(NPT) >= NX) 
                fprintf (F17,...
'IOBS, JOBS or KOBS values have errors. %d Process is stopped.\n',NPT);
            end
            if (IOBS(NPT) <= 1)
                fprintf (F17,...
'IOBS, JOBS or KOBS values have errors. %d Process is stopped.\n',NPT);
            end
            if (JOBS(NPT) >= NY)
                fprintf (F17,...
'IOBS, JOBS or KOBS values have errors. %d Process is stopped.\n',NPT);
            end
            if (JOBS(NPT) <= 1) 
                fprintf (F17,...
'IOBS, JOBS or KOBS values have errors. %d Process is stopped.\n',NPT);
            end
            if (KOBS(NPT) >= NZ) 
                fprintf (F17,...
'IOBS, JOBS or KOBS values have errors. %d Process is stopped.\n',NPT);
            end
            if (KOBS(NPT) <= 1) 
                fprintf (F17,...
'IOBS, JOBS or KOBS values have errors. %d Process is stopped.\n',NPT);
            end
        end
        for NPT=1:NTEST,
             if ((NTYPE(NPT) >= 10) || (NTYPE(NPTS) <= 0))
fprintf (F17,'NTYPE Error at sample point %d Process is stopped\n',NPT);
             end
        end
        
        
    end   
    for II=1:NPTS,
        
% SAVING RESULTS THE OUTPUT FILE
   
             fprintf (F10,'%d ',STORE(II));
    end
    fprintf (F10,'\n');
end
T=NSTOP*DT;

fprintf (F17,'EXIT TIME %14.7f SECONDS, AT TIME STEP %d',T,NSTOP);

fclose(F10);
fclose(F17);

toc
