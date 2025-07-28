function [Br0,Bz0,Bphi0] = MacInterp2D(R,Z,Br,Bz,Bphi,R0,Z0)
%% Interperate B-field at a given point (R0,Z0) 

R1 = R(1:end-1,1:end-1);
Z1 = Z(1:end-1,1:end-1);

R2 = R(2:end,1:end-1);
Z2 = Z(2:end,1:end-1);

R3 = R(2:end,2:end);
Z3 = Z(2:end,2:end);

R4 = R(1:end-1,2:end);
Z4 = Z(1:end-1,2:end);

A1 = atan2(Z1-Z0,R1-R0)*180/pi; II=find(A1<0); A1(II) = A1(II)+360;
A2 = atan2(Z2-Z0,R2-R0)*180/pi; II=find(A2<0); A2(II) = A2(II)+360;
A3 = atan2(Z3-Z0,R3-R0)*180/pi; II=find(A3<0); A3(II) = A3(II)+360;
A4 = atan2(Z4-Z0,R4-R0)*180/pi; II=find(A4<0); A4(II) = A4(II)+360;

B1 = A2-A1; II=find(B1<-180); B1(II)=B1(II)+360; II=find(B1>180); B1(II)=B1(II)-360;
B2 = A3-A2; II=find(B2<-180); B2(II)=B2(II)+360; II=find(B2>180); B2(II)=B2(II)-360;
B3 = A4-A3; II=find(B3<-180); B3(II)=B3(II)+360; II=find(B3>180); B3(II)=B3(II)-360;
B4 = A1-A4; II=find(B4<-180); B4(II)=B4(II)+360; II=find(B4>180); B4(II)=B4(II)-360;

C  = B1 + B2 + B3 + B4;
[I1,I2] = find(abs(abs(C)-360)<10);

R1 = R1(I1,I2);  R2 = R2(I1,I2);  R3 = R3(I1,I2);  R4 = R4(I1,I2);
Z1 = Z1(I1,I2);  Z2 = Z2(I1,I2);  Z3 = Z3(I1,I2);  Z4 = Z4(I1,I2);

A = [1.0 R1 Z1 R1*Z1
     1.0 R2 Z2 R2*Z2
     1.0 R3 Z3 R3*Z3
     1.0 R4 Z4 R4*Z4];

B = Br;  Y = [B(I1,I2); B(I1+1,I2); B(I1+1,I2+1); B(I1,I2+1)];
C = A\Y;
Br0 = C(1) + C(2)*R0 + C(3)*Z0 + C(4)*R0*Z0;

B = Bz;  Y = [B(I1,I2); B(I1+1,I2); B(I1+1,I2+1); B(I1,I2+1)];
C = A\Y;
Bz0 = C(1) + C(2)*R0 + C(3)*Z0 + C(4)*R0*Z0;

B = Bphi;  Y = [B(I1,I2); B(I1+1,I2); B(I1+1,I2+1); B(I1,I2+1)];
C = A\Y;
Bphi0 = C(1) + C(2)*R0 + C(3)*Z0 + C(4)*R0*Z0;


