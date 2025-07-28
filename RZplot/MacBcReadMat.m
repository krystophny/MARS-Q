function [B1,B2,B3,psiz] = MacBcReadMat(filename,psiz,kk)
%% Read input coupling matrices for backward coupling scheme

global Fc A1 A1a A2a A2b A3

data0 = load(filename); 

%read data
K = round(size(data0,1)/3);
data = data0(1:K,:);  data=data';
N = round(size(data,1)/2); M = size(data,2); 
B1 = data(1:N,:) + data(N+1:end,:)*i;
data = data0(K+1:2*K,:); data=data';
B2 = data(1:N,:) + data(N+1:end,:)*i;
data = data0(2*K+1:3*K,:); data=data';
B3 = data(1:N,:) + data(N+1:end,:)*i;

%coordinate transform
if kk==1
if Fc.cariddi>0
  B1inv = inv(B1)*A1;
  B2 = (A2a*B2-A2b*B1)*B1inv;
  B3 = (A3*B3)*B1inv;
  B1 = eye(size(B1));
  psiz = psiz*B1inv;
else
  B1inv = inv(B1);
  B2 = B2*B1inv;
  B3 = B3*B1inv;
  B1 = eye(size(B1));
end
end

if kk==2 & Fc.cariddi>0
B2 = (A2a*B2-A2b*B1)*Fc.fac;
B1 = A1a*B1*Fc.fac;
B3 = A3*B3*Fc.fac;
end

%replace B2 component by B3 for m=0 harmonic
M0=find(abs(Fc.mm)<1e-10);
B2(M0,:) = B3(M0,:); 

%select subset of harmonics
ii=Fc.ii;
if size(B1,1)==size(B1,2)
   B1 = B1(ii,ii); B2 = B2(ii,ii); B3 = B3(ii,ii); 
else
   B1 = B1(ii,:); B2 = B2(ii,:); B3 = B3(ii,:); 
end   

