%read in error field (BR,BZ) at one surface
%convert into B1
%and save to INPUT_BNM_EF if Mac.RunEF=1
 
function MacReadFLOW

% read in data
fid = fopen('../Temp/FLOW_MARS_input.txt','r');
data=fscanf(fid,'%f');
fclose(fid);

NPSI    = round(data(1));
NTET    = round(data(2));
ISHOCK  = round(data(3));
IP      = 3;
data2   = data(IP+1:IP+NPSI*5);
data2   = reshape(data2,5,NPSI);
data2   = data2';
PSI     = data2(:,1);
FF      = data2(:,2);
GG      = data2(:,3);
FLOW_TOR= data2(:,4);
FLOW_POL= data2(:,5);
IP      = IP+NPSI*5;
data    = data(IP+1:end);
data    = reshape(data,3,NPSI*NTET);
data    = data';
RR      = data(:,1);
ZZ      = data(:,2);
RHO     = data(:,3);
RR      = reshape(RR,NTET,NPSI);
ZZ      = reshape(ZZ,NTET,NPSI);
RHO     = reshape(RHO,NTET,NPSI);


% plotting
TET = linspace(0,2*pi,NTET+1);
TET = TET(1:end-1);

figure(1)
plot(PSI,FF),  xlabel('\psi_p'),   ylabel('FF')

figure(2)
plot(PSI,GG),  xlabel('\psi_p'),   ylabel('GG')

figure(3)
plot(PSI,FLOW_TOR),  xlabel('\psi_p'),   ylabel('FLOW-TOR')

figure(4)
plot(PSI,FLOW_POL),  xlabel('\psi_p'),   ylabel('FLOW-POL')

figure(5)
surf(PSI,TET,RR),  shading interp, xlabel('\psi_p'),   ylabel('\theta'),  zlabel('R')

figure(6)
surf(PSI,TET,ZZ),  shading interp, xlabel('\psi_p'),   ylabel('\theta'),  zlabel('Z')

figure(7)
surf(PSI,TET,RHO),  shading interp, xlabel('\psi_p'),   ylabel('\theta'),  zlabel('\rho')

figure(8)
surf(RR,ZZ,RHO),  shading interp, xlabel('R'),   ylabel('Z'),  zlabel('\rho')

figure(9),
plot(TET,RHO(:,end)), xlabel('\theta'),   ylabel('\rho(edge)')

figure(10),
plot(TET,RHO(:,2)), xlabel('\theta'),   ylabel('\rho(core)')
