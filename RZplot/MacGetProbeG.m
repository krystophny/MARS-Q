%read in probe-g vacuum field computed on efit grid 
%save to fiels for MARS-F REORBIT runs
 
function MacGetProbeG(filename,R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)

global Mac
global SDIR SFLD

if Mac.RunEF==4
%parameter setting
NSP = 24;
NR  = 65;
NZ  = 65;
n   = 6;
mm = [-33:33];  

%read in data from PROBE_G code output
d = load(filename);
phi = d(:,1)*pi/180;
RR  = d(:,2);
ZZ  = d(:,3);
BP  = d(:,4);
BR  = d(:,5);
BZ  = d(:,6);

if NR*NZ*NSP ~= length(phi), disp('probe_gb_efitgrid.out: wrong size'); pause; end

phi = reshape(phi,NR,NZ,NSP);
RR  = reshape(RR,NR,NZ,NSP);
ZZ  = reshape(ZZ,NR,NZ,NSP);
BR  = reshape(BR,NR,NZ,NSP);
BZ  = reshape(BZ,NR,NZ,NSP);
BP  = reshape(BP,NR,NZ,NSP);

x = zeros(1,NSP);  for k=1:NSP, x(k)=phi(1,1,k); end, phi = x;
RR  = RR(:,:,1);
ZZ  = ZZ(:,:,1);

if 1==0
%read in vacuum data on control surface to do comparison
d = load([SDIR '../Expt/probe_gb.out']);
RR1  = d(:,2);
ZZ1  = d(:,3);
BP1  = d(:,4);
BR1  = d(:,5);
BZ1  = d(:,6);

NST = floor(length(RR1)/NSP);
RR1  = reshape(RR1,NST,NSP); RR1 = RR1(:,1);
ZZ1  = reshape(ZZ1,NST,NSP); ZZ1 = ZZ1(:,1);
BR1  = reshape(BR1,NST,NSP);
BZ1  = reshape(BZ1,NST,NSP);
BP1  = reshape(BP1,NST,NSP);

%interpolate along control surface and compare with control surface data
kphi = 5;
BR2 = griddata(RR,ZZ,BR(:,:,kphi),RR1,ZZ1);
BZ2 = griddata(RR,ZZ,BZ(:,:,kphi),RR1,ZZ1);
BP2 = griddata(RR,ZZ,BP(:,:,kphi),RR1,ZZ1);

if Mac.plot_EF > 0
   tt = atan2(ZZ1,RR1-Mac.R0EXP); [tt,II]=sort(tt);

   hf=figure(10*Mac.plot_EF+1);
   plot(tt,BR1(II,kphi),'r-',tt,BR2(II),'b-','LineWidth',2), hold on
   xlabel('poloidal angle','FontSize',16,'FontWeight','Bold')
   ylabel('{\delta}B_R','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+2);
   plot(tt,BZ1(II,kphi),'r-',tt,BZ2(II),'b-','LineWidth',2), hold on
   xlabel('poloidal angle','FontSize',16,'FontWeight','Bold')
   ylabel('{\delta}B_Z','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+3);
   plot(tt,BP1(II,kphi),'r-',tt,BP2(II),'b-','LineWidth',2), hold on
   xlabel('poloidal angle','FontSize',16,'FontWeight','Bold')
   ylabel('{\delta}B_{\phi}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
end
end

%Fourier decomposition along phi-angle
BR = permute(BR,[3 1 2]);
BZ = permute(BZ,[3 1 2]);
BP = permute(BP,[3 1 2]);

dphi = phi(2)-phi(1);
ephi = exp(-1i*n*phi);

X = RR; Y = X; P = X;
for k=1:NR
for j=1:NZ
    X(k,j) = ephi*BR(:,k,j);
    Y(k,j) = ephi*BZ(:,k,j);
    P(k,j) = ephi*BP(:,k,j);
end
end

BR = X*dphi/2/pi;
BZ = Y*dphi/2/pi;
BP = P*dphi/2/pi;

%map into (s,chi)-mesh
BR = griddata(RR/Mac.R0EXP,ZZ/Mac.R0EXP,BR/Mac.B0EXP,R,Z);
BZ = griddata(RR/Mac.R0EXP,ZZ/Mac.R0EXP,BZ/Mac.B0EXP,R,Z);
BP = griddata(RR/Mac.R0EXP,ZZ/Mac.R0EXP,BP/Mac.B0EXP,R,Z);

BR(isnan(BR)) = 0;
BZ(isnan(BZ)) = 0;
BP(isnan(BP)) = 0;

%convert (BR,BZ,BP) to (B1,B2,B3)
B1 = R.*(dZdchi.*BR-dRdchi.*BZ);
B2 = R.*(dRds.*BZ-dZds.*BR);
B3 = jacobian.*BP./R;

%Fourier decomposition along poloidal angle
chi=Mac.chi(1:end-1); chi=chi(:); hchi=(Mac.chi(2)-Mac.chi(1))/2/pi;
expchi = exp(-1i*chi*mm);
B1m    = B1(:,1:end-1)*expchi*hchi;
B2m    = B2(:,1:end-1)*expchi*hchi;
B3m    = B3(:,1:end-1)*expchi*hchi;

%B2 and B3 defined at half-integer radial mesh
B2m(1:end-1,:) = (B2m(1:end-1,:)+B2m(2:end,:))/2;
B3m(1:end-1,:) = (B3m(1:end-1,:)+B3m(2:end,:))/2;

%cut the last radial point
B1m = B1m(1:end-1,:);
B2m = B2m(1:end-1,:);
B3m = B3m(1:end-1,:);
 
%save data 
res = [Mac.Ns1-1 Mac.Ns-Mac.Ns1 mm(1) mm(end) n 0; ...
       real(B1m(:)) imag(B1m(:)) real(B2m(:)) imag(B2m(:)) real(B3m(:)) imag(B3m(:))];

save(['FIELD_RE_PROBEG_n' int2str(n)],'res','-ascii','-double')

%plot field harmonics
if Mac.plot_EF > 0
   hf=figure(10*Mac.plot_EF+4);
   plot(Mac.s(1:end-1),real(B1m),'LineWidth',1), hold on
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Re(B^1)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+5);
   plot(Mac.s(1:end-1),imag(B1m),'LineWidth',1), hold on
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Im(B^1)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+6);
   plot(Mac.s(1:end-1),real(B2m),'LineWidth',1), hold on
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Re(B^2)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+7);
   plot(Mac.s(1:end-1),imag(B2m),'LineWidth',1), hold on
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Im(B^2)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+8);
   plot(Mac.s(1:end-1),real(B3m),'LineWidth',1), hold on
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Re(B^3)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

   hf=figure(10*Mac.plot_EF+9);
   plot(Mac.s(1:end-1),imag(B3m),'LineWidth',1), hold on
   xlabel('s','FontSize',16,'FontWeight','Bold')
   ylabel('Im(B^3)','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
end

end


