% read data BPLASMA_RECTRZ, plot, and save BnTOT.m for ASCOT 

SDIR = '/home/liuy/Work/MAST-U/46943/Data/'; BNORM=5.7660e-02; %==>dBnMaxPls=100G 

kplot2d = 1;
NZ      = 221;
plotBB  = 2;   %2D plots

dataB = load([SDIR 'BPLASMA_RECTRZ_Case1.OUT']);

R    = dataB(:,1);
Z    = dataB(:,2);
BR   = (dataB(:,3)+dataB(:,4)*i)*BNORM;
BZ   = (dataB(:,5)+dataB(:,6)*i)*BNORM;
BP   = (dataB(:,7)+dataB(:,8)*i)*BNORM;

NR   = round(length(R)/NZ);

R    = reshape(R,NZ,NR); 
Z    = reshape(Z,NZ,NR);
BR   = reshape(BR,NZ,NR);
BZ   = reshape(BZ,NZ,NR);
BP   = reshape(BP,NZ,NR);

B    = sqrt(abs(BR).^2 + abs(BZ).^2 + abs(BP).^2);

% plot |B| in the R-Z plane
if kplot2d == 1
     hf=figure(10*plotBB + 0);
     pcolor(R,Z,B*1e+4), hold on, shading interp
     axis equal
     xlabel('R [m]','FontSize',14)
     ylabel('Z [m]','FontSize',14)
     title('|B| (Gauss)','FontSize',14)
     colorbar,  colormap(jet)
     ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')
     maxAmpB = max(max(B))*1e+4
end

% save data for ASCOT
%note: need to change sign of BP to match ASCOT definition
BP = -BP;
Bn.r = R(1,:);
Bn.z=  Z(:,1);
Bn.Br = BR;
Bn.Bz = BZ;
Bn.Bphi =BP;
save BnTOT Bn

if kplot2d == 1

% plot BR in the R-Z plane
hf=figure(10*plotBB + 1);
pcolor(R,Z,real(BR)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Re[\delta{B}_R]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

hf=figure(10*plotBB + 2);
pcolor(R,Z,imag(BR)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Im[\delta{B}_R]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

% plot BZ in the R-Z plane
hf=figure(10*plotBB + 3);
pcolor(R,Z,real(BZ)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Re[\delta{B}_Z]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

hf=figure(10*plotBB + 4);
pcolor(R,Z,imag(BZ)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Im[\delta{B}_Z]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

% plot Re(BP) in the R-Z plane
hf=figure(10*plotBB + 5);
pcolor(R,Z,real(BP)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Re[\delta{B}_\phi]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

hf=figure(10*plotBB + 6);
pcolor(R,Z,imag(BP)), hold on, shading interp
axis equal
xlabel('R [m]','FontSize',14)
ylabel('Z [m]','FontSize',14)
title('Im[\delta{B}_\phi]','FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14,'FontWeight','Bold')

end

