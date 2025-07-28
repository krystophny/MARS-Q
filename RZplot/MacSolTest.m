% read in SOLTEST.OUT and plot results
% for relative error (RE) of force balance equations
% RE is calulated using mean value of solution vector as normalization factor
% note that V1-3 below is the relative errors for the momentum equations for V1-3

SDIR = '/cscratch/liuy/WorkIR_D/BU/';
kplot = 1; %1:linear scale; 2:log10 scale 

d = load([SDIR 'PROFEQ.OUT']);
s = d(2:end-1,1).^2;

d = load([SDIR 'SOLTEST.OUT']);
N = d(1,1);
M = d(1,2); 
V1 = d(2:end,1);
V2 = d(2:end,2);
V3 = d(2:end,3);

V1 = reshape(V1,M,N); 
V2 = reshape(V2,M,N);
V3 = reshape(V3,M,N);

M0 = floor(M/2) + 1;
mm = [1:M]-M0;

NFIT = 2;
for k=1:NFIT
    V1(:,k) = V1(:,NFIT+1);
    V2(:,k) = V2(:,NFIT+1);
end

NFIT=2;
for k=0:NFIT-1
    V1(:,end-k) = V1(:,end-NFIT);
end

RES_mean = [mean(mean(V1)) mean(mean(V2)) mean(mean(V3))];
RES_maxm = [max(max(V1)) max(max(V2)) max(max(V3))];
RES = [RES_mean RES_maxm]  

TS = '';
if kplot==2
   TS = 'log_{10}';
   V1 = log10(V1);
   V2 = log10(V2);
   V3 = log10(V3);
end

hf=figure(1);
pcolor(s,mm,V1), shading interp
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('m','FontSize',18,'FontWeight','Bold')
title([TS '|RE for V^1-equation|'],'FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
a=axis; axis([0 1 a(3) a(4)])

hf=figure(2);
pcolor(s,mm,V2), shading interp
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('m','FontSize',18,'FontWeight','Bold')
title([TS '|RE for V^2-equation|'],'FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
a=axis; axis([0 1 a(3) a(4)])

hf=figure(3);
pcolor(s,mm,V3), shading interp
xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('m','FontSize',18,'FontWeight','Bold')
title([TS '|RE for V^3-equation|'],'FontSize',14)
colorbar,  colormap(jet)
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16)
a=axis; axis([0 1 a(3) a(4)])

