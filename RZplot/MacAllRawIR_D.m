% go through database 
% generate and save raw data for Jpara(s,m) and X1U(s,m) 

global Mac

SDIR1 = '/cscratch/liuy/IR_D/Database/';
%SDIR1 = '/home/liuy/Work/IR_D/Database/DIII-D/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';

%spline all raw data on the same radial mesh s2
s2   = linspace(0,1,200);     
s2   = s2.^(1./(2+4*s2));     
s2   = [0 s2(2)/2 s2(2:end)]; 
s2   = s2(:);  %s-mesh for NN output raw data
CP2  = linspace(0,2*pi,9);    %coil phasing for two rows
CP3  = linspace(0,2*pi,9);    %coil phasing for three rows
EP2  = exp(i*CP2(1:8));
EP3  = exp(i*CP3(1:8));

cd(SDIR1)
D  = split(ls);
ND = length(D)-1;

D = {'DIII-D','ITER'}; ND = 2;

for kd=1:ND
DD = D{kd};
SDIR2 = [SDIR1 DD '/'];
if strcmp(DD,'MAST'),   C1='UL'; end
if strcmp(DD,'AUG'),    C1='UL'; end
if strcmp(DD,'DIII-D'), C1='ULC'; end
if strcmp(DD,'ITER'),   C1='ULM'; end

cd(SDIR2)    
F  = split(ls);
NF = length(F)-1;
for kf=1:NF
    FF = F{kf};
    SDIR3 = [SDIR2 FF '/3D/'];
    cd(SDIR3)
    G  = split(ls);
    NG = length(G)-1;
    for kn=1:NG
    GG = G{kn};
    SDIR4 = [SDIR3 GG '/'];
    cd(SDIR4)
    NC1 = length(C1);
    Jall = {};
    Xall = {};
    for kc=1:NC1
        CC = C1(kc);
        SDIR5 = [SDIR4 CC '/'];
        SDIR  = SDIR5;

        cd(SDIRM)
        if 1==0
        MacMainRawIR_D

        res  = [real(Mac.RAW_Jpara) imag(Mac.RAW_Jpara)];
        resn = zeros(length(s2),size(res,2));
        for k=1:size(res,2)
            resn(:,k) = spline(Mac.s(1:Mac.Ns1),res(:,k),s2);
        end
        save([SDIR5 'NNraw_Jpara.dat'],'resn','-ascii') 
        else
	resn=load([SDIR5 'NNraw_Jpara.dat']);
        end
        Jall{kc} = resn;

        if 1==0
        res  = [real(Mac.RAW_X1U) imag(Mac.RAW_X1U)];
        resn = zeros(length(s2),size(res,2));
        for k=1:size(res,2)
            resn(:,k) = spline(Mac.s(1:Mac.Ns1),res(:,k),s2);
        end
        save([SDIR5 'NNraw_X1U.dat'],'resn','-ascii') 
        else
        resn=load([SDIR5 'NNraw_X1U.dat']);
        end
        Xall{kc} = resn;
     end    

     %linear superposition: for two rows
     if length(C1)==2
     for kp=1:length(EP2)
         SDIR6 = [SDIR4 int2str(kp) '/'];
         mkdir(SDIR6);

         U  = Jall{1};
         L  = Jall{2};
         N  = round(size(U,2)/2);
         UC = U(:,1:N) + U(:,N+1:end)*i;
         LC = L(:,1:N) + L(:,N+1:end)*i;
         AC = UC*EP2(kp)+LC;
         res = [real(AC) imag(AC)];
         save([SDIR6 'NNraw_Jpara.dat'],'res','-ascii') 

         U  = Xall{1};
         L  = Xall{2};
         N  = round(size(U,2)/2);
         UC = U(:,1:N) + U(:,N+1:end)*i;
         LC = L(:,1:N) + L(:,N+1:end)*i;
         AC = UC*EP2(kp)+LC;
         res = [real(AC) imag(AC)];
         save([SDIR6 'NNraw_X1U.dat'],'res','-ascii') 
     end
     end
      
     %linear superposition: for three rows
     if length(C1)==3
     for kp1=1:length(EP3)
     for kp3=1:length(EP3)
         SDIR6 = [SDIR4 int2str(kp1) int2str(kp3) '/'];
         mkdir(SDIR6);

         U  = Jall{1};
         L  = Jall{2};
         M  = Jall{3};
         N  = round(size(U,2)/2);
         UC = U(:,1:N) + U(:,N+1:end)*i;
         LC = L(:,1:N) + L(:,N+1:end)*i;
         MC = M(:,1:N) + M(:,N+1:end)*i;
         AC = UC*EP3(kp1)+LC + MC*EP3(kp3);
         res = [real(AC) imag(AC)];
         save([SDIR6 'NNraw_Jpara.dat'],'res','-ascii') 

         U  = Xall{1};
         L  = Xall{2};
         M  = Xall{3};
         N  = round(size(U,2)/2);
         UC = U(:,1:N) + U(:,N+1:end)*i;
         LC = L(:,1:N) + L(:,N+1:end)*i;
         MC = M(:,1:N) + M(:,N+1:end)*i;
         AC = UC*EP3(kp1)+LC + MC*EP3(kp3);
         res = [real(AC) imag(AC)];
         save([SDIR6 'NNraw_X1U.dat'],'res','-ascii') 
     end
     end
     end
      
     end
end
end

cd(SDIRM)
