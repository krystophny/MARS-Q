% go through database 
% collect global databse matrices for Jpara(s,m) and X1U(s,m)
% then perform global SVD for the database matrices
% and save (S*V)-coefficients back to each datapoint for next use as NN output 
% database includes 3D fields from linear superposition of different rows

global Mac

SDIR1 = '/cscratch/liuy/IR_D/Database/';
SDIR0 = '/home/liuy/Work/IR_D/Database/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';

%collect NNraw_Jpara.dat and NNraw_X1U.dat and form database matrices
NK  = 5;  %truncation of SVD to first NK S-values and eigenvectors
NM  = 67; %number of poloidal Fourier harmonics (must be fixed at 67)
NN  = 13467; %size(MatJ,1)
NS_MAST = 35 ; %number of shots for MAST
NS_D3D  = 27;  %number of shots for DIII-D
NS_AUG  = 13;  %number of shots for AUG
NS_ITER = 32;  %number of shots for ITER

%DEV = {'MAST','DIII-D','AUG','ITER'};
%DEV = {'DIII-D','AUG','ITER'};
DEV = {'MAST'};

%NCASES = {'n1','n2','n3','n4'};
NCASES = {'n4'};

SOUT = 'nn';
if length(NCASES)==1, SOUT=NCASES{1}; end

NS = 0;
for kd=1:size(DEV,2)
    if strcmp(DEV{kd},'DIII-D'), NS=NS+NS_D3D*(3+64); end
    if strcmp(DEV{kd},'AUG'),    NS=NS+NS_AUG*(2+8); end
    if strcmp(DEV{kd},'ITER'),   NS=NS+NS_ITER*(3+64); end
    if strcmp(DEV{kd},'MAST'),   NS=NS+NS_MAST*(2+8); end
end

NS = NS*length(NCASES);

I1   = [1:NM];
I2   = [NM+1:2*NM];
MatJ = zeros(NN,NS);
MatX = MatJ;
ND   = size(DEV,2);

ka   = 0;
for kd=1:ND
DD = DEV{kd};
SDIR2 = [SDIR1 DD '/'];
cd(SDIR2)    
F  = split(ls);
NF = length(F)-1;
for kf=1:NF
    FF = F{kf};
    SDIR3 = [SDIR2 FF '/3D/'];
    cd(SDIR3)
    G   = NCASES; 
    NG = length(G);
    for kn=1:NG
    GG = G{kn};
    SDIR4 = [SDIR3 GG '/'];
    cd(SDIR4)
    C = split(ls); C1={};
    NC = length(C)-1;
    kc1 = 0;
    for kc=1:NC
        if length(C{kc})<=2
           kc1 = kc1+1;
           C1{kc1}=C{kc};
        end
    end
    NC1 = kc1;
    for kc=1:NC1
        CC = C1{kc};
        SDIR5 = [SDIR4 CC '/'];
        ka    = ka + 1;

        d = load([SDIR5 'NNraw_Jpara.dat']);
        d = d(:,I1)+d(:,I2)*i;
        MatJ(:,ka) = d(:);

        d = load([SDIR5 'NNraw_X1U.dat']);
        d = d(:,I1)+d(:,I2)*i;
        MatX(:,ka) = d(:);
     end     
     end
end
end

%global SVD for MatJ and MatX
%test accuracy of SVD
disp('start SVD for Jpara...')
[U,S,V] = svd(MatJ);
V_J = transpose(conj(V));
SV  = S*V_J;
hf=figure(1);
SS = diag(S);
plot(SS(1:10),'b-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b')
xlabel('SVD number','FontSize',18,'FontWeight','Bold')
ylabel('S-value for J^{||}','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes');
set(ha,'FontSize',16,'FontWeight','Bold')

MatJn = U(:,1:NK)*SV(1:NK,:);
UM = sum(sum(abs(MatJ).^2));
UD = sum(sum(abs(MatJ-MatJn).^2));
UE_J = UD/UM;

hf=figure(2);
plot(log10(abs(V_J(1:NK,:))),angle(V_J(1:NK,:)),'bo')
xlabel('log_{10}|V| for J_{||}')
ylabel('phase(V) for J_{||}')

disp('start SVD for X1U...')
[U,S,V] = svd(MatX);
V_X = transpose(conj(V));
SV  = S*V_X;
hf=figure(3);
SS = diag(S);
plot(SS(1:10),'b-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b')
xlabel('SVD number','FontSize',18,'FontWeight','Bold')
ylabel('S-value for \xi^1','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes');
set(ha,'FontSize',16,'FontWeight','Bold')

MatXn = U(:,1:NK)*SV(1:NK,:);
UM = sum(sum(abs(MatX).^2));
UD = sum(sum(abs(MatX-MatXn).^2));
UE_X = UD/UM;

hf=figure(4);
plot(log10(abs(V_X(1:NK,:))),angle(V_X(1:NK,:)),'bo')
xlabel('log_{10}|V| for \xi^1')
ylabel('phase(V) for \xi^1')

RES_MSE = [UE_J UE_X]

UN = sum(abs(U).^2,1);
figure(5)
plot(UN,'o')

MM = size(V,2);
II = 1:MM-1;
VN = abs(V);
figure(6)
plot([1:MM],VN(:,II),'-')

if 1==1
%save first NK V-values to each datapoint
disp('save V-values for each datapoint...')
ka   = 0;
for kd=1:ND
DD = DEV{kd};
SDIR2 = [SDIR1 DD '/'];
cd(SDIR2)    
F  = split(ls);
NF = length(F)-1;
for kf=1:NF
    FF = F{kf};
    SDIR3 = [SDIR2 FF '/3D/'];
    cd(SDIR3)
    G   = NCASES;
    NG = length(G);
    for kn=1:NG
    GG = G{kn};
    SDIR4 = [SDIR3 GG '/'];
    cd(SDIR4)
    C = split(ls); C1={};
    NC = length(C)-1;
    kc1 = 0;
    for kc=1:NC
        if length(C{kc})<=2
           kc1 = kc1+1;
           C1{kc1}=C{kc};
        end
    end
    NC1 = kc1;
    for kc=1:NC1
        CC = C1{kc};
        SDIR5 = [SDIR4 CC '/'];
        SDIR6 = [SDIR0(1:15) SDIR5(15:end)];
        ka    = ka + 1;

        d = V_J(1:NK,ka);
        res = [real(d); imag(d)];
        save([SDIR6 'SVD1_Jpara_' SOUT '.txt'],'res','-ascii') 

        d = V_X(1:NK,ka);
        res = [real(d); imag(d)];
        save([SDIR6 'SVD1_X1U_' SOUT '.txt'],'res','-ascii') 
     end     
     end
end
end

end

cd(SDIRM)

