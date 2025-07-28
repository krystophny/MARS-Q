%merge all RE_LOSS_* files into one file
%according to FIDIST* files
%also merge all FIDIST* files into one file

SDIRM = '/home/liuy/Codes/MarsQ/RZplot/';
%SDIRD = '/home/liuy/Work/MAST-U/46943/Data/FIDIST_FO_E0.1/';
%SDIRD = '/home/liuy/Work/KSTAR/Data/30306.7850_gpec/FIDIST_5D/3D/';
SDIRD = '/home/liuy/Work/DIII-D/D3D177028/Data/FIDIST/RootAR/RE_PERTURB0.01_G/';
%SDIRD = '/home/liuy/Work/ITER/ITER_RE/DataAE_HALO/68ms/FIDIST/Root-c/RE_PERTURB0.1/';

PSS   = 'P10';

%go through all RE_LOSS_* files
cd(SDIRD);
FF = split(ls);
NF = length(FF)-1;

FIDIST = [];
RELOSS = [];

for kf=1:NF
    F = FF{kf};

%if length(F)>10 & strcmp(F(1:8),'RE_LOSS_') & strcmp(F(end-1:end),'_1')
%if strcmp(F,'FIDIST_ALL_SAVE.IN') | strcmp(F,'FIDIST_7A_SAVE.IN')
if length(F)==9 & strcmp(F(1:8),'RE_LOSS_')
       %SF = ['FIDIST_' PSS '_' F(9:end)];
       SF = ['../../FIDIST_' PSS '_' F(9:end) '.IN'];

       FI  = load(SF,'-ascii');
       RE  = load(F,'-ascii');

       NRE = size(RE,1);
       if (RE(1,5)+2)==NRE
          RE  = RE(3:end,:);
	  NRE = NRE-2;
       end

       if size(FI,1)-1 ~= NRE
          disp(['WARNING: wrong data for ' F])
       end

       FIDIST = [FIDIST; FI(2:end,:)];
       RELOSS = [RELOSS; RE];
    end
end

N = size(FIDIST,1);
FIDIST = [N*ones(1,size(FIDIST,2)); FIDIST];

save([SDIRD 'FIDIST_' PSS '_ALL'], 'FIDIST','-ascii','-double') 
save([SDIRD 'RE_LOSS_ALL'],'RELOSS','-ascii','-double') 

cd(SDIRM);
                         
