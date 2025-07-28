% select and form NN database for 3D equilibria
% according to the input/output format defined in drawMemo.m 
% use DIII-D and AUG data for training, and ITER data for testing

SDIR1 = '/home/liuy/Work/IR_D/Database/'; 
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';
SDIRO = '/home/liuy/Work/IR_D/DataNN/MixedV1/'; 

%DEV = {'MAST','DIII-D','AUG','ITER'};
DEV = {'DIII-D','AUG','ITER'};
%DEV = {'ITER'};
SN  = '1234';
%SN = '4';
LN  = zeros(1,size(DEV,2));
kUV = 1;
koutput =-3; %-1:         select all S-values as output
             %-2:         select all U(kUV) as output
             %-3:         select all V(kUV) as output
             % 0:         select all [S, U, V] as output
             % k=1,...,5: select [S(k), U(k), V(k)] as output 
rfrac = 0.8;  %fraction of random selection for training data

if length(SN)==1
   SDIRO = [SDIRO 'n' SN '/'];
else
   SDIRO = [SDIRO 'nn/'];
end

NNdataJ = [];
NNdataX = [];
for kd=1:size(DEV,2)
    if strcmp(DEV{kd},'DIII-D'), SC='ULC'; end
    if strcmp(DEV{kd},'AUG'),    SC='UL';  end
    if strcmp(DEV{kd},'ITER'),   SC='ULM'; end
    if strcmp(DEV{kd},'MAST'),   SC='UL'; end

    SDIR2 = [SDIR1 DEV{kd} '/'];
    cd(SDIR2)
    F  = split(ls);
    NF = length(F)-1
    cd(SDIRM)

    for kf=1:NF
        FF = F{kf};
        for kn=1:length(SN)
        for kc=1:length(SC)
            SDIR4 = [SDIR2 FF '/3D/n' SN(kn) '/'];
            SDIR5 = [SDIR4 SC(kc) '/'];
            
            d = load([SDIR5 'NNdata_Jpara.dat']);
            if length(SN)>1, 
               d=[d(1:8); str2num(SN(kn)); d(9:end)]; 
               d(1) = d(1) + 1;
            end
            NNdataJ = [NNdataJ d];

            d = load([SDIR5 'NNdata_X1U.dat']);
            %ke = 345;
            %if log10(d(ke))>3.5 | log10(d(ke))<-0.5
            %   disp([SDIR5 ': S=' num2str(d(ke))]) 
            %else   
            if length(SN)>1, 
               d=[d(1:8); str2num(SN(kn)); d(9:end)]; 
               d(1) = d(1) + 1;
            end
            NNdataX = [NNdataX d];
            %end
         end
         end
    end
end

NIN = round(d(1));
NOU = round(d(2));
NRZ = round(d(3));
NPR = round(d(4));
NBV = round(d(5));
NSS = round(d(6));
NUS = round(d(7));
NVM = round(d(8));

NIN = NIN - NBV;
NNdataJ = [NNdataJ(9:NIN+8,:); NNdataJ(9+NIN+NBV:end,:)];
NNdataX = [NNdataX(9:NIN+8,:); NNdataX(9+NIN+NBV:end,:)];

%take log10(S) as NN output
if 1==1
NNdataJ(NIN+1:NIN+NSS,:) = log10(NNdataJ(NIN+1:NIN+NSS,:));
NNdataX(NIN+1:NIN+NSS,:) = log10(NNdataX(NIN+1:NIN+NSS,:));
end

hf=figure(1);
plot(NNdataJ(NIN+1,:),'bo')
title('S(1) for J_{||}')
hf=figure(2);
plot(NNdataX(NIN+1,:),'bo')
title('S(1) for \xi_n')

% select data with S(1)<?
if 1==1
II = find(NNdataJ(NIN+1,:)<200);
NNdataJ = NNdataJ(:,II);
II = find(NNdataX(NIN+1,:)<3.5 & NNdataX(NIN+1,:)>-0.5);
NNdataX = NNdataX(:,II);
end

%select output data for NN
if koutput==-1
   %select only S-values as NN output
   NNdataJ = NNdataJ(1:NIN+NSS,:);
   NNdataX = NNdataX(1:NIN+NSS,:);
   NOU = NSS;
elseif koutput==-2
   %select only U-values as NN output
   k   = kUV;
   II  = 1:NUS; II=II(:);
   IIU = NIN + NSS + [(k-1)*NUS+II; (NSS+k-1)*NUS+II];
   NNdataJ = [NNdataJ(1:NIN,:); NNdataJ(IIU,:)];
   NNdataX = [NNdataX(1:NIN,:); NNdataX(IIU,:)];
   NSS = 0;
   NVM = 0;
   NOU = NUS*2;
elseif koutput==-3
   %select only U-values as NN output
   k   = kUV;
   II  = 1:NVM; II=II(:);
   IIV = NIN + NSS + 2*NUS*NSS + [(k-1)*NVM+II; (NSS+k-1)*NVM+II];
   NNdataJ = [NNdataJ(1:NIN,:); NNdataJ(IIV,:)];
   NNdataX = [NNdataX(1:NIN,:); NNdataX(IIV,:)];
   NSS = 0;
   NUS = 0;
   NOU = NVM*2;
elseif koutput>0
   %select [S(k), U(k), V(k)], k=koutput, as NN output
   k   = koutput;
   IIS = NIN + k;
   II  = 1:NUS; II=II(:);
   IIU = NIN + NSS + [(k-1)*NUS+II; (NSS+k-1)*NUS+II];
   II  = 1:NVM; II=II(:);
   IIV = NIN + NSS + 2*NUS*NSS + [(k-1)*NVM+II; (NSS+k-1)*NVM+II];
   NNdataJ = [NNdataJ(1:NIN,:); NNdataJ([IIS;IIU;IIV],:)];
   NNdataX = [NNdataX(1:NIN,:); NNdataX([IIS;IIU;IIV],:)];
   NSS = 1;
   NOU = 1 + (NUS + NVM)*2;
end

% randomly select rfrac datapoints for training and (1-rfrac) for testing
NSJA = size(NNdataJ,2);
NSJ1 = round(NSJA*rfrac);
NSJ2 = NSJA - NSJ1;
NSXA = size(NNdataX,2);
NSX1 = round(NSXA*rfrac);
NSX2 = NSXA - NSX1;
II  = randperm(NSJA);
NNdataJ1 = NNdataJ(:,II(1:NSJ1));
NNdataJ2 = NNdataJ(:,II(NSJ1+1:end));
II  = randperm(NSXA);
NNdataX1 = NNdataX(:,II(1:NSX1));
NNdataX2 = NNdataX(:,II(NSX1+1:end));

NALL = [NIN; NOU; NRZ; NPR; NBV; NSS; NUS; NVM];
FJ1   = [repmat(' %15.7e',1,NSJ1),'\n'];
FJ2   = [repmat(' %15.7e',1,NSJ2),'\n'];
FX1   = [repmat(' %15.7e',1,NSX1),'\n'];
FX2   = [repmat(' %15.7e',1,NSX2),'\n'];

if NSJ1>0
fid = fopen([SDIRO 'NN_TRAIN_J.dat'],'w');
fprintf(fid,'%6i %6i %6i %6i %6i %6i %6i %6i %6i\n',[NSJ1; NALL]);
fprintf(fid,FJ1,transpose(NNdataJ1));
fclose(fid);
end

if NSX1>0
fid = fopen([SDIRO 'NN_TRAIN_X.dat'],'w');
fprintf(fid,'%6i %6i %6i %6i %6i %6i %6i %6i %6i\n',[NSX1; NALL]);
fprintf(fid,FX1,transpose(NNdataX1));
fclose(fid);
end

if NSJ2>0
fid = fopen([SDIRO 'NN_TEST_J.dat'],'w');
fprintf(fid,'%6i %6i %6i %6i %6i %6i %6i %6i %6i\n',[NSJ2; NALL]);
fprintf(fid,FJ2,transpose(NNdataJ2));
fclose(fid);
end

if NSX2>0
fid = fopen([SDIRO 'NN_TEST_X.dat'],'w');
fprintf(fid,'%6i %6i %6i %6i %6i %6i %6i %6i %6i\n',[NSX2; NALL]);
fprintf(fid,FX2,transpose(NNdataX2));
fclose(fid);
end


