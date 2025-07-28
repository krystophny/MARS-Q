% select and form NN database for 3D equilibria
% according to the input/output format defined in drawMemo.m 
% use V-values from global SVD (for whole database) as NN output
% database is expanded to include linear superposition of different rows

SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';
SDIR1 = '/home/liuy/Work/IR_D/Database/'; 
SDIRO = '/home/liuy/Work/IR_D/DataNN/Raw1/'; 

%DEV = {'MAST','DIII-D','AUG','ITER'};
%DEV = {'DIII-D','AUG','ITER'};
DEV = {'MAST'};
%SN  = '1234';
SN = '4';
koutput =0;  % 0:   select all V as output (as stored in database point)
             % 1-5: select first koutput V-values as output
             %(-1)-(-5) : select |koutput| V value as output
rfrac = 0.8;  %fraction of random selection for training data

if length(SN)==1
   if koutput<0
      SDIRO = [SDIRO 'n' SN '_S' int2str(-koutput) '_P/'];
      mkdir(SDIRO)
   else 
      SDIRO = [SDIRO 'n' SN '/'];
   end
   SOUT = ['n' SN];
else
   SDIRO = [SDIRO 'nn/'];
   SOUT = 'nn';
end

mkdir(SDIRO);

NNdataJ = [];
NNdataX = [];
for kd=1:size(DEV,2)
    SDIR2 = [SDIR1 DEV{kd} '/'];
    cd(SDIR2)
    F  = split(ls);
    NF = length(F)-1

    for kf=1:NF
        FF    = F{kf};
        SDIR3 = [SDIR2 FF '/3D/'];
        din1 = load([SDIR3 'NNdata_INPUT2D.dat']);
        NRZ  = round(din1(2));
        NPR  = round(din1(3));

        for kn=1:length(SN)
        SDIR4 = [SDIR3 'n' SN(kn) '/'];
        cd(SDIR4)
        C = split(ls); SC={};
        NC = length(C)-1;
        kc1 = 0;
        for kc=1:NC
            if length(C{kc})<=2
               kc1 = kc1+1;
               SC{kc1}=C{kc};
            end
        end
        NC1 = kc1;
        for kc=1:NC1
            SDIR5 = [SDIR4 SC{kc} '/'];
            
            din2 = load([SDIR5 'NNdata_INPUT3D.dat']);
            din  = [din1(4:end); din2(3:end)];
            if length(SN)>1, din = [str2num(SN(kn)); din]; end
            dout = load([SDIR5 'SVD1_Jpara_' SOUT '.txt']);
            NBV  = round(din2(1));
            NOU  = length(dout);
            NSS  = round(NOU/2);
            NIN  = length(din);

            d = [din; dout];
            NNdataJ = [NNdataJ d];

            dout = load([SDIR5 'SVD1_X1U_' SOUT '.txt']);
            d = [din; dout];
           NNdataX = [NNdataX d];
         end
         end
    end
end

NUS = 0;
NVM = 67;

%take log10(|V|) and angle(V) as NN output
if 1==1
   NOU2 = round(NOU/2); NOUb = NOU; NOU2b = NOU2;
   d = NNdataJ(NIN+1:NIN+NOU2,:) + NNdataJ(NIN+NOU2+1:end,:)*i;
   if koutput>0, d=d(1:koutput,:); NOUb=koutput*2; NOU2b=koutput; end
   if koutput<0, d=d(-koutput,:); NOUb=2; NOU2b=1; end
   %if koutput<0, d=d(-koutput,:); NOUb=1; NOU2b=1; end
   a = angle(d); %II=find(a>0); a(II)=a(II)-2*pi;
   %a = tan(a/2.5);
   %d = [log10(abs(d)); a];    
   %d = [a]; 
   d = [real(d); imag(d)]*1e+0;
   NNdataJ=[NNdataJ(1:NIN,:); d];

   d = NNdataX(NIN+1:NIN+NOU2,:) + NNdataX(NIN+NOU2+1:end,:)*i;
   if koutput>0, d=d(1:koutput,:); NOUb=koutput*2; NOU2b=koutput; end
   if koutput<0, d=d(-koutput,:); NOUb=2; NOU2b=1; end
   %if koutput<0, d=d(-koutput,:); NOUb=1; NOU2b=1; end
   a = angle(d); %II=find(a>0); a(II)=a(II)-2*pi;
   %a = tan(a/2.5);
   %d = [log10(abs(d)); a]; 
   %d = [a]; 
   d = [real(d); imag(d)]*1e+0;
   NNdataX=[NNdataX(1:NIN,:); d];

   NOU = NOUb; NOU2 = NOU2b; NSS=NOU2b;
end

if NOU>NOU2
hf=figure(1);
plot(NNdataJ(NIN+1:NIN+NOU2,:),NNdataJ(NIN+NOU2+1:end,:),'bo'), hold on
xlabel('log10(|V|) for J_{||}')
ylabel('phase(V) for J_{||}')

hf=figure(2);
plot(NNdataX(NIN+1:NIN+NOU2,:),NNdataX(NIN+NOU2+1:end,:),'bo'), hold on
xlabel('log10(|V|) for \xi^1')
ylabel('phase(V) for \xi^1')
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

disp(['number of training data points = ',int2str([NSJ1 NSX1])]),
disp(['number of testing data points = ',int2str([NSJ2 NSX2])]),

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


cd(SDIRM)
