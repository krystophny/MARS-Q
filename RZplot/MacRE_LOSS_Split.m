%read in RE_LOSS_*_RAW data
%find particles (by numbering) that have been traced, and exclude from FIDIST*.IN file
%thus split FIDSI*.IN into two parts (traced and not traced)
%also patch data with MF-MD<=10

SDIRM = '/home/liuy/Codes/MarsQ/RZplot/';
%SDIRD = '/home/liuy/Work/MAST-U/46943/Data/FIDIST_FO_E0.1/';
%SDIRD = '/home/liuy/Work/KSTAR/Data/30306.7850_gpec/FIDIST_5D/3D/';
SDIRD = '/home/liuy/Work/DIII-D/D3D177028/Data/FIDIST/RootAR/RE_PERTURB0.01_G/';
%SDIRD = '/home/liuy/Work/ITER/ITER_RE/DataAE_HALO/68ms/FIDIST/Root-c/RE_PERTURB0.01/';

PSS = 'P10';

%go through all RE_LOSS_* and find RAW files
cd(SDIRD);
FF = split(ls);
NF = length(FF)-1;

for kf=1:NF
    F = FF{kf};

    if strcmp(F(end-2:end),'RAW')
       disp(F)
       I1  = 9; 
       I2  = length(F); 
       SFI = F(I1:I2-4);

       eval(['!sed -i s/''LOG_MPI....\.OUT:RECORD:''/''' '''/g ' F]);

       d  = load(F);
       f  = load(['../../FIDIST_' PSS '_' SFI '.IN']);
       %f  = load(['FIDIST_' PSS '_' SFI]);

       MF = f(1,1);
       MD = size(d,1);
       if MF-MD <= 200 
          [I1,II]=sort(d(:,1));
          d = d(II,:);
          II = setdiff([1:MF],d(:,1));
          disp(['number of missing data in ' F ' =' int2str(length(II))])

          for k=1:length(II)
              L=II(k);
              d = [d(1:L-1,:); d(L-1,:); d(L:end,:)];
          end
          d = d(:,2:end);

          save(['RE_LOSS_' SFI],'d','-ascii','-double')
          delete(F);

       elseif size(d,1)>0

       I1 = d(:,1); d = d(:,2:end);
       save(['RE_LOSS_' SFI '_1'],'d','-ascii','-double') 

       N  = f(1,1); f = f(2:end,:);
       I2 = setdiff([1:N],I1); 
       I2 = I2(:);
       I2 = flipud(I2);
       N1 = length(I1);
       N2 = length(I2);
       M  = size(f,2);
       f1 = [N1*ones(1,M); f(I1,:)];
       f2 = [N2*ones(1,M); f(I2,:)];
       save(['FIDIST_' PSS '_' SFI '_1'],'f1','-ascii','-double') 
       save(['FIDIST_' PSS '_' SFI '_2'],'f2','-ascii','-double') 

       delete(F);

       else
    
          disp(['No data for ' F])

       end
    end
end

cd(SDIRM);
                         
