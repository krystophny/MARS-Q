% run MARS-F to produce perturbed 3D equilibria
% use input from the database 

global Mac Acad

kstep = 2; %1: quick MARS-F run for generating geometry from EXPEQ
           %2: MARS-F vacuum runs 
           %3: save vacuum results and MARS-F plasma response runs
           %4: save plasma response data

nn    = 1; %=NTOR
           %for each nn, repeat kstep=1-4

kprof = 2; %=1: generate PROF*.IN profiles from reading PROFEQ.OUT from the exisitng database
           %=2: copy files from existing PROF*.IN data

%NB: need to choose kdev for a device
kdev  = 3; %=1: D3D; 2: EAST; 3: KSTAR; 4=AUG

SDIR1 = '/home/liuy/Work/IR_D/Run_SAVE/';
SDIR2 = '/cscratch/liuy/WorkIR_D/';
SDIR4 = '/home/liuy/IRIS/mars-code/RZplot/';
%SDIR5 = '/cscratch/liuy/Database_More/DIII-D/';
SDIR5 = '/cscratch/liuy/Database_New/KSTAR/';
%SDIRP = '/home/liuy/Work/RMPthreshold/AUG/Data/';  %only needed when kprof=2
SDIRP = '/home/liuy/Work/KSTAR/Data/';              %only needed when kprof=2

%KSTAR
%SHOTS = {'11341','11465','16586','18451','19204','19211','19212','21575','22737','25607','28332'};
%RUNS  = 'ABCDEFGHIJK';
%SHOTS = {'30306.4850','30306.7850'};
%RUNS  = 'AB';
SHOTS = {'30306.7850_gpec'};
RUNS  = 'C';

%DIII-D
%SHOTS = {'158076','164330','164361','170005','174689','179328','179358'};
%RUNS  = 'ABCDEFG';
%SHOTS = {'164330'};
%RUNS  = 'D';

%EAST
%SHOTS = {'67578', '69635', '74184', '74198', '94048', '94088'};
%RUNS = 'ABCDEF';

%AUG
%SHOTS = {'34842', '34205', '35307'};
%RUNS = 'ABC';

if kstep==1
   for k=1:length(RUNS)
       SDIRR = [SDIR2 RUNS(k) 'U/'];
       copyfile([SDIR5 SHOTS{k} '/2D/EXPEQ'],[SDIRR 'EXPEQ_SAVE'],'f');
       ss = [SDIR5 SHOTS{k} '/3D/n' int2str(nn) '/'];

       if kprof==1
       d = load([ss 'PROFEQ.OUT']);
       s = d(:,1);
       den = d(:,5);
       rot = d(:,6);
       Ti  = d(:,10);
       Te  = d(:,11);
       RES = [[length(s) 1]; s den];
       save([SDIRR 'PROFDEN_SAVE'],'RES','-ascii');
       RES = [[length(s) 1]; s rot];
       save([SDIRR 'PROFROT_SAVE'],'RES','-ascii');
       RES = [[length(s) 1]; s Ti];
       save([SDIRR 'PROFTI_SAVE'],'RES','-ascii');
       RES = [[length(s) 1]; s Te];
       save([SDIRR 'PROFTE_SAVE'],'RES','-ascii');
       else
       copyfile([SDIRP SHOTS{k} '/PROFDEN_SAVE.IN'],[SDIRR 'PROFDEN_SAVE'],'f');
       copyfile([SDIRP SHOTS{k} '/PROFROT_SAVE.IN'],[SDIRR 'PROFROT_SAVE'],'f');
       copyfile([SDIRP SHOTS{k} '/PROFTI_SAVE.IN'], [SDIRR 'PROFTI_SAVE'], 'f');
       copyfile([SDIRP SHOTS{k} '/PROFTE_SAVE.IN'], [SDIRR 'PROFTE_SAVE'], 'f');
       end
       
       copyfile([SDIR1 'RmarsQ_quick'],[SDIRR 'RmarsQ_quick'],'f');
       copyfile([SDIR1 'all2.bat'],[SDIRR 'all.bat'],'f');
       copyfile([SDIR1 'runit.bat'],[SDIRR 'runit.bat'],'f');

       d  = load([SDIR5 SHOTS{k} '/2D/global.txt']);

       %prepare Rchease in *U/
       Schease = [SDIRR 'Rchease'];
       copyfile([SDIR1 'Rchease'], Schease, 'f');
       eval(['!sed -i s/"NTOR=1"/"NTOR=' int2str(nn) '"/g ' Schease])
       ss = sprintf('%10.4e',d(2));
       eval(['!sed -i s/"B0EXP=2.0000e-00"/"B0EXP=' ss '"/g ' Schease])
       ss = sprintf('%10.4e',d(1));
       eval(['!sed -i s/"R0EXP=1.8000e-00"/"R0EXP=' ss '"/g ' Schease])

       cd(SDIRR);
       eval('!sbatch runit.bat')
   end
elseif kstep==2
   for k=1:length(RUNS)
       copyfile([SDIR1 'all4.bat'],[SDIR2 RUNS(k) 'U/all.bat'],'f');
       copyfile([SDIR1 'all4.bat'],[SDIR2 RUNS(k) 'L/all.bat'],'f');
       copyfile([SDIR1 'runit.bat'],[SDIR2 RUNS(k) 'L/.'],'f');
       if kdev==3
       copyfile([SDIR1 'all4.bat'],[SDIR2 RUNS(k) 'M/all.bat'],'f');
       copyfile([SDIR1 'runit.bat'],[SDIR2 RUNS(k) 'M/.'],'f');
       end
       if 1==1
       cd([SDIR2 RUNS(k) 'L']);
       eval(['!ln -s ../' RUNS(k) 'U/OUT*MAR .'])
       eval(['!ln -s ../' RUNS(k) 'U/*_SAVE .'])
       end
       if 1==1 & kdev==3
       cd([SDIR2 RUNS(k) 'M']);
       eval(['!ln -s ../' RUNS(k) 'U/OUT*MAR .'])
       eval(['!ln -s ../' RUNS(k) 'U/*_SAVE .'])
       end

       d  = load([SDIR5 SHOTS{k} '/2D/global.txt']);

       %prepare Rmars in *U 
       Smars = [SDIR2 RUNS(k) 'U/RmarsQ'];
       copyfile([SDIR1 'RmarsQ'], Smars, 'f');
       eval(['!sed -i s/"RNTOR  =-1"/"RNTOR  =-' int2str(nn) '"/g ' Smars])
       ss = sprintf('%10.4e',d(10));
       eval(['!sed -i s/"ETA    = 1.0000e-08"/"ETA    = ' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',d(11));
       eval(['!sed -i s/"ROTE   = 1.0000e-02"/"ROTE   =' ss '"/g ' Smars]) 
       ss = sprintf('%10.4e',d(12));
       eval(['!sed -i s/"FEEDI  = (1.0000e-04"/"FEEDI  = (' ss '"/g ' Smars])

       copyfile(Smars, [SDIR2 RUNS(k) 'L/.'], 'f');
       if kdev==3, copyfile(Smars, [SDIR2 RUNS(k) 'M/.'], 'f'); end

       cd(SDIR4)
       SDIR = [SDIR2 RUNS(k) 'U/'];
       if kdev==1, Mac.coilN  = [2.184 1.012 2.394 0.504];     end %D3D Upper I-coils
       if kdev==2, Mac.coilN  = [2.092 0.759 2.278 0.577];     end %EAST Upper coils
       if kdev==3, Mac.coilN = [2.138 0.942 2.518 0.550];      end %KSTAR Upper coils
       if kdev==4, Mac.coilN  = [1.9423 0.8490 2.0982 0.5846]; end %AUG Upper coils
       MacMainKSTAR_G
       ss = sprintf('%11.4e',Acad.FCCHI);
       eval(['!sed -i s/"FCCHI  = 2.0000e-01"/"FCCHI  =' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',Acad.FWCHI);
       eval(['!sed -i s/"FWCHI  = 2.0000e-02"/"FWCHI  =' ss '"/g ' Smars]) 
       eval(['!sed -i s/"IFEED  = 50"/"IFEED  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       eval(['!sed -i s/"ISENS  = 50"/"ISENS  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
  
       %prepare Rmars in *M 
       if kdev==3
       cd(SDIR4)
       Smars = [SDIR2 RUNS(k) 'M/RmarsQ'];
       Mac.coilN = [2.482 0.550 2.482 -0.550]; %KSTAR Middle coils
       MacMainKSTAR_G
       ss = sprintf('%11.4e',Acad.FCCHI);
       eval(['!sed -i s/"FCCHI  = 2.0000e-01"/"FCCHI  =' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',Acad.FWCHI);
       eval(['!sed -i s/"FWCHI  = 2.0000e-02"/"FWCHI  =' ss '"/g ' Smars]) 
       eval(['!sed -i s/"IFEED  = 50"/"IFEED  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       eval(['!sed -i s/"ISENS  = 50"/"ISENS  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       end

       %prepare Rmars in *L 
       cd(SDIR4)
       Smars = [SDIR2 RUNS(k) 'L/RmarsQ'];
       if kdev==1, Mac.coilN  = [2.394 -0.504 2.184 -1.012];     end %D3D Lower I-coils
       if kdev==2, Mac.coilN  = [2.278 -0.577 2.092 -0.759];     end %EAST Lower I-coils
       if kdev==3, Mac.coilN  = [2.138 -0.942 2.518 -0.550];     end %KSTAR Lower coils
       if kdev==4, Mac.coilN  = [2.0102 -0.5867 1.8509 -0.8291]; end %AUG Lower coils
       MacMainKSTAR_G
       ss = sprintf('%11.4e',Acad.FCCHI);
       eval(['!sed -i s/"FCCHI  = 2.0000e-01"/"FCCHI  =' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',Acad.FWCHI);
       eval(['!sed -i s/"FWCHI  = 2.0000e-02"/"FWCHI  =' ss '"/g ' Smars]) 
       eval(['!sed -i s/"IFEED  = 50"/"IFEED  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       eval(['!sed -i s/"ISENS  = 50"/"ISENS  = ' int2str(Acad.IFEED) '"/g ' Smars]) 

       cd([SDIR2 RUNS(k) 'U']);
       eval('!sbatch runit.bat')

       if kdev==3
       cd([SDIR2 RUNS(k) 'M']);
       eval('!sbatch runit.bat')
       end

       cd([SDIR2 RUNS(k) 'L']);
       eval('!sbatch runit.bat')
   end
elseif kstep==3
   for k=1:length(RUNS)
       ss = [SDIR5 SHOTS{k} '/3D/n' int2str(nn) '/'];
       mkdir([ss 'U'])
       if kdev==3, mkdir([ss 'M']); end
       mkdir([ss 'L'])

       %save vacuum data
       copyfile([SDIR2 RUNS(k) 'U/PROFEQ.OUT'], [ss '.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/RMZM_F.OUT'], [ss '.'], 'f')
       cd(SDIR4)

       SDIR = [SDIR2 RUNS(k) 'U/'];
       MacMainIR_D 
       copyfile([SDIR4 'dBnormal.txt'], [ss 'U/.'], 'f')      

       if kdev==3
       SDIR = [SDIR2 RUNS(k) 'M/'];
       MacMainIR_D 
       copyfile([SDIR4 'dBnormal.txt'], [ss 'M/.'], 'f')      
       end

       SDIR = [SDIR2 RUNS(k) 'L/'];
       MacMainIR_D 
       copyfile([SDIR4 'dBnormal.txt'], [ss 'L/.'], 'f')      

       %prepare and run MARS-F plasma response
       Smars = [SDIR2 RUNS(k) 'U/RmarsQ'];
       eval(['!sed -i s/"INCFEED= 4"/"INCFEED= 8"/g ' Smars]) 
       cd([SDIR2 RUNS(k) 'U']);
       eval('!sbatch runit.bat')

       if kdev==3           
       Smars = [SDIR2 RUNS(k) 'M/RmarsQ'];
       eval(['!sed -i s/"INCFEED= 4"/"INCFEED= 8"/g ' Smars]) 
       cd([SDIR2 RUNS(k) 'M']);
       eval('!sbatch runit.bat')
       end
           
       Smars = [SDIR2 RUNS(k) 'L/RmarsQ'];
       eval(['!sed -i s/"INCFEED= 4"/"INCFEED= 8"/g ' Smars]) 
       cd([SDIR2 RUNS(k) 'L']);
       eval('!sbatch runit.bat')
   end
elseif kstep==4
   for k=1:length(RUNS)
       ss = [SDIR5 SHOTS{k} '/3D/n' int2str(nn) '/'];
       copyfile([SDIR2 RUNS(k) 'U/BPLASMA.OUT'], [ss 'U/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/JPLASMA.OUT'], [ss 'U/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/XPLASMA.OUT'], [ss 'U/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/PPLASMA.OUT'], [ss 'U/.'], 'f')

       if kdev==3
       copyfile([SDIR2 RUNS(k) 'M/BPLASMA.OUT'], [ss 'M/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'M/JPLASMA.OUT'], [ss 'M/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'M/XPLASMA.OUT'], [ss 'M/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'M/PPLASMA.OUT'], [ss 'M/.'], 'f')
       end

       copyfile([SDIR2 RUNS(k) 'L/BPLASMA.OUT'], [ss 'L/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'L/JPLASMA.OUT'], [ss 'L/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'L/XPLASMA.OUT'], [ss 'L/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'L/PPLASMA.OUT'], [ss 'L/.'], 'f')
   end
end

cd(SDIR4)

