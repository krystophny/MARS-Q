% run MARS-F to produce perturbed 3D equilibria
%
global Mac Acad

kstep = 5; %1: run CHEASE with g-file as input
           % :  generate EXPEQ with drawBound.m etc.
           %2: quick MARS-F run for generating geometry from EXPEQ
           % : prepare PROF*_SAVE
           % : prepare 2D/global.txt
           %3: MARS-F vacuum runs 
           %4: save vacuum results and MARS-F plasma response runs
           %5: save plasma response data

nn    = 4; %=NTOR
           %for each nn, repeat kstep=2-5

SDIR1 = '/home/liuy/Work/IR_D/Run_SAVE/';
SDIR2 = '/cscratch/liuy/WorkIR_D/';
SDIR3 = '/home/liuy/Work/IR_D/Work/';
SDIR4 = '/home/liuy/IRIS/mars-code/RZplot/';
SDIR5 = '/home/liuy/Work/DIII-D/';

%KSTAR
%SHOTS = {'11341','11465','16586','18451','19204','19211','19212','21575','22737','25607','28332'};
%RUNS  = 'ABCDEFGHIJK';

%DIII-D
SHOTS = {'158076','164330','164361','170005','174689','179328','179358'};
RUNS  = 'ABCDEFG';

if kstep==1
   for k=1:length(RUNS)
       mkdir([SDIR2 RUNS(k) 'U']);
       mkdir([SDIR2 RUNS(k) 'L']);
       mkdir([SDIR2 RUNS(k) 'M']);

       copyfile([SDIR5 SHOTS{k} '/2D/gfile'],[SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR1 'all1.bat'],[SDIR2 RUNS(k) 'U/all.bat'],'f');
       copyfile([SDIR1 'runit.bat'],[SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR1 'Rchease_g'],[SDIR2 RUNS(k) 'U/.'],'f');
       cd([SDIR2 RUNS(k) 'U']);
       eval('!sbatch runit.bat')
   end
elseif kstep==2
   for k=1:length(RUNS)
       copyfile([SDIR5 SHOTS{k} '/2D/EXPEQ'],[SDIR2 RUNS(k) 'U/EXPEQ_SAVE'],'f');
       copyfile([SDIR1 'RmarsQ_quick'],[SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR1 'all2.bat'],[SDIR2 RUNS(k) 'U/all.bat'],'f');

       d  = load([SDIR5 SHOTS{k} '/2D/global.txt']);

       %prepare Rchease in *U/
       Schease = [SDIR2 RUNS(k) 'U/Rchease'];
       copyfile([SDIR1 'Rchease'], Schease, 'f');
       eval(['!sed -i s/"NTOR=1"/"NTOR=' int2str(nn) '"/g ' Schease])
       ss = sprintf('%10.4e',d(2));
       eval(['!sed -i s/"B0EXP=2.0000e-00"/"B0EXP=' ss '"/g ' Schease])
       ss = sprintf('%10.4e',d(1));
       eval(['!sed -i s/"R0EXP=1.8000e-00"/"R0EXP=' ss '"/g ' Schease])

       cd([SDIR2 RUNS(k) 'U']);
       eval('!sbatch runit.bat')
   end
elseif kstep==3
   for k=1:length(RUNS)
       copyfile([SDIR5 SHOTS{k} '/PROFDEN_SAVE'],[SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR5 SHOTS{k} '/PROFROT_SAVE'],[SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR5 SHOTS{k} '/PROFTI_SAVE'], [SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR5 SHOTS{k} '/PROFTE_SAVE'], [SDIR2 RUNS(k) 'U/.'],'f');
       copyfile([SDIR1 'all4.bat'],[SDIR2 RUNS(k) 'U/all.bat'],'f');
       copyfile([SDIR1 'all4.bat'],[SDIR2 RUNS(k) 'M/all.bat'],'f');
       copyfile([SDIR1 'all4.bat'],[SDIR2 RUNS(k) 'L/all.bat'],'f');
       copyfile([SDIR1 'runit.bat'],[SDIR2 RUNS(k) 'M/.'],'f');
       copyfile([SDIR1 'runit.bat'],[SDIR2 RUNS(k) 'L/.'],'f');
       if 1==0
       cd([SDIR2 RUNS(k) 'M']);
       eval(['!ln -s ../' RUNS(k) 'U/OUT*MAR .'])
       eval(['!ln -s ../' RUNS(k) 'U/*_SAVE .'])
       cd([SDIR2 RUNS(k) 'L']);
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

       copyfile(Smars, [SDIR2 RUNS(k) 'M/.'], 'f');
       copyfile(Smars, [SDIR2 RUNS(k) 'L/.'], 'f');

       cd(SDIR4)
       SDIR = [SDIR2 RUNS(k) 'U/'];
       %Mac.coilN = [2.138 0.942 2.518 0.550]; %KSTAR Upper coils
       Mac.coilN  = [2.184 1.012 2.394 0.504]; %D3D Upper I-coil
       MacMainKSTAR_G
       ss = sprintf('%11.4e',Acad.FCCHI);
       eval(['!sed -i s/"FCCHI  = 2.0000e-01"/"FCCHI  =' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',Acad.FWCHI);
       eval(['!sed -i s/"FWCHI  = 2.0000e-02"/"FWCHI  =' ss '"/g ' Smars]) 
       eval(['!sed -i s/"IFEED  = 50"/"IFEED  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       eval(['!sed -i s/"ISENS  = 50"/"ISENS  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
  
       %prepare Rmars in *L 
       cd(SDIR4)
       Smars = [SDIR2 RUNS(k) 'L/RmarsQ'];
       %Mac.coilN = [2.138 -0.942 2.518 -0.550]; %KSTAR Lower coils
       Mac.coilN  = [2.394 -0.504 2.184 -1.012]; %D3D Lower I-coils
       MacMainKSTAR_G
       ss = sprintf('%11.4e',Acad.FCCHI);
       eval(['!sed -i s/"FCCHI  = 2.0000e-01"/"FCCHI  =' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',Acad.FWCHI);
       eval(['!sed -i s/"FWCHI  = 2.0000e-02"/"FWCHI  =' ss '"/g ' Smars]) 
       eval(['!sed -i s/"IFEED  = 50"/"IFEED  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       eval(['!sed -i s/"ISENS  = 50"/"ISENS  = ' int2str(Acad.IFEED) '"/g ' Smars]) 

       %prepare Rmars in *M 
       cd(SDIR4)
       Smars = [SDIR2 RUNS(k) 'M/RmarsQ'];
       %Mac.coilN = [2.482 0.550 2.482 -0.550]; %KSTAR Middle coils
       Mac.coilN  = [3.23 0.80 3.23 -0.80]; %D3D C-coils
       MacMainKSTAR_G
       ss = sprintf('%11.4e',Acad.FCCHI);
       eval(['!sed -i s/"FCCHI  = 2.0000e-01"/"FCCHI  =' ss '"/g ' Smars]) 
       ss = sprintf('%11.4e',Acad.FWCHI);
       eval(['!sed -i s/"FWCHI  = 2.0000e-02"/"FWCHI  =' ss '"/g ' Smars]) 
       eval(['!sed -i s/"IFEED  = 50"/"IFEED  = ' int2str(Acad.IFEED) '"/g ' Smars]) 
       eval(['!sed -i s/"ISENS  = 50"/"ISENS  = ' int2str(Acad.IFEED) '"/g ' Smars]) 

       cd([SDIR2 RUNS(k) 'U']);
       eval('!sbatch runit.bat')

       cd([SDIR2 RUNS(k) 'L']);
       eval('!sbatch runit.bat')
  
       cd([SDIR2 RUNS(k) 'M']);
       eval('!sbatch runit.bat')
   end
elseif kstep==4
   for k=1:length(RUNS)
       %mkdir([SDIR5 SHOTS{k} '/3D'])
       ss = [SDIR5 SHOTS{k} '/3D/n' int2str(nn) '/'];
       mkdir(ss)
       mkdir([ss 'U'])
       mkdir([ss 'L'])
       mkdir([ss 'C'])

       %save vacuum data
       copyfile([SDIR2 RUNS(k) 'U/PROFEQ.OUT'], [ss '.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/RMZM_F.OUT'], [ss '.'], 'f')
       cd(SDIR4)

       SDIR = [SDIR2 RUNS(k) 'U/'];
       MacMainIR_D 
       copyfile([SDIR4 'dBnormal.txt'], [ss 'U/.'], 'f')      

       SDIR = [SDIR2 RUNS(k) 'L/'];
       MacMainIR_D 
       copyfile([SDIR4 'dBnormal.txt'], [ss 'L/.'], 'f')      

       SDIR = [SDIR2 RUNS(k) 'M/'];
       MacMainIR_D 
       copyfile([SDIR4 'dBnormal.txt'], [ss 'C/.'], 'f')  

       %prepare and run MARS-F plasma response
       Smars = [SDIR2 RUNS(k) 'U/RmarsQ'];
       eval(['!sed -i s/"INCFEED= 4"/"INCFEED= 8"/g ' Smars]) 
       cd([SDIR2 RUNS(k) 'U']);
       eval('!sbatch runit.bat')
           
       Smars = [SDIR2 RUNS(k) 'L/RmarsQ'];
       eval(['!sed -i s/"INCFEED= 4"/"INCFEED= 8"/g ' Smars]) 
       cd([SDIR2 RUNS(k) 'L']);
       eval('!sbatch runit.bat')
           
       Smars = [SDIR2 RUNS(k) 'M/RmarsQ'];
       eval(['!sed -i s/"INCFEED= 4"/"INCFEED= 8"/g ' Smars]) 
       cd([SDIR2 RUNS(k) 'M']);
       eval('!sbatch runit.bat')
   end
elseif kstep==5
   for k=1:length(RUNS)
       ss = [SDIR5 SHOTS{k} '/3D/n' int2str(nn) '/'];
       copyfile([SDIR2 RUNS(k) 'U/BPLASMA.OUT'], [ss 'U/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/JPLASMA.OUT'], [ss 'U/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/XPLASMA.OUT'], [ss 'U/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'U/PPLASMA.OUT'], [ss 'U/.'], 'f')

       copyfile([SDIR2 RUNS(k) 'L/BPLASMA.OUT'], [ss 'L/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'L/JPLASMA.OUT'], [ss 'L/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'L/XPLASMA.OUT'], [ss 'L/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'L/PPLASMA.OUT'], [ss 'L/.'], 'f')

       copyfile([SDIR2 RUNS(k) 'M/BPLASMA.OUT'], [ss 'C/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'M/JPLASMA.OUT'], [ss 'C/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'M/XPLASMA.OUT'], [ss 'C/.'], 'f')
       copyfile([SDIR2 RUNS(k) 'M/PPLASMA.OUT'], [ss 'C/.'], 'f')
   end
end

cd(SDIR4)

