function iflag=MacACADprepareMAST(SDIR1,SDIR2,SDIR3,STASK)

global Acad

% read in EXPEQ-file from 24460@330 (CDN) and 25075@450 (LSN)
% continous scan of shapes from CDN to LSN, with parameter sfac
% sfac=0 ==> 24460; sfac=1 ==> 25075
% keep all other profiles the same as 24460
% scan Ip to make sure the same q95 while scanning sfac
% assume even parity coil phasing for MAST coils 
% prepare RmarsQ-file for MARS-* runs on STASK='RMP' or 'MHD' 
% SDIR1 = input folder where the original files are stored, existence assumed
% SDIR2 = output folder to store the new equilibrium file, existence not assumed
% SDIR3 = directory to run CHEASE and MARS-F, existence not assumed  

%% prepare folders for CHEASE run
mkdir(SDIR2);
mkdir(SDIR3);

copyfile('MacCOMMON/Rchease',[SDIR3 'Rchease'],'f')

SDIRM = pwd;  %%current directory to run Matlab 
SDIRM = [SDIRM '/'];

%% read in BOUNDP for 24460 and 25075 and do combination 
cd(SDIR1);

%%all data below are for 24460
Acad.R0EXP = 0.8740; 
Acad.B0EXP = 0.4661; 
Acad.I0EXP = 0.8423; 
Acad.QZEROG = 1.3533;
Acad.QPSI95G = 5.0412;
Acad.QEDGEG = 8.9923;
Acad.ASPECT = 1.4285;
Acad.NE0 = 4.6602e+19;
Acad.TI0 = 1.3975e+3;
Acad.TE0 = 1.1788e+3;
Acad.ROT0P =5.1939e+4;
Acad.WE0P = Acad.ROT0P;

d1=load('BOUNDP_24460')*Acad.R0EXP;
d2=load('BOUNDP_25075')*0.8606;
d2(:,2) = d2(:,2)-min(d2(:,2)) + min(d1(:,2));

d2(:,2)=-d2(:,2);  %to make up-down mirroring of the shape
t1 = atan2(d1(:,2),d1(:,1)-Acad.R0EXP);
[t1,II]=sort(t1); R1=d1(II,1); Z1=d1(II,2);

t2 = atan2(d2(:,2),d2(:,1)-Acad.R0EXP);
[t2,II]=sort(t2); R2=d2(II,1); Z2=d2(II,2);
R2 = spline(t2,R2,t1);
Z2 = spline(t2,Z2,t1);

ccfac = linspace(0,1,6);
SDIR2_INIT = SDIR2;

for kf=1:length(ccfac)
    cfac = ccfac(kf);
    sfac = sprintf('%4.2f',cfac);
    Acad.DRSEP =cfac;

    RN = ((1-cfac)*R1 + cfac*R2)/Acad.R0EXP;
    ZN = ((1-cfac)*Z1 + cfac*Z2)/Acad.R0EXP;

    SDIR2=[SDIR2_INIT 'cfac_' sfac '/'];

    mkdir(SDIR2);

    %%process EXPEQ_MAST24460 data 
    fid = fopen([SDIR1 'EXPEQ_MAST24460'],'r');
    d   = fscanf(fid,'%f');
    fclose(fid);
    NB = floor(d(4));
    B  = d(7:7+4*NB-1);
    B  = transpose(reshape(B,2,2*NB));
    B(1:NB,1) = RN; 
    B(1:NB,2) = ZN; 

    %%save new EXPEQ file
    fid = fopen([SDIR2 'EXPEQ_SAVE'],'w');
    fprintf(fid,'%16.8e\n',d(1:3));
    fprintf(fid,'%4i %3i %3i\n',[NB; 2; 2]);
    fprintf(fid,'%16.8e %16.8e\n',transpose(B));
    N  = 7+4*NB;
    NP = floor(d(N));
    fprintf(fid,'%4i\n',d(N:N+1));
    fprintf(fid,'%16.8e\n',d(N+2:N+2+NP*3-1));
    fclose(fid);  

    cd(SDIR3);

    %%repare Rchease-file and run it
    copyfile([SDIRM 'MacCOMMON/Rchease'],'Rchease','f');

    ss = ['NTOR=' int2str(Acad.n)];
    eval(['!sed -i s/''NTOR=1''/''' ss '''/g Rchease']);

    s1 = sprintf('%6.4f',Acad.R0EXP);
    s2 = ['R0EXP=' s1];
    eval(['!sed -i s/''R0EXP=1.6955''/''' s2 '''/g Rchease']);      

    s1 = sprintf('%6.4f',Acad.B0EXP);
    s2 = ['B0EXP=' s1];
    eval(['!sed -i s/''B0EXP=1.9000''/''' s2 '''/g Rchease']);  

    if 1==1 %run CHEASE with fixed q95
       s2 = ['QSPEC =' sprintf('%6.4f',Acad.QPSI95G)];
       eval(['!sed -i s/''QSPEC =1.140''/''' s2 '''/g Rchease']);
       eval(['!sed -i s/''CSSPEC=0.000''/''CSSPEC=0.97468''/g Rchease']);
       eval(['!sed -i s/''NCSCAL = 4''/''NCSCAL = 1''/g Rchease']);
    end

    copyfile([SDIR2 'EXPEQ_SAVE'],'EXPEQ','f');
    copyfile('Rchease',[SDIR2 'Rchease_SAVE'],'f');  
    disp('Running CHEASE with EXPEQ-file ...')
    eval('!./Rchease');

    %%get global parameters from log_chease
    eval('!grep Q_ZERO log_chease > tmp.txt')
    eval('!grep Q_EDGE log_chease >> tmp.txt')
    eval('!grep ''Q AT 95% FLUX SURFACE'' log_chease >> tmp.txt')
    eval('!grep GEXP log_chease >> tmp.txt')
    eval('!grep ''RW='' log_chease >> tmp.txt')
    fid = fopen('tmp.txt','r');
    ss  = fgetl(fid);  Acad.QZERO  = str2num(ss(5:18));
    ss  = fgetl(fid);  Acad.QEDGE  = str2num(ss(5:18));
    ss  = fgetl(fid);  Acad.QPSI95 = str2num(ss(26:end));
    ss  = fgetl(fid);  Acad.BETAN  = str2num(ss(28:end));
    ss  = fgetl(fid);  Acad.NW     = floor(str2num(ss(51:end)));
    fclose(fid);

    %%get dimensionless parameters
    cd(SDIRM)
    MacACADquantities;
    
    Acad.ROT0 = Acad.ROT0P*Acad.TAUA;
    Acad.WE0  = Acad.WE0P*Acad.TAUA;
    Acad.ALPHAP = Acad.TI0/(Acad.TI0+Acad.TE0);

    %%prepare RmarsQ-file
    cd(SDIR3)
    copyfile([SDIRM 'MacCOMMON/RmarsQ'],'RmarsQ','f');
    
    %%modify common portion
    %%valid for both STASK='MHD' and 'RMP'
    ss = ['ETA    =' sprintf('%11.4e',Acad.ETA)];
    eval(['!sed -i s/''ETA    = 1.0000e-10''/''' ss '''/g RmarsQ']);

    ss = ['RNTOR  = -' sprintf('%1i',Acad.n)];
    eval(['!sed -i s/''RNTOR  = -1''/''' ss '''/g RmarsQ']);

    eval(['!sed -i s/''NWALL  = 1''/''NWALL  = 0''/g RmarsQ']);

    ss = ['IFEED  = ' sprintf('%3i',Acad.NW) ','];
    eval(['!sed -i s/''IFEED  = 47,''/''' ss '''/g RmarsQ']);

    ss = ['ISENS  = ' sprintf('%3i',Acad.NW) ','];
    eval(['!sed -i s/''ISENS  = 47,''/''' ss '''/g RmarsQ']);

    ss = ['ROTE   =' sprintf('%11.4e',Acad.ROT0)];
    eval(['!sed -i s/''ROTE   = 1.0000e-02''/''' ss '''/g RmarsQ']);

    ss = ['ROTWE0 =' sprintf('%11.4e',Acad.WE0)];
    eval(['!sed -i s/''ROTWE0 = 1.0000e-02''/''' ss '''/g RmarsQ']);

    ss = ['ALPHAP =' sprintf('%11.4e',Acad.ALPHAP)];
    eval(['!sed -i s/''ALPHAP = 5.0000e-01''/''' ss '''/g RmarsQ']);

    ss = ['OMEGACI0 =' sprintf('%11.4e',Acad.OMEGACI0)];
    eval(['!sed -i s/''OMEGACI0 = 1.0000e+02''/''' ss '''/g RmarsQ']);

    if STASK=='MHD'
       copyfile('RmarsQ',[SDIR2 'RmarsQ_SAVE'],'f');
    end

    if STASK=='RMP'
       eval(['!sed -i s/''NCASE  = 1''/''NCASE  = 2''/g RmarsQ']);
       eval(['!sed -i s/''NITMAX = 100,''/''NITMAX = 1,''/g RmarsQ']);
       eval(['!sed -i s/''M1     =-9''/''M1     =-1''/g RmarsQ']);
       eval(['!sed -i s/''M2     = 33''/''M2     = 1''/g RmarsQ']);
       eval(['!sed -i s/''INCFEED= 0''/''INCFEED= 4''/g RmarsQ']);
       eval(['!sed -i s/''PVISC  = 0.1''/''PVISC  = 1.5''/g RmarsQ']);
       eval(['!sed -i s/''TALPHA1= (1.00000E-01''/''TALPHA1= (1.00000E-10''/g RmarsQ']);

       copyfile([SDIR1 'PROFDEN_MAST24460'],'PROFDEN.IN','f');
       copyfile([SDIR1 'PROFROT_MAST24460'],'PROFROT.IN','f');
       copyfile([SDIR1 'PROFTI_MAST24460'],'PROFTI.IN','f');
       copyfile([SDIR1 'PROFTE_MAST24460'],'PROFTE.IN','f');
       copyfile([SDIR1 'PROFROT_MAST24460'],'PROFWE.IN','f');
       disp('Running quick MARS-F ...')
       eval('!./RmarsQ');

       %%setup coil geometry and coil current (=1kAt)
       %%for DIII-D I-coils
       copyfile('RMZM_F.OUT',[SDIR2 'RMZM_F_EQAC'],'f');
       cd(SDIRM)
       
       MacACAD_ICOIL(SDIR2)

       Acad.FEEDI = MacACADcoilcurrent(Acad.n,Acad.R0EXP,Acad.B0EXP);

       cd(SDIR3)
       %%setup new RmarsQ file
       eval(['!sed -i s/''M1     =-1''/''M1     =-33''/g RmarsQ']);
       eval(['!sed -i s/''M2     = 1''/''M2     = 33''/g RmarsQ']);
       eval(['!sed -i s/''INCFEED= 4''/''INCFEED= 8''/g RmarsQ']);

       eval(['!sed -i s/''NCOIL  = 2''/''NCOIL  = 1''/g RmarsQ']);
       ss = ['FEEDI  = (' sprintf('%11.4e',Acad.FEEDI)];
       eval(['!sed -i s/''FEEDI  = ( 1.0000e-03''/''' ss '''/g RmarsQ']);

       %%upper I-coil geometry
       copyfile('RmarsQ','RmarsQ_IU','f');
       ss = ['FCCHI  =' sprintf('%9.6f',Acad.FCCHI(1))];
       eval(['!sed -i s/''FCCHI  =-0.200000''/''' ss '''/g RmarsQ_IU']);
       ss = ['FWCHI  =' sprintf('%9.6f',Acad.FWCHI(1))];
       eval(['!sed -i s/''FWCHI  = 0.100000''/''' ss '''/g RmarsQ_IU']);
       copyfile('RmarsQ_IU',[SDIR2 'RmarsQ_IU_SAVE'],'f');

       %%lower I-coil geometry
       copyfile('RmarsQ','RmarsQ_IL','f');
       ss = ['FCCHI  =' sprintf('%9.6f',Acad.FCCHI(2))];
       eval(['!sed -i s/''FCCHI  =-0.200000''/''' ss '''/g RmarsQ_IL']);
       ss = ['FWCHI  =' sprintf('%9.6f',Acad.FWCHI(2))];
       eval(['!sed -i s/''FWCHI  = 0.100000''/''' ss '''/g RmarsQ_IL']);
       copyfile('RmarsQ_IL',[SDIR2 'RmarsQ_IL_SAVE'],'f');
    end

    save([SDIR2 'Acad.mat'],'Acad')
end

cd(SDIRM);
