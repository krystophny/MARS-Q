function iflag=MacACADprepare(SDIR1,SDIR2,SDIR3,STASK)

global Acad

% read in g-file and generate the EXPEQ-file
% read in p-file and generate the PROF*.IN files for MARS-*
% prepare RmarsQ-file for MARS-* runs, depending on STASK='RMP' or 'MHD' 
% SDIR1 = input folder where the g-files are stored, existence assumed
% SDIR2 = output folder to store the EXPEQ-file, existence not assumed
% SDIR3 = directory to run CHEASE and MARS-F, existence not assumed  
% NB: if there are multiple g-files in SDIR1, each file will be processed and 
%     the corresponding EXPEQ-file generated

%% prepare folders for CHEASE run
mkdir(SDIR2);
mkdir(SDIR3);

copyfile('MacCOMMON/Rchease_g',[SDIR3 'Rchease_g'],'f')

SDIRM = pwd;  %%current directory to run Matlab 
SDIRM = [SDIRM '/'];

%% read in g-files from SDIR1
cd(SDIR1);
FF = split(ls);
NF = length(FF)-1;

SDIR2_INIT = SDIR2;

for kf=1:NF
    F = FF{kf};

    %%work on g-file
    if F(1)=='g'
       %SDIR2=[SDIR2_INIT F(2:end) '/'];
       SDIR2 = SDIR2_INIT;

       if isfield(Acad,'DRPW')
          ss = sprintf('%3.1f',Acad.DRPW);
          SDIR2 = [SDIR2 'DRPW_' ss '/'];
       end

       mkdir(SDIR2);

       copyfile([SDIR1 F],[SDIR3 'EXPEQ'],'f');
       cd(SDIR3);
       ss = ['NTOR=' int2str(Acad.n)];
       eval(['!sed -i s/''NTOR=1''/''' ss '''/g Rchease_g']);
       disp('Running CHEASE with g-file ...')
       eval('!./Rchease_g');

       %%get global parameters from log_chease
       eval('!grep ''EQDSK: R0EXP'' log_chease > tmp.txt')
       eval('!grep ''EQDSK: B0EXP'' log_chease >> tmp.txt')
       eval('!grep ''EQDSK: I0EXP'' log_chease >> tmp.txt')
       eval('!grep ''EQDSK: QZERO'' log_chease >> tmp.txt')
       eval('!grep ''EQDSK: QPSI95'' log_chease >> tmp.txt')
       eval('!grep ''EQDSK: QEDGE'' log_chease >> tmp.txt')
       eval('!grep ''ASPECT RATIO ;'' log_chease >> tmp.txt');
       fid = fopen('tmp.txt','r');
       ss  = fgetl(fid);  Acad.R0EXP = str2num(ss(20:end));
       ss  = fgetl(fid);  Acad.B0EXP = str2num(ss(20:end));
       ss  = fgetl(fid);  Acad.I0EXP = str2num(ss(20:end));
       ss  = fgetl(fid);  Acad.QZEROG = str2num(ss(20:end));
       ss  = fgetl(fid);  Acad.QPSI95G = str2num(ss(20:end));
       ss  = fgetl(fid);  Acad.QEDGEG = str2num(ss(20:end));
       ss  = fgetl(fid);  Acad.ASPECT = str2num(ss(5:18));
       fclose(fid);

       %%process EXPEQ.OUT data 
       fid = fopen('EXPEQ.OUT','r');
       d   = fscanf(fid,'%f');
       fclose(fid);
       NB = floor(d(4));
       B  = d(7:7+2*NB-1);
       B  = transpose(reshape(B,2,NB));
       save([SDIRM 'BOUNDP'],'B','-ascii','-double')

       cd(SDIRM)
       %%smooth plasma boundary
       B = MacACADbound(Acad.R0EXP);

       N  = 7+2*NB;
       NP = floor(d(N));
       %%modify equilibrium pressure or current profile
       P  = MacACADpj(d(N+2:N+2+NP*3-1));   

       %%save EXPEQ file
       fid = fopen('EXPEQ','w');
       fprintf(fid,'%16.8e\n',d(1:3));
       fprintf(fid,'%4i %3i %3i\n',[255; 2; 2]);
       fprintf(fid,'%16.8e %16.8e\n',transpose(B));
       fprintf(fid,'%4i\n',d(N:N+1));
       fprintf(fid,'%16.8e\n',P);
       fclose(fid);  

       copyfile('EXPEQ',[SDIR2 'EXPEQ_SAVE'],'f');

       %%work on p-file and generate PROF*.IN files 
       cd(SDIRM);
       RES=MacACADprof(SDIR1,SDIR2,['p' F(2:end)]);

       %%collect global parameters for running quantities.m
       Acad.NE0 = RES(1);
       Acad.TI0 = RES(2);
       Acad.TE0 = RES(3);
       Acad.ROT0 = RES(4);
       Acad.WE0 = RES(5);


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

       if 1==1 %run CHEASE with fied qa, need to turn off for general run
          s2 = ['QSPEC =' sprintf('%6.4f',Acad.QPSI95G)];
          eval(['!sed -i s/''QSPEC =1.140''/''' s2 '''/g Rchease']);
          eval(['!sed -i s/''CSSPEC=0.000''/''CSSPEC=0.9747''/g Rchease']);
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
       
       Acad.ROT0 = Acad.ROT0*Acad.TAUA;
       Acad.WE0  = Acad.WE0*Acad.TAUA;
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

       ss = ['IWALL  = ' sprintf('%3i',Acad.NW) ','];
       eval(['!sed -i s/''IWALL  = 48,''/''' ss '''/g RmarsQ']);

       ss = ['IFEED  = ' sprintf('%3i',Acad.NW-1) ','];
       eval(['!sed -i s/''IFEED  = 47,''/''' ss '''/g RmarsQ']);

       ss = ['ISENS  = ' sprintf('%3i',Acad.NW-1) ','];
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
          %ss = ['OMEGACI0 =' sprintf('%11.4e',Acad.OMEGACI0)];
          %eval(['!sed -i s/''OMEGACI0 = 1.0000e+02''/''' ss '''/g RmarsQ']);

          copyfile([SDIR2 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
          copyfile([SDIR2 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
          copyfile([SDIR2 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
          copyfile([SDIR2 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
          copyfile([SDIR2 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
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
end

cd(SDIRM);
