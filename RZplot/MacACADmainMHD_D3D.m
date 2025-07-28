% Automatic Chease&mars-* setup for D3D  (ACAD)
% starting from standard g-file and p-file for DIII-D equilibrium
global Acad

% note that the options below need to be executed in the same order as listed
kaction = 1; %1: prepare MARS-F runs 
             %2: quick MARS-F run to generate PROFEQ.OUT and save it
	     %3: plot all q-profiles, ne, Ti, Te profiles
             %4: setup and run MARS-F 
             %5: save MARS-F run data and plot eigenfunctions
	     %6: setup and run MARS-F with ...
 
Acad.n    = 1;
Acad.DEVICE = 'D3D';
%Acad.DRPW = 0.0;  %[cm]

Acad.ModJ    = 0;
Acad.ModJa1  = 0.65;
Acad.ModJa2  = 1.0 - Acad.ModJa1;
Acad.ModJa3  = 2.0;


SDIRM  = '~/IRIS/mars-code/RZplot/';
SDIR16 = '~/Work/D3D165370/';
SDIR3 = '/cscratch/liuy/WorkTEMP/';
SDIR1 = [SDIR16 'Expt/'];
SDIR2 = [SDIR16 'Data/'];

if kaction==1
   cd(SDIR1)
   FF  = split(ls);
   NF = length(FF);
   for kf=1:NF
       F = FF{kf};
       if length(F)>0 & strcmp(F(1),'g')
       SDIR4 = [SDIR2];

       cd(SDIRM)
       MacACADprepare([SDIR1],SDIR4,SDIR3,'MHD');
       end
   end
end

if kaction==2|kaction==3
   cd(SDIR2)
   FF  = split(ls);
   NF = length(FF);
  
   cd(SDIR3)
   TIMS = [];
   for kf=1:NF
       F = FF{kf};
       if length(F)>0
       SDIR4 = [SDIR2 F '/'];
       TIMS = [TIMS; F(1:6) '.' F(8:11)];

       if kaction==2
       copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
       copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
       eval('!./Rchease');

       copyfile([SDIR4 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
       copyfile([SDIR4 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
       copyfile([SDIR4 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
       copyfile([SDIR4 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
       copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
       copyfile([SDIR4 'RmarsQ_SAVE'],'RmarsQ','f');

       eval(['!sed -i s/''M1     =-9''/''M1     =-1''/g RmarsQ']);
       eval(['!sed -i s/''M2     = 33''/''M2     = 1''/g RmarsQ']);
       eval(['!sed -i s/''INCFEED= 0''/''INCFEED= 4''/g RmarsQ']);
       eval(['!sed -i s/''NITMAX = 100''/''NITMAX = 1''/g RmarsQ']);
       eval(['!sed -i s/''NCASE  = 1''/''NCASE  = 2''/g RmarsQ']);
       eval('!./RmarsQ')

       copyfile('PROFEQ.OUT',[SDIR4 'PROFEQ.OUT'],'f');
       end

       if kaction==3
          SCOL = [(kf-1)/(NF-1) 0 (NF-kf)/(NF-1)];
          load([SDIR4 'Acad.mat']);

          hf=figure(1);
          d = load([SDIR4 'PROFEQ.OUT']);
          plot(d(:,1).^2,d(:,2),'-','LineWidth',3,'Color',SCOL), hold on 
          xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
          ylabel('safety factor','FontSize',18,'FontWeight','Bold'),
          ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

          hf=figure(2);
          d = load([SDIR4 'PROFDEN_SAVE.IN']); d=d(2:end,:);
          plot(d(:,1).^2,d(:,2)*Acad.NE0/1e+19,'-','LineWidth',3,'Color',SCOL), hold on 
          xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
          ylabel('density [x10^{19}m^{-3}]','FontSize',18,'FontWeight','Bold'),
          ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

          hf=figure(3);
          d = load([SDIR4 'PROFTI_SAVE.IN']); d=d(2:end,:);
          plot(d(:,1).^2,d(:,2)/1000,'-','LineWidth',3,'Color',SCOL), hold on 
          xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
          ylabel('Ti [keV]','FontSize',18,'FontWeight','Bold'),
          ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

          hf=figure(4);
          d = load([SDIR4 'PROFTE_SAVE.IN']); d=d(2:end,:);
          plot(d(:,1).^2,d(:,2)/1000,'-','LineWidth',3,'Color',SCOL), hold on 
          xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
          ylabel('Te [keV]','FontSize',18,'FontWeight','Bold'),
          ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
       end
    end
    end
    if kaction==3, 
       figure(1),  plot([0 1],[1 1],'k--'), legend(TIMS), 
       figure(2),  legend(TIMS),        
       figure(3),  legend(TIMS), 
       figure(4),  legend(TIMS), 
     end
end
 
if kaction==4
   cd(SDIR2)
   FF  = split(ls);
   NF = length(FF);
  
   cd(SDIR3)
   TIMS = [];
   for kf=1:NF
       F = FF{kf};
       if length(F)>4
       %if strcmp(F,'179991_3100_7')
       SDIR4 = [SDIR2 F '/'];
   %SDIR4 = [SDIR2 '179824.03230/'];
   %SDIR4 = [SDIR2 '179993_1450_5/'];

   %ss = sprintf('%3.1f',Acad.DRPW);
   %SDIR4 = [SDIR4 'DRPW_' ss '/'];
   
   cd(SDIR3)
   copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
   copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
   copyfile([SDIR4 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
   copyfile([SDIR4 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
   copyfile([SDIR4 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
   copyfile([SDIR4 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
   copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
   copyfile([SDIR4 'RmarsQ_SAVE'],'RmarsQ','f');
   copyfile([SDIRM 'MacCOMMON/runit.bat'],'runit.bat','f');
   copyfile([SDIRM 'MacCOMMON/runit.sbatch'],'runit.sbatch','f');

   load([SDIR4 'Acad.mat']);
   if 1==1 %run MARS-F with zero rotation and fixed ETA
   s1 = ['ROTE   =' sprintf('%11.4e',Acad.ROT0)];
   s2 = ['ROTE   =' sprintf('%11.4e',0)];
   eval(['!sed -i s/''' s1 '''/''' s2 '''/g RmarsQ']);

   s1 = ['ETA    =' sprintf('%11.4e',Acad.ETA)];
   %s2 = ['ETA    =' sprintf('%11.4e',1.2106e-6)];
   s2 = ['ETA    =' sprintf('%11.4e',0)];
   eval(['!sed -i s/''' s1 '''/''' s2 '''/g RmarsQ']);

   eval(['!sed -i s/''PVISC  = 0.1''/''PVISC  = 0.0''/g RmarsQ']);
   eval(['!sed -i s/''NWALL  = 1''/''NWALL  = 0''/g RmarsQ']);
   end

   krunopt = 3;
   if krunopt==1
      eval('!sbatch runit.sbatch')
   elseif krunopt==2 
      copyfile('RMZM_F.OUT',[SDIR4 'RMZM_F_EQAC'],'f');
      copyfile('PROFEQ.OUT',[SDIR4 'PROFEQ.OUT'],'f');
      copyfile('BPLASMA.OUT',[SDIR4 'BPLASMA_.OUT'],'f');
      copyfile('JPLASMA.OUT',[SDIR4 'JPLASMA_.OUT'],'f');
      copyfile('XPLASMA.OUT',[SDIR4 'XPLASMA_.OUT'],'f');
   elseif krunopt==3  %beta scan
      load([SDIR4 'Acad.mat']);
      Acad.n = 3;
      s2 = ['QSPEC =' sprintf('%6.4f',Acad.QEDGE)];
      %s2 = ['QSPEC =' sprintf('%6.4f',2.1)];
      eval(['!sed -i s/''QSPEC =1.140''/''' s2 '''/g Rchease']);
      eval(['!sed -i s/''CSSPEC=0''/''CSSPEC=1''/g Rchease']);
      eval(['!sed -i s/''NCSCAL = 4''/''NCSCAL = 1''/g Rchease']);
      s3 = ['NTOR=' sprintf('%1i',Acad.n)];
      eval(['!sed -i s/''NTOR=1''/''' s3 '''/g Rchease']);

      %CFP = linspace(1.5,0.5,11);
      CFP = 0.9;
      GR  = 3e-2;
      s1 = ['CFBAL=' sprintf('%7.5f',1)];
      g1 = ['TALPHA1= (' sprintf('%11.5e',GR)];
      w1 = 'NV     = 200';
      eval(['!sed -i s/''TALPHA1= (1.00000E-01''/''' g1 '''/g RmarsQ']);
      s3 = ['RNTOR  = -' sprintf('%1i',Acad.n)];
      eval(['!sed -i s/''RNTOR  = -1''/''' s3 '''/g RmarsQ']);
      RES_ALL = []; 

      for kc=1:length(CFP)
          CF = CFP(kc);

          s2 = ['CFBAL=' sprintf('%7.5f',CF)];
          eval(['!sed -i s/''' s1 '''/''' s2 '''/g Rchease']);
          s1 = s2;
          eval('!./Rchease')

          %%get global parameters from log_chease
          eval('!grep Q_ZERO log_chease > tmp.txt')
          eval('!grep Q_EDGE log_chease >> tmp.txt')
          eval('!grep GEXP log_chease >> tmp.txt')
          eval('!grep ''RW='' log_chease >> tmp.txt')
          fid = fopen('tmp.txt','r');
          ss  = fgetl(fid);  QZERO  = str2num(ss(5:18));
          ss  = fgetl(fid);  QEDGE  = str2num(ss(5:18));
          ss  = fgetl(fid);  BETAN  = str2num(ss(28:end));
          ss  = fgetl(fid);  NW     = floor(str2num(ss(51:end)));
          fclose(fid);
          
          g2 = ['TALPHA1= (' sprintf('%11.5e',GR)];
          eval(['!sed -i s/''' g1 '''/''' g2 '''/g RmarsQ']);
          w2 = ['NV     = ' sprintf('%3i',NW-1)];
          %eval(['!sed -i s/''' w1 '''/''' w2 '''/g RmarsQ']);

          g1 = g2;
          w1 = w2;
          eval('!./RmarsQ')

          d=load('RESULT.OUT');
          GR = d(5);
          RES_ALL = [RES_ALL; CF BETAN QZERO QEDGE d(2) d(5:6)];
          eval('!cat RESULT.OUT >> RESULT_SAVE');
          
          if GR<1e-3, break; end
      end
      %save([SDIR4 'RES_ALL_n' num2str(Acad.n)],'RES_ALL','-ascii');  
      end
      end
   end
end

if kaction==5
   %SDIR4 = [SDIR2 '179990_4000_7/'];  
   %SDIR4 = [SDIR2 '900000_4000_1/'];  
   SDIR4 = [SDIR2 '179991_3100_M/'];  

   %ss = sprintf('%3.1f',Acad.DRPW);
   %SDIR4 = [SDIR4 'DRPW_' ss '/'];
   load([SDIR4 'Acad.mat']);

   Acad.axis=[1 2.3 -1 1];

   Acad.AMP = 1/5000/1.915/0.813;
   
   if 1==0
      if 1==0
      cd(SDIR3)
      copyfile('RMZM_F.OUT',[SDIR4 'RMZM_F_EQAC'],'f');
      copyfile('PROFEQ.OUT',[SDIR4 'PROFEQ.OUT'],'f');
      copyfile('BPLASMA.OUT',[SDIR4 'BPLASMA_.OUT'],'f');
      copyfile('JPLASMA.OUT',[SDIR4 'JPLASMA_.OUT'],'f');
      copyfile('XPLASMA.OUT',[SDIR4 'XPLASMA_.OUT'],'f');
      end
      cd(SDIRM)
      MacACADplot(SDIR4,'','',''); 
      Acad
   else
      cd(SDIR3)
      copyfile('RMZM_F.OUT', 'RMZM_F_EQAC','f');
      copyfile('BPLASMA.OUT','BPLASMA_.OUT','f');
      copyfile('JPLASMA.OUT','JPLASMA_.OUT','f');
      copyfile('XPLASMA.OUT','XPLASMA_.OUT','f');
      cd(SDIRM)
      MacACADplot(SDIR3,'','','');  
   end 

   %d=load([SDIRM 'MacCOMMON/BOUNDW_D3D'])*1.6995;
   %figure(70), plot(d(:,1),d(:,2),'b-','LineWidth',3), axis([-0.5 2.5 -1.5 1.5]) 
end 

cd(SDIRM)

