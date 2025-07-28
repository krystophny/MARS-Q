% Automatic Chease&mars-* setup for D3D  (ACAD)
% starting from standard g-file and p-file for DIII-D equilibrium
global Acad

% note that the options below need to be executed in the same as listed
kaction = 20; %1: prepare MARS-F runs 
              %2: run MARS-F for upper and lower row RMP response, fluid rotation
	      %3: save MARS-F response data, fluid rotation
              %12: run MARS-F for upper and lower row RMP response, ExB flow
              %13: save MARS-F response data, ExB flow
              %4: run MARS-F for vacuum field
              %5: save MARS-F vacuum field data
              %6: get pure plasma response dB data, fluid flow
              %16: get pure plasma response dB data, ExB flow
              %7: combine response with given coil phasing, fluid flow
              %17: combine response with given coil phasing, ExB flow
              %11: generate RMZM_F_PEST
	      %19: save and plot MARS-F data for surface displacement and dBp along wall, 
              %    for all cases. Note that need to repeat with SEXB='' and '_EXB'
              %20: plot results from kaction=19
              %21: plot all MARS-F run data for selected cases
  
Acad.n = 3;
Acad.DEVICE='MAST';  %'D3D','MAST','MAST-U','AUG'


if kaction==2|kaction==3|kaction==6|kaction==7,     SEXB='';     end
if kaction==12|kaction==13|kaction==16|kaction==17, SEXB='_EXB'; end

SDIRM  = '~/IRIS/mars-code/RZplot/';
SDIR0 = '~/Work/MAST_DRSEP/';
SDIR1 = [SDIR0 'Expt/'];
SDIR2 = [SDIR0 'DataN/'];
SDIR3 = '/cscratch/liuy/WorkTEMP/';

cd(SDIRM)

if kaction==19
   AcadA.DRSEP=[];
   AcadA.XedgeX=[];
   AcadA.XedgeY=[];
   AcadA.BpestX=[];
   AcadA.BpestY=[];
   AcadA.QprofX=[];
   AcadA.QprofY=[];
end

if kaction==1
   cd(SDIRM)
   MacACADprepareMAST(SDIR1,SDIR2,SDIR3,'RMP'); 
end

if kaction < 20 & kaction > 1
   cd(SDIR2)
   F  = split(ls);
   NF = length(F)-1;

   for kf=1:NF
   FF = F{kf};
   SDIR4=[SDIR2 FF '/'];
   SDIRU = [SDIR3 FF '/IU/'];
   SDIRL = [SDIR3 FF '/IL/'];
   
   if kaction==11
      cd(SDIRU)
      copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
      eval('!sed -i s/''NEGP=-1''/''NEGP= 0''/g Rchease');
      eval('!sed -i s/''NER=1''/''NER=2''/g Rchease');
      eval('!./Rchease');

      copyfile([SDIR4 'RmarsQ_IU_SAVE'],'RmarsQ','f');
      eval(['!sed -i s/''M1     =-33''/''M1     =-1''/g RmarsQ']);
      eval(['!sed -i s/''M2     = 33''/''M2     = 1''/g RmarsQ']);
      eval(['!sed -i s/''INCFEED= 8''/''INCFEED= 4''/g RmarsQ']);
      eval('!./RmarsQ')

      copyfile('RMZM_F.OUT',[SDIR4 'RMZM_F_PEST'],'f');
   end

   if kaction==2 | kaction==12
      mkdir(SDIRU)
      mkdir(SDIRL)
      
      cd(SDIRU)
      copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
      copyfile([SDIR1 'PROFDEN_MAST24460'],'PROFDEN.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFROT.IN','f');
      copyfile([SDIR1 'PROFTI_MAST24460'],'PROFTI.IN','f');
      copyfile([SDIR1 'PROFTE_MAST24460'],'PROFTE.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFWE.IN','f');
      copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
      copyfile([SDIR4 'RmarsQ_IU_SAVE'],'RmarsQ','f');
      copyfile([SDIRM 'MacCOMMON/runit.bat'],'runit.bat','f');
      copyfile([SDIRM 'MacCOMMON/runit.sbatch'],'runit.sbatch','f');

      if kaction==12
         load([SDIR4 'Acad.mat']);
         copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFROT.IN','f');
         s1 = ['ROTE   =' sprintf('%11.4e',Acad.ROT0)];
         s2 = ['ROTE   =' sprintf('%11.4e',Acad.WE0)];
         eval(['!sed -i s/''' s1 '''/''' s2 '''/g RmarsQ']);
      end
      eval('!sbatch runit.sbatch')

      cd(SDIRL)
      copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
      copyfile([SDIR1 'PROFDEN_MAST24460'],'PROFDEN.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFROT.IN','f');
      copyfile([SDIR1 'PROFTI_MAST24460'],'PROFTI.IN','f');
      copyfile([SDIR1 'PROFTE_MAST24460'],'PROFTE.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFWE.IN','f');
      copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
      copyfile([SDIR4 'RmarsQ_IL_SAVE'],'RmarsQ','f');
      copyfile([SDIRM 'MacCOMMON/runit.bat'],'runit.bat','f');
      copyfile([SDIRM 'MacCOMMON/runit.sbatch'],'runit.sbatch','f');
      if kaction==12
         copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFROT.IN','f');
         eval(['!sed -i s/''' s1 '''/''' s2 '''/g RmarsQ']); 
      end
      eval('!sbatch runit.sbatch')
   end

   if kaction==3 | kaction==13

      copyfile([SDIRU 'RMZM_F.OUT'],[SDIR4 'RMZM_F_EQAC'],'f');
      copyfile([SDIRU 'PROFEQ.OUT'],[SDIR4 'PROFEQ.OUT'],'f');
      copyfile([SDIRU 'XPLASMA.OUT'],[SDIR4 'XPLASMA' SEXB '_IU.OUT'],'f');
      copyfile([SDIRU 'BPLASMA.OUT'],[SDIR4 'BPLASMA' SEXB '_IU.OUT'],'f');
      copyfile([SDIRU 'JPLASMA.OUT'],[SDIR4 'JPLASMA' SEXB '_IU.OUT'],'f');
      copyfile([SDIRU 'PPLASMA.OUT'],[SDIR4 'PPLASMA' SEXB '_IU.OUT'],'f');
      copyfile([SDIRL 'XPLASMA.OUT'],[SDIR4 'XPLASMA' SEXB '_IL.OUT'],'f');
      copyfile([SDIRL 'BPLASMA.OUT'],[SDIR4 'BPLASMA' SEXB '_IL.OUT'],'f');
      copyfile([SDIRL 'JPLASMA.OUT'],[SDIR4 'JPLASMA' SEXB '_IL.OUT'],'f');
      copyfile([SDIRL 'PPLASMA.OUT'],[SDIR4 'PPLASMA' SEXB '_IL.OUT'],'f');
   end

   if kaction==4
      cd(SDIRU)
      copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
      copyfile([SDIR1 'PROFDEN_MAST24460'],'PROFDEN.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFROT.IN','f');
      copyfile([SDIR1 'PROFTI_MAST24460'],'PROFTI.IN','f');
      copyfile([SDIR1 'PROFTE_MAST24460'],'PROFTE.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFWE.IN','f');
      copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
      copyfile([SDIR4 'RmarsQ_IU_SAVE'],'RmarsQ','f');
      eval(['!sed -i s/''INCFEED= 8''/''INCFEED= 4''/g RmarsQ']);
      copyfile([SDIRM 'MacCOMMON/runit.bat'],'runit.bat','f');
      copyfile([SDIRM 'MacCOMMON/runit.sbatch'],'runit.sbatch','f');
      eval('!sbatch runit.sbatch')

      cd(SDIRL)
      copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
      copyfile([SDIR1 'PROFDEN_MAST24460'],'PROFDEN.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFROT.IN','f');
      copyfile([SDIR1 'PROFTI_MAST24460'],'PROFTI.IN','f');
      copyfile([SDIR1 'PROFTE_MAST24460'],'PROFTE.IN','f');
      copyfile([SDIR1 'PROFROT_MAST24460'],'PROFWE.IN','f');
      copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
      copyfile([SDIR4 'RmarsQ_IL_SAVE'],'RmarsQ','f');
      eval(['!sed -i s/''INCFEED= 8''/''INCFEED= 4''/g RmarsQ']);
      copyfile([SDIRM 'MacCOMMON/runit.bat'],'runit.bat','f');
      copyfile([SDIRM 'MacCOMMON/runit.sbatch'],'runit.sbatch','f');
      eval('!sbatch runit.sbatch')
   end
   if kaction==5
      copyfile([SDIRU 'BPLASMA.OUT'],[SDIR4 'BPLASMA_IUVAC.OUT'],'f');
      copyfile([SDIRL 'BPLASMA.OUT'],[SDIR4 'BPLASMA_ILVAC.OUT'],'f');
   end

   if kaction==6 | kaction==16
      BTOT=load([SDIR4 'BPLASMA' SEXB '_IU.OUT']);
      BVAC=load([SDIR4 'BPLASMA_IUVAC.OUT']);
      MSMAX=floor(BTOT(1,1));
      BPLS=[BTOT(1:MSMAX+1,:); BTOT(MSMAX+2:end,:)-BVAC(MSMAX+2:end,:)];
      save([SDIR4 'BPLASMA' SEXB '_IUPLS.OUT'],'BPLS','-ascii');

      BTOT=load([SDIR4 'BPLASMA' SEXB '_IL.OUT']);
      BVAC=load([SDIR4 'BPLASMA_ILVAC.OUT']);
      MSMAX=floor(BTOT(1,1));
      BPLS=[BTOT(1:MSMAX+1,:); BTOT(MSMAX+2:end,:)-BVAC(MSMAX+2:end,:)];
      save([SDIR4 'BPLASMA' SEXB '_ILPLS.OUT'],'BPLS','-ascii');
   end

   if kaction==7 | kaction==17
      BIU=load([SDIR4 'BPLASMA' SEXB '_IUPLS.OUT']);
      BIL=load([SDIR4 'BPLASMA' SEXB '_ILPLS.OUT']);
      MSMAX=floor(BIU(1,1));
      BIE=[BIU(1:MSMAX+1,:); BIU(MSMAX+2:end,:)+BIL(MSMAX+2:end,:)];
      save([SDIR4 'BPLASMA' SEXB '_IEPLS.OUT'],'BIE','-ascii');

      BIU=load([SDIR4 'BPLASMA' SEXB '_IU.OUT']);
      BIL=load([SDIR4 'BPLASMA' SEXB '_IL.OUT']);
      MSMAX=floor(BIU(1,1));
      BIE=[BIU(1:MSMAX+1,:); BIU(MSMAX+2:end,:)+BIL(MSMAX+2:end,:)];
      save([SDIR4 'BPLASMA' SEXB '_IE.OUT'],'BIE','-ascii');

      XIU=load([SDIR4 'XPLASMA' SEXB '_IU.OUT']);
      XIL=load([SDIR4 'XPLASMA' SEXB '_IL.OUT']);
      MSMAX=floor(XIU(1,1)+XIU(1,2)); 
      XIE=[XIU(1:MSMAX+1,:); XIU(MSMAX+2:end,:)+XIL(MSMAX+2:end,:)];
      save([SDIR4 'XPLASMA' SEXB '_IE.OUT'],'XIE','-ascii');

      JIU=load([SDIR4 'JPLASMA' SEXB '_IU.OUT']);
      JIL=load([SDIR4 'JPLASMA' SEXB '_IL.OUT']);
      MSMAX=floor(JIU(1,1)); 
      JIE=[JIU(1:MSMAX+1,:); JIU(MSMAX+2:end,:)+JIL(MSMAX+2:end,:)];
      save([SDIR4 'JPLASMA' SEXB '_IE.OUT'],'JIE','-ascii');

      PIU=load([SDIR4 'PPLASMA' SEXB '_IU.OUT']);
      PIL=load([SDIR4 'PPLASMA' SEXB '_IL.OUT']);
      MSMAX=floor(PIU(1,1)); 
      PIE=[PIU(1:MSMAX+1,:); PIU(MSMAX+2:end,:)+PIL(MSMAX+2:end,:)];
      save([SDIR4 'PPLASMA' SEXB '_IE.OUT'],'PIE','-ascii');
   end

   if kaction==19
      load([SDIR4 'Acad.mat']);
      %Acad.axis = [1 2.4 -1.2 1.2];
      SCOIL = 'IE';  %'IU','IL','IE' 
      SEXB  = '';    %''(fluid),'_EXB'
      SFLD  = '';    %''(TOT),'PLS','VAC'

      cd(SDIRM)
      MacACADplot(SDIR4,SCOIL,'PLS',SEXB); 
      MacACADpest(SDIR4,SCOIL,'',SEXB);  

      AcadA.DRSEP = [AcadA.DRSEP Acad.DRSEP];

      AcadA.XedgeX = [AcadA.XedgeX Acad.Xedge(:,1)];
      AcadA.XedgeY = [AcadA.XedgeY Acad.Xedge(:,2)];

      AcadA.BpestX = [AcadA.BpestX; 0; Acad.Bpest(:,1)];
      AcadA.BpestY = [AcadA.BpestY; 0; Acad.Bpest(:,2)];

      AcadA.QprofX = [AcadA.QprofX Acad.Qprof(:,1)];
      AcadA.QprofY = [AcadA.QprofY Acad.Qprof(:,2)];
   end
   end
end

if kaction==19
   SDIR8 = [SDIR2 'Resu/'];
   mkdir(SDIR8);
   save([SDIR8 'AcadA' SEXB '_' SCOIL '.mat'],'AcadA'); 
end

if kaction==20
   SDIR8 = [SDIR2 'Resu/'];

   %%plot results with fluid rotation and even parity n=3 I-coils
   load([SDIR8 'AcadA_IE.mat']);
   [DRSEP,II] = sort(AcadA.DRSEP);

   psip = AcadA.QprofX(:,II);
   qprof= AcadA.QprofY(:,II);

   chi = AcadA.XedgeX(:,II);
   Xn  = AcadA.XedgeY(:,II);

   JJ  = find(AcadA.BpestX==0); 
   NN  = [diff(JJ)-1; length(AcadA.BpestX)-JJ(end)]; 
   JJ  = JJ(II); NN = NN(II);
   Bpsi  = AcadA.BpestX;
   Bpest = AcadA.BpestY;
   
   ND = length(DRSEP);
   DMIN = min(DRSEP); 
   DMAX = max(DRSEP);
   XRATIO = DRSEP;
   B1EDGE = DRSEP;
   for k20=1:ND
       DR = DRSEP(k20);
       if DR>0, SCOL=[DR/DMAX 1-DR/DMAX 0]; end
       if DR<0, SCOL=[0 1-DR/DMIN DR/DMIN]; end
       if DR==0, SCOL=[0 1 0]; end

       hf=figure(1);
       plot(psip(:,k20),qprof(:,k20),'-','LineWidth',3,'Color',SCOL), hold on 
       xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
       ylabel('safety factor','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

       hf=figure(2);
       plot(chi(:,k20),Xn(:,k20),'-','LineWidth',3,'Color',SCOL), hold on 
       xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
       ylabel('{\xi}_n [mm]','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
       XRATIO(k20) = max(Xn(:,k20))/Xn(1,k20);

       hf=figure(3);
       KK=JJ(k20)+1:JJ(k20)+NN(k20);
       plot(Bpsi(KK),Bpest(KK),'-','LineWidth',3,'Color',SCOL), hold on 
       xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
       ylabel('|B^1_{res}|{\times}10^4','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
       [X,L]=min(abs(Bpsi(KK)-0.995)); 
       B1EDGE(k20) = Bpest(KK(L));
   end
   figure(3), a=axis; axis([0.9 1 a(3) a(4)]);

   figure(1), legend(num2str(DRSEP(:)))
   figure(2), legend(num2str(DRSEP(:)))
   figure(3), legend(num2str(DRSEP(:)))

   hf=figure(4);
   plot(DRSEP,XRATIO,'r--s','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','r'), hold on 
   xlabel('{\alpha}_S','FontSize',18,'FontWeight','Bold'),
   ylabel('{\xi}_{max}/{\xi}_{LFS}','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

   hf=figure(5);
   plot(DRSEP,B1EDGE,'b-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b'), hold on 
   xlabel('{\alpha}_S','FontSize',18,'FontWeight','Bold'),
   ylabel('|B^1_{0.995}|{\times}10^4','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
end

if kaction==21
   SDIR4 = [SDIR2 'cfac_-0.60/'];

   load([SDIR4 'Acad.mat']);
   Acad.axis = [0.2 1.5 -1.1 1.1];
   SCOIL = 'IE';  %'IU','IL','IE' 
   SEXB  = '';    %''(fluid),'_EXB'
   SFLD  = '';    %''(TOT),'PLS','VAC'

   cd(SDIRM)
   MacACADplot(SDIR4,SCOIL,'PLS',SEXB);  
   MacACADpest(SDIR4,SCOIL,'',SEXB);  

   if 1==0
   hf=figure;
   df = load([SDIR1 'PROFROT_MAST24460']);
   plot(df(2:end,1).^2,df(2:end,2)*51.939,'r-','LineWidth',3), hold on 
   xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
   ylabel('rotation frequency [krad/s]','FontSize',18,'FontWeight','Bold'),
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
   end
end

cd(SDIRM)

