% Automatic Chease&mars-* setup for D3D  (ACAD)
% starting from standard g-file and p-file for DIII-D equilibrium
global Acad

% note that the options below need to be executed in the same order as listed
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
  
Acad.n = 2;
Acad.DEVICE='D3D';  %'D3D','MAST','MAST-U','AUG'


if kaction==2|kaction==3|kaction==6|kaction==7,     SEXB='';     end
if kaction==12|kaction==13|kaction==16|kaction==17, SEXB='_EXB'; end

SDIRM  = '~/IRIS/mars-code/RZplot/';
SDIR3 = '/cscratch/liuy/WorkTEMP/';

SDIR12 = {'q95_-3.42/'};
SDIR13 = {'betan_1.80/','betan_2.00/'};
SDIR14 = {'neped_0.18/','neped_0.20/'};
SDIR16 = '~/Work/D3D_DRSEP/';
%SDIR16 = '~/Work/D3D164277/';

for k2=1:1
for k3=1:1 
for k4=1:1
    SDIR  = [SDIR12{k2} SDIR13{k3} SDIR14{k4}];
    %SDIR  = '';
    %SDIR21 = [SDIR16 'Expt/' SDIR];
    SDIR21 = [SDIR16 'Data/' SDIR];
    cd(SDIR21)
    F  = split(ls);
    NF = length(F)-1;

    cd(SDIRM);
    
    if kaction==19
       AcadA.DRSEP=zeros(NF,1);
       AcadA.XedgeX=[];
       AcadA.XedgeY=[];
       AcadA.BpestX=[];
       AcadA.BpestY=[];
       AcadA.QprofX=[];
       AcadA.QprofY=[];
    end

    if kaction < 20
    for k6=1:NF
        FF = F{k6};
        Acad.DRSEP = str2num(FF(7:end));
        SDIR1 = [SDIR21 '/'];
        SDIR2 = [SDIR16 'Data/' SDIR FF '/'];

        cd(SDIR2)
        G=split(ls);
        GG=G{1}; GG=GG(1:end);
        SDIR4=[SDIR2 GG '/'];
        SDIRU = [SDIR3 FF '/IU/'];
        SDIRL = [SDIR3 FF '/IL/'];
        
        if kaction==1
           cd(SDIRM)
           MacACADprepare(SDIR1,SDIR2,SDIR3,'RMP'); 
        end
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
           copyfile([SDIR4 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
           copyfile([SDIR4 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
           copyfile([SDIR4 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
           copyfile([SDIR4 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
           copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
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
           copyfile([SDIR4 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
           copyfile([SDIR4 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
           copyfile([SDIR4 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
           copyfile([SDIR4 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
           copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
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
           copyfile([SDIR4 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
           copyfile([SDIR4 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
           copyfile([SDIR4 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
           copyfile([SDIR4 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
           copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
           copyfile([SDIR4 'Rchease_SAVE'],'Rchease','f');
           copyfile([SDIR4 'RmarsQ_IU_SAVE'],'RmarsQ','f');
           eval(['!sed -i s/''INCFEED= 8''/''INCFEED= 4''/g RmarsQ']);
           copyfile([SDIRM 'MacCOMMON/runit.bat'],'runit.bat','f');
           copyfile([SDIRM 'MacCOMMON/runit.sbatch'],'runit.sbatch','f');
           eval('!sbatch runit.sbatch')

           cd(SDIRL)
           copyfile([SDIR4 'EXPEQ_SAVE'],'EXPEQ','f');
           copyfile([SDIR4 'PROFDEN_SAVE.IN'],'PROFDEN.IN','f');
           copyfile([SDIR4 'PROFROT_SAVE.IN'],'PROFROT.IN','f');
           copyfile([SDIR4 'PROFTI_SAVE.IN'],'PROFTI.IN','f');
           copyfile([SDIR4 'PROFTE_SAVE.IN'],'PROFTE.IN','f');
           copyfile([SDIR4 'PROFWE_SAVE.IN'],'PROFWE.IN','f');
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
           ang0=linspace(0,360,17); ang0=ang0(1:end-1);
           phi0=ang0*pi/180; 
           ep0 = exp(i*phi0);
           for kkp=1:length(ep0)
           SANG = ['_' num2str(ang0(kkp))];
           BIU=load([SDIR4 'BPLASMA' SEXB '_IUPLS.OUT']);
           BIL=load([SDIR4 'BPLASMA' SEXB '_ILPLS.OUT']);
           MSMAX=floor(BIU(1,1));
           TU = BIU; 
           TL = BIL;
           TU1 = TU(MSMAX+2:end,1)+i*TU(MSMAX+2:end,2);
           TU2 = TU(MSMAX+2:end,3)+i*TU(MSMAX+2:end,4);
           TU3 = TU(MSMAX+2:end,5)+i*TU(MSMAX+2:end,6);
           TL1 = TL(MSMAX+2:end,1)+i*TL(MSMAX+2:end,2);
           TL2 = TL(MSMAX+2:end,3)+i*TL(MSMAX+2:end,4);
           TL3 = TL(MSMAX+2:end,5)+i*TL(MSMAX+2:end,6);
           TE1 = TU1 + ep0(kkp)*TL1;
           TE2 = TU2 + ep0(kkp)*TL2;
           TE3 = TU3 + ep0(kkp)*TL3;
           TE  = [real(TE1) imag(TE1) real(TE2) imag(TE2) real(TE3) imag(TE3)];
           BIE=[BIU(1:MSMAX+1,:); TE];
           save([SDIR4 'BPLASMA' SEXB SANG '_IEPLS.OUT'],'BIE','-ascii');

           BIU=load([SDIR4 'BPLASMA' SEXB '_IU.OUT']);
           BIL=load([SDIR4 'BPLASMA' SEXB '_IL.OUT']);
           MSMAX=floor(BIU(1,1));
           TU = BIU; 
           TL = BIL;
           TU1 = TU(MSMAX+2:end,1)+i*TU(MSMAX+2:end,2);
           TU2 = TU(MSMAX+2:end,3)+i*TU(MSMAX+2:end,4);
           TU3 = TU(MSMAX+2:end,5)+i*TU(MSMAX+2:end,6);
           TL1 = TL(MSMAX+2:end,1)+i*TL(MSMAX+2:end,2);
           TL2 = TL(MSMAX+2:end,3)+i*TL(MSMAX+2:end,4);
           TL3 = TL(MSMAX+2:end,5)+i*TL(MSMAX+2:end,6);
           TE1 = TU1 + ep0(kkp)*TL1;
           TE2 = TU2 + ep0(kkp)*TL2;
           TE3 = TU3 + ep0(kkp)*TL3;
           TE  = [real(TE1) imag(TE1) real(TE2) imag(TE2) real(TE3) imag(TE3)];
           BIE=[BIU(1:MSMAX+1,:); TE];
           save([SDIR4 'BPLASMA' SEXB SANG '_IE.OUT'],'BIE','-ascii');

           XIU=load([SDIR4 'XPLASMA' SEXB '_IU.OUT']);
           XIL=load([SDIR4 'XPLASMA' SEXB '_IL.OUT']);
           MSMAX=floor(XIU(1,1)+XIU(1,2)); 
           TU = XIU; 
           TL = XIL;
           TU1 = TU(MSMAX+2:end,1)+i*TU(MSMAX+2:end,2);
           TU2 = TU(MSMAX+2:end,3)+i*TU(MSMAX+2:end,4);
           TU3 = TU(MSMAX+2:end,5)+i*TU(MSMAX+2:end,6);
           TL1 = TL(MSMAX+2:end,1)+i*TL(MSMAX+2:end,2);
           TL2 = TL(MSMAX+2:end,3)+i*TL(MSMAX+2:end,4);
           TL3 = TL(MSMAX+2:end,5)+i*TL(MSMAX+2:end,6);
           TE1 = TU1 + ep0(kkp)*TL1;
           TE2 = TU2 + ep0(kkp)*TL2;
           TE3 = TU3 + ep0(kkp)*TL3;
           TE  = [real(TE1) imag(TE1) real(TE2) imag(TE2) real(TE3) imag(TE3)];
           XIE=[XIU(1:MSMAX+1,:); TE];
           save([SDIR4 'XPLASMA' SEXB SANG '_IE.OUT'],'XIE','-ascii');

           JIU=load([SDIR4 'JPLASMA' SEXB '_IU.OUT']);
           JIL=load([SDIR4 'JPLASMA' SEXB '_IL.OUT']);
           MSMAX=floor(JIU(1,1)); 
           TU = JIU; 
           TL = JIL;
           TU1 = TU(MSMAX+2:end,1)+i*TU(MSMAX+2:end,2);
           TU2 = TU(MSMAX+2:end,3)+i*TU(MSMAX+2:end,4);
           TU3 = TU(MSMAX+2:end,5)+i*TU(MSMAX+2:end,6);
           TL1 = TL(MSMAX+2:end,1)+i*TL(MSMAX+2:end,2);
           TL2 = TL(MSMAX+2:end,3)+i*TL(MSMAX+2:end,4);
           TL3 = TL(MSMAX+2:end,5)+i*TL(MSMAX+2:end,6);
           TE1 = TU1 + ep0(kkp)*TL1;
           TE2 = TU2 + ep0(kkp)*TL2;
           TE3 = TU3 + ep0(kkp)*TL3;
           TE  = [real(TE1) imag(TE1) real(TE2) imag(TE2) real(TE3) imag(TE3)];
           JIE=[JIU(1:MSMAX+1,:); TE];
           save([SDIR4 'JPLASMA' SEXB SANG '_IE.OUT'],'JIE','-ascii');

           PIU=load([SDIR4 'PPLASMA' SEXB '_IU.OUT']);
           PIL=load([SDIR4 'PPLASMA' SEXB '_IL.OUT']);
           MSMAX=floor(PIU(1,1)); 
           PIE=[PIU(1:MSMAX+1,:); PIU(MSMAX+2:end,:)+PIL(MSMAX+2:end,:)];
           save([SDIR4 'PPLASMA' SEXB '_IE.OUT'],'PIE','-ascii');
           end
        end

        if kaction==19
           load([SDIR4 'Acad.mat']);
           Acad.axis = [1 2.4 -1.2 1.2];
           ang0=linspace(0,360,17); ang0=ang0(1:end-1);
           for kkp=1:length(ang0)
           SANG = [num2str(ang0(kkp))];
           SCOIL = [SANG '_IE'];  %'IU','IL','IE' 
           SEXB  = '';    %''(fluid),'_EXB'
           SFLD  = '';    %''(TOT),'PLS','VAC'

           cd(SDIRM)
           MacACADplot(SDIR4,SCOIL,'PLS',SEXB); 
           MacACADpest(SDIR4,SCOIL,'',SEXB);  

           AcadA.DRSEP(k6) = Acad.DRSEP;

           AcadA.XedgeX = [AcadA.XedgeX Acad.Xedge(:,1)];
           AcadA.XedgeY = [AcadA.XedgeY Acad.Xedge(:,2)];

           AcadA.BpestX = [AcadA.BpestX; 0; Acad.Bpest(:,1)];
           AcadA.BpestY = [AcadA.BpestY; 0; Acad.Bpest(:,2)];

           AcadA.QprofX = [AcadA.QprofX Acad.Qprof(:,1)];
           AcadA.QprofY = [AcadA.QprofY Acad.Qprof(:,2)];
           end
        end
    end
    end

    if kaction==19
       SDIR8 = [SDIR16 'Data/' SDIR 'Resu/'];
mkdir(SDIR8);
       save([SDIR8 'AcadA' SEXB '_' SCOIL '.mat'],'AcadA'); 
    end

    if kaction==20
       SDIR8 = [SDIR16 'Data/' SDIR 'Resu/'];

       %%plot results with fluid rotation and even parity n=3 I-coils
       load([SDIR8 'AcadA_IE.mat']);
       ang0 = linspace(0,360,17); ang0=ang0(1:end-1);
       AcadA.DRSEP = ang0;
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
       XDISP  = DRSEP;
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
           xlabel('poloidal angle [deg.]','FontSize',18,'FontWeight','Bold'),
           ylabel('{\xi}_n [mm]','FontSize',18,'FontWeight','Bold'),
           ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
           XRATIO(k20) = max(Xn(:,k20))/Xn(1,k20);
           [C,KK]=min(abs(chi(:,k20)-90));
           XDISP(k20) = Xn(KK,k20);

           hf=figure(3);
           KK=JJ(k20)+1:JJ(k20)+NN(k20);
           plot(Bpsi(KK),Bpest(KK),'-','LineWidth',3,'Color',SCOL), hold on 
           xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
           ylabel('|B^1_{res}|{\times}10^4','FontSize',18,'FontWeight','Bold'),
           ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
           [X,L]=min(abs(Bpsi(KK)-0.995)); 
           B1EDGE(k20) = Bpest(KK(L));
       end
       figure(3), a=axis; axis([0.95 1 a(3) a(4)]);

       figure(2), legend(num2str(DRSEP))
       figure(3), legend(num2str(DRSEP))

       hf=figure(4);
       plot(DRSEP,XRATIO,'b-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b'), hold on 
       xlabel('\Delta\Phi [deg.]','FontSize',18,'FontWeight','Bold'),
       ylabel('{\xi}_{max}/{\xi}_{LFS}','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

       hf=figure(5);
       plot(DRSEP,XDISP,'b-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b'), hold on 
       xlabel('\Delta\Phi [deg.]','FontSize',18,'FontWeight','Bold'),
       ylabel('{\xi}_X [mm]','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')

       hf=figure(6);
       plot(DRSEP,B1EDGE,'b-o','LineWidth',3,'MarkerSize',9,'MarkerFaceColor','b'), hold on 
       xlabel('\Delta\Phi [deg.]','FontSize',18,'FontWeight','Bold'),
       ylabel('|B^1_{0.995}|{\times}10^4','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
    end

    if kaction==21
       %SDIR4 = [SDIR16 'Data/' SDIR 'drsep_-0.03/146065.02510.6/'];
       %SDIR4 = [SDIR16 'Data/' SDIR 'drsep_-4.60/146065.02510.1/'];
       SDIR4 = [SDIR16 'Data/' SDIR 'drsep_4.60/146065.02510.11/'];

       load([SDIR4 'Acad.mat']);
       Acad.axis = [1 2.4 -1.2 1.2];
       SCOIL = 'IE';  %'IU','IL','IE' 
       SEXB  = '';    %''(fluid),'_EXB'
       SFLD  = '';    %''(TOT),'PLS','VAC'

       cd(SDIRM)
       MacACADplot(SDIR4,SCOIL,'PLS',SEXB);  
       MacACADpest(SDIR4,SCOIL,'',SEXB);  

       if 1==0
       hf=figure;
       df = load([SDIR4 'PROFROT_SAVE.IN']);
       de = load([SDIR4 'PROFWE_SAVE.IN']);
       plot(df(2:end,1).^2,df(2:end,2)/1e+3,'r-','LineWidth',3), hold on 
       plot(de(2:end,1).^2,de(2:end,2)/1e+3,'b-','LineWidth',3), hold on 
       xlabel('{\psi}_p','FontSize',18,'FontWeight','Bold'),
       ylabel('rotation frequency [krad/s]','FontSize',18,'FontWeight','Bold'),
       ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
       legend('FLUID','ExB')
       end
    end
end
end
end

cd(SDIRM)

