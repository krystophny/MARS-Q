function MacGetXIRZ(R,Z,filename)

global Mac SDIR

kavrg = 2; %0: no average, direct plot of (XI_R,XI_Z)
           %1: average over poloidal angle
           %2: average over both poloidal angle and minor radius

% read in XPLASMA_XIRZ-data
d  = load(filename);
Ns = Mac.Ns1-1;
Nc = round(size(d,1)/Ns);

s  = d(:,1);
chi= d(:,2);  II=find(chi>pi); chi(II)=chi(II)-2*pi;
XIR= d(:,3)+d(:,4)*i;
XIZ= d(:,5)+d(:,6)*i;

s   = transpose(reshape(s,Nc,Ns)); 
chi = transpose(reshape(chi,Nc,Ns)); 
XIR = transpose(reshape(XIR,Nc,Ns)); 
XIZ = transpose(reshape(XIZ,Nc,Ns)); 
s   = [s s(:,1)];
chi = [chi chi(:,1)];
XIR = [XIR XIR(:,1)];
XIZ = [XIZ XIZ(:,1)];

if kavrg==0
   XIR(1:10,:) = 0;
   XIZ(1:10,:) = 0;
elseif kavrg==1
   onev= ones(1,size(XIR,2));
   for k=1:size(XIR,1)
       XIR(k,:)=mean(XIR(k,:))*onev;
       XIZ(k,:)=mean(XIZ(k,:))*onev;
   end
elseif kavrg==2
   onev = ones(size(XIR)); 
   XIR = mean(mean(XIR))*onev;
   XIZ = mean(mean(XIZ))*onev;
end

[ss,cc] = meshgrid(Mac.s(1:Mac.Ns1),Mac.chi);
ss = transpose(ss);
cc = transpose(cc);
RR = griddata(ss,cc,R(1:Mac.Ns1,:),s,chi);
ZZ = griddata(ss,cc,Z(1:Mac.Ns1,:),s,chi);

if Mac.plot_XIRZ
  hf=figure(10*Mac.plot_XIRZ+1);
  pcolor(RR*Mac.R0EXP,ZZ*Mac.R0EXP,real(XIR)), hold on,
  shading interp
  axis equal
  colorbar
  xlabel('R [m]','FontSize',18,'FontWeight','Bold')
  ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
  title('\xi_R','FontSize',18,'FontWeight','Bold')
  ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
  axis([1 2.5 -1.25 1.25])


  %figure(10*Mac.plot_XIRZ+2)
  %pcolor(RR*Mac.R0EXP,ZZ*Mac.R0EXP,imag(XIR)), hold on,
  %shading interp
  %axis equal
  %colorbar

  hf=figure(10*Mac.plot_XIRZ+3);
  pcolor(RR*Mac.R0EXP,ZZ*Mac.R0EXP,real(XIZ)), hold on,
  shading interp
  axis equal
  colorbar
  xlabel('R [m]','FontSize',18,'FontWeight','Bold')
  ylabel('Z [m]','FontSize',18,'FontWeight','Bold')
  title('\xi_Z','FontSize',18,'FontWeight','Bold')
  ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
  axis([1 2.5 -1.25 1.25])

  %figure(10*Mac.plot_XIRZ+4)
  %pcolor(RR*Mac.R0EXP,ZZ*Mac.R0EXP,imag(XIZ)), hold on,
  %shading interp
  %axis equal
  %colorbar
end


