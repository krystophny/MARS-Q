%=====================================================================
%convert the perturbation to Boozzer coordinats, as the Input for ORBIC
%Here, only transfer \xi^\psi
%steps 

%1. get J_p/J_B; J_p: Jacobian in PEST,  J_B  Jacobian in Boozer
%2. get relatetion between thetaB thetaP; 
    %thetaB: Boozer poloidal angle,
    %thetaP: PEST poloidal angle,
%3. omega: periodic function determins the transformation 
           % from (thetaP, phiP) to (thetaB,phiB)
           % it is only function of $$\thetaB
%4. get \xi^\psi_mn in PEST
   %and caculate the transformation using the
   %Fourier method


   
%5 mmin_orbit, mmax_orbit, the minimum/maximum harmonic number used in ORBIT

%-----------------------------------
%Modify History, 
% 10/3 add r_grid, mode_orbit by G.Z.Hao
%created by G.Z.Hao 0819_2015

%======================================================

   
   global Mac
   global SDIR 
   format long 
   SDIR = '/cscratch/liuy/WorkReKink/';
   Mac.Nchi   = 1025;    % number of points along poloidal angle 'chi'
   Mac.Norm   = 2.4229e-3;  % normalization factor
   Mac.cut_RZ  = -1.0; %cut RZ-space beyond the value, not effective when <0 
   Mac.plot_RZM  = 0;
   Mac.plot_RZ  = 0;
   Mac.Nm2    = 80;     % number of poloidal harmonics for RZ-construction
   Mac.Nm0    = [];     % number of poloidal harmonics for equilibrium
   Mac.Nm1    = [];     % number of poloidal harmonics for stability
   Mac.rw     = 1.0; %rviw;
   Mac.rs     = Mac.rw; % the position of equivalent sheet current from CARIDDI
   Mac.plot_JACOB0=0;
   Mac.CHIS    = 0.0;
   
   
   Mac.plot_shape =0; 
   Mac.plot_thickness=0;
   Mac.plot_EQ1=0.0;
   Mac.plot_EQ2=0.0;
   Mac.VNORM=1.0;
   Mac.plot_VM = 0;
   %----------------------------------------------------------------------------
   Mac.chi = linspace(-pi,pi,Mac.Nchi); %The region of \chi is in[-pi, pi] in output of Chease
   %Mac.theta_psi_B = linspace(-pi,pi,Mac.Nchi); %Boozer poloidal angle
   %theta_B_psi=Mac.theta_psi_B;
    
   theta_B_psi=zeros(1,Mac.Nchi); 
   d_theta_B_psi=zeros(1,Mac.Nchi);
   omega_psi=zeros(1,Mac.Nchi);
   omega_psi1=zeros(1,Mac.Nchi);
   
   [RM,ZM] = MacReadRMZM([SDIR 'RMZM_F.OUT']); % RM ZM, the harmonics of the (R,Z) coordinates, from RMZM_F.OUT after running Chease
   R0=Mac.R0EXP;
   B0=Mac.B0EXP;

   [R,Z] = MacGetRZ(RM,ZM);
   
   d_theta_B=zeros(Mac.Ns1,Mac.Nchi);
   theta_B=zeros(Mac.Ns1,Mac.Nchi);%Mac.Ns1 from RMZM_F
   omega=zeros(Mac.Ns1,Mac.Nchi);%Mac.Ns1 from RMZM_F
   omega1=zeros(Mac.Ns1,Mac.Nchi);%Mac.Ns1 from RMZM_F
   C_psi=zeros(1,Mac.Ns1);% the normalization factor at each flux surface. 
   itheta_B_0_psi=zeros(1,Mac.Ns1);%record the element # for theta_B=0 for each flux
   omega_psi_0=zeros(1,Mac.Ns1);% omega integration constant
   
   
   R2=R(1:Mac.Ns1,:).^2;
   R2=R2*R0^2;
   %calculate the magnetic field B
   [dRds,dZds,dRdchi,dZdchi,jacobian] = MacGetUnitVec(R,Z);
   jacobian2=jacobian(1:Mac.Ns1,:).^2;%s1chi_1 s1chi_2 s1chi_3 ....]
                                      %s2chi_1 s2chi_2 s1chi_3 ....
   dRdchi2=dRdchi(1:Mac.Ns1,:).^2;
   dZdchi2=dZdchi(1:Mac.Ns1,:).^2;
   
   PROF=load([SDIR 'PROFEQ.OUT']);
   DPSIDS=PROF(:,12);
   DPSIDS_full=repmat(DPSIDS,1,Mac.Nchi);
   DPSIDS_full2 = DPSIDS_full.^2;
   Q=PROF(:,2); % the safety factor from PEST 
   cs=PROF(:,1); % poloidal flux
   
   %cs=cs.^2;%\psi_P
   %ics_plot=0.9;%to plot which radial position
   
  % ics_plot=0.9486;
  ics_plot=0.81;
   [~,ics]=min(abs(cs-ics_plot));
   
   
   T=PROF(:,13);%T value, similar to toroidal field
   T=repmat(T,1,Mac.Nchi);
   T2=T.^2;
   
   B2=T2./R2+DPSIDS_full2.*(dRdchi2+dZdchi2)./jacobian2;
   B2=B2*B0^2;
   
   Dchi=2*pi/(Mac.Nchi -1); %
   B=sqrt(B2);

   %1. get J_p/J_B; J_p= R^2 * B^2 : Jacobian in PEST,  J_B  Jacobian in Boozer
   RB2 = R2 .* B2; % R^2 * B^2  the first row [s1 \chi1,  s1 \chi2,..... s1 \chi_end]
  
   [~,ichi0]=min(abs(Mac.chi) - 0);% index for chi=0
   
   %2. get relatetion between thetaB thetaP; for each flux surface 
   
      for ipsi=1:1:Mac.Ns1
          temp=sum(RB2(ipsi,:),2)*Dchi;
          C_psi(ipsi)=2*pi/temp;
      end
   
      
   for ipsi=1:1:Mac.Ns1 %loop of poloidal flux
   for ii=1:1:Mac.Nchi  % loop of poloidal angle for each flux surface
     if ii<ichi0 
       theta_B_psi(ii)=sum(RB2(ipsi,ii:ichi0),2)*(-1.0*Dchi); %sum of column at same row, case <0. [-pi,pi]
      else
       theta_B_psi(ii)=sum(RB2(ipsi,ichi0+1:ii),2)*Dchi; %sum of column at same row, the cooresponding \theta_B at each flux surface.
     end
   end 
   %C_psi(ipsi)=pi/max(theta_B_psi);
   theta_B_psi=theta_B_psi*C_psi(ipsi);
   theta_B(ipsi,:)=theta_B_psi(:);% theta_B(\psi,\theta)
   
   temp=theta_B_psi(2:end)-theta_B_psi(1:end-1);
   d_theta_B_psi=[temp,temp(end)];% delta theta_B 
   d_theta_B(ipsi,:)=d_theta_B_psi(:);
   hold on;
   %plot(Mac.chi,theta_B_psi*C_psi(ii))
   end
   hold on;
   %plot(Mac.chi,Mac.chi,'-r');
   
   
   
  %3. omega: periodic function determins the transformation 
  
  %omega_0  the integration constant for the \omega
     
      for ipsi=1:1:Mac.Ns1
          temp=(1.0./RB2(ipsi,1:end-1)+1.0./RB2(ipsi,2:end))*0.5.*d_theta_B(ipsi,1:end-1);
         % temp1=(1.0./RB2(ipsi,1)+1.0./RB2(ipsi,2:end))*0.5.*d_theta_B(ipsi,1:end-1);
           temp=-sum(temp,2)/C_psi(ipsi)*Q(ipsi);
          omega_psi_0(ipsi)=2*pi*Q(ipsi)+temp;
         % omega_psi_0(ipsi)=0.0;
      end

  
  
  for ipsi=1:1:Mac.Ns1
   [~,itheta_B_0]=min(abs(theta_B(ipsi,:)) - 0);% index for theta_B=0
   itheta_B_0_psi(ipsi)=itheta_B_0;
   for ii=1:1:Mac.Nchi
       if ii<=itheta_B_0
           temp=(1.0./RB2(ipsi,ii:itheta_B_0)).*(-1.0*d_theta_B(ipsi,ii:itheta_B_0));% at ii=1, the integration is [0, -pi], d_theta_B should be minus
           if theta_B(ipsi,ii) >0
               tem_B=theta_B(ipsi,ii-1) ;
           else
               tem_B=theta_B(ipsi,ii) ;
           end
               
         omega_psi(ii)=-Q(ipsi)*sum(temp,2)/C_psi(ipsi)+Q(ipsi).*tem_B;%0828
         omega_psi1(ii)=-sum(temp,2)/C_psi(ipsi);
           
       else
           temp=(1.0./RB2(ipsi,itheta_B_0+1:ii-1)).*d_theta_B(ipsi,itheta_B_0+1:ii-1);
            if theta_B(ipsi,ii) <0
               tem_B=theta_B(ipsi,ii+1) ;
           else
               tem_B=theta_B(ipsi,ii) ;
           end
           omega_psi(ii)=-Q(ipsi)*sum(temp,2)/C_psi(ipsi)+Q(ipsi).*tem_B;%0828
           omega_psi1(ii)=-sum(temp,2)/C_psi(ipsi);
           
          
           
       end
       
   end
   omega_psi(:)=omega_psi(:)-omega_psi_0(ipsi);
   omega(ipsi,:)=omega_psi(:);
   omega1(ipsi,:)=omega_psi1(:);
   
  % omega(ipsi,:)=omega_psi(:)*2*cs(ipsi);% because s=sqrt(\psi_p)
   %omega1(ipsi,:)=omega_psi1(:)*2*cs(ipsi);
   end
   %omega(:,:)=0;

   
    %4 get X1^UP from PEST   
   [VM1,VM2,VM3,DPSIDS,T] = MacReadVPLASMA([SDIR 'XPLASMA.OUT']);

   
   MN = Mac.Mm; 
   MM = MN - Mac.Mm(1) + 1; % from m1 to m2; m1 the minimum harmonic number, m2 the maximum mode number
   
   s_mode=repmat(Mac.s(1:Mac.Ns1),1,MN(end)-MN(1)+1);
   
      VM1=VM1.*s_mode*2.0; % * 2.0 * sqrt(\psi_p)
      VM2=VM2.*s_mode*2.0;
      VM3=VM3.*s_mode*2.0;
% Note: in PEST(adoped in MARS), the flux label is sqrt of poloidal flux,
% but, in ORBIT, flux label is poloidal flux.Hence, the displacement should
% multiply 2*sqrt(\psi_p)
   
   
   
   %for PEST coordinates
   MN1=Mac.Mm(1);%real harmonic
   MN2=Mac.Mm(end);%harmonics freom PEST are used to make the full displacment
   
   %MN1=8;
  % MN2=8;
   
   MM1=MN1 - Mac.Mm(1) + 1; % label of harmoic 
   MM2=MN2 - Mac.Mm(1) + 1;
  
   %for boozer coordinates
   MN1=-10;
   MN2=40;
   MN1_B=MN1;%real harmonic lower limit in Boozer 
   MN2_B=MN2;%real harmonic ; upper limit in Boozer 
   VM1_B=zeros(Mac.Ns1,MN2_B-MN1_B+1);
   MN_B_real=zeros(1,MN2_B-MN1_B+1);
   % convert it to Boozer coordinates, have assumed the poloidal mode number are same =1. 
   for ipsi=1:1:Mac.Ns1
   for im_B=1:1:MN2_B-MN1_B+1%label of harmonic in Boozer
      Sum_1_s=0.0;
      
      MN_B=im_B+MN1_B-1;%real harmonic in Boozer
      for im=MM1:1:MM2
      %Sum_11=0.0;
      %Sum_12=0.0;
      
      MN_mode=im+Mac.Mm(1)-1;%at beginning= -21
     
      temp= 1.0i * (omega(ipsi,:)*(Q(ipsi)-MN_mode)/Q(ipsi)+MN_mode*theta_B(ipsi,:)-MN_B*theta_B(ipsi,:));
      temp=exp(temp);
      temp2=d_theta_B(ipsi,:);
      
      %integration along the poloidal number 
      
      Sum_1=sum(VM1(ipsi,im)*temp.*temp2,2);%only one m harmonics from PEST, one row, multi columns(Mac.Nchi)
    
      Sum_1_s=Sum_1_s+Sum_1; %sum all of harmonics in PEST coordinates
      end
      
   VM1_B(ipsi,im_B)=Sum_1_s;
   MN_B_real(im_B)=MN_B;
   end 
   end
%VM1_B(radial_coordinates, m number)
%5 choose the harmonices used in 
   r_grids=151;
   r_orbit=linspace(0,1,r_grids);
   
   mmin_orbit=-1;
   mmax_orbit=9;
   
   total_mode_orbit=mmax_orbit-mmin_orbit+1;
   
   mmin=mmin_orbit-MN1_B+1;
   mmax=mmax_orbit-MN1_B+1;
   
   mode_orbit=zeros(r_grids,mmax-mmin+1);
   imm=1;
   for im=mmin:1:mmax
       
       x=Mac.s(1:Mac.Ns1).^2;
       y=VM1_B(:,im)/(2*pi);
       %A=polyfit(x,y,11); 
       %temp=polyval(A,r_orbit);
       xi=r_orbit;
       temp=interp1(x,y,xi);
       mode_orbit(:,imm)=temp;
       imm=imm+1;
   end
   VM1_B_save=VM1_B/2.0/pi;
   
close all   
plot(r_orbit,real(mode_orbit(:,:)),'O')
hold on
plot(Mac.s(1:Mac.Ns1).^2,real( VM1_B_save(:,:)),'-r')
   
   
   save r_orbit r_orbit -ascii
   mode_orbit(1:3,:)=0; %all of perturbation is equal to zero, required by orbit
   [hang,lie]=size(mode_orbit);
   temp=hang*lie;
   temp=reshape(mode_orbit,1,temp);%to reshape to one row , which is needed by orbit 
   temp_real=real(temp);
   temp_real=temp_real*Mac.Norm*Mac.R0EXP;
   
   temp_imag=imag(temp);
   temp_imag=temp_imag*Mac.Norm*Mac.R0EXP;
   amp=abs(temp);
   amp=amp/max(amp);
   save mode_real.dat temp_real -ascii
   save mode_imag.dat temp_imag -ascii
   save amp.dat amp -ascii
   
   hf=figure(1);
   plot(cs, B(:,ichi0),'linewidth',3);
   xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   ylabel('|B|','FontSize',24,'FontWeight','Bold') 
   %legend('\gamma /\omega_A','\omega_r/\omega_A')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   
   
   hf=figure(2);
   m_plot=2;% to plot which harmonic
  
   MM_P= m_plot-Mac.Mm(1)+1; % the corresponding number of elements in perturbation in PEST
   MM_B= m_plot - MN1_B +1; 
   Max_P=max(real(VM1(:,MM_P)));
   Max_B=max(real(VM1_B(:,MM_B)));
   VM1_save=VM1;
   VM1_B_save=VM1_B/2.0/pi;
   VM1(:,MM_P)=VM1(:,MM_P)/Max_P;
   VM1_B(:,MM_B)=VM1_B(:,MM_B)/Max_B;
 
   
   subplot(2,1,1),plot(cs,real(VM1(:,MM_P)),'linewidth',3) ; hold on
   subplot(2,1,1),plot(cs,real(VM1_B(:,MM_B)),'--r','linewidth',3)
   xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   ylabel('Re(VM1)','FontSize',24,'FontWeight','Bold') 
   legend('in PEST','in Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   
   subplot(2,1,2),plot(cs,imag(VM1(:,MM_P)),'-','linewidth',3) ; hold on
   subplot(2,1,2),plot(cs,imag(VM1_B(:,MM_B)),'--r','linewidth',3)
   xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   ylabel('Im(VM1)','FontSize',24,'FontWeight','Bold') 
   legend('in PEST','in Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   
   
   hf=figure(4);
   plot(cs,omega(:,256),'linewidth',3)
   xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   ylabel('\omega','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
  
   hf=figure(401);
   plot(theta_B(ics,:),omega(ics,:),'linewidth',3)
   xlabel('\theta_B','FontSize',24,'FontWeight','Bold')
   ylabel('\omega','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   hf=figure(402);
   plot(theta_B(ics,:),omega1(ics,:),'linewidth',3)
   xlabel('\theta_B','FontSize',24,'FontWeight','Bold')
   ylabel('\omega integration part ','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   hf=figure(403);
   plot(Mac.chi(:),omega(ics,:),'linewidth',3)
   xlabel('\theta_P','FontSize',24,'FontWeight','Bold')
   ylabel('\omega','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   hf=figure(404);
   plot(Mac.chi(:),omega1(ics,:),'linewidth',3)
   xlabel('\theta_P','FontSize',24,'FontWeight','Bold')
   ylabel('\omega integration part ','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
  
   
   
   VM1=VM1_save;
   VM1_B=VM1_B_save;
   
   hf=figure(6);
   
   [~,ics]=min(abs(cs-ics_plot));
   subplot(1,2,1), plot(MN_B_real,real(VM1_B(ics,:)),'-O','linewidth',3)
   xlabel('Poloidal harmonic m','FontSize',24,'FontWeight','Bold')
   ylabel('Real(\xi^{\psi}_{n=1})','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   
   %cs_plot=0.81;%to plot which radial position
   
   subplot(1,2,2),plot(MN_B_real,imag(VM1_B(ics,:)),'-O','linewidth',3)
   %xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   xlabel('Poloidal harmonic m','FontSize',24,'FontWeight','Bold')
   ylabel('Imag(\xi^{\psi}_{n=1})','FontSize',24,'FontWeight','Bold')  
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
    hf=figure(7);
    %to plot which radial position
   [~,ics]=min(abs(cs-ics_plot));
   plot(MN_B_real,abs(VM1_B(ics,:))/max(abs(VM1_B(ics,:))),'-O','linewidth',3)
   %xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   xlabel('Poloidal harmonic m, in Boozer','FontSize',24,'FontWeight','Bold')
   ylabel('ABS(\xi^{\psi}_{n=1})','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
       hf=figure(8);
   plot(cs,C_psi,'-O','linewidth',3);
   xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
  %xlabel('Poloidal harmonic m','FontSize',24,'FontWeight','Bold')
   ylabel('C(\psi)','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
    
   hf=figure(9);
    %to plot which radial position
   plot(Mac.Mm,abs(VM1(ics,:))/max(abs(VM1(ics,:))),'-O','linewidth',3)
   %xlabel('$s=\sqrt{\psi_p}$','FontSize',24,'FontWeight','Bold','interpreter','Latex')
   xlabel('Poloidal harmonic m, PEST','FontSize',24,'FontWeight','Bold')
   ylabel('ABS(\xi^{\psi}_{n=1})','FontSize',24,'FontWeight','Bold') 
   %legend('in_PEST','in_Boozer')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',24,'FontWeight','Bold')
   
   
