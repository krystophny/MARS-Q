function [PM1,PM2,PM3,PM4,PM5] = MacReadPPLASMA(filename)

global Mac

fac_wk = 1.0;

PPLASMA = load(filename);

Mac.Nm1 = PPLASMA(1,1);

if Mac.Ns1 ~= PPLASMA(1,2),
    disp('Number of radial points is different for equilibrium and stability!')
end

Mac.n  = round(PPLASMA(1,3));
Mac.Mm = round(PPLASMA(2:Mac.Nm1+1,1));

if length(Mac.mm_plot)>0 
   Mac.mm_true = Mac.mm_plot-Mac.Mm(1)+1;
else
   Mac.mm_true = Mac.Mm-Mac.Mm(1)+1;
end

PM1 = PPLASMA(Mac.Nm1+2:end,1) + PPLASMA(Mac.Nm1+2:end,2)*i;  %PRE
PM2 = PPLASMA(Mac.Nm1+2:end,5) + PPLASMA(Mac.Nm1+2:end,6)*i;  %PPERP
PM3 = PPLASMA(Mac.Nm1+2:end,7) + PPLASMA(Mac.Nm1+2:end,8)*i;  %PPARA
PM4 = PPLASMA(Mac.Nm1+2:end,3) + PPLASMA(Mac.Nm1+2:end,4)*i;  %PEE

PM1 = reshape(PM1,Mac.Ns1,Mac.Nm1)*Mac.PNORM;
PM2 = reshape(PM2,Mac.Ns1,Mac.Nm1)*Mac.PNORM;
PM3 = reshape(PM3,Mac.Ns1,Mac.Nm1)*Mac.PNORM;
PM4 = reshape(PM4,Mac.Ns1,Mac.Nm1)*Mac.PNORM;

PM5 = PM4;
if size(PPLASMA,2)==10
   PM5 = PPLASMA(Mac.Nm1+2:end,9) + PPLASMA(Mac.Nm1+2:end,10)*i;  %DPHI
   PM5 = reshape(PM5,Mac.Ns1,Mac.Nm1)*Mac.FNORM;
end


% for KSTAR where deltaPe is measured, adding a factor
if 1==0
   Ti0  = 2.75069e+03; 
   Te0  = 2.90103e+03;
   facP = Te0/(Te0+Ti0);
   PM1  = PM1*facP;
end

% note that PM* are defined at half-points, recompute at integer-points
if 1==0
x = (Mac.s(1:Mac.Ns1-1) + Mac.s(2:Mac.Ns1))*0.5;
PM1new = PM1; PM1new(2:end-1,:) = transpose(spline(x',transpose(PM1(1:end-1,:)),Mac.s(2:Mac.Ns1-1)'));
PM1new(1,:) = PM1new(2,:);  PM1new(end,:) = PM1new(end-1,:);  PM1 = PM1new; clear PM1new

PM2new = PM2; PM2new(2:end-1,:) = transpose(spline(x',transpose(PM2(1:end-1,:)),Mac.s(2:Mac.Ns1-1)'));
PM2new(1,:) = PM2new(2,:);  PM2new(end,:) = PM2new(end-1,:);  PM2 = PM2new; clear PM2new

PM3new = PM3; PM3new(2:end-1,:) = transpose(spline(x',transpose(PM3(1:end-1,:)),Mac.s(2:Mac.Ns1-1)'));
PM3new(1,:) = PM3new(2,:);  PM3new(end,:) = PM3new(end-1,:);  PM3 = PM3new; clear PM3new
end

k0=find(Mac.s<=Mac.core); k0=k0(end);
for k=1:k0
  PM1(k,:) = PM1(k0+1,:);
  PM2(k,:) = PM2(k0+1,:);
  PM3(k,:) = PM3(k0+1,:);
  PM4(k,:) = PM4(k0+1,:);
  PM5(k,:) = PM5(k0+1,:);
end

if 1==1
I = find(Mac.s>Mac.edge & Mac.s<=1.0); 
PM1(I,:) = 0;
PM2(I,:) = 0;
PM3(I,:) = 0;
PM4(I,:) = 0;
PM5(I,:) = 0;
end

if 1==0
PM1(end,:) = 2*PM1(end-1,:)-PM1(end-2,:);
PM2(end,:) = 0;
PM3(end,:) = 0;
PM4(end,:) = 0;
PM5(end,:) = 0;
end

if Mac.plot_PM > 0
   if 1==1
   figure(10*Mac.plot_PM+0)
   SS = Mac.SS;
   PM = PM1; a = angle(PM); C = max(max(abs(PM)));
   I = find(abs(PM)/C<0.2); a(I)=0; I = find(a<0); a(I)=a(I)+pi;
   subplot(3,2,1), plot(Mac.s(1:Mac.Ns1),real(PM),SS), hold on,
                   ylabel('Re(P_{fluid})','FontSize',14)
   subplot(3,2,2), plot(Mac.s(1:Mac.Ns1),imag(PM),SS), hold on,
                   xlabel('s','FontSize',14), 
                   ylabel('Im(P_{fluid})','FontSize',14)

   PM = PM4; a = angle(PM); C = max(max(abs(PM)));
   I = find(abs(PM)/C<0.2); a(I)=0; I = find(a<0); a(I)=a(I)+pi;
   subplot(3,2,3), plot(Mac.s(1:Mac.Ns1),real(PM),SS), hold on,
                   ylabel('Re(P_e)','FontSize',14)
   subplot(3,2,4), plot(Mac.s(1:Mac.Ns1),imag(PM),SS), hold on,
                   xlabel('s','FontSize',14), 
                   ylabel('Im(P_e)','FontSize',14)

   PM = PM5; a = angle(PM); C = max(max(abs(PM)));
   %I = find(abs(PM)/C<0.2); a(I)=0; I = find(a<0); a(I)=a(I)+pi;
   subplot(3,2,5), plot(Mac.s(1:Mac.Ns1),real(PM),SS), hold on,
                   ylabel('Re(P_{||})','FontSize',14)
   subplot(3,2,6), plot(Mac.s(1:Mac.Ns1),imag(PM),SS), hold on,
                   xlabel('s','FontSize',14), 
                   ylabel('Im(P_{||})','FontSize',14)
   end

   
   kf = 10*Mac.plot_PM;
   hf=figure(kf+1);
   Cn=max(max(abs(PM1))); 
   Cn=1.0;
   plot(Mac.s(1:Mac.Ns1),real(PM1)/Cn,'-','LineWidth',1), hold on,
   ylabel('Re[p_f^{(m)}]','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight
   
   hf=figure(kf+2);
   plot(Mac.s(1:Mac.Ns1),imag(PM1)/Cn,'-','LineWidth',1), hold on,
   ylabel('Im[p_f^{(m)}]','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight

   if 1==0
   s0=Mac.s(1:Mac.Ns1);
   s1=s0(1:end-1);
   P20 = imag(PM2);
   P21 = imag(PM2(1:end-1,:));
   for k=1:size(PM2,2)
      P20(:,k) = spline(s1,P21(:,k),s0);
   end
   plot(s0,P20)
   end

   if 1==0
   hf=figure(kf+3);
   plot(Mac.s(1:Mac.Ns1),imag(PM3)*fac_wk,SS,'LineWidth',1), hold on,
   ylabel('Im[p_{||}^{(m)}]','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight

   hf=figure(kf+4);
   plot(Mac.s(1:Mac.Ns1),real(PM1+PM4),'r-','LineWidth',1), hold on,
   plot(Mac.s(1:Mac.Ns1),real(PM2),'b-','LineWidth',1), hold on,
   plot(Mac.s(1:Mac.Ns1),real(PM3),'k-','LineWidth',1), hold on,
   ylabel('Re[p^{(m)}]','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight

   hf=figure(kf+5);
   plot(Mac.s(1:Mac.Ns1),imag(PM1+PM4),'r-','LineWidth',1), hold on,
   plot(Mac.s(1:Mac.Ns1),imag(PM2),'b-','LineWidth',1), hold on,
   plot(Mac.s(1:Mac.Ns1),imag(PM3),'k-','LineWidth',1), hold on,
   ylabel('Im[p^{(m)}]','FontSize',16,'FontWeight','Bold')
   xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight


   hf=figure(kf+6);
   [Y,I]=min(abs(Mac.s-0.5));
   plot(Mac.Mm,imag(PM1(I,:)+PM4(I,:)),'r-o','LineWidth',1), hold on,
   plot(Mac.Mm,imag(PM2(I,:)),'b-s','LineWidth',1), hold on,
   plot(Mac.Mm,imag(PM3(I,:)),'k-d','LineWidth',1), hold on,
   ylabel('P_{\perp}','FontSize',16,'FontWeight','Bold')
   xlabel('m','FontSize',16,'FontWeight','Bold')
   ha=get(hf,'CurrentAxes');
   set(ha,'FontSize',16,'FontWeight','Bold')
   axis tight
   end
end

if Mac.plot_PMM > 0
   figure(Mac.plot_PMM)
   C=max(max(abs(imag(PM1)),[],1));
   C1=max(imag(PM1),[],1);
   C2=min(imag(PM1),[],1);
   C3=C1;
   I=find(abs(C1)<abs(C2));  C3(I)=C2(I);
   subplot(3,1,1), plot(Mac.Mm,C3,'b-o','LineWidth',2,'MarkerSize',6), hold on,
                   ylabel('harmonic amplitude','FontSize',14)
   C=max(max(abs(real(PM2)),[],1));
   C1=max(real(PM2),[],1);
   C2=min(real(PM2),[],1);
   C3=C1;
   I=find(abs(C1)<abs(C2));  C3(I)=C2(I);
   subplot(3,1,2), plot(Mac.Mm,C3,'b-o','LineWidth',2,'MarkerSize',6), hold on,
                   ylabel('harmonic amplitude','FontSize',14)
   C=max(max(abs(real(PM3)),[],1));
   C1=max(real(PM3),[],1);
   C2=min(real(PM3),[],1);
   C3=C1;
   I=find(abs(C1)<abs(C2));  C3(I)=C2(I);
   subplot(3,1,3), plot(Mac.Mm,C3,'b-o','LineWidth',2,'MarkerSize',6), hold on,
                   ylabel('harmonic amplitude','FontSize',14)

                   xlabel('m','FontSize',14), 
end

