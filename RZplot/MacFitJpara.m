function MacFitJpara

global Mac

N = Mac.Ns1;
M = length(Mac.Mm);
s = Mac.s(1:N); s=s(:);

% analytic fitting form
% J(s,m) = A(s)*B(m)
% A(s)   = sum sa(k)*exp[-(s-sn(k))^2/sg(k)^2] 
% B(m)   = sin[mn(1) + mn(2)*m + mn(3)*m^2]/(1+mn(4)*m);
sa = [1.3*exp(i*pi/4*2.5) 2.5*exp(-i*pi/4*1.6)];
sn = [0.97 0.99];
sg = [0.006 0.006];
mn = [pi/16 1.9 0 1];

A = 0;
for k=1:length(sa)
    A = A + sa(k)*exp(-(s-sn(k)).^2/sg(k)^2);
end

B = zeros(1,M);
for k=1:length(Mac.Mm)
    mm = Mac.Mm(k);
    if mm >= 0
       B(k) = sin(mn(1)+mn(2)*mm + mn(3)*mm^2)/(1+mn(4)*mm);
    end
end

Jm = A*B;

SS    = '--';
sleft = 0.95;
mleft = 0;
plot_JM_KD = 4; %=1: 1-D plot along s for Fourier harmonics; 
                %=2: 1-D plot along m for all s
                %=3: 2-D plot in (m,s)-space
                %=4: 2-D plot in (chi,s)-space

if plot_JM_KD==1
   hf=figure(10*Mac.plot_Jpara+1);
   plot(s,real(Jm),SS,'LineWidth',0.5,'Color','b'), hold on,
   xlabel('s','FontSize',18)
   ylabel('Re(J^{||}_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   %a=axis; axis([sleft 1 a(3) a(4)]);

   hf=figure(10*Mac.plot_Jpara+2);
   plot(s,imag(Jm),SS,'LineWidth',0.5,'Color','b'), hold on,
   xlabel('s','FontSize',18)
   ylabel('Im(J^{||}_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   %a=axis; axis([sleft 1 a(3) a(4)]);
end

if plot_JM_KD==2
   hf=figure(10*Mac.plot_Jpara+1);
   plot(Mac.Mm,real(Jm),SS,'LineWidth',0.5,'Color','b'), hold on,
   xlabel('m','FontSize',18)
   ylabel('Re(J^{||}_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   %a=axis; axis([mleft a(2) a(3) a(4)]);

   hf=figure(10*Mac.plot_Jpara+2);
   plot(Mac.Mm,imag(Jm),SS,'LineWidth',0.5,'Color','b'), hold on,
   xlabel('m','FontSize',18)
   ylabel('Im(J^{||}_{mn})','FontSize',18)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18)
   %a=axis; axis([mleft a(2) a(3) a(4)]);
end

if plot_JM_KD==3
   hf=figure(10*Mac.plot_Jpara+3);
   pcolor(Mac.Mm,s,real(Jm)), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('s','FontSize',16)
   title('Re(J^{||}_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_Jpara+4);
   pcolor(Mac.Mm,s,imag(Jm)), colorbar, shading interp,
   xlabel('m','FontSize',16)
   ylabel('s','FontSize',16)
   title('Im(J^{||}_{mn})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([mleft a(2) sleft 1]); %colormap(hot)
end

if plot_JM_KD==4
   expmchi = exp(Mac.Mm*Mac.chi*i);
   Jp = Jm*expmchi;
   hf=figure(10*Mac.plot_Jpara+3);
   pcolor(Mac.chi*180/pi,s,real(Jp)), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('s','FontSize',16)
   title('Re(J^{||})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)

   hf=figure(10*Mac.plot_Jpara+4);
   pcolor(Mac.chi*180/pi,s,imag(Jp)), colorbar, shading interp,
   xlabel('\chi [degree]','FontSize',16)
   ylabel('s','FontSize',16)
   title('Im(J^{||})','FontSize',14)
   ha=get(hf,'CurrentAxes'); set(ha,'FontSize',14)
   a=axis; axis([a(1) a(2) sleft 1]); %colormap(hot)
end





