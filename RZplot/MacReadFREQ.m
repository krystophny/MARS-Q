% plot all frequencies related to the drift kinetic resonances
function MacReadFREQ(filename)

data = load(filename);

s    = data(:,1);
we   = data(:,2);
wsni = data(:,3);
wsne = data(:,4);
wsti = data(:,5);
wste = data(:,6);
awbp = data(:,7);
awbt = data(:,8);
awdi = data(:,9);   awdi(1)=awdi(3); awdi(2)=awdi(3);
awde = data(:,10);
%awda = data(:,11);  awda(1)=awda(3); awda(2)=awda(3); 
awda = data(:,9);  awda(1)=awda(3); awda(2)=awda(3); 
fracp = data(:,12);
fract = data(:,13);

lam1  = data(:,14);
lam0  = data(:,15);
lam2  = data(:,16);
awda1 = data(:,17);
awda2 = data(:,18);

% smoothing precession frequencies
awdi = csaps(s,awdi,0.99999,s);
awda = csaps(s,awda,0.99999,s);

figure(1)
%plot(s,we,'k-','LineWidth',3), hold on,
plot(s,awbt,'b-','LineWidth',3), hold on,
plot(s,wsni+wsti,'b--','LineWidth',3), hold on,
plot(s,awdi,'b:','LineWidth',3), hold on,
%plot(s,awda,'r-','LineWidth',3), hold on,
xlabel('\psi_p^{1/2}','FontSize',18,'FontWeight','Bold'),
ylabel('frequency/\omega_A','FontSize',18,'FontWeight','Bold'),
ha=get(1,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),
%legend('\omega_E','\omega_*^{thi}','<\omega_b^{thi}>','<\omega_d^{thi}>','<\omega_d^{hot}>')
legend('<\omega_b^{thi}>','\omega_*^{thi}','<\omega_d^{thi}>')
%plot([0 1],[0 0],'k:'),

figure(2)
plot(s,fract,'r-','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold'),
ylabel('trapped particle fraction','FontSize',16,'FontWeight','Bold'),
ha=get(2,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),

figure(3)
semilogy(s,abs(we+awda),'r-','LineWidth',2), hold on,
xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold'),
ylabel('|\omega_E+<\omega_d^{hot}>|','FontSize',16,'FontWeight','Bold'),
ha=get(3,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),


figure(4)
plot(s,lam1,'b-','LineWidth',2), hold on,
plot(s,lam2,'b-','LineWidth',2), hold on,
plot(s,lam0,'b--','LineWidth',2), hold on,
s2 = 0.58275;  s0 = 0.8741; a = axis;
plot([s2 s2],[a(3) a(4)],'k--'), hold on,
plot([s0 s0],[a(3) a(4)],'k-.'), hold on,
plot([0 1],[1 1],'r--','LineWidth',2), hold on,
plot([0.2366 0.3914],[1 1],'r-','LineWidth',6), hold on,
plot([s0 1],[1 1],'r-','LineWidth',6), hold on,
xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold'),
ylabel('\Lambda','FontSize',16,'FontWeight','Bold'),
ha=get(4,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),

figure(5)
plot(s,awda1,'r-','LineWidth',2), hold on,
plot(s,awda2,'b-','LineWidth',2), hold on,
s2 = 0.58275;  s0 = 0.8741; a = axis;
plot([s2 s2],[a(3) a(4)],'k--'), hold on,
plot([s0 s0],[a(3) a(4)],'k-.'), hold on,
xlabel('\psi_p^{1/2}','FontSize',16,'FontWeight','Bold'),
ylabel('<\omega_D^\alpha>-omega_E','FontSize',16,'FontWeight','Bold'),
ha=get(5,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold'),


