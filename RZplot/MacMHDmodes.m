% plot sketches of MHD modes
% for normal displacement

x = linspace(0,1,101);

% internal kink
hf = figure(1);
y  = 1 - 1./(1+exp((0.3-x)*50));
plot(x,y,'b-','LineWidth',4),
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModeIK.eps','-depsc')

% external kink
hf = figure(2);
y  = x.^2.*(1.1-x);
plot(x,y,'b-','LineWidth',4),
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModeEK.eps','-depsc')

% tearing mode
hf = figure(3);
xx = x-0.5;
y  = xx.*exp(-xx.^2*1000);
plot(x,y,'b-','LineWidth',4),
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModeTM.eps','-depsc')

% interchange mode
hf = figure(4);
xx = x-0.5;
y  = exp(-xx.^2*1000);
plot(x,y,'b-','LineWidth',4),
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModeIC.eps','-depsc')

% peeling mode
hf = figure(5);
y  = x.^50;
plot(x,y,'b-','LineWidth',4),
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModePM.eps','-depsc')

% infernal mode
hf = figure(6);
y  = x.^3.*exp(-x.^4*4).*(1.2-x);
plot(x,y,'b-','LineWidth',4),
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModeIF.eps','-depsc')

% ballooning mode
hf = figure(7);
x0 = linspace(0.9,1,21); 
a1 = linspace(0.5,1,21);
a2 = a1.*(1.5-a1.^2);
for k=1:length(a1)
    y = x.*exp(-(x-x0(k)).^2*100)*a2(k);
    plot(x,y,'b-','LineWidth',1), hold on
end
xlabel('r','FontSize',18,'FontWeight','Bold'),
ylabel('\xi_n','FontSize',18,'FontWeight','Bold'),
ha=get(hf,'CurrentAxes'); set(ha,'XTick',[],'YTick',[])
print(hf,'MhdModeBM.eps','-depsc')
