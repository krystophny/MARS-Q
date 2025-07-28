% plot eigenvalues from MARS-K runs with outer loop (NSWEEP)
% and extrapolate the converged eigenvalue
format short e

SDIR = '/cscratch/liuy/D3D_IK/Run03/';
NITMAX = 100;

d = load([SDIR 'RESULT.OUT']);
N = d(:,2);
g = d(:,5);
w = d(:,6);

II = find(N==NITMAX);  %data with bad convergence
x  = 1:length(g); 

% find analytic fitting: y=y0 + c/x^p
p = 5; xx=1./x.^p;

hf = figure(1);
c = (g(end)-g(end-1))/(xx(end)-xx(end-1));
g0= g(end)-c*xx(end);
y = g0 + c*xx(end-5:end);
plot(x,g,'b-o','LineWidth',2,'MarkerSize',12), hold on,
plot(x(II),g(II),'bx','LineWidth',2,'MarkerSize',12)
plot(x(1),g(1),'ro','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r')
plot(x(end-5:end),y,'k--')
xlabel('iteration','FontSize',16)
ylabel('\gamma\tau_A','FontSize',16)
ha = get(hf,'CurrentAxes');
set(ha,'FontSize',16)

hf = figure(2);
c = (w(end)-w(end-1))/(xx(end)-xx(end-1));
w0= w(end)-c*xx(end);
y = w0 + c*xx(end-5:end);
plot(x,w,'b-o','LineWidth',2,'MarkerSize',12), hold on,
plot(x(II),w(II),'bx','LineWidth',2,'MarkerSize',12)
plot(x(1),w(1),'ro','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r')
plot(x(end-5:end),y,'k--')
xlabel('iteration','FontSize',16)
ylabel('\omega\tau_A','FontSize',16)
ha = get(hf,'CurrentAxes');
set(ha,'FontSize',16)

hf = figure(3);
plot(g,w,'b-o','LineWidth',2,'MarkerSize',12), hold on,
plot(g(II),w(II),'bx','LineWidth',2,'MarkerSize',12)
plot(g(1),w(1),'ro','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r')
plot(g0,w0,'r+','LineWidth',2,'MarkerSize',12)
xlabel('\gamma\tau_A','FontSize',16)
ylabel('\omega\tau_A','FontSize',16)
ha = get(hf,'CurrentAxes');
set(ha,'FontSize',16)

EigValConv = g0 + i*w0
