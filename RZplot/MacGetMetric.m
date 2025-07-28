function MacGetMetric(R,Z,dRds,dZds,dRdchi,dZdchi,jacobian)
%% plot metrics tensor elements

global Mac

G22 = (dRdchi.^2 + dZdchi.^2);
G33 = R.^2;

figure(Mac.plot_metric)
plot(R(:,1)*Mac.R0EXP,G22(:,1),'r-'), hold on,
N = round(Mac.Nchi/2);
plot(R(:,N)*Mac.R0EXP,G22(:,N),'r-'), hold on,
