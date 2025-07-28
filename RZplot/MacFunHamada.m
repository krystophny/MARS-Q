function [gm] = MacGetHamada(g,t_hamada,f_hamada,tc_hamada,mk)

%compute Fourier harmonics of function g(s,chi) in 
%Hamada coordinates (V,theta,zeta)
%defined only in the plasma region

global Mac

N   = Mac.Ns1;
ga  = g(1:N,:).*tc_hamada.*exp(-i*Mac.n*f_hamada);
gam = (ga(:,1:end-1)+ga(:,2:end))/2;
tm  = (t_hamada(:,1:end-1)+t_hamada(:,2:end))/2;

gm = zeros(N,length(mk));

for k=1:length(mk)
    gm(:,k) = sum(gam.*exp(-i*mk(k)*tm),2);
end
gm = gm*(Mac.chi(2)-Mac.chi(1))/2/pi;

