function [t_hamada,f_hamada,tc_hamada] = MacGetHamada(R,jacobian)

%compute functions theta(chi) and f(chi) for Hamada coordinates
%defined only in the plasma region

global Mac

N  = Mac.Ns1;
jm = (jacobian(1:N,1:end-1)+jacobian(1:N,2:end))/2;
jm = [zeros(N,1) jm];
jsum = cumsum(jm,2)*(Mac.chi(2)-Mac.chi(1));
j0   = jsum(:,end)*ones(1,Mac.Nchi)/2/pi;
t_hamada = -pi + jsum./j0;
tc_hamada = jacobian(1:N,:)./j0;

jr2 = jacobian(1:N,:)./R(1:N,:).^2;
jm = (jr2(:,1:end-1)+jr2(:,2:end))/2;
jm = [zeros(N,1) jm];
jsum = cumsum(jm,2)*(Mac.chi(2)-Mac.chi(1));
Mac.dpsi(1) = Mac.dpsi(2);
f_hamada = (Mac.q*ones(1,Mac.Nchi)).*(t_hamada+pi) - ((Mac.F./Mac.dpsi)*ones(1,Mac.Nchi)).*jsum;

