function MacCheckJacob(filename,jacobRZ)

global Mac

command = ['!cp ' filename ' JACOBIAN'];
eval(command), load JACOBIAN

Mac.Nm1 = JACOBIAN(1,1);

if Mac.Ns1 > JACOBIAN(1,2),
   disp('Number of radial points is different for equilibrium and stability!')
   Ns = JACOBIAN(1,2);
   Mac.Ns = Ns;
   Mac.s  = Mac.s(1:Ns);
else 
   Ns1 = JACOBIAN(1,2);
end

%Mac.Mm = round(JACOBIAN(2:Mac.Nm1+1,1));

JACM = JACOBIAN((2+Mac.Nm1*1):end,1) + JACOBIAN((2+Mac.Nm1*1):end,2)*i;

JACM = reshape(JACM,Ns1,Mac.Nm1);
JACM(:,2:end) = 2*JACM(:,2:end);

JACM = JACM(1:Mac.Ns1,:);

m = [0:1:Mac.Nm1-1]';
expmchi = exp(m*Mac.chi*i);

JAC = real(JACM*expmchi);

% compare JAC and jacobRZ
jacerr = abs(JAC(2:end,:)-jacobRZ(2:Mac.Ns1,:))./abs(JAC(2:end,:));
maxerr = max(max(jacerr));
meanerr = sum(sum(jacerr))/size(jacerr,1)/size(jacerr,2);

disp(['   Check jacobian: maximum error = ' num2str(maxerr) ', average error = ' num2str(meanerr)])

if Mac.plot_JACOB>0
   figure(Mac.plot_JACOB)
   subplot(3,1,1), plot(Mac.s(1:Mac.Ns1),real(JACM)), hold on,
                   ylabel('Re(JAC_m)','FontSize',14)
   subplot(3,1,2), plot(Mac.s(1:Mac.Ns1),imag(JACM)), hold on,
                   xlabel('s','FontSize',14), ylabel('Im(JAC_m)','FontSize',14)
   subplot(3,1,3),
   plot(Mac.s(2:Mac.Ns1),JAC(2:end,1)./jacobRZ(2:Mac.Ns1,1),'r-'), hold on
   axis([0 1 0.5 1.5])
   xlabel('s','FontSize',14)
   ylabel('JacobCHEASE/JacobRZ','FontSize',14)
end

if Mac.plot_JACOB0
   figure(Mac.plot_JACOB0)
   plot(Mac.chi,JAC(Mac.Ns1,:),'b--'), hold on,
end
