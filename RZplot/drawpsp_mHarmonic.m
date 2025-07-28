% fft for psp

% build psp(chi)
Mac.psp    = [0.12552176E+01 0.3 0.44; 0.12552176E+01 -0.44 -0.3];   % thin new (0.4,-0.4)

k=1;
rc   = Mac.psp(k,1);
chi1 = Mac.psp(k,2)*pi;
chi2 = Mac.psp(k,3)*pi;

k=2;
chi3 = Mac.psp(k,2)*pi+2*pi;
chi4 = Mac.psp(k,3)*pi+2*pi;

Nchi = 320;
current = zeros(Nchi,1);
current(:) = 1e-4;
chi = linspace(0,2*pi,Nchi);
[lchi,II1]=min(abs(chi-chi1));
[lchi,II2]=min(abs(chi-chi2)); 
[lchi,II3]=min(abs(chi-chi3));
[lchi,II4]=min(abs(chi-chi4)); 
current(II1:II2,1)=1;
current(II3:II4,1)=1;

figure(1);
plot(chi,current,'b-','LineWidth',2), hold on
xlabel('\chi','FontSize',16)
ylabel('current','FontSize',16)

Y=fft(current);

j_new=ifft(Y);
figure(2);
plot(chi,j_new,'r--','LineWidth',2), hold on
figure(3);
plot(chi,Y,'r--','LineWidth',2), hold on


fs = 50;   % sampling frequency
Ts = 1/fs;  % sample time
N  = 320;   % length of signal
%N  = 80;
n  = 0:N-1;
t  = n*Ts;  % time vector

Y1 = fft(current,N);
f = (0:N-1)'*fs/N;

figure(4);
plot(f,abs(Y1),'b-','LineWidth',2), hold on

current_new = ifft(Y1);
nn=length(current_new);
chi1 = linspace(0,2*pi,nn);

figure(1);
plot(chi1,current_new,'r--','LineWidth',2), hold on

