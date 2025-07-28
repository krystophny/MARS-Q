mark = '975134';
N    = length(mark);

h  = 1/(N-1);

for k=1:N
    c = [(k-1)*h 4*(k-1)*h*(1-(k-1)*h) 1-(k-1)*h];
    MacRfaCase
    MacRfaXnSurf
end
figure(47)
hl=legend('q_{95}=4.5085','q_{95}=4.7542','q_{95}=4.9334','q_{95}=5.1267','q_{95}=5.3359','q_{95}=5.4469');
set(hl,'FontSize',12)
