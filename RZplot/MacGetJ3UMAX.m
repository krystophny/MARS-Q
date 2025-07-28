function MacGetJ3UMAX(filename)

global Mac

command = ['!cp ' filename ' J3UMAX'];
eval(command), load J3UMAX,
data = J3UMAX;

q = data(:,1);
jr = data(:,2);
ji = data(:,3);

figure(Mac.plot_J3U)
plot(q,jr,'r+',q,ji,'b+'), hold on,

N = 20;
a = linspace(0,1,N+1);
qs= q(1:end-1);
for k=2:N
  qs = [qs; (1-a(k))*q(1:end-1)+a(k)*q(2:end)];
end
qs = sort(qs);

jrs = spline(q,jr,qs);
jis = spline(q,ji,qs);
II = find(jrs<=0); jrs(II) = min(jr); 
II = find(jis<=0); jis(II) = min(ji); 

plot(qs,jrs,'r-',qs,jis,'b-'), hold on,
grid on

% get local maxima
eps = 0.1;
qm = [2 3 4 5 6];
jm = max(jrs,jis);
qqm = [];
for k=1:length(qm)
  I1 = find(qs>qm(k)-eps & qs<qm(k));
  [y,I2] = max(jm(I1));
  qqm = [qqm; qs(I1(I2))];

  I1 = find(qs<qm(k)+eps & qs>qm(k));
  [y,I2] = max(jm(I1));
  qqm = [qqm; qs(I1(I2))];
end

q_res = qqm
