% check banana tip for trapped particles
% for REORBIT module

k = 1;

d = load('TEMP_01');
N = round(size(d,1)/5);

if k <= N
   II = 5*(k-1);
   t  = d(II+1,:);  t = t-t(6);
   s  = d(II+2,:);
   c  = d(II+3,:);
   p  = d(II+4,:);
   y  = d(II+5,:);
end


hf=figure(1);
II = 1:length(t);
subplot(2,2,1), plot(t(II),s(II),'b-+'), xlabel('t'), ylabel('s'), hold on,
                plot(t(6),s(6),'bo')
                plot(t(7),s(7),'bs')
                plot(t(5),s(5),'rx')
subplot(2,2,2), plot(t(II),c(II),'b-+'), xlabel('t'), ylabel('c'), hold on,
                plot(t(6),c(6),'bo')
                plot(t(7),c(7),'bs')
                plot(t(5),c(5),'rx')
subplot(2,2,3), plot(t(II),p(II),'b-+'), xlabel('t'), ylabel('p'), hold on,
                plot(t(6),p(6),'bo')
                plot(t(7),p(7),'bs')
                plot(t(5),p(5),'rx')
subplot(2,2,4), plot(s(II),c(II),'b-+'), xlabel('s'), ylabel('c'), hold on,
                plot(s(6),c(6),'bo')
                plot(s(7),c(7),'bs')
                plot(s(5),c(5),'rx')

hf=figure(2);
subplot(1,1,1), plot(t(II),y(II),'b-+'), xlabel('t'), ylabel('y'), hold on,
                plot(t(6),y(6),'bo')
                plot(t(7),y(7),'bs')
                plot(t(5),y(5),'rx')



