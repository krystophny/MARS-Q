x = [0 0.1 0.12 0.14];
y = [4.0 4.3 4.7 5.2];

xx=linspace(0,1,101);
yy=spline(x,y,xx);

plot(x,y,'o',xx,yy,'-')
