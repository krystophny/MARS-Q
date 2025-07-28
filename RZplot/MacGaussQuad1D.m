function [x,w] = MacGaussQuad1D(n)

%% Compute rule for 1D for interval [-1,1], with the order = n

eps = 2.0e-15; 
xstart = -1.0;  xstop = 1.0;

x = zeros(1,n);
w = zeros(1,n);

m  = round((n+1)/2);
xm = 0.5*(xstop+xstart);
xl = 0.5*(xstop-xstart);
 
if n>1
for j1=1:m
    z = cos(pi*(j1-0.25)/(n+0.5));
    z1 = 2*z;
    while (abs(z-z1)>eps)
      p1=1.0;
      p2=0.0;
      for j2=1:n
        p3=p2;  
        p2=p1;
        p1=((2.0*j2-1.0)*z*p2-(j2-1.0)*p3)/j2;
      end
      pp=n*(z*p1-p2)/(z*z-1);
      z1=z;
      z=z1-p1/pp;
    end
    k=n+1-j1;
    x(j1)=xm-xl*z;
    x(k)=xm+xl*z;
    w(j1)=2*xl/((1-z*z)*pp*pp);
    w(k)=w(j1);
end
end

if n==1
   x(1) = 0.0;
   w(1) = 2.0;
end

%disp(['   GaussQuad1D: x=[' num2str(x) '], w=[',num2str(w) ']'])


