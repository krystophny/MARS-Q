function [xx,yy,ww] = MacGaussQuad2D(n)

%% Define Gauss quadrature rule for a rectangle (-1,1)x(-1,1)

[x,w] = MacGaussQuad1D(n);

% compute Gauss quadrature rule for 2D
[xx,yy] = meshgrid(x,x); 
ww = w'*w;
xx = xx(:);  yy = yy(:);  ww = ww(:);

