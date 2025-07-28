%% find x-values from a given data function y(x), 
%% at given array of y-values yy
%% if y(x) is a non-monotonic function, find all the x-values
function [xn,yn]=MacFindX(x,y,yy)

xn = [];
yn = [];

%first do a spline
xs = linspace(x(1),x(end),1001);
ys = spline(x,y,xs);

for k=1:length(yy)
    I = find((ys(1:end-1)-yy(k)).*(ys(2:end)-yy(k)) <= 0);
    for m=1:length(I)
       J = I(m);
       xn = [xn xs(J) + (xs(J+1)-xs(J))*(yy(k)-ys(J))/(ys(J+1)-ys(J))];
       yn = [yn yy(k)];
    end
end

