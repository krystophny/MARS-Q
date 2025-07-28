load ../TestCirc/ResMatrix_p.asc
A0 = ResMatrix_p;
load ../TestCirc/VacMatrix_p.asc
A1 = VacMatrix_p;

A = A0 - A1;
M = size(A,1);
A = A(:,1:M) + i*A(:,M+1:end);

figure(1)
x = [-9:19]; y=x; [xx,yy]=meshgrid(x,y);
surfl(xx,yy,abs(A)),
shading flat
