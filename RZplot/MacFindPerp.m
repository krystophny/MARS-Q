% given three points O=(R0,Z0), A=(R1,Z1), B=(R2,Z2)
% find point X=(R3,Z3), such that OX is perpendicular to AB

R0=-1;     
Z0=2.5;
R1=0;   
Z1=0;
R2=3;     
Z2=4;

DEL = (R2-R1)^2+(Z2-Z1)^2;
R3  = (R2-R1)^2*R0+(Z2-Z1)^2*R1-(R2-R1)*(Z2-Z1)*(Z1-Z0);
Z3  = (R2-R1)^2*Z1+(Z2-Z1)^2*Z0-(R2-R1)*(Z2-Z1)*(R1-R0);
R3  = R3/DEL;
Z3  = Z3/DEL;

figure
plot([R1 R2],[Z1 Z2],'b-','LineWidth',2), hold on,
plot([R0 R3],[Z0 Z3],'r-','LineWidth',2), hold on,
axis equal


 
