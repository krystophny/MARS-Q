load SURFS
R = SURFS(:,1);  
Z = SURFS(:,2);  

N = length(R)/3;
R = reshape(R,N,3);
Z = reshape(Z,N,3);

figure
plot(R,Z,'b-','LineWidth',2), hold on
axis equal
xlabel('R [m]', 'FontSize',16)
ylabel('Z [m]', 'FontSize',16)
