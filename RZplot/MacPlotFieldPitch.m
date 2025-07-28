function MacPlotFieldPitch(TEQ,DPSIEQ,R,jacobian,Bn)

global Mac

%find magnetic surface inside plasma, where to plot field lines pitch 
[Y,II] = min(abs(Mac.s-Mac.spitch));

hf = figure(Mac.plot_pitch);

%plot perturbed normal field Bn at Mac.spitch surface
Bn1 = Bn(II,:);
Bn2 = real(exp(i*Mac.n*Mac.phi)*Bn1); Bn2 = transpose(Bn2);
[pp,cc] = meshgrid(Mac.phi,Mac.chi);
pp = pp*180/pi;
cc = cc*180/pi;
pcolor(pp,cc,Bn2), shading interp, hold on 
%colorbar,  
colormap(jet)

%plot equilibrium field line pitch
qloc = (jacobian(II,:)*TEQ(II))./(DPSIEQ(II)*R(II,:).^2);
phi  = cumsum(qloc)*(Mac.chi(2)-Mac.chi(1));

phi0 = linspace(-2*pi,2*pi,21);

for k=1:length(phi0)
    plot(Mac.chi*180/pi,(phi0(k)+phi)*180/pi,'w-','LineWidth',3), hold on,
end
x1 = Mac.chi(1)*180/pi;
x2 = Mac.chi(end)*180/pi;
axis([0 360 x1 x2]);
xlabel('toroidal angle [deg]','FontSize',16,'FontWeight','Bold')
ylabel('poloidal angle [deg]','FontSize',16,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')
