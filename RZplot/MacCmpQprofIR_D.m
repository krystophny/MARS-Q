% compare q-profiles for 4 cases from MAST, D3D, AUG and ITER

SDIR = {};
SDIR{1} = '/home/liuy/Work/IR_D/Database/MAST/25075/3D/n1/PROFEQ.OUT';
SDIR{2} = '/home/liuy/Work/IR_D/Database/DIII-D/157376/3D/n1/PROFEQ.OUT';
SDIR{3} = '/home/liuy/Work/IR_D/Database/AUG/31131.6400/3D/n1/PROFEQ.OUT';
SDIR{4} = '/home/liuy/Work/IR_D/Database/ITER/131025.0001/3D/n1/PROFEQ.OUT';
SC = 'rkbg';

hf = figure(1);
for k=1:size(SDIR,2)
    d = load(SDIR{k});
    plot(d(:,1).^2,d(:,2),'LineWidth',3,'color',SC(k)), hold on,
end

xlabel('\psi_p','FontSize',18,'FontWeight','Bold')
ylabel('safety factor','FontSize',18,'FontWeight','Bold')
ha=get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
legend('MAST#25075','DIII-D#157376','AUG#31131','ITER#131025')
