% plot final results from EFC optimization

SDIRR = '/cscratch/liuy/WorkEFC/Results/';
ICcut = 100;

SS = ['01';'02';'03';'04'];
ICoptA0 = [];
ICoptB0 = [];
ICoptC0 = [];
for k=1:size(SS,1)
    load([SDIRR  'RES_ICoptA_' SS(k,:) '.mat']);
    ICoptA0 = [ICoptA0; ICoptA]; 
    load([SDIRR  'RES_ICoptB_' SS(k,:) '.mat']);
    ICoptB0 = [ICoptB0; ICoptB]; 
    load([SDIRR  'RES_ICoptC_' SS(k,:) '.mat']);
    ICoptC0 = [ICoptC0; ICoptC]; 
end

%II = 1;     %middle row
%II = 2:5;   %U+L
%II = 6:21;  %U+M+L, 1:1:1
%II = 22:37; %U+M+L, 3:1:3
 II = 1:37;

ICoptA0 = ICoptA0(:,II);
ICoptB0 = ICoptB0(:,II);
ICoptC0 = ICoptC0(:,II);

% results for Criterion-A
ICopt = abs(ICoptA0);
II = find(ICopt>ICcut); ICopt(II) = ICcut;
hf = figure(11);
histogram(ICopt);
xlabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
ylabel('numbers of cases','FontSize',16,'FontWeight','Bold'),
title('Criterion-A','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf = figure(12);
plot([1:size(ICopt,1)],ICopt,'b+','LineWidth',0.5,'MarkerSize',8),
xlabel('EF number','FontSize',16,'FontWeight','Bold'),
ylabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
title('Criterion-A','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

ICopt = min(abs(ICoptA0),[],2); ICoptMax(1) = max(ICopt); 
hf = figure(13);
histogram(ICopt);
xlabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
ylabel('number of random EFs','FontSize',16,'FontWeight','Bold'),
title('Criterion-A','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% results for Criterion-B
ICopt = abs(ICoptB0);
II = find(ICopt>ICcut); ICopt(II) = ICcut;
hf = figure(21);
histogram(ICopt);
xlabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
ylabel('numbers of cases','FontSize',16,'FontWeight','Bold'),
title('Criterion-B','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf = figure(22);
plot([1:size(ICopt,1)],ICopt,'b+','LineWidth',0.5,'MarkerSize',8),
xlabel('EF number','FontSize',16,'FontWeight','Bold'),
ylabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
title('Criterion-B','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

ICopt = min(abs(ICoptB0),[],2); ICoptMax(2) = max(ICopt); 
hf = figure(23);
histogram(ICopt);
xlabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
ylabel('number of random EFs','FontSize',16,'FontWeight','Bold'),
title('Criterion-B','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

% results for Criterion-C
ICopt = abs(ICoptC0);
II = find(ICopt>ICcut); ICopt(II) = ICcut;
hf = figure(31);
histogram(ICopt);
xlabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
ylabel('numbers of cases','FontSize',16,'FontWeight','Bold'),
title('Criterion-C','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

hf = figure(32);
plot([1:size(ICopt,1)],ICopt,'b+','LineWidth',0.5,'MarkerSize',8),
xlabel('EF number','FontSize',16,'FontWeight','Bold'),
ylabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
title('Criterion-C','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

ICopt = min(abs(ICoptC0),[],2); ICoptMax(3) = max(ICopt); 
hf = figure(33);
histogram(ICopt);
xlabel('optimal EFCC current [kAt]','FontSize',16,'FontWeight','Bold'),
ylabel('number of random EFs','FontSize',16,'FontWeight','Bold'),
title('Criterion-C','FontSize',16,'FontWeight','Bold'),
ha = get(hf,'CurrentAxes'); set(ha,'FontSize',16,'FontWeight','Bold')

RES = ICoptMax
