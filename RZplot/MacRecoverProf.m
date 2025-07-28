% recover PROF*.IN MARS-F input files from PROFEQ.OUT

SDIR = '/cscratch/liuy/WorkIR_D/BU/';

d = load([SDIR 'PROFEQ.OUT']);
s = d(:,1);
II = [5 6 10 11];

N = 101;
ss = linspace(0,1,N); ss=ss(:);
Y  = spline(s,transpose(d(:,II)),ss); Y=transpose(Y);

RES = [N 1; [ss Y(:,1)]];
save([SDIR 'PROFDEN_SAVE'],'RES','-ascii')

RES = [N 1; [ss Y(:,2)]];
save([SDIR 'PROFROT_SAVE'],'RES','-ascii')

RES = [N 1; [ss Y(:,3)]];
save([SDIR 'PROFTI_SAVE'],'RES','-ascii')

RES = [N 1; [ss Y(:,4)]];
save([SDIR 'PROFTE_SAVE'],'RES','-ascii')
