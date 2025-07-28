% generate PROF*.IN files from PROFEQ.OUT stored in IR_D database

SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';
SDIR1 = '/cscratch/liuy/Database_More/';
SDIRR = '/cscratch/liuy/WorkIR_D/AU/';
SSHOT = 'EAST/94048/';
SCASE = '3D/n4/';

SDIR  = [SDIR1 SSHOT SCASE];

d = load([SDIR 'PROFEQ.OUT']);
s = d(:,1);
den = d(:,5);
rot = d(:,6);
Ti  = d(:,10);
Te  = d(:,11);

RES = [[length(s) 1]; s den];
save([SDIRR 'PROFDEN_SAVE'],'RES','-ascii');

RES = [[length(s) 1]; s rot];
save([SDIRR 'PROFROT_SAVE'],'RES','-ascii');

RES = [[length(s) 1]; s Ti];
save([SDIRR 'PROFTI_SAVE'],'RES','-ascii');

RES = [[length(s) 1]; s Te];
save([SDIRR 'PROFTE_SAVE'],'RES','-ascii');

copyfile([SDIR1 SSHOT '2D/EXPEQ'],[SDIRR 'EXPEQ_SAVE'],'f');
