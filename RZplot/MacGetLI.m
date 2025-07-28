% compute LI for a given paralellel current density profile
% note that this script involves running CHEASE code

SDIR1 = '~/IRIS/mars-code/RZplot/';
SDIR2 = '~/Work/RE_AVALANCHE/Data/';
SDIR3 = '/cscratch/liuy/WorkTEMP/';

% generate EXPEQ-file using current density profile stored in SDIR1
% the other profiles and plasma shape from EXPEQ_PART_I and 
% EXPEQ_PART_II stored in SDIR2
% store new EXPEQ to EXPEQ_NEW in SDIR3
d = load([SDIR2 'EXPEQ_PART_II']);
N = round(length(d)/3);
s1= d(1:N);
p1= d(N+1:2*N);

d = load([SDIR1 'JparaSurfRE.txt']);
s2= d(:,1);
j2= d(:,2);
j1= spline(s2,j2,s1);

fid = fopen([SDIR3 'EXPEQ_PART_II'],'w');
fprintf(fid,'%4i\n',[N; 3]);
fprintf(fid,'%16.8e\n',[s1; p1; j1]);
fclose(fid);  

eval(['!cat ' SDIR2 'EXPEQ_PART_I ' SDIR3 'EXPEQ_PART_II > ' SDIR3 'EXPEQ_NEW']);

% run CHEASE with EXPEQ_NEW to calculate LI
copyfile([SDIR2 'Rchease_TEMP'],[SDIR3 'Rchease_TEMP'],'f');
copyfile([SDIR3 'EXPEQ_NEW'],[SDIR3 'EXPEQ'],'f');
cd(SDIR3);
disp('Running CHEASE with EXPEQ-file ...')
eval('!./Rchease_TEMP');

% get LI from log_chease
eval('!grep ''LI       '' log_chease > tmp.txt')
fid = fopen('tmp.txt','r');
ss  = fgetl(fid);  LI = str2num(ss(30:end));
fclose(fid);
disp(['LI = ',num2str(LI)])

cd(SDIR1)


