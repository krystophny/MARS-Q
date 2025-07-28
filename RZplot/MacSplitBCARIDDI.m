function MacSplitBCARIDDI(filename)

command = ['!cp ' filename ' MacDataBs'];
eval(command), load MacDataBs,

N = size(MacDataBs,1);
N = floor(N/2);

tmp = MacDataBs(1:N,:);
save MacDataBs_TMP1 tmp -ascii -double

tmp = MacDataBs(N+1:end,:);
save MacDataBs_TMP2 tmp -ascii -double

