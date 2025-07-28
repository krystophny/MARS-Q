% go through database 
% generate and save MoR/SVD data

global Mac

SDIR1 = '/cscratch/liuy/WorkIR_D/Database/DIII-D/';
%SDIR1 = '/home/liuy/Work/IR_D/Database/DIII-D/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';

cd(SDIR1)
F  = split(ls);
NF = length(F)-1;
SC = 'C';

cd(SDIRM)
for kf=1:NF
    FF = F{kf};
    cd([SDIR1 FF '/3D/'])
    G  = split(ls);
    cd(SDIRM)
    NG = length(G)-1;
    for kn=1:NG
    GG = G{kn};
    for kc=1:length(SC)
        SDIR  = [SDIR1 FF '/3D/' GG '/'];
        SDIR2 = [SC(kc) '/']; 

        MacMainSvdIR_D

        res1 = [1 size(Mac.SVD_JU,1) size(Mac.SVD_JV,1)];
        res2 = [res1 zeros(1,length(Mac.SVD_JS)*2-3); transpose([Mac.SVD_JS; Mac.SVD_JS(1:end-2); Mac.SVD_JRER; Mac.SVD_JRE]); real(Mac.SVD_JU) imag(Mac.SVD_JU); real(Mac.SVD_JV) imag(Mac.SVD_JV)];

        save tmp.txt res2 -ascii
        copyfile('tmp.txt',[SDIR SDIR2 'SVD_Jpara5.txt'],'f')   

        res1 = [1 size(Mac.SVD_XU,1) size(Mac.SVD_XV,1)];
        res2 = [res1 zeros(1,length(Mac.SVD_XS)*2-3); transpose([Mac.SVD_XS; Mac.SVD_XS(1:end-2); Mac.SVD_XRER; Mac.SVD_XRE]); real(Mac.SVD_XU) imag(Mac.SVD_XU); real(Mac.SVD_XV) imag(Mac.SVD_XV)];

        save tmp.txt res2 -ascii
        copyfile('tmp.txt',[SDIR SDIR2 'SVD_X1U5.txt'],'f')   
     end     
     end
end
