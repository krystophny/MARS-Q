% go through database 
% and generate RMZM_F_PEST.OUT file
% also copy RMZM_F.OUT to RMZM_F_EQAC.OUT 

SDIR1 = '/cscratch/liuy/Database_More/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';
SDIRR = '/cscratch/liuy/WorkIR_D/AU/';

cd(SDIR1)
DV  = split(ls);
NDV = length(DV)-1;

for kdv=1:NDV
    SDEV = DV{kdv};
    SDIR2 = [SDIR1 SDEV '/'];
    cd(SDIR2)
    SH  = split(ls);
    NSH = length(SH)-1;
    for ksh=1:NSH
        SSHOT = SH{ksh};
        SDIR3 = [SDIR2 SSHOT '/3D/'];
        cd(SDIR3)
        NN  = split(ls);
        NNN = length(NN)-1;
 
        copyfile([SDIR2 SSHOT '/2D/EXPEQ'],[SDIRR 'EXPEQ_SAVE'],'f');        
        
        for knn=1:NNN
            SNUM = NN{knn};
            SDIR4 = [SDIR3 SNUM '/'];
            SDIR = SDIR4;

            copyfile([SDIR4 'RMZM_F.OUT'],[SDIR4 'RMZM_F_EQAC.OUT'],'f');

            d  = load([SDIR2 SSHOT '/2D/global.txt']);
            copyfile([SDIRR 'Rchease_PEST'],[SDIRR 'Rchease'],'f');
            Schease = [SDIRR 'Rchease'];
            eval(['!sed -i s/"NTOR=1"/"NTOR=' SNUM(2) '"/g ' Schease])
            ss = sprintf('%10.4e',d(2));
            eval(['!sed -i s/"B0EXP=2.0000e-00"/"B0EXP=' ss '"/g ' Schease])
            ss = sprintf('%10.4e',d(1));
            eval(['!sed -i s/"R0EXP=1.8000e-00"/"R0EXP=' ss '"/g ' Schease])

            cd(SDIRR)
            eval('!./Rchease');
            eval('!./RmarsQ_quick');
            copyfile([SDIRR 'RMZM_F.OUT'],[SDIR4 'RMZM_F_PEST.OUT'],'f');
        end
    end
end

cd(SDIRM)

