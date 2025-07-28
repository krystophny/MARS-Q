% go through database 
% obtain RMP ELM suppression criteria given the threshold coil current as measured in experiments

global SDIR SDIR5 k_plot

SDIR1 = '/cscratch/liuy/Database_More/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/';
k_plot= 0;

cd(SDIR1)
DV  = split(ls);
NDV = length(DV)-1;

DV  = {'KSTAR'};
NDV = 1;

for kdv=1:NDV
    SS = DV{kdv};
    SDIR2 = [SDIR1 SS '/'];
    cd(SDIR2)
    SH  = split(ls);
    NSH = length(SH)-1;
    for ksh=1:NSH
        SS = SH{ksh};
        SDIR3 = [SDIR2 SS '/3D/'];
        cd(SDIR3)
        NN  = split(ls);
        NNN = length(NN)-1;
        cd(SDIRM)
        for knn=1:NNN
            SS = NN{knn};
            SDIR4 = [SDIR3 SS '/'];
            SDIR = SDIR4;
            if strcmp(DV{kdv},'DIII-D'), CL = 'UL';  end  
            if strcmp(DV{kdv},'EAST'),   CL = 'UL';  end  
            if strcmp(DV{kdv},'KSTAR'),  CL = 'ULM'; end
            if strcmp(DV{kdv},'AUG'),    CL = 'UL';  end  
            NCL=length(CL); 
            for kcl=1:NCL
                SS = CL(kcl);
                SDIR5 = [SDIR4 SS '/'];

                MacMainRMPthreshold
                MacRfaCtRMPthreshold            
            end

            %linear superposition of criterion quantities
            %first need to run Work/drawCurrThreshold.m to save 'current_threshold.txt'
            CC = load([SDIR4 'current_threshold.txt']);
            NCC = round(size(CC,1)/2);
            RES1 = [];
            RES2 = [];
            RES1M = [];
            RES2M = [];
            for kcc=1:NCC
            kc2  = (kcc-1)*2;
            CRIT1U = 0;
            CRIT1L = 0;
            CRIT2 = 0;
            CRIT1UM = 0;
            CRIT1LM = 0;
            CRIT2M = 0;
            for kcl=1:NCL
                SS    = CL(kcl);
                SDIR5 = [SDIR4 SS '/'];
 
                Y       = load([SDIR5 'CRIT_Xn']); 
                Z       = Y(1) + i*Y(2);
                CRIT1U  = CRIT1U + Z*(CC(1+kc2,kcl)+i*CC(2+kc2,kcl));
                Z       = Y(3) + i*Y(4);
                CRIT1L  = CRIT1L + Z*(CC(1+kc2,kcl)+i*CC(2+kc2,kcl));

                Y       = load([SDIR5 'CRIT_Xn_mean']); 
                Z       = Y(1) + i*Y(2);
                CRIT1UM = CRIT1UM + Z*(CC(1+kc2,kcl)+i*CC(2+kc2,kcl));
                Z       = Y(3) + i*Y(4);
                CRIT1LM = CRIT1LM + Z*(CC(1+kc2,kcl)+i*CC(2+kc2,kcl));

                Y       = load([SDIR5 'CRIT_b1res']); 
                Z       = Y(1) + i*Y(2);
                CRIT2 = CRIT2 + Z*(CC(1+kc2,kcl)+i*CC(2+kc2,kcl));

                Y       = load([SDIR5 'CRIT_b1res_mean']); 
                Z       = Y(1) + i*Y(2);
                CRIT2M  = CRIT2M + Z*(CC(1+kc2,kcl)+i*CC(2+kc2,kcl));
            end
            CRIT1  = CRIT1U; if abs(CRIT1U)<abs(CRIT1L); CRIT1=CRIT1L; end
            CRIT1M = CRIT1UM; if abs(CRIT1UM)<abs(CRIT1LM); CRIT1M=CRIT1LM; end
            RES1   = [RES1; real(CRIT1) imag(CRIT1)];
            RES2   = [RES2; real(CRIT2) imag(CRIT2)];
            RES1M  = [RES1M; real(CRIT1M) imag(CRIT1M)];
            RES2M  = [RES2M; real(CRIT2M) imag(CRIT2M)];
            end
            save([SDIR4 'CRIT_Xn.txt'],'RES1','-ascii') 
            save([SDIR4 'CRIT_b1res.txt'],'RES2','-ascii') 
            save([SDIR4 'CRIT_Xn_mean.txt'],'RES1M','-ascii') 
            save([SDIR4 'CRIT_b1res_mean.txt'],'RES2M','-ascii') 
        end
    end
end

cd(SDIRM)

