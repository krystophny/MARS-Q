function mars_adaptive_packing_routines (new_NRP1,old_equlocate,highresol_equlocate,new_equlocate, eigenfile_dir)
    equ_old = adaptive_load_equ(0,old_equlocate);
    equ_highres = adaptive_load_equ(0,highresol_equlocate);
    cequ_highres = adaptive_mergeEquGrid (equ_highres,0);
    cse_new = adaptive_generate_adaptive_grid(eigenfile_dir,equ_old,cequ_highres,new_NRP1);
    cequ_new = adaptive_interpolate_equ (cequ_highres,cse_new,new_NRP1,0);
    adaptive_output_equ(0,new_equlocate,cequ_new);
end
function component = adaptive_loadeigen (typeid,dir)

    close all;    
    vn = 3;
    if exist ('typeid') == 0 
        filename='XPLASMA.OUT';
    elseif typeid == 0
            filename = 'XPLASMA.OUT';
    elseif typeid == 1
            filename = 'VPLASMA.OUT';
    elseif typeid == 2
            filename ='BPLASMA.OUT';
    elseif typeid == 3
            filename='JPLASMA.OUT';
    elseif typeid == 4
            filename='PPLASMA.OUT';
            vn = 4;
    end

    pfilename=[dir filename];
    fid = fopen (pfilename,'r');
    [name,totCol] = fscanf (fid,'%e',6); 
    MAXM = name(1);
    NRP = name(2);
    for i=1: MAXM
        for k=1:vn
            component(k).Mnum(i,:) = fscanf (fid,'%e',2);
        end
    end
   if ( typeid == 0 || typeid == 1 )
        for j=1:NRP
            DPSIDS(j,:) = fscanf (fid,'%e',3);
            T(j,:) = fscanf (fid,'%e',3);            
        end
    end    
    for i=1:MAXM
        for j=1:NRP
            for k=1:vn
                a = fscanf (fid,'%e',1);
                b = fscanf (fid,'%e',1);
                component(k).val(i,j) = complex(a,b);
            end
        end
    end
    
    fclose (fid);
end


function EQU = adaptive_load_equ(type_id,pfilename)
%    typeid = 1;
    if exist ('type_id') == 0 
        type_id = 0;
    end
    % import plasma equilibrium file
    if (type_id == 0)
        fid=fopen(pfilename,'r');
        EQU.NRP1 = fscanf (fid,'%i',1);
        EQU.NCHI = fscanf (fid,'%i',1);
        EQU.ASPECT = fscanf (fid,'%e',1);
        EQU.R0EXP = fscanf (fid,'%e',1);
        EQU.B0EXP = fscanf (fid,'%e',1);
        %read plasma equilibrium
        for J = 1:(EQU.NRP1 - 1)
                I = J;
                EQUDATAM.CSEM(I)=fscanf (fid,'%e',1);
                PSIISO(2*I) = fscanf (fid,'%e',1);
                EQUDATAM.PEQM(I) = fscanf (fid,'%e',1);
                EQUDATAM.TM(I) = fscanf (fid,'%e',1);
                EQUDATAM.TTPM(I)= fscanf (fid,'%e',1);
                EQUDATAM.PPEQM(I) = fscanf (fid,'%e',1);
                EQUDATAM.DPSIDSM(I) = fscanf (fid,'%e',1);
                EQUDATAM.REQM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZEQM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.RJAM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.G11LM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.G22LM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.G33LM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.G12LM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.RDCDZM(I,:) = fscanf (fid,'%e',EQU.NCHI); 
                EQUDATAM.RDSDZM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.RBZM(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.TPM(I) = EQUDATAM.TTPM(I)/EQUDATAM.TM(I);
                I = J + 1 ;
                EQUDATA.CSE(I)=fscanf (fid,'%e',1);
                PSIISO(2*I-1) = fscanf (fid,'%e',1);
                EQUDATA.PEQ(I) = fscanf (fid,'%e',1);
                EQUDATA.T(I) = fscanf (fid,'%e',1);
                EQUDATA.TTP(I)= fscanf (fid,'%e',1);
                EQUDATA.PPEQ(I) = fscanf (fid,'%e',1);
                EQUDATA.DPSIDS(I) = fscanf (fid,'%e',1);
                EQUDATA.REQ(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.ZEQ(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.RJA(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.G11L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.G22L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.G33L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.G12L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.RDCDZ(I,:) = fscanf (fid,'%e',EQU.NCHI); 
                EQUDATA.RDSDZ(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.RBZ(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATA.TP(I) = EQUDATA.TTP(I)/EQUDATA.T(I);

        end
%         EQU.NRP1 = NRP1;
%         EQU.NCHI = NCHI;
        EQU.PSIISO = PSIISO;
        EQU.EQUDATAM = EQUDATAM;
        EQU.EQUDATA = EQUDATA;
        fclose (fid);
    end    
    if (type_id > 0)
        fid=fopen(pfilename,'r');  % open the vacuum equilibrium file
        EQU.NV = fscanf (fid,'%i',1);
        EQU.NCHI = fscanf (fid,'%i',1);
        EQU.NWBPS = fscanf (fid,'%e',1);
        
        for I=1:EQU.NV
            EQUDATA.ZS1(I) = fscanf (fid,'%e',1);
            EQUDATA.ZRJA(I,:) = fscanf (fid,'%e',EQU.NCHI);
            EQUDATA.ZG11L(I,:) = fscanf (fid,'%e',EQU.NCHI);
            EQUDATA.ZG22L(I,:) = fscanf (fid,'%e',EQU.NCHI);
            EQUDATA.ZG33L(I,:) = fscanf (fid,'%e',EQU.NCHI);
            EQUDATA.ZG12L(I,:) = fscanf (fid,'%e',EQU.NCHI);
            EQUDATA.ZRAN(I,:) = fscanf (fid,'%e',EQU.NCHI);
            EQUDATA.ZZAN(I,:) = fscanf (fid,'%e',EQU.NCHI);
        
            if I < EQU.NV
                EQUDATAM.ZS1(I) = fscanf (fid,'%e',1);
                EQUDATAM.ZRJA(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZG11L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZG22L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZG33L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZG12L(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZRAN(I,:) = fscanf (fid,'%e',EQU.NCHI);
                EQUDATAM.ZZAN(I,:) = fscanf (fid,'%e',EQU.NCHI);
            end
        end
        EQU.EQUDATAM = EQUDATAM;
        EQU.EQUDATA = EQUDATA;
        fclose (fid); 
    end
end         

function adaptive_output_equ(type_id,pfilename,cequ)
%    typeid = 1;
    if exist ('type_id') == 0 
        type_id = 0;
    end
    % export plasma equilibrium file
    if (type_id == 0)
        fid=fopen(pfilename,'w');
        fprintf (fid, '%30i %30i\n',cequ.NRP1,cequ.NCHI);
        fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.ASPECT,cequ.R0EXP,cequ.B0EXP);
        %output plasma equilibrium
        for I = 2:(2*cequ.NRP1 - 1)
                fprintf (fid, '%30.20e %30.20e %30.20e %30.20e\n',cequ.CSE(I),cequ.PSIISO(I),cequ.PEQ(I),cequ.T(I));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.TTP(I),cequ.PPEQ(I),cequ.DPSIDS(I));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.REQ(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.ZEQ(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.RJA(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.G11L(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.G22L(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.G33L(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.G12L(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.RDCDZ(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.RDSDZ(I,:));
                fprintf (fid, '%30.20e %30.20e %30.20e\n',cequ.RBZ(I,:));
        end
        fclose (fid);
    end    
    if (type_id > 0)
        fid=fopen(pfilename,'w');  % open the vacuum equilibrium file
        fprintf (fid, '%30i %30i %30i',cequ.NV,cequ.NCHI,cequ.NWBPS);
        
        for I=1:2*cequ.NV-1
            fprintf (fid, '%30.20e\n', cequ.ZS1(I));
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZRJA(I,:));
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZG11L(I,:));
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZG22L(I,:) );
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZG33L(I,:));
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZG12L(I,:));
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZRAN(I,:));
            fprintf (fid, '%30.20e %30.20e %30.20e\n', cequ.ZZAN(I,:));
        end
        fclose (fid); 
    end
end         

function cequ = adaptive_mergeEquGrid (equ,typeid)
%   typeid == 0 is for the plasma region
%   typeid == 1 is for the vacuum region
% merge the grids of plasma region   
    if (typeid == 0) 
        cequ.NCHI   = equ.NCHI;
        cequ.NRP1   = equ.NRP1;
        cequ.ASPECT = equ.ASPECT;
        cequ.R0EXP  = equ.R0EXP;
        cequ.B0EXP  = equ.B0EXP;
        cequ.CSE    = zeros (1,cequ.NRP1);
        cequ.PEQ    = zeros (1,cequ.NRP1);
        cequ.T      = zeros (1,cequ.NRP1);
        cequ.TTP    = zeros (1,cequ.NRP1);
        cequ.PPEQ   = zeros (1,cequ.NRP1);
        cequ.DPSIDS = zeros (1,cequ.NRP1);
        cequ.TP     = zeros (1,cequ.NRP1);
        cequ.RJA    = zeros (cequ.NRP1,cequ.NCHI);
        cequ.G11L   = zeros (cequ.NRP1,cequ.NCHI);
        cequ.G22L   = zeros (cequ.NRP1,cequ.NCHI);
        cequ.G33L   = zeros (cequ.NRP1,cequ.NCHI);
        cequ.G12L   = zeros (cequ.NRP1,cequ.NCHI);
%        cequ.G11U   = zeros (cequ.NRP1,cequ.NCHI);
%        cequ.G22U   = zeros (cequ.NRP1,cequ.NCHI);
%        cequ.G33U   = zeros (cequ.NRP1,cequ.NCHI);
%        cequ.G12U   = zeros (cequ.NRP1,cequ.NCHI);
        cequ.REQ    = zeros (cequ.NRP1,cequ.NCHI);
        cequ.ZEQ    = zeros (cequ.NRP1,cequ.NCHI);
        cequ.RDCDZ  = zeros (cequ.NRP1,cequ.NCHI);
        cequ.RDSDZ  = zeros (cequ.NRP1,cequ.NCHI);
        cequ.RBZ    = zeros (cequ.NRP1,cequ.NCHI);

        for i=1:equ.NRP1
           k = 2*i-1;
            cequ.CSE(k)    = equ.EQUDATA.CSE(i);
            cequ.PEQ(k)    = equ.EQUDATA.PEQ(i);
            cequ.T(k)      = equ.EQUDATA.T(i);
            cequ.TTP(k)     = equ.EQUDATA.TTP(i);
            cequ.PPEQ(k)   = equ.EQUDATA.PPEQ(i);
            cequ.DPSIDS(k) = equ.EQUDATA.DPSIDS(i);
            cequ.TP(k)     = equ.EQUDATA.TP(i);
            cequ.RJA(k,:)  = equ.EQUDATA.RJA(i,:);
            cequ.G11L(k,:) = equ.EQUDATA.G11L(i,:);
            cequ.G22L(k,:) = equ.EQUDATA.G22L(i,:);
            cequ.G33L(k,:) = equ.EQUDATA.G33L(i,:);
            cequ.G12L(k,:) = equ.EQUDATA.G12L(i,:);
            cequ.REQ(k,:)  = equ.EQUDATA.REQ(i,:);
            cequ.ZEQ(k,:)  = equ.EQUDATA.ZEQ(i,:);
            cequ.RDCDZ(k,:)= equ.EQUDATA.RDCDZ(i,:);
            cequ.RDSDZ(k,:)= equ.EQUDATA.RDSDZ(i,:);
            cequ.RBZ(k,:)  = equ.EQUDATA.RBZ(i,:);
%             if (i == 1)
%                 cequ.G11U(k,:) = zeros(1,equ.NCHI);
%                 cequ.G22U(k,:) = zeros(1,equ.NCHI);
%                 cequ.G33U(k,:) = zeros(1,equ.NCHI);
%                 cequ.G12U(k,:) = zeros(1,equ.NCHI);
%             else
%                 for j=1:equ.NCHI
%                     G11L = cequ.G11L(k,j);
%                     G22L = cequ.G22L(k,j);
%                     G33L = cequ.G33L(k,j);
%                     G12L = cequ.G12L(k,j);
%                    
%                     cequ.G11U(k,j) = G22L / (G11L * G22L - G12L * G12L);
%                     cequ.G22U(k,j) = G11L / (G11L * G22L - G12L * G12L);
%                     cequ.G33U(k,j) = 1.0 / G33L;
%                     cequ.G12U(k,j) = -G12L / (G11L * G22L - G12L * G12L);
%                end
%            end
            if (i < equ.NRP1)
                k = 2*i;
                cequ.CSE(k)    = equ.EQUDATAM.CSEM(i);
                cequ.PEQ(k)    = equ.EQUDATAM.PEQM(i);
                cequ.T(k)      = equ.EQUDATAM.TM(i);
                cequ.TTP(k)    = equ.EQUDATAM.TTPM(i);                
                cequ.PPEQ(k)   = equ.EQUDATAM.PPEQM(i);
                cequ.DPSIDS(k) = equ.EQUDATAM.DPSIDSM(i);
                cequ.TP(k)     = equ.EQUDATAM.TPM(i);
                cequ.RJA(k,:)  = equ.EQUDATAM.RJAM(i,:);
                cequ.G11L(k,:) = equ.EQUDATAM.G11LM(i,:);
                cequ.G22L(k,:) = equ.EQUDATAM.G22LM(i,:);
                cequ.G33L(k,:) = equ.EQUDATAM.G33LM(i,:);
                cequ.G12L(k,:) = equ.EQUDATAM.G12LM(i,:);
                cequ.REQ(k,:)  = equ.EQUDATAM.REQM(i,:);
                cequ.ZEQ(k,:)  = equ.EQUDATAM.ZEQM(i,:);
                cequ.RDCDZ(k,:)= equ.EQUDATAM.RDCDZM(i,:);
                cequ.RDSDZ(k,:)= equ.EQUDATAM.RDSDZM(i,:);
                cequ.RBZ(k,:)  = equ.EQUDATAM.RBZM(i,:);
%                 for j=1:equ.NCHI
%                     G11L = cequ.G11L(k,j);
%                     G22L = cequ.G22L(k,j);
%                     G33L = cequ.G33L(k,j);
%                     G12L = cequ.G12L(k,j);
%                     
%                     cequ.G11U(k,j) = G22L / (G11L * G22L - G12L * G12L);
%                     cequ.G22U(k,j) = G11L / (G11L * G22L - G12L * G12L);
%                     cequ.G33U(k,j) = 1.0 / G33L;
%                     cequ.G12U(k,j) = -G12L / (G11L * G22L - G12L * G12L);
%                 end
            end
        end
        cequ.PSIISO=equ.PSIISO;
    end
 % merge the grids of vacuum region   
    if (typeid == 1)
        cequ.NCHI   = equ.NCHI;
        cequ.NV     = equ.NV * 2 - 1;
        cequ.CSE    = zeros (1,cequ.NV);
        cequ.RJA    = zeros (cequ.NV,cequ.NCHI);
        cequ.G11L   = zeros (cequ.NV,cequ.NCHI);
        cequ.G22L   = zeros (cequ.NV,cequ.NCHI);
        cequ.G33L   = zeros (cequ.NV,cequ.NCHI);
        cequ.G12L   = zeros (cequ.NV,cequ.NCHI);
%        cequ.G11U   = zeros (cequ.NV,cequ.NCHI);
%        cequ.G22U   = zeros (cequ.NV,cequ.NCHI);
%        cequ.G33U   = zeros (cequ.NV,cequ.NCHI);
%        cequ.G12U   = zeros (cequ.NV,cequ.NCHI);
        cequ.REQ    = zeros (cequ.NV,cequ.NCHI);
        cequ.ZEQ    = zeros (cequ.NV,cequ.NCHI);
        
        for i=1:equ.NV
            k = 2*i-1;
            cequ.CSE(k)    = equ.EQUDATA.ZS1(i);
            cequ.RJA(k,:)  = equ.EQUDATA.ZRJA(i,:);
            cequ.G11L(k,:) = equ.EQUDATA.ZG11L(i,:);
            cequ.G22L(k,:) = equ.EQUDATA.ZG22L(i,:);
            cequ.G33L(k,:) = equ.EQUDATA.ZG33L(i,:);
            cequ.G12L(k,:) = equ.EQUDATA.ZG12L(i,:);
            cequ.REQ(k,:)  = equ.EQUDATA.ZRAN(i,:);
            cequ.REZ(k,:)  = equ.EQUDATA.ZZAN(i,:);
            
%             for j=1:equ.NCHI
%                 G11L = cequ.G11L(k,j);
%                 G22L = cequ.G22L(k,j);
%                 G33L = cequ.G33L(k,j);
%                 G12L = cequ.G12L(k,j);
%                     
%                 cequ.G11U(k,j) = G22L / (G11L * G22L - G12L * G12L);
%                 cequ.G22U(k,j) = G11L / (G11L * G22L - G12L * G12L);
%                 cequ.G33U(k,j) = 1.0 / G33L;
%                 cequ.G12U(k,j) = -G12L / (G11L * G22L - G12L * G12L);
%             end
    
            if (i < equ.NV)
                k = 2*i;
                cequ.CSE(k)    = equ.EQUDATAM.ZS1(i);
                cequ.RJA(k,:)  = equ.EQUDATA.ZRJA(i,:);
                cequ.G11L(k,:) = equ.EQUDATAM.ZG11L(i,:);
                cequ.G22L(k,:) = equ.EQUDATAM.ZG22L(i,:);
                cequ.G33L(k,:) = equ.EQUDATAM.ZG33L(i,:);
                cequ.G12L(k,:) = equ.EQUDATAM.ZG12L(i,:);
                cequ.REQ(k,:)  = equ.EQUDATAM.ZRAN(i,:);
                cequ.REZ(k,:)  = equ.EQUDATAM.ZZAN(i,:);               
%                 for j=1:equ.NCHI
%                     G11L = cequ.G11L(k,j);
%                     G22L = cequ.G22L(k,j);
%                     G33L = cequ.G33L(k,j);
%                     G12L = cequ.G12L(k,j);
%                     
%                     cequ.G11U(k,j) = G22L / (G11L * G22L - G12L * G12L);
%                     cequ.G22U(k,j) = G11L / (G11L * G22L - G12L * G12L);
%                     cequ.G33U(k,j) = 1.0 / G33L;
%                     cequ.G12U(k,j) = -G12L / (G11L * G22L - G12L * G12L);
%                 end
            end
        end
    end
end

function cse = adaptive_generate_adaptive_grid(dir,equ,cequ_highres,new_NRP1)
    global sp;
    V_comp = adaptive_loadeigen(1,dir);
    Q_comp  = adaptive_loadeigen(2,dir);
    J_comp  = adaptive_loadeigen(3,dir);
    cs = equ.EQUDATA.CSE;
    csm = equ.EQUDATAM.CSEM;
    
%   integer grid    
    sp.V1 = spline (cs(2:equ.NRP1),V_comp(1).val(:,2:equ.NRP1));
    sp.Q1 = spline (cs(2:equ.NRP1),Q_comp(1).val(:,2:equ.NRP1));
    sp.J2 = spline (cs(2:equ.NRP1),J_comp(2).val(:,2:equ.NRP1));
    sp.J3 = spline (cs(2:equ.NRP1),J_comp(3).val(:,2:equ.NRP1));
%   half integer grid
    sp.V2 = spline (csm(1:equ.NRP1-1),V_comp(2).val(:,1:equ.NRP1-1));
    sp.J1 = spline (csm(1:equ.NRP1-1),J_comp(1).val(:,1:equ.NRP1-1));

    totharmo = size(V_comp(1).val,1)*1;
    y0 = zeros(totharmo,6);
    options = odeset('RelTol',1e-7,'AbsTol',1e-8);
    s0 = cequ_highres.CSE(3);
    s1 = 1;
    [T,Y] = ode45 (@adaptive_integrand,[s0 s1],y0,options);
    T=T';
    
    if (new_NRP1 > size(T,2) + 5)
        disp('Use integral point directly or re-run with higher integral tolerance.');
        cse_new = T;
        cse_new = [0 cse_new];
        new_NRP1 = size(cse_new,2);
    else
%   generate new grid
        cse=zeros(1,new_NRP1);
        T_size = size(T,2);
        step = floor((T_size-2)/(new_NRP1-3));
        smod = mod (T_size-2, new_NRP1-3);
        cse_new(1)=0.0;
        cse_new(2)=s0;
        cse_new(new_NRP1) = 1;
        p1 = 2;
        p2 = step+1;
        for i=3:new_NRP1-2
            smod = smod - 1;
            if (smod <=0)
                smod = 0;
                pmod = 0;
            else
                pmod = 1;
            end
            p2 = p2 + pmod;
            cse_new(i)=mean(T(p1:p2));
            p1 = p2 + 1;
            p2 = p2 + step;
        end
        cse_new(new_NRP1-1) = mean(T(p1:T_size-1));
    end
    cse_new_size=size(cse_new,2);
    cse=zeros(1,2*cse_new_size-1);
    for i=1:cse_new_size
        cse(2*i-1)=cse_new(i);
        if (i <= cse_new_size -1)
            cse(2*i)=(cse_new(i)+cse_new(i+1))*0.5;
        end
    end
    if (cse(2) < cequ_highres.CSE(2))
        disp('Reset half integer point cse(2).');
        cse(2) = cequ_highres.CSE(2);
    end
%     newv1=ppval(sp.V1,cse);
%     figure; plot(cse,newv1,'-+');
%     newv2=ppval(sp.V2,cse);
%     figure; plot(cse,newv2,'-+');
    

end

function dy = adaptive_integrand (t,y)
    global sp;
    V1 = ppval(sp.V1,t);
    V2 = ppval(sp.V2,t);
    Q1 = ppval(sp.Q1,t);
    J1 = ppval(sp.J1,t);
    J2 = ppval(sp.J2,t);
    J3 = ppval(sp.J3,t);
    dy = [V1.' V2.' Q1.' J1.' J2.' J3.'].';

end
function equ_new = adaptive_interpolate_equ (cequ_highres,cse_new,new_NRP1,typeid)

    if (typeid==0)
        msize = 2*new_NRP1-1;
        equ_new.NCHI   = cequ_highres.NCHI;
        equ_new.NRP1   = new_NRP1;
        equ_new.ASPECT = cequ_highres.ASPECT;
        equ_new.R0EXP  = cequ_highres.R0EXP;
        equ_new.B0EXP  = cequ_highres.B0EXP;
        
        equ_new.CSE    = cse_new;
        equ_new.PSIISO = zeros (1,msize);
        equ_new.PEQ    = zeros (1,msize);
        equ_new.T      = zeros (1,msize);
        equ_new.TTP    = zeros (1,msize);
        equ_new.PPEQ   = zeros (1,msize);
        equ_new.DPSIDS = zeros (1,msize);
        equ_new.TP     = zeros (1,msize);
        
        equ_new.PSIISO(2:end) = spline (cequ_highres.CSE(2:end),cequ_highres.PSIISO(2:end),cse_new(2:end));
        equ_new.PEQ(2:end)    = spline (cequ_highres.CSE(2:end),cequ_highres.PEQ(2:end),   cse_new(2:end));
        equ_new.T(2:end)      = spline (cequ_highres.CSE(2:end),cequ_highres.T(2:end),     cse_new(2:end));        
        equ_new.TTP(2:end)    = spline (cequ_highres.CSE(2:end),cequ_highres.TTP(2:end),   cse_new(2:end));        
        equ_new.PPEQ(2:end)   = spline (cequ_highres.CSE(2:end),cequ_highres.PPEQ(2:end),  cse_new(2:end));        
        equ_new.DPSIDS(2:end) = spline (cequ_highres.CSE(2:end),cequ_highres.DPSIDS(2:end),cse_new(2:end));        
        equ_new.TP(2:end)     = spline (cequ_highres.CSE(2:end),cequ_highres.TP(2:end),    cse_new(2:end));        
        
        equ_new.RJA    = zeros (msize,equ_new.NCHI);
        equ_new.G11L   = zeros (msize,equ_new.NCHI);
        equ_new.G22L   = zeros (msize,equ_new.NCHI);
        equ_new.G33L   = zeros (msize,equ_new.NCHI);
        equ_new.G12L   = zeros (msize,equ_new.NCHI);
%        equ_new.G11U   = zeros (msize,equ_new.NCHI);
%        equ_new.G22U   = zeros (msize,equ_new.NCHI);
%        equ_new.G33U   = zeros (msize,equ_new.NCHI);
%        equ_new.G12U   = zeros (msize,equ_new.NCHI);
        equ_new.REQ    = zeros (msize,equ_new.NCHI);
        equ_new.ZEQ    = zeros (msize,equ_new.NCHI);
        equ_new.RDCDZ  = zeros (msize,equ_new.NCHI);
        equ_new.RDSDZ  = zeros (msize,equ_new.NCHI);
        equ_new.RBZ    = zeros (msize,equ_new.NCHI);   
        
        chi = linspace(0,2*pi,cequ_highres.NCHI+1);
        chi = chi(1:end-1);
        [chi_new s_new ] = meshgrid( chi, cse_new );        
        [chi_high s_high] = meshgrid (chi, cequ_highres.CSE);
        equ_new.RJA(2:end,1:end)   = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.RJA(2:end,1:end),   chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.G11L(2:end,1:end)  = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.G11L(2:end,1:end),  chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.G22L(2:end,1:end)  = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.G22L(2:end,1:end),  chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.G33L(2:end,1:end)  = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.G33L(2:end,1:end),  chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.G12L(2:end,1:end)  = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.G12L(2:end,1:end),  chi_new(2:end,1:end), s_new(2:end,1:end),  'spline');                 
        equ_new.REQ(2:end,1:end)   = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.REQ(2:end,1:end),   chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.ZEQ(2:end,1:end)   = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.ZEQ(2:end,1:end),   chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.RDCDZ(2:end,1:end) = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.RDCDZ(2:end,1:end), chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.RDSDZ(2:end,1:end) = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.RDSDZ(2:end,1:end), chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
        equ_new.RBZ(2:end,1:end)   = interp2 ( chi_high(2:end,1:end), s_high(2:end,1:end), cequ_highres.RBZ(2:end,1:end),   chi_new(2:end,1:end), s_new(2:end,1:end),  'spline'); 
%        check
%        figure; plot(equ_new.CSE,equ_new.RJA(:,1));hold on;
%        plot(cequ_highres.CSE,cequ_highres.RJA(:,1),':r');hold on;        
    end
    if (typeid==1)
%       vacuum interpolation will be deveopled later.        
    end

end
