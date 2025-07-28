% generate NN input database for 3D equilibria
% according to the input/output format defined in drawMemo.m 
% also do linear superposition of various rows of RMP coils while varying coil phasing

SDIR1 = '/home/liuy/Work/IR_D/Database/';
SDIRM = '/home/liuy/IRIS/mars-code/RZplot/'; 

%DEV = {'MAST','DIII-D','AUG','ITER'};
DEV = {'DIII-D','ITER'};
SN  = '1234';
LN  = zeros(1,size(DEV,2));

chi1 = linspace(-pi,pi,301);  %internal poloidal mesh for plasma boundary
s1   = linspace(0,1,20);      
s1   = s1.^(1./(2+4*s1));     
s1   = [0 s1(2)/2 s1(2:end)]; %s-mesh for NN input equilibrium profiles 
s2   = linspace(0,1,200);     
s2   = s2.^(1./(2+4*s2));     
s2   = [0 s2(2)/2 s2(2:end)]; %s-mesh for NN output U-vectors 
CP2  = linspace(0,2*pi,9);    %coil phasing for two rows
CP3  = linspace(0,2*pi,9);    %coil phasing for three rows
EP2  = exp(i*CP2(1:8));
EP3  = exp(i*CP3(1:8));


for kd=1:size(DEV,2)
    if strcmp(DEV{kd},'MAST') | strcmp(DEV{kd},'AUG'), SC='UL';  end
    if strcmp(DEV{kd},'DIII-D'),                       SC='ULC'; end
    if strcmp(DEV{kd},'ITER'),                         SC='ULM'; end

    SDIR2 = [SDIR1 DEV{kd} '/'];
    cd(SDIR2)
    F  = split(ls);
    NF = length(F)-1
    cd(SDIRM)

    for kf=1:NF
        FF = F{kf};
        SDIRA = [SDIR2 FF '/3D/'];
        SDIR3 = [SDIRA 'n1/'];

        if 1==0
        %input:RZ: plasma boundary
        d  = load([SDIR3 'RMZM_F.OUT']);
        M  = round(d(1,1));
        N1 = round(d(1,2));
        N2 = round(d(1,3));
        N  = N1 + N2;
        RM = d(N+2:end,1) + d(N+2:end,2)*1i;
        ZM = d(N+2:end,3) + d(N+2:end,4)*1i;

        RM = reshape(RM,N,M);
        ZM = reshape(ZM,N,M);
        RM = RM(N1,:);
        ZM = ZM(N1,:);
        RM(2:end) = 2*RM(2:end);
        ZM(2:end) = 2*ZM(2:end);

        mm = transpose([0:1:M-1]); 
        expmchi = exp(mm*chi1*1i);
        R  = real(RM*expmchi);
        Z  = real(ZM*expmchi);
        
        Rmid = (max(R)+min(R))/2;
        R  = R/Rmid;
        Z  = Z/Rmid;
        Zmid = (max(Z)+min(Z))/2;

        t  = atan2(Z-Zmid,R-1);
        [t,I] = sort(t); R=R(I); Z=Z(I);
        [t,I] = unique(t); R=R(I); Z=Z(I);
        I1 = find(t>pi/2);
        I2 = find(t<-pi/2);
        tt = [t(I1)-2*pi t t(I2)+2*pi];
        RR = [R(I1) R R(I2)];
        ZZ = [Z(I1) Z Z(I2)];
        [tt,I] = unique(tt); RR=RR(I); ZZ=ZZ(I);
	t1 = linspace(-pi,pi,33); t1 = t1(1:end-1);
        R1 = spline(tt,RR,t1);
        Z1 = spline(tt,ZZ,t1);
        
        RZ     = R1;
        RZ(9)  = Z1(9);
        RZ(25) = Z1(25);
        RZ     = [RZ Zmid];  %=Zmid/Rmid
 
        d  = load([SDIR3 'PROFEQ.OUT']);
        s  = d(:,1);

        %input:q: safety factor profile
        q = spline(s,d(:,2),s1(1:end-1));

        %input:P: pressure profile
        P = spline(s,d(:,4),s1);

        %input:rho: density profile
        rho = spline(s,d(:,5),s1);

        %input:w: toroidal rotation frequency profile
        w = spline(s,d(:,6),s1);

        %input:eta: resistivity profile
        eta = spline(s,d(:,7),s1(1:end-1));

        AI = [RZ(:); q(:); P(:); rho(:); w(:); eta(:)];
        AI = [length(AI); length(RZ); length(P); AI];
        save([SDIRA 'NNdata_INPUT2D.dat'],'AI','-ascii') 
        end

        for kn=1:length(SN)
            SDIR4 = [SDIRA 'n' SN(kn) '/'];

            %input dBvac: vacuum normal field along poloidal angle at s=1.01
            %             for separate rows of RMP coils
            dBA = [];
            for kc=1:length(SC)
            SDIR5 = [SDIR4 SC(kc) '/'];
            d  = load([SDIR5 'dBnormal.txt']);
            t  = d(:,1);
            y  = d(:,2) + d(:,3)*1i;

            [t,I] = sort(t);   y=y(I);
            [t,I] = unique(t); y=y(I);
            I1 = find(t>pi/2);
            I2 = find(t<-pi/2);
            tt = [t(I1)-2*pi; t; t(I2)+2*pi];
            yy = [y(I1); y; y(I2)];
            [tt,I] = unique(tt); yy=yy(I);
	    t1 = linspace(-pi,pi,101); t1 = t1(1:end-1);
            dBvac = spline(tt,yy,t1); dBvac = dBvac(:);
            dBA = [dBA dBvac];

            %input:AI: dBvac
            AI = [real(dBvac); imag(dBvac)];
            AI = [length(AI); length(dBvac); AI];
            save([SDIR5 'NNdata_INPUT3D.dat'],'AI','-ascii') 
            end

            %input dBvac: vacuum normal field along poloidal angle at s=1.01
            %             for linear superposition of two rows
            if length(SC)==2
            for kp=1:length(EP2)
                SDIR6 = [SDIR4 int2str(kp) '/'];
                mkdir(SDIR6);
                dBN = dBA(:,1)*EP2(kp)+dBA(:,2);

                %input:AI: dBN
                AI = [real(dBN); imag(dBN)];
                AI = [length(AI); size(dBN,1); AI];
                save([SDIR6 'NNdata_INPUT3D.dat'],'AI','-ascii') 
            end
            end

            %input dBvac: vacuum normal field along poloidal angle at s=1.01
            %             for linear superposition of three rows
            if length(SC)==3
            for kp1=1:length(EP3)
            for kp3=1:length(EP3)
                SDIR6 = [SDIR4 int2str(kp1) int2str(kp3) '/'];
                mkdir(SDIR6);
                dBN = dBA(:,1)*EP3(kp1)+dBA(:,2)+dBA(:,3)*EP3(kp3);

                %input:AI: dBN
                AI = [real(dBN); imag(dBN)];
                AI = [length(AI); size(dBN,1); AI];
                save([SDIR6 'NNdata_INPUT3D.dat'],'AI','-ascii') 
            end
            end
            end
         end
     end
end

