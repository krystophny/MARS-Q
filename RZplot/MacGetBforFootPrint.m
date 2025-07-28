function MacGetBforFootPrint(R,Z,Br,Bz)

global Mac
global SDIR

% define flux surfaces to do the conversion
% from MARS-F poloidal angle to geometric angle
sfp_max = 2.6;  %for MAST25056 to cover bounding box [0.2 1.7 -1.9 1.9]
sfp = linspace(0,sfp_max,201);
sfp = [sfp(2:end) 1];
sfp = unique(sfp);
sfp = sort(sfp);

% find nearby surfaces in Mac.s for sfp
Ifp = ones(size(sfp));
for k=1:length(sfp)
    [Y,J] = min(abs(Mac.s-sfp(k)));
    Ifp(k) = J;
end
Ifp = unique(Ifp);

[Y,J] = min(abs(Mac.s-sfp_max));
Ifp = [2:J];

% new (R,Z,Br,Bz) in Foureir space in
% geometric poloidal angle
MfpE = 100;
MfpP = 50;
Rfp = zeros(length(Ifp),MfpE+1);
Zfp = Rfp;
Brfp = zeros(length(Ifp),2*MfpP+1);
Bzfp = Brfp;

% find magnetic axis, and use it as the origin
% for defining the gepmetric theta angle (= tet)
R0 = R(1,1);
Z0 = Z(1,1);

% for each surface, compute (R,Z,Br,Bz) 
% in geometric angle Fourier decomposition
for kfp = 1:length(Ifp)
    J = Ifp(kfp);
    Rk = R(J,:);
    Zk = Z(J,:);
    Brk = Br(J,:);
    Bzk = Bz(J,:);

    % compute geometric angle
    % and define dtet, taking into account possible jumps
    tet = atan2(Zk-Z0,Rk-R0);
    dtet = diff(tet);
    II = find(abs(dtet)>pi);
    dtet(II) = dtet(II) - sign(dtet(II))*2*pi;
    tetm = tet(1:end-1) + dtet/2;

    % Fourier decompose (R,Z) in geometric angle
    % first compute all quantities at half-points in theta
    expt = exp(-i*tetm(:)*[0:MfpE]);
    Rk = (Rk(1:end-1)+Rk(2:end))/2;
    Zk = (Zk(1:end-1)+Zk(2:end))/2;
    Rfp(kfp,:) = Rk.*dtet*expt/2/pi;
    Zfp(kfp,:) = Zk.*dtet*expt/2/pi;

    % Fourier decompose (Br,Bz) in geometric angle
    % first compute all quantities at half-points in theta
    expt = exp(-i*tetm(:)*[-MfpP:MfpP]);
    Brk = (Brk(1:end-1)+Brk(2:end))/2;
    Bzk = (Bzk(1:end-1)+Bzk(2:end))/2;
    Brfp(kfp,:) = Brk.*dtet*expt/2/pi;
    Bzfp(kfp,:) = Bzk.*dtet*expt/2/pi;
end

% re-normalize (R,Z)
Rfp = Rfp*Mac.R0EXP;
Zfp = Zfp*Mac.R0EXP;

% save data to a file
sfpn = Mac.s(Ifp);
res = [sfpn(:) real(Rfp) imag(Rfp) real(Zfp) imag(Zfp)];
save FootPrintData_RZ res -ascii
res = [sfpn(:) real(Brfp) imag(Brfp) real(Bzfp) imag(Bzfp)];
save FootPrintData_B res -ascii

% plot data
kplot = 1;
if kplot > 0
   figure, plot(sfpn,real(Rfp))
   figure, plot(sfpn,imag(Rfp))
   figure, plot(sfpn,real(Zfp))
   figure, plot(sfpn,imag(Zfp))
   figure, plot(sfpn,real(Brfp))
   figure, plot(sfpn,imag(Brfp))
   figure, plot(sfpn,real(Bzfp))
   figure, plot(sfpn,imag(Bzfp))
end


