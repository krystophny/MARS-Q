function MacGetJatGaussQuad(R, Z, Jr, Jz, Jphi)

global Mac

[xx,yy,ww] = MacGaussQuad2D(Mac.NGauss);

N = Mac.Ns1;

% reduce the number of radial points by factor of 2
if (Mac.ReduceNs) 
Nn = floor(N/2);
II = [2*[1:Nn]-1 N];
R = R(II,:);  Z = Z(II,:);
Jr = Jr(II,:); Jz = Jz(II,:); Jphi = Jphi(II,:); 
N = length(II);
end

R1 = R(1:N-1,1:end-1)*Mac.R0EXP;  Z1 = Z(1:N-1,1:end-1)*Mac.R0EXP;
R2 = R(2:N,1:end-1)*Mac.R0EXP;    Z2 = Z(2:N,1:end-1)*Mac.R0EXP;
R3 = R(2:N,2:end)*Mac.R0EXP;      Z3 = Z(2:N,2:end)*Mac.R0EXP;
R4 = R(1:N-1,2:end)*Mac.R0EXP;    Z4 = Z(1:N-1,2:end)*Mac.R0EXP;

R1 = R1(:); Z1 = Z1(:);
R2 = R2(:); Z2 = Z2(:);
R3 = R3(:); Z3 = Z3(:);
R4 = R4(:); Z4 = Z4(:);

ar = 0.25*(+ R1 + R2 + R3 + R4);  az = 0.25*(+ Z1 + Z2 + Z3 + Z4);
br = 0.25*(- R1 + R2 + R3 - R4);  bz = 0.25*(- Z1 + Z2 + Z3 - Z4);
cr = 0.25*(- R1 - R2 + R3 + R4);  cz = 0.25*(- Z1 - Z2 + Z3 + Z4);
dr = 0.25*(+ R1 - R2 + R3 - R4);  dz = 0.25*(+ Z1 - Z2 + Z3 - Z4);

RR = 0.25*R1*((1-xx').*(1-yy')) + ...
     0.25*R2*((1+xx').*(1-yy')) + ...
     0.25*R3*((1+xx').*(1+yy')) + ...
     0.25*R4*((1-xx').*(1+yy'));
ZZ = 0.25*Z1*((1-xx').*(1-yy')) + ...
     0.25*Z2*((1+xx').*(1-yy')) + ...
     0.25*Z3*((1+xx').*(1+yy')) + ...
     0.25*Z4*((1-xx').*(1+yy'));

weight = (br.*cz - cr.*bz)*ww' + (br.*dz - dr.*bz)*(xx'.*ww') + (dr.*cz - cr.*dz)*(yy'.*ww');

J1 = Jr(1:N-1,1:end-1); J2 = Jr(2:N,1:end-1); J3 = Jr(2:N,2:end); J4 = Jr(1:N-1,2:end);
JJr = 0.25*J1(:)*((1-xx').*(1-yy')) + ...
      0.25*J2(:)*  ((1+xx').*(1-yy')) + ...
      0.25*J3(:)*    ((1+xx').*(1+yy')) + ...
      0.25*J4(:)*  ((1-xx').*(1+yy'));

J1 = Jz(1:N-1,1:end-1); J2 = Jz(2:N,1:end-1); J3 = Jz(2:N,2:end); J4 = Jz(1:N-1,2:end);
JJz = 0.25*J1(:)*((1-xx').*(1-yy')) + ...
      0.25*J2(:)*  ((1+xx').*(1-yy')) + ...
      0.25*J3(:)*    ((1+xx').*(1+yy')) + ...
      0.25*J4(:)*  ((1-xx').*(1+yy'));

J1 = Jphi(1:N-1,1:end-1); J2 = Jphi(2:N,1:end-1); J3 = Jphi(2:N,2:end); J4 = Jphi(1:N-1,2:end);
JJphi = 0.25*J1(:)*((1-xx').*(1-yy')) + ...
        0.25*J2(:)*  ((1+xx').*(1-yy')) + ...
        0.25*J3(:)*    ((1+xx').*(1+yy')) + ...
        0.25*J4(:)*  ((1-xx').*(1+yy'));

FinalRes = [RR(:) ZZ(:) weight(:) real(JJr(:)) imag(JJr(:)) real(JJz(:)) imag(JJz(:)) real(JJphi(:)) imag(JJphi(:))];

save MacDataJp FinalRes -ascii

disp(['   Check total area = ' num2str(sum(weight(:)))])
disp(['   Check total Jr = ' num2str(sum(weight(:).*JJr(:)))])
disp(['   Check total Jz = ' num2str(sum(weight(:).*JJz(:)))])
disp(['   Check total Jphi = ' num2str(sum(weight(:).*JJphi(:)))])


