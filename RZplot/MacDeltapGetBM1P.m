function [BM1P,BM1A] = MacDeltapGetBM1P(BM1,BM2,BM3)

global Mac

KDELTAP = 0;

JJ = Mac.Iratsurf;
J  = JJ(Mac.Kratsurf)+1;

if KDELTAP==0
   BM1P = (BM1(J+1,:)-BM1(J,:))/(Mac.s(J+1)-Mac.s(J));
elseif KDELTAP==1
   BM1P = -i*(Mac.Mm'.*BM2(J+1,:) + Mac.n*BM3(J+1,:));
elseif KDELTAP==2
   II=[J+2:Mac.Ns1];
   xo = Mac.s(II)-Mac.s(J);
   yo = real(BM1(II,:)); 
   nu = (log(abs(yo(2,:)))-log(abs(yo(1,:))))./(log(xo(2))-log(xo(1))) - 1;
   logDelta = log(abs(yo(1,:))) - (1+nu)*log(xo(1));
   Delta = exp(logDelta);
   IK = find(yo(1,:)<0); Delta(IK)=-Delta(IK);
   BM1P = Delta;
elseif KDELTAP==5
   xs  = Mac.s(J);
   nu0 = 1.3381e-02;
   II=[J+1:Mac.Ns1];
   xo  = Mac.s(II)-xs; yo = real(BM1(II,:));
   BM1P = MacDeltapGetDs(xo,yo,nu0,2,[1:8],0);
end

BM1P = BM1P(:);
BM1A = BM1(JJ,:);
BM1A = transpose(BM1A);

