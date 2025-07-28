% Compare BnEDGE abd XnEDGE

load TMP_XnE
Tgx = Tg(:);
XnE = BnEDGE(:);
REx = Rs(:);

load TMP_BnE
Tgb = Tg(:);
BnE = (BnEDGE(:)-BnEDGE_VAC(:));
BnV = BnEDGE_VAC(:).*REx.^2;

res = [Tgx abs(BnE.*REx.^2./XnE) abs(BnV.*REx.^2./XnE)];

figure
plot(Tgx,abs(XnE)/max(abs(XnE)),'r-','LineWidth',2), hold on,
plot(Tgb,abs(BnE)/max(abs(BnE)),'b--','LineWidth',2), hold on,
plot(Tgb,abs(BnV)/max(abs(BnV)),'k-.','LineWidth',2), hold on,
