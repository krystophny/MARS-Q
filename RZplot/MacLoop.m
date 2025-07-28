global flux_mars
global KFEED

fluxR =[];
fluxI =[];

KK = [6];
for k=1:length(KK)
  KFEED=KK(k);
  MacMain
  fluxR = [fluxR; flux_mars(1,:)];  
  fluxI = [fluxI; flux_mars(2,:)];
end
  
