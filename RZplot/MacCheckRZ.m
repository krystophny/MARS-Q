function MacCheckRZ(filename)

global Mac

command = ['!cp ' filename ' RMZM_R'];
eval(command), load RMZM_R

MM = RMZM_R(1,1); NN=RMZM_R(1,2);
s = RMZM_R(2:NN+1,1);
chi = RMZM_R(NN+2:NN+MM+1,1);
R = RMZM_R(NN+MM+2:end,1); Z=RMZM_R(NN+MM+2:end,2);

R = reshape(R,MM,NN);
Z = reshape(Z,MM,NN);

% plot RZ coordinates
if Mac.plot
   plot(R*Mac.R0EXP,Z*Mac.R0EXP,'y-'), hold on,
   plot(R(:,end)*Mac.R0EXP,Z(:,end)*Mac.R0EXP,'r-'), hold on
end

