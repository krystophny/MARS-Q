function MacReadGAMMA(filename)

global Mac

command = ['!cp ' filename ' GAMMA'];
eval(command), load GAMMA,

%Mac.CNORM = GAMMA(2);
Mac.CNORM = abs(GAMMA(3)-2.0e-4);
