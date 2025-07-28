function MacInputPROFEQ(filename)

global Mac

data = load(filename);

Mac.q    = data(:,2);
Mac.P    = data(:,4);
Mac.F    = data(:,13);
Mac.dpsi = data(:,12);
Mac.omega= data(:,6);
Mac.den  = data(:,5);
