%plot plasma displacement vs coil phasing
%for D3D shot 161199

%input data from MARS-F
dataR = [
%                n=2                               n=4
%s     Re(Xi_U) Im(xi_U)  Re(Xi_L) Im(Xi_L)  Re(Xi_U) Im(xi_U)  Re(Xi_L) Im(Xi_L)
0.985  -0.1642  0.8706   -0.5362   1.0152    0.0821  -0.4353   -0.3055  -0.4860
0.995  -0.3256  0.8615   -0.7434   1.1507    0.1628  -0.4308   -0.3124  -0.6096
];
dataZ = [
%                n=2                               n=4
%s     Re(Xi_U) Im(xi_U)  Re(Xi_L) Im(Xi_L)  Re(Xi_U) Im(xi_U)  Re(Xi_L) Im(Xi_L)
0.985 -2.9191  -3.1236   -4.5383  -2.1950    1.4596   1.5618    2.0851	-1.4164
0.995 -3.6027  -5.0358   -8.3709  -1.9585    1.8014   2.5179    2.9408	-3.1351
];

%input data from Expt.
dataExp=[
%row1=x; row2=y_low; row3=y_high
29.411762505941752, 60.00002243939401, 89.41178494533574, 119.41177372563874, 149.999988780303, 179.41175128624482, 210.00001121969706, 240.58822627436126, 269.999988780303, 300.58820383496726, 330.58823749405826, 359.41172884685074
-0.004316546028779191, -0.004244601940711101, -0.0033812948389847957, -0.0034532416714392316, 0.0009352516471921191, 0.005467627653186957, 0.0061151079794816866, 0.0038848933927114852, 0.0005755394400107101, -0.0006474798688970058, -0.0023740995611223102, -0.0043884873724609344
-0.0027338117683037193, -0.002877697200053553, -0.0005755385252152623, 0.0005035980963289666, 0.003237410322030409, 0.006834532393844505, 0.007553956808207323, 0.00625899341123152, 0.0025179859076675908, 0.0012230225106917847, -0.0009352507323966714, -0.0030215826318033857
];

%set up phase parameters
N0 = 37; 
ND = 37;
%p0 = linspace(0,360,N0); %toroidal angle of measurement
p0 = 0; N0=1;
dp = linspace(0,360,ND); %coil phasing
thet = 18*pi/180;         %"tilting" angle 

%superposition and plotting
for k=1:size(dataR,1);
    XiU2 = (dataR(k,2)+dataR(k,3)*1i)*cos(thet) + (dataZ(k,2)+dataZ(k,3)*1i)*sin(thet);
    XiL2 = (dataR(k,4)+dataR(k,5)*1i)*cos(thet) + (dataZ(k,4)+dataZ(k,5)*1i)*sin(thet);
    XiU4 = (dataR(k,6)+dataR(k,7)*1i)*cos(thet) + (dataZ(k,6)+dataZ(k,7)*1i)*sin(thet);
    XiL4 = (dataR(k,8)+dataR(k,9)*1i)*cos(thet) + (dataZ(k,8)+dataZ(k,9)*1i)*sin(thet);

    hf = figure(k);
    for j=1:N0
        XiU = XiU2*exp(i*2*p0(j)*pi/180) + XiU4*exp(i*4*p0(j)*pi/180)*1;
        XiL = XiL2*exp(i*2*p0(j)*pi/180) + XiL4*exp(i*4*p0(j)*pi/180)*1;
        Xi  = 2*real( XiU + XiL*exp(i*dp*pi/180) );
        col = [(j-1)/N0 0 (N0-j)/N0];
        plot(dp,Xi,'LineWidth',2,'Color',col), hold on
    end
    for j=1:size(dataExp,2)
        plot([dataExp(1,j) dataExp(1,j)],[dataExp(2,j) dataExp(3,j)]*1000,'k-','LineWidth',3)
    end
    xlabel('\Delta\Phi [deg]','FontSize',18,'FontWeight','Bold')
    ylabel('{\xi}_R [mm]','FontSize',18,'FontWeight','Bold')
    ha = get(hf,'CurrentAxes'); set(ha,'FontSize',18,'FontWeight','Bold')
    a = axis; axis([0 360 a(3) a(4)]) 
end

 
