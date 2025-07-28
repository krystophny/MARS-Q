% Set basic parameters for MacRead* procedures
% kscale = 0: direct MARS-Q output
%          1: in physical units
%          2: normalised profiles to 1

%kdev    = 67; 
             %1=MAST25075; 2=ANATOR; 3=D3D138593; 4=JET77329; 5=D3D147626; 6=DELTAPANA
             %7=MAST28118; 8=JET40542; 9=ITER; 10=MAST27873; 11=MAST25056; 12=NTV_CHECK
	     %13=ITER_ABT4ZL, 14=MAST30090; 15=D3D147131; 16=EXAMPLE; 17=MAST31128.3500_n=4; 
             %18=MAST31128.3500_n=2; %19=MAST31131.6400_n=2; %20=MAST30684.6950_n=2; 
             %21=AUG31021; 22=F4E_15MAFT; 23=F4E_9MA; 24=D3D161243; 25=D3D157376
	     %26=D3D158104; 27=IK; 28=D3D156908; 29=156898; 30=158115; 
             %31=ITER_15MA/5.3T/Q10_ISF; 32=MAST22264@200; 33=MAST30128@250;
             %34=VDE_TEST; 35=MAST28912@290; 36=MAST29222@290;
             %37=Dong; 38=AUG33353; 39=D3D174045(avalanche); 
             %40=ITER 5MA/1.8T_ISF; 41=ITER 7.5MA/2.65T_ISF
             %42=ITER 10MA/5.3T_ISF; 43=ITER 5MA/1.8T_ISF_neg8; 
             %44=D3D177322(avalanche); 45=D3D157376(n=6)
	     %46=D3D184003.2750(NTV); 47=D3D184003.3880(NTV);
             %48=ITER 7.5MA/5.3T_ISF; 49=ITER 15MA/5.3T/Q5_ISF;
             %55=MAST-U45272; 56=STEP_CASE1; 57=F4E_9MA
             %58=ITER_5MA_GAP3 (=Standard)
             %59=ITER_5MA_GAP4 (=Clearance)
             %60=ITER_5MA_GAP5 (=Outergap)
	     %61=STEP_CASE2; 62=STEP_CASE3; 63=STEP_CASE4; 
             %64=ITER_10MA_RIPPLE; 65=KSTAR30306.3900
	     %66=STEP_CASE5; 67=EAST52340
             %68=D3D163519.1750; 69=D3D163519.2650; 70=D3D163520.1750; 71=D3D163520.2650; 			       

if kdev==1
   SDIR    = '../MAST/MAST25075/';
   ASPECT  = 1.6789; %aspect ratio
   mm      = [-29:29];
   n       = 3;      %toroidal mode number
   B0      = 0.4299; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8938; %major radius in [m]
   TAUA    = 7.9461e-7;  %Alfven time [s]
elseif kdev==2
   SDIR    = '~/Work/ANATOR04/Data/';
   B0      = 2.0;    %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 3.0;    %major radius in [m]
   ASPECT  = 3.0;    %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   TAUA    = 5.2198e-07;  %Alfven time [s]
elseif kdev==3
   SDIR    = '/.automount/funsrv1/root/home/yliu/D3D138593/Work/';
   ASPECT  = 2.8896; %aspect ratio
   mm      = [-29:29];
   n       = 3;      %toroidal mode number
   B0      = 1.8990; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.0522e-7;  %Alfven time [s]
elseif kdev==4
   SDIR    = '/.automount/funsrv1/root/home/yliu/JET77329/';
   ASPECT  = 3.2031; %aspect ratio
   mm      = [-29:29];
   n       = 1;      %toroidal mode number
   B0      = 1.8653; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 2.9600; %major radius in [m]
   TAUA    = 5.730454441699772E-007;  %Alfven time [s]
   kscale  = 1;      %=1: plot in physical units
   SDIR    = '/.automount/funsrv4.ccfe.ac.uk/root/home1/yliu/NTV_CHECK_CASE/Liu/';
elseif kdev==5
   SDIR    = '/.automount/funsrv1/root/home/yliu/D3D147626/';
   ASPECT  = 2.8609; %aspect ratio
   mm      = [-29:29];
   n       = 1;      %toroidal mode number
   B0      = 1.7550; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.8821e-7;  %Alfven time [s]
elseif kdev==6
   SDIR    = '/.automount/funsrv1/root/home/yliu/DeltapAna/Work/';
   ASPECT  = 10.0; %aspect ratio
   mm      = [-5:5];
   n       = 1;      %toroidal mode number
   B0      = 2.5000; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 10.0; %major radius in [m]
   TAUA    = 1.0e-6;  %Alfven time [s]
elseif kdev==7
   SDIR    = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';
   ASPECT  = 1.4970; %aspect ratio
   mm      = [-29:29];
   n       = 3;      %toroidal mode number
   B0      = 0.4634; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8780; %major radius in [m]
   TAUA    = 6.5293e-7;  %Alfven time [s]
elseif kdev==8
   SDIR    = '/.automount/funsrv4.ccfe.ac.uk/root/home1/yliu/JET40542/Temp05/';
   ASPECT  = 3.0794; %aspect ratio
   mm      = [-29:29];
   n       = 1;      %toroidal mode number
   B0      = 1.2180; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 2.8960; %major radius in [m]
   TAUA    = 8.5710e-7;  %Alfven time [s]
elseif kdev==9
   SDIR    = '/.automount/funsrv1/root/home/yliu/Scen2_V02/';
   ASPECT  = 3.0947; %aspect ratio
   mm      = [-29:29];
   n       = 4;      %toroidal mode number
   B0      = 5.3000; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2000; %major radius in [m]
   TAUA    = 9.0100e-7;  %Alfven time [s]
elseif kdev==10
   SDIR    = '/.automount/funsrv1/root/home/yliu/MAST_ADD/';
   ASPECT  = 1.4757; %aspect ratio
   mm      = [-35:35];
   n       = 3;      %toroidal mode number
   B0      = 0.4831; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8419; %major radius in [m]
   TAUA    = 6.3183e-7;  %Alfven time [s]
elseif kdev==11
   %SDIR    = '../MAST/MAST25056/';   
   SDIR    = '~/Work/MAST/DataNQ/';
   ASPECT  = 1.4298; %aspect ratio
   mm      = [-35:35];
   n       = 3;      %toroidal mode number
   B0      = 0.4814; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8473; %major radius in [m]
   TAUA    = 6.2697e-7;  %Alfven time [s]
elseif kdev==12
   SDIR    = '/home/liuy/Work/IK/';
   ASPECT  = 10.0; %aspect ratio
   mm      = [-5:25];
   n       = 1;      %toroidal mode number
   B0      = 1.0; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 2.0; %major radius in [m]
   TAUA    = 2.8995E-07;  %Alfven time [s]
elseif kdev==13
   SDIR    = '/home/yliu/ITER/Scen2_V03/ResDATA_v1.2/n=4D/NCASE6KPB_PLS_v1.2/';
   ASPECT  = 3.0908; %aspect ratio
   mm         = [-35:35];
   n             = 4;      %toroidal mode number
   B0          = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0          = 6.2; %major radius in [m]
   TAUA     = 8.5642e-07;  %Alfven time [s]
elseif kdev==14
   SDIR    = '/home/yliu/MAST_ADD/';
   ASPECT  = 1.7150; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 0.4311; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.9457; %major radius in [m]
   TAUA    = 8.4482e-7;  %Alfven time [s]
elseif kdev==15
   SDIR    = '/home/yliu/D3D147131/Work/';
   ASPECT  = 2.7935; %aspect ratio
   mm      = [-29:29];
   n       = 1;      %toroidal mode number
   B0      = 1.7228; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.6158e-7;  %Alfven time [s]
elseif kdev==16
   SDIR    = '/home/yliu/Temp/';
   ASPECT  = 3.0000; %aspect ratio
   mm      = [-29:29];
   n       = 1;      %toroidal mode number
   B0      = 1.5000; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 3.0000; %major radius in [m]
   TAUA    = 5.0000e-7;  %Alfven time [s]
elseif kdev==17
   SDIR    = '/home/yliu/ASDEX/Work/';
   ASPECT  = 3.3208; %aspect ratio
   mm      = [-29:29];
   n       = 4;      %toroidal mode number
   B0      = 1.7897; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.7060; %major radius in [m]
   TAUA    = 3.7334e-7;  %Alfven time [s]
elseif kdev==18
   SDIR    = '/home/yliu/Temp/Work/';
   ASPECT  = 3.3208; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 1.7897; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.7060; %major radius in [m]
   TAUA    = 3.7334e-7;  %Alfven time [s]
elseif kdev==19
   SDIR    = '/home/yliu/ASDEX/Work/';
   ASPECT  = 3.3747; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 1.7390; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.7357; %major radius in [m]
   TAUA    = 4.2795e-7;  %Alfven time [s]
elseif kdev==20
   SDIR    = '/home/yliu/ASDEX/Work/';
   ASPECT  = 3.3764; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 2.3886; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.7117; %major radius in [m]
   TAUA    = 3.5943e-07;  %Alfven time [s]
elseif kdev==21
   SDIR    = '/home/yliu/Temp/';
   ASPECT  = 3.4630; %aspect ratio
   mm      = [-29:29];
   n       = 1;      %toroidal mode number
   B0      = 2.4741; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.7141; %major radius in [m]
   TAUA    = 4.2212e-07;  %Alfven time [s]
elseif kdev==22
   SDIR    = '/home/yliu/F4E_2014/DataRes/10470/n=1/';
   %SDIR    = '/home/yliu/F4E_2014/DataRes/10470/n=1/15MAFTFITBM/';
   ASPECT  = 3.1166; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 5.3000; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2000; %major radius in [m]
   TAUA    = 8.9397e-07;  %Alfven time [s]
elseif kdev==23
   %SDIR    = '/home/yliu/F4E_2014/DataRes/10100/n=1/';
   SDIR    = '../F4E_2014/DataRes/10100/n=1/9MAFTFITBMELM/';
   ASPECT  = 3.1034; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 5.3000; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2000; %major radius in [m]
   TAUA    = 7.7025e-07;  %Alfven time [s]
elseif kdev==24
   %SDIR    = '../../YQLiu/Collaborators/OkabayashiMichio/ModeLockingFeedback/Work/';
   SDIR    = '../Temp/';
   ASPECT  = 2.7651; %aspect ratio
   mm      = [-26:26];
   n       = 1;      %toroidal mode number
   B0      = 1.9967; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 1.1226e-06;  %Alfven time [s]
elseif kdev==25
   SDIR    = '~/Work/D3D157376/Data/';
   %SDIR    = '/cscratch/liuy/WorkD3D157376e/';
   ASPECT  = 2.9219; %aspect ratio
   mm      = [-33:33];
   n       = 3;      %toroidal mode number
   B0      = 1.9104; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.6541e-07;  %Alfven time [s]
elseif kdev==26
   SDIR    = '/cscratch/liuy/D3D158104/Data/';
   ASPECT  = 2.8141; %aspect ratio
   mm      = [-33:33];
   n       = 2;      %toroidal mode number
   B0      = 1.8534; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 3.8355e-07;  %Alfven time [s]
elseif kdev==27
   SDIR    = '~/IRIS/WorkIK/';
   ASPECT  = 10; %aspect ratio
   mm      = [-5:11];
   n       = 1;      %toroidal mode number
   B0      = 1.5; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 10; %major radius in [m]
   TAUA    = 1.0e-06;  %Alfven time [s]
elseif kdev==28
   SDIR    = '~/Work/D3D156908/Data/';
   ASPECT  = 2.9118; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 1.7690; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.6026e-07;  %Alfven time [s]
elseif kdev==29
   SDIR    = '~/Work/D3D156898/Data/';
   ASPECT  = 2.9187; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 1.7673; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.5848e-07;  %Alfven time [s]
elseif kdev==30
   SDIR    = '/cscratch/liuy/Temp/Data/';
   ASPECT  = 2.8113; %aspect ratio
   mm      = [-29:29];
   n       = 2;      %toroidal mode number
   B0      = 1.8951; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.1490e-07;  %Alfven time [s]
elseif kdev==31
   %SDIR    = '/cscratch/liuy/WorkITERa/';
   SDIR    = '/home/liuy/Work/ITER/Data_QL/15MAQ10/';
   ASPECT  = 3.0908; %aspect ratio
   mm      = [-35:35];
   n       = 3;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 7.7812e-07;  %Alfven time [s]
elseif kdev==32
   %SDIR    = '/cscratch/liuy/WorkMASTa4/';
   SDIR    = '/home/liuy/Work/MAST/DataBQ/';
   ASPECT  = 1.4226; %aspect ratio
   mm      = [-33:33];
   n       = 2;      %toroidal mode number
   B0      = 0.4642; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8779; %major radius in [m]
   TAUA    = 7.2120e-07;  %Alfven time [s]
elseif kdev==33
   %SDIR    = '/cscratch/liuy/WorkMASTb4/';
   SDIR    = '/home/liuy/Work/MAST/DataBQ/';
   ASPECT  = 1.5190; %aspect ratio
   mm      = [-33:33];
   n       = 2;      %toroidal mode number
   B0      = 0.4505; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8519; %major radius in [m]
   TAUA    = 6.1841e-07;  %Alfven time [s]
elseif kdev==34
   SDIR    = '/cscratch/liuy/WorkVDE/';
   ASPECT  = 3; %aspect ratio
   mm      = [-1:1];
   n       = 0;      %toroidal mode number
   B0      = 2.0; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.7; %major radius in [m]
   TAUA    = 5.0e-07;  %Alfven time [s]
elseif kdev==35
   %SDIR    = '/cscratch/liuy/WorkMASTa/';
   SDIR    = '/home/liuy/Work/MAST/DataNQ/';
   ASPECT  = 1.4045; %aspect ratio
   mm      = [-33:33];
   n       = 2;      %toroidal mode number
   B0      = 0.4322; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8435; %major radius in [m]
   TAUA    = 9.1141e-07;  %Alfven time [s]
elseif kdev==36
   %SDIR    = '/cscratch/liuy/WorkMASTb2/';
   SDIR    = '/home/liuy/Work/MAST/DataNQ/';
   ASPECT  = 1.3815; %aspect ratio
   mm      = [-33:33];
   n       = 2;      %toroidal mode number
   B0      = 0.5038; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8088; %major radius in [m]
   TAUA    = 6.3546e-07;  %Alfven time [s]
elseif kdev==37
   SDIR    = '/home/liuy/Work/Dong/';
   ASPECT  = 1.3815; %aspect ratio
   mm      = [-9:49];
   n       = 1;      %toroidal mode number
   B0      = 0.5038; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8088; %major radius in [m]
   TAUA    = 6.3546e-07;  %Alfven time [s]
elseif kdev==38
   SDIR    = '/home/liuy/Work/AUG33353/DataQL/';
   ASPECT  = 3.2652; %aspect ratio
   mm      = [-33:33];
   n       = 2;      %toroidal mode number
   B0      = 1.7573; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6991; %major radius in [m]
   TAUA    = 4.6956e-07;  %Alfven time [s]
elseif kdev==39
   %SDIR    = '/home/liuy/Work/RE_AVALANCHE/TimeNew/';
   %SDIR    = '/cscratch/liuy/Work_VS/';
   SDIR    = '/cscratch/liuy/WorkTEMP_03/';
   %SDIR     = '/home/liuy/Work/RE_AVALANCHE/TimeVDE/';
   ASPECT  = 2.5606; %aspect ratio
   mm      = [-15:15];
   n       = 0;      %toroidal mode number
   B0      = 2.1330; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6995; %major radius in [m]
   TAUA    = 2.0008e-07;  %Alfven time [s]
elseif kdev==40
   %SDIR    = '/home/liuy/Work/ITER/Data_QL/5MA/';
   %SDIR    = '/cscratch/liuy/WorkTEMP/';
   SDIR     = '/home/liuy/Work/ITER_NTV/A/rot1/n4/';
   ASPECT  = 3.0913; %aspect ratio
   mm      = [-33:33];
   n       = 4;      %toroidal mode number
   B0      = 1.8; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 6.3284e-07;  %Alfven time [s]
elseif kdev==41
   %SDIR    = '/home/liuy/Work/ITER/Data_QL/7d5MAHalfB/';
   SDIR     = '/home/liuy/Work/ITER_NTV/B/rot3/n3/';
   ASPECT  = 3.0921; %aspect ratio
   mm      = [-33:33];
   n       = 3;      %toroidal mode number
   B0      = 2.65; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 9.8103e-07;  %Alfven time [s]
elseif kdev==42
   SDIR    = '/home/liuy/Work/ITER/Data_QL/10MA/';
   ASPECT  = 3.0921; %aspect ratio
   mm      = [-49:49];
   n       = 3;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 6.5542e-07;  %Alfven time [s]
elseif kdev==43
   SDIR    = '/home/liuy/Work/ITER/Data_QL/5MA_neg8/';
   ASPECT  = 3.0913; %aspect ratio
   mm      = [-42:42];
   n       = 3;      %toroidal mode number
   B0      = 1.8; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 6.3284e-07;  %Alfven time [s]
elseif kdev==44
   SDIR    = '/home/liuy/Work/RE_AVALANCHE/Time2/';
   ASPECT  = 2.7756; %aspect ratio
   mm      = [-25:25];
   n       = 0;      %toroidal mode number
   B0      = 2.1025; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6995; %major radius in [m]
   TAUA    = 4.4366e-07;  %Alfven time [s]
elseif kdev==45
   SDIR    = '~/Work/D3D157376/Data_Mn6/';
   %SDIR    = '/cscratch/liuy/WorkD3D157376e/';
   ASPECT  = 2.9219; %aspect ratio
   mm      = [-35:35];
   n       = 6;      %toroidal mode number
   B0      = 1.9104; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.6541e-07;  %Alfven time [s]
elseif kdev==46
   SDIR    = '/cscratch/liuy/WorkD3D184003/';
   ASPECT  = 2.7682; %aspect ratio
   mm      = [-33:33];
   n       = 3;      %toroidal mode number
   B0      = 1.9169; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 3.1662e-07;  %Alfven time [s]
elseif kdev==47
   %SDIR    = '/cscratch/liuy/WorkD3D184003/';
   %SDIR    = '/cscratch/liuy/WorkTEMP_08/';
   SDIR     = '/home/liuy/Work/D3D184003/Data/';
   ASPECT  = 2.7610; %aspect ratio
   mm      = [-33:33];
   n       = 3;      %toroidal mode number
   B0      = 1.9098; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 3.1397e-07;  %Alfven time [s]
   psipp   = [0.9361579      0.9469605     0.9675397     0.9765310      0.9912874]; %2740
   %psipp   = [0.9370174      0.9455787     0.9686599     0.9767946      1]; %3880
elseif kdev==48
   SDIR    = '/home/liuy/Work/ITER/Data_QL/New/case4_7d5MA5d3T_qa8/';
   ASPECT  = 3; %aspect ratio
   mm      = [-33:33];
   n       = 3;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 2.8897e-07;  %Alfven time [s]
elseif kdev==49
   SDIR    = '/home/liuy/Work/ITER/Data_QL/New/case6_15MA5d3T_qa3d9/';
   ASPECT  = 3; %aspect ratio
   mm      = [-33:33];
   n       = 3;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 2.8897e-07;  %Alfven time [s]
elseif kdev==50
   SDIR    = '/cscratch/liuy/WorkIR_D/AT/';
   ASPECT  = 3.0913; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 1.8; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 8.7744e-07;  %Alfven time [s]
elseif kdev==51
   SDIR    = '/cscratch/liuy/WorkIR_D/BT/';
   ASPECT  = 3.0921; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 2.65; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 8.8677e-07;  %Alfven time [s]
elseif kdev==52
   SDIR    = '/cscratch/liuy/WorkIR_D/CT/';
   ASPECT  = 3.0935; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 4.1391e-07;  %Alfven time [s]
elseif kdev==53
   SDIR    = '/cscratch/liuy/D3D136642_VDE/';
   ASPECT  = 3.1219; %aspect ratio
   mm      = [-23:23];
   n       = 0;      %toroidal mode number
   B0      = 1.7; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.5227e-07;  %Alfven time [s]
elseif kdev==54
   SDIR    = '/cscratch/liuy/D3D088806_VDE/';
   ASPECT  = 2.8524; %aspect ratio
   mm      = [-23:23];
   n       = 0;      %toroidal mode number
   B0      = 2.0874; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.4961e-07;  %Alfven time [s]
elseif kdev==55
   SDIR    = '/home/liuy/Work/MAST-U45272/Data/';
   ASPECT  = 1.5643; %aspect ratio
   mm      = [-9:39];
   n       = 1;      %toroidal mode number
   B0      = 0.5885; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 0.8000; %major radius in [m]
   TAUA    = 6.3274e-07;  %Alfven time [s]
elseif kdev==56
   SDIR    = '/home/liuy/Work/STEP-P/Case1/DataN/';
   %SDIR    = '/cscratch/liuy/STEP-P/Case1_n2A/'
   ASPECT  = 1.8031; %aspect ratio
   mm      = [-31:31];
   n       = 2;      %toroidal mode number
   B0      = 3.2004; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 3.5995; %major radius in [m]
   TAUA    = 1.1433e-06;  %Alfven time [s]

   if n==10
      copyfile([SDIR 'PROFEQ_n1'],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'qlin_rot_5e-3_fld/TIMEEVOL_n1_THR1e-8_53'],[SDIR 'TIMEEVOL.OUT'],'f')
   elseif n==2
      copyfile([SDIR 'PROFEQ_n2'],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'TIMEEVOL_n2A_01'],[SDIR 'TIMEEVOL.OUT'],'f')
    end
elseif kdev==57
   SDIR    = '/home/liuy/Work/ITER_RIPPLE/Data/9MA/RES_ALL_EP_NDB0/';
   ASPECT  = 3.1022; %aspect ratio
   mm      = [-80:80];
   n       = 18;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 7.7025e-07;  %Alfven time [s]
elseif kdev==58
   SDIR    = '/home/liuy/Work/ITER_RIPPLE/Data/5MA_GAP3/RES_ALL_NDB0/';
   ASPECT  = 3.1163; %aspect ratio
   mm      = [-80:80];
   n       = 18;      %toroidal mode number
   B0      = 2.65; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 6.7469e-07;  %Alfven time [s]
elseif kdev==59
   SDIR    = '/home/liuy/Work/ITER_RIPPLE/Data/5MA_GAP4/RES_ALL_NDB0/';
   ASPECT  = 3.2644; %aspect ratio
   mm      = [-80:80];
   n       = 18;      %toroidal mode number
   B0      = 2.65; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 6.7341e-07;  %Alfven time [s]
elseif kdev==60
   SDIR    = '/home/liuy/Work/ITER_RIPPLE/Data/5MA_GAP5/RES_ALL_NDB0/';
   ASPECT  = 3.2743; %aspect ratio
   mm      = [-80:80];
   n       = 18;      %toroidal mode number
   B0      = 2.65; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 7.8767e-07;  %Alfven time [s]
elseif kdev==61
   SDIR    = '/home/liuy/Work/STEP-P/Case2/DataN/';
   %SDIR     = '/home/liuy/Work/STEP-P/Case2_q10.05_n1/';
   %SDIR    = '/cscratch/liuy/STEP-P/Case2_n2A/'
   ASPECT  = 1.8031; %aspect ratio
   mm      = [-31:31];
   n       = 2;      %toroidal mode number
   B0      = 3.2004; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 3.5995; %major radius in [m]
   TAUA    = 1.1433e-06;  %Alfven time [s]

   if 1==1
      copyfile([SDIR 'PROFEQ_n' int2str(n)],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'TIMEEVOL_n' int2str(n) 'B_04'],[SDIR 'TIMEEVOL.OUT'],'f')
    end
elseif kdev==62
   SDIR    = '/home/liuy/Work/STEP-P/Case3/DataN_q95_11.9/';
   %SDIR     = '/home/liuy/Work/STEP-P/Case3_q11.4_n2/';
   %SDIR    = '/cscratch/liuy/STEP-P/Case3_n2A/'
   ASPECT  = 1.8033; %aspect ratio
   mm      = [-31:31];
   n       = 2;      %toroidal mode number
   B0      = 3.1926; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 3.6083; %major radius in [m]
   TAUA    = 1.1477e-06;  %Alfven time [s]

   if 1==1
      copyfile([SDIR 'PROFEQ_n' int2str(n)],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'TIMEEVOL_n' int2str(n) 'A_01'],[SDIR 'TIMEEVOL.OUT'],'f')
    end
elseif kdev==63
   SDIR    = '/home/liuy/Work/STEP-P/Case4/DataN_q95_9.8/';
   %SDIR    = '/cscratch/liuy/STEP-P/Case4_n2A/'
   ASPECT  = 1.7965; %aspect ratio
   mm      = [-31:31];
   n       = 2;      %toroidal mode number
   B0      = 3.1968; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 3.6036; %major radius in [m]
   TAUA    = 9.9311e-07;  %Alfven time [s]

   if n==10
      copyfile([SDIR 'PROFEQ_n1'],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'TIMEEVOL_n1_04'],[SDIR 'TIMEEVOL.OUT'],'f')
   elseif n==2
      copyfile([SDIR 'PROFEQ_n2'],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'TIMEEVOL_n2B_01'],[SDIR 'TIMEEVOL.OUT'],'f')
    end
elseif kdev==64
   SDIR    = '/home/liuy/Work/ITER_RIPPLE/Data/10MA/RES_ALL_EP_NDB0/';
   ASPECT  = 3.0921; %aspect ratio
   mm      = [-80:80];
   n       = 18;      %toroidal mode number
   B0      = 5.3; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 6.2; %major radius in [m]
   TAUA    = 6.5542e-07;  %Alfven time [s]
elseif kdev==65
   SDIR    = '/home/liuy/Work/KSTAR/Data/30306.3900/';
   ASPECT  = 3.8399; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 3.0; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 4.3; %major radius in [m]
   TAUA    = 4.2824e-07;  %Alfven time [s]
elseif kdev==66
   SDIR    = '/home/liuy/Work/STEP-P/Case5/DataN/';
   %SDIR    = '/cscratch/liuy/STEP-P/Case5_q95_10.969_n2/';
   %SDIR    = '/cscratch/liuy/STEP-P/Case5_q95_11.1_n2/'
   %SDIR    = '/cscratch/liuy/STEP-P/Case5_q95_10.969_n2/'
   %SDIR    = '/cscratch/liuy/STEP-P/Case5_n3A/'
   ASPECT  = 1.8000; %aspect ratio
   mm      = [-31:31];
   n       = 3;      %toroidal mode number
   B0      = 3.0; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 4.3; %major radius in [m]
   TAUA    = 1.1887e-06;  %Alfven time [s]

   if n==3 | n==4
     mm = [-39:39];
   end
   
   if n>1
      copyfile([SDIR 'PROFEQ_n' int2str(n)],[SDIR 'PROFEQ.OUT'],'f')
      copyfile([SDIR 'TIMEEVOL_n' int2str(n) 'B_06'],[SDIR 'TIMEEVOL.OUT'],'f')
    end
elseif kdev==67
   SDIR    = '/cscratch/liuy/ZhangHY/EAST52340_3150/';
   ASPECT  = 4.2294; %aspect ratio
   mm      = [-33:33];
   n       = 1;      %toroidal mode number
   B0      = 2.3010; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.8; %major radius in [m]
   TAUA    = 4.2824e-07;  %Alfven time [s]
elseif kdev==68
   SDIR    = '/home/liuy/Work/DIII-D/D3D_QH_CETOP/Data/163519.1750/';
   ASPECT  = 2.8495; %aspect ratio
   mm      = [-9:33];
   n       = 1;      %toroidal mode number
   B0      = 2.0083; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 3.3210e-07;  %Alfven time [s]
   
   SDIR    = [SDIR 'n' int2str(n) '/'];
   copyfile([SDIR 'TIMEEVOL_QL_05b'],[SDIR 'TIMEEVOL.OUT'],'f')
elseif kdev==69
   SDIR    = '/home/liuy/Work/DIII-D/D3D_QH_CETOP/Data/163519.2650/';
   ASPECT  = 2.8719; %aspect ratio
   mm      = [-9:33];
   n       = 2;      %toroidal mode number
   B0      = 2.0089; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 3.9370e-07;  %Alfven time [s]
   
   SDIR    = [SDIR 'n' int2str(n) '/'];
   copyfile([SDIR 'TIMEEVOL_QL_02b'],[SDIR 'TIMEEVOL.OUT'],'f')
elseif kdev==70
   SDIR    = '/home/liuy/Work/DIII-D/D3D_QH_CETOP/Data/163520.1750/';
   ASPECT  = 2.8848; %aspect ratio
   mm      = [-9:33];
   n       = 2;      %toroidal mode number
   B0      = 2.0112; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 3.3785e-07;  %Alfven time [s]
   
   SDIR    = [SDIR 'n' int2str(n) '/'];
   copyfile([SDIR 'TIMEEVOL_QL_03'],[SDIR 'TIMEEVOL.OUT'],'f')
elseif kdev==71
   SDIR    = '/home/liuy/Work/DIII-D/D3D_QH_CETOP/Data/163520.2650/';
   ASPECT  = 2.8585; %aspect ratio
   mm      = [-9:33];
   n       = 2;      %toroidal mode number
   B0      = 2.0126; %toroidal on axis vacuum magnetic field in [Tesla]
   R0      = 1.6955; %major radius in [m]
   TAUA    = 4.0317e-07;  %Alfven time [s]
   
   SDIR    = [SDIR 'n' int2str(n) '/'];
   copyfile([SDIR 'TIMEEVOL_QL_03b'],[SDIR 'TIMEEVOL.OUT'],'f')
end

