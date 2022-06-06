% 5MW Turbine Parameters

global  rho Rrad Ngear Jtot TSRopt Geff
global PratedM WratedRPM WratedRad Trated K 

TSRopt=7.55;
Ngear=97;
Rrad=63;
rho=1.225;
Geff=.944;
PratedM=5*10^6/Geff;
WratedRPM=12.1;
WratedRad=WratedRPM*pi/30;
Trated=PratedM/(WratedRad*Ngear);
K = 0.025576386;
Jtot=38759228; % From Fast.fsm file

%% Copy to Parameters Struct for Completeness

Parameters.Turbine.TSRopt=7.55;
Parameters.Turbine.Ngear=97;
Parameters.Turbine.Rrad=63;
Parameters.Turbine.rho=1.225;
Parameters.Turbine.Geff=.944;
Parameters.Turbine.PratedM=5*10^6/Geff;
Parameters.Turbine.WratedRPM=12.1;
Parameters.Turbine.WratedRad=WratedRPM*pi/30;
Parameters.Turbine.Trated=PratedM/(WratedRad*Ngear);
Parameters.Turbine.K = 0.025576386;
Parameters.Turbine.Jtot=38759228; % From Fast.fsm file

%% Yaw Controller Paramters

% Parameters.Turbine.Yaw = 0;     %for initialization


% global  WT_DATA WT_XDATA WT_YDATA
% load('5MW_LU_TSR2.mat')

% global  OP
% OP=csvread('OpPts_New_UW_CP_SE_OP.csv');
