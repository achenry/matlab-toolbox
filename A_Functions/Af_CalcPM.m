function [PM,DELdat,PM_Desc] = Af_CalcPM(OutData,OutList,tout,varargin)
% J.Aho 9/22/11- Some code from perf99bynumlateY by F.Dunne
% deleted unused PM calcs
% edited so side-to-side IMU isn't calc'd till 300s
% PM contains up to 50 performance measures, listed in measurenumbers.xls

% This code will calculate the DELS for various components
% DEL Equivalant saved in PM
% DEL rainflow count vectors stored in cell array DELdat

%% Handle Varargin
if ~isempty(varargin)
    ts = varargin{1}(1);
    tend = varargin{1}(2);
    
else
    ts=10;     %Time to start calculating DELS in seconds
    tend = tout(end);
end

%% User Inputs
ts2=30;   %Time to start calculating Nacelle Accelerations
WhC=10;   %Wohler curve exponent for material (slope of log(S) vs log(N))
WhS=4;    %Eg. C=composite: 10-12, S=steel: 3
WhShft=5;


%%

PM=zeros(50,1);
TS2 = round(1+ts2/(tout(2)-tout(1)));        % # steps equiv to 300 sec

dt=tout(2)-tout(1);
TS=round(ts/dt);
nend=round(tend/dt);
tout_trim=tout(TS:nend,:);                % ignore the initial ts seconds of data  


DELdat=cell(12); 
PM_Desc_Overview{1}=['DEL calculated from time ',num2str(ts),' to time ' num2str(nend/dt)];


PM_Desc{1,1}='RootMyb1';
[PM(1), DELdat{1}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch(PM_Desc{1,1},OutList)) , WhC);
PM_Desc{2,1}='RootMyb2';
[PM(2), DELdat{2}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch(PM_Desc{2,1},OutList)) , WhC);
PM_Desc{3,1}='RootMyb3';
[PM(3), DELdat{3}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch(PM_Desc{3,1},OutList)) , WhC);
%PM_Desc{4,1}='Spn2MLyb1';
%[PM(4), DELdat{4}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('Spn2MLyb1',OutList)) , WhC);
%PM_Desc{5,1}='Spn2MLyb2';
%[PM(5), DELdat{5}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('Spn2MLyb2',OutList)) , WhC);
%PM_Desc{6,1}='Spn2MLyb3';
%[PM(6), DELdat{6}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('Spn2MLyb3',OutList)) , WhC);
%PM_Desc{7,1}='YawBrMyp';
%[PM(7), DELdat{7}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('YawBrMyp',OutList)) , WhS);
PM_Desc{8,1}='TwrBsMxt';
[PM(8), DELdat{8}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('TwrBsMxt',OutList)) , WhS);
PM_Desc{8,1}='TwrBsMyt';
[PM(8), DELdat{8}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('TwrBsMyt',OutList)) , WhS);
PM_Desc{9,1}='LSShftMxa';
[PM(9), DELdat{9}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('LSShftMxa',OutList)) , WhShft);

[blip1, temp11] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('RootMyc1',OutList)),WhC); 
[blip2, temp12] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('RootMyc2',OutList)),WhC); 
[blip3, temp13] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('RootMyc3',OutList)),WhC); 


[blop1, tempo11] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('RootMxc1',OutList)),WhC);
[blop2, tempo12] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('RootMxc2',OutList)),WhC);
[blop3, tempo13] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('RootMxc3',OutList)),WhC); 

DELdat{11}=[temp11,tempo11];
DELdat{12}=[temp12,tempo12];
DELdat{13}=[temp13,tempo13];

bl1=sqrt((blip1^2+blop1^2)/2);
bl2=sqrt((blip2^2+blop2^2)/2);
bl3=sqrt((blip3^2+blop3^2)/2);

PM_Desc{10,1}='Combined IP/OP BRBM';
PM(10)=(bl1+bl2+bl3)/3;            %combined in-plane and oop blade root loads

%XX Change this number later?
PM_Desc{14,1}='TwrBsMxt';
[PM(14), DELdat{14}] = mcdel(tout_trim, OutData.signals.values(TS:nend,strmatch('TwrBsMxt',OutList)) , WhS);


%power performance
PM_Desc{11,1}='Mean GenPwr';
PM(11) = mean(OutData.signals.values(TS:nend,strmatch('GenPwr',OutList)));   %GenPwr is in kW

%rotor speed
PM_Desc{17,1}='Max RotSpeed';
PM(17) = max(OutData.signals.values(TS:nend,strmatch('RotSpeed',OutList)));

%actuator effort
PM_Desc{22,1}='Mean RMS Pitch Rate';
allpitch = [OutData.signals.values(TS:nend,strmatch('BldPitch1',OutList)); OutData.signals.values(TS:nend,strmatch('BldPitch2',OutList)); OutData.signals.values(TS:nend,strmatch('BldPitch3',OutList))];
pitchrate = [diff(OutData.signals.values(TS:nend,strmatch('BldPitch1',OutList)))./diff(tout_trim); diff(OutData.signals.values(TS:nend,strmatch('BldPitch2',OutList)))./diff(tout_trim); diff(OutData.signals.values(TS:nend,strmatch('BldPitch3',OutList)))./diff(tout_trim)];
    %allpitch and pitchrate concatenate values for the three blades [1;2;3]
PM(22) = sqrt(mean(pitchrate.^2));



%IMUs (inertial measurement units)
PM_Desc{34,1}='RMS X Nacelle Accel';
PM(34) = sqrt(mean(OutData.signals.values(TS2:nend,strmatch('NcIMUTAxs',OutList)).^2));
PM_Desc{36,1}='RMS Y Nacelle Accel';
PM(36) = sqrt(mean(OutData.signals.values(TS2:nend,strmatch('NcIMUTAys',OutList)).^2));
PM_Desc{38,1}='RMS Z Nacelle Accel';
PM(38) = sqrt(mean(OutData.signals.values(TS2:nend,strmatch('NcIMUTAzs',OutList)).^2));

