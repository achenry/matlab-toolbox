%% Af_StateAso

%% Internal State Description from 'ReadFastInput.m'
% This is the order of states that come out of FAST Simulation
%Useful for doing analysis on the [Q,Qdot] outputs, such as the frequency
%tests with states as the output

State_Desc=cell(22,2);

State_Desc{1,1}='DOF_Sg';       State_Desc{1,2}='DOF index for platform surge.';                State_Desc{1,3}='DOF ptfm. surge.';
State_Desc{2,1}='DOF_Sw';       State_Desc{2,2}='DOF index for platform sway.';                 State_Desc{2,3}='DOF ptfm. sway.';
State_Desc{3,1}='DOF_Hv';       State_Desc{3,2}='DOF index for platform heave.';                State_Desc{3,3}='DOF ptfm. heave.';
State_Desc{4,1}='DOF_R';        State_Desc{4,2}='DOF index for platform roll.';                 State_Desc{4,3}='DOF ptfm. roll.';
State_Desc{5,1}='DOF_P';        State_Desc{5,2}='DOF index for platform pitch.';                State_Desc{5,3}='DOF ptfm. pitch.';
State_Desc{6,1}='DOF_Y';        State_Desc{6,2}='DOF index for platform yaw.';                  State_Desc{6,3}='DOF ptfm. yaw.';
State_Desc{7,1}='DOF_TFA1';     State_Desc{7,2}='DOF index for 1st tower fore-aft mode.';       State_Desc{7,3}='DOF 1st tower f-a mode.';
State_Desc{8,1}='DOF_TSS1';     State_Desc{8,2}='DOF index for 1st tower side-to-side mode.';   State_Desc{8,3}='DOF 1st tower s-s mode.'; 
State_Desc{9,1}='DOF_TFA2';     State_Desc{9,2}='DOF index for 2nd tower fore-aft mode.';       State_Desc{9,3}='DOF 2nd tower f-a mode.';
State_Desc{10,1}='DOF_TSS2';    State_Desc{10,2}='DOF index for 2nd tower side-to-side mode.';  State_Desc{10,3}='DOF 2nd tower s-s mode.';
State_Desc{11,1}='DOF_Yaw';     State_Desc{11,2}='DOF index for nacelle-yaw.';                  State_Desc{11,3}='DOF nacelle-yaw.';
State_Desc{12,1}='DOF_RFrl';    State_Desc{12,2}='DOF index for rotor-furl.';                   State_Desc{12,3}='DOF rotor-furl.';
State_Desc{13,1}='DOF_GeAz';    State_Desc{13,2}='DOF index for the generator azimuth.';        State_Desc{13,3}='DOF gen azimuth.';
State_Desc{14,1}='DOF_DrTr';    State_Desc{14,2}='DOF index for drivetrain rotational-flexibility.'; State_Desc{14,3}='DOF drivetrain rot.-flex.';
State_Desc{15,1}='DOF_TFrl';    State_Desc{15,2}='DOF index for tail-furl.';                         State_Desc{15,3}='DOF tail-furl.';
State_Desc{16,1}='DOF_BF(1,1)';   State_Desc{16,2}='DOF index for 1st blade flap mode for blade 1';  State_Desc{16,3}='DOF 1st flap mode-blade 1';
State_Desc{17,1}='DOF_BE(1,1)';   State_Desc{17,2}='DOF index for 1st blade edge mode for blade 1';  State_Desc{17,3}='DOF 1st edge mode-blade 1';
State_Desc{18,1}='DOF_BF(1,2)';   State_Desc{18,2}='DOF index for 2nd blade flap mode for blade 1';  State_Desc{18,3}='DOF 2nd flap mode-blade 1';
State_Desc{19,1}='DOF_BF(2,1)';   State_Desc{19,2}='DOF index for 1st blade flap mode for blade 2';  State_Desc{19,3}='DOF 1st flap mode-blade 2';
State_Desc{20,1}='DOF_BE(2,1)';   State_Desc{20,2}='DOF index for 1st blade edge mode for blade 2';  State_Desc{20,3}='DOF 1st edge mode-blade 2';
State_Desc{21,1}='DOF_BF(2,2)';   State_Desc{21,2}='DOF index for 2nd blade flap mode for blade 2';  State_Desc{21,3}='DOF 2nd flap mode-blade 2';
State_Desc{22,1}='DOF_BF(3,1)';   State_Desc{22,2}='DOF index for 1st blade flap mode for blade 3';  State_Desc{22,3}='DOF 1st flap mode-blade 3';
State_Desc{23,1}='DOF_BE(3,1)';   State_Desc{23,2}='DOF index for 1st blade edge mode for blade 3';  State_Desc{23,3}='DOF 1st edge mode-blade 3';
State_Desc{24,1}='DOF_BF(3,2)';   State_Desc{24,2}='DOF index for 2nd blade flap mode for blade 3';  State_Desc{24,3}='DOF 2nd flap mode-blade 3';
State_Desc{25,1}='DOF_Teet';      State_Desc{25,2}='DOF index for rotor-teeter.'; State_Desc{25,3}='DOF rotor-teeter.';





