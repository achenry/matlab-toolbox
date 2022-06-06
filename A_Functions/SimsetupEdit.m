% This script creates necessary workspace variables to run a Simulink model
% using the FAST dynamics and aerodynamic S-function block.  Before running
% a simulation, the character array, input_fast, must contain the FAST
% input file name, and the script Read_FAST_Input.m must run.
%clc;
%clear all;
%option instead of clear all (FD 3/15/09): clear functions;
%clear('-regexp', '^[^A]'); %clears all vars that don't start with A (to
%save AllPM)

% Prompt the user for the input file name.
%commented out FD 03/15/09:[input_fast,PathName] = uigetfile('*.fst','Select FAST model');
% input_fast = 'FAST_IF/FASTinputfile.fst';
% Name_FAST_New
% Read FAST input file and set initial conditions
Read_FAST_Input

% -------------------------------------------------------------------------
%c'd out FD 3/16/09: addpath Subs;   % access to functions located in the ./Subs folder
% Adding the path of the GUI files for easy access to the GUI by just
% calling FASTDataPlotter from the command window
%c'd out FD 3/16/09: addpath GUI
% JPA- Took Filter variables and combined with PIDPitchController and
% renamed as BaselineController.m

% % Global variables needed for the azimuth_unwrap.m function
% global az_prev counter;
% az_prev = 0;
% counter = 0;
% 
% % Rated properties of the turbine
% TgRated = 43.09355;         % rated generator torque [kNm]
% Prated = 5000;              % rated power [kW]
% w_rated_rpm = 12.1;         % rated rotor speed [rpm]
