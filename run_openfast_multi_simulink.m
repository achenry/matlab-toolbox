%% Run FAST from Simulink in parallel

%% TODO
% use generateCases function X
% create edit FASTIF function using struct and templates X
% adapt filepaths for mac and unix X

%% Initialize workspace and add paths
close all;
clc;
clearvars;

TMax = 700;
n_seeds = 100;

if ismac
    
    home_dir = '/Users/aoifework/Documents';
    project_dir = fullfile(home_dir, 'Research/ipc_tuning/');

    addpath(fullfile(home_dir, 'toolboxes/matlab-toolbox'));
    addpath(fullfile(home_dir, 'toolboxes/matlab-toolbox/A_Functions'));
    addpath(fullfile(home_dir, 'toolboxes/matlab-toolbox/Utilities'));
    addpath(fullfile(home_dir, 'toolboxes/turbsim-toolbox/A_Functions/'));
    addpath(project_dir);

    libext = '.dylib';
    
    % fast_install_dir = fullfile(home_dir, 'dev/WEIS/OpenFAST/install');
    fast_install_dir = fullfile(home_dir, 'usflowt_src/openfast/install');

    fast_models_dir = project_dir;
    FAST_runDirectory = fullfile(home_dir, 'Research/ipc_tuning/simulations');

    windfiles_dir = fullfile(project_dir, 'WindFiles', 'rated_turbulent');

    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean';
    addpath(fullfile(project_dir)); % sl model
    
    % FAST_SFunc location
    addpath(fullfile(fast_install_dir, 'lib'));

elseif isunix
    home_projects_dir = '/projects/aohe7145';
    home_storage_dir = '/scratch/summit/aohe7145';
    project_dir = fullfile(home_projects_dir, 'projects/ipc_tuning');
    
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox'));
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox/A_Functions'));
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox/Utilities'));
    addpath(fullfile(home_projects_dir, 'toolboxes/turbsim-toolbox/A_Functions/'));
    addpath(project_dir);

    libext = '.so';
    
    fast_install_dir = fullfile(home_projects_dir, 'toolboxes/usflowt/openfast/install');
    fast_models_dir = project_dir;
    FAST_runDirectory = fullfile(home_storage_dir, 'OpenfastSimulations/SOAR_rated_turbulent');
    windfiles_dir = fullfile(project_dir, 'WindFiles/rated_turbulent');

    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean';
    addpath(fullfile(project_dir)); % sl model

    % FAST_SFunc location
    addpath(fullfile(fast_install_dir, 'lib'));

end

fastRunner.FAST_exe = fullfile(fast_install_dir, 'bin/openfast');
fastRunner.FAST_lib = fullfile(fast_install_dir, ['lib/libopenfastlib', libext]);
fastRunner.FAST_directory = fullfile(fast_models_dir, 'SOAR-25_V2f_IF');
fastRunner.FAST_InputFile = 'weis_job_00';

if ~exist(FAST_runDirectory, 'dir')
    mkdir(FAST_runDirectory);
end

%% Generate simulation cases

Simulation.OpenFAST         = 1;
Simulation.OpenFAST_Date    = '042821';
Simulation.AeroDynVer       = '15';
Simulation.DT = 0.0125;
DT = Simulation.DT;
cut_transients              = 100;

Name_Model                  = 'AD_SOAR_c7_V2f_c73_Clean';
Name_Control                = 'C_BL_SOAR25_V2f_c73_Clean';
Parameters.Turbine.String   = 'SOAR-25-V2f';
Parameters.Turbine.fine_pitch = 0.303;
CpCtCqFile = fullfile(fastRunner.FAST_directory, 'weis_job_00_Cp_Ct_Cq.txt');
Parameters.Turbine.ConeAngle = 3.6;%2.5;
Parameters.Turbine.ShaftTilt = 8.4;%7.0; % deg
Parameters.Tower.Height  = 193.287;

if ~exist('OutList')
    OutList = manualOutList([fullfile(fastRunner.FAST_directory, ...
        fastRunner.FAST_InputFile), '.SFunc.sum']);
end

ux_mean = 11.4;

CaseGen.dir_matrix = FAST_runDirectory;
CaseGen.namebase = 'AD_SOAR_c7_V2f_c73_Clean';
[~, CaseGen.model_name, ~] = fileparts(fastRunner.FAST_InputFile);

case_basis.InflowWind.FileName_BTS = {};
case_basis.InflowWind.WindType = {'3'};
case_basis.InflowWind.HWindSpeed = {num2str(ux_mean)};
case_basis.Fst.TMax = {num2str(TMax)};

for bts_idx = 1:n_seeds
    case_basis.InflowWind.FileName_BTS{bts_idx} = ['"' fullfile(windfiles_dir, ...
        ['B_', replace(num2str(ux_mean), '.', '-'), '_', num2str(bts_idx), '.bts']) '"'];
end

[case_list, case_name_list, n_cases] = generateCases(case_basis, CaseGen.namebase, true);

%% Generate OpenFAST input files for each case

input_mode = 2;
def_fst_file = fullfile(fastRunner.FAST_directory, fastRunner.FAST_InputFile);
def_infw_file = fullfile(fastRunner.FAST_directory, [fastRunner.FAST_InputFile, '_InflowFile']);
templates_dir = fullfile(fastRunner.FAST_directory, 'templates');
template_fst_dir = fullfile(templates_dir, 'Fst');
template_infw_dir = fullfile(templates_dir, 'InflowWind');

copyfile(fullfile(fastRunner.FAST_directory, 'Airfoils'), fullfile(FAST_runDirectory, 'Airfoils'));
copyfile(fullfile(fastRunner.FAST_directory, '*.txt'), FAST_runDirectory);
copyfile(fullfile(fastRunner.FAST_directory, '*.dat'), FAST_runDirectory);
copyfile(fullfile(fastRunner.FAST_directory, '*.fst'), FAST_runDirectory);


parfor case_idx=1:n_cases
    
    % case_list(case_idx).FAST_directory = FAST_runDirectory;
    % case_list(case_idx).FAST_runDirectory = fullfile(FAST_runDirectory, ['case_' num2str(case_idx)]);
    % case_list(case_idx).FAST_runDirectory = FAST_runDirectory;
%     if ~exist(FAST_runDirectory, 'dir')
%         mkdir(FAST_runDirectory);
%     end
    %case_list(case_idx).FAST_inputFilename = ...
    %    fullfile(case_list(case_idx).FAST_runDirectory, [case_name_list{case_idx}, '.fst']);

    new_fst_name = fullfile(FAST_runDirectory, ...
        case_name_list{case_idx});
    new_infw_name = fullfile(FAST_runDirectory, ...
        [case_name_list{case_idx}, '_InflowWind']);

    fst_lines = fields(case_list(case_idx).Fst);
    fst_edits = {};

    for l = 1:length(fst_lines)
        fst_edits{l} = case_list(case_idx).Fst.(fst_lines{l});
    end

    if isfield(case_list, 'InflowWind')
        fst_lines{l + 1} = 'InflowFile';
        fst_edits{l + 1} = ['"' new_infw_name '.dat' '"'];
    end
    
    infw_lines = fields(case_list(case_idx).InflowWind);
    infw_edits = {};

    for l = 1:length(infw_lines)
        infw_edits{l} = case_list(case_idx).InflowWind.(infw_lines{l});
    end
    
    Af_EditFast(fst_lines, fst_edits, new_fst_name, def_fst_file, template_fst_dir, input_mode);
    Af_EditInflow(infw_lines, infw_edits, new_infw_name, def_infw_file, template_infw_dir, input_mode);
end

% Af_EditSub;
% Af_EditServo;
% Af_EditLin;
% Af_EditHydro;
% Af_EditElast;
% Af_EditDriver;
% Af_EditBeam;
% Af_EditAero15;
% Af_EditADriver;

%% Prepare parallel SL Simulations
cd(project_dir);
for case_idx = 1:n_cases

    C_BL_SOAR25_V2f_c73_Clean;

    

    FAST_InputFileName = fullfile(FAST_runDirectory, ...
        [case_name_list{case_idx}, '.fst']);

    % Populate thread parameters
    sim_inputs(case_idx) = Simulink.SimulationInput(FAST_SimulinkModel);
%     sim_inputs(case_idx).thread_id = case_idx;
    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('TMax', TMax);
    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('FAST_InputFileName', FAST_InputFileName);
    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('DT', Simulation.DT);
    sim_inputs(case_idx) = sim_inputs(case_idx).setModelParameter('StartTime','0','StopTime','TMax');
end

% Return to starting directory
%cd(start_dir);


%% Run simulations in multiple parallel threads

RUN_SIM = true;

if RUN_SIM
    Future = parsim(sim_inputs, ...
                   'TransferBaseWorkspaceVariables', true, ... % Run simulation
                   'RunInBackground', 'off', ...
                   'ShowProgress',  'on');
        
     %Future = sim(FAST_SimulinkModel);
end
% wait(Future);
% test_idx = 1;
% sim_data = Future(test_idx);
% time = sim_data.OutData.time;
% BldPitch1 = sim_data.Cyc_BldPitch1;
% GenTorq = sim_data.Out_Tg;
% % Torque = sim_data.sigsOut{14}.Values.Data;
% RotSpeed = sim_data.sigsOut{9}.Values.Data;
% Uest = sim_data.v_hat_out.Data;
% 
% figure;
% tiledlayout(4, 1);
% ax1 = nexttile;
% plot(time, BldPitch1); title('BldPitch1');
% ax2 = nexttile
% plot(time, GenTorq); title('GenTorq');
% ax3 = nexttile;
% plot(time, RotSpeed); title('RotSpeed');
% ax4 = nexttile;
% plot(time, Uest); title('Uest');
% linkaxes([ax1, ax2, ax3, ax4], 'x');
% xlabel('time');