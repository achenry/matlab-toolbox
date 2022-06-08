%% Run FAST from Simulink in parallel

%% TODO
% use generateCases function X
% create edit FASTIF function using struct and templates X
% adapt filepaths for mac and unix X

%% Initialize workspace and add paths
close all;
clc;
clearvars;

start_dir = pwd;

if ismac
    
    home_dir = '/Users/aoifework/Documents';
    project_dir = fullfile(home_dir, 'Research/ipc_tuning/');

    addpath(fullfile(home_dir, 'toolboxes/matlab-toolbox/A_Functions'));
    addpath(fullfile(home_dir, 'toolboxes/matlab-toolbox/Utilities'));
    addpath(fullfile(home_dir, 'toolboxes/turbsim-toolbox/A_Functions/'));
    addpath(project_dir);

    libext = '.dylib';
    
    % fast_install_dir = fullfile(home_dir, 'dev/WEIS/OpenFAST/install');
    fast_install_dir = fullfile(home_dir, 'usflowt_src/openfast/install');

    fast_models_dir = fullfile(home_dir, 'Research/ipc_tuning');
    fastRunner.FAST_runDirectory = fullfile(home_dir, 'Research/ipc_tuning/simulations');

    windfiles_dir = fullfile(project_dir, 'WindFiles', 'rated_turbulent');

    FAST_SimulinkModel = fullfile(project_dir, 'AD_SOAR_c7_V2f_c73_Clean');
    
    % FAST_SFunc location
    addpath(fullfile(fast_install_dir, 'lib'));

elseif isunix
    home_projects_dir = '/projects/aohe7145';
    home_storage_dir = '/scratch/summit/aohe7145';
    project_dir = fullfile(home_projects_dir, 'projects/ipc_tuning');

    addpath(fullfile(home_projects_dir, 'toolboxes/dev/matlab-toolbox/A_Functions'));
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox/Utilities'));
    addpath(fullfile(home_projects_dir, 'toolboxes/turbsim-toolbox/A_Functions/'));
    addpath(project_dir);

    libext = '.so';
    
    fast_install_dir = fullfile(home_projects_dir, 'toolboxes/dev/WEIS/OpenFAST/install');
    fast_models_dir = fullfile(home_projects_dir, 'models');
    fastRunner.FAST_runDirectory = fullfile(home_storage_dir, 'OpenfastSimulations/SOAR_rated_turbulent');
    windfiles_dir = fullfile(home_storage_dir, 'WindFiles/rated_turbulent');

    FAST_SimulinkModel = fullfile(home_dir, fast_models_dir, 'AD_SOAR_c7_V2f_c73_Clean');

    % FAST_SFunc location
    addpath(fullfile(fast_install_dir, 'lib'));

end

fastRunner.FAST_exe = fullfile(fast_install_dir, 'bin/openfast');
fastRunner.FAST_lib = fullfile(fast_install_dir, ['lib/libopenfastlib', libext]);
fastRunner.FAST_directory = fullfile(fast_models_dir, 'SOAR-25_V2f_IF');
fastRunner.FAST_InputFile = 'weis_job_00';

if ~exist(fastRunner.FAST_runDirectory, 'dir')
    mkdir(fastRunner.FAST_runDirectory);
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
    OutList = manualOutList([fullfile(fastRunner.FAST_directory, fastRunner.FAST_InputFile), '.SFunc.sum']);
end

TMax = 150;
n_seeds = 100;
ux_mean = 11.4;

CaseGen.dir_matrix = fastRunner.FAST_runDirectory;
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

[fastRunner.case_list, fastRunner.case_name_list, n_cases] = generateCases(case_basis, CaseGen.namebase, true);

%% Generate OpenFAST input files for each case

input_mode = 2;
def_fst_file = fullfile(fastRunner.FAST_directory, fastRunner.FAST_InputFile);
def_infw_file = fullfile(fastRunner.FAST_directory, [fastRunner.FAST_InputFile, '_InflowFile']);
templates_dir = fullfile(fastRunner.FAST_directory, 'templates');
template_fst_dir = fullfile(templates_dir, 'Fst');
template_infw_dir = fullfile(templates_dir, 'InflowWind');

copyfile(fullfile(fastRunner.FAST_directory, 'Airfoils'), fullfile(fastRunner.FAST_runDirectory, 'Airfoils'));
copyfile(fullfile(fastRunner.FAST_directory, '*.txt'), fastRunner.FAST_runDirectory);
copyfile(fullfile(fastRunner.FAST_directory, '*.dat'), fastRunner.FAST_runDirectory);
copyfile(fullfile(fastRunner.FAST_directory, '*.fst'), fastRunner.FAST_runDirectory);

for case_idx=1:n_cases
    
    fastRunner.case_list(case_idx).FAST_directory = fastRunner.FAST_runDirectory;
    % fastRunner.case_list(case_idx).FAST_runDirectory = fullfile(fastRunner.FAST_runDirectory, ['case_' num2str(case_idx)]);
    fastRunner.case_list(case_idx).FAST_runDirectory = fastRunner.case_list(case_idx).FAST_directory;
    if ~exist(fastRunner.case_list(case_idx).FAST_runDirectory, 'dir')
        mkdir(fastRunner.case_list(case_idx).FAST_runDirectory);
    end
    fastRunner.case_list(case_idx).FAST_inputFilename = ...
        fullfile(fastRunner.case_list(case_idx).FAST_runDirectory, [fastRunner.case_name_list{case_idx}, '.fst']);

    new_fst_name = fullfile(fastRunner.case_list(case_idx).FAST_runDirectory, ...
        fastRunner.case_name_list{case_idx});
    new_infw_name = fullfile(fastRunner.case_list(case_idx).FAST_runDirectory, ...
        [fastRunner.case_name_list{case_idx}, '_InflowWind']);

    fst_lines = fields(fastRunner.case_list(case_idx).Fst);
    fst_edits = {};

    for l = 1:length(fst_lines)
        fst_edits{l} = fastRunner.case_list(case_idx).Fst.(fst_lines{l});
    end

    if isfield(fastRunner.case_list, 'InflowWind')
        fst_lines{l + 1} = 'InflowFile';
        fst_edits{l + 1} = ['"' new_infw_name '.dat' '"'];
    end
    
    infw_lines = fields(fastRunner.case_list(case_idx).InflowWind);
    infw_edits = {};

    for l = 1:length(infw_lines)
        infw_edits{l} = fastRunner.case_list(case_idx).InflowWind.(infw_lines{l});
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

for case_idx = 1:length(n_cases)

    C_BL_SOAR25_V2f_c73_Clean;

    cd(fastRunner.case_list(case_idx).FAST_runDirectory);

    FAST_InputFileName = fastRunner.case_list(case_idx).FAST_inputFilename;

    % Populate thread parameters
    sim_inputs(case_idx) = Simulink.SimulationInput(FAST_SimulinkModel);
%     sim_inputs(case_idx).thread_id = case_idx;
    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('TMax', TMax);
    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('FAST_InputFileName', FAST_InputFileName);
    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('DT', Simulation.DT);
    sim_inputs(case_idx) = sim_inputs(case_idx).setModelParameter('StartTime','0','StopTime','TMax');
end

% Return to starting directory
cd(start_dir);


%% Run simulations in multiple parallel threads

RUN_SIM = true;

if RUN_SIM
%     Future = parsim(sim_inputs, ...
%          'RunInBackground', 'off', ...
%         'ShowProgress', 'on', ...
%         'TransferBaseWorkspaceVariables', true); % Run simulation
     Future = sim(FAST_SimulinkModel);
end

fprintf('\nOpenFAST Simulations Completed\n\n\n');

