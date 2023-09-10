%% Run FAST from Simulink in parallel
% set(findall(figureHandle,'type','text'),'fontSize',14,'fontWeight','bold')

%% Initialize workspace and add paths
close all;
clc;
% clearvars;
clear all;
restoredefaultpath;

TMax = 700;
n_seeds = 100;
RUN_OL_DQ = 0; % ['false']
RUN_OL_BLADE = 0;
RUN_CL = 1;
WIND_TYPE = 'steady';
% WIND_SPEEDS = 10:0.5:22;
% WIND_SPEEDS = 0:0.5:24;
% WIND_TYPE = 'turbsim';
WIND_SPEEDS = [14];
RUN_SIMS_PAR = 0;
RUN_SIMS_SINGLE = 1;
COMPUTE_SS_VALS = 1;
SWEEP_ANALYSIS = 0;
LIN_SENS_ANALYSIS = 0;
NONLIN_SENS_ANALYSIS = 0;
SETUP_NONLIN_HINF = 0;
SETUP_LIN_HINF = 0;

% if sum([RUN_OL_DQ, RUN_OL_BLADE, RUN_CL]) ~= 1
%     error('Choose run type');
% end

% set default plot settings
set(groot, 'defaultAxesFontSize', 16);
% get(groot, 'default'), get(groot, 'factory')

% blade_ops = [35; 30; 30];
% theta = 90;
% c_comp = (1/3) .* sum(blade_ops, 2);
% d_comp = (2/3) .* [cosd(theta) cosd(theta+120) cosd(theta+240)] * blade_ops;
% q_comp =  (2/3) .* [sind(theta) sind(theta+120) sind(theta+240)] * blade_ops;


if RUN_OL_DQ
    DistAmp = deg2rad(15);
    RampRate = DistAmp; %1e-3;
    SettleTime = 100;

    BaselineSteadyState = SettleTime:(2 * SettleTime);

    dRampStart = 2 * SettleTime + 1; % leave 100 seconds to settle, then another 100 to get mean and max baseline
    dRampStop = dRampStart + (DistAmp / RampRate); % constant by 110
    dSteadyState = dRampStop + SettleTime:dRampStop + (2 * SettleTime);
    
    qRampStart = dRampStop + (2 * SettleTime) + 1; % leave 100 seconds to settle, then another 100 to get mean and max
    qRampStop = qRampStart + (DistAmp / RampRate); % constant by 220
    qSteadyState = qRampStop + SettleTime:qRampStop + (2 * SettleTime);
    
    TMax = qSteadyState(end); % leave 100 seconds to settle, then another 100 to get mean and max
elseif RUN_OL_BLADE
    DistAmp = deg2rad(15); %1e-3; % rad
    RampRate = DistAmp; %1e-3;
    SettleTime = 100;

    BaselineSteadyState = SettleTime:(2 * SettleTime);

    b1RampStart = 2 * SettleTime + 1; % leave 100 seconds to settle, then another 100 to get mean and max baseline
    b1RampStop = b1RampStart + (DistAmp / RampRate); % constant by 110
    b1SteadyState = b1RampStop + SettleTime:b1RampStop + (2 * SettleTime);
    
    TMax = b1SteadyState(end); % leave 100 seconds to settle, then another 100 to get mean and max
end

if ismac % TODO setup for Manuel to run LIN_SYS_ANALYSIS
    
    home_dir = '/Users/aoifework/Documents';
    project_dir = fullfile(home_dir, 'Research', 'ipc_tuning');
    % project_dir = fullfile(home_dir, 'Research', 'learning_actuation');
    simulink_model_dir = fullfile(project_dir, 'simulink_models');
    fig_dir = fullfile(project_dir, 'figs');
    toolbox_dir = fullfile(home_dir, 'toolboxes');
    
    addpath(simulink_model_dir);

    addpath(fullfile(toolbox_dir, 'matlab-toolbox'));
    addpath(fullfile(toolbox_dir, 'matlab-toolbox/A_Functions'));
    addpath(fullfile(toolbox_dir, 'matlab-toolbox/Utilities'));
    addpath(fullfile(toolbox_dir, 'matlab-toolbox', 'MBC', 'Source'));
    addpath(fullfile(toolbox_dir, 'turbsim-toolbox/A_Functions/'));
    addpath(fullfile(toolbox_dir, 'PMtools/'));
    autrun;

    addpath(project_dir);

    libext = '.dylib';
    
    % fast_install_dir = fullfile(home_dir, 'dev/WEIS/OpenFAST/install');
    fast_install_dir = fullfile(toolbox_dir, 'openfast/install');
    % fast_install_dir = fullfile(home_dir, 'src/openfast/install');
    % fast_install_dir = fullfile(home_dir, 'usflowt_src/openfast/sl_install');
    addpath(fast_install_dir);
    % addpath(fullfile(fast_install_dir, 'include'));
    % addpath(fullfile(fast_install_dir, 'lib'));

    FAST_runDirectory = fullfile(project_dir, 'ss_simulations');
    if ~exist(FAST_runDirectory)
        mkdir(FAST_runDirectory);
    end
    addpath(FAST_runDirectory);

    windfiles_dir = fullfile(project_dir, 'WindFiles', 'rated_turbulent');
    
    addpath(fullfile(project_dir)); % sl model
    
    % FAST_SFunc location
    addpath(fullfile(fast_install_dir, 'lib'));
    % addpath(fullfile(toolbox_dir, 'openfast/build/glue-codes/simulink'));

elseif isunix
    home_projects_dir = '/projects/aohe7145';
    home_storage_dir = '/scratch/alpine/aohe7145/ipc_tuning';
    project_dir = fullfile(home_projects_dir, 'projects/ipc_tuning');
    simulink_model_dir = fullfile(project_dir, 'simulink_models');
    addpath(simulink_model_dir);
    % plant_setup_dir = fullfile(project_dir, 'plant_setup_package/');
    
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox'));
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox/A_Functions'));
    addpath(fullfile(home_projects_dir, 'toolboxes/matlab-toolbox/Utilities'));
    addpath(fullfile(home_projects_dir, 'toolboxes', 'matlab-toolbox', 'MBC', 'Source'));
    addpath(fullfile(home_projects_dir, 'toolboxes/turbsim-toolbox/A_Functions/'));
    addpath(project_dir);
    autrun;
    libext = '.so';
    
    fast_install_dir = fullfile(home_projects_dir, 'toolboxes/openfast_dev/install');
    FAST_runDirectory = fullfile(home_storage_dir, 'ss_simulations');
    windfiles_dir = fullfile(project_dir, 'WindFiles/rated_turbulent');
    
    if ~exist(FAST_runDirectory)
        mkdir(FAST_runDirectory);
    end
    addpath(FAST_runDirectory);

    addpath(fullfile(project_dir)); % sl model

    % FAST_SFunc location
    addpath(fast_install_dir);
    addpath(fullfile(fast_install_dir, 'lib'));
end

FAST_SimulinkModel_dir = simulink_model_dir;
if RUN_OL_DQ
    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean_OL_DQ';
elseif RUN_OL_BLADE
    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean_OL_BLADE';
elseif RUN_CL && SWEEP_ANALYSIS
    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean_newIPC';
elseif RUN_CL && COMPUTE_SS_VALS
    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean';
else
    FAST_SimulinkModel = 'AD_SOAR_c7_V2f_c73_Clean_newIPC';
end

fastRunner.FAST_exe = fullfile(fast_install_dir, 'bin/openfast');
fastRunner.FAST_lib = fullfile(fast_install_dir, ['lib/libopenfastlib', libext]);
fastRunner.FAST_InputFile = 'weis_job_00';
% fastRunner.FAST_InputFile = 'IEA-15-240-RWT-UMaineSemi';

if ~exist(FAST_runDirectory, 'dir')
    mkdir(FAST_runDirectory);
end


%% Generate simulation cases

Simulation.OpenFAST         = 1;
Simulation.OpenFAST_Date    = '042821';
Simulation.AeroDynVer       = '15';
Simulation.DT = 0.0125;
% Simulation.DT = 0.025;
DT = Simulation.DT;
cut_transients              = 100;

% Name_Model                  = 'AD_SOAR_c7_V2f_c73_Clean';
% Name_Control                = 'C_BL_SOAR25_V2f_c73_Clean';
Parameters.Turbine.String   = 'SOAR-25-V2f';
% Parameters.Turbine.String = 'IEA-15-240-RWT-UMaineSemi';
Parameters.Turbine.fine_pitch = 0.303;  % deg

Parameters.Turbine.ConeAngle = 3.6;%2.5;% deg
Parameters.Turbine.ShaftTilt = 8.4;%7.0; % deg
Parameters.Tower.Height  = 193.287;  % meters

fastRunner.FAST_directory = fullfile(project_dir, [Parameters.Turbine.String '_IF']);
FAST_InputFileName = fullfile(fastRunner.FAST_directory, [fastRunner.FAST_InputFile '.fst']);
CpCtCqFile = fullfile(fastRunner.FAST_directory, 'weis_job_00_Cp_Ct_Cq.txt');
% CpCtCqFile = fullfile(fastRunner.FAST_directory, 'Cp_Ct_Cq.IEA15MW.txt');

addpath(fastRunner.FAST_directory);

lin_models_dir = fullfile(fastRunner.FAST_directory, 'linearization', 'steady_wind-CL', 'linfiles', 'excGenDOF_incSecOrd');
% lin_models_dir = fullfile(fastRunner.FAST_directory, 'linearization', 'steady_wind-CL', 'linfiles');

C_BL_SOAR25_V2f_c73_Clean;
if COMPUTE_SS_VALS
    Parameters.Control.IPCDQ.Enable = 0;
end

% ux_mean = 11.4;
ux_mean = 14;

if sum([RUN_OL_DQ, RUN_OL_BLADE, RUN_CL]) == 1
    CaseGen.dir_matrix = FAST_runDirectory;
    CaseGen.namebase = FAST_SimulinkModel;
    [~, CaseGen.model_name, ~] = fileparts(fastRunner.FAST_InputFile);
    
    if strcmp(WIND_TYPE, 'turbsim')
        case_basis.InflowWind.WindType = {'3'};
        case_basis.InflowWind.FileName_BTS = {};
        for bts_idx = 1:n_seeds
            case_basis.InflowWind.FileName_BTS{bts_idx} = ['"' fullfile(windfiles_dir, ...
                ['B_', replace(num2str(ux_mean), '.', '-'), '_', num2str(bts_idx), '.bts']) '"'];
      
        end
        case_basis.InflowWind.HWindSpeed = {num2str(ux_mean)};
        
    %     case_basis.InflowWind.HWindSpeed = num2str(ux_mean);
    elseif strcmp(WIND_TYPE, 'steady')
        case_basis.InflowWind.WindType = {'1'};
        case_basis.InflowWind.HWindSpeed = split(num2str(WIND_SPEEDS));
    end
    
    case_basis.Fst.TMax = {num2str(TMax)};
    
    % Integrator Gain Sweep TODO
%     case_basis.Parameters = Parameters;
%     case_basis.Parameters.cIPC.DQ_Ki_1P = logspace(-8, 2, 10);
%     Parameters.cIPC.DQ_Ki_1P = kmul*30e-8;

    [case_list, case_name_list, n_cases] = generateCases(case_basis, CaseGen.namebase, true);
    
    % if OutList doesn't exist, run
    % /Users/aoifework/Documents/dev/WEIS/OpenFAST/install/bin/openfast
    % weis_job_00.fst from command line with all Modes in ServoDyn file set
    % to 0 for TMax=10 to output OutList
    if ~exist('OutList.mat') || true
        % OutList = manualOutList([fullfile(FAST_runDirectory, ...
        %     case_name_list{1}), '.sum']);
        OutList = manualOutList(fullfile(fastRunner.FAST_directory, 'weis_job_00_cmdline.sum'));
        save(fullfile(project_dir, 'OutList.mat'), 'OutList');
    else
        load OutList.mat
    
    %     OutList = manualOutList([fullfile(fastRunner.FAST_directory, ...
    %         fastRunner.FAST_InputFile), '.SFunc.sum']);
    end
end
if exist('OutList.mat')
    load OutList.mat
end

dqOutList = OutList;
for op_label = OutList'
    % if rotating quantity
    if strcmp(op_label{1}(end), '1')
        % get all corresponding quantities
        blade_op_labels = cellfun(@(b) [op_label{1}(1:end-1) b], {'1', '2', '3'}, 'UniformOutput', false);
        cdq_op_labels = cellfun(@(b) [op_label{1}(1:end-1) b], {'C', 'D', 'Q'}, 'UniformOutput', false);
        
        % replace in transformed signal matrix
        dqOutList(ismember(dqOutList, blade_op_labels)) = cdq_op_labels;
    end
end
save(fullfile(project_dir, 'dqOutList.mat'), 'dqOutList');

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

if RUN_SIMS_PAR
    % load(fullfile(fastRunner.FAST_directory, 'ss_vals'));
    parfor case_idx=1:n_cases

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
elseif RUN_SIMS_SINGLE
    % load(fullfile(fastRunner.FAST_directory, 'ss_vals'));
    case_idx = 1;%find(WIND_SPEEDS == 14);
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
% TODO make this more systematic, only works if single Fst set of input
% files
if SWEEP_ANALYSIS
origParameters = Parameters;
sweepParameters = origParameters;
IPCDQEnable = [0 1];
IPC3DQEnable = [0 1];
load(fullfile(project_dir, 'controller_setup_package', 'i_ctrl_range.mat')) % K_1PD_range, K_1PQ_range, K_3PD_range, K_3PQ_range
n_gains = 5;
Ki_1PD_gains = [mean([K_1PD_first_s(1), K_1PD_range(1, 1)])]; 
Ki_1PQ_gains = [mean([K_1PQ_first_s(1), K_1PQ_range(1, 1)])]; 
Ki_3PD_gains = [mean([K_3PD_first_s(1), K_3PD_range(1, 1)])]; 
Ki_3PQ_gains = [mean([K_3PQ_first_s(1), K_3PQ_range(1, 1)])];
% Ki_1PD_gains = [1e-6];
% Ki_1PQ_gains = [1e-6];
% Ki_3PD_gains = [1];
% Ki_3PQ_gains = [1];
% K_1PD_range = [1e-8, 1e-7, 1e-5, 1e-4];
% K_1PQ_range = [1e-8, 1e-7, 1e-5, 1e-4];
% K_3PD_range = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
% K_3PQ_range = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
if true
for r = 1:size(K_1PD_range, 1)
    Ki_1PD_gains = [Ki_1PD_gains, linspace(K_1PD_range(r, 1), K_1PD_range(r, 2), n_gains)];
end
for r = 1:size(K_1PQ_range, 1)
    Ki_1PQ_gains = [Ki_1PQ_gains, linspace(K_1PQ_range(r, 1), K_1PQ_range(r, 2), n_gains)];
end
for r = 1:size(K_3PD_range, 1)
    Ki_3PD_gains = [Ki_3PD_gains, linspace(K_3PD_range(r, 1), K_3PD_range(r, 2), n_gains)];
end
for r = 1:size(K_3PQ_range, 1)
    Ki_3PQ_gains = [Ki_3PQ_gains, linspace(K_3PQ_range(r, 1), K_3PQ_range(r, 2), n_gains)];
end
end
% gains = logspace(-10, 4, 14);
sweepParameters.cIPC.D_Ki_1P = []; % repelem(gains, length(repelem));
sweepParameters.cIPC.Q_Ki_1P = [];
sweepParameters.cIPC.D_Ki_3P = [];
sweepParameters.cIPC.Q_Ki_3P = [];
sweepParameters.Control.IPCDQ.Enable = [];
sweepParameters.Control.IPC3DQ.Enable = [];
for e1p = IPCDQEnable
    for e3p = IPC3DQEnable
        if (e1p == 0) && (e3p == 0) % if noIPC => consider only first gain for both IPC controllers
            sweepParameters.cIPC.D_Ki_1P = [sweepParameters.cIPC.D_Ki_1P 0]; %repelem(gains, length(repelem));
            sweepParameters.cIPC.Q_Ki_1P = [sweepParameters.cIPC.Q_Ki_1P 0];
            sweepParameters.cIPC.D_Ki_3P = [sweepParameters.cIPC.D_Ki_3P 0];
            sweepParameters.cIPC.Q_Ki_3P = [sweepParameters.cIPC.Q_Ki_3P 0];
            sweepParameters.Control.IPCDQ.Enable = [sweepParameters.Control.IPCDQ.Enable e1p];
            sweepParameters.Control.IPC3DQ.Enable = [sweepParameters.Control.IPC3DQ.Enable e3p];
        elseif (e1p == 1) && (e3p == 0) % if 1PIPC => consider only first gain for 3P IPC controller and all 1P IPC gains
            % test full range of 1PD gains for 0 3P gains and stable 1PQ
            % gain
            for k1pd = Ki_1PD_gains
                sweepParameters.cIPC.D_Ki_1P = [sweepParameters.cIPC.D_Ki_1P k1pd]; %repelem(gains, length(repelem));
                sweepParameters.cIPC.Q_Ki_1P = [sweepParameters.cIPC.Q_Ki_1P Ki_1PQ_gains(1)]; 
                sweepParameters.cIPC.D_Ki_3P = [sweepParameters.cIPC.D_Ki_3P 0];
                sweepParameters.cIPC.Q_Ki_3P = [sweepParameters.cIPC.Q_Ki_3P 0]; 
                sweepParameters.Control.IPCDQ.Enable = [sweepParameters.Control.IPCDQ.Enable e1p];
                sweepParameters.Control.IPC3DQ.Enable = [sweepParameters.Control.IPC3DQ.Enable e3p];
            end
            % test full range of 1PQ gains for 0 3P gains and stable 1PD
            % gain
            for k1pq = Ki_1PQ_gains(2:end)
                sweepParameters.cIPC.Q_Ki_1P = [sweepParameters.cIPC.Q_Ki_1P k1pq]; %repelem(gains, length(repelem));
                sweepParameters.cIPC.D_Ki_1P = [sweepParameters.cIPC.D_Ki_1P Ki_1PD_gains(1)];
                sweepParameters.cIPC.D_Ki_3P = [sweepParameters.cIPC.D_Ki_3P 0];
                sweepParameters.cIPC.Q_Ki_3P = [sweepParameters.cIPC.Q_Ki_3P 0];
                sweepParameters.Control.IPCDQ.Enable = [sweepParameters.Control.IPCDQ.Enable e1p];
                sweepParameters.Control.IPC3DQ.Enable = [sweepParameters.Control.IPC3DQ.Enable e3p];
            end
        elseif (e1p == 0) && (e3p == 1) % if only 3PIPC => consider only first gain for 1P IPC controller and all 3P IPC gains
            % test full range of 3PD gains for 0 1P gains and stable 3PQ
            % gain
            for k3pd = Ki_3PD_gains
                sweepParameters.cIPC.D_Ki_3P = [sweepParameters.cIPC.D_Ki_3P k3pd]; %repelem(gains, length(repelem));
                sweepParameters.cIPC.Q_Ki_3P = [sweepParameters.cIPC.Q_Ki_3P Ki_3PQ_gains(1)]; 
                sweepParameters.cIPC.D_Ki_1P = [sweepParameters.cIPC.D_Ki_1P 0];
                sweepParameters.cIPC.Q_Ki_1P = [sweepParameters.cIPC.Q_Ki_1P 0]; 
                sweepParameters.Control.IPCDQ.Enable = [sweepParameters.Control.IPCDQ.Enable e1p];
                sweepParameters.Control.IPC3DQ.Enable = [sweepParameters.Control.IPC3DQ.Enable e3p];
            end
            % test full range of 3PQ gains for 0 1P gains and stable 3PD
            % gain
            for k3pq = Ki_3PQ_gains(2:end)
                sweepParameters.cIPC.Q_Ki_3P = [sweepParameters.cIPC.Q_Ki_3P k3pq]; %repelem(gains, length(repelem));
                sweepParameters.cIPC.D_Ki_3P = [sweepParameters.cIPC.D_Ki_3P Ki_3PD_gains(1)];
                sweepParameters.cIPC.D_Ki_1P = [sweepParameters.cIPC.D_Ki_1P 0];
                sweepParameters.cIPC.Q_Ki_1P = [sweepParameters.cIPC.Q_Ki_1P 0];
                sweepParameters.Control.IPCDQ.Enable = [sweepParameters.Control.IPCDQ.Enable e1p];
                sweepParameters.Control.IPC3DQ.Enable = [sweepParameters.Control.IPC3DQ.Enable e3p];
            end
        end
    end
end

n_param_cases = length(sweepParameters.cIPC.D_Ki_1P);
param_idx = 1;
end

if RUN_SIMS_PAR
    cd(project_dir);
    if SWEEP_ANALYSIS
        for case_idx = 1:n_cases
            for param_idx = 1:n_param_cases
        
                FAST_InputFileName = fullfile(FAST_runDirectory, ...
                    [case_name_list{case_idx}, '.fst']);
            
                % Populate thread parameters
                sim_inputs(case_idx * param_idx) = Simulink.SimulationInput(FAST_SimulinkModel);
                sim_inputs(case_idx * param_idx) = sim_inputs(case_idx * param_idx).setVariable('TMax', TMax);
                sim_inputs(case_idx * param_idx) = sim_inputs(case_idx * param_idx).setVariable('FAST_InputFileName', FAST_InputFileName);
                sim_inputs(case_idx * param_idx) = sim_inputs(case_idx * param_idx).setVariable('DT', Simulation.DT);
                
                Parameters = origParameters;
                Parameters.cIPC.D_Ki_1P = sweepParameters.cIPC.D_Ki_1P(param_idx);
                Parameters.cIPC.Q_Ki_1P = sweepParameters.cIPC.Q_Ki_1P(param_idx);
                Parameters.cIPC.D_Ki_3P = sweepParameters.cIPC.D_Ki_3P(param_idx);
                Parameters.cIPC.Q_Ki_3P = sweepParameters.cIPC.Q_Ki_3P(param_idx);
                Parameters.Control.IPCDQ.Enable = sweepParameters.Control.IPCDQ.Enable(param_idx);
                Parameters.Control.IPC3DQ.Enable = sweepParameters.Control.IPC3DQ.Enable(param_idx);
                
                sim_inputs(case_idx * param_idx) = sim_inputs(case_idx * param_idx).setVariable('Parameters', Parameters);
        
                if strcmp(WIND_TYPE, 'steady')
                    sim_inputs(case_idx * param_idx) = sim_inputs(case_idx * param_idx).setVariable('HWindSpeed', str2num(case_list(case_idx).InflowWind.HWindSpeed));
                end
            
            end
        end
    else %elseif COMPUTE_SS_VALS
%         Parameters = origParameters;
        for case_idx = 1:n_cases
        
                FAST_InputFileName = fullfile(FAST_runDirectory, ...
                    [case_name_list{case_idx}, '.fst']);
            
                % Populate thread parameters
                sim_inputs(case_idx) = Simulink.SimulationInput(FAST_SimulinkModel);
                sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('TMax', TMax);
                sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('FAST_InputFileName', FAST_InputFileName);
                sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('DT', Simulation.DT);
                
                sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('Parameters', Parameters);
        
%                 if strcmp(WIND_TYPE, 'steady')
                    sim_inputs(case_idx) = sim_inputs(case_idx).setVariable('HWindSpeed', str2num(case_list(case_idx).InflowWind.HWindSpeed));
%                 end
        end
    end



    %% Run simulations in multiple parallel threads
    
    sim_out_list = parsim(sim_inputs, ...
                   'TransferBaseWorkspaceVariables', true, ... % Run simulation
                   'RunInBackground', 'off', ...
                   'ShowProgress',  'on', ...
                   'ShowSimulationManager', 'on');

    %% Perform Transformations on OutData

%     load(fullfile(project_dir, 'sim_out_list_ol_dq.mat'));
    if RUN_OL_BLADE || RUN_OL_DQ
        load(fullfile(fastRunner.FAST_directory, 'op_absmax'));
    end
    for c = 1:length(sim_out_list)
        sim_out_list(c).OutData.signals.dqValues = mbcTransformOutData(sim_out_list(c).OutData.signals.values, OutList);
        if RUN_OL_BLADE || RUN_OL_DQ
            sim_out_list(c).OutData.signals.normalizedValues = normalizeOutData(sim_out_list(c).OutData.signals.values, table2array(op_absmax.blade(c, :)));
            sim_out_list(c).OutData.signals.dqNormalizedValues = normalizeOutData(sim_out_list(c).OutData.signals.dqValues, table2array(op_absmax.dq(c, :)));
    
        end
    end
elseif RUN_SIMS_SINGLE
    % run single case
    case_idx = 1;
    FAST_InputFileName = fullfile(FAST_runDirectory, ...
        [case_name_list{case_idx}, '.fst']);
    DT = Simulation.DT;
    HWindSpeed = str2num(case_list(case_idx).InflowWind.HWindSpeed);
    clear OutList;
    clear DT;
%     Parameters = origParameters;
%     Parameters.cIPC.DQ_Ki_1P = sweepParameters.cIPC.DQ_Ki_1P(case_idx); % TODO change to be automated
    sim_out_list = sim(fullfile(FAST_SimulinkModel_dir, FAST_SimulinkModel), [0, TMax]);%, 'StartTime', '0', 'StopTime', 'TMax');
    
end

%% Save/Load Simulation Data
if RUN_SIMS_PAR || RUN_SIMS_SINGLE
    if RUN_CL && COMPUTE_SS_VALS% strcmp(WIND_TYPE, 'steady') && ~(RUN_OL_DQ || RUN_OL_BLADE)
        save(fullfile(project_dir, 'sim_out_list_ss.mat'), 'sim_out_list', '-v7.3');
    elseif RUN_CL && SWEEP_ANALYSIS && strcmp(WIND_TYPE, 'steady')
        save(fullfile(project_dir, 'sim_out_list_sweep_steady.mat'), 'sim_out_list', '-v7.3');
    elseif RUN_CL && SWEEP_ANALYSIS && strcmp(WIND_TYPE, 'turbsim')
        save(fullfile(project_dir, 'sim_out_list_sweep_turbsim.mat'), 'sim_out_list', '-v7.3');
    elseif RUN_OL_DQ
        save(fullfile(project_dir, 'sim_out_list_ol_dq.mat'), 'sim_out_list', '-v7.3');
    elseif RUN_OL_BLADE
        save(fullfile(project_dir, 'sim_out_list_ol_blade.mat'), 'sim_out_list', '-v7.3');
    elseif RUN_CL && strcmp(WIND_TYPE, 'turbsim')
        save(fullfile(project_dir, 'sim_out_list_turbsim.mat'), 'sim_out_list', '-v7.3');
    end
elseif COMPUTE_SS_VALS
    % get closed-loop to steady-state simulations from memory
    load(fullfile(project_dir, 'sim_out_list_ss.mat'));
elseif SWEEP_ANALYSIS && strcmp(WIND_TYPE, 'steady')
    load(fullfile(project_dir, 'sim_out_list_sweep_steady.mat'));
elseif SWEEP_ANALYSIS && strcmp(WIND_TYPE, 'turbsim')
    load(fullfile(project_dir, 'sim_out_list_sweep_turbsim.mat'));
elseif NONLIN_SENS_ANALYSIS && RUN_OL_DQ
    % get open-loop sensitivity to IPC simulations from memory
    load(fullfile(project_dir, 'sim_out_list_ol_dq.mat'));
elseif NONLIN_SENS_ANALYSIS && RUN_OL_BLADE
    % get open-loop sensitivity to IPC simulations from memory
    load(fullfile(project_dir, 'sim_out_list_ol_blade.mat'));
end


%% Sweep Data Analysis QUESTION what else to analyze?
if SWEEP_ANALYSIS
    
    
    % QUESTION MANUEL: are these the loads Daniel was referring to?
    % 'LSSGagMya', 'LSSGagMza' - Rotating low-speed shaft bending moment at the shaft's strain gage (shaft strain gage located by input ShftGagL)
    % 'YawBrMxp', 'YawBrMyp', 'YawBrMzp' - Nonrotating tower-top / yaw bearing roll moment, ..., Tower-top / yaw bearing yaw moment
    
    absmax_Myc_blade_vals = zeros(3, length(sim_out_list));
    absmax_Myc_cdq_vals = zeros(3, length(sim_out_list));
    mean_Myc_blade_vals = zeros(3, length(sim_out_list));
    mean_Myc_cdq_vals = zeros(3, length(sim_out_list));

    if false
%         time = sim_out_list(1).OutData.time;
%         sim_out_list = [struct];
%         sim_out_list(1).OutData = OutData;
%         sim_out_list(1).ErrorMessage = {};
        time_length = (TMax / DT) + 1 - (50 / DT);
        Myc_blade_vals = zeros(time_length, 3, length(sim_out_list));
        Myc_cdq_vals = zeros(time_length, 3, length(sim_out_list));
        fft_Myc_blade_vals = zeros(time_length, 3, length(sim_out_list));
        fft_Myc_cdq_vals = zeros(time_length, 3, length(sim_out_list));
        
%         fft_
        parfor c = 1:length(sim_out_list)
        %     op_data = [getData(sim_out_list(c).OutData.signals.values, OutList, 'RootMyc1'), ...
        %     getData(sim_out_list(c).OutData.signals.values, OutList, 'RootMyc2'), ...
        %     getData(sim_out_list(c).OutData.signals.values, OutList, 'RootMyc3')];
            
    %         time = sim_out_list(c).OutData.signals.values(:, strmatch('Time', OutList));
            sprintf(['Processing Simulation ' num2str(c)])

            op_data_blade = [sim_out_list(c).OutData.signals.values(50/DT + 1:end, strmatch('RootMyc1', OutList)), ...
                sim_out_list(c).OutData.signals.values(50/DT + 1:end, strmatch('RootMyc2', OutList)), ...
                sim_out_list(c).OutData.signals.values(50/DT + 1:end, strmatch('RootMyc3', OutList))];
            
    %         if ~fields()
            all_data_cdq = mbcTransformOutData(sim_out_list(c).OutData.signals.values, OutList);
    
            op_data_cdq = [all_data_cdq(50/DT + 1:end, strmatch('RootMycC', dqOutList)), ...
                all_data_cdq(50/DT + 1:end, strmatch('RootMycD', dqOutList)), ...
                all_data_cdq(50/DT + 1:end, strmatch('RootMycQ', dqOutList))];
    
    %         absmax_Myc_blade_vals(:, c) = abs(max(op_data_blade, [], 1, 'ComparisonMethod', 'abs'));
    %         absmax_Myc_cdq_vals(:, c) = abs(max(op_data_cdq, [], 1, 'ComparisonMethod', 'abs'));
    %         
    %         mean_Myc_blade_vals(:, c) = mean(op_data_blade, 1);
    %         mean_Myc_cdq_vals(:, c) = mean(op_data_cdq, 1);
            op_data_blade = [op_data_blade; zeros(max(0, time_length - length(op_data_blade)), size(op_data_blade, 2))];
            op_data_cdq = [op_data_cdq; zeros(max(0, time_length - length(op_data_cdq)), size(op_data_cdq, 2))];

            Myc_blade_vals(:, :, c) = op_data_blade;
            Myc_cdq_vals(:, :, c) = op_data_cdq;
            fft_Myc_blade_vals(:, :, c) = fft(op_data_blade);
            fft_Myc_cdq_vals(:, :, c) = fft(op_data_cdq);
        
        end
        figure(1);
        n_uns = 0;
        for c = 1:length(sim_out_list)
            if length(sim_out_list(c).ErrorMessage)
                n_uns = n_uns + 1;
            end
        end

        uns_idx = 0
        for c = 1:length(sim_out_list)
            time = 50/DT:1:(TMax / DT);
            if length(sim_out_list(c).ErrorMessage) == 0
                subplot(n_uns + 1, 1, 1);
                plot(time, Myc_blade_vals(:, 1, c)); 
            else
                c
                uns_idx = uns_idx + 1;
                subplot(n_uns + 1, 1, uns_idx + 1);
                plot(time, Myc_blade_vals(:, 1, c));
            end
        end
        xlabel('Tims [s]'); ylabel('Root Blade 1 Bending Mode [kNm]')
        hold off;

        if strcmp(WIND_TYPE, 'turb')
            save(fullfile(project_dir, 'fft_Myc_blade_vals_turb.mat'), 'fft_Myc_blade_vals', '-v7.3');
            save(fullfile(project_dir, 'fft_Myc_cdq_vals_turb.mat'), 'fft_Myc_cdq_vals', '-v7.3');
        elseif strcmp(WIND_TYPE, 'steady')
            save(fullfile(project_dir, 'fft_Myc_blade_vals_steady.mat'), 'fft_Myc_blade_vals', '-v7.3');
            save(fullfile(project_dir, 'fft_Myc_cdq_vals_steady.mat'), 'fft_Myc_cdq_vals', '-v7.3');
        end
    else
        if strcmp(WIND_TYPE, 'turb')
            load(fullfile(project_dir, 'fft_Myc_blade_vals_turb'));
            load(fullfile(project_dir, 'fft_Myc_cdq_vals_turb'));
        elseif strcmp(WIND_TYPE, 'steady')
            load(fullfile(project_dir, 'fft_Myc_blade_vals_steady'));
            load(fullfile(project_dir, 'fft_Myc_cdq_vals_steady'));
        end
    end

%     figure(2)
%     tiledlayout(2, 2)
%     nexttile
%     semilogx(unique(sweepParameters.cIPC.DQ_Ki_1P), absmax_Myc_blade_vals)
%     legend('1', '2', '3');
%     title('Blade Absmax')
%     nexttile
%     semilogx(unique(sweepParameters.cIPC.DQ_Ki_1P), absmax_Myc_cdq_vals)
%     legend('C', 'D', 'Q');
%     title('CDQ Abxmax')
%     nexttile
%     semilogx(unique(sweepParameters.cIPC.DQ_Ki_1P), mean_Myc_blade_vals)
%     legend('1', '2', '3');
%     title('Blade Mean')
%     nexttile
%     semilogx(unique(sweepParameters.cIPC.DQ_Ki_1P), mean_Myc_cdq_vals)
%     legend('C', 'D', 'Q');
%     title('CDQ Mean')


    % Plot FFT of RootMyb with no IPC, IPC1P and IPC 3P for different
    % gains
    % MANUEL: 3P transform of 3P filtered Mdq with 3P content has frequency (2
    % * omega_3P_Hz ?) 0.532 Hz
    % 1P transform of 1P filtered Mdq has frequency 0.344 Hz
    figure(3);
    tiledlayout(3, 1);
    omega_1P_Hz = Parameters.Turbine.wr_rated * (2*pi/60) * (1/(2*pi));
    omega_2P_Hz = omega_1P_Hz * 2;
    omega_3P_Hz = omega_1P_Hz * 3;

    nexttile; % noIPC case
    sim_indices_noIPC = (sweepParameters.Control.IPCDQ.Enable == 0) .* (sweepParameters.Control.IPC3DQ.Enable == 0);
    [f, P1, Peaks_noIPC] = getFFT(fft_Myc_blade_vals, find(sim_indices_noIPC), DT, omega_1P_Hz);
    plot(f, P1);
    hold on;
    plot([omega_1P_Hz omega_1P_Hz], [0, max(P1)]);
    plot([omega_2P_Hz omega_2P_Hz], [0, max(P1)]);
    plot([omega_3P_Hz omega_3P_Hz], [0, max(P1)]);
    hold off;
    title("Single-Sided Amplitude Spectrum of MyC for no IPC")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    xlim([0, 0.3]);
    ylim([0, 2e4]);

    nexttile; % only 1P IPC controller case
    sim_indices_1P = (sweepParameters.Control.IPCDQ.Enable == 1) .* (sweepParameters.Control.IPC3DQ.Enable == 0);
    [f, P1, Peaks_IPC1] = getFFT(fft_Myc_blade_vals, find(sim_indices_1P), DT, omega_1P_Hz);
%     labels = {};
%     for c = 1:size(P1, 3)
%         plot(f, P1(:, :, c));
%         labels{c} = ['Ki1P ' num2str(sweepParameters.cIPC.DQ_Ki_1P(c))]
%         hold on;
%     end
%     legend(labels);
%     idx = find((sweepParameters.cIPC.D_Ki_1P(find(sim_indices_1P)) == K_1PD_first_s(1)) .* (sweepParameters.cIPC.Q_Ki_1P(find(sim_indices_1P)) == K_1PQ_first_s(1)));
    idx = find((sweepParameters.cIPC.D_Ki_1P(find(sim_indices_1P)) == Ki_1PD_gains(1)) .* (sweepParameters.cIPC.Q_Ki_1P(find(sim_indices_1P)) == Ki_1PQ_gains(1)));
    idx = idx(1);
    plot(f, P1(:, :, idx));
    hold on;
    plot([omega_1P_Hz omega_1P_Hz], [0, max(P1, [], 'all')]);
    plot([omega_2P_Hz omega_2P_Hz], [0, max(P1, [], 'all')]);
    plot([omega_3P_Hz omega_3P_Hz], [0, max(P1, [], 'all')]);
    hold off;
    title("Single-Sided Amplitude Spectrum of MyC for 1P IPC")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    xlim([0, 0.3]);
    ylim([0, 2e4]);
    
    nexttile(3);
    sim_indices_3P = (sweepParameters.Control.IPCDQ.Enable == 0) .* (sweepParameters.Control.IPC3DQ.Enable == 1);
    [f, P1, Peaks_IPC3] = getFFT(fft_Myc_blade_vals, find(sim_indices_3P), DT, omega_1P_Hz);
%     labels = {};
%     for c = 1:size(P1, 3)
%         plot(f, P1(:, :, c));
%         labels{c} = ['Ki1P=' num2str(sweepParameters.cIPC.DQ_Ki_1P(c)) ', Ki3P=' num2str(sweepParameters.cIPC.DQ_Ki_3P(c))]
%         hold on;
%     end
%     legend(labels);
     
    % get index of simulations with first stable 3P gains in D and Q
    idx = find((sweepParameters.cIPC.D_Ki_3P(find(sim_indices_3P)) == Ki_3PD_gains(1)) .* (sweepParameters.cIPC.Q_Ki_3P(find(sim_indices_3P)) == Ki_3PQ_gains(1)));
    plot(f, P1(:, :, idx));
    hold on;
    plot([omega_1P_Hz omega_1P_Hz], [0, max(P1, [], 'all')]);
    plot([omega_2P_Hz omega_2P_Hz], [0, max(P1, [], 'all')]);
    plot([omega_3P_Hz omega_3P_Hz], [0, max(P1, [], 'all')]);
    hold off;
    title("Single-Sided Amplitude Spectrum of MyC for 1P+3P IPC")
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    xlim([0, 0.3]);
    ylim([0, 2e4]);
    legend('FFT', '1P', '2P', '3P');
    
    set(gcf, 'Position', get(0, 'Screensize'));
    savefig(fullfile(project_dir, 'nonlin_figs', 'fft.fig'));

    % Plot FFT values at 1P, 2P and 3P for different DQ gains
    figure(2);
    tiledlayout(2, 2); % Varying D, Q integrator gains along rows; IPC1P, IPC3P activated controllers along columns

    ax1 = nexttile(1); % 1P/2P/3P peaks, Varying D1P gain, IPC1P
    idx = find((sweepParameters.cIPC.Q_Ki_1P(find(sim_indices_1P)) == Ki_1PQ_gains(1)));
    peaks_1P = squeeze(Peaks_IPC1.P1(:, :, idx));
    peaks_2P = squeeze(Peaks_IPC1.P2(:, :, idx));
    peaks_3P = squeeze(Peaks_IPC1.P3(:, :, idx));
    semilogx([Ki_1PD_gains(1) Ki_1PD_gains(end)], [Peaks_noIPC.P1, Peaks_noIPC.P1]);
    hold on;
    semilogx([Ki_1PD_gains(1) Ki_1PD_gains(end)], [Peaks_noIPC.P2, Peaks_noIPC.P2]);
    semilogx([Ki_1PD_gains(1) Ki_1PD_gains(end)], [Peaks_noIPC.P3, Peaks_noIPC.P3]);
    semilogx(Ki_1PD_gains, peaks_1P);
    semilogx(Ki_1PD_gains, peaks_2P);
    semilogx(Ki_1PD_gains, peaks_3P);
    legend('noIPC 1P Peak', 'noIPC 2P Peak', 'noIPC 3P Peak', 'IPC 1P 1P Peak', 'IPC 1P 2P Peak', 'IPC 1P 3P Peak');
    title('1P/2P/3P Amp for IPC 1P vs KiD Gain');
    hold off;

    ax2 = nexttile(2); % 1P/2P/3P peaks, Varying D3P gain, IPC3P
    idx = find((sweepParameters.cIPC.Q_Ki_3P(find(sim_indices_3P)) == Ki_3PQ_gains(1)));
    peaks_1P = squeeze(Peaks_IPC3.P1(:, :, idx));
    peaks_2P = squeeze(Peaks_IPC3.P2(:, :, idx));
    peaks_3P = squeeze(Peaks_IPC3.P3(:, :, idx));
    semilogx([Ki_3PD_gains(1) Ki_3PD_gains(end)], [Peaks_noIPC.P1, Peaks_noIPC.P1]);
    hold on;
    semilogx([Ki_3PD_gains(1) Ki_3PD_gains(end)], [Peaks_noIPC.P2, Peaks_noIPC.P2]);
    semilogx([Ki_3PD_gains(1) Ki_3PD_gains(end)], [Peaks_noIPC.P3, Peaks_noIPC.P3]);
    semilogx(Ki_3PD_gains, peaks_1P);
    semilogx(Ki_3PD_gains, peaks_2P);
    semilogx(Ki_3PD_gains, peaks_3P);
    legend('noIPC 1P Peak', 'noIPC 2P Peak', 'noIPC 3P Peak', 'IPC 3P 1P Peak', 'IPC 3P 2P Peak', 'IPC 3P 3P Peak');
    title('1P/2P/3P Amp for IPC 3P vs KiD Gain');
    hold off;

    ax3 = nexttile(3); % 1P/2P/3P peaks, Varying Q1P gain, IPC1P
    idx = find((sweepParameters.cIPC.D_Ki_1P(find(sim_indices_1P)) == Ki_1PD_gains(1)));
    peaks_1P = squeeze(Peaks_IPC1.P1(:, :, idx));
    peaks_2P = squeeze(Peaks_IPC1.P2(:, :, idx));
    peaks_3P = squeeze(Peaks_IPC1.P3(:, :, idx));
    semilogx([Ki_1PQ_gains(1) Ki_1PQ_gains(end)], [Peaks_noIPC.P1, Peaks_noIPC.P1]);
    hold on;
    semilogx([Ki_1PQ_gains(1) Ki_1PQ_gains(end)], [Peaks_noIPC.P2, Peaks_noIPC.P2]);
    semilogx([Ki_1PQ_gains(1) Ki_1PQ_gains(end)], [Peaks_noIPC.P3, Peaks_noIPC.P3]);
    semilogx(Ki_1PQ_gains, peaks_1P);
    semilogx(Ki_1PQ_gains, peaks_2P);
    semilogx(Ki_1PQ_gains, peaks_3P);
    legend('noIPC 1P Peak', 'noIPC 2P Peak', 'noIPC 3P Peak', 'IPC 1P 1P Peak', 'IPC 1P 2P Peak', 'IPC 1P 3P Peak');
    title('1P/2P/3P Amp for IPC 1P vs KiQ Gain');
    hold off;

    ax4 = nexttile(4); % 1P/2P/3P peaks, Varying Q3P gain, IPC3P
    idx = find((sweepParameters.cIPC.D_Ki_3P(find(sim_indices_3P)) == Ki_3PD_gains(1)));
    peaks_1P = squeeze(Peaks_IPC3.P1(:, :, idx));
    peaks_2P = squeeze(Peaks_IPC3.P2(:, :, idx));
    peaks_3P = squeeze(Peaks_IPC3.P3(:, :, idx));
    semilogx([Ki_3PQ_gains(1) Ki_3PQ_gains(end)], [Peaks_noIPC.P1, Peaks_noIPC.P1]);
    hold on;
    semilogx([Ki_3PQ_gains(1) Ki_3PQ_gains(end)], [Peaks_noIPC.P2, Peaks_noIPC.P2]);
    semilogx([Ki_3PQ_gains(1) Ki_3PQ_gains(end)], [Peaks_noIPC.P3, Peaks_noIPC.P3]);
    semilogx(Ki_3PQ_gains, peaks_1P);
    semilogx(Ki_3PQ_gains, peaks_2P);
    semilogx(Ki_3PQ_gains, peaks_3P);
    legend('noIPC 1P Peak', 'noIPC 2P Peak', 'noIPC 3P Peak', 'IPC 3P 1P Peak', 'IPC 3P 2P Peak', 'IPC 3P 3P Peak');
    title('1P/2P/3P Amp for IPC 3P vs KiD Gain');
    hold off;

%     linkaxes([ax1, ax2], 'y')
%     linkaxes([ax3, ax4], 'y')

    set(gcf, 'Position', get(0, 'Screensize'));
    savefig(fullfile(project_dir, 'nonlin_figs', 'fft_peaks.fig'));


end


%% Get SS values

if COMPUTE_SS_VALS
    [ss_vals, op_absmax] = compute_ss_vals(sim_out_list, OutList, dqOutList, Parameters);
    save(fullfile(fastRunner.FAST_directory, 'ss_vals'), 'ss_vals');
    save(fullfile(fastRunner.FAST_directory, 'op_absmax'), 'op_absmax');
end
% plot blade-pitch, root blade bending moment (RootMyb1, RootMyb2, RootMyb3), 
% tower foreaft (YawBrFxp, TwrBsMyt)
% ss_vals.BlPitch1(WIND_SPEEDS == str2num(case_list(1).InflowWind.HWindSpeed))

%% Investigate Sensitivity of Outputs to Individual BldPitch1 Variation for RUN_OL Nonlinear model
if NONLIN_SENS_ANALYSIS
    Delta_op_nonlin = [];
    % for each case, compute output mean and absmax before and after
    % RampStart
    OutList_op = OutList;
    exc_outdata_fields = {'Time', 'BldPitch1', 'BldPitch2', 'BldPitch3', ...
        'Azimuth', 'RotSpeed', 'Wind1VelX', 'Wind1VelY', 'Wind1VelZ', ...
        '-React', 'Alpha', 'Rt', 'M2N', 'M8N', 'GenTq'};
    inc_outdata_indices = 1:length(OutList);
    for field = exc_outdata_fields
        field_idx = cellfun(@(a) ~isempty(a) && a > 0, strfind(OutList_op, field));
        OutList_op(field_idx) = [];
        inc_outdata_indices(field_idx) = [];
    end
    
    absmax_varnames = {};
    mean_varnames = {};
    abs_absmax_varnames = {};
    abs_mean_varnames = {};

    load(fullfile(fastRunner.FAST_directory, 'op_absmax'));

    for c = 1:length(sim_out_list)
        
        ws = mean(getData(sim_out_list(c).OutData.signals.values, OutList, 'Wind1VelX'));

        nanvals = zeros(length(OutList_op), 1); % ismember(OutList_op, {'GenPwr'});%all(pre_dist.data == 0);

        if RUN_OL_DQ
            pre_dist.data = sim_out_list(c).OutData.signals.dqNormalizedValues(...
                BaselineSteadyState / DT, inc_outdata_indices);

            post_dist.data.d = sim_out_list(c).OutData.signals.dqNormalizedValues(...
            dSteadyState / DT, inc_outdata_indices);
            post_dist.data.q = sim_out_list(c).OutData.signals.dqNormalizedValues(...
                qSteadyState / DT, inc_outdata_indices);

            % compute mean over all rows (time) for each column (output) before/after RampStart    
            pre_dist.mean = mean(pre_dist.data, 1);
            post_dist.mean.d = mean(post_dist.data.d, 1);
            post_dist.mean.q = mean(post_dist.data.q, 1);

            % compute absmax over all rows (time) for each column (output) before/after RampStart
            pre_dist.absmax = abs(max(pre_dist.data, [], 1, 'ComparisonMethod', 'abs'));

            post_dist.absmax.d = abs(max(post_dist.data.d, [], 1, 'ComparisonMethod', 'abs'));
            post_dist.absmax.q = abs(max(post_dist.data.q, [], 1, 'ComparisonMethod', 'abs'));
            
            rel_mean.d = 100 * (post_dist.mean.d - pre_dist.mean) ./ pre_dist.mean;
            rel_mean.q = 100 * (post_dist.mean.q - pre_dist.mean) ./ pre_dist.mean;
    
            rel_absmax.d = 100 * (post_dist.absmax.d - pre_dist.absmax) ./ pre_dist.absmax;
            rel_absmax.q = 100 * (post_dist.absmax.q - pre_dist.absmax) ./ pre_dist.absmax;
    
            abs_rel_mean.d = abs(rel_mean.d);
            abs_rel_mean.q = abs(rel_mean.q);
    
            abs_rel_absmax.d = abs(rel_absmax.d);
            abs_rel_absmax.q = abs(rel_absmax.q);

            absmax_varnames.d{c} = ['Relative Change in AbsMax WS = ' num2str(ws) 'm/s' ' D'];
            mean_varnames.d{c} = ['Relative Change in Mean WS = ' num2str(ws) 'm/s' ' D'];
            absmax_varnames.q{c} = ['Relative Change in AbsMax WS = ' num2str(ws) 'm/s' ' Q'];
            mean_varnames.q{c} = ['Relative Change in Mean WS = ' num2str(ws) 'm/s' ' Q'];
    
            abs_absmax_varnames.d{c} = ['Absolute ' absmax_varnames.d{c}];
            abs_mean_varnames.d{c} = ['Absolute ' mean_varnames.d{c}];
            abs_absmax_varnames.q{c} = ['Absolute ' absmax_varnames.q{c}];
            abs_mean_varnames.q{c} = ['Absolute ' mean_varnames.q{c}];
    
            base_varnames = {mean_varnames.d{c}, absmax_varnames.d{c}, ...
            abs_mean_varnames.d{c}, abs_absmax_varnames.d{c}, ...
            mean_varnames.q{c}, absmax_varnames.q{c}, ...
            abs_mean_varnames.q{c}, abs_absmax_varnames.q{c}};

            Delta_op_nonlin = [Delta_op_nonlin, ...
            table(...
            rel_mean.d(:, ~nanvals)', rel_absmax.d(:, ~nanvals)', ...
            abs_rel_mean.d(:, ~nanvals)', abs_rel_absmax.d(:, ~nanvals)', ...
            rel_mean.q(:, ~nanvals)', rel_absmax.q(:, ~nanvals)', ...
            abs_rel_mean.q(:, ~nanvals)', abs_rel_absmax.q(:, ~nanvals)', ...
            'VariableNames', base_varnames, ...
            'RowNames', OutList_op(~nanvals))];

        elseif RUN_OL_BLADE
            pre_dist.data = sim_out_list(c).OutData.signals.normalizedValues(...
                BaselineSteadyState / DT, inc_outdata_indices);

            post_dist.data = sim_out_list(c).OutData.signals.normalizedValues(...
            b1SteadyState / DT, inc_outdata_indices);

            % compute mean over all rows (time) for each column (output) before/after RampStart    
            pre_dist.mean = mean(pre_dist.data, 1);
            post_dist.mean = mean(post_dist.data, 1);

            % compute absmax over all rows (time) for each column (output) before/after RampStart
            pre_dist.absmax = abs(max(pre_dist.data, [], 1, 'ComparisonMethod', 'abs'));

            post_dist.absmax = abs(max(post_dist.data, [], 1, 'ComparisonMethod', 'abs'));
            
            rel_mean = 100 * (post_dist.mean - pre_dist.mean) ./ pre_dist.mean;
    
            rel_absmax = 100 * (post_dist.absmax - pre_dist.absmax) ./ pre_dist.absmax;
    
            abs_rel_mean = abs(rel_mean);
    
            abs_rel_absmax = abs(rel_absmax);

            absmax_varnames{c} = ['Relative Change in AbsMax WS = ' num2str(ws) 'm/s'];
            mean_varnames{c} = ['Relative Change in Mean WS = ' num2str(ws) 'm/s'];
    
            abs_absmax_varnames{c} = ['Absolute ' absmax_varnames{c}];
            abs_mean_varnames{c} = ['Absolute ' mean_varnames{c}];
    
            base_varnames = {mean_varnames{c}, absmax_varnames{c}, ...
            abs_mean_varnames{c}, abs_absmax_varnames{c}};

            Delta_op_nonlin = [Delta_op_nonlin, ...
            table(...
            rel_mean(:, ~nanvals)', rel_absmax(:, ~nanvals)', ...
            abs_rel_mean(:, ~nanvals)', abs_rel_absmax(:, ~nanvals)', ...
            'VariableNames', base_varnames, ...
            'RowNames', OutList_op(~nanvals))];
        end
        

    end
    
    if RUN_OL_DQ
        for comp = {'d', 'q'}
    
            Delta_op_nonlin = [Delta_op_nonlin, ...
                table(mean(table2array(Delta_op_nonlin(:, ...
                abs_absmax_varnames.(comp{1}))), 2), ...
                'VariableNames', {['Mean Relative Change in AbsMax ' comp{1}]})];
            Delta_op_nonlin = [Delta_op_nonlin, ...
                table(median(table2array(Delta_op_nonlin(:, ...
                abs_absmax_varnames.(comp{1}))), 2), ...
                'VariableNames', {['Median Relative Change in AbsMax ' comp{1}]})];
            Delta_op_nonlin = [Delta_op_nonlin, ...
                table(max(table2array(Delta_op_nonlin(:, ...
                abs_absmax_varnames.(comp{1}))), [], 2), ...
                'VariableNames', {['Max Relative Change in AbsMax ' comp{1}]})];
    
            Delta_op_nonlin = [Delta_op_nonlin, ...
                table(mean(table2array(Delta_op_nonlin(:, ...
                abs_mean_varnames.(comp{1}))), 2), ...
                'VariableNames', {['Mean Relative Change in Mean ' comp{1}]})];
            Delta_op_nonlin = [Delta_op_nonlin, ...
                table(median(table2array(Delta_op_nonlin(:, ...
                abs_mean_varnames.(comp{1}))), 2), ...
                'VariableNames', {['Median Relative Change in Mean ' comp{1}]})];
            Delta_op_nonlin = [Delta_op_nonlin, ...
                table(max(table2array(Delta_op_nonlin(:, ...
                abs_mean_varnames.(comp{1}))), [], 2), ...
                'VariableNames', {['Max Relative Change in Mean ' comp{1}]})];
        end
        Delta_op_nonlin = sortrows(Delta_op_nonlin, 'Max Relative Change in Mean d', 'descend');
        writetable(Delta_op_nonlin, fullfile(fig_dir, 'ol_dq_rel_change.csv'));
    elseif RUN_OL_BLADE
        Delta_op_nonlin = [Delta_op_nonlin, ...
            table(mean(table2array(Delta_op_nonlin(:, ...
            abs_absmax_varnames)), 2), ...
            'VariableNames', {['Mean Relative Change in AbsMax']})];
        Delta_op_nonlin = [Delta_op_nonlin, ...
            table(median(table2array(Delta_op_nonlin(:, ...
            abs_absmax_varnames)), 2), ...
            'VariableNames', {['Median Relative Change in AbsMax']})];
        Delta_op_nonlin = [Delta_op_nonlin, ...
            table(max(table2array(Delta_op_nonlin(:, ...
            abs_absmax_varnames)), [], 2), ...
            'VariableNames', {['Max Relative Change in AbsMax']})];

        Delta_op_nonlin = [Delta_op_nonlin, ...
            table(mean(table2array(Delta_op_nonlin(:, ...
            abs_mean_varnames)), 2), ...
            'VariableNames', {['Mean Relative Change in Mean']})];
        Delta_op_nonlin = [Delta_op_nonlin, ...
            table(median(table2array(Delta_op_nonlin(:, ...
            abs_mean_varnames)), 2), ...
            'VariableNames', {['Median Relative Change in Mean']})];
        Delta_op_nonlin = [Delta_op_nonlin, ...
            table(max(table2array(Delta_op_nonlin(:, ...
            abs_mean_varnames)), [], 2), ...
            'VariableNames', {['Max Relative Change in Mean']})];

        Delta_op_nonlin = sortrows(Delta_op_nonlin, 'Max Relative Change in Mean', 'descend');
        writetable(Delta_op_nonlin, fullfile(fig_dir, 'ol_blade_rel_change.csv'));
    end

    % plot barchart of relative change in mean/max from before RampStart to
    % after for 'big movers'
    figure(1);
    tiledlayout(2, 1);
    title('Relative Change in Mean/Absolute Maximum of Outputs over Disturbance[%]')
    test_case_idx = find(WIND_SPEEDS == 14);
    ws_suffix = num2str(WIND_SPEEDS(test_case_idx));
    axes = [];
    
    if RUN_OL_DQ
        % barchart data on rel change of DQ outputs
        mean_data = table2array(Delta_op_nonlin(:, ...
            {mean_varnames.d{test_case_idx} mean_varnames.q{test_case_idx}}));
        absmax_data = table2array(Delta_op_nonlin(:, ...
            {absmax_varnames.d{test_case_idx} absmax_varnames.q{test_case_idx}}));

        % remove outputs which have nan values in mean or absmax, d or q
        nan_idx = all((~isnan(mean_data)) & (~isnan(absmax_data)), 2);
        OutList_op = OutList_op(nan_idx);
        mean_data = mean_data(nan_idx, :);
        absmax_data = absmax_data(nan_idx, :);

        % only plot outliers
        n_stds = 2;
        outlier_idx = ...
            (any(mean_data > mean(mean_data, 1) + (n_stds * std(mean_data, 1)), 2)) ...
            | (any(mean_data < mean(mean_data, 1) - (n_stds * std(mean_data, 1)), 2)) ...
            | (any(absmax_data > mean(absmax_data, 1) + (n_stds * std(absmax_data, 1)), 2)) ...
            | (any(absmax_data < mean(absmax_data, 1) - (n_stds * std(absmax_data, 1)), 2));
        OutList_op = OutList_op(outlier_idx);
        mean_data = mean_data(outlier_idx, :);
        absmax_data = absmax_data(outlier_idx, :);
        
        % plot mean before/after RampStart
        axes = [axes, nexttile(1)];
        bar(mean_data);
        title('Mean');
        xticks(1:length(OutList_op));
        xticklabels(OutList_op);
        legend('d', 'q');
    
        % plot max before/after RampStart
        axes = [axes, nexttile(2)];
        bar(absmax_data);
        title('Absolute Maximum');
        xticks(1:length(OutList_op));
        xticklabels(OutList_op);
        legend('d', 'q');
        linkaxes(axes, 'x');

        savefig(fullfile(fig_dir, ['ol_dq_rel_change_' ws_suffix '.fig']));
        saveas(gcf, fullfile(fig_dir, ['ol_dq_rel_change_' ws_suffix '.png']));
    elseif RUN_OL_BLADE
        
        % barchart data on rel change of BLADE outputs
        mean_data = table2array(Delta_op_nonlin(:, ...
            {mean_varnames{test_case_idx}}));
        absmax_data = table2array(Delta_op_nonlin(:, ...
            {absmax_varnames{test_case_idx}}));

        % remove outputs which have nan values in mean or absmax, d or q
        nan_idx = all((~isnan(mean_data)) & (~isnan(absmax_data)), 2);
        OutList_op = OutList_op(nan_idx);
        mean_data = mean_data(nan_idx, :);
        absmax_data = absmax_data(nan_idx, :);

        % only plot outliers
        n_stds = 0.75;
        outlier_idx = ...
            (any(mean_data > mean(mean_data, 1) + (n_stds * std(mean_data, 1)), 2)) ...
            | (any(mean_data < mean(mean_data, 1) - (n_stds * std(mean_data, 1)), 2)) ...
            | (any(absmax_data > mean(absmax_data, 1) + (n_stds * std(absmax_data, 1)), 2)) ...
            | (any(absmax_data < mean(absmax_data, 1) - (n_stds * std(absmax_data, 1)), 2));
        OutList_op = OutList_op(outlier_idx);
        mean_data = mean_data(outlier_idx, :);
        absmax_data = absmax_data(outlier_idx, :);

        % plot mean before/after RampStart
        axes = [axes, nexttile(1)];
        bar(mean_data);
        title('Mean');
        xticks(1:length(OutList_op));
        xticklabels(OutList_op);
    
        % plot max before/after RampStart
        axes = [axes, nexttile(2)];
        bar(absmax_data);
        title('Absolute Maximum');
        xticks(1:length(OutList_op));
        xticklabels(OutList_op);
        linkaxes(axes, 'x');

        savefig(fullfile(fig_dir, ['ol_blade_rel_change_' ws_suffix '.fig']));
        saveas(gcf, fullfile(fig_dir, ['ol_blade_rel_change_' ws_suffix '.png']))
    end

    
%     n_output_plots = 0;
%     for o = 1:length(OutList)
%         op = OutList{o};
%         if ~strcmp(op, 'Time') && ~strcmp(op(end), '2') && ~strcmp(op(end), '3')
%             n_output_plots = n_output_plots + 1;
%         end
%     end
    
    % Blade root bending moment: out of plane, RootMyc1,2,3
    % Shaft bending moment (My), LSSTipMya. LSSTipMys
    % Yaw bearing yaw moment (Mz&My), YawBrMyp, YawBrMzp, YawBrMyn, YawBrMzn
    % Deflection: 
    % Tower Base
    if RUN_OL_DQ
        OutList_plot = {
            'BldPitchC', 'BldPitchD', 'BldPitchQ', ...
            'RootMycC', 'RootMycD', 'RootMycqQ', ... % Blade 3 out-of-plane moment (i.e., the moment caused by out-of-plane forces) at the blade root
            'OoPDeflC', 'OoPDeflD', 'OoPDeflQ'};
%             'GenSpeed', ...
%             'TwrBsMyt', ... % Tower base pitching (or fore-aft) moment (i.e., the moment caused by fore-aft forces)
%             'YawBrMyp', ... % Nonrotating tower-top / yaw bearing pitch moment
%             'YawBrMzp', ... % Tower-top / yaw bearing yaw moment
%             'LSSGagMya' % Nonrotating low-speed shaft bending moment at the shaft tip (teeter pin for 2-blader, apex of rotation for 3-blader)
%         };
%             'IPDeflC', 'IPDeflD', 'IPDeflQ'};
        n_output_plots = 3;

        time_data = getData(sim_out_list(1).OutData.signals.values, OutList, 'Time');
        
        figure(2);
        tiledlayout(n_output_plots, 1);
        axs = [];
       
        
        for op_label = OutList_plot
            if strcmp(op_label{1}(end), 'C')
                axs = [axs, nexttile];
                
                cdq_op_labels = cellfun(@(q) [op_label{1}(1:end-1) q], ...
                    {'C', 'D', 'Q'}, 'UniformOutput', false);
                for cdq_label = cdq_op_labels
                    plot(time_data, getData(sim_out_list(test_case_idx).OutData.signals.dqValues, dqOutList, cdq_label{1})); hold on;
                end
     
                plot([dRampStart, dRampStart], ylim, 'k--');
                plot([dRampStop, dRampStop], ylim, 'k--');
                plot([qRampStart, qRampStart], ylim, 'k--');
                plot([qRampStop, qRampStop], ylim, 'k--');
                hold off;
                ylabel(op_label{1}(1:end-1));
                legend('c', 'd', 'q');
            elseif ~strcmp(op_label{1}(end), 'D') && ~strcmp(op_label{1}(end), 'Q')
                axs = [axs, nexttile];
                plot(time_data, getData(sim_out_list(test_case_idx).OutData.signals.dqValues, dqOutList, op_label{1}));
                hold on;
                plot([dRampStart, dRampStart], ylim, 'k--');
                plot([dRampStop, dRampStop], ylim, 'k--');
                plot([qRampStart, qRampStart], ylim, 'k--');
                plot([qRampStop, qRampStop], ylim, 'k--');
                hold off;
                ylabel(op_label{1});
            end
        end
        linkaxes(axs, 'x');
        xlabel('Time [s]');
        xlim([0, max(time_data, [], 1)]);
    
        savefig(fullfile(fig_dir, ['ol_dq_ts_' ws_suffix '.fig']));
        saveas(gcf, fullfile(fig_dir, ['ol_dq_ts_' ws_suffix '.png']));
    elseif RUN_OL_BLADE
        OutList_plot = {
            'BldPitch1', 'BldPitch2', 'BldPitch3', ...
            'RootMyc1', 'RootMyc2', 'RootMycq3', ... % Blade 3 out-of-plane moment (i.e., the moment caused by out-of-plane forces) at the blade root
            'OoPDefl1', 'OoPDefl2', 'OoPDefl3'};
            %'GenSpeed', ...
            %'TwrBsMyt', ... % Tower base pitching (or fore-aft) moment (i.e., the moment caused by fore-aft forces)
            %'YawBrMyp', ... % Nonrotating tower-top / yaw bearing pitch moment
            %'YawBrMzp', ... % Tower-top / yaw bearing yaw moment
            %'LSSGagMya' % Nonrotating low-speed shaft bending moment at the shaft tip (teeter pin for 2-blader, apex of rotation for 3-blader)
%             };
        %'TwrBsFxt', TwrBsMxt, 
    %         'OoPDefl1', 'OoPDefl2', 'OoPDefl3', ...
    %         'IPDefl1', 'IPDefl2', 'IPDefl3'};
        n_output_plots = 3;

        time_data = getData(sim_out_list(1).OutData.signals.values, OutList, 'Time');
        
        figure(2);
        tiledlayout(n_output_plots, 1);
        axs = [];
        
        for op_label = OutList_plot
            if strcmp(op_label{1}(end), '1')
                axs = [axs, nexttile];
                
                blade_op_labels = cellfun(@(q) [op_label{1}(1:end-1) q], ...
                    {'1', '2', '3'}, 'UniformOutput', false);
                for blade_label = blade_op_labels
                    plot(time_data, getData(sim_out_list(test_case_idx).OutData.signals.values, OutList, blade_label{1})); hold on;
                end
                plot([b1RampStart, b1RampStart], ylim, 'k--');
                plot([b1RampStop, b1RampStop], ylim, 'k--');
                hold off;
                ylabel(op_label{1}(1:end-1));
                legend('1', '2', '3');
            elseif ~strcmp(op_label{1}(end), '2') && ~strcmp(op_label{1}(end), '3')
                axs = [axs, nexttile];
                plot(time_data, getData(sim_out_list(test_case_idx).OutData.signals.values, OutList, op_label{1}));
                hold on;
                plot([b1RampStart, b1RampStart], ylim, 'k--');
                plot([b1RampStop, b1RampStop], ylim, 'k--');
                hold off;
                ylabel(op_label{1});
            end
        end
        linkaxes(axs, 'x');
        xlabel('Time [s]')
        xlim([0, max(time_data, [], 1)]);
    
        savefig(fullfile(fig_dir, ['ol_blade_ts_' ws_suffix '.fig']));
        saveas(gcf, fullfile(fig_dir, ['ol_blade_ts_' ws_suffix '.png']));
    end

end

%% Investigate Sensitivity of Outputs to Individual BldPitch1 Variation for Linear model

if LIN_SENS_ANALYSIS

    % Import and MBC-transform linear models
    all_linfiles = dir(lin_models_dir);

    input_arr = {'Blade C pitch command', ...
            'Blade D pitch command', ...
            'Blade Q pitch command'};
    h2_norms = table([], [], [], 'VariableNames', input_arr, 'RowNames', {}); % arrayfun(@(ws) num2str(ws), WIND_SPEEDS, 'UniformOutput', false));
    hinf_norms = table([], [], [], 'VariableNames', input_arr, 'RowNames', {});
    plotting_idx = find(WIND_SPEEDS == 14);

    load(fullfile(fastRunner.FAST_directory, 'op_absmax'));

    parfor ws_idx = 1:length(WIND_SPEEDS)

        ws = WIND_SPEEDS(ws_idx);

%         if ws_idx ~= plotting_idx
%             continue
%         end

        sprintf(['Processing linfiles for ' num2str(ws) 'm/s'])

        linfile_prefix = ['lin_' num2str(ws_idx - 1, '%02u') '.'];
        linfiles = {};
        linfile_idx = 1;
        for f_idx = 1:length(all_linfiles)
            filepath = fullfile(lin_models_dir, all_linfiles(f_idx).name);
            if any(strfind(all_linfiles(f_idx).name, linfile_prefix)) && isfile(filepath)
                linfiles{linfile_idx} = filepath;
                linfile_idx = linfile_idx + 1;
            end
        end
        
        % 1,2,3 -> c(0), d(cos), q(sin)
        [MBC, matData, FAST_linData, VTK] = fx_mbc3(linfiles);
        [MBC_3P, matData_3P, FAST_linData_3P, VTK_3P] = fx_mbc3_3P(linfiles);

        
        if matData.WindSpeed(1) ~= ws
            error('*** Linfiles do not correspond to provided WIND_SPEEDS.');
        end

        % Rename labels
        MBC.DescStates = transformLabels(MBC.DescStates);
        MBC_3P.DescStates = transformLabels(MBC_3P.DescStates);
        matData.DescCntrlInpt = transformLabels(matData.DescCntrlInpt);
        matData_3P.DescCntrlInpt = transformLabels(matData_3P.DescCntrlInpt);
        matData.DescOutput = transformLabels(matData.DescOutput);
        matData_3P.DescOutput = transformLabels(matData_3P.DescOutput);

        sys = ss(MBC.AvgA, MBC.AvgB, MBC.AvgC, MBC.AvgD, ...
            'StateName', MBC.DescStates, ...
            'InputName', matData.DescCntrlInpt, ...
            'OutputName', matData.DescOutput);
        sys_3P = ss(MBC_3P.AvgA, MBC_3P.AvgB, MBC_3P.AvgC, MBC_3P.AvgD, ...
            'StateName', MBC_3P.DescStates, ...
            'InputName', matData_3P.DescCntrlInpt, ...
            'OutputName', matData_3P.DescOutput);

%         sum(sys_3P.A - sys.A, 'all')
%         sum(sys_3P.B - sys.B, 'all')
%         sum(sys_3P.C - sys.C, 'all')
%         sum(sys_3P.D - sys.D, 'all')
        % generate indices of OutList outputs corresponding to sys
        % outputs
        if ~exist('dqOutList2Lin_idx', 'var')
            dqOutList2Lin_idx = -1 * ones(length(dqOutList), 1);
            for lin_op_idx = 1:length(matData.DescOutput)
                for nonlin_op_idx = 1:length(dqOutList)
                    if contains(matData.DescOutput(lin_op_idx), dqOutList(nonlin_op_idx))
                        dqOutList2Lin_idx(lin_op_idx) = nonlin_op_idx;%find(matData.DescOutput(lin_op_idx), OutList);
                    end
                end
            end
        end

        % remove outputs in sys not present in OutList
        y_drop_idx = dqOutList2Lin_idx == -1;
        sys = sys(~y_drop_idx, :);
        sys_3P = sys_3P(~y_drop_idx, :);
        matData.DescOutput = matData.DescOutput(~y_drop_idx);
        matData_3P.DescOutput = matData_3P.DescOutput(~y_drop_idx);
        dqOutList2Lin_idx = dqOutList2Lin_idx(~y_drop_idx);
        
        % normalize TODO this does not work, different wind speeds ?? Ask
        % Manuel
%         sys = sys * diag(1 ./ table2array(op_absmax.dq(ws_idx, dqOutList2Lin_idx)))';
%         if ws_idx == plotting_idx
%             bodemag(sys(ismember(matData.DescOutput, ...
%                 {'ED RootMycD, (kN-m)', 'ED RootMycQ, (kN-m)'}), ...
%                 ismember(matData.DescCntrlInpt, ...
%                 {'ED Blade D pitch command, rad', 'ED Blade Q pitch command, rad'})), ...
%                 sys_3P(ismember(matData_3P.DescOutput, ...
%                 {'ED RootMycD, (kN-m)', 'ED RootMycQ, (kN-m)'}), ...
%                 ismember(matData_3P.DescCntrlInpt, ...
%                 {'ED Blade D pitch command, rad', 'ED Blade Q pitch command, rad'})))
%             legend('1P', '3P')
%         end

        sys = sysclean(sys); % zero almost zero elements
        sys_3P = sysclean(sys_3P); % zero almost zero elements
        
        az_idx = ismember(MBC.DescStates, 'ED Variable speed generator DOF (internal DOF index = DOF_GeAz), rad');
        sys = modred(sys, az_idx, 'truncate'); % remove azimuth state
        az_idx_3P = ismember(MBC.DescStates, 'ED Variable speed generator DOF (internal DOF index = DOF_GeAz), rad');
        sys_3P = modred(sys_3P, az_idx_3P, 'truncate'); % remove azimuth state
        
        exc_ops = {'IfW Wind1VelX, (m/s)', 'IfW Wind1VelY, (m/s)', ...
            'SrvD GenPwr, (kW)', 'SrvD GenTq, (kN-m)', 'ED RotTorq, (kN-m)', ...
            'ED BldPitchC, (deg)', 'ED BldPitchD, (deg)', 'ED BldPitchQ, (deg)'};
        
        u_zero_idx = sum(abs(sys.B), 1) == 0; % remove cols from B corresponding to zero sum, rows from C, D corresponding to zero
        y_zero_idx = (sum(abs(sys.C), 2) == 0) .* (sum(abs(sys.D), 2) == 0);
        y_zero_idx(ismember(matData.DescOutput, exc_ops)) = 1;
        sys = sys(~y_zero_idx, ~u_zero_idx);
        op_arr = matData.DescOutput(~y_zero_idx);
        ip_arr = matData.DescCntrlInpt(~u_zero_idx);

        u_zero_idx_3P = sum(abs(sys_3P.B), 1) == 0; % remove cols from B corresponding to zero sum, rows from C, D corresponding to zero
        y_zero_idx_3P = (sum(abs(sys_3P.C), 2) == 0) .* (sum(abs(sys_3P.D), 2) == 0); % QUESTION MANUEL is this correct, should it be AND rather than OR
        y_zero_idx_3P(ismember(matData_3P.DescOutput, exc_ops)) = 1;
        sys_3P = sys_3P(~y_zero_idx_3P, ~u_zero_idx_3P);
        op_arr_3P = matData_3P.DescOutput(~y_zero_idx_3P); % QUESTIION MANUEL this removes RootMyCDQ
        ip_arr_3P = matData_3P.DescCntrlInpt(~u_zero_idx_3P);

        % spy(sys.D) % nonzero elements

        state_arr = MBC.DescStates(~az_idx);
        short_op_arr = cellfun(@(op) simplifyOpName(op), op_arr', 'UniformOutput', false);
        short_ip_arr = cellfun(@(ip) simplifyOpName(ip), ip_arr', 'UniformOutput', false);
        short_state_arr = cellfun(@(state) simplifyOpName(state), state_arr', 'UniformOutput', false);

        state_arr_3P = MBC_3P.DescStates(~az_idx_3P);
        short_op_arr_3P = cellfun(@(op) simplifyOpName(op), op_arr_3P', 'UniformOutput', false);
        short_ip_arr_3P = cellfun(@(ip) simplifyOpName(ip), ip_arr_3P', 'UniformOutput', false);
        short_state_arr_3P = cellfun(@(state) simplifyOpName(state), state_arr_3P', 'UniformOutput', false);
        
        % remove any outputs not output in nonlinear model
        y_zero_idx = zeros(length(op_arr), 1);
        for op_idx = 1:length(short_op_arr)
            if sum(ismember(dqOutList, short_op_arr(op_idx))) == 0
                y_zero_idx(op_idx) = 1;
            end
        end
        sys = sys(~y_zero_idx, :);
        op_arr = op_arr(~y_zero_idx);
        short_op_arr = short_op_arr(~y_zero_idx);

        y_zero_idx_3P = zeros(length(op_arr_3P), 1);
        for op_idx = 1:length(short_op_arr_3P)
            if sum(ismember(dqOutList, short_op_arr_3P(op_idx))) == 0
                y_zero_idx_3P(op_idx) = 1;
            end
        end

        sys_3P = sys_3P(~y_zero_idx_3P, :);
        op_arr_3P = op_arr_3P(~y_zero_idx_3P);
        short_op_arr_3P = short_op_arr_3P(~y_zero_idx_3P);

        sys.StateName = short_state_arr';
        sys.InputName = short_ip_arr';
        sys.OutputName = short_op_arr';
        sys_3P.StateName = short_state_arr_3P';
        sys_3P.InputName = short_ip_arr_3P';
        sys_3P.OutputName = short_op_arr_3P';

        sys = sys(short_op_arr, short_ip_arr);
        sys_3P = sys_3P(short_op_arr_3P, short_ip_arr_3P);
        
        % remove 2nd bending modes: 
        % Options 
        % a) modred, truncate (removing rows/cols in A, B matrix) 
        % b) modred, matchdc ie residualization (assume states to remove are
        % settled, Enforce matching DC gains)
        second_bm_states = contains(sys.StateName, '2nd');
        sys_red = modred(sys, second_bm_states, 'MatchDC');
        second_bm_states_3P = contains(sys_3P.StateName, '2nd');
        sys_3P_red = modred(sys_3P, second_bm_states_3P, 'MatchDC');
        
        % LPF
        max_freq = max(abs(eig(sys))) * 2;
        lpf = tf([max_freq], [1 max_freq]);
        sys = lpf * sys; % zero out D matrix

        max_freq = max(abs(eig(sys_red))) * 2;
        lpf = tf([max_freq], [1 max_freq]);
        sys_red = lpf * sys_red; % zero out D matrix

        max_freq = max(abs(eig(sys_3P))) * 2;
        lpf = tf([max_freq], [1 max_freq]);
        sys_3P = lpf * sys_3P; % zero out D matrix

        max_freq = max(abs(eig(sys_3P_red))) * 2;
        lpf = tf([max_freq], [1 max_freq]);
        sys_3P_red = lpf * sys_3P_red; % zero out D matrix
        
        xop = matData.Avgxop(~az_idx);
        xop_3P = matData_3P.Avgxop(~az_idx);
        % additional state added for each output due to low-pass filter
        xop_red = [zeros(size(sys_red.OutputName, 1), 1); xop(~second_bm_states)];
        xop_3P_red = [zeros(size(sys_3P_red.OutputName, 1), 1); xop_3P(~second_bm_states_3P)];
%         BldPitch2Op_filt2 = BldPitch2Op; % Manuel: Don't do this, D is
%         relevant at lower frequencies as it adds zeros
%         BldPitch2Op_filt2.D = 0;
        
        xop_arr(:, ws_idx) = xop;
        xop_red_arr(:, ws_idx) = xop_red;
        sys_arr(:, :, ws_idx) = sys;
        sys_red_arr(:, :, ws_idx) = sys_red;

        xop_3P_arr(:, ws_idx) = xop_3P;
        xop_3P_red_arr(:, ws_idx) = xop_3P_red;
        sys_3P_arr(:, :, ws_idx) = sys_3P;
        sys_3P_red_arr(:, :, ws_idx) = sys_3P_red;

        if 0

        h2_norm = singleNorm(sys);
        hinf_norm = singleNorm(sys, Inf);
        
        % compute H2 norm for each channel seperately
        h2_norm = h2_norm ./ table2array(op_absmax.dq(num2str(WIND_SPEEDS(ws_idx)), short_op_arr))';
        h2_norms((ws_idx - 1) * length(op_arr) + 1:(ws_idx) * length(op_arr), :) = table(...
            h2_norm(:, 1), h2_norm(:, 2), h2_norm(:, 3), ...
            'VariableNames', input_arr', ...
            'RowNames', cellfun(@(op) ...
            [op '_' replace(num2str(ws), '.', '_')], ...
            op_arr', 'UniformOutput', false));
        
        % compute Hinf norm
        hinf_norm = hinf_norm ./ table2array(op_absmax.dq(num2str(WIND_SPEEDS(ws_idx)), short_op_arr))';
        hinf_norms((ws_idx - 1) * length(op_arr) + 1:(ws_idx) * length(op_arr), :) = table( ...
            hinf_norm(:, 1), hinf_norm(:, 2), hinf_norm(:, 3), ...
            'VariableNames', input_arr, ...
            'RowNames', cellfun(@(op) ...
            [op '_' replace(num2str(ws), '.', '_')], ...
            op_arr', 'UniformOutput', false));
        
        end

        if false && ws_idx == plotting_idx
            
            plotting_op_arr = {
                'RootMycC', ...
                'RootMycD', ...
                'RootMycQ'
                };

            ws_suffix = replace(num2str(ws), '.', '_');
            
            figure(3);
            spy(sys);
            savefig(fullfile(fig_dir, ['sys_spy_' ws_suffix '.fig']));
            saveas(gcf, fullfile(fig_dir, ['sys_spy_' ws_suffix '.png']));

%             spy(h2_norm);
%             savefig(fullfile(fig_dir, ['h2_norm_spy_' ws_suffix '.fig']));
%             saveas(gcf, fullfile(fig_dir, ['h2_norm_spy_' ws_suffix '.png']));
% 
%             spy(hinf_norm);
%             savefig(fullfile(fig_dir, ['hinf_norm_spy_' ws_suffix '.fig']));
%             saveas(gcf, fullfile(fig_dir, ['hinf_norm_spy_' ws_suffix '.png']));
            
            % TODO 
            % tighter freq range
            figure(4);
            p = bodeoptions;
            p.MagUnits = 'abs';
            p.Title.String = 'Frequency Response: Blade-Pitch to Blade Root Bending Moment';
            p.Title.FontSize = 14;
            p.XLabel.FontSize = 14;
            p.YLabel.FontSize = 14;
            p.TickLabel.FontSize = 14;
            p.InputLabels.FontSize = 14;
            p.OutputLabels.FontSize = 14;
            %p.XLim = [];
            % p.XLimMode = 'manual';
            p.PhaseVisible = 'off';
%             bodeplot(BldPitch2Op(plotting_op_arr, input_arr), ...
%                 BldPitch2Op_filt(plotting_op_arr, input_arr), p);
            bodeplot(sys(plotting_op_arr, input_arr), p);
%             legend('Unfiltered', 'LPF Filtered');
%             gcf.se
            savefig(fullfile(fig_dir, ['bodemag_' ws_suffix '.fig']));
            saveas(gcf, fullfile(fig_dir, ['bodemag_' ws_suffix '.png']));
            
            figure(5);
            h2_norm_data = [h2_norm(:, 2) h2_norm(:, 3)];
            hinf_norm_data = [hinf_norm(:, 2) hinf_norm(:, 3)];

            % only plot outliers
            n_stds = 2;
            outlier_idx = ...
                (any(h2_norm_data > mean(h2_norm_data, 1) + (n_stds * std(h2_norm_data, 1)), 2)) ...
                | (any(h2_norm_data < mean(h2_norm_data, 1) - (n_stds * std(h2_norm_data, 1)), 2)) ...
                | (any(hinf_norm_data > mean(hinf_norm_data, 1) + (n_stds * std(hinf_norm_data, 1)), 2)) ...
                | (any(hinf_norm_data < mean(hinf_norm_data, 1) - (n_stds * std(hinf_norm_data, 1)), 2));
            
            h2_norm_data = h2_norm_data(outlier_idx, :);
            hinf_norm_data = hinf_norm_data(outlier_idx, :);

            tiledlayout(2, 2); % h2 and hinf for each bldpitch c/d/q
%             axes(1) = nexttile(1); bar(h2_norm(:, 1)); title('H2 Norm BldPitchC'); xticks(1:length(short_op_arr)); xticklabels(short_op_arr);
            axes(1) = nexttile(1); bar(h2_norm_data(:, 1)); title('H2 Norm BldPitchD'); xticks(1:length(short_op_arr)); xticklabels(short_op_arr);
            axes(2) = nexttile(2); bar(h2_norm_data(:, 2)); title('H2 Norm BldPitchQ'); xticks(1:length(short_op_arr)); xticklabels(short_op_arr);
%             axes(4) = nexttile(4); bar(hinf_norm(:, 1)); title('Hinf Norm BldPitchC'); xticks(1:length(short_op_arr)); xticklabels(short_op_arr);
            axes(3) = nexttile(3); bar(hinf_norm_data(:, 1)); title('Hinf Norm BldPitchD'); xticks(1:length(short_op_arr)); xticklabels(short_op_arr);
            axes(4) = nexttile(4); bar(hinf_norm_data(:, 2)); title('Hinf Norm BldPitchQ'); xticks(1:length(short_op_arr)); xticklabels(short_op_arr);
            linkaxes([axes(1), axes(2)], 'xy');
            linkaxes([axes(3), axes(4)], 'xy');
            savefig(fullfile(fig_dir, ['Hnorm_' ws_suffix '.fig']));
            saveas(gcf, fullfile(fig_dir, ['Hnorm_' ws_suffix '.png']));

        end
    end
    
    save(fullfile(plant_setup_dir, 'xop_arr'), 'xop_arr');
    save(fullfile(plant_setup_dir, 'sys_arr'), 'sys_arr');
    save(fullfile(plant_setup_dir, 'xop_red_arr'), 'xop_red_arr');
    save(fullfile(plant_setup_dir, 'sys_red_arr'), 'sys_red_arr');
    save(fullfile(plant_setup_dir, 'xop_3P_arr'), 'xop_3P_arr');
    save(fullfile(plant_setup_dir, 'sys_3P_arr'), 'sys_3P_arr');
    save(fullfile(plant_setup_dir, 'xop_3P_red_arr'), 'xop_3P_red_arr');
    save(fullfile(plant_setup_dir, 'sys_3P_red_arr'), 'sys_3P_red_arr');

    save(fullfile(plant_setup_dir, 'h2_norms'), 'h2_norms');
    save(fullfile(plant_setup_dir, 'hinf_norms'), 'hinf_norms');
    
%     writematrix(sys_arr, fullfile(fig_dir, 'sys_arr.csv'));
    writetable(h2_norms, fullfile(fig_dir, 'h2_norms.csv'));
    writetable(hinf_norms, fullfile(fig_dir, 'hinf_norms.csv'));
% else
%     load(fullfile(project_dir, 'xop_arr'));
%     load(fullfile(project_dir, 'sys_arr'));
%     load(fullfile(project_dir, 'h2_norms'));
%     load(fullfile(project_dir, 'hinf_norms'));
end

%% Hinf Controller
% TODO
% 1) choose e1, e2, d1, d2: 
% e1 = d/q BldPitch command
% e2 = d/q RootMyc measurement/error
% d1 = d/q BldPitch disturbance
% d2 = d/q RootMyc noise
% 
% 2) normalize units w/ maximum time-domain magnitude of signal
% 
% 3) Choose TuningGoal as Margins to set Weighting Matrices
% 
% 4) Define controller K(s) as blkdiag(PI, PI)
% 
% 5) Tune and test w/ single wind speed = 14 m/s
%
% ...
%
% -) Tune and test for individual wind speeds, all wind speeds
% -) Use MIMO PI (vs 2 along diagonal) to consider dq coupling
% -) Consider additional outputs in e2

if SETUP_LIN_HINF
    % define plant, choosing only d/q blade-pitches as inputs and d/q
    % RootMyc as outputs
    ws = 14;
    input_arr = {
        'Blade D pitch command', ...
          'Blade Q pitch command'
            };
    min_op_arr = {
            'RootMycD', ...
            'RootMycQ'
            };
    G = sys_arr(:, :, WIND_SPEEDS == ws); % plant at 14m/s
    G = G(min_op_arr, input_arr);

    % define filter TODO, check simulink
    p3_filter = blkdiag(tunableTF('P3_filt_D', 0, 1), tunableTF('P3_filt_Q', 0, 1));
    p6_filter = blkdiag(tunableTF('P6_filt_D', 0, 1), tunableTF('P6_filt_Q', 0, 1));
%     lp_filter = tunableTF(); % first or second order, before 6p
    
    % define controller
    PI_D = tunablePID('PI_D', 'pi'); % tunable PI for d components
    PI_D.Kp.Value = Parameters.cIPC.DQ_Kp_1P;
    PI_D.Ki.Value = Parameters.cIPC.DQ_Ki_1P;
    PI_Q = tunablePID('PI_Q', 'pi'); % tunable PI for q components
    PI_Q.Kp.Value = Parameters.cIPC.DQ_Kp_1P;
    PI_Q.Ki.Value = Parameters.cIPC.DQ_Ki_1P;
    PI_DQ = blkdiag(PI_D, PI_Q);

    % define analysis points
    output_ap = blkdiag(AnalysisPoint('RootMycD'), AnalysisPoint('RootMycQ'));
    input_ap = blkdiag(AnalysisPoint('BldPitchD'), AnalysisPoint('BldPitchQ'));

    % build closed-loop transfer function
    T0 = feedback(output_ap * G * input_ap * p3_filter * p6_filter * PI_DQ, eye(2));
    T0.InputName = {'r2', 'r3'};
    T0.OutputName = {'RootMycD', 'RootMycQ'};

    % Define Tuning Goals
    % TODO Start with 1, 2 etc Tuning Goals, remove redundant goals, add
    % complexity
    % margins are always important for robustness, location of complex gain
    % matters for MIMO systems, useful to have at most critical location
    % one margin location at a time - loop-at-a-time margin
    Req1 = TuningGoal.Margins({'RootMycD', 'RootMycQ'}, 2, 10); % 6dB of Gain Margin and 45deg of PM
    Req2 = TuningGoal.Margins({'BldPitchD', 'BldPitchQ'}, 2, 10);
    Req3 = TuningGoal.Tracking({'BldPitchD', 'BldPitchQ'}, {'RootMycD', 'RootMycQ'}, 1); 
    
    % TODO feed with filtered signal, set inputs to BldPitch and outputs to
    % bending moments

    rng('default')
    Options = systuneOptions('RandomStart', 3   );
    [T, fSoft] = systune(T0, [Req1, Req2, Req3], Options);

    getBlockValue(T,'PI_D')
    getBlockValue(T,'PI_Q')
    
    showTunable(T0)
    showTunable(T)  % tuned values of all tunable elements
    
    figure(1);
%     tiledlayout(1, 2);
%     nexttile;
    viewGoal([Req1, Req2, Req3],T0); title('initial'); 
%     nexttile;
    figure(2);
    viewGoal([Req1, Req2, Req3],T); title('tuned');
%     legend('initial', 'tuned')
end

%% TODO plug in LPV model instead of OpenFAST block
if SETUP_NONLIN_HINF
    % open SL controller
    open_system(FAST_SimulinkModel);

    % create slTuner instance and list tunable blocks
    opts = slTunerOptions('RateConversionMethod', 'tustin');
    ipc_mdl = 'Baseline Controller/Cyclic Pitch controller/1P Cyclic Pitch Controller1/';
    ST0 = slTuner(FAST_SimulinkModel, opts);
    addBlock(ST0, {[ipc_mdl, 'PI - D/Kp'], [ipc_mdl, 'PI - D/Ki'], ...
        [ipc_mdl, 'PI - Q/Kp'], [ipc_mdl, 'PI - Q/Ki']});
    
    % Set tunable parameters
%     Kp_D = tunableGain('Kp_D', Parameters.cIPC.DQ_Kp_1P);
%     Ki_D = tunableGain('Ki_D', Parameters.cIPC.DQ_Ki_1P);
%     Kp_Q = tunableGain('Kp_Q', Parameters.cIPC.DQ_Kp_1P);
%     Ki_Q = tunableGain('Ki_Q', Parameters.cIPC.DQ_Ki_1P);
%     setBlockParam(ST0, 'PI - D/Kp', Kp_D, 'PI - D/Ki', Ki_D, 'PI - Q/Kp', Kp_Q, 'PI - Q/Ki', Ki_Q)
%     getBlockValue(ST0, [ipc_mdl, 'PI - D/Kp'])

    % add list of analysis points
    base_mdl = [FAST_SimulinkModel, '/', ipc_mdl];
    % inpt_ap = {[base_mdl, 'PI - D/b_d'], [base_mdl, 'PI - Q/b_q']};
    ref_ap = {[base_mdl, 'zero']};
%     op_ap = {[base_mdl, 'PI - D/Md'], [base_mdl, 'PI - Q/Mq']};
    op_ap = {'RootMyc1', 'RootMyc2', 'RootMyc3', [base_mdl, 'PI - D/Md'], [base_mdl, 'PI - Q/Mq']};
    addPoint(ST0, [ref_ap, op_ap]);
    getPoints(ST0)

     % Define Tuning Goals TODO how to define
    MarginReq = TuningGoal.Margins('AD_SOAR_c7_V2f_c73_Clean/Baseline Controller/Cyclic Pitch controller/1P Cyclic Pitch Controller1/1P IPCDQ Filtering/1[Md]', 6, 45); % 6dB of Gain Margin and 45deg of PM
    resp_time = 10; % in seconds, ST0.TimeUnit
    StepTrackingReq = TuningGoal.StepTracking(...
        {'zero', 'zero', 'zero'}, ...
        {'AD_SOAR_c7_V2f_c73_Clean/extract Myc1/1[RootMyc1]', 'AD_SOAR_c7_V2f_c73_Clean/extract Myc2/1[RootMyc2]', 'AD_SOAR_c7_V2f_c73_Clean/extract Myc3/1[RootMyc3]'}, 1);
        %{[ipc_mdl, 'PI - D/Md'], [ipc_mdl, 'PI - Q/Mq']}, resp_time);

    % Tune initial control system
    ST1 = systune(ST0, MarginReq);
    y = getBlockValue(ST1, 'PI - D');
    
    % get initial block value
    x = getBlockValue(ST0, 'PI - D/Kp');

    % compare step response of initial and tuned controllers
    T0 = getIOTransfer(ST0, {[base_mdl, 'PI - D/1[b_d]'], [base_mdl, 'PI - Q/1[b_q]']}, ...
        {[base_mdl, '1P IPCDQ Filtering/1[Md]'], [base_mdl, '1P IPCDQ Filtering/2[Mq]']});
    T1 = getIOTransfer(ST1, {[base_mdl, 'PI - D/1[b_d]'], [base_mdl, 'PI - Q/1[b_q]']}, ...
        {[base_mdl, '1P IPCDQ Filtering/1[Md]'], [base_mdl, '1P IPCDQ Filtering/2[Mq]']});
    step(T0, T1);
    legend('initial', 'tuned')
end


function [f, P1, Peaks] = getFFT(fft_vals, sim_indices, DT, Omega_1P) 
     Y = fft_vals(:, :, sim_indices);
     L = length(Y);
     P2 = abs(Y / L); % two-sided spectrum
     P1 = P2(1:L/2+1, :, :); % single-sided spectrum
     P1(2:end-1, :, :) = 2 * P1(2:end-1, :, :);
     P1 = P1(:, 1, :); % choose first blade only
     f = (1 / DT) * (0:(L/2))/L;

     Peaks = struct;
     Peaks.P1 = interp1(f, P1, Omega_1P);
     Peaks.P2 = interp1(f, P1, 2*Omega_1P);
     Peaks.P3 = interp1(f, P1, 3*Omega_1P);
end
