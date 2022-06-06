%% Run FAST from Simulink in parallel

%% TODO
% use generateCases function
% adapt for mac and unix
% adapt for most recent fast version

close all;
clc;
clearvars;

if ismac
    addpath('/Users/aoifework/Documents/dev/WEIS/OpenFAST/install/bin');
    addpath('/Users/aoifework/Documents/toolboxes/matlab-toolbox/A_Functions');
elseif isunix
    addpath('/projects/aohe7145/toolboxes/dev/WEIS/OpenFAST/install/bin');
    addpath('/projects/aohe7145/toolboxes/matlab-toolbox/A_Functions');
end

save_dir = '../controller_verf/model3.0_Dist_v1_DLC11_SL';
[~, ~, ~] = mkdir(save_dir);

fast_IF_dir = 'FAST8_IF';
Name_Model = 'USFL_SL_MDL_c6_r2018a';
% Name_Model = 'USFL_SL_MDL_c5_elenya_fullBaseline_attempt06_newFF';
c4_dir = '../controller_dev/';
addpath(c4_dir);

% Metocean environment conditions
load_metocean;
outlist_cache = load('outlist_cache');


main_wd = cd(fast_IF_dir);

%% Start embedded python environment
pe = pyenv;

%% Set up parameters for simulations to run

% Dynamic case
enable_inflow = 1;
enable_aero = 2;
enable_servo = 1;

wind_type = 3;
wave_mod = 2;
rot_speed_init = 9.6;

% Wind cases
% urefs = 4:24;
urefs = 4:2:24;
% urefs = [10 11:2:23];
% urefs = 12:24;
% urefs = 24;
% urefs = 14;
% urefs = [4 10 12 14 16 24];
urefs = fliplr(urefs);

% Add controllers to evaluate
% controller_cb    = true;
controller_cb    = false;
% cb_forward_tilt  = 5; % deg
cb_forward_tilt  = 0; % deg
controller_names = {};
controller_files = {};

controller_names = cat(2, controller_names, 'c6');
controller_files = cat(2, controller_files, 'simulink_rosco_c4_redone');

% controller_names = cat(2, controller_names, 'c6-fl-dual-comp');
% controller_files = cat(2, controller_files, 'simulink_rosco_c4_fl_gentq_shared_sched_sat');

controllers = 1:length(controller_names);
if controller_cb
    for ii = 1:length(controller_names)
        controller_names{ii} = [controller_names{ii} '-CB1'];
%         controller_names{ii} = [controller_names{ii} '-CB2'];
    end
end

% Random seeds
% num_seeds = 6;
% num_seeds = 2;
% seeds = 1:num_seeds;
% seeds = 1:6;
seeds = 1;

TMax = 800;

% Some empirical data
model_data = load('../v5_mean_data');


% Parameter file names
FAST_template_filename = 'Model.fst.t';
FAST_filename_format = 'U_%d_C_%s_S_%d_T_%d';

inflow_template_filename = 'InflowWind.dat.t';
inflow_filename_format = 'InflowWind_U_%d_S_%d.dat';

hydro_template_filename = 'HydroDyn.dat.t';
hydro_filename_format = 'HydroDyn_U_%d.dat';

servo_template_filename = 'ServoDyn.dat.t';
servo_filename_format = 'ServoDyn_C_%s.dat';

elasto_template_filename = 'ElastoDyn.dat.t';
elasto_filename_format = 'ElastoDyn_U_%d.dat';

moor_template_filename = 'MoorDyn.dat.t';
moor_filename_format = 'MoorDyn_U_%d.dat';

% Load template files
% FAST_template_data = fileread(['templates/' FAST_template_filename]);
FAST_template_file = py.open(['templates/' FAST_template_filename]);
FAST_template_data = FAST_template_file.read();
FAST_template_file.close();
FAST_template_object = py.string.Template(FAST_template_data);

inflow_template_data = fileread(['templates/' inflow_template_filename]);
inflow_template_object = py.string.Template(inflow_template_data);

hydro_template_data = fileread(['templates/' hydro_template_filename]);
hydro_template_object = py.string.Template(hydro_template_data);

% servo_template_data = fileread(['templates/' servo_template_filename]);
% servo_template_object = py.string.Template(servo_template_data);

elasto_template_data = fileread(['templates/' elasto_template_filename]);
elasto_template_object = py.string.Template(elasto_template_data);

moor_template_data = fileread(['templates/' moor_template_filename]);
moor_template_object = py.string.Template(moor_template_data);


% Load parameters for each running instance (slow because structure is not
% preallocated, but only needs to be done once at init)
num_runs = length(seeds) * length(controllers) * length(urefs);
run_index = 1;
for u_index = 1:length(urefs)
    wind_select = urefs(u_index);

    % Populate case-specific template files

    % ElastoDyn (rotor speed)
%     rot_speed_select = rot_speed_init; % Could replace with lookup table for region 2
    rot_speed_select = interp1(model_data.wind, model_data.rotspd_data, wind_select);
    blpitch_select = interp1(model_data.wind, model_data.blpitch_data, wind_select);
    elasto_filename = sprintf(['ED/' elasto_filename_format], urefs(u_index));
    elasto_output_file = fopen(elasto_filename, 'w');
    fwrite(elasto_output_file, elasto_template_object.substitute(py.dict(pyargs(...
        'INITROTSPD', rot_speed_select, ...
        'INITBLPITCH', blpitch_select))).char());
    fclose(elasto_output_file);
    
    % MoorDyn (can ballast)
    if controller_cb
        fine_ballast = 1e-5;
        base_ballast = 2;
        cb_level = 0.63186 * (interp1(model_data.wind, model_data.pitch_data, wind_select) + cb_forward_tilt);
        moor_filename = sprintf(['MD/' moor_filename_format], urefs(u_index));
        moor_output_file = fopen(moor_filename, 'w');
        fwrite(moor_output_file, moor_template_object.substitute(py.dict(pyargs(...
            'BAL1LEN', base_ballast + cb_level, 'BAL2LEN', base_ballast+fine_ballast, 'BAL3LEN', base_ballast+fine_ballast))).char());
        fclose(moor_output_file);
    else
        moor_filename = 'MD/MoorDyn.dat';
    end

    for s_index = 1:length(seeds)
        seed_select = seeds(s_index);

        % HydroDyn (metocean conditions)
        wavehs_select = interp1(met_wind_table, met_wavehs_table, wind_select);
        wavetp_select = interp1(met_wind_table, met_wavetp_table, wind_select);
        wavepkshp_select = interp1(met_wind_table, met_wavepkshp_table, wind_select);
        wave_s1 = round(wind_select * seed_select * wavehs_select * wavetp_select * wavepkshp_select) + 123;
        wave_s2 = wave_s1 * wave_s1 + 456;
        hydro_filename = sprintf(['HD/' hydro_filename_format], wind_select);
        hydro_output_file = fopen(hydro_filename, 'w');
        fwrite(hydro_output_file, hydro_template_object.substitute(py.dict(pyargs(...
            'WAVEMOD', num2str(wave_mod), 'WAVETMAX', num2str(TMax + 100),...
            'WAVEHS', wavehs_select, 'WAVETP', wavetp_select, 'WAVEPKSHP', wavepkshp_select,...
            'WAVES1', num2str(wave_s1), 'WAVES2', num2str(wave_s2)))).char());
        fclose(hydro_output_file);

        % InflowWind (wind speed, random seed)
        inflow_filename = sprintf(['IfW/' inflow_filename_format], wind_select, seed_select);
        inflow_output_file = fopen(inflow_filename, 'w');
        fwrite(inflow_output_file, inflow_template_object.substitute(py.dict(pyargs(...
                'WINDTYPE', num2str(wind_type), 'UREF', num2str(wind_select), 'SEED', num2str(seed_select), 'ETM', ''))).char());
        fclose(inflow_output_file);

        for c_index = 1:length(controllers)
            c_name = controller_names{c_index};
            c_file = controller_files{c_index};

            % FAST input file (need separate model file for each controller so output files don't collide)
            FAST_filename_base = sprintf(FAST_filename_format, ...
                    wind_select, c_name, seed_select, TMax);
            FAST_filename = ['FASTInput/' FAST_filename_base '.fst'];
            FAST_output_file = fopen(FAST_filename, 'w');
            fwrite(FAST_output_file, FAST_template_object.substitute(py.dict(pyargs('TMAX', num2str(TMax), ...
                    'COMP_INFLOW', num2str(enable_inflow), 'COMP_AERO', num2str(enable_aero), 'COMP_SERVO', num2str(enable_servo), ...
                    'ELASTODYN_DAT', ['../' elasto_filename], ...
                    'INFLOWWIND_DAT', ['../' inflow_filename], ...
                    'HYDRODYN_DAT', ['../' hydro_filename], ...
                    'SERVODYN_DAT', '../SrvD/ServoDyn_SL.dat', ...
                    'SUBDYN_DAT', '../SD/SubDyn.dat', ...
                    'MOORDYN_DAT', ['../' moor_filename]...
                    ))).char());
            fclose(FAST_output_file);


            % Populate thread parameters
            sim_params(run_index).thread_id = run_index;
            sim_params(run_index).TMax = TMax;
            sim_params(run_index).model = Name_Model;
            sim_params(run_index).outlist_cache = outlist_cache;
            
            sim_params(run_index).wind = wind_select;
            sim_params(run_index).c_name = c_name;
            sim_params(run_index).c_file = c_file;
            sim_params(run_index).seed = seed_select;
            sim_params(run_index).init_rot_speed = rot_speed_select;
            sim_params(run_index).init_blpitch = blpitch_select;

            % Input filenames
            sim_params(run_index).save_dir = save_dir;
            sim_params(run_index).FAST_filename = [fast_IF_dir '/' FAST_filename];
            sim_params(run_index).inflow_filename = inflow_filename;
            sim_params(run_index).hydro_filename = hydro_filename;

            run_index = run_index + 1;
        end
    end
end

% Return to starting directory
cd(main_wd);


%% Run simulations in multiple parallel threads

num_threads = 10; % Just shy of the 12 present in machine

% Automatically identify runs that haven't finished (e.g if the sim was interrupted)
run_select = 1%:num_runs;
% find_done;
num_thread_runs = length(run_select);
num_batches = ceil(num_thread_runs / num_threads);
fprintf('Executing %d runs of %d (expect approx %.1f hours to completion)\n', ...
    num_thread_runs, num_runs, 1.5*num_batches);
%parfor (run_index = 1:num_thread_runs, num_threads)
for run_index = 1:num_thread_runs
    thread_index = run_select(run_index);
    run_openfast_simulink(sim_params(thread_index));
end

fprintf('\nOpenFAST Simulations Completed\n\n\n');

