function run_openfast_simulink(params)
% Executes a single simulation run with input parameters

wind = params.wind; % Windspeed reference
seed = params.seed;
init_rot_speed = params.init_rot_speed;
init_blpitch = deg2rad(params.init_blpitch);
c_name = params.c_name;

%% Init.
Name_Model                  = params.model;
Name_Control                = params.c_file;
FAST_filename               = params.FAST_filename;
OutList                     = params.outlist_cache.OutList;

TMax                  = params.TMax;
Simulation.TMax       = TMax;
DT                    = 0.005 ;
Simulation.DT         = DT;
Control.DT            = DT;
Simulation.Floating   = 1; % ALWAYS 1 FOR USFLOWT
Simulation.AeroDynVer = '15';

Simulation.Name_Control = Name_Control;
Simulation.Name_Model   = Name_Model;

Parameters.Turbine.gen_eff = 96;


% Evaluate controller init function
eval(Simulation.Name_Control);

% Model initial conditions
IC_rotspd   = init_rot_speed; %in rpm
IC_genspd   = (R.GBRatio*(pi/30)*IC_rotspd); %in rad/s
IC.genTq    = R.VS_Rgn2K*IC_genspd^2 ;%in Nm
IC.pwr      = IC.genTq*IC_genspd; %in W
IC.bldPitch = init_blpitch; %rad
R.WE_v0     = wind;
R.WE_om0    = IC_rotspd;


% Create hidden worker directory
worker_dir = ['.worker_' num2str(params.thread_id)];
[~, ~, ~] = mkdir(worker_dir);
[~, ~, ~] = copyfile('Cp_Ct_Cq.DTU10MW_Nautilus.txt', worker_dir);
old_dir = cd(worker_dir);

% Run simulation
FAST_InputFileName = ['../' FAST_filename];

try
    sim([old_dir '/' Name_Model], [], simset('SrcWorkspace','current'));
    success = 1;
catch
    fprintf(['Simulation ' FAST_filename ' failed - continuing\n']);
    success = 0;
end

% Return to starting directory
cd(old_dir);

if success
    % Delete worker directory
    [~, ~, ~] = rmdir(worker_dir, 's');
    
    % Move output file
    [outdir, basename, ~] = fileparts(FAST_filename);
    movefile([outdir '/' basename '.SFunc.out'], [params.save_dir '/' basename '.out']);
end

fprintf('Thread %d finished\n', params.thread_id);

end
