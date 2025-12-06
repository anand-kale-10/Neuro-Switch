function ns_generate_pid_rollouts()
%NS_GENERATE_PID_ROLLOUTS
% Generate a small grid of PID-controlled rollouts for NeuroSwitch.
% Run from anywhere; this function assumes it lives in NeuroSwitch_Lab/scripts.

clc;
fprintf('=== NeuroSwitch PID rollout generator ===\n');

%% 0) Find project root from this file location
thisFile   = mfilename('fullpath');
[scriptsDir,~] = fileparts(thisFile);           % .../NeuroSwitch_Lab/scripts
projRoot     = fileparts(scriptsDir);           % .../NeuroSwitch_Lab

fprintf('Project root detected: %s\n', projRoot);

cd(projRoot);                                   % work from project root

%% 1) Paths
dataDir = fullfile(projRoot,'dataset','raw');
if ~isfolder(dataDir)
    mkdir(dataDir);
    fprintf('Created dataDir: %s\n', dataDir);
end

addpath(fullfile(projRoot,'scripts'));          % for extractor + init script

%% 2) Load Simulink model
modelName = 'inverter_base';
modelPath = fullfile(projRoot,'simulink_models',modelName);

if ~exist([modelPath '.slx'],'file')
    error('Cannot find model file: %s.slx', modelPath);
end

load_system(modelPath);                         % loads under name 'inverter_base'
fprintf('Loaded model: %s\n', modelName);

%% 3) Rollout grid (same as your A_meas table)
R_values   = [5 10 20];                         % ohm
Vdc_values = [350 400 450];                     % V

totalRuns = numel(R_values)*numel(Vdc_values);
runIdx = 0;

%% 4) Main loop
for i = 1:numel(R_values)
    for j = 1:numel(Vdc_values)
        R  = R_values(i);
        Vd = Vdc_values(j);
        runIdx = runIdx + 1;

        fprintf('\n=== Run %d/%d: R=%.1f, Vdc=%.0f ===\n', ...
                runIdx, totalRuns, R, Vd);

        % 4a) Initialise model parameters
        run(fullfile(projRoot,'scripts','ns_init_params.m'));  % sets defaults

        % override R_load & Vdc_nom for this run
        ws = get_param(modelName,'ModelWorkspace');
        assignin(ws,'R_load',  R);
        assignin(ws,'Vdc_nom', Vd);

        % 4b) Simulate
        simOut = sim(modelName, ...
            'StopTime','0.10', ...
            'ReturnWorkspaceOutputs','on');

        % 4c) Save rollout
        fname  = sprintf('rollout_pid_R%.1f_Vdc%.0f_%s.h5', ...
                         R, Vd, datestr(now,'HHMMSSFFF'));
        outH5  = fullfile(dataDir, fname);

        extract_and_save_rollout(simOut, outH5);

        fprintf('Saved rollout: %s\n', outH5);
    end
end

fprintf('\n=== All %d PID rollouts generated. ===\n', totalRuns);
end
