function outPath = ns_run_rollout()
    projRoot = fileparts(mfilename('fullpath'));
    projRoot = fileparts(projRoot);              % go up from scripts/
    addpath(fullfile(projRoot,'scripts'));

    cd(projRoot);
    dataDir = fullfile(projRoot,'dataset','raw');
    if ~isfolder(dataDir), mkdir(dataDir); end

    simOut = sim('simulink_models/inverter_base', ...
        'StopTime','0.10', ...
        'ReturnWorkspaceOutputs','on');

    outPath = fullfile(dataDir, ...
        sprintf('rollout_%s.h5', datestr(now,'yyyymmdd_HHMMSS')));
    extract_and_save_rollout(simOut, outPath);
    fprintf('Saved rollout: %s\n', outPath);
end
