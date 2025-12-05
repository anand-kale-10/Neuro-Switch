%% inspect_rollouts.m
% Inspect HDF5 rollouts in dataset/raw and print a concise summary.

clear; clc;

%% === Locate dataset/raw folder ===
scriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptPath);            % one level up
dataDir = fullfile(projectRoot, 'dataset', 'raw');

fprintf("Project root detected: %s\n", projectRoot);

if ~isfolder(dataDir)
    fprintf("ERROR: dataset/raw folder not found at: %s\n", dataDir);
    return;
end

files = dir(fullfile(dataDir, '*.h5'));
if isempty(files)
    fprintf("No .h5 files found in dataset/raw.\n");
    return;
end

fprintf("\n=== Rollout Inspection (%d files) ===\n", numel(files));

%% === Loop over all HDF5 files ===
for i = 1:numel(files)

    fname = fullfile(files(i).folder, files(i).name);
    fprintf("\n--- %s ---\n", files(i).name);

    % Get dataset names
    try
        info = h5info(fname);
    catch ME
        fprintf("  ERROR reading file: %s\n", ME.message);
        continue;
    end

    dnames = {info.Datasets.Name};

    fprintf("  Datasets found: %s\n", strjoin(dnames, ", "));

    % Summaries
    dataset_list = {"t", "Vout", "Iout", "Vdc", "m_oracle"};

    for k = 1:numel(dataset_list)
        dsname = dataset_list{k};

        if ismember(dsname, dnames)
            try
                data = h5read(fname, "/" + dsname);
                data = double(data(:));
                fprintf("  %s: len=%d, min=%.4g, max=%.4g, mean=%.4g\n", ...
                    dsname, numel(data), min(data), max(data), mean(data));
            catch
                fprintf("  %s: ERROR reading dataset\n", dsname);
            end
        else
            fprintf("  %s: [MISSING]\n", dsname);
        end
    end

    % File size
    f = dir(fname);
    fprintf("  File size: %.1f KB\n", f.bytes/1024);

end

fprintf("\n=== END OF INSPECTION ===\n");
