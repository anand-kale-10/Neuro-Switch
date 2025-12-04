function extract_and_save_rollout()
% NeuroSwitch - Extract and Save Rollout to HDF5
% Works with any simOut structure produced by inverter_base model
% No dependencies. Clean, robust, and reliable.

fprintf('\n=== NeuroSwitch Rollout Extractor ===\n');

%% Configure paths
baseRoot = fullfile('C:','Users','user','OneDrive','Documents','MATLAB','NeuroSwitch_Lab');
outFolder = fullfile(baseRoot,'dataset','raw');
if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% Require simOut
if ~evalin('base','exist(''simOut'',''var'')')
    error('simOut does not exist in base workspace. Run simulation with ReturnWorkspaceOutputs on.');
end

simOut = evalin('base','simOut');

%% Get list of element names
try 
    names = simOut.getElementNames();
catch
    names = simOut.who();
end
disp('simOut elements:'); disp(names);

%% ------------------ FIND MASTER TIME VECTOR ------------------
t_num = [];

% Prefer element 'tout'
if any(strcmp(names,'tout'))
    tmp = simOut.get('tout');
    if isnumeric(tmp)
        t_num = double(tmp(:));
    end
end

% If not found, search for fields Time/time in elements
if isempty(t_num)
    for ii = 1:numel(names)
        E = simOut.get(names{ii});
        try
            if isstruct(E) && isfield(E,'Time')
                t_num = double(E.Time(:));
                break;
            end
            if isstruct(E) && isfield(E,'time')
                t_num = double(E.time(:));
                break;
            end
            if isa(E,'timeseries')
                t_num = double(E.Time(:));
                break;
            end
        catch
        end
    end
end

if isempty(t_num)
    error('No time vector found in simOut.');
end

Tlen = numel(t_num);
fprintf('Master time length = %d\n', Tlen);

%% ------------------ EXTRACT NUMERIC ARRAYS ------------------
Vnum = []; Inum = []; Dnum = []; Mnum = [];

for ii = 1:numel(names)
    nm = names{ii};
    E = simOut.get(nm);

    extracted = [];

    % Case 1: timeseries
    try
        if isa(E,'timeseries')
            extracted = double(E.Data(:));
        end
    catch; end

    % Case 2: Struct with Data/Time
    if isempty(extracted)
        try
            if isstruct(E) && isfield(E,'Data') && isfield(E,'Time')
                extracted = double(E.Data(:));
            end
        catch; end
    end

    % Case 3: Legacy struct: signals.values or signals(1).values
    if isempty(extracted)
        try
            if isstruct(E) && isfield(E,'signals')
                S = E.signals;

                if isstruct(S) && isfield(S,'values')
                    extracted = double(S.values(:));
                elseif iscell(S) && ~isempty(S)
                    c = S{1};
                    if isstruct(c) && isfield(c,'values')
                        extracted = double(c.values(:));
                    end
                elseif isstruct(S) && numel(S)>=1 && isfield(S(1),'values')
                    extracted = double(S(1).values(:));
                end
            end
        catch; end
    end

    % Case 4: first numeric field
    if isempty(extracted)
        try
            if isstruct(E)
                fns = fieldnames(E);
                for jj = 1:numel(fns)
                    f = fns{jj};
                    if any(strcmpi(f,{'time','Time','blockName'})), continue; end
                    val = E.(f);
                    if isnumeric(val) && ~isempty(val)
                        extracted = double(val(:));
                        break;
                    end
                end
            end
        catch; end
    end

    % Case 5: raw numeric
    if isempty(extracted)
        try
            if isnumeric(E)
                extracted = double(E(:));
            end
        catch; end
    end

    % Assign to correct output variable
    if strcmp(nm,'Vout'),     Vnum = extracted; end
    if strcmp(nm,'Iout'),     Inum = extracted; end
    if strcmp(nm,'Vdc'),      Dnum = extracted; end
    if strcmp(nm,'m_oracle'), Mnum = extracted; end
end

%% ------------------ ALIGN TO TIME LENGTH ------------------
alignOne = @(arr) local_align(arr, Tlen);

V_al = alignOne(Vnum);
I_al = alignOne(Inum);
D_al = alignOne(Dnum);
M_al = alignOne(Mnum);

fprintf('\nAligned sizes:\n');
fprintf(' t       : %s\n', mat2str(size(t_num)));
fprintf(' Vout    : %s\n', mat2str(size(V_al)));
fprintf(' Iout    : %s\n', mat2str(size(I_al)));
fprintf(' Vdc     : %s\n', mat2str(size(D_al)));
fprintf(' m_oracle: %s\n', mat2str(size(M_al)));

%% ------------------ WRITE HDF5 ------------------
outFile = fullfile(outFolder, 'rollout_manual_fixed.h5');
if exist(outFile,'file'), delete(outFile); end

h5create(outFile,'/t',size(t_num)); h5write(outFile,'/t',t_num);

if ~isempty(V_al), h5create(outFile,'/Vout',size(V_al)); h5write(outFile,'/Vout',V_al); end
if ~isempty(I_al), h5create(outFile,'/Iout',size(I_al)); h5write(outFile,'/Iout',I_al); end
if ~isempty(D_al), h5create(outFile,'/Vdc',size(D_al)); h5write(outFile,'/Vdc',D_al); end
if ~isempty(M_al), h5create(outFile,'/m_oracle',size(M_al)); h5write(outFile,'/m_oracle',M_al); end

fprintf('\nSaved HDF5 → %s\n', outFile);
info = h5info(outFile);
for k = 1:numel(info.Datasets)
    fprintf(' - /%s  size=%s\n', info.Datasets(k).Name, mat2str(info.Datasets(k).Dataspace.Size));
end

fprintf('\n=== Rollout extraction complete ===\n');

end

%% --------- Local helper: alignment function ---------
function out = local_align(arr, Tlen)
    out = [];
    if isempty(arr), return; end
    if isscalar(arr)
        out = repmat(double(arr), Tlen, 1); return
    end
    if size(arr,1) == Tlen
        out = double(arr); return
    end
    if size(arr,2) == Tlen
        out = double(arr)'; return
    end
    if numel(arr) == Tlen
        out = reshape(double(arr),[],1); return
    end
    % otherwise can't align: leave empty
end
