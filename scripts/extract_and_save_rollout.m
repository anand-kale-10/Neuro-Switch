function extract_and_save_rollout(simOut, out_h5_path)
% Robustly extract signals from a Simulink simOut and save to HDF5.
% Usage: extract_and_save_rollout(simOut, out_h5_path)
%
% Writes datasets:
%   /t, /Vout, /Iout, /Vdc, /m_oracle, /I_inv, /V_filter, /A_meas, /A_ref

if nargin < 2
    error('Usage: extract_and_save_rollout(simOut, out_h5_path)');
end

fprintf('Extractor: starting...\n');

% ----- helper: read one element safely -----
    function [t, y, ok] = read_element(name)
        t = []; y = []; ok = false;
        try
            val = simOut.get(name);
        catch
            try val = simOut.(name); catch val = []; end
        end
        if isempty(val), return; end

        % timeseries
        if isa(val,'timeseries')
            t = double(val.Time(:));
            y = double(val.Data(:));
            ok = true; return;
        end

        % StructureWithTime from To Workspace
        if isstruct(val) && isfield(val,'time')
            try
                t = double(val.time(:));
                if isfield(val,'signals') && isfield(val.signals,'values')
                    y = double(val.signals.values(:));
                elseif isfield(val,'values')
                    y = double(val.values(:));
                else
                    fn = fieldnames(val);
                    for ii = 1:numel(fn)
                        f = val.(fn{ii});
                        if isnumeric(f)
                            y = double(f(:)); break;
                        elseif isstruct(f) && isfield(f,'values') && isnumeric(f.values)
                            y = double(f.values(:)); break;
                        end
                    end
                end
                ok = ~isempty(y);
                return;
            catch
            end
        end

        % SimulationData.Signal with Time/Data
        try
            if isprop(val,'Time') && isprop(val,'Data')
                t = double(val.Time(:)); y = double(val.Data(:)); ok = true; return;
            end
        catch
        end

        % numeric array fallback
        if isnumeric(val)
            y = double(val(:));
            ok = true;
            return;
        end
    end

% ----- 1) master time vector (prefer tout) -----
master_t = [];
try
    tt = simOut.get('tout');
    if isa(tt,'timeseries'), master_t = double(tt.Time(:)); end
    if isempty(master_t) && isstruct(tt) && isfield(tt,'time')
        master_t = double(tt.time(:));
    end
catch
end

% signals we care about
names = {'Vout','Iout','Vdc','m_oracle','I_inv','V_filter','A_meas','A_ref'};
data = struct();

for k = 1:numel(names)
    nm = names{k};
    [t,y,ok] = read_element(nm);
    if ok
        data.(nm) = y;
        data.([nm '_t']) = t;
        if isempty(master_t)
            master_t = t;
        end
    else
        data.(nm) = [];
        data.([nm '_t']) = [];
    end
end

% fallback master_t if still empty
if isempty(master_t)
    fns = fieldnames(data);
    for i = 1:numel(fns)
        if endsWith(fns{i},'_t') && ~isempty(data.(fns{i}))
            master_t = data.(fns{i}); break;
        end
    end
end
if isempty(master_t)
    N = 0;
    for k = 1:numel(names)
        if ~isempty(data.(names{k})), N = numel(data.(names{k})); break; end
    end
    if N == 0, error('No numeric signals found in simOut to save.'); end
    master_t = (0:(N-1))';
end
master_t = master_t(:);

% ----- helper: align any signal to master_t -----
    function y_out = align_to_master(y_in, t_in)
        if isempty(y_in)
            y_out = nan(numel(master_t),1); return;
        end
        y_in = double(y_in(:));
        if isempty(t_in)
            if numel(y_in) == 1
                y_out = repmat(y_in, numel(master_t),1); return;
            elseif numel(y_in) == numel(master_t)
                y_out = y_in; return;
            else
                xq = linspace(1,numel(y_in), numel(master_t));
                y_out = interp1(1:numel(y_in), y_in', xq, 'linear', 'extrap')';
                return;
            end
        end
        t_in = double(t_in(:));
        if numel(t_in) == 1
            y_out = repmat(y_in(end), numel(master_t),1); return;
        end
        try
            y_out = interp1(t_in, y_in, master_t, 'linear', 'extrap'); y_out = y_out(:);
        catch
            y_out = interp1(t_in, y_in, master_t, 'nearest', 'extrap'); y_out = y_out(:);
        end
    end

% ----- 2) build arrays -----
T          = master_t;
Vout_arr   = align_to_master(data.Vout,      data.Vout_t);
Iout_arr   = align_to_master(data.Iout,      data.Iout_t);
Vdc_arr    = align_to_master(data.Vdc,       data.Vdc_t);
m_arr      = align_to_master(data.m_oracle,  data.m_oracle_t);
Iinv_arr   = align_to_master(data.I_inv,     data.I_inv_t);
Vfilter_arr= align_to_master(data.V_filter,  data.V_filter_t);
Ameas_arr  = align_to_master(data.A_meas,    data.A_meas_t);
Aref_arr   = align_to_master(data.A_ref,     data.A_ref_t);

% ----- 3) ensure output folder exists -----
[outDir,~,~] = fileparts(out_h5_path);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
    fprintf('Created output folder: %s\n', outDir);
end

% overwrite file if exists
if exist(out_h5_path,'file'), delete(out_h5_path); end

% helper: create & write 1D dataset
    function create_and_write(dset_path, arr)
        arr = double(arr(:));
        if isempty(arr), arr = zeros(numel(T),1); end
        h5create(out_h5_path, dset_path, size(arr), 'Datatype','double');
        h5write(out_h5_path, dset_path, arr);
    end

fprintf('Writing HDF5 to: %s\n', out_h5_path);
create_and_write('/t',        T);
create_and_write('/Vout',     Vout_arr);
create_and_write('/Iout',     Iout_arr);
create_and_write('/Vdc',      Vdc_arr);
create_and_write('/m_oracle', m_arr);
create_and_write('/I_inv',    Iinv_arr);
create_and_write('/V_filter', Vfilter_arr);
create_and_write('/A_meas',   Ameas_arr);
create_and_write('/A_ref',    Aref_arr);

fprintf('Done: wrote datasets: /t, /Vout, /Iout, /Vdc, /m_oracle, /I_inv, /V_filter, /A_meas, /A_ref\n');
end
