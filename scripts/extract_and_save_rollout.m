function extract_and_save_rollout(simOut, out_h5_path)
% Robustly extract signals from a Simulink simOut and save to HDF5.
% Usage: extract_and_save_rollout(simOut, out_h5_path)
%
% Writes datasets: /t, /Vout, /Iout, /Vdc, /m_oracle, /I_inv, /V_filter

if nargin < 2
    error('Usage: extract_and_save_rollout(simOut, out_h5_path)');
end

fprintf('Extractor: starting...\n');

% Helper: robust reader for a named element from simOut
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
                % possible places for values
                if isfield(val,'signals') && isfield(val.signals,'values')
                    y = double(val.signals.values(:));
                elseif isfield(val,'values')
                    y = double(val.values(:));
                else
                    % attempt to find numeric field
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
                % fallthrough
            end
        end

        % SimulationData.Signal or other object with Time/Data or values
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

% 1) Build master time vector: prefer 'tout' if present
master_t = [];
try
    tt = simOut.get('tout');
    if isa(tt,'timeseries'), master_t = double(tt.Time(:)); end
    if isempty(master_t) && isstruct(tt) && isfield(tt,'time'), master_t = double(tt.time(:)); end
catch
    % nothing
end

% Names to extract (desired)
names = {'Vout','Iout','Vdc','m_oracle','I_inv','V_filter'};

data = struct();

% read each element; collect times if different
for k=1:numel(names)
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

% final fallback: if master_t still empty, try any *_t field
if isempty(master_t)
    fns = fieldnames(data);
    for i=1:numel(fns)
        if endsWith(fns{i},'_t') && ~isempty(data.(fns{i}))
            master_t = data.(fns{i}); break;
        end
    end
end

% if still empty, create generic index time
if isempty(master_t)
    warning('No time vector found; creating default 0..N-1 sample index as time.');
    % find first non-empty data vector length
    N = 0;
    for k=1:numel(names)
        if ~isempty(data.(names{k})), N = numel(data.(names{k})); break; end
    end
    if N == 0, error('No numeric signals found in simOut to save.'); end
    master_t = (0:(N-1))';
end

% Ensure master_t is column
master_t = master_t(:);

% Helper to expand/align y to master_t length
    function y_out = align_to_master(y_in, t_in)
        % returns column vector same length as master_t
        if isempty(y_in)
            y_out = nan(numel(master_t),1);
            return;
        end
        y_in = double(y_in(:));

        % If t_in empty or single sample, handle scalar/single-sample cases
        if isempty(t_in)
            if isscalar(y_in)
                y_out = repmat(y_in, numel(master_t),1);
                return;
            else
                % no time base but vector: try to align by length
                if numel(y_in) == numel(master_t)
                    y_out = y_in;
                    return;
                else
                    % resample by linear interpolation over index
                    xq = linspace(1,numel(y_in), numel(master_t));
                    y_out = interp1(1:numel(y_in), y_in', xq, 'linear', 'extrap')';
                    return;
                end
            end
        end

        t_in = double(t_in(:));
        % If only one time sample, repeat value
        if numel(t_in) == 1
            y_out = repmat(y_in(end), numel(master_t),1);
            return;
        end

        % now both t_in and master_t have >= 2 samples => safe to interp
        try
            y_out = interp1(t_in, y_in, master_t, 'linear', 'extrap');
            y_out = double(y_out(:));
        catch
            % fallback to nearest
            y_out = interp1(t_in, y_in, master_t, 'nearest', 'extrap');
            y_out = double(y_out(:));
        end
    end

% Build final arrays to write
T = master_t;
Vout_arr = align_to_master(data.Vout, data.Vout_t);
Iout_arr = align_to_master(data.Iout, data.Iout_t);
Vdc_arr  = align_to_master(data.Vdc,  data.Vdc_t);
m_arr    = align_to_master(data.m_oracle, data.m_oracle_t);
Iinv_arr = align_to_master(data.I_inv, data.I_inv_t);
Vfilter_arr = align_to_master(data.V_filter, data.V_filter_t);

% Create / overwrite HDF5 file
try
    if exist(out_h5_path,'file')
        delete(out_h5_path);
    end
catch
end

% Helper to create & write datasets (handles 1D arrays)
    function create_and_write(dset_path, arr)
        arr = double(arr(:));
        if isempty(arr)
            arr = zeros(numel(T),1);
        end
        % Ensure dataset dims
        if exist(out_h5_path,'file')==2
            % ok
        end
        h5create(out_h5_path, dset_path, size(arr), 'Datatype', 'double');
        h5write(out_h5_path, dset_path, arr);
    end

fprintf('Writing HDF5 to: %s\n', out_h5_path);
create_and_write('/t', T);
create_and_write('/Vout', Vout_arr);
create_and_write('/Iout', Iout_arr);
create_and_write('/Vdc', Vdc_arr);
create_and_write('/m_oracle', m_arr);
create_and_write('/I_inv', Iinv_arr);
create_and_write('/V_filter', Vfilter_arr);

fprintf('Done: wrote datasets: /t, /Vout, /Iout, /Vdc, /m_oracle, /I_inv, /V_filter\n');

end
