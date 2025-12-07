function results = ns_sweep_Vdc_Rload()
%NS_SWEEP_VDC_RLOAD  Sweep R_load and Vdc_nom and measure amplitude tracking.
%
% Uses inverter_base.slx with PID amplitude controller.
% Assumes ns_init_params.m sets defaults from workspace variables:
%   R_load, Vdc_nom, etc.

    % ---- sweep grids ----
    R_vals   = [5, 10, 20];      % ohm
    Vdc_vals = [400, 800];       % V

    % ---- simulation options ----
    Tstop       = 0.10;          % total sim time [s]
    Twindow_ss  = 0.02;          % use last 20 ms for steady state average

    % ---- preallocate results table ----
    results = table('Size',[0 4], ...
                    'VariableTypes', {'double','double','double','double'}, ...
                    'VariableNames', {'R_load','Vdc_nom','A_meas_mean','error_percent'});

    row = 1;

    for R = R_vals
        for Vdc = Vdc_vals

            % --- set parameters in workspace and re-init ---
            R_load  = R;          %#ok<NASGU>
            Vdc_nom = Vdc;        %#ok<NASGU>
            ns_init_params();     % uses R_load, Vdc_nom from base workspace

            fprintf('Running: R_load = %.1f ohm, Vdc_nom = %.0f V\n', R, Vdc);

            % --- run simulation ---
            simOut = sim('simulink_models/inverter_base', ...
                         'StopTime', num2str(Tstop), ...
                         'ReturnWorkspaceOutputs', 'on');

            % --- pull A_meas and A_ref ---
            A_meas_sig = simOut.A_meas;
            A_ref_sig  = simOut.A_ref;

            % Robust extraction (works for timeseries OR structure-with-time)
            [t, y_meas] = extract_sig(A_meas_sig);
            [~, y_ref]  = extract_sig(A_ref_sig);

            % --- steady-state window: last Twindow_ss seconds of A_meas ---
            t_end = t(end);
            t_ss_start = max(t(1), t_end - Twindow_ss);
            idx_ss = (t >= t_ss_start) & (t <= t_end);

            if ~any(idx_ss)
                warning('No samples in steady-state window! Check sample times.');
                A_meas_mean = NaN;
                A_ref_mean  = NaN;
                err_pct     = NaN;
            else
                A_meas_mean = mean(y_meas(idx_ss));

                % If A_ref has same length as t, average on same window,
                % otherwise just treat it as constant (common case).
                if numel(y_ref) == numel(t)
                    A_ref_mean = mean(y_ref(idx_ss));
                else
                    A_ref_mean = mean(y_ref);   % constant reference
                end

                err_pct = 100 * (A_meas_mean - A_ref_mean) / A_ref_mean;
            end

            % --- store into table ---
            results(row,:) = {R, Vdc, A_meas_mean, err_pct};
            row = row + 1;
        end
    end

    disp(results);
end


function [t, y] = extract_sig(sig)
%EXTRACT_SIG  Extract time and data from either timeseries or struct-with-time

    if isa(sig, 'timeseries')
        t = sig.Time;
        y = sig.Data;
    elseif isstruct(sig) && isfield(sig, 'time') && isfield(sig, 'signals')
        t = sig.time;
        y = sig.signals.values;
    else
        error('Unsupported signal type in extract_sig()');
    end
end
