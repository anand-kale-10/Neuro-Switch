function results = ns_sweep_Rload_400V()
%NS_SWEEP_RLOAD_400V  Sweep R_load at Vdc_nom = 400 V and measure amplitude tracking.
%
% Uses inverter_base.slx with PID amplitude controller.
% Assumes ns_init_params.m sets:
%   R_load, Vdc_nom, Vgrid_peak (A_ref) in base workspace.
%
% A_meas is expected as "Structure with time" from To Workspace:
%   .time, .signals.values

    % ---- sweep grid ----
    R_vals = [5, 10, 20];   % ohm, or whatever you want

    % ---- simulation options ----
    Tstop      = 0.10;      % total sim time [s]
    Twindow_ss = 0.02;      % use last 20 ms for steady-state average

    % ---- preallocate results table ----
    results = table('Size',[0 4], ...
        'VariableTypes', {'double','double','double','double'}, ...
        'VariableNames', {'R_load','Vdc_nom','A_meas_mean','error_percent'});

    row = 1;

    for R = R_vals

        % --- set parameters in base workspace and re-init ---
        R_load  = R;        %#ok<NASGU>
        Vdc_nom = 400;      %#ok<NASGU>

        ns_init_params();   % uses R_load, Vdc_nom from base workspace

        % Get reference amplitude (constant) from workspace
        A_ref_nom = evalin('base', 'Vgrid_peak');  % 325.3 V in your log

        fprintf('Running: R_load = %.1f ohm, Vdc_nom = %.0f V\n', R, Vdc_nom);

        % --- run simulation ---
        simOut = sim('simulink_models/inverter_base', ...
                     'StopTime', num2str(Tstop), ...
                     'ReturnWorkspaceOutputs', 'on');

        % --- pull A_meas as Structure-with-time ---
        A_meas_ts = simOut.A_meas;

        % Structure with time: fields .time and .signals.values
        t      = A_meas_ts.time;
        y_meas = A_meas_ts.signals.values;

        % --- steady-state window: last Twindow_ss seconds of A_meas ---
        t_end      = t(end);
        t_ss_start = t_end - Twindow_ss;
        idx_ss     = (t >= t_ss_start) & (t <= t_end);

        if ~any(idx_ss)
            warning('No samples in steady-state window! Check sample times.');
            A_meas_mean = NaN;
            err_pct     = NaN;
        else
            A_meas_mean = mean(y_meas(idx_ss));

            % A_ref is constant = A_ref_nom
            A_ref_mean  = A_ref_nom;

            err_pct = 100 * (A_meas_mean - A_ref_mean) / A_ref_mean;
        end

        % --- store into table ---
        results(row,:) = {R, Vdc_nom, A_meas_mean, err_pct};
        row = row + 1;
    end

    disp(results);
end
