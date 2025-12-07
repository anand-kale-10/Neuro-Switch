function ns_init_params()
% NS_INIT_PARAMS  Central parameter init for NeuroSwitch inverter model.
% This is called from inverter_base InitFcn.

%% ==== Plant rating (15 kW grid-tied) ====

P_rated      = 15e3;      % 15 kW
Vgrid_rms    = 230;       % per-phase rms (assuming single-phase model)
Vgrid_peak   = Vgrid_rms * sqrt(2);  % ~325 V
A_ref       = Vgrid_peak;

%% ==== LCL FILTER PARAMETERS ====
% Inductors
L1 = 5.066e-3;    % H
L2 = 5.066e-3;    % H

% Capacitor
C_filter = 10e-6; % F

% Damping / series resistances
R_damp    = 5.0;   % ohm  (in series with C or damping branch)
R_series1 = 0.05;  % ohm  (L1 series)
R_series2 = 0.05;  % ohm  (L2 series)

%% ==== LOAD & DC BUS ====
R_load   = 10.0;          % ohm (baseline)
Vdc_nom      = 400;       % example nominal DC bus (adjust if needed)
I_nom_rms    = P_rated / Vgrid_rms;  % ~65 A
I_nom        = I_nom_rms;           % use rms as nominal scale for now

%% ==== AC FUNDAMENTAL & ENVELOPE DETECTOR ====
f_ac   = 50;                % Hz (output fundamental)
w_ac   = 2*pi*f_ac;

% Envelope lowpass (A_meas)
f_env  = 16.0;                % Hz cutoff (envelope filter)
tau_env = 1/(2*pi*f_env);   % s
% In Simulink Transfer Fcn: num = 1, den = [tau_env 1]

%% ==== VOLTAGE CONTROLLER (PI ON AMPLITUDE) ====
% We designed PID assuming approx gain K ~ Vdc_nom and target tau_cl ~ 0.02 s
tau_cl = 0.02;              % desired closed-loop time constant (s)
K_plant = Vdc_nom;          % rough small-signal m -> Vout amplitude gain

Ts_ctrl = 1e-4;

Kp_V = 1 / (K_plant * tau_cl);  % proportional gain
Ki_V = Kp_V / tau_cl;           % integral gain (continuous-time)

m_max = 1;

%% ==== PUSH VARIABLES TO BASE OR MODEL WORKSPACE ====
% This function is called from the model InitFcn, so using assignin('base',...)
% makes them visible to blocks that use workspace evaluation.

vars = who;  % all local variables
for k = 1:numel(vars)
    if strcmp(vars{k}, 'ans'), continue; end
    assignin('base', vars{k}, eval(vars{k}));
end

fprintf('[ns_init_params] Parameters initialized: LCL, load, Vdc, PI gains, Vgrid_peak=%.1f V\n', Vgrid_peak);
end