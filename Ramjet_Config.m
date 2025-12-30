% RAMJET_CONFIG
% global constants, flight conditions,  derived variables.
% call it before running generation or solvers.

clear; clc;

%% Inputs
% -------------------------------------------------------------------------
inputs.Mach      = input('Enter Freestream Mach Number (Default 4.0): ');
if isempty(inputs.Mach), inputs.Mach = 4.0; end

inputs.Theta_deg = input('Enter Cone Half-Angle [deg] (Default 15): ');
if isempty(inputs.Theta_deg), inputs.Theta_deg = 15; end

inputs.Altitude  = 20000; %m 

% grid
inputs.IMax = 200; % axial
inputs.JMax = 100;  % Radial

%% constants of air
% -------------------------------------------------------------------------
phys.gamma = 1.4;
phys.R     = 287.05;
phys.Cv    = phys.R / (phys.gamma - 1);
phys.Cp    = phys.gamma * phys.R / (phys.gamma - 1);

%% 3. FS conditions, standard atmosphere @ 20km
% -------------------------------------------------------------------------
T_inf = 216.65;     % [K] Temperature
P_inf = 5529.0;     % [Pa] Pressure
rho_inf = P_inf / (phys.R * T_inf);

% speed and velocity of sound
a_inf = sqrt(phys.gamma * phys.R * T_inf);
U_inf = inputs.Mach * a_inf;
prim.rho_inf = rho_inf;
prim.u_inf   = U_inf;
prim.v_inf   = 0.0;
prim.p_inf   = P_inf;

%% conservative variables
% -------------------------------------------------------------------------
E_inf = rho_inf * (phys.Cv * T_inf + 0.5 * U_inf^2);
IC.U1 = rho_inf;
IC.U2 = rho_inf * U_inf;
IC.U3 = 0.0;
IC.U4 = E_inf;

fprintf('Configuration Loaded.\n');
fprintf('  > Mach: %.2f\n', inputs.Mach);
fprintf('  > Velocity: %.2f m/s\n', U_inf);
fprintf('  > Initial Energy Density: %.2f J/m^3\n', E_inf);
% --- COMBUSTION PARAMS ---
% modelled heat addition as a volumetric source, Delta_T: approx 1000K-1500K.
% Q_vol approx rho * Cp * Delta_T * (Velocity / Length)
% Estimate: 1.0 kg/m3 * 1000 J/kgK * 1200 K * (1000 m/s / 1.0 m) = 1.2e9 W/m^3
% USER NOTE:REDUCE inputs.Heat_Release_Max WHEN NEEDED TO PRVENENT CRASHING
% --- HEAT  ---
inputs.Heat_Release_Max = 4.0e7; 
inputs.Burn_Start_Iter  = 1000;
inputs.Burn_Ramp_Length = 1000;