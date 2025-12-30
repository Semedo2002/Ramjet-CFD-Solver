%PERFORMANCE_ANALYSIS
%mass balance and  aligning cell centers with faces.

fprintf('\n=== RAMJET PERFORMANCE REPORT ===\n');

% geometry
% -------------------------------------------------------------------------
% --- i=1---
dy = GRID.Y(2:end, 1) - GRID.Y(1:end-1, 1);
Y_face = 0.5 * (GRID.Y(2:end, 1) + GRID.Y(1:end-1, 1));
A_inlet_vec = sqrt(dy.^2) .* (2 * pi * Y_face);

% ---i=end---
dy = GRID.Y(2:end, end) - GRID.Y(1:end-1, end);
Y_face = 0.5 * (GRID.Y(2:end, end) + GRID.Y(1:end-1, end));
A_exit_vec = sqrt(dy.^2) .* (2 * pi * Y_face);

% ---j=end ---
dx_top = GRID.X(end, 2:end) - GRID.X(end, 1:end-1);
dy_top = GRID.Y(end, 2:end) - GRID.Y(end, 1:end-1);
S_top  = sqrt(dx_top.^2 + dy_top.^2);
Y_top  = 0.5 * (GRID.Y(end, 2:end) + GRID.Y(end, 1:end-1));

A_top_vec = S_top .* (2 * pi * Y_top);

% n.v
nx_top = -dy_top ./ S_top;
ny_top =  dx_top ./ S_top;

% mask before cowl 
L_inlet = 1.2; 
x_centers = 0.5 * (GRID.X(end, 1:end-1) + GRID.X(end, 2:end));
mask_open = x_centers < L_inlet;

%flow integration
% -------------------------------------------------------------------------

% ---inlet fluxes ---
rho_in = U(:, 1, 1);
u_in   = U(:, 1, 2) ./ (rho_in+1e-16);
mdot_inlet = sum( rho_in .* u_in .* A_inlet_vec );

% ---entrainment---
rho_t = U(end, :, 1);
u_t   = U(end, :, 2) ./ (rho_t+1e-16);
v_t   = U(end, :, 3) ./ (rho_t+1e-16);
v_dot_n = u_t .* nx_top + v_t .* ny_top;
flux_top = rho_t .* v_dot_n .* A_top_vec;

% sum ONLY open part
mdot_top_loss = sum( flux_top(mask_open) );

% --- exit fluxes---
rho_ex = U(:, end, 1);
u_ex   = U(:, end, 2) ./ (rho_ex+1e-16);
mdot_exit = sum( rho_ex .* u_ex .* A_exit_vec );

% balance?
% -------------------------------------------------------------------------
Total_Mass_In = mdot_inlet - mdot_top_loss; 

fprintf('Flight Conditions: Mach %.2f\n', inputs.Mach);
fprintf('\nMass Balance:\n');
fprintf('  Inlet Face:       %8.4f kg/s\n', mdot_inlet);
fprintf('  Top Entrainment:  %8.4f kg/s\n', -mdot_top_loss);
fprintf('  ----------------------------\n');
fprintf('  TOTAL INFLOW:     %8.4f kg/s\n', Total_Mass_In);
fprintf('  TOTAL EXIT:       %8.4f kg/s\n', mdot_exit);

err = abs(mdot_exit - Total_Mass_In) / Total_Mass_In * 100;
fprintf('  ERROR:            %.2f%%\n', err);

%% propulsion
p_ex = (phys.gamma-1)*(U(:,end,4) - 0.5*rho_ex.*u_ex.^2);
p_in = (phys.gamma-1)*(U(:,1,4) - 0.5*rho_in.*u_in.^2);

term_exit = (rho_ex .* u_ex.^2 + p_ex) .* A_exit_vec;
term_inlet = (rho_in .* u_in.^2 + p_in) .* A_inlet_vec;

F_gross = sum(term_exit) - sum(term_inlet);

fprintf('\nPropulsion:\n');
fprintf('  Gross Thrust: %.2f N\n', F_gross);

% added heat
Q_vol_input = inputs.Heat_Release_Max;
Combustor_Start = GRID.L_inlet + GRID.L_isolator;
mask_burn = METRICS.Xc > Combustor_Start;
Q_total_Watts = sum(sum( mask_burn .* Q_vol_input .* METRICS.Vol ));

fprintf('  Total Heat Release:   %.2f MW\n', Q_total_Watts/1e6);

mdot_fuel = Q_total_Watts / 42e6; 
Isp = F_gross / (mdot_fuel * 9.81);

if F_gross > 0
    fprintf('  Specific Impulse (Isp): %.1f s\n', Isp);
else
    fprintf('  STATUS: DRAG DOMINATED\n');
end