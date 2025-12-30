%SOLVER_MAIN 
%runs the Ramjet simulation with dynmic plotting.

clear; clc;
addpath('./'); 

%% Init
fprintf('Initializing Simulation...\n');
Ramjet_Config;                  
GRID = Grid_Gen_Combined(inputs); 
METRICS = Grid_Metrics(GRID);
U = Initialize_Flow(inputs, METRICS, IC);

%% Params
Target_CFL = 0.4;       
CFL = 0.05;             

MaxIter = 2500;         
SaveInterval = 10;      
Residual_History = [];

figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 1000 800]);

fprintf('Starting Solver Loop (%d iterations)...\n', MaxIter);

%% time method
for iter = 1:MaxIter
    
    % CFL ramping
    if iter < 800
        CFL = 0.05 + (Target_CFL - 0.05) * (iter / 800);
    else
        CFL = Target_CFL;
    end
    
    % --- RK2---
    dt = Compute_TimeStep(U, METRICS, phys, CFL);
    U_n = U;
    
    % predictor
    [UL_i, UR_i] = Reconstruct_MUSCL(U_n, phys, 'i');
    [UL_j, UR_j] = Reconstruct_MUSCL(U_n, phys, 'j');
    F_i = Flux_Roe(UL_i, UR_i, METRICS.nxi, METRICS.nyi, phys);
    F_j = Flux_Roe(UL_j, UR_j, METRICS.nxj, METRICS.nyj, phys);
    
    dF_i = zeros(size(U_n));
    dF_i(:, 2:end-1, :) = (F_i(:, 1:end-1, :) .* METRICS.Si(:, 1:end-1) - ...
                           F_i(:, 2:end, :)   .* METRICS.Si(:, 2:end));              
    dF_j = zeros(size(U_n));
    dF_j(2:end-1, :, :) = (F_j(1:end-1, :, :) .* METRICS.Sj(1:end-1, :) - ...
                           F_j(2:end, :, :)   .* METRICS.Sj(2:end, :));
    
    S_total = Compute_Sources(U_n, METRICS, phys, inputs, GRID, iter);
    
    RHS = zeros(size(U_n));
    for k=1:4
        RHS(:,:,k) = (dF_i(:,:,k) + dF_j(:,:,k)) ./ METRICS.Vol + S_total(:,:,k);
    end
    U_star = U_n + dt * RHS;
    U_star = Boundary_Conditions(U_star, inputs, IC, METRICS);
    U_star = Enforce_Positivity(U_star, phys);

    % corrector
    [UL_i, UR_i] = Reconstruct_MUSCL(U_star, phys, 'i');
    [UL_j, UR_j] = Reconstruct_MUSCL(U_star, phys, 'j');
    F_i = Flux_Roe(UL_i, UR_i, METRICS.nxi, METRICS.nyi, phys);
    F_j = Flux_Roe(UL_j, UR_j, METRICS.nxj, METRICS.nyj, phys);
    
    dF_i = zeros(size(U_n));
    dF_i(:, 2:end-1, :) = (F_i(:, 1:end-1, :) .* METRICS.Si(:, 1:end-1) - ...
                           F_i(:, 2:end, :)   .* METRICS.Si(:, 2:end));       
    dF_j = zeros(size(U_n));
    dF_j(2:end-1, :, :) = (F_j(1:end-1, :, :) .* METRICS.Sj(1:end-1, :) - ...
                           F_j(2:end, :, :)   .* METRICS.Sj(2:end, :));
    
    S_total = Compute_Sources(U_star, METRICS, phys, inputs, GRID, iter);
    
    RHS_star = zeros(size(U_n));
    for k=1:4
        RHS_star(:,:,k) = (dF_i(:,:,k) + dF_j(:,:,k)) ./ METRICS.Vol + S_total(:,:,k);
    end
    U = 0.5 * (U_n + U_star + dt * RHS_star);
    U = Boundary_Conditions(U, inputs, IC, METRICS);
    U = Enforce_Positivity(U, phys);
    
    % ---check---
    if any(isnan(U(:)))
        fprintf('\n!!! CRITICAL ERROR: NaNs detected at Iter %d !!!\n', iter);
        U = U_n; break;
    end
    dRho = rms(rms(real(U(:,:,1) - U_n(:,:,1))));
    Residual_History = [Residual_History; dRho];
    
    %% --- Viz---
    if mod(iter, SaveInterval) == 0
        
        ramp = (iter - inputs.Burn_Start_Iter) / inputs.Burn_Ramp_Length;
        heat_pct = max(0, min(100, ramp*100));
        
        fprintf('Iter: %4d | CFL: %.2f | Heat: %5.1f%% | dRho: %.2e\n', ...
                iter, CFL, heat_pct, dRho);
        
        % density
        rho_plot = real(U(:,:,1));
        
        subplot(2,1,1);
        pcolor(METRICS.Xc, METRICS.Yc, rho_plot);
        shading interp; 
        
        colormap(jet); colorbar; axis equal tight;
        caxis([0.0, 0.5]);
        title(sprintf('Density (Iter %d) [Smooth]', iter));
        xlabel('Axial (x)'); ylabel('Radial (r)');
        
        %residuals
        subplot(2,1,2);
        semilogy(Residual_History, 'LineWidth', 1.5);
        grid on; title('Convergence');
        xlabel('Iteration'); ylabel('dRho');
        xlim([0, MaxIter]);
        
        drawnow;
    end
end

fprintf('Simulation Complete.\n');

function U = Enforce_Positivity(U, phys)
    min_rho = 1e-4; min_p = 1e-4;
    rho = U(:,:,1);
    u = U(:,:,2) ./ (rho + 1e-16);
    v = U(:,:,3) ./ (rho + 1e-16);
    E = U(:,:,4);
    
    bad_rho = rho < min_rho;
    if any(bad_rho(:)), rho(bad_rho) = min_rho; u(bad_rho)=0; v(bad_rho)=0; end
    
    p = (phys.gamma - 1) * (E - 0.5 * rho .* (u.^2 + v.^2));
    bad_p = p < min_p;
    if any(bad_p(:))
        p(bad_p) = min_p;
        E(bad_p) = p(bad_p) ./ (phys.gamma - 1) + 0.5 * rho(bad_p) .* (u(bad_p).^2 + v(bad_p).^2);
    end
    U(:,:,1) = rho; U(:,:,2) = rho .* u; U(:,:,3) = rho .* v; U(:,:,4) = E;
end