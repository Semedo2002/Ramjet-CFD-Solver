function [dt] = Compute_TimeStep(U, METRICS, phys, CFL)
% TIMESTEP
% maximum stable time step using cd.

    gamma = phys.gamma;
    
    rho = U(:,:,1);
    u   = U(:,:,2) ./ rho;
    v   = U(:,:,3) ./ rho;
    p   = (gamma - 1) * (U(:,:,4) - 0.5 * rho .* (u.^2 + v.^2));
    
    % Sound speed
    a = sqrt(gamma * p ./ rho);
 
    % characteristic speeds
    lambda_x = abs(u) + a;
    lambda_r = abs(v) + a;
    
    % inv ts per cell
    inv_dt = (lambda_x ./ METRICS.dx) + (lambda_r ./ METRICS.dr);
    inv_dt(inv_dt < 1e-6) = 1e-6;
    
    % local ts
    dt_local = CFL ./ inv_dt;
    
    % global ts
    dt = min(dt_local(:));

end