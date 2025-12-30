function [S] = Compute_Sources(U, METRICS, phys, inputs, GRID, iter)
%only pressure source for radial mom + Heat Addition.
%!!!!!!!!!!!:continuity& axial momentum sources HAVE TO BE 0.

    gamma = phys.gamma;
    Xc = METRICS.Xc;
    
    rho = U(:,:,1);
    u   = U(:,:,2) ./ (rho + 1e-16);
    v   = U(:,:,3) ./ (rho + 1e-16);
    E   = U(:,:,4);
    p   = (gamma - 1) * (E - 0.5 * rho .* (u.^2 + v.^2));
    
    S = zeros(size(U));
    
    %radial momentum ONLY
    inv_r = 1 ./ (METRICS.Yc + 1e-6);
    S(:,:,3) = p .* inv_r; 
    
    % heat
    current_ramp = (iter - inputs.Burn_Start_Iter) / inputs.Burn_Ramp_Length;
    heat_scale = max(0, min(1.0, current_ramp));
    
    if heat_scale > 0
        Combustor_Start = GRID.L_inlet + GRID.L_isolator;
        mask_burn = Xc > Combustor_Start;
        Q_applied = inputs.Heat_Release_Max * heat_scale;
        
        S(:,:,4) = S(:,:,4) + (mask_burn * Q_applied);
    end
end