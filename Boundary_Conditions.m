function [U] = Boundary_Conditions(U, inputs, IC, METRICS)
%BOUNDARY CONDITIONS
%physics at the domain edges handles  cowl lip transition from Freestream to wall.

    [JMax, IMax, ~] = size(U);
    
    % i = 1 -> hold fixed at IC
    % ------------------------------------
    U(:, 1, 1) = IC.U1;
    U(:, 1, 2) = IC.U2;
    U(:, 1, 3) = IC.U3;
    U(:, 1, 4) = IC.U4;
    
    %i = IMax -> Extrapolation s.sonic exit
    % -------------------------------------------------------
    U(:, end, :) = U(:, end-1, :);
    
    %j = 1 -> Ramp / centerbody wall
    % ----------------------------------------------------
    % tangency cond to inviscid wall
    rho_bot = U(1, :, 1);
    u_bot   = U(1, :, 2) ./ rho_bot;
    v_bot   = U(1, :, 3) ./ rho_bot;

    x_c = METRICS.Xc(1, :);
    theta_inlet = deg2rad(inputs.Theta_deg);
    theta_combustor = deg2rad(2.0); 

    L_inlet = 1.2; 
    L_isolator = 0.5;
    
    wall_angle_bot = zeros(size(x_c));
    mask_ramp = x_c <= L_inlet;
    mask_iso  = (x_c > L_inlet) & (x_c <= L_inlet + L_isolator);
    mask_comb = x_c > L_inlet + L_isolator;
    
    wall_angle_bot(mask_ramp) = theta_inlet;
    wall_angle_bot(mask_iso)  = 0.0;
    wall_angle_bot(mask_comb) = -theta_combustor;%drop -> -ve
    
    %tangency bot
    Vel_bot = sqrt(u_bot.^2 + v_bot.^2);
    u_new_b = Vel_bot .* cos(wall_angle_bot);
    v_new_b = Vel_bot .* sin(wall_angle_bot);
    
    U(1, :, 2) = rho_bot .* u_new_b;
    U(1, :, 3) = rho_bot .* v_new_b;
    
    % j = JMax fs & cowl wall
    % ------------------------------------------------------------
    rho_top = U(end, :, 1);
    u_top   = U(end, :, 2) ./ rho_top;
    v_top   = U(end, :, 3) ./ rho_top;
    x_c_top = METRICS.Xc(end, :);
    
    % mask cowl
    mask_cowl = x_c_top >= L_inlet;
    mask_open = x_c_top < L_inlet;
    
    % --- dirichlet fs ---
    U(end, mask_open, 1) = IC.U1;
    U(end, mask_open, 2) = IC.U2;
    U(end, mask_open, 3) = IC.U3;
    U(end, mask_open, 4) = IC.U4;
    
    % ---cowl wall---

    if any(mask_cowl)
        u_cowl = u_top(mask_cowl);
        v_cowl = v_top(mask_cowl);
        rho_cowl = rho_top(mask_cowl);

        v_cowl = 0;
        
        U(end, mask_cowl, 3) = rho_cowl .* v_cowl;
        % slip wall , axial momentum n.c
    end

end