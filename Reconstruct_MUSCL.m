function [L_State, R_State] = Reconstruct_MUSCL(U, phys, direction)
%===============================================================
% RECON MUSCL cell centered interpolation to cell interfaces by 2nd order MUSCL with VL limiting.
%
% Inputs:
%   U:[JMax, IMax, 4] 
%  'i' (axial) or 'j' (radial) 
%=============================================================
    gamma = phys.gamma;
    [JMax, IMax, ~] = size(U);
   
    % ---------------------------------------------------------------------
    W = zeros(size(U));
    
    rho = U(:,:,1);
    u   = U(:,:,2) ./ rho;
    v   = U(:,:,3) ./ rho;
    E   = U(:,:,4);
    p   = (gamma - 1) * (E - 0.5 * rho .* (u.^2 + v.^2));
    
    W(:,:,1) = rho;
    W(:,:,2) = u;
    W(:,:,3) = v;
    W(:,:,4) = p;

    % recon loop
    % ---------------------------------------------------------------------
    
    if strcmp(direction, 'i')
        % ---axial recon---

        WL = zeros(JMax, IMax-1, 4); 
        WR = zeros(JMax, IMax-1, 4);
        
        for k = 1:4
            Q = W(:,:,k); 
            
            % pad
            Q_pad = [Q(:,1), Q, Q(:,end)];
            
            % diff
            d_minus = Q_pad(:, 2:end-1) - Q_pad(:, 1:end-2); 
            d_plus  = Q_pad(:, 3:end)   - Q_pad(:, 2:end-1); 
            
            epsilon = 1e-12; 
            
            % ---left at i+1/2 ---
  
            delta_i_minus = d_minus(:, 1:end-1);
            delta_i_plus  = d_plus(:,  1:end-1);
            
            r_i = delta_i_minus ./ (delta_i_plus + epsilon);
            phi_i = Slope_Limiter(r_i);

            WL(:,:,k) = Q(:, 1:end-1) + 0.5 * phi_i .* delta_i_plus;
            
            % ---right at i+1/2 ---
            % cells i+1 neeeded -> iterate cells 2 to IMax.

            delta_ip1_minus = d_minus(:, 2:end);
            delta_ip1_plus  = d_plus(:,  2:end);
            
            r_ip1 = delta_ip1_plus ./ (delta_ip1_minus + epsilon);
            phi_ip1 = Slope_Limiter(r_ip1);

            WR(:,:,k) = Q(:, 2:end) - 0.5 * phi_ip1 .* delta_ip1_minus;
            
        end
        
    elseif strcmp(direction, 'j')
        % ---radial recon ---
       
        WL = zeros(JMax-1, IMax, 4); 
        WR = zeros(JMax-1, IMax, 4);
        
        for k = 1:4
            Q = W(:,:,k);
            Q_pad = [Q(1,:); Q; Q(end,:)];
            
            % diff
            d_minus = Q_pad(2:end-1, :) - Q_pad(1:end-2, :);
            d_plus  = Q_pad(3:end, :)   - Q_pad(2:end-1, :);
            
            epsilon = 1e-12;
            
            % --- j = 1 -> JMax-1) ---
            delta_j_minus = d_minus(1:end-1, :);
            delta_j_plus  = d_plus(1:end-1, :);
            
            r_j = delta_j_minus ./ (delta_j_plus + epsilon);
            phi_j = Slope_Limiter(r_j);
            
            WL(:,:,k) = Q(1:end-1, :) + 0.5 * phi_j .* delta_j_plus;
            
            % ---j+1 = 2 -> JMax---
            delta_jp1_minus = d_minus(2:end, :);
            delta_jp1_plus  = d_plus(2:end, :);
            
            r_jp1 = delta_jp1_plus ./ (delta_jp1_minus + epsilon);
            phi_jp1 = Slope_Limiter(r_jp1);
            
            WR(:,:,k) = Q(2:end, :) - 0.5 * phi_jp1 .* delta_jp1_minus;
        end
        
    else
        error('Direction must be i or j');
    end
    % ---------------------------------------------------------------------
    L_State = PrimToCons(WL, gamma);
    R_State = PrimToCons(WR, gamma);

end

function phi = Slope_Limiter(r)
    % VL limiter
    r_abs = abs(r);
    phi = (r + r_abs) ./ (1 + r_abs);
end

function U = PrimToCons(W, gamma)
    rho = W(:,:,1);
    u   = W(:,:,2);
    v   = W(:,:,3);
    p   = W(:,:,4);
    
    U = zeros(size(W));
    U(:,:,1) = rho;
    U(:,:,2) = rho .* u;
    U(:,:,3) = rho .* v;
    U(:,:,4) = p ./ (gamma - 1) + 0.5 * rho .* (u.^2 + v.^2);
end