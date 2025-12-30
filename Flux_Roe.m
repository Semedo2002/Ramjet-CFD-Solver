function [F] = Flux_Roe(UL, UR, nx, ny, phys)
% ROE FLUXES ,numerical fluxes across  interface by roe FDS, handles [J, I, 4] structured inputs by reshaping.
    gamma = phys.gamma;
    
    % dimensions & flatten
    % ---------------------------------------------------------------------
    [Rows, Cols, Vars] = size(UL);
    NumCells = Rows * Cols;
    
    % flatten 
    UL_flat = reshape(UL, NumCells, 4);
    UR_flat = reshape(UR, NumCells, 4);
    
    % flatten normals to [N x 1]
    nx_flat = reshape(nx, NumCells, 1);
    ny_flat = reshape(ny, NumCells, 1);
    
    %flux Output [N x 4]
    F_flat = zeros(NumCells, 4);

    % ---------------------------------------------------------------------
    % ls
    rhoL = UL_flat(:,1);
    uL   = UL_flat(:,2) ./ (rhoL + 1e-16);
    vL   = UL_flat(:,3) ./ (rhoL + 1e-16);
    E_L  = UL_flat(:,4);
    pL   = (gamma - 1) * (E_L - 0.5 * rhoL .* (uL.^2 + vL.^2));
    HL   = (E_L + pL) ./ (rhoL + 1e-16);
    
    % rs
    rhoR = UR_flat(:,1);
    uR   = UR_flat(:,2) ./ (rhoR + 1e-16);
    vR   = UR_flat(:,3) ./ (rhoR + 1e-16);
    E_R  = UR_flat(:,4);
    pR   = (gamma - 1) * (E_R - 0.5 * rhoR .* (uR.^2 + vR.^2));
    HR   = (E_R + pR) ./ (rhoR + 1e-16);
    
    % roe avgs
    % ---------------------------------------------------------------------
    R = sqrt(rhoR ./ (rhoL + 1e-16));
    
    rho_roe = R .* rhoL; 
    u_roe   = (R .* uR + uL) ./ (R + 1);
    v_roe   = (R .* vR + vL) ./ (R + 1);
    H_roe   = (R .* HR + HL) ./ (R + 1);
    
    V_sq_roe = u_roe.^2 + v_roe.^2;
    a_roe    = sqrt((gamma - 1) * (H_roe - 0.5 * V_sq_roe));
    
    % projected vel
    U_cont_roe = u_roe .* nx_flat + v_roe .* ny_flat;
    
    % eigenvalues
    % ---------------------------------------------------------------------
    lambda1 = U_cont_roe - a_roe;
    lambda2 = U_cont_roe;
    lambda3 = U_cont_roe;
    lambda4 = U_cont_roe + a_roe;
    
    % harten
    epsilon = 0.1 * a_roe; 
    lambda1 = abs(lambda1);
    idx1 = lambda1 < epsilon;
    lambda1(idx1) = 0.5 * (lambda1(idx1).^2 ./ epsilon(idx1) + epsilon(idx1));
    
    lambda4 = abs(lambda4);
    idx4 = lambda4 < epsilon;
    lambda4(idx4) = 0.5 * (lambda4(idx4).^2 ./ epsilon(idx4) + epsilon(idx4));
    
    lambda2 = abs(lambda2);
    lambda3 = abs(lambda3);
    
    % 5. wwave strength
    % ---------------------------------------------------------------------
    dp = pR - pL;
    drho = rhoR - rhoL;
    du = uR - uL;
    dv = vR - vL;
    
    dU_perp = du .* nx_flat + dv .* ny_flat;
    
    alpha4 = 0.5 * (dp ./ (rho_roe .* a_roe) + dU_perp) ./ a_roe;
    alpha1 = 0.5 * (dp ./ (rho_roe .* a_roe) - dU_perp) ./ a_roe;
    alpha2 = drho - (dp ./ a_roe.^2);
    alpha3 = du .* ny_flat - dv .* nx_flat; 
    
    % 6. dissipation
    % ---------------------------------------------------------------------
    % u-a
    D1_1 = alpha1 .* lambda1;
    D1_2 = alpha1 .* lambda1 .* (u_roe - a_roe .* nx_flat);
    D1_3 = alpha1 .* lambda1 .* (v_roe - a_roe .* ny_flat);
    D1_4 = alpha1 .* lambda1 .* (H_roe - u_roe .* a_roe .* nx_flat - v_roe .* a_roe .* ny_flat);
    
    % entropy
    D2_1 = alpha2 .* lambda2;
    D2_2 = alpha2 .* lambda2 .* u_roe;
    D2_3 = alpha2 .* lambda2 .* v_roe;
    D2_4 = alpha2 .* lambda2 .* 0.5 .* V_sq_roe;
    
    % shear
    D3_1 = zeros(NumCells, 1);
    D3_2 = alpha3 .* lambda3 .* ny_flat;
    D3_3 = alpha3 .* lambda3 .* (-nx_flat);
    D3_4 = alpha3 .* lambda3 .* (u_roe .* ny_flat - v_roe .* nx_flat);
    
    % u+a
    D4_1 = alpha4 .* lambda4;
    D4_2 = alpha4 .* lambda4 .* (u_roe + a_roe .* nx_flat);
    D4_3 = alpha4 .* lambda4 .* (v_roe + a_roe .* ny_flat);
    D4_4 = alpha4 .* lambda4 .* (H_roe + u_roe .* a_roe .* nx_flat + v_roe .* a_roe .* ny_flat);

    % total
    Diss_1 = D1_1 + D2_1 + D3_1 + D4_1;
    Diss_2 = D1_2 + D2_2 + D3_2 + D4_2;
    Diss_3 = D1_3 + D2_3 + D3_3 + D4_3;
    Diss_4 = D1_4 + D2_4 + D3_4 + D4_4;
    
    % physical fluxes
    % ---------------------------------------------------------------------
    U_perp_L = uL .* nx_flat + vL .* ny_flat;
    FL_1 = rhoL .* U_perp_L;
    FL_2 = rhoL .* uL .* U_perp_L + pL .* nx_flat;
    FL_3 = rhoL .* vL .* U_perp_L + pL .* ny_flat;
    FL_4 = rhoL .* HL .* U_perp_L;
    
    U_perp_R = uR .* nx_flat + vR .* ny_flat;
    FR_1 = rhoR .* U_perp_R;
    FR_2 = rhoR .* uR .* U_perp_R + pR .* nx_flat;
    FR_3 = rhoR .* vR .* U_perp_R + pR .* ny_flat;
    FR_4 = rhoR .* HR .* U_perp_R;
    
    % end roe flux
    % ---------------------------------------------------------------------
    F_flat(:,1) = 0.5 * (FL_1 + FR_1 - Diss_1);
    F_flat(:,2) = 0.5 * (FL_2 + FR_2 - Diss_2);
    F_flat(:,3) = 0.5 * (FL_3 + FR_3 - Diss_3);
    F_flat(:,4) = 0.5 * (FL_4 + FR_4 - Diss_4);
    % ---------------------------------------------------------------------
    F = reshape(F_flat, Rows, Cols, 4);

end