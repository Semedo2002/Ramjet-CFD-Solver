function [U] = Initialize_Flow(inputs, METRICS, IC)
% INITIALIZE_FLOW
% sets ic for the  domain, sets everything to freestream atm

    [JMax, IMax] = size(METRICS.Vol);
    
    U = zeros(JMax, IMax, 4);
    
    % fs IC
    U(:,:,1) = IC.U1;
    U(:,:,2) = IC.U2;
    U(:,:,3) = IC.U3;
    U(:,:,4) = IC.U4;
    
    fprintf('Flow Field Initialized to Mach %.2f Freestream.\n', inputs.Mach);
end