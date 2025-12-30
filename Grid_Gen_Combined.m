function [GRID] = Grid_Gen_Combined(inputs)
% large domain for shock capture to prevent drag orcrashes.

    IMax = inputs.IMax; 
    JMax = inputs.JMax;
    theta_inlet = deg2rad(inputs.Theta_deg);
    
    L_inlet = 1.2;       
    L_isolator = 0.5;    
    L_combustor = 1.0;   
    
    H_ramp_throat = L_inlet * tan(theta_inlet); 
    theta_combustor = deg2rad(2.0); 
    
    % --- stabilility ---
    H_outer_start = 1.5; 

    Total_L = L_inlet + L_isolator + L_combustor;
    N1 = floor(IMax * (L_inlet/Total_L));
    N2 = floor(IMax * (L_isolator/Total_L));
    N3 = IMax - N1 - N2;
    
    x1 = linspace(0, L_inlet, N1);
    x2 = linspace(L_inlet, L_inlet + L_isolator, N2+1);
    x3 = linspace(L_inlet + L_isolator, Total_L, N3+1);
    x_1d = [x1, x2(2:end), x3(2:end)];
    IMax = length(x_1d);
    
    X = zeros(JMax, IMax);
    Y = zeros(JMax, IMax); 

    for i = 1:IMax
        x_loc = x_1d(i);
        
        % bottom wall
        if x_loc <= L_inlet
            y_bot = x_loc * tan(theta_inlet);
        elseif x_loc <= (L_inlet + L_isolator)
            y_bot = H_ramp_throat; 
        else
            dx_comb = x_loc - (L_inlet + L_isolator);
            y_bot = H_ramp_throat - dx_comb * tan(theta_combustor); 
        end
        
        %top
        if x_loc < L_inlet
             y_top = H_outer_start + x_loc * 0.2; 
        else
             y_top = H_outer_start + L_inlet * 0.2; 
        end

        y_col = linspace(y_bot, y_top, JMax);
        X(:, i) = x_loc;
        Y(:, i) = y_col;
    end

    GRID.X = X; GRID.Y = Y;
    GRID.IMax = IMax; GRID.JMax = JMax;
    GRID.L_inlet = L_inlet; GRID.L_isolator = L_isolator;
    
    fprintf('Grid Restored to Stable Configuration (High Ceiling).\n');
end