function [METRICS] = Grid_Metrics(GRID)
    X = GRID.X;
    Y = GRID.Y; % rad
    [JMax_Nodes, IMax_Nodes] = size(X);

    % cell centers
    Xc = 0.25*(X(1:end-1, 1:end-1)+X(2:end, 1:end-1)+X(1:end-1, 2:end)+X(2:end, 2:end));
    Yc = 0.25*(Y(1:end-1, 1:end-1)+Y(2:end, 1:end-1)+Y(1:end-1, 2:end)+Y(2:end, 2:end));
    
    % ---vertical ---
    x_i_bot = X(1:end-1, 2:end-1); y_i_bot = Y(1:end-1, 2:end-1);
    x_i_top = X(2:end,   2:end-1); y_i_top = Y(2:end,   2:end-1);
    
    dx_i = x_i_top - x_i_bot;
    dy_i = y_i_top - y_i_bot;
    S_i_2D = sqrt(dx_i.^2 + dy_i.^2);
    Y_face_i = 0.5 * (y_i_bot + y_i_top); 
    S_i_3D = S_i_2D .* (2 * pi * Y_face_i); % rotational area
    
    nx_i =  dy_i ./ S_i_2D;
    ny_i = -dx_i ./ S_i_2D;
    
    % --- horizontal ---
    x_j_left = X(2:end-1, 1:end-1); y_j_left = Y(2:end-1, 1:end-1);
    x_j_right= X(2:end-1, 2:end);   y_j_right= Y(2:end-1, 2:end);
    
    dx_j = x_j_right - x_j_left;
    dy_j = y_j_right - y_j_left;
    S_j_2D = sqrt(dx_j.^2 + dy_j.^2);
    Y_face_j = 0.5 * (y_j_left + y_j_right);
    S_j_3D = S_j_2D .* (2 * pi * Y_face_j);
    
    nx_j = -dy_j ./ S_j_2D;
    ny_j =  dx_j ./ S_j_2D;

    %pappus
    x1=X(1:end-1,1:end-1); y1=Y(1:end-1,1:end-1);
    x2=X(1:end-1,2:end);   y2=Y(1:end-1,2:end);
    x3=X(2:end,2:end);     y3=Y(2:end,2:end);
    x4=X(2:end,1:end-1);   y4=Y(2:end,1:end-1);
    
    Area_2D = 0.5 * abs((x1.*y2 - y1.*x2) + (x2.*y3 - y2.*x3) + ...
                        (x3.*y4 - y3.*x4) + (x4.*y1 - y4.*x1));
                    
    Vol_3D = Area_2D .* (2 * pi * Yc); % rotational volume
    
    dx_cell = abs(0.5*((x2+x3)-(x1+x4)));
    dr_cell = abs(0.5*((y3+y4)-(y1+y2)));

    METRICS.Xc = Xc; METRICS.Yc = Yc;
    METRICS.Vol = Vol_3D; 
    METRICS.dx = dx_cell; METRICS.dr = dr_cell;
    
    METRICS.Si = S_i_3D; METRICS.nxi = nx_i; METRICS.nyi = ny_i;
    METRICS.Sj = S_j_3D; METRICS.nxj = nx_j; METRICS.nyj = ny_j;
end