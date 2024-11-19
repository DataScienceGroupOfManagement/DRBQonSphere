function [x_groups, yp_groups] = sdesign_data_generate(x0, q, fun)

% rotation matrix 
A = @(theta) [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

x_groups = cell(1,q);
yp_groups = cell(1,q);

if q == 1
    x_groups{1} = x0;
    yp_groups{1} = fun(x0);
else
    for i = 1:q
        theta_i = pi*i/q;
        x_i = (A(theta_i)*x0')'; % the rotated x_i,   Nx3
        yp_i = fun(x_i); % Nx1
        x_groups{i} = x_i;
        yp_groups{i} = yp_i;
    end
end

