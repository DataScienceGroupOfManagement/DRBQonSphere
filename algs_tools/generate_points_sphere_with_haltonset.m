function [points_sphere] = generate_points_sphere_with_haltonset(n_points)
%GENERATE_POINTS_SPHERE_WITH_HALTONSET 
%==========================================================================
%
% Generate quadrature nodes by mapping points from the 2 dimensional
% halton set (in [0,1]X[0,1]) onto the sphere of radius Sphere_Radius.
%
%==========================================================================

% Generate the quadrature nodes
points_sphere=zeros(n_points,3);
p  = haltonset(2);
tz = 2*p(1:n_points,:)-1;
points_sphere(:,3)  = 1.*tz(:,2);
r  = sqrt(1.^2-points_sphere(:,3).^2);
points_sphere(:,1)  = r.*cos(pi*tz(:,1));
points_sphere(:,2)  = r.*sin(pi*tz(:,1));

end