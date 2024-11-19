function show_s2_points(points_x, co, r)

% SHOW_S2_POINT_SET 3D illustration of a point set
% points_x: each column is a point
% co: color for points
% r: the size value of point, i.e., the radius of points


N = size(points_x,2);
surf_jet;
[X, Y, Z] = sphere;
rX = r * X; 
rY = r * Y; 
rZ = r * Z;

CO(:,:,1) = co(1) * ones(21); % red
CO(:,:,2) = co(2) * ones(21); % green
CO(:,:,3) = co(3) * ones(21); % blue

for n = 1:N
   surf(points_x(1,n)+rX, points_x(2,n)+rY, points_x(3,n)+rZ,CO, ...
   'FaceColor','interp','FaceLighting','phong','EdgeColor','none')
end
%
axis equal
axis off
grid off
hold off