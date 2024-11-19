function show_s2_sphere(N, co)
% Draws a surf plot of an n-by-n sphere with the color of 'co'
%
% Syntax
% show_s2_sphere

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.


[X,Y,Z] = sphere(N);
CO(:,:,1) = co(1) * ones(N+1); % red
CO(:,:,2) = co(2) * ones(N+1); % green
CO(:,:,3) = co(3) * ones(N+1); % blue

surf(X,Y,Z,CO,...
     'FaceColor','interp','FaceLighting','phong',...
     'EdgeColor','none','SpecularStrength',0)
 
light('Position', [3 1 3]);

axis equal
grid off
axis off
