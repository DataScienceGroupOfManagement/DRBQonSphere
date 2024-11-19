function [fig1,fig2] = pltfunc(F, X, lengendticks, tstr, ifig, scale, ipr)
% [K, fh1, fh2, fig1,fig2] = pltfunc(F, X, tstr, ifig, scale, ipr)
% Plot the function F at the m points X on the unit sphere S^2 in R^3
% Title string tstr using figure windows ifig and ifig+1 and scaling scale
% are all optional

% Default arguments
if nargin < 7
    ipr = 0;
end
if nargin < 6
    scale = 0.5;
end
if nargin < 5
    ifig = 1;
end
if nargin < 4
    tstr = [];
end
if nargin < 3
    lengendticks = [];
end

K = convhulln(X');

X1=[X(1,K(:,1)); X(1,K(:,2)); X(1,K(:,3))];
Y1=[X(2,K(:,1)); X(2,K(:,2)); X(2,K(:,3))];
Z1=[X(3,K(:,1)); X(3,K(:,2)); X(3,K(:,3))];
C = F(K');

if ifig > 0
    fig1 = figure(ifig); clf;
    fh1 = patch(X1, Y1, Z1, C);
    colormap(jet(255));
    view(90,0);
    axis vis3d
    axis equal tight
    view([1 1 1]);
    grid off
    set(gca, 'Visible', 'off')
    cbh1 = colorbar('Location','EastOutside');
    cbh1.Position(1) = 0.85;
    if ipr > 10
        hold on
        plot3(X(1,:), X(2,:), X(3,:), 'k.', 'MarkerSize', 16);
        hold off
    else
        set(fh1, 'EdgeColor', 'none')
    end
    title(tstr);
    % view(165,-10);
    view(45,135);
end

[Fmax, ~] = max(F);
[Fmin, ~] = min(F);

fprintf('Minimum function value = %.6f, Maximum function value = %.6f\n', Fmin, Fmax);

FS = 1 + (scale/(Fmax-Fmin))*(C-Fmin);

figureUnits = 'centimeters';
figureWidth = 12;
figureHeight = 14;

% figures 
fig2 = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight])

% fig2 = figure(abs(ifig)+1);
clf;
fh2 = patch(X1.*FS, Y1.*FS, Z1.*FS, C);

colormap(jet(1000));
view(90,0);
axis vis3d
axis equal tight
if isempty(lengendticks)
    colorbar('Location','north','FontSize',18, ...
        'Position', [0.136563876651982 0.886170764602137 0.759911894273128 0.0444591416813621]);
else
    colorbar('Location','north','FontSize',18, 'Ticks', lengendticks, ...
        'Position', [0.136563876651982 0.886170764602137 0.759911894273128 0.0444591416813621]);
end

% cbh2.Position(1) = 0.85;
set(fh2, 'EdgeColor', 'none');

axis off
title(tstr);
% view(165,-10);
view(45,135);


