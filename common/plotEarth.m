%% Plots 3d Earth
% dependencies
%   earthmap.jpg
%--------------------------------------------------------%
% Author : Siddharth Deore
% email  : siddharthdeore@gmail.com
%%

% Plot Earth with texture
function plotEarth(gmst)
% Turn off the normal axes

%set(gca, 'NextPlot','add', 'Visible','off');
axis equal; axis auto; axis vis3d;
image_file='earthmap.jpg';

GMST0 = gmst;
radius = 6371; %km
[xs,ys,zs]=sphere(32);

xs = xs*radius;
ys = ys*radius;
zs = zs*radius;

globe=surf(xs,ys,-zs,'FaceAlpha',1.0,'EdgeColor',0.5*[1 1 1],'FaceColor', 'none','FaceLighting','flat');
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end
cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');

end
