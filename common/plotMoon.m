
% Plot Moon with texture
function plotMoon()
% Turn off the normal axes

set(gca, 'NextPlot','add', 'Visible','off');
axis equal; axis auto; axis vis3d;
image_file='moonmap.jpg';

GMST0 = 0.1;
radius = 1737.1; %km
[xs,ys,zs]=sphere(32);

xs = xs*radius;
ys = ys*radius;
zs = zs*radius;

globe=surf(xs,ys,zs,'FaceAlpha',1.0,'EdgeColor',0.5*[1 1 1],'FaceColor', 'none','FaceLighting','flat');
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end
cdata = imread(image_file);
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
end
