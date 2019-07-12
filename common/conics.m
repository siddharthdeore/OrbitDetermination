%% returns 3d elipse cordinates for given kepler elements
% Inputs
%   a - semimajor axis [km]
%   e - eccentricity
%   i - inclination [rad]
%   RA - Right Accertion [rad]
%   omega - Argument of periapsis [rad]
%   M0 - Mean anomaly [rad]
%   GM - Gravity parameter [Km^3/s^2]
% Outputs
%   x - vector [km]
%   y - vector [km]
%   z - vector [km]
% dependencies
%   invKepler(M,e)
%	plotEarth
%--------------------------------------------------------%
% Author : Siddharth Deore
% email  : siddharthdeore@gmail.com
%%
function [x,y,z] = conics(a,e,i,RA,omega,M0,GM)
% Mean motion
n = sqrt(GM/(a*a*a));
period=round(2*pi*sqrt(a*a*a/GM));
Npoints = 256;
t=linspace(-period/2,period/2,Npoints);
M=n*t+M0;
if e == 1
    % parabola
    p=2*a;
    focus = a*e-a;
elseif e<1
    % ellipse
    p=a*(1-e*e);
    focus = a*e;
else
    % hyperbola
    p=a*(e*e-1);
    focus = 2*a-a*e;
end
E=M;
v=M;
for k=1:length(M)
[E(k),v(k)]=invKepler(M(k),e); % this line gives evenly spaced orbit. converts mean anomaly to Ecentric anomaly
end
%% periapsis
[xp,yp,zp] = getXYZ(a,e,i,0,0,omega,RA,p); % perigee cordinates

%% trajectory
[x,y,z] = getXYZ(a,e,i,E,v,omega,RA,p);
%% plot
plot3(x',y',z','-','LineWidth',1);
hold on
plot3(0,0,0,'r*'); % planet position (planet focus)

plot3(xp,yp,zp,'ko'); % Perigee postion
if(e>1)
    %% Apoapsis
    xa=-10*xp;
    ya=-10*yp;
    za=-10*zp;
else
    [xa,ya,za] = getXYZ(a,e,i,pi,pi,omega,RA,p); % apogee cordinates
    plot3([xp,xa],[yp,ya],[zp,za],'k-.');  % line of major axis
    plot3(xa+xp,ya+yp,za+zp,'k*'); % empty focus
end
hold on
axis equal;
end

%%
function [x,y,z] = getXYZ(a,e,i,E,v,omega,RA,p)
r=p./(1+e*cos(v));
%r = a*(1-e*cos(E));

theta=v+omega;
x1=r.*cos(theta);%+ 0*focus;
y1=r.*sin(theta);

x = x1 * cos(RA) - y1 * cos(i) * sin(RA); 
y = x1 * sin(RA) + y1 * cos(i) * cos(RA); 
z = y1 * sin(i);
end