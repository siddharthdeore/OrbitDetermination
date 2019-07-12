% 
% Input parameters:
% a: Semi major axis (meters)
% e: Eccentricity
% i: Inclination angle of the orbit [rad]
% RA: Right ascension [rad]
% omega: Argument of perigee [rad]
% M0: Mean anomaly at time t = 0 [rad]
% GM :  Gravitational Constant of planet [km^3/s^2]
% t: time vector for which the orbit should be computed
% 
% Output parameters:
% P = [x,y,z]: ECI position of the satellite
%  Dependencies
%  invKepler
function [P] = getKeplerEllipse(a, e, i, RA, omega, M0,GM, t)
% Mean motion
n = sqrt(GM/(a*a*a));
% Mean anomaly
M = M0 + n*t; 
% Eccentric anomaly evaluation using Kepler equation
[E,v] = invKepler(M, e); 

[x, y, z] = getXYZ(a,e,i,E,v,omega,RA);
P = [x', y', z'];

%figure(1),plot3(P(:,1),P(:,2),P(:,3)),hold on, disEarth, hold off, axis equal
orbit_name = strcat('a = ',num2str(a),'e = ',num2str(e),'i = ',num2str(i));

[peri_x,peri_y,peri_z] = getXYZ(a,e,i,0,0,omega,RA); % perigee cordinates
[apo_x,apo_y,apo_z] = getXYZ(a,e,i,pi,pi,omega,RA);  % apogee cordinates
plot3(P(:,1),P(:,2),P(:,3),'DisplayName',orbit_name);
hold on;
plot3(peri_x,peri_y,peri_z,'ro','DisplayName','');
plot3(apo_x,apo_y,apo_z,'bo','DisplayName','');
plot3([peri_x,apo_x],[peri_y,apo_y],[peri_z,apo_z],'--','Color',0.5*[1 1 1]);
axis equal; 
view(0,90);
end

function [x,y,z] = getXYZ(a,e,i,E,v,omega,RA)
% Actual distance of the satellite from the center of the earth
r = a*(1-e*cos(E));

% Rotation of the perigee
theta = v+omega;
% X and Y in the orbital plane
xp = r .* cos(theta); 
yp = r .* sin(theta);     

% Rotation of the orbital plane due to the inclination and right ascension
x = xp*cos(RA) - yp*cos(i)*sin(RA); 
y = xp*sin(RA) + yp*cos(i)*cos(RA); 
z = yp*sin(i); 
end