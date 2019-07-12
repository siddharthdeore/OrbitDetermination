%% Compute topocentric RA and DEC from Azimuth Elevation and station cordinates
% Inputs
%   az - Azimuth [rad]
%   el - Elevation [rad]
%   stationLat - site/station Lattitude [rad]
%   stationLon - site/station Longitude [rad]
%   lst - local sidereal time [rad]

% Outputs
% author : Siddharth Deore

function [RA,DEC] = azel2radec(az,el,stationLat,stationLon,lst)
DEC = asin(sin(el)*sin(stationLat)+cos(el)*cos(stationLat)*cos(az));
local_hour_angle = atan2(-sin(az)*cos(el)/cos(DEC),(sin(el)-sin(DEC)*sin(stationLat))/(cos(DEC)*cos(stationLat)));
RA=mod(lst-local_hour_angle,2*pi);
end