%orbit radius velocity latitude longitude ECEF
% Richard Rieberexit

% October 17, 2006
% rrieber@gmail.com
%
% [Lat, Long] = RVtoLatLong(ECEF)
%
% Revision 8/21/07: Added H1 line for lookfor functionality.
%
% Revision 9/25/07 - Fixed typo referring to velocity in comments.
%                    Velocity not needed in this function.
%
% Revision 9/20/09: Removed GMST since it isn't used
%
% Purpose:  This fuction convertes ECEF coordinates to Geocentric latitude
%           and longitude given ECEF radius in km. Valid for any planetary
%           body.
%
% Inputs:  o ECI  - A 3x1 vector of Earth-Centered Inertial (IJK)
%                   coordinates in km.
%
% Outputs: o Lat  - Geocentic latitude of spacecraft in radians
%          o Long - Longitude of spacecraft in radians

function [Lat, Long] = RVtoLatLong(ECEF)

if nargin < 1
    error('Too few inputs.  See help RVtoLatLong')
elseif nargin > 1
    error('Too many inputs.  See help RVtoLatLong')
end

if length(ECEF) ~= 3
    error('ECI length incorrect.  See help RVtoLatLong')
end
if nargout ~= 2
    error('Incorrect number of outputs.  See help RVtoLatLong')
end

r_delta = norm(ECEF(1:2));

sinA = ECEF(2)/r_delta;
cosA = ECEF(1)/r_delta;

Long = atan2(sinA,cosA);

if Long < -pi
    Long = Long + 2*pi;
end

Lat = asin(ECEF(3)/norm(ECEF));
end