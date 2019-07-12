function LLA = ECEFm2LLA(ECEF)
% ECI2LLA.m
% usage: ENU = ECEF2LLA(ECEF)
% Inputs:
% 			ECEF is a matrix of row vectors of form [ time x y z]
% 					x, y, and z are in meters
%              time in seconds from GMT (Zulu time)
% Output:
%			LLA is a matrix of row vectors of form [time lat lon HAE] 
% 				lat, lon are in degrees
%           HAE in meters
%				time = input_time - 68280 (time of launch)
%

% constants
%WGS84=[6378.137 0.081819190842622]
ecc    = 0.081819190842622;   % first eccentricity of the earth
omegae = 7.2921159e-5 ;     % earth rotation rate in radians per second
% f      = 0.003352782 ;      % earth flattening factor
f = ecc2flat([6378.137 0.081819190842622]);
ae     = 6378137;  % length of earth semimajor axis in meters
ft2m = 12*.0254;


% loop to convert from ECI to ENU
for k=1:size(ECEF,1)
   
% rotate earth from GMT    
%	theta = omegae*ECI(k,1);
    
% coordinate transform (rotation matrix) from ECI to ECEF   
%	R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
   eci_pos = [ECEF(k,2) ECEF(k,3) ECEF(k,4)]';
      
   Pose = eci_pos;
   xe = Pose(1);        % x component in ECEF (kilometers)
   ye = Pose(2);			% y component in ECEF (meters)
   ze = Pose(3);			% z component in ECEF (meters)
   
% routines taken from WS-21231/1 Appendix 1
	rxy = sqrt(xe^2+ye^2);
	lat1 = atan(ze/((1-ecc^2)*rxy));
	if (ye<0)
    	long = -acos(xe/rxy);
	else
    	long = acos(xe/rxy);
	end
	r = sqrt(xe^2+ye^2+ze^2);
	phi = asin(ze/r);
	lambda = acos(rxy/r);
	gamma = f*(2-f);
	alpha = ae/sqrt(1+gamma*(sin(phi))^2/(1-gamma));
	del = gamma*sin(phi)*cos(lambda)/(1-gamma*(cos(lambda))^2);
	del1 = alpha*del/r;
	H = (r-alpha)*(1-del*del1/2);

	Rn1 = ae/sqrt(1-ecc^2*(sin(lat1))^2);
	lat2 = atan(ze*(Rn1+H)/((Rn1*(1-ecc^2)+H)*rxy));
	Rn2 = ae/sqrt(1-ecc^2*(sin(lat2))^2);
	lat = atan(ze*(Rn2+H)/((Rn2*(1-ecc^2)+H)*rxy));

	Long = long*180/pi;
   Lat = lat*180/pi;
   
% Form output matrix
	LLA(k,1) = ECEF(k,1);
	LLA(k,2) = Lat;
	LLA(k,3) = Long;
	LLA(k,4) = H;

end


