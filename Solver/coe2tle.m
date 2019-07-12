%%
function [tle1,tle2] = coe2tle(tle1,tle2,ecc,inc,rtasc,argp,m,no)
% Update tle string using kepler orbital elemets
% input
% tle1 - tle1 string
% tle2 - tle2 string
% inc - inclination [rad]
% rtasc - Right ascension of asccending node [rad]
% ecc - eccentricity  (0 < ecc < 1)
% argp - argument of perigee [rad]
% m - mean annomaly [rad]
% no - mean motion rad per min
%% siddharth@deore.in
inc=mod(inc*180/pi,360);
rtasc=mod(rtasc*180/pi,360);
argp=mod(argp*180/pi,360);
m=mod(m*180/pi,360);
%no=n;%no*1440/(2*pi); % rad per min to rev per day
str_Inc=pad(num2str(inc),8,'left');
str_RA=pad(num2str(rtasc),8,'left');
str_Ecc=num2str(fix((ecc-fix(ecc))*1e7),'%07d');	% Eccentricity (only fraction part)
str_argp=pad(num2str(argp),8,'left');               % argument of perigee [deg]
str_m=pad(num2str(m,'%8.4f'),8,'left');             % mean anomaly [deg]
str_no=pad(num2str(no,'%11.8f'),11,'left');              % mean motion rev per day
tle2(9:16)=str_Inc;
tle2(18:25)=str_RA;
tle2(27:33)=str_Ecc;
tle2(35:42)=str_argp;
tle2(44:51)=str_m;
tle2(53:63)=str_no;
%fprintf('%s\n',tle2);
end
