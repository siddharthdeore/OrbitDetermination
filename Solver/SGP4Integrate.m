clc; clear all; close all;
addpath("celestrack");
%%  TLE
%ISS (ZARYA)
tle1='1 25544U 98067A   19117.72128472  .00001222  00000-0  26993-4 0  9992';
tle2='2 25544  51.6413 252.0397 0001001 239.0409 324.7377 15.52596785167414';

%tle1='1 40029U PLANET   19101.50005787  .00000000  00000+0  19102-3 0    07'; tle2='2 40029 097.9238 024.3277 0012381 282.5665 309.0415 14.89736175    02';
%INSAT-4A                
%tle1='1 28911U 05049A   19117.56999302 -.00000176  00000-0  00000+0 0  9995';
%tle2='2 28911   0.0529 271.8348 0003047  69.4313 162.2070  1.00271614 48947';


%INTELSAT 5 (IS-5)       
%tle1='1 24916U 97046A   19101.08666139  .00000087  00000-0  00000+0 0  9993'; tle2='2 24916   4.8431  67.8322 0004355 320.5247  64.8378  1.00272374 79303';
%% Ground Station
stationLat=41.890251*pi/180; % rad
stationLon=12.483373*pi/180; % rad
stationAlt=12; % m
%% SGP4 configration
% input arguments
opsmode= 'a'; whichconst = 84; eqeterms = 1; typerun='c'; typeinput = 'e';
% convert TLE to satellite record
[startmfe, stopmfe, deltamin, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
period=2*pi/satrec.no; % period of satellite in min
%% Some Constants
lod=0;xp=0;yp=0;
dut1 = -0.4399619; dat  = 32; ddpsi = -0.052195 * pi / (180*3600); ddeps = -0.003875 * pi / (180*3600);


%% Date time
current_utc_time   = datetime('now','TimeZone','UTC'); % current date time
current_local_time = datetime('now','TimeZone','Europe/London'); % current date time
current_jd = juliandate(current_utc_time); % current julian date
epoch_jd = satrec.jdsatepoch+satrec.jdsatepochf; 

min_since_epoch = (current_jd-epoch_jd)*1440; % one JD = 24*60 = 1440 min;
%for i=0:min_since_epoch
%[satrec, ~, ~] = sgp4(satrec, i);    
%end 
t= linspace(min_since_epoch,min_since_epoch+period*1,200); % get time vector current time to one orbit 

[year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf );
fprintf("Time\nEpoch  %d %2d %d  %d:%2d:%2.3f \n",day,mon,year,hr,min,sec);
fprintf("UTC    %s\n",datestr(current_utc_time,  'dd mm yyyy HH:MM:SS.FFF'));
fprintf("Local  %s\n",datestr(current_local_time,'dd mm yyyy HH:MM:SS.FFF'));

% declaration of empty vectors
Lon=zeros(size(t));
Lat=Lon;


%% Propogator
for i=1:length(t);
[satrec, r, v] = sgp4(satrec, t(i)); % r v is in TEME
jd=satrec.jdsatepoch+satrec.jdsatepochf + t(i)/1440;
[rsun,rtasc,decl] = sun(jd);
[year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf + t(i)/1440 );
% convert time
[ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] = convtime ( year, mon, day, hr, min, sec, 0, dut1, dat );
% ECI
[reci, veci, aeci] = teme2eci  ( r', v', [1,1,1]', ttt, ddpsi, ddeps);
% Aximuth Elevation
[~,az,el,~,~,~] = rv2razel ( reci,veci, stationLat,stationLon,stationAlt,ttt,jdut1,0,0,0,2,ddpsi,ddeps );
gst = gstime(jdut1+jdut1frac); %Greanwitch sideral time
U=R_z(gst); % U is rotation matrix rotated by angle gst
[lat,lon]=RVtoLatLong(rot3(reci,gst));
if i==1
    [rsun,rtasc,decl] = sun(satrec.jdsatepoch+satrec.jdsatepochf + t(1)/(24));
    [sun_lat,sun_lon]=RVtoLatLong(rot3(rsun,gst));
end
Lat(i)=lat;
Lon(i)=lon;
if Lon(i)<=-pi
 Lon(i)=Lon(i)+2*pi;
elseif Lon(i)>=pi
 Lon(i)=Lon(i)-2*pi;
end
% RV(i,:)=reci; V(i,:)=veci; Stored vector or %t_R_V_Lat_lon(i,:) =[t(i),reci',veci',Lat(i),Lon(i)];
%pause(0.001);
end
%figure
%plot3(R(:,1),R(:,2),R(:,3));
plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon);
fprintf("\n Lon %2.2f Lat %2.2f",Lon(1)*180/pi,Lat(1)*180/pi);

rmpath('celestrack');

%% Returns rotation matirx of Rotate around z axis by angle
function rotmat = R_z(angle)
  C = cos(angle);
  S = sin(angle);
  rotmat = zeros(3,3);
  rotmat(1,1) =    C; rotmat(1,2) = S; rotmat(1,3) = 0;
  rotmat(2,1) = -1*S; rotmat(2,2) = C; rotmat(2,3) = 0;
  rotmat(3,1) =    0; rotmat(3,2) = 0; rotmat(3,3) = 1;
end

%% Plot ground track from set of lat aand lon Vector
function plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon)
img = imread('pixmaps/earthmap.jpg'); % import imgage
image('CData',img,'XData',[-180 180],'YData',[90 -90]); % scale image xy axis
hold on;
plot(Lon*180/pi,Lat*180/pi,'w.','LineWidth',1); %plots
axis([-180  180  -90 90]);
plot(stationLon(1)*180/pi,stationLat(1)*180/pi,'r.','LineWidth',1); % Station cordinates
plot(sun_lon*180/pi,sun_lat*180/pi,'yO','LineWidth',8); % Station cordinates
plot(Lon(1)*180/pi,Lat(1)*180/pi,'ks','LineWidth',4); axis([-180  180  -90 90]); xlabel('Longitude');ylabel('Lattitude');
end


function [Lat, Long] = RVtoLatLong(ECEF)
%orbit radius velocity latitude longitude ECEF
%
% [Lat, Long] = RVtoLatLong(ECEF)
% Purpose:  This fuction convertes ECEF coordinates to Geocentric latitude
%           and longitude given ECEF radius in km. Valid for any planetary
%           body.
%
% Inputs:  o ECI  - A 3x1 vector of Earth-Centered Inertial (IJK)
%                   coordinates in km.
%
% Outputs: o Lat  - Geocentic latitude of spacecraft in radians
%          o Long - Longitude of spacecraft in radians

r_delta = norm(ECEF(1:2));

sinA = ECEF(2)/r_delta;
cosA = ECEF(1)/r_delta;

Long = atan2(sinA,cosA);

if Long < -pi
    Long = Long + 2*pi;
end

Lat = asin(ECEF(3)/norm(ECEF));
end
