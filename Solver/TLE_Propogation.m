clc; clear all; close all;
addpath("../common");
addpath("vallado");
%% Fetch satellite TLE from Name
search_satellite = 'ENV';
[name,tle1,tle2]=fetchTLE(search_satellite,'visual.txt');

%%  TLE
%ISS (ZARYA)
%tle1='1 25544U 98067A   19117.72128472  .00001222  00000-0  26993-4 0  9992';
%tle2='2 25544  51.6413 252.0397 0001001 239.0409 324.7377 15.52596785167414';

%tle1='1 40029U PLANET   19101.50005787  .00000000  00000+0  19102-3 0    07'; tle2='2 40029 097.9238 024.3277 0012381 282.5665 309.0415 14.89736175    02';
%INSAT-4A                
%tle1='1 28911U 05049A   19117.56999302 -.00000176  00000-0  00000+0 0  9995';
%tle2='2 28911   0.0529 271.8348 0003047  69.4313 162.2070  1.00271614 48947';


%INTELSAT 5 (IS-5)       
%tle1='1 24916U 97046A   19101.08666139  .00000087  00000-0  00000+0 0  9993'; tle2='2 24916   4.8431  67.8322 0004355 320.5247  64.8378  1.00272374 79303';
%% Ground Station
stationLat=41.890251*pi/180; % rad
stationLon=12.483373*pi/180; % rad
stationAlt=0; % m
timezone=0;
%% Earth Orientation parameters
current_utc_time   = datetime('now','TimeZone','Europe/Rome'); % current date time
[xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(current_utc_time,'yyyy mm dd'));

%% SGP4 configration input arguments
opsmode= 'a'; whichconst = 84; terms = 2; typerun='c'; typeinput = 'e'; min_in_day=1440;
% convert TLE to satellite record
[startmfe, stopmfe, deltamin, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
period=2*pi/satrec.no; % period of satellite in min

%% Date time
current_utc_time   = datetime('now','TimeZone','Europe/London'); % current date time
current_jd = juliandate(current_utc_time); % current julian date
epoch_jd = satrec.jdsatepoch+satrec.jdsatepochf; 

min_since_epoch = (current_jd-epoch_jd)*min_in_day; % one JD = 24*60 = 1440 min;

t= linspace(min_since_epoch,min_since_epoch+period*1,200); % get time vector current time to one orbit 

[year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf );
fprintf("Time\nEpoch  %d %2d %d  %d:%2d:%2.3f \n",day,mon,year,hr,min,sec);
fprintf("Local    %s\n",datestr(current_utc_time,  'dd mm yyyy HH:MM:SS.FFF'));

% declaration of empty vectors
Lon=zeros(size(t));
Lat=Lon;


%% Propogator
for i=1:length(t);
[satrec, r, v] = sgp4(satrec, t(i)); % r v is in TEME
jd=satrec.jdsatepoch+satrec.jdsatepochf + t(i)/min_in_day;
[rsun,rtasc,decl] = sun(jd);
[year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf + t(i)/min_in_day );
%% convert time
[ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
jd=jdut1+ jdut1frac;
a=[0, 0, 0];
%% ECEF
[recef, vecef, aecef] = teme2ecef( r',v',a',ttt,jdut1,lod,xp, yp,terms);
%% ECI
[reci, veci, aeci] = teme2eci  ( r', v', a', ttt, ddpsi, ddeps);
%% Geodectic Lat Lon  -- Corrected
[lst,gst] = lstime ( stationLon, jd);

[lat,lon]=RVtoLatLong(rot3(r,gst)); % varified with orbitron 
%[latgc,latgd,loni,hellp] = ijk2lle ( r',jd); % direct lat lon without gst earth rotation
%% Aximuth Elevation -- Corrected
[rho,az,el,drho,daz,del] = rv2razel ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
%% RA DEC -- Corected
%tropocentric
%[rho,trtasc,tdecl,drho,dtrtasc,dtdecl] = rv2tradc ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
%geocentric
%[~,rtasc,decl,~,~,~] = rv2radec( reci,veci);
dec = asin(sin(el)*sin(stationLat)+cos(el)*cos(stationLat)*cos(az));
local_hour_angle = atan2(-sin(az)*cos(el)/cos(dec),(sin(el)-sin(dec)*sin(stationLat))/(cos(dec)*cos(stationLat)));
rtasc=lst-local_hour_angle;

Lat(i)=lat;
Lon(i)=lon;
Az(i)=az;
El(i)=el;
RA(i)=rtasc;
DEC(i)=dec;
Rho(i)=rho;
dRho(i)=drho;

if i==1
    gst = gstime(satrec.jdsatepoch+satrec.jdsatepochf + min_since_epoch/min_in_day); %Greanwitch sideral time
    [rsun,~,~] = sun(satrec.jdsatepoch+satrec.jdsatepochf + min_since_epoch/min_in_day);
    [sun_lat,sun_lon]=RVtoLatLong(rot3(rsun,gst));
end
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

rmpath('../common');
rmpath('vallado');

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
function plotGruondTrack_mapbox(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon)
worldmap('World')
load coastlines;
plotm(coastlat,coastlon);

hold on;
geoshow(Lat*180/pi,Lon*180/pi);

geoshow(stationLat(1)*180/pi,stationLon(1)*180/pi,'DisplayType','point'); % Station cordinates
geoshow(sun_lat*180/pi,sun_lon*180/pi,'DisplayType','point'); % Sun cordinates
geoshow(Lat(1)*180/pi,Lon(1)*180/pi,'DisplayType','point');
%axis([-180  180  -90 90]);
end
%% Plot ground track from set of lat aand lon Vector
function plotGruondTrack_old(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon)
img = imread('pixmaps/earthmap2k.png'); % import imgage
image('CData',img,'XData',[-180 180],'YData',[90 -90]); % scale image xy axis
hold on;
plot(Lon*180/pi,Lat*180/pi,'w-','LineWidth',1); %plots
axis([-180  180  -90 90]);
plot(stationLon(1)*180/pi,stationLat(1)*180/pi,'r.','LineWidth',1); % Station cordinates
plot(sun_lon*180/pi,sun_lat*180/pi,'yO','LineWidth',8); % Station cordinates
plot(Lon(1)*180/pi,Lat(1)*180/pi,'ks','LineWidth',4); axis([-180  180  -90 90]); xlabel('Longitude');ylabel('Lattitude');
end

%% Plot ground track from set of lat aand lon Vector
function plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon)
img = imread('pixmaps/earthmap2k.png'); % import imgage
image('CData',img,'XData',[-180 180],'YData',[90 -90]); % scale image xy axis
hold on;
tol = 5; % distance > tol indicates discontinuity
dl = diff([Lat;Lon],1,2); % look up what that command does if you don't know it
euler_dist = sqrt((dl(1,:)+dl(2,:)).^2); % distance between data points

jumpind = [0 euler_dist>tol]; % now if jumpind(i) = true, we know that the 
                 %   point [lat(i) lon(i)] is the first after a jump
blocks = cumsum(jumpind); % points that belong to the same continuous part
%disp(euler_dist)
% Now just loop over the continuous blocks to draw a separate line for each one
for i=0:blocks(end)
    plot(Lon(blocks==i)*180/pi,Lat(blocks==i)*180/pi,'w-','LineWidth',1);
    hold on;
end
axis([-180  180  -90 90]);
plot(stationLon(1)*180/pi,stationLat(1)*180/pi,'r.','LineWidth',1); % Station cordinates
plot(sun_lon*180/pi,sun_lat*180/pi,'yO','LineWidth',8); % Station cordinates
plot(Lon(1)*180/pi,Lat(1)*180/pi,'ks','LineWidth',4);
xlabel('Longitude');ylabel('Lattitude');
end
