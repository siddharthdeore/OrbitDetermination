%% propogate and predict if satellite is visible or not
clc; clear all; close all;
addpath("../common","vallado");

%% Fetch satellite name
search_satellite = 'ENV';
[name,tle1,tle2]=fetchTLE(search_satellite,'visual.txt');
%number of orbits
ipArg=6;
min_in_day=1440;
%% Ground Station
stationLat=41.89*pi/180; % rad
stationLon=12.49*pi/180; % rad
stationAlt=0; % m
%% SGP4 configration
% input arguments
opsmode= 'a'; whichconst = 84; terms = 2; typerun='c'; typeinput = 'e';
min_in_day=1440;
timezone=0;
% convert TLE to satellite record
[startmfe, stopmfe, deltamin, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
period=2*pi/satrec.no; % period of satellite in min

%% Date time
current_utc_time=datetime('today');
%current_utc_time   = datetime('now','TimeZone','Europe/Rome'); % current date time
current_jd = juliandate(current_utc_time);          % current julian date
epoch_jd = satrec.jdsatepoch+satrec.jdsatepochf;
min_since_epoch = (current_jd-epoch_jd)*min_in_day; % one JD = 24*60 = min_in_day min;

%% Earth Orientation parameters
[xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(current_utc_time,'yyyy mm dd'));

switch ipArg
    case 1
        t= linspace(min_since_epoch,min_since_epoch+period*1,100); % get time vector current time to one orbit
    case 2
        t= linspace(min_since_epoch,min_since_epoch+3*period,300); % get time vector current time to 5 orbits
    otherwise
        t= min_since_epoch:1:min_since_epoch+1440; % get time vector current time to one solar day
end


%[year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf );
%fprintf("Time\nEpoch  %d %2d %d %d:%2d:%2.3f \n",day,mon,year,hr,min,sec);
fprintf('Epoch %s\n',datestr(datetime(satrec.jdsatepoch+satrec.jdsatepochf,'ConvertFrom','juliandate','TimeZone','Europe/Rome')));

% declaration of empty vectors
Lon=zeros(size(t));
Lat=Lon;

prev_index=0;
%% Propogator
for i=1:length(t)
    [satrec, r, v] = sgp4(satrec, t(i)); % r v is in TEME
    jd=satrec.jdsatepoch+satrec.jdsatepochf + t(i)/min_in_day;
    %[rsun,rtasc,decl] = sun(jd);
    [year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf + t(i)/min_in_day );
    %% convert time
    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
    jd=jdut1+ jdut1frac;
    %a=[0, 0, 0];
    %% ECI
    [reci, veci, aeci] = teme2eci  ( r', v', [0, 0, 0]', ttt, ddpsi, ddeps);
    %% Geodectic Lat Lon  -- Corrected
    [lst,gst] = lstime ( stationLon, jd);  % local and Greanwitch sideral time
    [lat,lon]=RVtoLatLong(rot3(r,gst)); % varified with orbitron
    %% Aximuth Elevation -- Corrected
    [rho,az,el,drho,daz,del] = rv2razel ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
    az=mod(az,2*pi);
    %% RA DEC -- Corected
    [rtasc,dec] = azel2radec(az,el,stationLat,stationLon,lst);
    
    %Lat(i)=lat;
    %Lon(i)=lon;
    Az(i)=az;
    El(i)=el;
    %RA(i)=rtasc;
    %DEC(i)=dec;
    %Rho(i)=rho;
    %dRho(i)=drho;
    
    %R(i,:)=reci; % V(i,:)=veci; Stored vector or % RA_DEC(i,:)=[rtasc,decl];
    %azel(i,:)=[az*180/pi,el*180/pi];
    
    strLight=light_1(reci,jd,'e');
    if(el > 5*pi/180 ) % chcek if elevation is  10 deg
        strLight=light_1(reci,jd,'e');
        if(strLight=="yes") % check if satellite is under daylight
            if((i-prev_index)>10) % check if there is gap between visible passes 
                fprintf('-------------------------------------------------------\n')
            end
            
            fprintf("%2d/%d/%d %02d:%02d:%02.2f \t ",day,mon,year,hr,min,sec);
            fprintf("RA:% 8.2f, Dec:% 8.2f %i\n",rtasc*180/pi,dec*180/pi,i);
            prev_index=i;            
            %fprintf("Az:%4.2f, El:%4.2f, RA:%f, Dec:%f, Rho%f dRho %f \n",az*180/pi,el*180/pi,rtasc*180/pi,dec*180/pi,rho, drho)
        end
    end
    
end


%{
if (El(1)>10*pi/180)
    strVisible='yes';
else
    strVisible='no';
end

fprintf(datestr(current_utc_time));
fprintf(strcat('\t\tIllum status : ',light_1( r,jd,'e'),'\t Visible : ',strVisible));
fprintf("\nLon : %2.2f Lat %2.2f | ",Lon(1)*180/pi,Lat(1)*180/pi)
fprintf("Az : %4.2f El: %4.2f | ",Az(1)*180/pi,El(1)*180/pi);
fprintf("RA: %4.2f DEC: %4.2f  | ",RA(1)*180/pi,DEC(1)*180/pi);
fprintf("Range: %4.2f RangeRate: %4.2f\n",Rho(1),dRho(1));

%% Sun Position
%gst = gstime(satrec.jdsatepoch+satrec.jdsatepochf + min_since_epoch/min_in_day); %Greanwitch sideral time
%[rsun,~,~] = sun(satrec.jdsatepoch+satrec.jdsatepochf + min_since_epoch/min_in_day);
%[sun_lat,sun_lon]=RVtoLatLong(rot3(rsun,gst));

%% plot track
%figure('Name','Ground Track');
%plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon);
%}

%% plot 3d trajectory ove globe with satellite
%{
figure('Name','Globe'); plotEarth(gst); hold on;
plot3(R(:,1),R(:,2),R(:,3),R(1,1),R(1,2),R(1,3),'rs');
axis equal;
%}
rmpath('../common');
rmpath('vallado');




%% Functions
%% Plot ground track from set of lat aand lon Vector
% plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon);
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
