%% realtime propogation for TLE
clc; clear all; close all;
addpath("../common","vallado");

%% Fetch satellite TLE from Name
search_satellite = 'ENV';
[name,tle1,tle2]=fetchTLE(search_satellite,'visual.txt');
%% Ground Station
stationLat=41.89*pi/180; % rad
stationLon=12.49*pi/180; % rad
stationAlt=0; % m
timezone=0;
%% Earth Orientation parameters
current_utc_time   = datetime('now','TimeZone','Europe/Rome'); % current date time
[xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(current_utc_time,'yyyy mm dd'));
prev_utc_time=current_utc_time; % for loop timer

%% SGP4 configration input arguments
opsmode= 'a'; whichconst = 84; terms = 2; typerun='c'; typeinput = 'e'; min_in_day=1440;
%% convert TLE to satellite record
[startmfe, stopmfe, deltamin, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
% period of satellite in min
period=2*pi/satrec.no;

figureHandle=figure('Name','Ground Track');
while(ishghandle(figureHandle))    
    current_utc_time   = datetime('now','TimeZone','Europe/Rome'); % current date time
    current_utc_time = dateshift(current_utc_time, 'start', 'second'); % round to nearest seconds
    if(seconds(current_utc_time-prev_utc_time)>=1) % loop only after 1 sec
        %% Date time
        current_jd = juliandate(current_utc_time); % current julian date
        epoch_jd = satrec.jdsatepoch+satrec.jdsatepochf;
        
        min_since_epoch = (current_jd-epoch_jd)*min_in_day; % one JD = 24*60 = 1440 min; this is min elapsed since epoch
        t= linspace(min_since_epoch,min_since_epoch+period,50); % get time vector current time to one orbit
        
        %t=min_since_epoch; % only one time
        
        % declaration of empty vectors
        Lon=zeros(size(t)); Lat=Lon;Az=Lon;El=Lon;RA=Lon;DEC=Lon; R=zeros(length(t),3); Rho=Lon;dRho=Lon;
        %% Propogate to find ground track
        for i=1:length(t)
            [satrec, r, v] = sgp4(satrec, t(i)); % r v is in TEME
            [year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch , satrec.jdsatepochf  + t(i)/min_in_day);
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
            %gst = gstime(jd);      % local and Greanwitch sideral time
            %lst=gst+stationLon; % local sidereal time angle
            [lat,lon]=RVtoLatLong(rot3(r,gst)); % varified with orbitron
            %[latgc,lat,lon,hellp] = ijk2lle ( r',jd); % direct lat lon without gst earth rotation
            %% Aximuth Elevation -- Corrected
            [rho,az,el,drho,daz,del] = rv2razel ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
            az=mod(az,2*pi);
            %% RA DEC -- Corected
            [rtasc,dec] = azel2radec(az,el,stationLat,stationLon,lst); % my correct
            %[rtasc,dec] = azl2radc(az, el, stationLat, lst); % valladoo correct
            % store the variables
            Lat(i)=lat;
            Lon(i)=lon;
            Az(i)=az;
            El(i)=el;
            RA(i)=mod(rtasc,2*pi);
            DEC(i)=dec;
            Rho(i)=rho;
            dRho(i)=drho;
            
            R(i,:)=r;
        end
        %% Check satellite visibility over observer horizon
        if el>10*pi/180
            strVisible='yes';
        else
            strVisible='no';
        end
        
        fprintf(datestr(current_utc_time,'mm-dd-yyyy HH:MM:SS'));
        %fprintf(strcat('\t\tIllum : ',light_1( r,jd,'e'),'\t Visible : ',strVisible));
        %fprintf("\nLon : %4.4f Lat %4.4f \t",Lon(1)*180/pi,Lat(1)*180/pi)
        %fprintf("Az : %4.1f El: %4.1f \t",Az(1)*180/pi,El(1)*180/pi);
        fprintf("\t RA: %4.4f DEC: %4.4f \n",RA(1)*180/pi,DEC(1)*180/pi);
        
        
        %% Sun Position
        gst = gstime(satrec.jdsatepoch+satrec.jdsatepochf + t(i)/min_in_day); %Greanwitch sideral time
        [rsun,~,~] = sun(satrec.jdsatepoch+satrec.jdsatepochf + t(i)/min_in_day);
        [sun_lat,sun_lon]=RVtoLatLong(rot3(rsun,gst));
        
        %% plot track
        plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon);
        drawnow;
        
        prev_utc_time=current_utc_time;
    end
end
%% plot 3d trajectory ove globe with satellite
%figure('Name','Globe'); plotEarth(gst); hold on;
%plot3(R(:,1),R(:,2),R(:,3),R(1,1),R(1,2),R(1,3),'rs');
%axis equal;

rmpath('../common');
rmpath('vallado');


%% Functions
%% Plot ground track from set of lat aand lon Vector
% plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon);
function plotGruondTrack(Lat,Lon,stationLat,stationLon,sun_lat,sun_lon)
img = imread('pixmaps/earthmap2k.png'); % import imgage
image('CData',img,'XData',[-180 180],'YData',[90 -90]); % scale image xy axis
hold on;
tol = 3; % distance > tol indicates discontinuity
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
grid minor;
xlabel('Longitude');ylabel('Lattitude');

end
