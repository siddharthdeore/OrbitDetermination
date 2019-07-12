close all;
clear all;
clc;
%% TLE of ENVSAT
% Telescope is located on the ground in ROME.
% Station co-ordinates in terns of Lattidue , longitude, Altitude (These are in ECEF  
% we will conver it in to ECI)

lat  = 41.890251; % Degree
lon  = 12.492373; % Degree
alt  = 13;      % in M
lla  = [lat,lon,alt]; % Vector
%% EOP Paramters
xp   = 0.171384; 
yp   = 0.280135;
lod = 0.0005471;
terms = 0;
dut1 = -0.2233684; % UT1- UTC  from EOP
dat  = 37; % Atomic time - UT1  from EOP
ddpsi=-0.113590;
ddeps=-0.008416;
timezone= 0;
UTC  = [2019 05 05 17 20 36];
p    = lla2eci(lla,UTC);  % ECI Co-ordinates in Meters

%% ENVSAT TLE from celestrack 
longstr1 = '1 27386U 02009A   04366.94033539  .00000036  00000-0  29592-4 0   386';
longstr2 = '2 27386  98.5414  70.7168 0001032  98.0426 262.0872 14.32251628148423';

% Convert TLE into position amd velocity
[~, ~, ~, satrec] = twoline2rv(longstr1, longstr2,'c', [],[], 84);
time = satrec.jdsatepoch:1:satrec.jdsatepoch+100;

% Time conversion90
%[~, ~, jdut1, ~, ~, ~, ~, ttt, ~, ~, ~, ~,~,~,~,~,~,~,~,~] = convtime( year, mon, day, hr, min, sec,0, dut1, dat);
    

for i = 1:length(time);   
    [year,mon,day,hr,min,sec] = invjday (time(i),0);
    [~, r, v] = sgp4(satrec,time(i));   % Propagator    
    R(i,:)=r;
    V(i,:)=v;
    [~, ~, jdut1, ~, ~, ~, ~, ttt, ~, ~, ~, ~,~,~,~,~,~,~,~,~] = convtime( year, mon, day, hr, min, sec,0, dut1, dat);   
    [~,rtasc,decl,~,~,~] = rv2radec(r,v); % Converting rv to Ra, Dec in TEME
    
    [reci, veci, ~] = teme2eci( r', v', 0, ttt, ddpsi, ddeps); 
    
   %[rho,az,el,~,~,~] = rv2razel ( reci,veci, lat,lon,alt,ttt',jdut1',lod,xp,yp,terms,ddpsi,ddeps);
    % Visibility
    [lit] = light_1 ( r, time(i), 'e');
    
    % RHO calulations
    Rho(i,:) = [cos(rtasc)*cos(decl); 
                cos(decl)*sin(rtasc);
                sin(decl)];

end















%% Plotting 
 plot3(R(:,1),R(:,2),R(:,3))
