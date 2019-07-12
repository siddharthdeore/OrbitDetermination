clc;
clear all;
close all;

%% Tracking station Rome
lat  = 41.890251;
lon  = 12.492373;
alt  = 12;
year = 2019;
mon = 05;
day = 03;
hr  = 15;
min = 0;
sec = 0;
UTC  = [2019 05 03 14 20 36]; % you dont need this if you are propogating from satrec epoch
lla  = [lat,lon,alt];
p    = lla2eci(lla,UTC);    %% postion of loaction in J2000
% EOP = '2019 10 30 58786  0.172470  0.281573 -0.2224960  0.0010547 -0.113826 -0.008542  0.000090  0.000212  37';
%%
dat =  37;
dut1 = -0.2224960;
lod = 0.0010547;
xp = 0.172470;
yp = 0.281573;
terms = 0;
ddpsi = -0.113826;
ddeps = -0.008542;
timezone = 0;
%% TIme Coversiion

% # ----------------------------------------------------------------------------------------------------
% #   Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% # (0h UTC)           "         "          s          s          "        "          "         "     s 
% # ----------------------------------------------------------------------------------------------------
% 2019 11 01 58788  0.171384  0.280135 -0.2241006  0.0005471 -0.113339 -0.008425  0.000130  0.000236  37
% # ----------------------------------------------------------------------------------------------------

%% TLE of ENVSAT
longstr1 = '1 27386U 02009A   02060.19034371 -.00000044  00000-0  00000+0 0    26';
longstr2 = '2 27386  98.5327 128.8079 0013202 261.4580  98.5110 14.34878788    25';
% longstr1 = '1 27386U 02009A   02060.12060913 -.00000044  00000-0  00000+0 0';
% longstr2 = '2 27386  98.5327 128.7390 0013202 261.6647  98.2927 14.34878788';

[~, ~, ~, satrec] = twoline2rv(longstr1, longstr2,'c', [],[], 84);

time = satrec.jdsatepoch:1:satrec.jdsatepoch+100;

for i = 1:length(time)
     % take note when function returns a value you shoud store it in variable of same size
     % you can not allocate new memory
 [~, r, v] = sgp4(satrec,time(i));
 R(i,:)=r; % this line allocates new memory and saves returned variable from fucntiron
 V(i,:)=v;
 [rr,rtasc,decl,drr,drtasc,ddecl] = rv2radec(r,v);
 RA(i,:)=rtasc;
 DEC(i,:)=decl;
 [~, ~, jdut1,~, ~, ~, ~, ttt, ~, ~, ~, ~, ~, ~, ~, ~,~, ~, ~,~] = convtime ( year,mon, day, hr, min, sec , 0, dut1, dat);
 [reci, veci, ~] = teme2eci  ( r', v', 0, ttt, ddpsi, ddeps);
 [rho,az,el,~,~,~] = rv2razel ( reci,veci, lat,lon,alt,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);
 Rho(i,:)=rho;
 Az(i,:)=az;
 El(i,:)=el;
%% Checking if the satellite in in Light when it is above the ROME
if abs(el) > 15*pi/180 % you check elevation not declination
  [lit] = Light ( r, time(i), 'e'); %% change it to Light
else 
   display('Not above Rome');
end

end
%% Computed obervalbes this is not necessory since youa are already getting rho from rv2razel (ref: line 57)
rho = [cos(decl).*cos(rtasc);cos(decl).*sin(rtasc);sin(decl)];


%% Plotting 
% Spehere Plotting
figure
[xs,ys,zs]=sphere(16);
xs=xs*6371;
ys=ys*6371;
zs=zs*6371;
%cdata=imread('Earth.jpg');
%globe=surf(xs,ys,-zs,'FaceColor','none');
hold on;
set(globe,'FaceColor','texturemap','CData',cdata,'FaceAlpha',1,'EdgeColor','none');
plot3(R(:,1),R(:,2),R(:,3))
%plot3(Rho(1,:),Rho(2,:),Rho(3,:))
figure
 plot(DEC,RA)
    