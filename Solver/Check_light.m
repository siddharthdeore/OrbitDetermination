clc; clear all; close all;
addpath('celestrack');
%Two line element set from celestrack.com
longstr1='1 03230U 68040B   19079.51173572  .00000107  00000-0  14911-4 0  9991';
longstr2='2 03230  74.0367  39.3621 0030513 183.2419 176.8570 14.90451929731055';

% input arguments
opsmode= 'a';  % afspc approach
whichconst = 84;
eqeterms = 1;
typerun='c';
typeinput = 'e';

%convert TLE to satellite record
[startmfe, stopmfe, deltamin, satrec] = twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst);

period=2*pi/satrec.no; %period of satellite
% from jullian day to readable date time
[year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf );

fprintf("Epoch UTC Time %d - %d - %d %d:%d:%2.2f \n",year,mon,day,hr,min,sec);
utc =[year,mon,day,hr,min,sec];
figure();
%% propogate satrec from epoch to time you want
to_time = period; %[min]
k=0;
for i=0:1:to_time;
[satrec, r, v] = sgp4(satrec, i);
r_save(i+1,:)=r;
v_save(i+1,:)=v;
title(strcat('epoch time utc + ', num2str(i) , ' min'));
[rr,ecllon,ecllat,drr,decllon,decllat] = rv2ell (r,v);
elonlat_save(i+1,:)=[rr,ecllon,ecllat,drr,decllon,decllat];
re=r/6371;
[year,mon,day,hr,min,sec]=invjday(satrec.jdsatepoch,satrec.jdsatepochf);

[jd,jf]=jday(year,mon,day,hr,min+i,sec);

    res=light_1(r/6371,jd+jf,'s');
    if(res=='yes');
        k=k+1; % minute array increments
        d(k)=jd+jf;
    end

%pause(0.001);

end
subplot(2,1,1);
plot3(r_save(:,1),r_save(:,2),r_save(:,3),r(1),r(2),r(3),'r*');
%axis([-8000 8000 -8000 8000 -8000 8000]);
axis equal
subplot(2,1,2);
axis([-180 180 -90 90]);
plot((elonlat_save(:,2)*180/pi),(elonlat_save(:,3)*180/pi),'k.');

for i=1:length(d)
        [year,mon,day,hr,min,sec] = invjday ( d(i), jf );
        fprintf("visble at UTC Time %d - %d - %d %d:%d:%2.2f \n",year,mon,day,hr,min,sec);
end


figure('Name','Position')
plot(r_save(:,1))
hold on;
plot(r_save(:,2))
plot(r_save(:,3))
legend('x','y','z');
xlabel('time from epoch [min]');
ylabel('position [km]');

figure('Name','Velocity')
plot(v_save(:,1));
hold on;
plot(v_save(:,2));
plot(v_save(:,3));
legend('vx','vy','vz');
xlabel('time from epoch [min]');
ylabel('velocity [km/s]');

rmpath('celestrack');
