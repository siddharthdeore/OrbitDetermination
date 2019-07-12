close; clear; clc;
addpath("celestrack");
ob=load('obsF.txt');
jd=ob(:,1);
RA=ob(:,2);
DEC=ob(:,3);

%{
%ob=load('obs.txt'); 
%ob=load('obs_delay.txt');
jdLocal=ob(:,1);
jd=ob(:,14);
ra1d=ob(:,2);
ra1m=ob(:,3);
ra1s=ob(:,4);
dec1d=ob(:,5);
dec1m=ob(:,6);
dec1s=ob(:,7);
ra2d=ob(:,8);
ra2m=ob(:,9);
ra2s=ob(:,10);
dec2d=ob(:,11);
dec2m=ob(:,12);
dec2s=ob(:,13);
for i=1:length(ra1d)
ra1(i,1)=dms2rad(ra1d(i),ra1m(i),ra1s(i));
ra2(i,1)=dms2rad(ra2d(i),ra2m(i),ra2s(i));
dec1(i,1)=dms2rad(dec1d(i),dec1m(i),dec1s(i));
dec2(i,1)=dms2rad(dec2d(i),dec2m(i),dec2s(i));
end
RA=(ra1+ra2)/2;
DEC=(dec1+dec2)/2;
%}

dd=[jd,RA,DEC];
%dd=[jdLocal,RA,DEC]

save observationsAll.txt dd -ascii -double
