%%
function [cost] = computeCost(state)
global terms min_in_day timezone rho_obs stationLat stationLon stationAlt jd_obs t_obs
global  xp  yp  dut1  lod  ddpsi  ddeps  dx  dy dat
global xp1 yp1 dut11 lod1 ddpsi1 ddeps1 dx1 dy1 dat1
global xp2 yp2 dut12 lod2 ddpsi2 ddeps2 dx2 dy2 dat2

    %[~, ~, ~, satrec] = twoline2rvDouble1(state); % with bstar
    [~, ~, ~, satrec] = twoline2rvDouble2(state); % without bstar
    
    epoch_jd=satrec.jdsatepoch + satrec.jdsatepochf;
    min_since_epoch = (jd_obs(1)-epoch_jd)*min_in_day;
%    [year,mon,day,hr,min,sec] = invjday ( satrec.jdsatepoch, satrec.jdsatepochf );
    
    t_mins=(jd_obs-jd_obs(1))*min_in_day;
    %% Propogator
    NN=length(t_mins);
    rho_comp=zeros(size(rho_obs));
    RA_comp = zeros(size(t_obs));
    DEC_comp = RA_comp;
    theta = RA_comp;
    for i=1:NN
        [satrec, r, v] = sgp4(satrec, min_since_epoch+t_mins(i));	% r v is in TEME
        [year,mon,day,hr,min,sec] = invjday ( epoch_jd,(min_since_epoch+t_mins(i))/min_in_day);
        if day==22
            xp=xp1; yp=yp1; dut1=dut11; lod=lod1; ddpsi=ddpsi1; ddeps=ddeps1; dx=dx1; dy=dy1; dat=dat1;
        elseif day==23
            xp=xp2; yp=yp2; dut1=dut12; lod=lod2; ddpsi=ddpsi2; ddeps=ddeps2; dx=dx2; dy=dy2; dat=dat2;
        end

        %% convert time
        [~, ~, jdut1, jdut1frac, ~, ~, ~, ttt, ~,~, ~, ~, ~,~, ~, ~,~, ~, ~,~ ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        jd=jdut1+ jdut1frac;
        a=[0, 0, 0];
        %% ECI
        [reci, veci, ~] = teme2eci  ( r', v', a', ttt, ddpsi, ddeps);
        %% Geodectic Lat Lon  -- Corrected
        [lst,~] = lstime ( stationLon, jd);  % local and Greanwitch sideral time
        %% Aximuth Elevation -- Corrected
        [~,az,el,~,~,~] = rv2razel ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
        z=mod(az,2*pi);
        %% RA DEC -- Corected for orbitron
        [rtasc,dec] = azel2radec(az,el,stationLat,stationLon,lst);
        %% RA DEC -- valadoo
        %[rr,rtasc,decl,drr,drtasc,ddecl] = rv2radec( reci,veci );
        
        RA_comp(i,:)=rtasc;
        DEC_comp(i,:)=dec;
                        
        rho_comp(i,:)=getRhoVector(RA_comp(i,:),DEC_comp(i,:));
        theta(i) = angl(rho_comp(i,:),rho_obs(i,:));
    end
    cost = 0.5/NN*(theta*theta'); % least square
end
%% some own functions
function [rho] = getRhoVector(rtasc,dec)
% retruns a vector from RA and DEC
rho=[cos(rtasc).*cos(dec), cos(dec).*sin(rtasc), sin(dec)];
end
