%%
function [RA_comp,DEC_comp] = getComputedRADEC(state)
global terms min_in_day timezone rho_obs stationLat stationLon stationAlt xp yp dut1 lod ddpsi ddeps dat jd_obs t_obs
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
        az=mod(az,2*pi);
        
        %% RA DEC -- Corected for orbitron
       [rtasc,dec] = azel2radec(az,el,stationLat,stationLon,lst);
        %% RA DEC -- valadoo
        %[rr,rtasc,dec,drr,drtasc,ddecl] = rv2radec( reci,veci );
        %% RA DEC -- valadoo Topocentric
        %[~,rtasc,dec,~,~,~] = rv2tradc ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
        
        RA_comp(i,:)=rtasc;
        DEC_comp(i,:)=dec;
                        
    end
end
