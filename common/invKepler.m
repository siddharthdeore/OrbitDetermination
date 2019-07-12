function [E,v] = invKepler(M,e)
En  = M;
e=abs(e);
%% Circular
    E=M;
    v=M;
%% Parabolic
    %
if e == 1
    E = En - (En+En.^3/3- M)./(En.*En+1);
    while ( abs(E-En) > 1e-12 )
        En = E;
        E = En - (En+En.^3/3- M)./(En.*En+1);
    end
    v=2*atan(E);
end
%% Eliptic
if e>0 && e<1
    E = En - (En-e*sin(En)- M)/(1 - e*cos(En));
    while ( abs(E-En) > 1e-12 )
        En = E;
        E = En - (En - e*sin(En) - M)/(1 - e*cos(En));
    end
    % elliptic orbit
    sv = sqrt(1 - e * e) * sin(E);
    cv = cos(E) - e;
    v = atan2(sv, cv);
end
%% Hyperbolic
if e>1
    E = En - (e*sinh(En)-En- M)/(e*cosh(En)-1);
    while ( abs(E-En) > 1e-12 )
        En = E;
        E = En - (e*sinh(En)-En- M)/(e*cosh(En)-1);
    end
    % hyperbolic orbit
    sv = sqrt(e * e - 1) * sinh(E);
    cv = e - cosh(E);
    v = atan2(sv, cv);
end
end