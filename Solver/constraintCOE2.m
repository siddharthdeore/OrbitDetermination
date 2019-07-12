% constraint the orbital elements
function Xest = constraintCOE2(Xest)
if Xest(1) >=1
    Xest(1) = 0.99999; %ecc
elseif Xest(1) <=0
    Xest(1)=0.00001;
end
Xest(2)=rem(Xest(2),2*pi);  % incl
Xest(3)=rem(Xest(3),2*pi);  % argp
Xest(4)=rem(Xest(4),2*pi);  % omega
Xest(5)=rem(Xest(5),2*pi);  % omega
%if Xest(6) >=25
%    Xest(6) = 25;
%elseif Xest(6) <=0
%    Xest(6)=0.0000001;
%end
%if Xest(7) <0
%    Xest(7) = 0;
%end
end