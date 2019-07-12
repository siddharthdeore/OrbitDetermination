function [T] = orbital_period(a,GM)
    T = 2*pi*(a.*a.*a./GM).^0.5;
end