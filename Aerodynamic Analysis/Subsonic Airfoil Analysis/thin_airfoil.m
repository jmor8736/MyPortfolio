% Jackson Morgan
% Thin Airfoil Theory Function

function [cl, cmc4] = thin_airfoil(NACA,alpha)
% Outputs:
%   cl = lift coefficient
%   cd = drag coefficient
% Inputs:
%   NACA = 4-digit NACA designation
%   alpha = angle of attack in degrees
%   n = number of differential sections for calculations

% obtain m and p from NACA number
NACA = num2str(NACA);
m = str2double(NACA(1))/100;
p = str2double(NACA(2))/10;

% convert alpha to radians
alph = alpha*(pi/180);

% calculate A0
lim = acos(1-2*p);
dzdx1 = @(theta) (m/p^2)*(2*p-(1-cos(theta)));
dzdx2 = @(theta) (m/(1-p)^2)*(2*p-(1-cos(theta)));
A0 = alph - (1/pi)*(integral(dzdx1,0,lim) + ...
    integral(dzdx2,lim,pi));

% calculate A1
integrand11 = @(theta) (m/p^2)*(2*p-(1-cos(theta))).*cos(theta);
integrand12 = @(theta) (m/(1-p)^2)*(2*p-(1-cos(theta))).*cos(theta);
A1 = (2/pi)*(integral(integrand11,0,lim) + integral(integrand12,lim,pi));

% Calculate A2
integrand21 = @(theta) (m/p^2)*(2*p-(1-cos(theta))).*cos(2*theta);
integrand22 = @(theta) (m/(1-p)^2)*(2*p-(1-cos(theta))).*cos(2*theta);
A2 = (2/pi)*(integral(integrand21,0,lim) + integral(integrand22,lim,pi));

% use A0 values to find cl and cm
cl = 2*pi*(A0 + A1/2);
cmc4 = (pi/4)*(A2-A1);

end



