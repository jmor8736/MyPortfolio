%  Standard & Off Standard Atmosphere Function [International Units]
%  for given values of pressure height

function [T, P, rho_std, rho] = atm_h(h, delta_T)
% inputs:
%   h = pressure height [meters]
%   delta_T = temperature offset [K]

% outputs:
%   T = temperature [K]
%   P = pressure [N/m^2]
%   rho_std = standard density [kg/m^3]
%   rho = density off standard [kg/m^3]

if h <= 11000
    A = 288.15;
    B = -6.5*10^(-3);
    C = 8.9619638;
    D = -.20216125*10^(-3);
    E = 5.2558797;
    I = 1.048840;
    J = -23.659414*10^(-6);
    L = 4.2558797;

    P = (C + D*h)^E;
    rho_std = (I + J*h)^L;

elseif h < 20000
    A = 216.65;
    B = 0;
    F = 128244.5;
    G = -.15768852*10^(-3);
    M = 2.0621400;
    N = -.15768852*10^(-3);

    P = F*exp(G*h);
    rho_std = M*exp(N*h);

elseif h <= 32000
    A = 196.65;
    B = 10^(-3);
    C = .34926867;
    D = 7.033098*10^(-6);
    E = -12.201149;
    I = .9726309;
    J = 4.94600*10^(-6);
    L = -35.163218;

    P = (C + D*h)^E;
    rho_std = (I + J*h)^L;

elseif h <= 47000
    A = 139.05;
    B = 2.8*10^(-3);
    C = .34936867;
    D = 7.0330980*10^(-6);
    E = -12.201149;
    I = .84392929;
    J = 16.993902*10^(-6);
    L = -13.201149;

    P = (C+D*h)^E;
    rho_std = (I + J*h)^L;

elseif h <= 50000
    A = 270.65;
    B = 0;
    F = 41828.42;
    G = -.12622656*10^(-3);
    M = .53839563;
    N = -.12622656*10^(-3);

    P = F*exp(G*h);
    rho_std = M*exp(N*h);
else
    disp('Error: Altitude out of range')
end

T = A + B*h;
rho = rho_std/(1 + (delta_T/T));

end
