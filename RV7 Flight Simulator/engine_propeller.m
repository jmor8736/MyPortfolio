function y = engine_propeller(u)

MAP = u(1);  % in Hg
TAS = u(2);  % kts
h   = u(3);  % ft
temp = u(4); % F

% y = [Power Thrust Eta RPM FuelConsumption BSFC];
% Power - hp
% Thrust - lbf
% Eta - 0 to 1
% FuelConsumption - US gal/ hour
% BSFC - lb/hp/hr
% MAP - inHg - minimum value - 13.5inHg, maximum value, 32inHg
% TAS - kts
% h - ft
% temp - F


%**************************************************************************
% Standard Atmosphere
%**************************************************************************
t_std = 288.15-6.5e-3*(h*0.3048);                                           % external temperature [K]
rho = ((1.048840 - 23.659414e-6*(h*0.3048))^4.2558797)/...                  % air density [kg/m3]
    (1+ ( ( (((temp-32)*5/9)+273.15) -t_std) /t_std) );

error = 10;
rpm = 2500;
k = 1;
while(error > 1 && k < 35)
    [power_eng,consumption] = engine_O360(rpm,MAP,h,temp,0);
    [Ct, Cp, thrust, power_prop, eta] = propeller(rpm,TAS,h,temp,0);
    error = power_eng-power_prop;
    rpm = rpm + 1*error;
%     rpm = min([rpm 3000]);
%     figure (1000)
%     plot(k,rpm,'o');
%     hold on
    k = k + 1;
    error = abs(error);
%     figure (100)
%     plot(power_eng,rpm,'.r'); hold on
%     plot(power_prop,rpm,'.b');
end


Power = power_eng;
Thrust = thrust;
Eta = eta/100;
RPM = rpm;
FuelConsumption = consumption;
BSFC = FuelConsumption*6.008/power_eng;

y = [Power Thrust Eta RPM FuelConsumption BSFC];
end