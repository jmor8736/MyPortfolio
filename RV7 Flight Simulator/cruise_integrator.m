% Jackson Morgan
% Cruise Function

function [Power, state2, vdot, SR] = cruise_integrator(state1, t_step, delta_T)
% INPUTS:
%   t_step - time step [s]
%   delta_T - temp offset [K]
%   cruise_type: 1-constant altitude, 2-???
%   state1 = [t x h w v MAP]
%       t - [s]
%       x - [nm]
%       h - [ft]
%       w - [lbs]
%       v - [kts]
%       MAP - [in hg]

% OUTPUTS:
%   vdot - acceleration [nm/s^2]
%   power = [Cl Cd Drag Thrust Pava Preq]
%       Cl - []
%       Cd - []
%       Drag - [lbf]
%       Thrust - [lbf]
%       Pava - [hp]
%       Preq - [hp]
%   state2 = [t x h w v MAP]
%       t - [s]
%       x - [nm]
%       h - [ft]
%       w - [lbs]
%       v - [kts]
%       MAP - [in hg]


% unit conversion factors:
ft2m = 0.305; % feet to meters
pa2hg = 0.00029529980164712; % pascals to lbs of mercury
kgm2slgft = .0019; % kg/m^3 to slug/ft^3
kt2ft = 6076.12; % nautical miles to feet
kt2fps = 1.688; % knots to feet per second
lb2slg = 0.031081; % lb to slug

% Extracting Values from State Vector:
time = state1(1); % [s]
x = state1(2)*kt2ft; % [ft]
h = state1(3); % [ft]
w = state1(4); % [lbs]
v = state1(5); % [kts]
MAP = state1(6); % [in hg]

% Atmospheric Conditions:
[T, P, ~, rho] = atm_h(h*ft2m, delta_T);

% RV-7 Specs:
S = 121; % [sq ft]

% Convert to Freedom Units
rho = rho/(16.018*32.2);% [slug/ft^3]
temp = (T - 273.15) * 9/5 + 32; % [deg F]

% Check that MAP is less than atmospheric pressure:
if MAP > P*pa2hg
    MAP = P*pa2hg; % Set MAP equal to Patm if it was higher
    disp(['MAP Adjusted to ', num2str(MAP),' in hg'])
end


    % Getting Flight Deck Info:
    u = [MAP,v,h,temp]; % Iscold's function setup
    y = engine_propeller(u); % call on Iscold function
    THRUST = y(2); % [lbf]
    Power_available = THRUST*(v*kt2fps); % [lb*ft/s]

    CL = 2*w/(rho*S*(v*kt2fps)^2); 
    CD = 0.024292 + 1.667*10^-9*CL + 0.071696*CL^2;
    Drag = 0.5*CD*rho*((v*kt2fps)^2)*S; % [lbf]
    Power_Required = Drag*(v*kt2fps);
    
    xdot = v;
    wdot = (-y(6)*y(1))/3600; % [lb/s]
    vdot = (THRUST - Drag)/(w*lb2slg); % [ft/s^2]
    hdot = 0;


% Update State:
x2 = x + xdot*t_step; % [ft]
x2 = x2/kt2ft; % [nm]
w2 = w + wdot*t_step;
v2 = v + vdot*t_step;
h2 = h + hdot*t_step;
MAP2 = MAP;
time2 = time + t_step;

SR = (xdot*kt2fps)/wdot; % specific range [ft/lb]
SR = -SR/kt2ft; % specific range [nm/lb]
    
% Output Vectors:
state2 = [time2; x2; h2; w2; v2; MAP2];
Power = [Power_available/550; Power_Required/550; THRUST; Drag; CL; CD ];

end