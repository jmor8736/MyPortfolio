% Jackson Morgan
% AERO 425
% Climb Euler Integrator

function [state2,ROC,SR] = climb_integrator(state1, t_step, delta_T)
% INPUTS:
%   t_step - time interval [s]
%   delta_T - temp offset [K]
%   state = [t x h w v MAP]
%       t - seconds
%       x - nautical miles
%       h - feet
%       w - lbs
%       v - knots
%       MAP - in hg

% OUTPUTS:
%   ROC - vertical speed [ft/s]
%   state2 = [t x h w v MAP]
%       t - [s]
%       x - [nm]
%       h - [ft]
%       w - [lbs]
%       v - [kts]
%       MAP - [in hg]
%   SR - specific range [nm/lb]

% unit conversion factors:
ft2m = 0.305; % feet to meters
pa2hg = 0.00029529980164712; % pascals to lbs of mercury
kgm2slgft = .0019; % kg/m^3 to slug/ft^3
kt2ft = 6076.12; % nautical miles to feet
kt2fps = 1.688; % knots to feet per second
lb2slg = 0.031081; % lb to slug

% RV-7 specs:
S = 121; % wing surface area [sqft]
Cdo = .024292; % parasitic drag coefficient
k = .071696;

% Extract State Vector Values:
t = state1(1); % time [s]
x = state1(2); % initial distance [nautical miles]
h = state1(3); % initial altitude [ft]
w = state1(4); % initial weight [lbs]
v = state1(5); % TAS [kts]
MAP = state1(6); % initial Manifold pressure

% check MAP and adjust accordingly:
[Temp, Press, ~, rho] = atm_h(h*ft2m, delta_T); % Temp [K], pressure [Pa]
P_hg = Press*pa2hg; % pressure [in hg]
if MAP > P_hg
    MAP = P_hg;
    disp(['MAP adjusted to ',num2str(MAP),' in hg'])
end

% atmospheric conditions:
Temp = (Temp-273.15)*(9/5) + 32; % temp [F]
rho = rho*kgm2slgft; % density [slugs/ft^3]

% ISCOLD FUNCTION:
u = [MAP v h Temp]; % input vector
y = engine_propeller(u);
P = y(1); % power [hp]
T = y(2); % thrust [lbf]
eta = y(3); % propeller efficiency
PSFC = y(6); % fuel consumption [lb/hp/hr]


%%

v_fps = v*kt2fps; % TAS [ft/s]

% calculate Cl, Cd, D
Cl = 2*w/(rho*S*v_fps^2);
Cd = .024292 + 1.6647e-9*Cl + .071696*Cl^2; % Cd from drag polar
L = 0.5*rho*(v_fps^2)*S*Cl; % [lbf]
D = 0.5*rho*(v_fps^2)*S*Cd; % [lbf]

% calculate climb angle:
singamma = P*550*eta/(v_fps*w) - D/w;
gamma = asin( singamma ); % climb angle [rad]

% Extract Vx and Vy from V:
vx = v_fps*cos(gamma); % horizantal speed [ft/s]
vy = v_fps*sin(gamma); % vertical speed [ft/s]

% find state derivitives for Vx and Vy:
vx_dot = ( T-D-L*sin(gamma) )/(w*lb2slg); % horizantal acceleration [ft/s^2]
vy_dot = (L-w)/(w*lb2slg); % vertical acceleration [ft/s^2]

% Update Vx and Vy:
vx2fps =  vx + vx_dot*t_step; % new horizantal speed [ft/s]
vx2 = vx2fps/kt2fps; % new horizantal speed [nm/s]
vy2 = vy + vy_dot*t_step; % new vertical speed [ft/s]

x_dot = vx2; % horizantal speed [nm/s]
x_dot = x_dot/3600; % horizantal speed [nm/s]
h_dot = vy2; % ROC [ft/s]
w_dot = -PSFC*P; % [lb/hr]
w_dot = w_dot/3600; % [lb/s]

ROC = h_dot;

SR = -vx/w_dot; % specific range [ft/lb]
SR = SR/kt2ft; % specific range [nm/lb]

% update states:
t2 = t + t_step; % time [s]
x2 = x + x_dot*t_step; % new distance [nautical miles]
h2 = h + h_dot*t_step; % new altitude [ft]
w2 = w + w_dot*t_step; % new weight [lbs]
v2 = sqrt(vx2fps^2 + vy2^2); % airspeed [ft/s]
v2 = v2/kt2fps; % airspeed [kts]
MAP2 = MAP; % initial Manifold pressure

state2 = [t2; x2; h2; w2; v2; MAP2];
end




