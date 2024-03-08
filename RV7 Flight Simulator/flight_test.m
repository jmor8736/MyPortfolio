% Jackson Morgan
% RV7 Simulator Flight Test Validation

clc
clear

% INITIALIZE STATE VECTOR:
t0 = 0;
x0 = 0;
h0 = 212; 
w0 = 1800;
v0 = 95; % climb velocity as defined by problem statement
MAP0 = 30; % wide-open throttle
state = [t0; x0; h0; w0; v0; MAP0];

i = 1;
t_step = 5;
delta_T = 0;

ROC(i) = 0;
SR(i) = 0;

% Initial Climb WOT for 2 minutes
while state(1,i) < 120
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

MAP = 23; % Normal climb MAP as defined by Problem
v = 115; % Normal Climb initial V as defined by problem
state(6,i) = MAP;
state(5,i) = v;


% Normal Climb until DOBRA
while state(2,i) < 6.4
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

DOBRA = i; % Index of DOBRA waypoint

% Climb to cruise height on way to KPRB
while state(3,i) < 8500
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

TOC1 = i; % Index for top of first climb

MAP = 20; % Cruise MAP as defined by Problem
state(6,i) = MAP;

% Cruise at 8500 ft until KPRB
while state(2,i) < 6.4+23.2
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end

KPRB = i; % Index at KPRB waypoint

% Cruise at 8500 ft until BONIT
while state(2,i) < 6.4+23.2+7
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end

BONIT = i; % Index at BONIT waypoint


%% CALCULATE DESCENT DISTANCE WITH TEMP ITERATOR AND STATE VECTOR 11111111

Poopoo = state; % temp state vector
j = i; % temp iterator
MAP = 17; % Descent MAP as defined by Problem
Poopoo(6,j) = MAP; % set MAP to descent MAP

% Measure Distance to Descend to Altitude of 
while Poopoo(3,j) > 374
    [Poopoo(:,j+1),~,~] = climb_integrator(Poopoo(:,j), t_step, delta_T);
    j = j + 1;
end
delta_x_descent1 = Poopoo(2,j)-state(2,i);

%% BACK TO ACTUAL FLIGHT SIMULATION


% Cruise at 8500 ft On way to KKIC until initial descent point
while state(2,i) < 6.4+23.2+7+35.3 - delta_x_descent1
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end


%% CALCULATE DESCENT DISTANCE WITH TEMP ITERATOR AND STATE VECTOR 22222222

Poopoo = state; % temp state vector
j = i; % temp iterator
MAP = 17; % Descent MAP as defined by Problem
Poopoo(6,j) = MAP; % set MAP to descent MAP

% Measure Distance to Descend to Altitude of 
while Poopoo(3,j) > 374
    [Poopoo(:,j+1),~,~] = climb_integrator(Poopoo(:,j), t_step, delta_T);
    j = j + 1;
end
delta_x_descent1 = Poopoo(2,j)-state(2,i);

%% BACK TO ACTUAL FLIGHT SIMULATION


% Cruise at 8500 ft On way to KKIC until initial descent point
while state(2,i) < 6.4+23.2+7+35.3 - delta_x_descent1
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end

initial_descent1 = i; % Index for first initial descent point

MAP = 17; % descent MAP as defined by Problem Statement 
state(6,i) = MAP;

% Descend to KKIC
while state(3,i) > 374
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

KKIC = i; % Index at KKIC touch and go waypoint

MAP = 30; % wide open throttle for initial climb
v = 95; % initial V for WOP climb
state(6,i) = MAP;
state(5,i) = v;

endclimb = state(1,i) + 120; % two minutes after climb starts
% Climb at WOT on way to KPRB
while state(1,i) < endclimb
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

MAP = 23; % Normal climb MAP as defined by Problem
v = 115; % initial V for normal climb
state(6,i) = MAP;

% Normal Climb until cruise altitude
while state(3,i) < 7500
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

TOC2 = i; % index for top of climb 2

MAP = 20; % Cruise MAP as defined in problem
state(6,i) = MAP;

% Cruise at 7500 ft to BONIT
while state(2,i) < 6.4+23.2+7+35.3+35.3
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end

BONIT2 = i; % Index at second BONIT waypoint

% Cruise at 7500 ft to KPRB
while state(2,i) < 6.4+23.2+7+35.3+35.3+7
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end

KPRB2 = i; % Index at second KPRB waypoint


%% CALCULATE DESCENT DISTANCE WITH TEMP ITERATOR AND STATE VECTOR 11111111

Peepee = state;
k = i;

MAP = 17; % descent MAP
Peepee(6,k) = MAP;

while Peepee(3,k) > 3500
    [Peepee(:,k+1),~,~] = climb_integrator(Peepee(:,k), t_step, delta_T);
    k = k + 1;
end

delta_x_descent2 = Peepee(2,k)-state(2,i); % x-distance for descent to 3500 ft

%% BACK TO ACTUAL FLIGHT PLAN


% Cruise to Initial descent Point
while state(2,i) < 6.4+23.2+7+35.3+35.3+7+21.8 - delta_x_descent2
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end


%% CALCULATE DESCENT DISTANCE WITH TEMP ITERATOR AND STATE VECTOR 22222222

Peepee = state;
k = i;

MAP = 17; % descent MAP
Peepee(6,k) = MAP;

while Peepee(3,k) > 3500
    [Peepee(:,k+1),~,~] = climb_integrator(Peepee(:,k), t_step, delta_T);
    k = k + 1;
end

delta_x_descent2 = Peepee(2,k)-state(2,i); % x-distance for descent to 3500 ft

%% BACK TO ACTUAL FLIGHT PROFILE


% Cruise to Initial descent Point
while state(2,i) < 6.4+23.2+7+35.3+35.3+7+21.8 - delta_x_descent2
    [Power, state(:,i+1), ~,SR(i+1)] = cruise_integrator(state(:,i), t_step, delta_T);
    ROC(i+1) = 0;
    i = i + 1;
end

initial_descent2 = i; % Index for second initial descent point

MAP = 17; % descent MAP as defined by Problem
state(6,i) = MAP;

% Descend to Final Point
while state(3,i) > 3500
    [state(:,i+1),ROC(i+1),SR(i+1)] = climb_integrator(state(:,i), t_step, delta_T);
    i = i + 1;
end

MORRO = i; % Index for end of flight (Morro Rock)

state(1,:) = state(1,:)/60; % convert t from s to min

indecies = [DOBRA,TOC1,KPRB,BONIT,initial_descent1,KKIC,TOC2,BONIT2,KPRB2...
    ,initial_descent2,MORRO];

save('State.mat','state')
save('SR.mat','SR')
save('indecies.mat','indecies')


% make table to display 

% Plotting Profiles:
figure(1)

subplot(4,1,1)
plot(state(1,:),state(2,:),'linewidth',1.5)
hold on
scatter(state(1,DOBRA),state(2,DOBRA),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB),state(2,KPRB),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT),state(2,BONIT),'filled','MarkerFaceColor','k');
scatter(state(1,KKIC),state(2,KKIC),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT2),state(2,BONIT2),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB2),state(2,KPRB2),'filled','MarkerFaceColor','k');
scatter(state(1,MORRO),state(2,MORRO),'filled','MarkerFaceColor','k');
ylabel('Distance [nm]')
grid on

subplot(4,1,2)
plot(state(1,:),state(3,:),'linewidth',1.5)
hold on
scatter(state(1,DOBRA),state(3,DOBRA),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB),state(3,KPRB),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT),state(3,BONIT),'filled','MarkerFaceColor','k');
scatter(state(1,KKIC),state(3,KKIC),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT2),state(3,BONIT2),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB2),state(3,KPRB2),'filled','MarkerFaceColor','k');
scatter(state(1,MORRO),state(3,MORRO),'filled','MarkerFaceColor','k');
ylabel('Altitude [ft]')
grid on

subplot(4,1,3)
plot(state(1,:),state(5,:),'linewidth',1.5)
hold on
scatter(state(1,DOBRA),state(5,DOBRA),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB),state(5,KPRB),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT),state(5,BONIT),'filled','MarkerFaceColor','k');
scatter(state(1,KKIC),state(5,KKIC),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT2),state(5,BONIT2),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB2),state(5,KPRB2),'filled','MarkerFaceColor','k');
scatter(state(1,MORRO),state(5,MORRO),'filled','MarkerFaceColor','k');
ylabel('Airspeed [kts]')
grid on

subplot(4,1,4)
plot(state(1,:),state(4,:),'linewidth',1.5)
hold on
scatter(state(1,DOBRA),state(4,DOBRA),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB),state(4,KPRB),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT),state(4,BONIT),'filled','MarkerFaceColor','k');
scatter(state(1,KKIC),state(4,KKIC),'filled','MarkerFaceColor','k');
scatter(state(1,BONIT2),state(4,BONIT2),'filled','MarkerFaceColor','k');
scatter(state(1,KPRB2),state(4,KPRB2),'filled','MarkerFaceColor','k');
scatter(state(1,MORRO),state(4,MORRO),'filled','MarkerFaceColor','k');
ylabel('Weight [lbs]')
xlabel('Time [min]')
grid on

spec_range = [SR(indecies(1)), SR(indecies(2)), SR(indecies(3)), ...
    SR(indecies(4)), SR(indecies(5)), SR(indecies(6)), ...
    SR(indecies(7)), SR(indecies(8)), SR(indecies(9)), ...
    SR(indecies(10)), SR(indecies(11))];

figure(2)
plot(state(1,:),SR,'linewidth',1.5,'color',[0 0.4470 0.7410])
grid on
hold on
xlabel('Time [min]')
ylabel('Specific Range [nm/lb]')
x1 = [state(1,(indecies(1))) state(1,(indecies(1)))];
y1 = [0 4.5];
line(x1,y1,'linestyle','--','color','k')
x2 = [state(1,(indecies(2))) state(1,(indecies(2)))];
y2 = [0 4.5];
line(x2,y2,'linestyle','--','color','k')
x3 = [state(1,(indecies(3))) state(1,(indecies(3)))];
y3 = [0 4.5];
line(x3,y3,'linestyle','--','color','k')
x4 = [state(1,(indecies(4))) state(1,(indecies(4)))];
y4 = [0 4.5];
line(x4,y4,'linestyle','--','color','k')
x5 = [state(1,(indecies(5))) state(1,(indecies(5)))];
y5 = [0 4.5];
line(x5,y5,'linestyle','--','color','k')
x6 = [state(1,(indecies(6))) state(1,(indecies(6)))];
y6 = [0 4.5];
line(x6,y6,'linestyle','--','color','k')
x7 = [state(1,(indecies(7))) state(1,(indecies(7)))];
y7 = [0 4.5];
line(x7,y7,'linestyle','--','color','k')
x8 = [state(1,(indecies(8))) state(1,(indecies(8)))];
y8 = [0 4.5];
line(x8,y8,'linestyle','--','color','k')
x9 = [state(1,(indecies(9))) state(1,(indecies(9)))];
y9 = [0 4.5];
line(x9,y9,'linestyle','--','color','k')
x10 = [state(1,(indecies(10))) state(1,(indecies(10)))];
y10 = [0 4.5];
line(x10,y10,'linestyle','--','color','k')
x11 = [state(1,(indecies(11))) state(1,(indecies(11)))];
y11 = [0 4.5];
line(x11,y11,'linestyle','--','color','k')

figure(3)
plot(state(3,:),SR,'linewidth',1.5,'color',[0 0.4470 0.7410])
grid on
hold on
xlabel('Altitude [ft]')
ylabel('Specific Range [nm/lb]')
% x1 = [state(3,(indecies(1))) state(3,(indecies(1)))];
% y1 = [0 4.5];
% line(x1,y1,'linestyle','--','color','k')
% x2 = [state(3,(indecies(2))) state(3,(indecies(2)))];
% y2 = [0 4.5];
% line(x2,y2,'linestyle','--','color','k')
% x3 = [state(3,(indecies(3))) state(3,(indecies(3)))];
% y3 = [0 4.5];
% line(x3,y3,'linestyle','--','color','k')
% x4 = [state(3,(indecies(4))) state(3,(indecies(4)))];
% y4 = [0 4.5];
% line(x4,y4,'linestyle','--','color','k')
% x5 = [state(3,(indecies(5))) state(3,(indecies(5)))];
% y5 = [0 4.5];
% line(x5,y5,'linestyle','--','color','k')
% x6 = [state(3,(indecies(6))) state(3,(indecies(6)))];
% y6 = [0 4.5];
% line(x6,y6,'linestyle','--','color','k')
% x7 = [state(3,(indecies(7))) state(3,(indecies(7)))];
% y7 = [0 4.5];
% line(x7,y7,'linestyle','--','color','k')
% x8 = [state(3,(indecies(8))) state(3,(indecies(8)))];
% y8 = [0 4.5];
% line(x8,y8,'linestyle','--','color','k')
% x9 = [state(3,(indecies(9))) state(3,(indecies(9)))];
% y9 = [0 4.5];
% line(x9,y9,'linestyle','--','color','k')
% x10 = [state(3,(indecies(10))) state(3,(indecies(10)))];
% y10 = [0 4.5];
% line(x10,y10,'linestyle','--','color','k')
% x11 = [state(3,(indecies(11))) state(3,(indecies(11)))];
% y11 = [0 4.5];
% line(x11,y11,'linestyle','--','color','k')













