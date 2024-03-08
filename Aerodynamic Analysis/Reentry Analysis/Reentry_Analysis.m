% Jackson Morgan
% AERO 405
% Final Project

clc, clear

%% Hard Coding Experimental Data
data = [968, 1045, 1431, 1505, 1538, 1553, 1572, 1602, 1637, 1678, 1880;...
    275000, 180000, 150000, 140000, 110000, 100000, 88000, 75000, 48000,...
    25000, 0; 25.6, 29.2, 12.6, 6, 4, 3, 2, 1.1, .7, .42, .05];
% row 1 - time [s], row 2 - altitude [ft], row 3 - Mach number

% temperature and pressure from 1976 standard atmosphere model:
T1 = [189.010, 259.831, 267.066, 258.532, 232.928, 227.130, 223.472,...
    219.510, 216.650, 238.620, 288.150]; % K
P1 = [0.450413, 40.6917, 130.495, 193.950, 692.301, 1090.16, 1898.28,...
    3497.80, 12767.4, 37600.9, 101325]; % Pa


%% Thermally perfect calculations

for i = 1:length(T1)
    if data(1,i) > 1 % points where normal shock occurs

        % get downstream conditions after shock:
        [M2(i), P2(i), T2(i)] = tpg_normal_shock(data(3,i),P1(i),T1(i));

        % isentropic relation after shock to get stagnation state:
        [~, Tstag_tpg(i)] = tpg_isentropic_flow(M2(i), P2(i), T2(i), 0);

    else % points where no normal shock occurs

        % just use isentropic relation for when there is no shock:
        [~, Tstag_tpg(i)] = tpg_isentropic_flow(data(3,i), P1(i), T1(i), 0);
    end
end


%% Calorically perfect calculations

gamma = 1.4; % constant

for i = 1:length(T1)
    if data(3,i) > 1 % points where normal shock occurs

        % obtain temperature ratio across shock:
        Traysh_norm(i) = ( (2*gamma*data(3,i)^2 - (gamma-1))*...
            ((gamma-1)*data(3,i)^2 + 2) )/( (gamma+1)^2 *data(3,i)^2 );

        % temperature downstream of shock:
        T2_cpg(i) = T1(i)*Traysh_norm(i);

        % Mach number downstream of shock (to be used for isen relation)
        M2_cpg(i) = sqrt( ( (gamma-1)*data(3,i)^2 + 2 )/...
            ( 2*gamma*data(3,i)^2 - (gamma-1) ) );

        % Isentropic temperature ratio to stagnation state:
        Traysh_isen(i) = ( 1+((gamma-1)/2)*M2_cpg(i)^2 )^(-1);

        % stagnation temperature:
        Tstag_cpg(i) = T2_cpg(i)/Traysh_isen(i);

    else % points where no normal shock occurs

        % isentropic temperature ratio to get to stagnation state:
        Traysh_isen(i) = ( 1+((gamma-1)/2)*data(3,i)^2 )^(-1);

        % stagnation temperature:
        Tstag_cpg(i) = T1(i)/Traysh_isen(i);
    end
end


%% Generate Deliverables

% plot data:
figure(1), hold on, grid on
plot(data(1,:),Tstag_cpg,'linewidth',1.5,'color',[0 0.4470 0.7410],'marker','x')
plot(data(1,:),Tstag_tpg,'linewidth',1.5,'color',[0.3010 0.7450 0.9330],'marker','x')
title('Temperature vs Time'), xlabel('Mission Time [s]'), ylabel('Temperature [K]')
legend('Calorically Perfect Gas Model','Thermally Perfect Gas Model')

% calculate difference between model results:
delta = Tstag_cpg - Tstag_tpg;

% tabulate data:
T = table(data(1,:)', Tstag_tpg', Tstag_cpg', delta','VariableNames',{'Mission Time [s]',...
    'TPG Stag Temp [K]','CPG Stag Temp [K]','Difference in Results'});







