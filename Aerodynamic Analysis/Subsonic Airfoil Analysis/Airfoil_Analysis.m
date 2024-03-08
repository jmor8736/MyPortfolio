% Jackson Morgan
% Airfoil Analysis Project

clear
clc

% Hard code experimental alpha data
alphax = [-18 -16 -14 -12 -10 -8 -6 -4 -2 0 ...
    2 4 6 8 10 12 14 16 18 20];

% Import xfoil data
xfoil = readtable('NACA4415data.txt');

% define NACA
NACA = 4415;

%% Part 1 - Grid Independence Analysis

% define constants for analysis
alpha = 2;
Q = 1.25;

% create loop until GCI is lower than desired value
i = 0;
GCI2 = 1;
while GCI2 >= 0.01
    % define current iteration's numbers of panels and spacing
    N1 = 10*2^i;
    N2 = 10*2^(i+1);
    N3 = 10*2^(i+2);
    h1 = 1/N1;
    h2 = 1/N2;
    h3 = 1/N3;
    
    % evaluate dependent terms (cl) at each N
    [phi1, m1] = lumped_vortex(NACA,N1,alpha);
    [phi2, m2] = lumped_vortex(NACA,N2,alpha);
    [phi3, m3] = lumped_vortex(NACA,N3,alpha);
    
    % obtain ratios of approximate error
    r12 = h1/h2;
    r23 = h2/h3;
    r = r12;
    
    % solve for n (order of accuracy)
    n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
    
    % calculate GCI
    GCI3 = (Q/(r23^n-1))*((phi3-phi2)/phi3);
    GCI2 = (r23^n)*GCI3*(phi3/phi2);
   
    % increase i for next iteration
    i = i+1;
end

% calculate relative error
relerr = abs(phi2-phi1);

NGCI = N2;
% Display Results
disp(['Grid Independent Solution at N = ',num2str(N2),' Panels']);
disp(['Expected Relative Error = ',num2str(relerr)]);

% Generate Values for plot from N=40 to N=1280
for i = 1:6
    % define current iteration's numbers of panels and grid spacing
    N1 = 10*2^(i-1);
    N2 = 10*2^(i);
    N3 = 10*2^(i+1);
    h1 = 1/N1;
    h2 = 1/N2;
    h3(i) = 1/N3;
    
    % evaluate dependent terms (cl) for each N
    [phi1, m1] = lumped_vortex(NACA,N1,alpha);
    [phi2, m2] = lumped_vortex(NACA,N2,alpha);
    [phi3, m3] = lumped_vortex(NACA,N3,alpha);
    
    % obtain ratios of approximate error
    r12 = h1/h2;
    r23 = h2/h3(i);
    r = r12;
    
    % solve for n (order of accuracy)
    n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
    
    % calculate GCI
    GCI3(i) = (Q/(r23^n-1))*((phi3-phi2)/phi3);
end

% plot results
figure(1)
loglog(h3,GCI3);
grid on
title('Log-Log of Grid Conversion Index vs Panel Length')
xlabel('Panel Length')
ylabel('Grid Convergence Index')


%% Part 2 - Lift Coefficients

% extract xfoil data
alphaxf = table2array(xfoil(:,1));
clxf = table2array(xfoil(:,2));

% call on Vortex Panel function
NACA = 4415;
alpha2 = -18:0.1:20;
N = NGCI;
cl = zeros(length(alpha2),1);
for i = 1:length(alpha2)
    [cl(i), cmle] = lumped_vortex(NACA,N,alpha2(i));
end

% call thin airfoil theory function
clt = zeros(length(alpha));
for i = 1:length(alpha2)
    [clt(i), cmc4t] = thin_airfoil(NACA,alpha2(i));
end

% hard code experimental cl data
clx = [-.74 -.86 -.9 -.82 -.6 -.4 -.2 0 .2 .4 .6 ... 
    .8 1.05 1.2 1.38 1.5 1.56 1.62 1.58 1.31];

% Plot results
figure(2)
plot(alphaxf,clxf,'k',alpha2,cl,'r',alpha2,clt,'m',alphax,clx,'bx')
grid on
title('NACA 4415 Lift Coefficient vs Angle of Attack');
xlabel('Angle of Attack');
ylabel('Lift Coefficient');
legend('xfoil','Vortex Panel','Thin Airfoil','Experimental Data', ...
    'location','southeast');


%% Part 3 - Moment Coefficients

% extract xfoil cm data
cmxf = table2array(xfoil(:,5));

% call vortex panel function
cmc4 = zeros(length(alpha2),1);
for i = 1:length(alpha2)
    [cl(i), cmle] = lumped_vortex(NACA,N,alpha2(i));
    cmc4(i) = cmle+(1/4)*cl(i); % obtain moment about c/4 from leading edge moment
end

% call thin airfoil theory function
for i = 1:length(alpha2)
    [clt, cmc4t(i)] = thin_airfoil(NACA,alpha2(i));
end

% hard coding experimental data
cmx = [0 -.025 -.05 -.08 -.1 -.095 -.093 -.092 -.091 -.09 ...
    -.085 -.084 -.083 -.081 -.08 -.07 -.065 -.065 -.07 -.075];

% Plotting and labeling
figure(3)
plot(alphaxf,cmxf,'k',alpha2,cmc4,'r',alpha2,cmc4t,'m',alphax,cmx,'bx');
grid on
ylim([-0.2, 0])
title('NACA 4415 Moment Coefficient vs Angle of Attack');
xlabel('Angle of Attack');
ylabel('Quarter Chord Moment Coefficient');
legend('xfoil','Vortex Panel','Thin Airfoil','Experimental Data', ...
    'location','southeast');


%% Part 4 - Drag Polar

% extract xfoil cd data
cdxf = table2array(xfoil(:,3));

% Viscous Analysis for lumped vortex and thin airfoil theory
Retr = 5e5;
Rec = 9e6;
cdl = 1.328; % cd for laminar boundary layer
cdt = 0.0315; % cd for turbulent boundary layer
cdlv = zeros(length(alpha2),1);
for i = 1:length(alpha2) 
    % cd doesn't change with i, I'm just building an array for plotting
    cdlv(i) = 2*((cdl/sqrt(Retr))*(Retr/Rec) + (cdt/Rec^(1/7)) ...
        - (cdt/Retr^(1/7))*(Retr/Rec));
end
cdta = cdlv;

% hard code experimental drag polar data
cls = [-1 -.9 -.8 -.7 -.6 -.5 -.4 -.3 -.2 -.1 0 .1 .2 ...
    .3 .4 .5 .6 .7 .8 .9 1 1.1 1.2 1.3 1.4 1.5];
cds = [.011 .01 .0091 .0088 .0082 .0078 .0075 .0072 .007 ...
    .0068 .0065 .0064 .0062 .0063 .0064 .0068 .0069 .0073 .008 .0088 ...
    .0098 .0108 .012 .0134 .015 .0172];

% plot results
figure(4)
plot(cdxf,clxf,'k',cdlv,cl,'r',cdta,cl,'m',cds,cls,'bx');
grid on
xlabel('Drag Coefficient');
ylabel('Lift Coefficient');
title('NACA 4415 Drag Polar');
legend('xfoil','Lumped Vortex','Thin Airfoil','Experimental Data', ...
    'location','southeast');





