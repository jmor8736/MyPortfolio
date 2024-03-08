% Jackson Morgan
% AERO 405
% Biconvex Airfoil Analysis Master File

clc
clear

% airfoil geometry information:
c = 1; % chord length
t = 0.1; % thickness as % of chord

%% Grid Convergence Analysis

% Perform a grid convergence analysis for each method by observing the
% drag coefficient at an angle of attack of 3 degrees:

alpha = deg2rad(3);
M = 2.13;
gamma = 1.4;



% linear potential method Grid Convergence Analysis:

Q = 1.25; % using three spacings to Q = 1.25
r = 2; % r will be constant, at 2

% initial three spacings before we start making them smaller:
h1 = 0.1;
[x,y] = biconvex(c,h1,t);
[~,phi1,~] = lin_pot(x,y,alpha,M);
h2 = h1/r;
[x,y] = biconvex(c,h2,t);
[~,phi2,~] = lin_pot(x,y,alpha,M);
h3 = h2/r;
[x,y] = biconvex(c,h3,t);
[~,phi3,~] = lin_pot(x,y,alpha,M);

% initial order of accuracy and GCIs
n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
eps3 = (phi3-phi2)/(r^n-1);
GCI3 = Q*abs(eps3/phi3);
GCI2 = (Q/(r^n-1))*abs((phi2-phi1)/phi2);

% decreasing spacing size until answers are grid independent:
while abs(GCI2-GCI3) > 10e-4
    h1 = h2;
    [x,y] = biconvex(c,h1,t);
    [~,phi1,~] = lin_pot(x,y,alpha,M);
    h2 = h3;
    [x,y] = biconvex(c,h2,t);
    [~,phi2,~] = lin_pot(x,y,alpha,M);
    h3 = h2/r;
    [x,y] = biconvex(c,h3,t);
    [~,phi3,~] = lin_pot(x,y,alpha,M);

    n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
    eps3 = (phi3-phi2)/(r^n-1);
    GCI3 = Q*abs(eps3/phi3);
    GCI2 = (Q/(r^n-1))*abs((phi2-phi1)/phi2);
end
h_lin = h2; % spacing required for a grid independent solution



% Busemann's method Grid Convergence Analysis

% initial three spacings before we start making them smaller:
h1 = 0.1;
[x,y] = biconvex(c,h1,t);
[phi1,~,~] = buse(x,y,alpha,M,gamma);
h2 = h1/r;
[x,y] = biconvex(c,h2,t);
[phi2,~,~] = buse(x,y,alpha,M,gamma);
h3 = h2/r;
[x,y] = biconvex(c,h3,t);
[phi3,~,~] = buse(x,y,alpha,M,gamma);

% initial order of accuracy and GCIs
n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
eps3 = (phi3-phi2)/(r^n-1);
GCI3 = Q*abs(eps3/phi3);
GCI2 = (Q/(r^n-1))*abs((phi2-phi1)/phi2);

% decreasing spacing size until answers are grid independent:
while abs(GCI2-GCI3) > 10e-4
    h1 = h2;
    [x,y] = biconvex(c,h1,t);
    [phi1,~,~] = buse(x,y,alpha,M,gamma);
    h2 = h3;
    [x,y] = biconvex(c,h2,t);
    [phi2,~,~] = buse(x,y,alpha,M,gamma);
    h3 = h2/r;
    [x,y] = biconvex(c,h3,t);
    [phi3,~,~] = buse(x,y,alpha,M,gamma);

    n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
    eps3 = (phi3-phi2)/(r^n-1);
    GCI3 = Q*abs(eps3/phi3);
    GCI2 = (Q/(r^n-1))*abs((phi2-phi1)/phi2);
end
h_bus = h2; % spacing required for a grid independent solution



% Shock and Expansion method Grid Convergence Analysis

% initial three spacings before we start making them smaller:
h1 = 0.1;
[x,y] = biconvex(c,h1,t);
[phi1,~,~,~,~] = shock_expan(x,y,alpha,M,gamma);
h2 = h1/r;
[x,y] = biconvex(c,h2,t);
[phi2,~,~,~,~] = shock_expan(x,y,alpha,M,gamma);
h3 = h2/r;
[x,y] = biconvex(c,h3,t);
[phi3,~,~,~,~] = shock_expan(x,y,alpha,M,gamma);

% initial order of accuracy and GCIs
n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
eps3 = (phi3-phi2)/(r^n-1);
GCI3 = Q*abs(eps3/phi3);
GCI2 = (Q/(r^n-1))*abs((phi2-phi1)/phi2);

% decreasing spacing size until answers are grid independent:
while abs(GCI2-GCI3) > 10e-4
    h1 = h2;
    [x,y] = biconvex(c,h1,t);
    [phi1,~,~] = shock_expan(x,y,alpha,M,gamma);
    h2 = h3;
    [x,y] = biconvex(c,h2,t);
    [phi2,~,~] = shock_expan(x,y,alpha,M,gamma);
    h3 = h2/r;
    [x,y] = biconvex(c,h3,t);
    [phi3,~,~] = shock_expan(x,y,alpha,M,gamma);

    n = (1/log(r))*log((phi2-phi1)/(phi3-phi2));
    eps3 = (phi3-phi2)/(r^n-1);
    GCI3 = Q*abs(eps3/phi3);
    GCI2 = (Q/(r^n-1))*abs((phi2-phi1)/phi2);
end
h_exp = h2; % spacing required for a grid independent solution


%% Obtaining Coefficient Values using Necessary spacing for each method

% neccessary constants:
M = 2.13;
gamma = 1.4;

%alphad = 5;
alphad = -10:28;
alpha = deg2rad(alphad);

% Linear Potential:
[x,y] = biconvex(c,h_lin,t);
for i = 1:length(alpha)
    [Cll(i),Cdl(i),Cml(i)] = lin_pot(x,y,alpha(i),M);
end

% Busemann's:
[x,y] = biconvex(c,h_bus,t);
for i = 1:length(alpha)
    [Clb(i),Cdb(i),Cmb(i)] = buse(x,y,alpha(i),M,gamma);
end

% Oblique Shock and Expansion Theory:
alphads = -10:13; % can't use shock and expansion theory past 13 degrees 
                  % because we start getting complex numbers
alphas = deg2rad(alphads);
[x,y] = biconvex(c,h_exp,t);
for i = 1:length(alphas)
    [Cls(i),Cds(i),Cms(i),prayshu,prayshl] = shock_expan(x,y,alphas(i),M,gamma);
end


% Experimental Results

% Cm = Cm - leading edge moment coefficient
% Cp = Cl - Lift coeficient
% Cr = Cd - drag coefficient

% array for angles of attack for every two degrees from -10 to 28
alphae = -10:2:28;

% hard coding data of 100Cp, 100Cr, and 100Cm (average values for angles of
% attack with two data points)
Cp = [-16.25, -13, -9.8, -6.05, -3.25, -.05, 2.95, 6.25, 9.65, 13.2,...
    16.65, 20.4, 23.8, 27.15, 30.6, 34.05, 37.5, 39.85, 44.4, 48.1];
Cr = [4.91, 3.725, 2.83, 2.2, 1.865, 1.755, 1.89, 2.16, 2.84, 3.74,...
    4.91, 6.65, 8.35, 10.44, 12.955, 15.32, 17.81, 20.725, 23.375, 26.3];
Cmm = [-5.4, -4.85, -3.2, -1.9, -0.7, .05, 0.9, 2.05, 3.65, 5, 6.5, 7.4...
    9.7, 10.7, 13.25, 14.8, 17.5, 19.7, 22.45, 25.1];

% dividing arrays by 100 and relabeling to modern nomenclature for plotting
Cle = 2*Cp/100; 
Cde = 2*Cr/100;
% also multiplied by two because they divided by (rho*v^2) instead of
% (1/2rho*v^2)

% making Cm negative because modern practice is that a positive moment
% causes the airfoil to pitch up:
Cme = -2*Cmm/100;


%% Plotting Curves

figure(1)
plot(alphad,Cll,'b',alphad,Clb,'r',alphads,Cls,'g','linewidth',1.5)
hold on
plot(alphae,Cle,'k','marker','x','linewidth',1)
grid on
legend('Linear Potential','Busemanns','Shocks and Expansions*',...
    'Experimental','location','southeast')
xlabel('Angle of Attack (degrees)')
ylabel('Cl')
title('Lift Curve')

figure(2)
plot(alphad,Cdl,'b',alphad,Cdb,'r',alphads,Cds,'g','linewidth',1.5)
hold on
plot(alphae,Cde,'k','marker','x','linewidth',1)
grid on
legend('Linear Potential','Busemanns','Shocks and Expansions*',...
    'Experimental','location','north')
xlabel('Angle of Attack (degrees)')
ylabel('Cd')
title('Drag Curve')

figure(3)
plot(alphad,Cml,'b',alphad,Cmb,'r',alphads,Cms,'g','linewidth',1.5)
hold on
plot(alphae,Cme,'k','marker','x','linewidth',1)
grid on
legend('Linear Potential','Busemanns','Shocks and Expansions*',...
    'Experimental','location','northeast')
xlabel('Angle of Attack (degrees)')
ylabel('Cmle')
title('Leading Edge Moment Curve')

figure(4)
plot(x,y,'b',x,-y,'b','linewidth',1)
axis equal
grid on
title('Biconvex Airfoil Profile')
xlabel('Chord (x/c)')
ylabel('Thickness (y/c)')


%% Tabulating Data:


for i = 1:length(alphad)/2+1
    alph(i,1) = alphad(2*i-1);
    Cl(i,1) = Cll(2*i-1);
    Cl(i,2) = Clb(2*i-1);
    
    Cd(i,1) = Cdl(2*i-1);
    Cd(i,2) = Cdb(2*i-1);
   

    Cm(i,1) = Cml(2*i-1);
    Cm(i,2) = Cmb(2*i-1);
    
    if i <= length(alphas)/2
        Cl(i,3) = Cls(2*i-1);
        Cd(i,3) = Cds(2*i-1);
        Cm(i,3) = Cms(2*i-1);
    else
        Cl(i,3) = NaN;
        Cd(i,3) = NaN;
        Cm(i,3) = NaN;
    end
end
lift = round([alph,Cle',Cl],3);
lift_table = array2table(lift,"VariableNames",["AoA","Experimental",...
    "Linear","Busemann","Shock/Expansion"]);
drag = round([alph,Cde',Cd],3);
drag_table = array2table(drag,"VariableNames",["AoA","Experimental",...
    "Linear","Busemann","Shock/Expansion"]);
mom = round([alph,Cme',Cm],2);
mom_table = array2table(mom,"VariableNames",["AoA","Experimental",...
    "Linear","Busemann","Shock/Expansion"]);
% writetable(lift_table,'Lift_table.csv')
% writetable(drag_table,'Drag_table.csv')
% writetable(mom_table,'Moment_table.csv')

%% Functions:



% Biconvex Airfoil Geometry:
function [x,y] = biconvex(c,n,t)
% Inputs:
%   c - chord length
%   n - spacing
%   t - thickness as a % of the chord
% Outputs:
%   x - x coordinates of airfoil
%   y - y coordinates of upper surface of airfoil

x = 0:n:c;

% rest of geometry values calculated from chord length:
R = 2.5*c; % Radius of circle for surfaces
h = c/2; % x-coord of circle center
k = -R+t/2; % y-coord of circle center

% (h,k) is the center of the circle
for i = 1:length(x)
    y(i) = sqrt((R^2)-(x(i)-h)^2) + k;
end
end



% Linear Supersonic Small Disturbance Potential Theory:
function [Cl,Cd,Cm] = lin_pot(x,y,alpha,M)
% Inputs:
%   x - x-coordinates of airfoil
%   y - y-coordinates of airfoil
%   alpha - angle of attack
%   M - mach number
% Outputs:
%   Cl - lift coefficient
%   Cd - drag coefficient
%   Cm - leading edge moment coefficient

    % lift coefficient only dependent on AoA using this theory so panels
    % aren't required for calculation:
    Cl = 4*alpha/sqrt(M^2-1);

    % Cd element due to angle of attack:
    Cdalpha = alpha^2;

    % Cm in this case is only dependent on angle of attack:
    Cm = (-2*alpha/sqrt(M^2-1));

    for i = 1:length(x)-1
        % dx and dy from i-1 to be used later:
        dx = x(i+1) - x(i);
        dy = y(i+1) - y(i);

        % xm to be used for trapz
        xm(i) = (x(i+1)+x(i))/2;

        % differential Cd element due to thickness:
        h_prime = dy/dx; % rate of change of thickness 
        Cdtdx(i) = (h_prime^2); 
    end

    % Cd element due to thickness
      Cdt = trapz(xm,Cdtdx);

    % total Cd:
    Cd = (4/sqrt(M^2-1))*(Cdalpha + Cdt);
end



% Busemann's Method:
function [Cl,Cd,Cm] = buse(x,y,alpha,M,gamma)
% Inputs:
%   x - x-coordinates of airfoil
%   y - y-coordinates of airfoil
%   alpha - angle of attack
%   M - mach number
%   gamma - specific heat ratio
% Outputs:
%   Cl - lift coefficient
%   Cd - drag coefficient
%   Cm - leading edge moment coefficient

for i = 1:length(x)-1
    % array of x-location of the center of each panel (to be used for trapz
    % function later on):
    xm(i) = (x(i)+x(i+1))/2;
    ym(i) = (y(i)+y(i+1))/2;
end

    for i = 1:length(x)-1
    
        % dx and dy from i-1 to i for calculating thetas
        dx = x(i+1) - x(i);
        dy = y(i+1) - y(i);

        % calculate theta: to be used for delta and for integrations of
        % pressure later on
        thetau(i) = atan(dy/dx);
        thetal(i) = atan(dy/dx);

        % calculate upper and lower deltas
        deltau(i) = thetau(i) - alpha;
        deltal(i) = thetal(i) + alpha;

        % upper surface pressure coefficient:
        Cpu(i) = (2*deltau(i))/sqrt(M^2-1) + (((gamma+1)*M^4-4*M^2+4)/...
            (2*(M^2-1)^2))*deltau(i)^2;

        % lower surface pressure coefficient:
        Cpl(i) = (2*deltal(i))/sqrt(M^2-1) + (((gamma+1)*M^4-4*M^2+4)/...
            (2*(M^2-1)^2))*deltal(i)^2;

        % differential length of each portion
        % (top and bottom are the same)
        dl(i) = sqrt(dx^2 + dy^2);

        % differential normal force coefficient:
        dCnu(i) = Cpu(i)*cos(thetau(i));
        dCnl(i) = Cpl(i)*cos(thetal(i));

        % differential axial force coefficients:
        dCau(i) = Cpu(i)*sin(thetau(i));
        dCal(i) = Cpl(i)*sin(thetal(i));

         % differential moment about the leading edge:
        dCmlex(i) = (dCnu(i)-dCnl(i))*xm(i);
        dCmley(i) = dCau(i)*sin(abs(thetau(i)))*ym(i) - ...
            dCal(i)*sin(abs(thetal(i)))*ym(i);
    end
    
    % summing differential elements of normal snd axial force
    Cn = trapz(xm,dCnl)-trapz(xm,dCnu);
    Ca = trapz(ym,dCal)+trapz(ym,dCau);

    % Calculating total lift, drag, and moment
    Cl = Cn*cos(alpha) - Ca*sin(alpha);
    Cd = Ca*cos(alpha) + Cn*sin(alpha);
    Cm = trapz(xm,dCmlex)+trapz(ym,dCmley);
end



% Oblique Shock and Expansion Wave Theory:
function [Cl, Cd, Cm,prayshu,prayshl] = shock_expan(x,y,alpha,M,gamma)
% Inputs:
%   x - x-coordinates of airfoil
%   y - y-coordinates of airfoil
%   alpha - angle of attack
%   M - mach number
%   gamma - specific heat ratio
% Outputs:
%   Cl - lift coefficient
%   Cd - drag coefficient
%   Cm - leading edge moment coefficient

for i = 1:length(x)-1 % loop to iterate through panels

        % dx and dy from i-1 to i for calculating thetas
        dx = x(i+1) - x(i);
        dy = y(i+1) - y(i);
        
        % thetas of each panel on the upper and lower surfaces
        thetau(i) = atan(dy/dx);
        thetal(i) = atan(dy/dx);
            
        % calculating the shock angle to determine whether expansion or
        % oblique shock
        if i == 1
            % angle of first panel determined wrt freestream
            shocktypeu(i) = thetau(i) - alpha;
            shocktypel(i) = thetal(i) + alpha;
        else
            % anlges of following panels determined wrt previous panel
            shocktypeu(i) = thetau(i) - thetau(i-1);
            shocktypel(i) = thetal(i) - thetal(i-1);
        end
        
        % Upper Surface:

        % first determine pressures along upper surface:
        % if statement to determine whether to do oblique shock analysis or
        % expansion fan analysis
        if shocktypeu(i) > 0
            % in this scenario, an oblique shock occurs
            
            % ensure that theta is positive
            theta_up(i) = abs(shocktypeu(i));
            
            % M1 value is the Mach number before oblique shock (mach
            % number at the previous panel)
            if i == 1
                M1u(i) = M;
            else
                M1u(i) = M2u(i-1);
            end
            
            mu = asin(1/M1u(i));

            % calculate special angles:
            betamax = 0.5*acos((1/gamma)*( ((gamma+1)/2 - cos(2*mu))...
                - sqrt(gamma+1)*sqrt(((gamma+1)/2 - cos(2*mu))^2 + ...
                (gamma/4)*(3-gamma))));

            betasonic = asin(sqrt((1/(2*gamma))*((((gamma-3)/2)*...
                sin(mu)^2 + (gamma+1)/2) + sqrt(4*gamma*sin(mu)^4 +...
                (((gamma-3)/2)*sin(mu)^2 + (gamma+1)/2)^2))));

            % calculate beta
            delta = 1; % must be a weak shock
            lambda = sqrt( (M1u(i)^2-1)^2 - 3*(1+((gamma-1)/2)*M1u(i)^2)*...
                (1+((gamma+1)/2)*M1u(i)^2)*tan(theta_up(i))^2);
            chi = ( (M1u(i)^2-1)^3 - 9*(1+((gamma-1)/2)*M1u(i)^2)* ...
                (1+((gamma-1)/2)*M1u(i)^2+((gamma+1)/4)*M1u(i)^4)* ...
                tan(theta_up(i))^2 )/(lambda^3);
            num = M1u(i)^2 - 1 + 2*lambda*cos((4*pi*delta+acos(chi))/3);            
            den = 3*(1+((gamma-1)/2)*M1u(i)^2)*tan(theta_up(i));
            beta = atan(num/den);

            if beta > betasonic
                error('Alpha too high, flow detached')
            end

            % calculate normal shock equivalent
            Mn1u(i) = M1u(i)*sin(beta);
            % calculate normal shock equiv M2:
            Mn2u(i) = sqrt((2+(gamma-1)*Mn1u(i)^2)/(2*gamma*Mn1u(i)^2-...
                (gamma-1)));
            % calculate M2 from normal shock equiv M2:
            M2u(i) = Mn2u(i)/sin(beta-theta_up(i));

            % calculate pressure ratio across shock (P2/P1)
            prayshu(i) = ((2*gamma)/(gamma+1))*Mn1u(i)^2 -...
                (gamma-1)/(gamma+1);
        else
            % in this scenario, expansion waves form
            
            % create theta for flow angle change and ensure it is positive
            theta_up(i) = abs(shocktypeu(i));
            
            % for first panel, M1 is freestream M, for all other panels, M1
            % is M2 from the previous panel
            if i == 1
                M1u(i) = M;
            else
                M1u(i) = M2u(i-1);
            end
            
            % needed values to find nu1:
            lambda = sqrt((gamma-1)/(gamma+1));
            beta = sqrt(M1u(i)^2 - 1);
            
            % calculate nu1 and obtain nu2 from that
            nu_M1 = (1/lambda)*atan(lambda*beta) - atan(beta);
            nu_M2 = theta_up(i) + nu_M1;
            
            % define function for fzero:
            fun = @(B) (1/lambda)*atan(lambda*B) - atan(B) - nu_M2;

            % find initial guess for fzero using I M Hall approx:
                lambda = sqrt((gamma-1)/(gamma+1));

                k = (4/(3*pi))*(1 + 1/lambda);

                nu_inf = (pi/2)*(1/lambda - 1);

                eta_inf = (3*nu_inf/(1-lambda^2))^(2/3);

                a1 = eta_inf/2;
                a2 = ((3+8*lambda^2)/40)*eta_inf^2;
                a3 = ((-1+328*lambda^2+104*lambda^4)/2800)*eta_inf^3;

                d1 = a1-1- (a3-k)/(a2-k);
                d2 = a2-a1- ((a1-1)*(a3-k))/(a2-k);
                d3 = ((a1-k)*(a3-k))/(a2-k) - a2 + k;

                e1 = -1 - (a3-k)/(a2-k);
                e2 = (a3-k)/(a2-k);

                Y = (nu_M1/nu_inf)^(2/3);

                M_est = (1 + d1*Y + d2*Y^2 + d3*Y^3)/(1 + e1*Y + e2*Y^2);

                B_est = sqrt(M_est^2 - 1);

            % root finder to obtain beta and M2
            beta = fzero(fun,B_est);
            M2u(i) = sqrt(beta^2 + 1);

            % calculate pressure ratio across expansion using isentropic
            % relations (P2/P1)
            prayshu(i) = ( (2+((gamma-1)/2)*M1u(i)^2)/(2+((gamma-1)/2)*...
                M2u(i)^2) )^(gamma/(gamma-1));
        end

        % LOWER SURFACE:

        % next determine pressures along lower surface:
        % if statement to determine whether to do oblique shock analysis or
        % expansion fan analysis
        if shocktypel(i) > 0
            % in this scenario, an oblique shock occurs
            
            % ensure that theta is positive
            theta_low(i) = abs(shocktypel(i));
            
            % M1 value is the Mach number before oblique shock (mach
            % number at the previous panel)
            if i == 1
                M1l(i) = M;
            else
                M1l(i) = M2l(i-1);
            end
            
            mu = asin(1/M1l(i));

            % calculate special angles:
            betamax = 0.5*acos((1/gamma)*( ((gamma+1)/2 - cos(2*mu))...
                - sqrt(gamma+1)*sqrt(((gamma+1)/2 - cos(2*mu))^2 + ...
                (gamma/4)*(3-gamma))));

            betasonic = asin(sqrt((1/(2*gamma))*((((gamma-3)/2)*...
                sin(mu)^2 + (gamma+1)/2) + sqrt(4*gamma*sin(mu)^4 +...
                (((gamma-3)/2)*sin(mu)^2 + (gamma+1)/2)^2))));

            % calculate beta
            delta = 1; % must be a weak shock
            lambda = sqrt( (M1l(i)^2-1)^2 - 3*(1+((gamma-1)/2)*M1l(i)^2)*...
                (1+((gamma+1)/2)*M1l(i)^2)*tan(theta_low(i))^2);
            chi = ( (M1l(i)^2-1)^3 - 9*(1+((gamma-1)/2)*M1l(i)^2)* ...
                (1+((gamma-1)/2)*M1l(i)^2+((gamma+1)/4)*M1l(i)^4)* ...
                tan(theta_low(i))^2 )/(lambda^3);
            num = M1l(i)^2 - 1 + 2*lambda*cos((4*pi*delta+acos(chi))/3);
            
            den = 3*(1+((gamma-1)/2)*M1l(i)^2)*tan(theta_low(i));
            beta = atan(num/den);

            if beta > betasonic
                error('Alpha too high, flow detached')
            end

            % calculate normal shock equivalent
            Mn1l(i) = M1l(i)*sin(beta);
            % calculate normal shock equiv M2:
            Mn2l(i) = sqrt((2+(gamma-1)*Mn1l(i)^2)/(2*gamma*Mn1l(i)^2-...
                (gamma-1)));
            % calculate M2 from normal shock equiv M2:
            M2l(i) = Mn2l(i)/sin(beta-theta_low(i));

            % calculate pressure ratio across shock (P2/P1)
            prayshl(i) = ((2*gamma)/(gamma+1))*Mn1l(i)^2 -...
                (gamma-1)/(gamma+1);
        else
            % in this scenario, expansion waves form
            
            % create theta for flow angle change and ensure it is positive
            theta_low(i) = abs(shocktypel(i));
            
            % for first panel, M1 is freestream M, for all other panels, M1
            % is M2 from the previous panel
            if i == 1
                M1l(i) = M;
            else
                M1l(i) = M2l(i-1);
            end
            
            % needed values to find nu1:
            lambda = sqrt((gamma-1)/(gamma+1));
            beta = sqrt(M1l(i)^2 - 1);
            
            % calculate nu1 and obtain nu2 from that
            nu_M1 = (1/lambda)*atan(lambda*beta) - atan(beta);
            nu_M2 = theta_low(i) + nu_M1;
            
            % define function for fzero:
            fun = @(B) (1/lambda)*atan(lambda*B) - atan(B) - nu_M2;

            % find initial guess for fzero using I M Hall approx:
                lambda = sqrt((gamma-1)/(gamma+1));

                k = (4/(3*pi))*(1 + 1/lambda);

                nu_inf = (pi/2)*(1/lambda - 1);

                eta_inf = (3*nu_inf/(1-lambda^2))^(2/3);

                a1 = eta_inf/2;
                a2 = ((3+8*lambda^2)/40)*eta_inf^2;
                a3 = ((-1+328*lambda^2+104*lambda^4)/2800)*eta_inf^3;

                d1 = a1-1- (a3-k)/(a2-k);
                d2 = a2-a1- ((a1-1)*(a3-k))/(a2-k);
                d3 = ((a1-k)*(a3-k))/(a2-k) - a2 + k;

                e1 = -1 - (a3-k)/(a2-k);
                e2 = (a3-k)/(a2-k);

                Y = (nu_M2/nu_inf)^(2/3);

                M_est = (1 + d1*Y + d2*Y^2 + d3*Y^3)/(1 + e1*Y + e2*Y^2);

                B_est = sqrt(M_est^2 - 1);

            beta = fzero(fun,B_est);
            M2l(i) = sqrt(beta^2 + 1);

            % calculate pressure ratio across expansion using isentropic
            % relations praysh = (P2/P1)
            prayshl(i) = ( (2+((gamma-1)/2)*M1l(i)^2)/(2+((gamma-1)/2)*...
                M2l(i)^2) )^(gamma/(gamma-1));
        end
end

% Pressure Coefficients:
for i = 1:length(x)-1
    if i == 1
        free_rayshu(i) = prayshu(i);
        free_rayshl(i) = prayshl(i);
    else
        % free_raysh is P/P_inf
        free_rayshu(i) = free_rayshu(i-1)*prayshu(i);
        free_rayshl(i) = free_rayshl(i-1)*prayshl(i);
    end

    % upper and lower pressure coefficients:
    Cpu(i) = 2*(free_rayshu(i)-1)/((gamma)*M2u(i)^2);
    Cpl(i) = 2*(free_rayshl(i)-1)/((gamma)*M2l(i)^2);
%     Cpu(i) = (prayshu(i)-1)/((gamma)*M2u(i)^2);
%     Cpl(i) = (prayshl(i)-1)/((gamma)*M2l(i)^2);

    % array of x-location of the center of each panel (to be used for trapz
    % function later on):
    xm(i) = (x(i)+x(i+1))/2;
    ym(i) = (y(i)+y(i+1))/2;

    % differential normal force coefficient:
    dCnu(i) = Cpu(i)*cos(thetau(i));
    dCnl(i) = Cpl(i)*cos(thetal(i));

    % differential axial force coefficients:
    dCau(i) = Cpu(i)*sin(thetau(i));
    dCal(i) = Cpl(i)*sin(thetal(i));

    % differential moment about the leading edge:
    dCmlex(i) = (dCnu(i)-dCnl(i))*xm(i);
    dCmley(i) = dCau(i)*sin(abs(thetau(i)))*ym(i) - ...
        dCal(i)*sin(abs(thetal(i)))*ym(i);
end

% summing differential elements of normal snd axial force
% Cn = trapz(xm,dCnl)-trapz(xm,dCnu);
% Ca = trapz(ym,dCal)+trapz(ym,dCau);
Cn = trapz(xm,Cpl-dCnu);
Ca = trapz(ym,Cpl-Cpu);

% Calculating total lift, drag, and moment
Cl = Cn*cos(alpha) - Ca*sin(alpha);
Cd = Ca*cos(alpha) + Cn*sin(alpha);
Cm = trapz(xm,dCmlex)+trapz(ym,dCmley);
end

