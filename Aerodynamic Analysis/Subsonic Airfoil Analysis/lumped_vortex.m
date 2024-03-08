% Jackson Morgan
% Lumped Vortex Panel Code

function [cl, cmle] = lumped_vortex(NACA,N,alpha)
% Inputs:
%   NACA = 4 digit NACA Designation
%   N = number of panels/vorticies/collocation points
%   alpha = angle of attack
% 
% Outputs:
%   cl = lift coefficient
%   mle = moment about the leading edge

%% Obtaining Panel Coordinates

Uinf = 100;

% setting n so that we will get correct number of panels
n = N + 1; % n is number of points, N is number of panels

x = linspace(0,1,n); % x (dependant on n number of points)

% extracting values from the NACA number
NACA = num2str(NACA);
m = str2double(NACA(1))/100;
p = str2double(NACA(2))/10;

% Iterate through x and calculate camber (z values)
z = zeros(1,n);
for i = 1:n
    if x(i) < p
        z(i) = (m/p^2)*(2*p*x(i) - x(i)^2);
    else
        z(i) = (m/(1-p)^2)*(1 - 2*p + 2*p*x(i) - x(i)^2);
    end
end


%% Obtaining Geometry Terms

for i = 1:N
    % getting normal vectors
    dx(i) = x(i+1) - x(i);
    dz(i) = z(i+1) - z(i);
    dl(i) = sqrt((x(i+1)-x(i))^2 + (z(i+1)-z(i))^2);
    normx(i) = -dz(i)/dl(i);
    normz(i) = dx(i)/dl(i);

    % getting vortex and collocation coordinates
    midx(i) = 0.5*(x(i)+x(i+1)); % midpoint holders
    midz(i) = 0.5*(z(i)+z(i+1)); 
    xv(i) = 0.5*(x(i)+midx(i)); % vortex locations
    zv(i) = 0.5*(z(i)+midz(i));
    xc(i) = 0.5*(midx(i)+x(i+1)); % collocation points
    zc(i) = 0.5*(midz(i)+z(i+1));

end

% obtaining freestream velocity components
Ufree = Uinf*[cosd(alpha); sind(alpha)];


%% Solving for Induced Velocities at Collocation Points

for i = 1:N
    for j = 1:N
        r = sqrt((xc(i)-xv(j))^2 + (zc(i)-zv(j))^2);
        u(i,j) = (1/(2*pi))*((zc(i)-zv(j))/(r^2));
        v(i,j) = -(1/(2*pi))*((xc(i)-xv(j))/(r^2));
    end
end


%% Building and Solving System of Equations

for i = 1:N
    norm = [normx(i); normz(i)];
    B(i,1) = dot(-Ufree,norm);
    for j = 1:N
        U = [u(i,j); v(i,j)];
        A(i,j) = dot(U,norm);
    end
end

% solve system of equations for gamma values
gamma = A\B;
Gamma = sum(gamma);


%% Obtaining Outputs (cd & cmle)

% use gamma values to get lift coefficient
cl = 2*Gamma/Uinf;

% for loop to find moment from each panel
for i = 1:N
    xg(i) = -xv(i)*gamma(i);
end

% sum moments to get total moment and find the leading edge coefficient
XG = sum(xg);
cmle = 2*XG*cosd(alpha)/Uinf;

end



