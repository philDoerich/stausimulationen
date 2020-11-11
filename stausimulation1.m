clear all;
close all;

% löse die Transportgleichung: u_t + u_u_x = 0
L = 5;              % 5 Kilometer Streckenlänge
T = 0.02;           % Time of interest
c = 130;          % Velocity in u_t + u*u_x = 0
CFL = 0.9;
Imax = 601;         % Anzahl Stützstellen
VMax = 130;         % Maximale Geschwindigkeit (bei RhoMin)   
RhoMax = 200;       % Maximale Anzahl an Autos innerhalb eines Kilometers   
RhoMin = 14.28;     % Minimale Anzahl an Autos innerhalb eines Kilometers (Maximaler Wert bei Einhaltung des Mindestabstands) 

x = linspace(0, L, Imax);
deltaX = L/(Imax - 1);
deltaT = CFL*deltaX/abs(c);
Nmax = ceil(T/deltaT) + 1;

u = zeros(Imax, Nmax);
t = linspace(0, T, Nmax);

% set initial conditions
% Gauss pulse centered at 4
for i = 1: Imax
    if(x(i)<3)
        u(i,1) = RhoMin;
    elseif(x(i)<4)
        u(i,1)=(RhoMax - RhoMin)*sin(2*pi*(x(i) - 2)) + RhoMin;
    else
        u(i,1)= RhoMin;
    end
end

for n = 2: Nmax
    %backward differences
    %periodic boundary conditions
    u(1, n) = RhoMin;
    u(Imax, n) = RhoMin;
    
    for i = 2: Imax-1
        rho = ((-2)*max(u(i, n - 1), 170)*VMax)/RhoMax + VMax;
        
        cplus = max(rho, 0);
        cminus = min(rho, 0);
        
        u(i, n) =  u(i, n - 1) - deltaT/deltaX*cplus*(u(i, n - 1) - u(i - 1, n - 1)) - deltaT/deltaX*cminus*(u(i + 1, n - 1) - u(i, n - 1));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'numerical Solution u','NumberTitle','off')

plot(x, u(1:Imax,1), 'g','LineWidth',2);
axis tight
set(gca,'nextplot','replacechildren');
% Record the movie
for j = 1:Nmax
    plot(x, u(1:Imax,j), 'g','LineWidth',2);
    F(j) = getframe;
end
movie(F,0) % play x+1 times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Numerical vs. exact solution','NumberTitle','off')
hold on;

plot(x, u(1:Imax,1), x, u(1:Imax,  ceil(Nmax/6)),'-.r', x, u(1:Imax, 2* ceil(Nmax/6)), '--m', x, u(1:Imax, 3* ceil(Nmax/6)), ':b', x, u(1:Imax, 4* ceil(Nmax/6)), '-.g', x, u(1:Imax, 5* ceil(Nmax/6)), ':r', x, u(1:Imax, Nmax), '-.m');
%plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','u','NumberTitle','off')
C = contourf( u');
xlabel('x')
ylabel('t')

figure('Name','numerical Solution u','NumberTitle','off')
mesh( u');
view (33, 12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Geschwindigkeit','NumberTitle','off')
plot([RhoMax;0],[0;VMax])
%plot(x, u(1:Imax,1), 'g','LineWidth',2);
%axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%