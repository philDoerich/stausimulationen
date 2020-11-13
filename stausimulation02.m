clear all;
close all;

% löse die Transportgleichung: u_t + u_u_x = 0
L = 230;            % 5 Kilometer Streckenlänge (Meter)
T = 40;             % Time of interest (Sekunden) 
CFL = 0.9;
Imax = 401;         % Anzahl Stützstellen
VMax = 36;         % Maximale Geschwindigkeit (bei RhoMin)  (m/s) 
RhoMax = 1.85;      % Maximale Anzahl an Autos innerhalb eines Kilometers   Auto / Meter
RhoStart = 0.81;     % Minimale Anzahl an Autos innerhalb eines Kilometers (Maximaler Wert bei Einhaltung des Mindestabstands) Auto/Meter

x = linspace(0, L, Imax);
deltaX = L/(Imax - 1);

u0 = zeros(Imax, 1);

% set initial conditions
% Gauss pulse centered at 4
for i = 1: Imax
    if(x(i)<2)
        u0(i,1) = RhoStart;
    elseif(x(i)<30)
        u0(i,1)= (1 + 0.01*sin(0.1*pi*(x(i) - 2)))* RhoStart;
    else
        u0(i,1)= RhoStart;
    end
end

maxU0 = -inf;

for i =1: Imax
    maxU0 = max(u0(i, 1), maxU0);
end

c = ((-2)*maxU0)*VMax/RhoMax + VMax;

deltaT = CFL*deltaX/abs(c);
Nmax = ceil(T/deltaT) + 1;

u = zeros(Imax, Nmax);
t = linspace(0, T, Nmax);


for i = 1: Imax
   u(i,1) = u0(i,1);
end


for n = 2: Nmax
    %backward differences
    %periodic boundary conditions
    
    u(1, n) = u(1, n - 1) + deltaT/deltaX*c*u(Imax, n - 1)*(u(1, n - 1)-u(Imax, n - 1));
    %u(1, n) = 1;
    u(Imax, n) = u(1, n - 1) + deltaT/deltaX*c*u(1, n - 1)*(u(1, n - 1)-u(Imax, n - 1));
    %u(Imax, n) = 1;
    for i = 2: Imax-1
        c = ((-2)*u(i, n - 1)*VMax)/RhoMax + VMax;
        
        cplus = max(c, 0);
        cminus = min(c, 0);
        
        u(i, n) = u(i, n - 1) - deltaT/deltaX*cplus*(u(i, n - 1) - u(i - 1, n - 1)) - deltaT/deltaX*cminus*(u(i + 1, n - 1) - u(i, n - 1));
        u(i, n) = min(u(i, n), RhoMax);
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