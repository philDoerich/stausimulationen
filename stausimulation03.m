clear all;
close all;

% löse die Transportgleichung: u_t + u_u_x = 0
L = 800;            % Streckenlänge (Meter)
T = 60;             % Time of interest (Sekunden) 
CFL = 0.9;
Imax = 500;         % Anzahl Stützstellen
VMax = 36;          % Maximale Geschwindigkeit (bei RhoMin) (m/s) 
RhoMax = 1.85;      % Maximale Anzahl an Autos innerhalb eines Kilometers   Auto / Meter
RhoStart = 1.05;    % Minimale Anzahl an Autos innerhalb eines Kilometers (Maximaler Wert bei Einhaltung des Mindestabstands) Auto/Meter

x = linspace(0, L, Imax);
deltaX = L/(Imax - 1);

u0 = zeros(Imax, 1);

% set initial conditions
% Sinus so eingestellt dass er eine Schwingung im bereich (120, 240) macht
% (Freuquenz)
AnfangDichtestoerung = 120;
EndeDichtestoerung = 240;
freq = 2*pi/(EndeDichtestoerung-AnfangDichtestoerung);

for i = 1: Imax
    if(x(i)<AnfangDichtestoerung)
        u0(i,1) = RhoStart;
    elseif(x(i)<EndeDichtestoerung)
        u0(i,1)= (1 + 0.01*sin(freq*(x(i) - 120)))* RhoStart;
    else
        u0(i,1)= RhoStart;
    end
end

maxU0 = -inf;

for i =1: Imax
    maxU0 = max(u0(i, 1), maxU0);
end

%c = -(((-2)*maxU0)*VMax/RhoMax + VMax);
c = (((-2)*RhoMax)*VMax/RhoMax + VMax);
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
    c = (((-2)*u(Imax, n - 1)*VMax)/RhoMax + VMax);
    cplus = max(c, 0);
    
    c = (((-2)*u(2, n - 1)*VMax)/RhoMax + VMax);
    cminus = min(c, 0);
    
    u(1, n) = u(1, n - 1) - deltaT/deltaX*cplus*(u(1, n - 1) - u(Imax, n - 1)) - deltaT/deltaX*cminus*(u(2, n - 1) - u(1, n - 1));
    
    
    c = (((-2)*u(Imax - 1, n - 1)*VMax)/RhoMax + VMax);   
    cplus = max(c, 0);
        
    c = (((-2)*u(1, n - 1)*VMax)/RhoMax + VMax);
    cminus = min(c, 0);

    u(Imax, n) = u(Imax, n - 1) - deltaT/deltaX*cplus*(u(Imax, n - 1) - u(Imax - 1, n - 1)) - deltaT/deltaX*cminus*(u(1, n - 1) - u(Imax, n - 1));
    for i = 2: Imax-1
        c = (((-2)*u(i -1, n - 1)*VMax)/RhoMax + VMax);
        cplus = max(c, 0);
        
        c = (((-2)*u(i + 1, n - 1)*VMax)/RhoMax + VMax);
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

%Achsenbeschriftungen
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Vekehrsdichte_{A/m}','FontAngle','italic');

% Record the movie
for j = 1:Nmax
    plot(x, u(1:Imax,j), 'g','LineWidth',2);
    F(j) = getframe;
end
movie(F,0) % play x+1 times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Numerical vs. exact solution','NumberTitle','off')
hold on;
%Achsenbeschriftungen
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Vekehrsdichte_{A/m}','FontAngle','italic');

plot(x, u(1:Imax,1), x, u(1:Imax,  ceil(Nmax/6)),'-.r', x, u(1:Imax, 2* ceil(Nmax/6)), '--m', x, u(1:Imax, 3* ceil(Nmax/6)), ':b', x, u(1:Imax, 4* ceil(Nmax/6)), '-.g', x, u(1:Imax, 5* ceil(Nmax/6)), ':r', x, u(1:Imax, Nmax), '-.m');
%plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
%lgd = legend('Hallo','Tim','und','Phillip','wie','gehts','euch?');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','u','NumberTitle','off')
C = contourf( u');
CB = colorbar;
CB.Label.String = 'Verkehrsdichte_{A/m}';
%Achsenbeschriftungen
xlabel('x')
ylabel('t')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','numerical Solution u','NumberTitle','off')
mesh( u');
view (33, 12);
CB = colorbar;
CB.Label.String = 'Verkehrsdichte_{A/m}';
%Achsenbeschriftungen
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Zeit_{t}','FontAngle','italic');
zlabel('Verkehrsdichte_{A/m}','FontAngle','italic');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Geschwindigkeit','NumberTitle','off')
plot([RhoMax;0],[0;VMax])

%Achsenbeschriftungen
xlabel('Verkehrsdichte_{Rho}','FontAngle','italic');
ylabel('Maximalgeschwindigkeit_{m/s}','FontAngle','italic');
%plot(x, u(1:Imax,1), 'g','LineWidth',2);
%axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Zusätzlicher_Plot

figure('Name', 'Rho und f(Rho)','NumberTitle','off')
hold on;

%Achsenbeschriftungen
xlabel('Rho_{A/m}','FontAngle','italic');
ylabel('f(Rho)','FontAngle','italic');
p = 0:1/100:RhoMax;
fRho = -(VMax/RhoMax).*p.*p+VMax.*p;
plot(p, fRho);
lgd = legend('Verkehrsfluss');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Numerical vs. exact solution','NumberTitle','off')
hold on;
%Achsenbeschriftungen
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Vekehrsdichte_{A/m}','FontAngle','italic');

plot(x, u(1:Imax,1), x, u(1:Imax,  ceil(Nmax/6)),'-.r', x, u(1:Imax, 2* ceil(Nmax/6)), '--m', x, u(1:Imax, 3* ceil(Nmax/6)), ':b', x, u(1:Imax, 4* ceil(Nmax/6)), '-.g', x, u(1:Imax, 5* ceil(Nmax/6)), ':r', x, u(1:Imax, Nmax), '-.m');
%plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
%lgd = legend('Hallo','Tim','und','Phillip','wie','gehts','euch?');