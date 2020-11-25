clear all;
close all;

% löse die Transportgleichung: u_t + u_u_x = 0
L = 230;            % Streckenlänge (Meter)
T = 60;             % Time of interest (Sekunden) 
CFL = 0.9;
Imax = 500;         % Anzahl Stützstellen
VMax = 11.1;        % Maximale Geschwindigkeit (bei RhoMin) (m/s) 
RhoMax = 0.17;      % Maximale Anzahl an Autos innerhalb eines Kilometers   Auto / Meter
RhoStart = 0.095;   % Minimale Anzahl an Autos innerhalb eines Kilometers (Maximaler Wert bei Einhaltung des Mindestabstands) Auto/Meter

x = linspace(0, L, Imax);
deltaX = L/(Imax - 1);

u0 = zeros(Imax, 1);

% set initial conditions
% Sinus so eingestellt dass er eine Schwingung im bereich (120, 240) macht
% (Freuquenz)
AnfangDichtestoerung = 80;
EndeDichtestoerung = 200;
freq = 2*pi/(EndeDichtestoerung-AnfangDichtestoerung);

for i = 1: Imax
    if(x(i)<AnfangDichtestoerung)
        u0(i,1) = RhoStart;
    elseif(x(i)<140)
        u0(i,1)= (1 + 0.3*sin(freq*(x(i) - 80)))* RhoStart;
    else
        u0(i,1)= RhoStart;
    end
end

%maxU0 = -inf;

%for i =1: Imax
%    maxU0 = max(u0(i, 1), maxU0);
%end

%c = -(((-2)*maxU0)*VMax/RhoMax + VMax);
a_Rho = (((-2)*RhoMax)*VMax/RhoMax + VMax);
deltaT = CFL*deltaX/abs(a_Rho);
Nmax = ceil(T/deltaT) + 1;

u = zeros(Imax, Nmax);
t = linspace(0, T, Nmax);
v = zeros(Imax, Nmax);
wg = zeros(Imax, Nmax); %wellengeschwindigkeiten
pos = zeros(Imax, Nmax);

for i = 1: Imax
   u(i,1) = u0(i,1);
   v(i, 1) = -(VMax/RhoMax)*u0(i, 1) + VMax;
   wg(i, 1) = (((-2)*u0(i, 1)*VMax)/RhoMax + VMax);
end



for n = 2: Nmax
    %backward differences
    %periodic boundary conditions
    a_Rho = (((-2)*u(Imax, n - 1)*VMax)/RhoMax + VMax);
    a_RhoPlus = max(a_Rho, 0);
    
    a_Rho = (((-2)*u(2, n - 1)*VMax)/RhoMax + VMax);
    a_RhoMinus = min(a_Rho, 0);
    
    u(1, n) = u(1, n - 1) - deltaT/deltaX*a_RhoPlus*(u(1, n - 1) - u(Imax, n - 1)) - deltaT/deltaX*a_RhoMinus*(u(2, n - 1) - u(1, n - 1));
    v(1, n) = -(VMax/RhoMax)*u(1, n) + VMax;
    if (a_RhoPlus ~= 0)
            wg(1, n) = a_RhoPlus;
        elseif (a_RhoMinus ~= 0)
            wg(1, n) = a_RhoMinus;
        else
            wg(1, n) = 0;
    end
    pos(1, n) = (a_Rho/2)*(deltaT^2) + v(1, n)*deltaT + pos(1, n-1);
    
    
    a_Rho = (((-2)*u(Imax - 1, n - 1)*VMax)/RhoMax + VMax);   
    a_RhoPlus = max(a_Rho, 0);
        
    a_Rho = (((-2)*u(1, n - 1)*VMax)/RhoMax + VMax);
    a_RhoMinus = min(a_Rho, 0);
    
    u(Imax, n) = u(Imax, n - 1) - deltaT/deltaX*a_RhoPlus*(u(Imax, n - 1) - u(Imax - 1, n - 1)) - deltaT/deltaX*a_RhoMinus*(u(1, n - 1) - u(Imax, n - 1));
    v(Imax, n) = -(VMax/RhoMax)*u(Imax, n) + VMax;
    if (a_RhoPlus ~= 0)
            wg(Imax, n) = a_RhoPlus;
        elseif (a_RhoMinus ~= 0)
            wg(Imax, n) = a_RhoMinus;
        else
            wg(Imax, n) = 0;
    end
    pos(Imax, n) = (a_Rho/2)*(deltaT^2) + v(Imax, n)*deltaT + pos(Imax, n-1);
    
    for i = 2: Imax-1
        a_Rho = (((-2)*u(i - 1, n - 1)*VMax)/RhoMax + VMax);
        a_RhoPlus = max(a_Rho, 0);
        
        a_Rho = (((-2)*u(i + 1, n - 1)*VMax)/RhoMax + VMax);
        a_RhoMinus = min(a_Rho, 0);
        
        u(i, n) = u(i, n - 1) - deltaT/deltaX*a_RhoPlus*(u(i, n - 1) - u(i - 1, n - 1)) - deltaT/deltaX*a_RhoMinus*(u(i + 1, n - 1) - u(i, n - 1));
        u(i, n) = min(u(i, n), RhoMax);
        v(i, n) = -(VMax/RhoMax)*u(i, n) + VMax;
        pos(i, n) = (a_Rho/2)*(deltaT^2) + v(i, n)*deltaT + pos(i-1, n-1);
        if (a_RhoPlus ~= 0)
            wg(i, n) = a_RhoPlus;
        elseif (a_RhoMinus ~= 0)
            wg(i, n) = a_RhoMinus;
        else
            wg(i, n) = 0;
        end
        
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
ylabel('Vekehrsdichte_{\frac{A}{m}}','FontAngle','italic');

%plot(x, u);
plot(x, u(1:Imax,1), x, u(1:Imax,  ceil(Nmax/6)),'-.r', x, u(1:Imax, 2* ceil(Nmax/6)), '--m', x, u(1:Imax, 3* ceil(Nmax/6)), ':b', x, u(1:Imax, 4* ceil(Nmax/6)), '-.g', x, u(1:Imax, 5* ceil(Nmax/6)), ':r', x, u(1:Imax, Nmax), '-.m');
%plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
%lgd = legend('Hallo','\lfloor Nmax/6 \rfloor','und','Phillip','wie','gehts','euch?');
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

%figure('Name', 'Geschwindigkeit','NumberTitle','off')
%plot([RhoMax;0],[0;VMax])

%Achsenbeschriftungen
%xlabel('Verkehrsdichte_{Rho}','FontAngle','italic');
%ylabel('Maximalgeschwindigkeit_{m/s}','FontAngle','italic');
%plot(x, u(1:Imax,1), 'g','LineWidth',2);
%axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Zusätzlicher_Plot

%figure('Name', 'Rho und f(Rho)','NumberTitle','off')
%hold on;

%Achsenbeschriftungen
%xlabel('Rho_{A/m}','FontAngle','italic');
%ylabel('f(Rho)','FontAngle','italic');
%p = 0:1/100:RhoMax;
%fRho = -(VMax/RhoMax).*p.*p+VMax.*p;
%plot(p, fRho);
%lgd = legend('Verkehrsfluss');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Strömungsgeschwindigkeit der Welle','NumberTitle','off')
hold on;
%Achsenbeschriftungen
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Strömungsgeschwindigkeit_{m/s}','FontAngle','italic');

plot(x, wg(1:Imax, 1), x, wg(1:Imax,  ceil(Nmax/6)),'-.r', x, wg(1:Imax, 2* ceil(Nmax/6)), '--m', x, wg(1:Imax, 3* ceil(Nmax/6)), ':b', x, wg(1:Imax, 4* ceil(Nmax/6)), '-.g', x, wg(1:Imax, 5* ceil(Nmax/6)), ':r', x, wg(1:Imax, Nmax), '-.m');
%plot(x, wellengeschwindigkeiten(1:Imax,1), x, wellengeschwindigkeiten(1:Imax,  ceil(Nmax/6)),'-.r', x, wellengeschwindigkeiten(1:Imax, 2* ceil(Nmax/6)), '--m', x, wellengeschwindigkeiten(1:Imax, 3* ceil(Nmax/6)), ':b', x, wellengeschwindigkeiten(1:Imax, 4* ceil(Nmax/6)), '-.g', x, wellengeschwindigkeiten(1:Imax, 5* ceil(Nmax/6)), ':r', x, wellengeschwindigkeiten(1:Imax, Nmax), '-.m');
%plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
%lgd = legend('Hallo','Tim','und','Phillip','wie','gehts','euch?');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Geschwindigkeiten der Fahrzeuge','NumberTitle','off')
hold on;
%Achsenbeschriftungen
xlabel('Streckenlänge_{m}','FontAngle','italic');
ylabel('Geschwindigkeiten_{m/s}','FontAngle','italic');

%plot(x, v);
plot(x, v(1:Imax,1), x, v(1:Imax,  ceil(Nmax/6)),'-.r', x, v(1:Imax, 2* ceil(Nmax/6)), '--m', x, v(1:Imax, 3* ceil(Nmax/6)), ':b', x, v(1:Imax, 4* ceil(Nmax/6)), '-.g', x, v(1:Imax, 5* ceil(Nmax/6)), ':r', x, v(1:Imax, Nmax), '-.m');
%plot(x, v(ceil(Imax/6), 1), x, v(ceil(Imax/5), 1),'-.r', x, v(ceil(Imax/4), 1), '--m', x, v(ceil(Imax/3), 1), ':b', x, v(ceil(Imax/2), 1), '-.g', x, v(ceil(Imax), 1));
lgd = legend('Hallo','Tim','und','Phillip','wie','gehts','euch?');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'Position im Kreis','NumberTitle','off')
hold on;
%Achsenbeschriftungen
xlabel('Position_{m}','FontAngle','italic');
ylabel('Zeit_{sek}','FontAngle','italic');

plot(pos(1:Imax, 1:25:Nmax));