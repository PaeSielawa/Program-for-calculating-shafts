clear all
close all
clc

%% Dane wału
a = 0.32; % w metrach
b = 0.18;
c = 0.075;
d = 0.005;
L = a+b+c+d;
r1 = 0.12; 
r2 = 0.3;
F1y = 10479.17; % w Netwonach
F1z = 1200;
F2z = 1100;
alfa = pi;
mi = 0.5;
F2z2 = F2z*exp(alfa*mi);
x = 0:0.001:(a+b+c+d);
xp = 0:0.001:(a+b+c+d);
Qb = 3.9*9.81; % N

%% Siły reakcji w płaszczyźnie XZ
%Rax = 0;
Vcz = ((a+b+c)*(F2z+F2z2) + a*F1z)/(a+b);
Vaz = F1z + F2z + F2z2 - Vcz;

%% Siły reakcji w płaszczyźnie XY
%Rax = 0;
Vcy = a*F1y/(a+b);
Vay = F1y - Vcy;

%% Wyliczanie wartości całkowitej sił reakcji
Va=(Vaz^2+Vay^2)^(1/2);
Vc=(Vcz^2+Vcy^2)^(1/2);

%% Wykresy sił tnących XZ
% Re1=200e6;

for i=1:1:581
    if i<=321
            Tz(i)=-Vaz;
        else if i >321
                if i<501
                Tz(i)=-Vaz+F1z;
                else 
                Tz(i)=-Vaz+F1z-Vcz;
                end
        end
    end
end

figure(1)
plot(x,Tz)
title('Siła tnąca w płaszczyźnie XZ')
xlabel('Odległość od początku wału [m]')
ylabel('Siła [N]')
xlim([0 0.581])
grid on

%% Wykres sił tnących XY
for i=1:1:581
    if i<=321
            Ty(i)=-Vay;
        else if i >321
                if i<501
                Ty(i)=-Vay+F1y;
                else 
                Ty(i)=-Vay+F1y-Vcy;
                end
        end
    end
end

figure(2)
plot(x,Ty)
title('Siła tnąca w płaszczyźnie XY')
xlabel('Odległość od początku wału [m]')
ylabel('Siła [N]')
xlim([0 0.581])
grid on

%% Wykres momentów zginających XY

for i=1:1:581
    if i<=321
            mgy(i)=x(i)*Vaz;
        else if i >321
                if i<501
                mgy(i)=x(i)*Vaz - (x(i)-a)*F1z;
                else 
                mgy(i)=x(i)*Vaz - (x(i)-a)*F1z + (x(i)-(a+b))*Vcz;
                end
        end
    end
end

figure(3)
plot(x,mgy)
title('Moment gnący z')
xlabel('Odległość od początku wału [m]')
ylabel('Moment zginający [Nm]')
xlim([0 0.575])
grid on

%% Wykres momentów zginających XZ
for i=1:1:581
    if i<=321
            mgz(i)=x(i)*Vay;
        else if i >321
                if i<501
                mgz(i)=x(i)*Vay - (x(i)-a)*F1y;
                else 
                mgz(i)=x(i)*Vay - (x(i)-a)*F1y + (x(i)-(a+b))*Vcy;
                end
        end
    end
end

figure(4)
plot(x,mgz)
title('Moment gnący y')
xlabel('Odległość od początku wału [m]')
ylabel('Moment zginający [Nm]')
xlim([0 0.581])
grid on

%% Wykres wypadkowego momentu gnącego
mg=sqrt(mgz.^2+mgy.^2);

figure(5)
plot(x,mg)
title('Wypadkowy moment gnący')
xlabel('Odległość od początku wału [m]')
ylabel('Moment zginający [Nm]')
xlim([0 0.575])
grid on

%% Wykres momentu skręcającego
for i=1:1:581
    if i<=321
        ms(i)=0;
    else if i >321
            ms(i) = F1y*r1;
         end
    end
end

figure(6)
plot(x,ms)
title('Moment skręcający')
xlabel('Odległość od początku wału [m]')
ylabel('Moment skręcający [Nm]')
xlim([0 0.581])
grid on

%% Wykres momentu zredukowanego (HMH)
for i=1:1:581             
     mzr(i)=sqrt(mg(i)^2 + 3/4*(ms(i)^2));            
end

figure(7)
plot(x,mzr)
title('Moment zredukowany')
xlabel('Odległość od początku wału [m]')
ylabel('Moment zredukowany [Nm]')
xlim([0 0.581])
grid on

%% Obliczanie średnicy wałów
% Rzeczywista średnica wału
E = 210e9;
Re = 235e6;
Rm = 340e6;
kg = 0.6*Re;

for i=1:1:581
    if i<=321
            mzrAB(i)=mzr(i);
        else if i >321
                if i<501
                mzrBC(i)=mzr(i);
                else 
                mzrCD(i)=mzr(i);
                end
            end
    end
end

for i=1:581
    D(i)=(32*mzr(i)/(pi*kg))^(1/3);
end

figure(8)
plot(x,D)
xlim([0 0.581])
title('Rzeczywista min średnica wału')
xlabel('Odległość od początku wału [m]')
ylabel('Średnica wału [m]')
grid on

% Minimalna średnica wału
dAB=(32*max(mzrAB)/(pi*kg))^(1/3); % min d na odcinku AB
dBC=(32*max(mzrBC)/(pi*kg))^(1/3); % min d na odcinku BC
dCD=(32*max(mzrCD)/(pi*kg))^(1/3); % min d na odcinku CD

disp('min d [m] na odcinku AB to') 
disp(dAB)
disp('min d [m] na odcinku BC to')
disp(dBC)
disp('min d [m] na odcinku CD to')
disp(dCD)

for i=1:1:581
    if i<=321
            d(i)=dAB;
        else if i >321
                if i<501
                d(i)=dBC;
                else 
                d(i)=dCD;
                end
            end
    end
end

figure(9)
plot(x,d);
title('Minimalna Średnica wału')
xlabel('Odległość od początku wału [m]')
ylabel('Średnica wału [m]')
xlim([0 0.581])
ylim([0 0.581/10])
grid on

disp('Dobrane d [m] na odcinku AA1 to') % na łożysko
dAA1d=0.05;
disp(dAA1d)
disp('Dobrane d [m] na odcinku A1B1 to')
dA1B1d=0.052;
disp(dA1B1d)
disp('Dobrane d [m] na odcinku B1B to')
dB1B2d=0.06;
disp(dB1B2d)
disp('Dobrane d [m] na odcinku BC to')
dB2Cd=0.055;
disp(dB2Cd)
disp('Dobrane d [m] na odcinku CD to') % na koło pasowe
dCDd=0.05;
disp(dCDd)

for i=1:34
d(i)=dAA1d;
end
for i=35:249
d(i)=dA1B1d;
end
for i=250:483
d(i)=dB1B2d;
end
for i=484:527
d(i)=dB2Cd;
end
for i=528:581
d(i)=dCDd;
end

figure(10)
plot(x,d);
title('Dobrane średnice wału')
xlabel('Odległość od początku wału [m]')
ylabel('Średnica wału [m]')
xlim([0 0.581])
ylim([0 0.07])
grid on

%% Strzałka ugięcia metodą Clebsh'a 
IAA1=(pi*(dAA1d/2)^4)/4;
IA1B1=(pi*(dA1B1d/2)^4)/4;
IB1B2=(pi*(dB1B2d/2)^4)/4;
IB2C=(pi*(dB2Cd/2)^4)/4;
ICD=(pi*(dCDd/2)^4)/4;

% XZ
D1 = 0; %wA(0) = 0
C1 = -Vaz*(0.5)^3/3 + (F1z+Qb)*(0.5-a)^3/3 - Vcz*(0.5-a-b)^3/3; %wC(0.5) = 0

wzB = -(D1 + C1*0.32 + Vaz*(0.32)^3/6)/(E*IB1B2); %wzB(0.32)
wzD = -(D1 + C1*0.575 + Vaz*(0.575)^3/6 - (F1z+Qb)*(0.575-a)^3/6 + Vcz*(0.575-a-b)^3/6)/(E*IB1B2); %wcD(0.575)

% XY
D2 = 0; 
C2 = -Vay*(0.5)^3/3 + F1y*(0.5-a)^3/3 - Vcy*(0.5-a-b)^3/3; %wC(0.5) = 0

wyB = -(D2 + C2*0.32 + Vay*(0.32)^3/6)/(E*IB1B2); %wzB(0.32)
wyD = -(D2 + C2*0.575 + Vay*(0.575)^3/6 - F1y*(0.575-a)^3/6 + Vcy*(0.575-a-b)^3/6)/(E*IB1B2); %wcD(0.575)

% Wypadkowe strzałki ugięcia
wB = sqrt(wyB^2 + wzB^2);
wD = sqrt(wyD^2 + wzD^2);

%Kąt ugięcia w punktach a i c

%XZ
wz1A = -C1/(E*IB2C); %kąt ugięcia wz'(0)
wz1C = -(C1 + Vaz*(0.5)^2/2 - F1z*(0.5-a)^2/2 + Vcz*(0.5-a-b)^2/2)/(E*IB2C); %kąt ugięcia wz'(0.5)

%XY
wy1A = -C2/(E*IB2C); %kąt ugięcia wy'(0)
wy1C = -(C2 + Vay*(0.5)^2/2 - F1y*(0.5-a)^2/2 + Vcy*(0.5-a-b)^2/2)/(E*IB2C); %kąt ugięcia wy'(0.5)

W1A = sqrt(wy1A^2 + wz1A^2);
W1C = sqrt(wy1C^2 + wz1C^2);

% Ugięcie dopuszczalne
fdop = 0.0004*L;

%% Kąt skręcenia
% G=80000000000; %Pa
% 
% for i=1:40
% phi(i) = ms(i)*(x(i)-0.8)/IAA1/G;
% end
% for i=41:250
% phi(i) = ms(i)*(x(i)-0.8)/IA1B1/G;
% end
% for i=251:480
% phi(i) = ms(i)*(x(i)-0.8)/IB1B2/G;
% end
% for i=481:527
% phi(i) = ms(i)*(x(i)-0.8)/IB2C/G;
% end
% for i=528:581
% phi(i) = ms(i)*(x(i)-0.8)/ICD/G;
% end
% 
% figure(11)
% plot(x,phi)
% title('Kąt skręcenia wału')
% xlabel('Długość wału [m]')
% ylabel('Kąt skręcenia [rad]')
% xlim([0 0.581])
% grid on

%% Dobór łożysk (kulkowe)
Pw = 0; % obciążenie wzdłużne
V = 1; % dla obracających się wałów = 1
X = 1; % dla poprzecznych = 1
Y = 0;
n = 500; % obr/min
epsilon = 3; % dla łożysk kulkowych
wsp_obc = 1.2; % współczynnik obciążenia
L10h = 10000; % oczekiwany czas pracy

L10=(60*n*L10h)/(10^6);

PzA = (X*V*Va + Y*Pw)*wsp_obc; % obciążenie zastępcze w punkcie A
PzC = (X*V*Vc + Y*Pw)*wsp_obc; % obciążenie zastępcze w punkcie C

Ca = PzA*(L10)^(1/epsilon); % nośność dynamiczna łożyska w punkcie A
% Dobieramy łożysko: 62210-2RS1  C = 35.1 kN

Cc = PzC*(L10)^(1/epsilon); % nośność dynamiczna łożyska w punkcie C
% Dobieramy łożysko: 6411  C = 99.5 kN

%% Dobieranie wpustów
Rew = 900; % MPa
ktw = 500;
kdw = 360;
Ms = 1257.5; % Nm

szerB=0.019;
szerD=0.016;

wysB=0.011;
wysD=0.010;

ldopB=2*Ms/(szerB*dB1B2d*ktw*10^(6)); % warunek na scinanie dla B
ldopB2=4*Ms/(wysB*dB1B2d*kdw*10^(6)); % warunek na doscik dla B

ldopD=2*Ms/(szerD*dCDd*ktw*10^(6)); %warunek na scinanie dla D
ldopD2=4*Ms/(wysD*dCDd*kdw*10^(6)); %warunek an doscik dla D

%% Obliczanie prędkości krytycznej

g = 9.81;
omegaZB = sqrt(g/abs(wzB));
omegaZD = sqrt(g/abs(wzD));
omegaZo = 1/omegaZB^2 + 1/omegaZD^2;
omegaZ = sqrt(1/omegaZo);

omegaYB = sqrt(g/abs(wyB));
omegaYD = sqrt(g/abs(wyD));
omegaYo = 1/omegaYB^2 + 1/omegaYD^2;
omegaY = sqrt(1/omegaYo);

omegaC = sqrt(omegaZ^2 + omegaY^2);
nkr = 30/pi * omegaC;

%% Obliczenia zmęczeniowe
% zginanie
WG = (pi*0.06^4)/32;
Sigma = 1218.9/WG;
yps = 1.5; 
Zgo = 250000000; % MPa
beta = 2.06;
delta = Zgo/(Sigma*beta*yps);

% skręcanie
WG2 = (pi*0.06^3)/16;
Zso = 90000000;
beta2 = 2;
yps2 = 1.54;
Tau = 1257.5/(4*WG2);

delta2 = Zso/(Tau*beta2*yps2);

Delta = delta*delta2/(sqrt(delta^2 + delta^2));
