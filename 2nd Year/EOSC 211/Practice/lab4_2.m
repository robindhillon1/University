
%constants
Tc1 = 25;
Tc2 = [-30:40]

Tk1 = Tc1+273.15;
Tk2= Tc2+273.15;
Rd = 287;       %J/(K*kg)
Lvap = 2.5*10^6; %J/kg
Rv = 461;      %J/(K*kg)
T0 = 273;      %K
e0 = 611;      %Pa

%Changeable constants
RH = 50;        % '% humidity'
RH2 = 100;      % '% humidity'
P = 102000 ;   %Pa
P2 = 90000 ;   %Pa

format long

%question 1
es = e0*exp((Lvap/Rv)*(1/T0-1/Tk1));        % es is in pascals

%question 2
e = es*(RH/100);
ee = Rd/Rv;
rr = (ee*e)/(P-e);

rho = P/(Rd*Tk1*(1+0.61*rr));     %kg/m^3

%question 3

esat = e0*exp((Lvap/Rv)*(1/T0-1./Tk2));

subplot(4,1,1);
plot(Tc2,esat)
xlabel('Temperature in Celsius')
ylabel('Saturation Vapor Pressure')

%question 4

e2 = esat.*(RH2/100);
rr2 = (ee.*e2)./(P2-e2);

subplot(4,1,2);
plot(Tc2,rr2)
xlabel('Temperature in Celsius')
ylabel('Mixing Ratio')

%question 5

[RH3,Tc3]=meshgrid([0:20:100],[-30:40])

Tk3= Tc3+273.15;
esat2 = e0*exp((Lvap/Rv)*(1/T0-1./Tk3));

e3 = esat2.*RH3./100;
mixrat = (ee.*e3)./(P2-e3);

subplot(4,1,3);
plot(Tc3,mixrat)
xlabel('Temperature (^oC)')
ylabel('r (g/g)')


%question 6

rho2 = P2./(Rd.*Tk3.*(1+0.61.*mixrat));

subplot(4,1,4);
plot(Tc3,rho2)
xlabel('Temperature (^oC)')
ylabel('rho kg/m^3')


