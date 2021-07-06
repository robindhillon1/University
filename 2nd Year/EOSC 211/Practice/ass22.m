%% Kyle Zell 50209121 -----------------------------------------------------------------------
clear all
clf
clc

% ------------------------------------------------
A = load('temperature_1880_2016.dat');
temp = gettemp(A);

fid = fopen('monthly_in_situ_co2_mlo.csv');
co2 = getco2(fid);
fclose(fid);


%% Figure 1 -------------------------------------------------

figure(1)
subplot(2,1,1)
plot(temp.time,temp.data,'r.')
text(7*10^5,1.5,'Average Temp Change of 1^oC/yr')
title('NCEI Temp 1880-Present')
xlabel('Time/yr')
ylabel('Temperature/^oC')
datetick('x')

subplot(2,1,2)
plot(co2.time,co2.data,'g.')
text(7.2*10^5,425,'Average Change of CO_2 ~100ppm/yr')
title('CO^2 Record 1958-Present')
xlabel('Time/yr')
ylabel('CO_2/ppm')
datetick('x')

% I chose the fifth colmn because that is the location of all of the
% meaningful data.  The column i chose has all of the raw monthly CO2 data

%% Figure 2 ------------------------------------------------

stats = getstats(co2.data,7); % mean and SD
T = co2.time;
X = co2.month==1; % logical index
Y = stats.mean(X); % logical indexing
Y = [Y(2:end)-Y(1:58) nan]; % annual difference


figure(2)
subplot(3,1,1)
plot(T,co2.data,'g.')
title('CO^2 Concentration 1958-Present')
xlabel('Time/yr')
ylabel('CO^2/ppm')
datetick('x')
ylim([300 450])
hold on
plot(T,stats.mean)

subplot(3,1,2)
plot(co2.time,stats.sd)
title('Standard Deviation 1958-Present')
ylabel('SD')
xlabel('Time (yr)')
datetick('x')

subplot(3,1,3)
plot(co2.time(X),Y,'r')
title('Annual Change in CO^2 1958-Present')
ylabel('Change (ppm)')
xlabel('Time (yr)')
datetick('x')


num1 = nansum(stats.sd)/length(stats.sd);
disp('Average Annual Standard Deviation is 1.2981')
num2 = num2str(nansum(Y)/length(Y));
disp('Average Annual Change is 1.4496')

% Average Annual Standard Deviation is 1.2981'
% Average Annual Change is 1.4496

% 6A. According to the figure the CO2 levels are not in a constant increase.
%  The general trend is a gradual increase in CO2 since the industrial
%  revolution however the levels are oscillating between increasing and
%  decreasing due to seasonal variability.

% 6B. The amount of years to detect a sufficiently large increase in CO2
% would be according to the level of the sufficiently large increase.  If a
% large increase was adequately calculated through the world it would take
% 2 years(2 cycles through each season) to judge whether the CO2 levels
% were increased due to the small standard deviation in the plot.

%% Figure 3 ----------------------------------------------------
figure(3)

stats = getstats(temp.data,7); % mean and SD
T = temp.time;
X = temp.month==1; % logical index
Y = stats.mean(X); % logical indexing
Y = [Y(2:end)-Y(1:136) nan];

subplot(3,1,1)
plot(T,temp.data,'g.')
title('NCEI Temp 1880-Present')
xlabel('Time/yr')
ylabel('Temperature/^oC')
datetick('x')
hold on
plot(T,stats.mean,'g')

subplot(3,1,2)
plot(T,stats.sd)
title('Standard Deviation 1880-Present')
ylabel('SD')
xlabel('Time/yr')
datetick('x')

subplot(3,1,3)
plot(T(X),Y,'r')
title('Annual Change in Temp 1880-Present')
ylabel('Change (^oC)')
xlabel('Time (yr)')
datetick('x')


% The problem with using annual temperature chages is the annual
% temperature change has so much variability and that the seasonality of
% the regions can change so much year to year.  Especially when you account
% for changes like el nino and la nina.

%% Figure 4 -------------------------------------------------
figure(4)

stats = getstats(temp.data,21); 
T = temp.time;
X = temp.month==1; 
Y = stats.mean(X); 
Y = [Y(2:end)-Y(1:136) nan];

subplot(3,1,1)
plot(T,temp.data,'r.')
title('NCEI Temp 1880-Present')
xlabel('Time/yr')
ylabel('Temperature/^oC')
datetick('x')
hold on
plot(T,stats.mean,'g')

subplot(3,1,2)
plot(T,stats.sd)
title('Standard Deviation 1880-Present')
ylabel('SD')
xlabel('Time/yr')
datetick('x')

subplot(3,1,3)
plot(T(X),Y,'r')
title('Annual Change in Temp 1880-Present')
ylabel('Change of Temp/^oC')
xlabel('Time/yr')
datetick('x')

% The overall temperature rise due to anthropogenic climate change can be
% best explained by the plots in figure 1.  The main effects of climate
% change can be shown in the plot from 1958 to present.  As levels of CO2
% rise from 1958 to present we see a steady increase in the temperature.
% From 1958 to present there is an increase in CO2 concentrations that
% parallel the increase in temperature.  


partner.name='Kinnon'; % Last name
partner.collab='Worked on the whole thing together';
Time_spent= 6; %hrs

