clear all

A = load('temperature_1880_2016.dat');
temp = gettemp(A);

fid = fopen('monthly_in_situ_co2_mlo.csv');
co2 = getco2(fid);
fclose(fid);

subplot(2,1,1)
plot(temp.time,temp.data,'r.')
text(7*10^5,1,'~Temp Change of 1 Degree')
title('NCEI Temp 1880-Present')
xlabel('Time (yr)')
ylabel('Temperature (^oC)')
datetick('x')

subplot(2,1,2)
plot(co2.time,co2.data,'g.')
text(7.2*10^5,400,'Co2 Change of 100ppm')
title('~CO^2 Concentration 1958-Present')
xlabel('Time (yr)')
ylabel('CO^2 (ppm)')
datetick('x')

% Using Column 5 because it is the unadjusted data


figure(2)

stats = getstats(co2.data,7); % mean and SD
t = co2.time;
X = co2.month==1; % logical index
Y = stats.mean(X); % logical indexing
Y = [Y(2:end)-Y(1:58) nan]; % annual difference

subplot(3,1,1)
plot(t,co2.data,'g.')
title('CO^2 Concentration 1958-Present')
xlabel('Time (yr)')
ylabel('CO^2 (ppm)')
datetick('x')
ylim([300 450])
hold on
plot(t,stats.mean)

subplot(3,1,2)
plot(t,stats.sd)
title('Standard Deviation 1958-Present')
ylabel('SD')
xlabel('Time (yr)')
datetick('x')

subplot(3,1,3)
plot(t(X),Y,'r')
title('Annual Change in CO^2 1958-Present')
ylabel('Change (ppm)')
xlabel('Time (yr)')
datetick('x')

num1 = nansum(stats.sd)/length(stats.sd);
disp('Average Annual Sandard Deviation is 1.2981')
num2 = num2str(nansum(Y)/length(Y));
disp('Average Annual Change is 1.4496 Since 1958')

% Pt 6.6
% a)
% Yes, co2 is increasing every year b/c the change is never a negative
% number.
% b)
% I think the time observed so far is enough data to assume that it is not
% seasonal changes is the levels of co2. Weather patterns usually show
% natural swings thoughout the course of a year and not between multiple
% years. Steady increase across multiple years is more likely to be a sign
% of a overall co2 production than a seasonal process. It is visible in
% fig.2 plot 3 that the increase is also accelerating which could not be
% due to anything but Co2 producion

figure(3)

stats = getstats(temp.data,7); % mean and SD
t = temp.time;
X = temp.month==1; % logical index
Y = stats.mean(X); % logical indexing
Y = [Y(2:end)-Y(1:136) nan];

subplot(3,1,1)
plot(t,temp.data,'r.')
title('NCEI Temp 1880-Present')
xlabel('Time (yr)')
ylabel('Temperature (^oC)')
datetick('x')
hold on
plot(t,stats.mean,'g')

subplot(3,1,2)
plot(t,stats.sd)
title('Standard Deviation 1880-Present')
ylabel('SD')
xlabel('Time (yr)')
datetick('x')

subplot(3,1,3)
plot(t(X),Y,'r')
title('Annual Change in Temp 1880-Present')
ylabel('Change (^oC)')
xlabel('Time (yr)')
datetick('x')

figure(4)

stats = getstats(temp.data,21); % mean and SD
t = temp.time;
X = temp.month==1; % logical index
Y = stats.mean(X); % logical indexing
Y = [Y(2:end)-Y(1:136) nan];

subplot(3,1,1)
plot(t,temp.data,'r.')
title('NCEI Temp 1880-Present')
xlabel('Time (yr)')
ylabel('Temperature (^oC)')
datetick('x')
hold on
plot(t,stats.mean,'g')

subplot(3,1,2)
plot(t,stats.sd)
title('Standard Deviation 1880-Present')
ylabel('SD')
xlabel('Time (yr)')
datetick('x')

subplot(3,1,3)
plot(t(X),Y,'r')
title('Annual Change in Temp 1880-Present')
ylabel('Change (^oC)')
xlabel('Time (yr)')
datetick('x')
