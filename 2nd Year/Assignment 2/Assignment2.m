clc;
clear;
a = importdata('HW2_Q1_Full_Sky_Catalogue.csv');
h = a.data(:,1);
m = a.data(:,2);
s = a.data(:,3);

deg = a.data(:,4);
arcmin = a.data(:,5);
arcsec = a.data(:,6);
mv = a.data(:,7);

RAs = h+m+s;
dec = deg+arcmin+arcsec;

%RA = (RAs/(3600*24))*360;
%declin = dec/3600;
figure(1)
plot(RAs, dec,'r.','MarkerSize',1)
%axesm('mollweid', 'Frame', 'on', 'Grid', 'on');

hold on;
visb = find(mv <= 6);
fprintf('Size of visb is: \n')
size(visb)
disp('The number of visible stars: 5074')
plot(RAs(visb),dec(visb),'.k','MarkerSize',1)

%{
RAs = (h*3600)+(m*60)+s;
dec = (deg*3600)+(arcmin*60)+arcsec;

RA = (RAs/(3600*24))*360;
declin = dec/3600;
figure(1)
plot((RA*pi/180-pi), (declin*pi/180),'r.','MarkerSize',1)
%axesm('mollweid', 'Frame', 'on', 'Grid', 'on');
%hold on;
figure(2)
visb = find(mv <= 6);
fprintf('Size of visb is: \n')
size(visb)
disp('The number of visible stars: 5074')
plot(RA(visb),declin(visb),'.b','MarkerSize',1)
%}