clc;
clear;
figure(1)
a = importdata('clusterUBV.txt');
id = a.data(:,1);
source = a.data(:,2);
V = a.data(:,3);
BV = a.data(:,4);
UB = a.data(:,5);
N = a.data(:,6);
plot(BV,V,'.r')
set(gca,'Ydir','reverse')

figure(2)
plot(BV,UB,'.r')
set(gca,'Ydir','reverse')
%%
b = importdata('UBV_intrinsic_ms.txt');
Mv = b.data(:,1);
BVo = b.data(:,2);
figure(3)
plot(BVo,Mv)
set(gca,'Xdir','reverse')
