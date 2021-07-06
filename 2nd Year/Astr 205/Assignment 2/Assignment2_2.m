clc;
clear;
figure(1)
a = importdata('HW2_Q2_data.txt');
vmag = a.data(:,1);
imag = a.data(:,2);
dx = a.data(:,3);
dy = a.data(:,4);
plot(vmag,(vmag-imag),'.r')
set(gca,'Ydir','reverse')
%HR Diagram
figure(2)
plot(dx,dy,'.r')
%Majority of the stars belong to the cluster. The CMD of the proper motions
%shows that the majority lie in between x = (-4:4) and y =
%(-4:4), approximately. Using these x and y values for the proper
%motion, we can estimate and isolate the 47 Tucanae cluster. 
figure(3)
rangx = find(dx>=-4 & dx<= 4);
rangy = find(dy>=-4 & dy<= 4);
plot(vmag(rangy),imag(rangx),'.r')
set(gca,'Ydir','reverse')