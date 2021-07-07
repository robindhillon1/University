clear all

b1 = 0.5E-3
b2 = 1E-3
b3 = 2E-3
bmin = 0.14E-3      %Becomes singular here. 
bminfull = 0.20E-3

xplot = linspace(-1E-3, 3E-3, 5000)

y = @(x, b)abs((sin(10000*pi*x)/(10000*pi*x) + sin(10000*pi*(x-b))/(10000*pi*(x-b))));

for i=1:5000
    y1(1,i) = y(xplot(i), b1);
    y2(1,i) = y(xplot(i), b2);
    y3(1,i) = y(xplot(i), b3);
    ymin(1,i) = y(xplot(i), bmin);
    yminfull(1,i) = y(xplot(i), bminfull);
end

figure()
hold on
plot(xplot*1000, y1)
plot(xplot*1000, y2)
plot(xplot*1000, y3)
xlabel('x (mm)')
ylabel('|g(x,0)|')
title('Magnitude of image as function of x(mm)')
legend('b = 0.5mm', 'b = 1mm', 'b = 2mm')
grid on

figure(2)
hold on
%plot(xplot*1000, ymin)
plot(xplot*1000, yminfull)
xlabel('x (mm)')
ylabel('|g(x,0)|')
title('Image has 2 discernible peaks when b=0.2mm')
grid on
    