clear all

b1 = 0.5E-3
b2 = 1E-3
b3 = 2E-3
bmin = 0.14E-3
bminfull = 0.22E-3

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
plot(xplot, y1)
plot(xplot, y2)
plot(xplot, y3)

    