clear all
w0 = 632.8E-9;
n0 = 120*3.14159;
ai = 1/n0;
as = 1.52/n0;
a1 = 1.65/n0;
a2 = 2.1/n0;

delta = @(x)3.1414*w0/(2*x);
F1 = @(delta)[cos(delta) -i*sin(delta)/a1 ; -i*a1*sin(delta) cos(delta)];
F2 = @(delta)[cos(delta) -i*sin(delta)/a2 ; -i*a2*sin(delta) cos(delta)];

F = @(F1,F2)(F1*F2);

FN = @(F)F^19;
R = @(A, B, C, D)abs((A*ai+ai*as*B-C-as*D)/(A*ai+ai*as*B+C+as*D));

x = linspace(580E-9, 696E-9, 100);

for xx = 1:100
    d = delta(x(xx));
    FF = FN(F(F1(d), F2(d)));
    A = FF(1, 1);
    B = FF(1, 2);
    C = FF(2, 1);
    D = FF(2, 2);
    
    R_vec(1,xx) = R(A, B, C, D);
end
figure
plot(x./w0, 100.*R_vec)
xlim([0.9 1.12])
%ylim([99.8 100])
xlabel('Normalized Wavelength')
ylabel('Reflectance(100%)')
title('Reflectivity at N=19 as function of wavelength')
grid on