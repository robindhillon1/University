clear all
a1 = 1;
a2 = .1;
a3 = .5;
x=linspace(-1, 6, 5000);

E1 = fresnelc(a1*x)+i*fresnels(a1*x)+fresnelc(inf)+i*fresnels(inf);
E2 = fresnelc(a2*x)+i*fresnels(a2*x)+fresnelc(inf)+i*fresnels(inf);
E3 = fresnelc(a3*x)+i*fresnels(a3*x)+fresnelc(inf)+i*fresnels(inf);

I1 = abs(E1).^2/abs(max(E1(:)).^2);
I2 = abs(E2).^2/abs(max(E2(:)).^2);
I3 = abs(E3).^2/abs(max(E3(:)).^2);
%I4 = abs(E4).^2/abs(max(E4(:)).^2);

figure(1)
plot(x, I1)
figure(2)
plot(x, I2)
figure(3)
plot(x, I3)
