% Setup
clc;
clearvars;
%% Calculate theoretical Fresnel diffraction pattern
% Parameters
w = 8.5e-3; % Slit width, m
lambda = 632.8e-9; % Laser wavelength, m
R = 383e-2; % Slit-to-screen distance, m
% Calculate ?v
delta_v = w*(2/(lambda*R))^(1/2);
% Initialize plot ranges
x_range = 3e-2; % One-sided x plot range, m
x_res = 1000; % Number of points between the two x plot range limits
x = linspace(-x_range, x_range, x_res); % x plot points
z_range = x_range/w; % One-sided z plot range
z_res = x_res; % Number of points between the two z plot range limits
z = linspace(-z_range, z_range, z_res); % z plot points
% Calculate intensity pattern (Fraunhofer)
beta = (z.*pi/2).*(delta_v).^2;
I = ((delta_v).^2).*((sin(beta)).^2)./(beta.^2);
I_norm_fraun = I./max(I);
% Calculate intensity pattern (Fresnel)
v_1 = -(z+0.5).*delta_v;
v_2 = -(z-0.5).*delta_v;
C = fresnelc(v_2) - fresnelc(v_1); % Fresnel cosine integral
S = fresnels(v_2) - fresnels(v_1); % Fresnel sine integral
I = C.^2 + S.^2;
I_norm_fresnel = I./max(I);
% Plot I vs. z
% figure(1);
% plot(z, I_norm);
% xlabel('z = x/w');
% ylabel('Normalized Intensity');
% title(['Fresnel Diffraction Pattern, \Deltav = ', num2str(delta_v, 2), ' ']);
% Plot I vs. x
figure;
plot(x*1000, I_norm_fraun);
xlabel('x (mm)');
ylabel('Normalized Intensity');
title(['Fresnel and Fraunhofer Diffraction Pattern, \Deltav = ', num2str(delta_v, 2), ' ']);
hold on;
plot(x*1000, I_norm_fresnel);
legend('Theoretical Fraunhofer Intensity','Theoretical Fresnel Intensity')