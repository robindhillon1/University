% Importing the relevant diffraction images
I_010 = readAndRotate('11.60mm.tif');
%I_030 = readAndRotate('0-3_1_polar1.tif');
%I_060 = readAndRotate('0-60_polar1.tif');
%I_100 = readAndRotate('1-0_1_polar1.tif');
%I_114 = readAndRotate('1-14_1_polar1.tif');
%I_150 = readAndRotate('1-50_polar1.tif');
%I_200 = readAndRotate('2-0_1_polar1.tif');
%I_300 = readAndRotate('3-0_polar1.tif');

% Combining the diffraction images for visualization
%combinedIm = cat(1,I_010,I_030,I_060,I_114,I_200,I_300);
%combinedIm8bit = uint8(combinedIm / 256);
%imwrite(combinedIm8bit, 'Diffraction.jpg')

% Calculating deltaV values ~ albeit inefficiently
v010 = deltaV(1.88*10^-3);
%v030 = deltaV(0.3*10^-3);
%v060 = deltaV(0.6*10^-3);
%v114 = deltaV(1.14*10^-3);
%v200 = deltaV(2*10^-3);
%v300 = deltaV(3*10^-3);

% Plotting the Image
figure
% z = -5:0.01:5;
%plotDif(I_010);
%hold on;
%plotDiff(I_030);
%hold on;
%plotDiff(I_060);
%imshow(I2)
%improfile;
%plot([-10:0.01:10],fraunhofer(1))

% Plotting the image data following some normalization
dpi = 600;
width = 1.88*10^-3;
[x,y,profile] = plotDiff(I_010);
x = ((x-382)*0.0254/dpi)/width; %1305 is the hardcoded position of the peak
plot(x, profile/max(profile))
hold on;

F = fresnel(v010,x);
plot(x,F/max(F))
hold on;

%Fr = fraunhofer(v114,x);
%plot(x,Fr/max(Fr))
%hold on;

% 
% hold on;
% plot(x,log(fresnel(0.088,x)/max(fresnel(0.088,x)))
% hold on;
% plot(x,fraunhofer(0.088,x)/max(fraunhofer(0.088,x)))
% set(gca, 'YScale', 'log')
% xlim([-150, 150])


function I = readAndRotate(path)
    I = imread(path);
    I = imrotate(I,90);
end

function [x,y,output] = plotDiff(I)
    x = [0 size(I,2)];
    y = [size(I,1)/3 size(I,1)/3];
    
    % Gaussian Filter to reduce noise
    I = imgaussfilt(I,2);

    [x,y,output] = improfile(I,x,y);
end

function output = fraunhofer(v, z)
    b = z*(pi/2);
    output = (v^2).*(sin(b).^2)./(b.^2);
end

function output = fresnel(v,z)
    output =arrayfun(@(x) fresnelInt(v,x),z);
end


function output = fresnelInt(v, z)
    func1 = @(x) cos((pi/2)*x.^2);
    func2 = @(x) sin((pi/2)*x.^2);
    c = integral(func1, -(z+0.5)*v, -(z-0.5)*v);
    s = integral(func2, -(z+0.5)*v, -(z-0.5)*v);
    output = c^2+s^2;
end

function v = deltaV(w)
    v = w*sqrt(2/(2.7*632.8*10^-9));
end
