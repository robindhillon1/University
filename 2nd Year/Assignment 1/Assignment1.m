clc;
clear;
figure(1)
fid = fopen('assign1_1.dat','r');
datacell = textscan(fid, '%f%f%f', 'HeaderLines', 1, 'Collect', 1);
fclose(fid);
A.data = datacell{1};
plot(A.data(:,1),A.data(:,2),'.r')
hold on;
grid on;

%The following model uses the given exponential function, and will be
%printed in a dotted line. 
x=A.data(:,1);

y = 9.4*exp(x/-21.7)+0.796;%+0.1*randn(size(x));
%Here, I used a = 10, b = -27, c = 0. The uncertainties in the values are
%less than 5%, as the values are within 95% cofidence bounds. 
f = fit(x,y,'exp1')
plot(f,'--k')
xlabel('x-data')
ylabel('y-data')
title('Data: y vs. x')
legend('Data','Fitted Curve')
%numbers: 9.402,-21.7,0.796
%uncertainties: 0.0273,0.1437,0.01347

%{
Matlab outputs the following:

General model Exp1:
     f(x) = a*exp(b*x)
     Coefficients (with 95% confidence bounds):
       a =          10  (9.999, 10)
       b =    -0.03705  (-0.03706, -0.03703)

Since I used c=0, just computes a and b. The output of the computer and my
choice was 10, which is good. However, the b value is too small, which
seems incorrect. Upon using the b proved by the computer, the fitted line
is completely off. So it's good not to straight trust the numbers given by
the computer. a looks reasonable, but b does not. 
%}

%% Part 2
figure(2)
F = importdata('assign1_2.dat');

c1=0.29594921;
s1=2.02564261;
u1=2.20713994;
c2=0.70392872;
s2=2.24227454;
u2=-1.87586983;
m = -((F-u1).^2)/(2*s1^2);
n = -((F-u2).^2)/(2*s2^2);
yy = ((abs(c1).*(exp(m))/sqrt(2*pi*s1.^2)))+((abs(c2).*(exp(n)))/sqrt(2*pi*s2.^2));
ff = fit(F,yy,'exp2')

HH = histfit(F,100); %Chose 100 bins
xlabel('x-axis')
ylabel('y-axis')
title('assign1.2 data Histogram')
legend('Data','Best-Fit')
grid on;
%At a glance, yes, the histogram looks like a Gaussian-distribution. But
%upon plotting the best fit through the computer, we can see that it isn't.
%Furthermore, applied more checks:

%{
I applied the jbtest, Jarque-Bera hypothesis test of composite normality,on
data. "jbtest(X) performs the Jarque-Bera goodness-of-fit test of composite
normality, i.e., that the data in the vector X came from an unspecified
normal distribution, and returns the result of the test in H. H=0 indicates
that the null hypothesis ('the data are normally distributed') cannot be
rejected at the 5% significance level. H=1 indicates that the null 
hypothesis can be rejected at the 5% level". 
%}

H = jbtest(F);

%Upon applying the jbtest, it returns a value of 1 for the data, which
%shows that the data is not Gaussian-distributed. Another test called
%the Lilliefors test(function name is lillietest), and Anderson-Darling test,adtest,
%undergo the same procedure and output 1 as well, confirming non-Gaussian-
%distributed data.

%Integrates to 1. 
%{
c1=0.29594921;
s1=2.02564261;
u1=2.20713994;
c2=0.70392872;
s2=2.24227454;
u2=-1.87586983;
%}