
clear all
close all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segments.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      

features_raw=num(1:end,8:end);  
age=num(1:end,3:5);

for i=1:size(features_raw,1) 
Fs = 1000;                    % Sampling frequency (Hz)
T = 1/Fs;                     % Sampling period (s)
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
Y = fft(features_raw(i,:));
% Compute the two-sided spectrum
P2 = abs(Y/L);
% Compute the single-sided spectrum based on the two-sided spectrum
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
% Define the frequency domain f
f = Fs*(0:(L/2))/L;

% Find the frequency range
threshold = 0.1;  % Define a threshold to identify significant frequencies
significant_freqs=f(P1 > threshold);
a(i)=max(significant_freqs);
b(i)=min(significant_freqs);
c(i)=mode(significant_freqs);
end

figure(1)
scatter(age,a)
figure(2)
scatter(age,b)
figure(3)
scatter(age,c)
