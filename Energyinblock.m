clear all
close all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segment.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      
% age=xlsread('age.xlsx')
% table=num;

features_raw=num(1:end,9:1008);  
age=num(1:end,5:5);

Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      
L = 1000; 
features_raw=num(1:end,7:end);  
% plot(features_raw)
for i=1:size(features_raw,1)
c=features_raw(i,:);
% F = fft(c); 
% pow(i,:) = F.*conj(F);
% power(i)=sum(pow(i,:));
%  min(i)=min(features(i
psdestx = psd(spectrum.periodogram,c,'Fs',1e3,'NFFT',length(c));
power(i) = avgpower(psdestx);
end

