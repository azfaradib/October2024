clear all
close all
% [num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segment.xlsx');
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segment.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      
table=num;
coeff = pca(table);
