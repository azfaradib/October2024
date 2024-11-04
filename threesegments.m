
clear all
close all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segments.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      

features_raw=num(1:end,8:end);  
age=num(1:end,5:5);


for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=2;   %order of filter
fcutlow=1;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'sym4');
% x = waverec(C,L,'db10');
% features(i,:)=x;
features(i,:)=appcoef(C,L,'sym4',4);
end

table=features_raw;

% XVal = transpose(table(21:30,1:end));
% YVal = transpose(age(21:30,1:end));

XTrain = transpose(table(2388:end,1:end));
YTrain = transpose(age(2388:end,1:end));
XTest= transpose(table(1:2387,1:end));
YTest = transpose(age(1:2387,1:end));

numFeatures = size(table,2);
numResponses = 1;
numHiddenUnits = 10000;


layers = [ ...
    sequenceInputLayer(numFeatures)
    %flattenLayer
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
     regressionLayer];

options = trainingOptions('sgdm', ...
    'MaxEpochs',15, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.05, ...
    'Verbose',false, ...
    'Plots','training-progress');


net= trainNetwork(XTrain,YTrain,layers,options);
test_outcome=round(predict(net,XTest));

% for i = 1:size(test_outcome, 2)
%     if (ceil(test_outcome(3, i) / 5) ~= test_outcome(2, i)) || (ceil(test_outcome(3, i) / 30) ~= test_outcome(1, i))
%         index(i) = 0;
%     else
%         index(i) = 1;
%     end
% end

% correct = sum(index== 1);
% accuracy= (correct/size(test_outcome,2))*100

multiplier=1
predicted=test_outcome*multiplier;
actual=YTest;
scatter(actual,predicted);
E=predicted-actual;
MAE= mae(E,predicted,actual);


