close all
clear all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segments.xlsx');

features_raw=num(1:end,7:1006);  
age=num(1:end,3:3);
max=6;

clear all
close all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segments.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      

features_raw=num(1:end,8:end);  
age=num(1:end,5:5);

input=70;
ratio=0.4;
threshold=90;


label=zeros(size(age,1),1)
[row, col] = find(age>input);

    for i=1:size(row,1)
       label(row(i))=1;
    end

[row, col]=find(label==input);

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
end

% Find the frequency range
% threshold = 7;  % Define a threshold to identify significant frequencies
significant_freqs=f(P1 > threshold);
a(i)=max(significant_freqs);
b(i)=min(significant_freqs);
c(i)=mode(significant_freqs);
featurefreq(i,:)=significant_freqs(1:300);
end

for i=1:size(features_raw,1) 
%     for i=1:15
signal=features_raw(i,:);
fs=500; %sample rate in kHz
order=2;   %order of filter
fcutlow=3;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,signal);  %filtered signal
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'sym4');
% x = waverec(C,L,'db10');
% features(i,:)=x;
featuressym(i,:)=appcoef(C,L,'sym4',4);
end


for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=2;   %order of filter
fcutlow=3;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db10');
% x = waverec(C,L,'db10');
% features(i,:)=x;
featuresdb(i,:)=appcoef(C,L,'db10',4);
end
% 

table=cat(2,featuressym,featuresdb,featurefreq);



table=featuressym;

% kernel1=[-512:512];
% increment1=size(kernel1,2);
% for s=1:size(table,1)
%     input1=table(s,:) 
% res1(:,:)=conv(kernel1,input1);
% i=1;
% j=1;
% k=1;
% for k=1:(size(res1,2)/increment1)
%     A1=res1(i,j:j+increment1-1);
% maxpool1(k)=max(A1);
% k=k+1;
% j=j+increment1;
% end
% output1(s,:)=(maxpool1-std(maxpool1))/mean(maxpool1);
% end

% table=output1;

numFeatures = size(table,2);
numResponses = 1;
numHiddenUnits = 100;
% 
% 

XTest = (table(1:1910,1:end));
YTest = (label(1:1910,1:end));
XTrain= (table(1911:end,1:end));
YTrain = (label(1911:end,1:end));

% index=transpose(age(1911:end,1:end));
% lower= find(index < input);
% upper=find(index > input);
% 
% if (size(lower,2)/size(XTrain,2))<=ratio
%     loweramp=1;
% else
%  loweramp=0;
% end
% 
% if (size(upper,2)/size(XTrain,2))<=ratio
%     upperamp=1;
% else
%  upperamp=0;
% end
% 
% accuracy=0;


for count=1:5

% if loweramp==1
%     XTrain=repmat(XTrain,1,5);
%      YTrain=repmat(YTrain,1,5);
%     % for i=1:size(lower,1)
%     %     k=lower(i,1);
%     %    [temp1]=XTrain(k,:);
%     %    [temp2]=YTrain(k,:);
%     %    XTrain=[XTrain; temp1];
%     %    YTrain=[YTrain;temp2];
%     % 
%     % end
% end
% 
%     if upperamp==1
%          XTrain=repmat(XTrain,1,5);
%      YTrain=repmat(YTrain,1,5);
%     %     for i=1:size(upper,2)
%     %           k=upper(i,1);
%     %    [temp1]=XTrain(k,:);
%     %    [temp2]=YTrain(k,:);
%     %    XTrain=[XTrain; temp1];
%     %    YTrain=[YTrain;temp2];
%     % 
%     % end
% end
% 

XTes = transpose(XTest);
YTes = transpose(YTest );
XTrai= transpose(XTrain);
YTrai = transpose(YTrain);

layers = [ ...
  sequenceInputLayer(numFeatures)
flattenLayer
 lstmLayer(numHiddenUnits)
fullyConnectedLayer(numResponses)
 regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',1, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.01, ...
    'Verbose',false, ...
    'Plots','training-progress');

net= trainNetwork(XTrai,YTrai,layers,options);

test_outcome=round(predict(net,XTes));

% for i = 1:size(test_outcome, 2)
%     if (ceil(test_outcome(3, i) / 5) ~= test_outcome(2, i)) || (ceil(test_outcome(3, i) / 30) ~= test_outcome(1, i))
%         index(i) = 0;

%     else
%         index(i) = 1;
%     end
% end

correct = sum(test_outcome==YTes);
accuracy= (correct/size(test_outcome,2))*100;
count=count+1;
  end

% row = find(index < val)