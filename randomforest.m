close all

clear all
close all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segment.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      

features_raw=num(1:end,9:1008);  
age=num(1:end,5:5);

input=70;
ratio=0.4;
threshold=90;

label=age;

% label=zeros(size(age,1),1)
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

%Find the frequency range
threshold = 7;  % Define a threshold to identify significant frequencies
significant_freqs=f(P1 > threshold);
a(i)=max(significant_freqs);
b(i)=min(significant_freqs);
c(i)=mode(significant_freqs);
featuresfreq(i,:)=significant_freqs(1:10);
end
mdl = fitlm(featuresfreq,age)
f=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
signal=features_raw(i,:);
fs=500; %sample rate in kHz
order=2;   %order of filter
fcutlow=1;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,signal);  %filtered signal
% Y(i,:) = real(fft(filtsig,4));
% filtsignal=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'sym4');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'sym4');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuressym(i,:)=[E,d1,d2,d3,d4];

% featuressym(i,:)=appcoef(C,L,'sym4',4);
end
mdl = fitlm(featuressym,age)
f5=mdl.Rsquared.ordinary*1e4;

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
% filtsignal=features_raw(i,:);
[C,L]=wavedec(filtsig,4,'db1');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db1');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb1(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb1,age)
f1=mdl.Rsquared.ordinary*1e4;

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
% filtsignal=features_raw(i,:);
[C,L]=wavedec(filtsig,4,'db2');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db2');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb2(i,:)=[E,d1,d2,d3,d4];
% featuresdb9(i,:)=appcoef(C,L,'db9',4);
end
mdl = fitlm(featuresdb2,age)
f2=mdl.Rsquared.ordinary*1e4;


for i=1:size(features_raw,1) 
%     for i=1:15
% filtsignal=features_raw(i,:);
fs=500; %sample rate in kHz
order=2;   %order of filter
fcutlow=3;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,signal);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,10,'db3');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db3');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb3(i,:)=[E,d1,d2,d3,d4];
% featuresdb10(i,:)=appcoef(C,L,'db10',4);
end
mdl = fitlm(featuresdb3,age)
f3=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
signal=features_raw(i,:);
fs=500; %sample rate in kHz
order=2;   %order of filter
fcutlow=3;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,signal);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db4');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db4');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb4(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb7(i,:)=appcoef(C,L,'db7',4);
end
mdl = fitlm(featuresdb4,age)
f4=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=1;   %order of filter
fcutlow=4;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db5');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db5');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb5(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb4(i,:)=appcoef(C,L,'db4',4);
end
mdl = fitlm(featuresdb5,age)
f5=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=1;   %order of filter
fcutlow=4;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db6');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db6');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb6(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb4(i,:)=appcoef(C,L,'db4',4);
end
mdl = fitlm(featuresdb6,age)
f6=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=1;   %order of filter
fcutlow=4;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db7');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db7');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb7(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb4(i,:)=appcoef(C,L,'db4',4);
end
mdl = fitlm(featuresdb7,age)
f7=mdl.Rsquared.ordinary*1e4;


for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=1;   %order of filter
fcutlow=4;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db8');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db8');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb8(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb4(i,:)=appcoef(C,L,'db4',4);
end
mdl = fitlm(featuresdb8,age)
f8=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=1;   %order of filter
fcutlow=4;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db9');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db9');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb9(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb4(i,:)=appcoef(C,L,'db4',4);
end
mdl = fitlm(featuresdb9,age)
f9=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
%     for i=1:15
input=features_raw(i,:);
fs=500; %sample rate in kHz
order=1;   %order of filter
fcutlow=4;   %low cut frequency in kHz
fcuthigh=4;   %high cut frequency in kHz
[b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
filtsig=filter(b,a,input);  %filtered signal
% filtsignal=features_raw(i,:);
% Y(i,:) = real(fft(filtsig,4));
[C,L]= wavedec(filtsig,4,'db10');
% x = waverec(C,L,'db10');
% features(i,:)=x;
E=appcoef(C,L,'db10');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb10(i,:)=[E,d1,d2,d3,d4];
filtsignal=features_raw(i,:);
% featuresdb4(i,:)=appcoef(C,L,'db4',4);
end
mdl = fitlm(featuresdb10,age);
f10=mdl.Rsquared.ordinary*1e4;

w1=f1/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w2=f2/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w3=f3/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);;
w4=f4/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w5=f5/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w6=f6/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w7=f7/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w8=f8/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w9=f9/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w10=f10/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w=f/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);

table=cat(2,featuresdb1*w1,featuresdb2*w2,featuresdb3*w3,featuresdb4*w4,featuresdb5*w5,featuresfreq*w,featuresdb6*w6,featuresdb7*w7,featuresdb8*w8,featuresdb9*w9,featuresdb10*w10);
%table=cat(2,featuresdb4*w1,featuresdb7*w4,featuresdb10*w3,featuresdb9*w2,featuresfreq*w6,featuresfreq*w5);
% table=cat(2,featuresdb4*w1,featuresfreq*w6,featuressym*w5);

% table=awgn(table,20);

% while(count<2)
 
close all
  % if(count==1);
  %     table=featuressym;
  % elseif(count==2) 
  %     table=featuresdb;
  %      elseif(count==3) 
  %     table=featuresfreq;
  %      elseif(count==4) 
  %     table=cat(2,featuressym,featuresdb);
  %      elseif(count==5) 
  %     table=cat(2,featurefreq,featuresdb);
  %        elseif(count==6) 
  %     table=cat(2,featurefreq,featuressym);
  %     elseif(count==7) 
% 
% R=rand(1,5);
% w1=R(1);
% w2=R(2);
% w3=R(3);
% w4=R(4);
% w5=R(5);
% table=cat(2,featuresdb4*w1,featuresdb7*w2,featuresdb10*w3,featuresdb9*w4,featuressym*w5);
 
% table=featuressym;

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
table=features_raw;
numFeatures = size(table,2);
numResponses = 1;
numHiddenUnits = 1000;

XTest = (table(1:2388,1:end));
YTest = (age(1:2388,1:end));
XTrain= (table(2389:end,1:end));
YTrain = (age(2389:end,1:end));

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

numTrees = 1000; % Number of trees
rfModel = TreeBagger(numTrees, XTrain, YTrain, 'OOBPrediction', 'On', 'Method', 'regression');

% View the out-of-bag error
figure;
oobErrorBaggedEnsemble = oobError(rfModel);
plot(oobErrorBaggedEnsemble);
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');

% Predict the test set
[Y_pred, scores] = predict(rfModel, XTest);
scatter(YTest,Ypred);
E=YTest-Ypred;
MAE= mae(E,YTest,Ypred);

E=YTest-Y_pred;
MAE= mae(E,YTest,Y_pred);
% 
% % Display the confusion matrix
% figure;
% confusionchart(Y_test, Y_pred);
% title('Confusion Matrix');

% Display feature importance
% figure;
% bar(rfModel.OOBPermutedPredictorDeltaError);
% xlabel('Feature Index');
% ylabel('Out-Of-Bag Feature Importance');
% title('Feature Importance');

% layers = [ ...
%   sequenceInputLayer(numFeatures)
% flattenLayer
%  lstmLayer(numHiddenUnits)
% fullyConnectedLayer(numResponses)
%  regressionLayer];
% 
% options = trainingOptions('adam', ...
%     'MaxEpochs',20, ...
%     'GradientThreshold',1, ...
%     'InitialLearnRate',0.01, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropPeriod',125, ...
%     'LearnRateDropFactor',0.01, ...
%     'Verbose',false, ...
%     'Plots','training-progress');
% 
% net= trainNetwork(XTrai,YTrai,layers,options);
% 
% test_outcome=round(predict(net,XTes));
% 
% multiplier=1;
% predicted=test_outcome*multiplier;
% actual=YTest;
% scatter(actual,predicted);
E=predicted-actual;
MAE= mae(E,predicted,actual);
% count=count+1;
% for i = 1:size(test_outcome, 2)
%     if (ceil(test_outcome(3, i) / 5) ~= test_outcome(2, i)) || (ceil(test_outcome(3, i) / 30) ~= test_outcome(1, i))
%         index(i) = 0;

%     else
%         index(i) = 1;
%     end
% end
% 
% correct = sum(test_outcome==YTes);
% accuracy= (correct/size(test_outcome,2))*100;



  % end
% row = find(index < val)