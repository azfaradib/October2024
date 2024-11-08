clear all
close all
[num,txt,raw] = xlsread('100 Hz 1000 Samples 3 segment.xlsx');
Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period      

% table=num;
%  for i=1:size(table,1)
%     for j=1:size(table,2)
%         if (table(i,j)==NaN)
%             table(i,j)=0;
%         end
%     end
%  end

features=num(1:end,9:1008);  
age=num(1:end,5:5);


% for i=1:size(features_raw,1) 
% Fs = 1000;                    % Sampling frequency (Hz)
% T = 1/Fs;                     % Sampling period (s)
% L = 1000;                     % Length of signal
% t = (0:L-1)*T;                % Time vector
% Y = fft(features_raw(i,:));
% % Compute the two-sided spectrum
% P2 = abs(Y/L);
% % Compute the single-sided spectrum based on the two-sided spectrum
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% % Define the frequency domain f
% f = Fs*(0:(L/2))/L;
% 
% % Find the frequency range
% threshold = 7;  % Define a threshold to identify significant frequencies
% significant_freqs=f(P1 > threshold);
% a(i)=max(significant_freqs);
% b(i)=min(significant_freqs);
% c(i)=mode(significant_freqs);
% featuresfreq(i,:)=significant_freqs(1:10);
% end

features_raw=[];

for i=1:size(features,1) 
%     for i=1:15
signal=features(i,:);
fs=1000; %sample rate in kHz
order=4;   %order of filter
% Define the cutoff frequencies (in Hz)
f_low = 1.1;    % Lower cutoff frequency
f_high = 2.9;    % Higher cutoff frequency
% Normalize the cutoff frequencies with respect to Nyquist frequency
Wn = [f_low f_high] / (fs/2);
% Design a 4th order Butterworth filter
[b, a] = butter(order, Wn, 'bandpass');
filtsig=filter(b,a,signal);  %filtered signal
features_raw(i,:)=filtsig;
end


for i=1:size(features_raw,1) 
 signal=features_raw(i,:);   
[C,L]= wavedec(signal,4,'sym4');
E=appcoef(C,L,'sym4');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% featuressym(i,:)=E;
featuressym(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuressym,age);
f5=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
signal=features_raw(i,:);
[C,L]=wavedec(filtsig,4,'db1');
E=appcoef(C,L,'db1');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% featuresdb1(i,:)=E;
featuresdb1(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb1,age);
f1=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
signal=features_raw(i,:);
filtsig=filter(b,a,signal);  %filtered signal
[C,L]=wavedec(filtsig,4,'db2');
E=appcoef(C,L,'db2');
% featuresdb2(i,:)=E;
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb2(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb2,age);
f2=mdl.Rsquared.ordinary*1e4;


for i=1:size(features_raw,1) 
filtsig=filter(b,a,signal);  %filtered signal
[C,L]= wavedec(filtsig,4,'db3');
E=appcoef(C,L,'db3');
% featuresdb3(i,:)=E;
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
featuresdb3(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb3,age);
f3=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
signal=features_raw(i,:);
filtsig=filter(b,a,signal);  %filtered signal
[C,L]= wavedec(filtsig,4,'db4');
E=appcoef(C,L,'db4');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% featuresdb4(i,:)=E;
featuresdb4(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb4,age);
f4=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
filtsig=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'db5');
E=appcoef(C,L,'db5');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% featuresdb5(i,:)=E;
featuresdb5(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb5,age);
f5=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
filtsig=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'db6');
E=appcoef(C,L,'db6');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% featuresdb6(i,:)=E;
featuresdb6(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb6,age);
f6=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
filtsig=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'db7');
E=appcoef(C,L,'db7');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% featuresdb7(i,:)=E;
featuresdb7(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb7,age);
f7=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
filtsig=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'db8');
E=appcoef(C,L,'db8');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
%featuresdb8(i,:)=E;
featuresdb8(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb8,age);
f8=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
filtsig=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'db9');
E=appcoef(C,L,'db9');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
%featuresdb9(i,:)=E;
featuresdb9(i,:)=[E,d1,d2,d3,d4];
end
mdl = fitlm(featuresdb9,age);
f9=mdl.Rsquared.ordinary*1e4;

for i=1:size(features_raw,1) 
filtsignal=features_raw(i,:);
[C,L]= wavedec(filtsig,4,'db10');
E=appcoef(C,L,'db10');
[d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
%featuresdb10(i,:)=E;
featuresdb10(i,:)=[E,d1,d2,d3,d4];
end
%transpose(mscohere(base,features_raw(i,:)));mdl = fitlm(featuresdb10,age);
f10=mdl.Rsquared.ordinary*1e4;

% w1=f1/(f1+f2+f3+f4);
% w2=f2/(f1+f2+f3+f4);
% w3=f3/(f1+f2+f3+f4);;
% w4=f4/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
f=0;
w1=f1/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w2=f2/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
w3=f3/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);;
w4=f4/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w5=f5/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w6=f6/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w7=f7/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w8=f8/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w9=f9/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w10=f10/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w=f/(f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f);
% w1=w1/4;
% w2=w2/4;
% w3=w3/4;
% w4=w4/4;

% % w=1;
% w1=1;
% w2=1;
% w3=1;
% w4=1;
% w5=1;
% w6=1;
% w7=1;
% w8=1;
% w9=1;
% w10=1;


% for i=1:size(features_raw,1)
%     table(:,:,i)= [
%     featuresdb1(1:950);
%     featuresdb2(1:950);
%     featuresdb3(1:950);
%     featuresdb4(1:950);
%     featuresdb5(1:950);
%     featuresdb6(1:950);
%     featuresdb7(1:950);
%     featuresdb8(1:950);
%     featuresdb9(1:950);
%     featuresdb10(1:950);
%     ]
% end

%table=cat(2,featuresdb1*w1,featuresdb2*w2,featuresdb3*w3,featuresdb4*w4,featuresdb5*w5,featuresdb6*w6,featuresdb7*w7,featuresdb8*w8,featuresdb9*w9,featuresdb10*w10);
% table=featuresdb1*w1;
% kernel1=[-16:16];
% increment1=size(kernel1,2);
% for s=1:size(table,1)
%     input1=table(s,:) ;
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

% % table=cat(2,featuresdb4*w1,featuresdb7*w4,featuresdb10*w3,featuresdb9*w2,featuresfreq*w6,featuresfreq*w5);
 table=cat(2,featuresdb1*w1,featuresdb2*w2,featuresdb3*w3,featuresdb4*w4);


%  for i=1:size(table,1)
%  table(i)=awgn(table(i),5);
%  end
% 
%  %table=featuresdb4*w4;
% % 
% % % table=awgn(table,20);
% % % 
% % % while(count<2)
% % 
% % close all
% 
numFeatures = size(table,2);
numResponses = 1;
numHiddenUnits = 1000+


































































+0;
% % 
% for i=1:size(age,1)
%     if age(i)<n 
%         index(i)=1;
%     else
%         index(i)=2;
%     end


XTrain = (table(1:2388,1:end));
YTrain = (age(1:2388,1:end));
XTest= (table(2389:end,1:end));
YTest = (age(2389:end,1:end));


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
    'MaxEpochs',5, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.1, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.01, ...
    'Verbose',false, ...
    'Plots','training-progress');

net= trainNetwork(XTrai,YTrai,layers,options);

test_outcome=round(predict(net,XTes));

multiplier=1;
predicted=test_outcome*multiplier;
actual=YTest;
scatter(actual,predicted);
E=predicted-actual;
MAE= mae(E,predicted,actual);


