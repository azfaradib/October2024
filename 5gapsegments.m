
close all
clear all
[num,txt,raw] = xlsread('100 Hz 1000 Samples.xlsx');

features_raw=num(1:end,8:1006);  
age=num(1:end,3:3);
% 
% segmentsize=200;
% 
% for k=1:size(features_raw,1) 
% input=features_raw(k,:);
% fs=500; %sample rate in kHz
% order=4;   %order of filter
% fcutlow=1.1;   %low cut frequency in kHz
% fcuthigh=2.9;   %high cut frequency in kHz
% [b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
% filtsig=filter(b,a,input);  %filtered signal
% features_filtered(k,:)=filtsig;
% end

% for i=1:size(features_raw,1)
%  input=features_raw (i,:);
% [C,L] = wavedec(input,4,'db10');
% E=appcoef(C,L,'db10');
% [d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% features_wavelet(i,:)=[E,d1,d2,d3,d4];
% end

%Apply k-means clustering with k=3 

kernel1=[-16:16];
increment1=size(kernel1,2);
for s=1:size(features_raw,1)
    input1=features_raw(s,:) ;
res1(s,:)=conv(kernel1,input1);
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
end

% kernel2=[-16:16];
% increment2=size(kernel2,2);
% for s=1:size(output1,1)
%     input2=output1(s,:);
% 
% res2(:,:)=conv(kernel2,input2);
% i=1;
% j=1;
% k=1;
% for k=1:(size(res2,2)/increment2)
%     A2=res2(i,j:j+increment2-1);
% maxpool2(k)=max(A2);
% k=k+1;
% j=j+increment2;
% end
% output2(s,:)=(maxpool2-std(maxpool2))/mean(maxpool2);
% end
% 
% 
% kernel3=[-8:8];
% increment3=size(kernel3,2);
% for s=1:size(output2,1)
%     input3=output2(s,:);
% 
% res3(:,:)=conv(kernel3,input3);
% i=1;
% j=1;
% k=1;
% for k=1:(size(res3,2)/increment3)
%     A3=res3(i,j:j+increment3-1);
% maxpool3(k)=max(A3);
% k=k+1;
% j=j+increment3;
% end
% output3(s,:)=(maxpool3-std(maxpool3))/mean(maxpool3);
% end
% 
% features_conv=res1;
% 
% [idx,C]=kmeans(features_conv,3,'Distance','cityblock');
% plot(age,idx);

table=features_conv;
numFeatures = size(table,2);
numResponses = 1;
numHiddenUnits = 10000;

XTest= (table(1:1910,1:end));
YTest =(age(1:1910,1:end));
XTrain= (table(1911:end,1:end));
YTrain =(age(1911:end,1:end));


% XTrain = transpose(table(1:1849,1:end));
% YTrain = transpose(age(1:1849,1:end));
% % XTrain=[XTrain1 XTrain1];
% % YTrain=[YTrain1 YTrain1];;
% XTest= transpose(table(1850:end,1:end));
% YTest = transpose(age(1850:end,1:end));

% 
% tree = fitrtree(XTrain,YTrain);
% test_outcome= predict(tree,XTest);
% 
% layers = [ ...
%   sequenceInputLayer(numFeatures)
% flattenLayer
%  lstmLayer(numHiddenUnits)
% fullyConnectedLayer(numResponses)
%  regressionLayer];
% 
% options = trainingOptions('adam', ...
%     'MaxEpochs',10, ...
%     'GradientThreshold',1, ...
%     'InitialLearnRate',0.1, ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropPeriod',125, ...
%     'LearnRateDropFactor',0.1, ...
%     'Verbose',false, ...
%     'Plots','training-progress');
% 
% net= trainNetwork(XTrain,YTrain,layers,options);
% 
% test_outcome=round(predict(net,XTest));
% % test_outcome(isnan(test_outcome))=2;
% 
% 
% test_accuracy=(size((find(test_outcome==YTest)),2)/size(YTest,2))*100;
% E=test_outcome-YTest;
% % 
% % % for i=1:size(YTest,2)
% % %    if (YTest(i)<18)
% % %           actual(i)=1;
% % %    else
% % %        actual(i)=2;
% % %    end
% % % end
% % % for i=1:size(test_outcome,2)
% % %    if (test_outcome(i)<18)
% % %           predicted(i)=1;
% % %    else
% % %        predicted(i)=2;
% % %    end
% % % end
% % % accuracy=(size((find(actual==predicted)),2)/size(predicted,2))*100;
% plot(YTest,test_outcome);
% E=test_outcome-YTest;
% MAE= mae(E,YTest,test_outcome);
% 
% 
