close all
clear all
[num,txt,raw] = xlsread('wavenet.xlsx');

features_raw=num(1:end,8:1075);  
age=num(1:end,3:3);


features=features_raw;
kernel1=[-32:32];
increment1=size(kernel1,2);
for s=1:size(features,1)
    input1=features(s,:) 
res1(:,:)=conv(kernel1,input1);
i=1;
j=1;
k=1;
for k=1:(size(res1,2)/increment1)
    A1=res1(i,j:j+increment1-1);
maxpool1(k)=max(A1);
k=k+1;
j=j+increment1;
end
output1(s,:)=(maxpool1-std(maxpool1))/mean(maxpool1);
end

% kernel2=[-256:256];
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
% kernel3=[-128:128];
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

table=output1;
numFeatures = size(table,2);
numResponses = 1;
numHiddenUnits = 1000;
% 
% 


XTest = transpose(table(2388:end,1:end));
YTest = transpose(age(2388:end,1:end));
XTrain= transpose(table(1:2387,1:end));
YTrain = transpose(age(1:2387,1:end));
layers = [ ...
  sequenceInputLayer(numFeatures)
flattenLayer
lstmLayer(numHiddenUnits)
fullyConnectedLayer(numResponses)
 regressionLayer];
options = trainingOptions('adam', ...
    'MaxEpochs',10, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.001, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',125, ...
    'LearnRateDropFactor',0.01, ...
    'Verbose',false, ...
    'Plots','training-progress');

net= trainNetwork(XTrain,YTrain,layers,options);
test_outcome=round(predict(net,XTest));
test_outcome(isnan(test_outcome))=2;
test_accuracy=(size((find(test_outcome==YTest)),2)/size(YTest,2))*100;
E=test_outcome-YTest;

% for i=1:size(YTest,2)
%    if (YTest(i)<18)
%           actual(i)=1;
%    else
%        actual(i)=2;
%    end
% end
% for i=1:size(test_outcome,2)
%    if (test_outcome(i)<18)
%           predicted(i)=1;
%    else
%        predicted(i)=2;
%    end
% end
% accuracy=(size((find(actual==predicted)),2)/size(predicted,2))*100;
plot(YTest,test_outcome);
% E=test_outcome-YTest;
% MAE= mae(E,YTest,test_outcome);


