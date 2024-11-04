close all
clear all
[num,txt,raw] = xlsread('100 Hz 1000 Samples.xlsx');

features_raw=num(1:end,7:1006);  
age=num(1:end,3:3);
histogram(age);


% for k=1:size(features_raw,1) 
% input(k,:)=features_raw(k,:);
% fs=5; %sample rate in kHz
% order=4;   %order of filter
% fcutlow=1;   %low cut frequency in kHz
% fcuthigh=4;   %high cut frequency in kHz
% [b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
% filtsig(k,:)=filter(b,a,input);  %filtered signal
% end 
% 
% input=filtsig;
% % end
% for i=1:size(input,1)
%     for j=1:size(input,2)
%         if (isnan(input(i,j)))
%             input(i,j)=1;
%         end
%     end
% end
% % % end
% % % end
% % % end
% % 
% % for i=1:size(input,1)
% %     for j=1:(size(input,2)-1)
% %     c(i,j)=(input(i,j+1)-input(i,j))/input(i,j);
% %     end
% % end
% % 
% % for i=1:size(c,1)
% %     for j=1:size(c,2)
% %         if (isnan(c(i,j)))
% %             c(i,j)=1;
% %         end
% %     end
% % end
% 
% % for i=1:size(c,1)
% %     for j=1:size(c,2)
% %         if (c(i,j)==inf)
% %             c(i,j)=1;
% %         end
% %     end
% % end
% % [C,L] = wavedec(input,4,'db10');
% % E=appcoef(C,L,'db10');
% % [d1,d2,d3,d4] = detcoef(C,L,[1 2 3 4]);
% % features_refresh(k,:)=E;
% % features_refresh(k,:)=horzcat(E,[d1,d2,d3,d4]);
% % input1 = awgn(input,10);
% % features(i,:) = filtsig;
% % features_refresh(K,:)=wprec(appcoef(C,L,'db10'));
% 
% 
% % for i=1:size(features,1) 
% %     for j=1:size(data,1)
% % % features_converted(i,j)= sum(features(i,:) .* features(j,:) ) ;
% % features_converted(i,:)= conv(features(i,:),data(j,:) ) ;
% %     end
% % end
% 
% 
% % for i=1:size(features_converted,1)
% %     for j=1:size(features_converted,2)
% %         if (isnan(features_converted(i,j)))
% %             features_converted(i,j)=1;
% %         end
% %     end
% % end
% 
% % table=c;
% % numFeatures = size(table,2);
% % numResponses = 1;
% % numHiddenUnits = 44;
% % 
% % 
% % XTrain = transpose(table(1:1910,1:end));
% % YTrain = transpose(age(1:1910,1:end));
% % XTest= transpose(table(1911:end,1:end));
% % YTest = transpose(age(1911:end,1:end));
% % layers = [ ...
% %   sequenceInputLayer(numFeatures)
% % flattenLayer
% %  lstmLayer(numHiddenUnits)
% % fullyConnectedLayer(numResponses)
% %  regressionLayer];
% % 
% % 
% % 
% % options = trainingOptions('sgdm', ...
% %     'MaxEpochs',1, ...
% %     'GradientThreshold',1, ...
% %     'InitialLearnRate',0.05, ...
% %     'LearnRateSchedule','piecewise', ...
% %     'LearnRateDropPeriod',125, ...
% %     'LearnRateDropFactor',0.1, ...
% %     'Verbose',false, ...
% %     'Plots','training-progress');
% % 
% % net= trainNetwork(XTrain,YTrain,layers,options);
% % 
% % test_outcome=round(predict(net,XTest));
% % test_outcome(isnan(test_outcome))=2;
% 
% 
% % test_accuracy=(size((find(test_outcome==YTest)),2)/size(YTest,2))*100;
% % E=test_outcome-YTest;
% % 
% % for i=1:size(YTest,2)
% %    if (YTest(i)<18)
% %           actual(i)=1;
% %    else
% %        actual(i)=2;
% %    end
% % end
% % for i=1:size(test_outcome,2)
% %    if (test_outcome(i)<18)
% %           predicted(i)=1;
% %    else
% %        predicted(i)=2;
% %    end
% % end
% % accuracy=(size((find(actual==predicted)),2)/size(predicted,2))*100;
% 
% predicted=test_outcome(3,:)*3;
% actual=YTest(3,:)'
% plot(actual,predicted);
% E=predicted-actual;
% MAE= mae(E,predicted,actual);
% 
% 
