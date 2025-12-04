% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

addpath(genpath('functions'))
fnames = dir('Dataset1_CNN_features/*.mat');
ns = length(fnames);
rng(42, 'twister') % for reproducibility

%% load data
for subj = 1:ns
    load([fnames(subj).folder '/' fnames(subj).name])

    if subj == 1
        subjIDs = [feat_struct.subjID]';
        X_subj = unique([feat_struct.subjID]);
        y_subj = unique([feat_struct.y]);
        y = [feat_struct.y]';
        X = {feat_struct.X}';
    else
        subjIDs_tmp = cat(1,subjIDs,[feat_struct.subjID]');
        subjIDs = subjIDs_tmp;

        X_subj_tmp = cat(1,X_subj,unique([feat_struct.subjID]));
        X_subj = X_subj_tmp;

        y_subj_tmp = cat(1,y_subj,unique([feat_struct.y]));
        y_subj = y_subj_tmp;

        y_tmp = cat(1,y,[feat_struct.y]');
        y = y_tmp;

        X_tmp = cat(1,X,{feat_struct.X}');
        X = X_tmp;
    end

end

y_labels = cell(length(y),1);
y_labels(y==0) = {'HC'};
y_labels(y==1) = {'SZ'};
X_mat = zeros(size(X{1},1),size(X{1},2),size(X{1},3),length(X));
for i = 1:length(X)
    X_mat(:,:,:,i) = X{i};
end

subjID_full = arrayDatastore(categorical(subjIDs));
X_subj_full = arrayDatastore(categorical(X_subj));
y_subj_full = arrayDatastore(categorical(y_subj));
X_full = arrayDatastore(X_mat,IterationDimension=4);
y_full = arrayDatastore(categorical(y));

XT_full = combine(X_full,y_full);

%% test ResNET-18

net = resnet18;
lgraph = layerGraph(net);
new_classNames = categorical({'HC','SZ'}');

% freeze first 10 convolution layers
cnt = 0;
for i = 1:length(lgraph.Layers)
    layer = lgraph.Layers(i);
    if isa(layer,'nnet.cnn.layer.Convolution2DLayer')
        newLayer = convolution2dLayer( ...
            layer.FilterSize, ...
            layer.NumFilters, ...
            'Stride', layer.Stride, ...
            'Padding', layer.PaddingSize, ...
            'Name', layer.Name, ...
            'Weights', layer.Weights, ...
            'Bias', layer.Bias, ...
            'WeightLearnRateFactor', 0, ...
            'BiasLearnRateFactor', 0);

        lgraph = replaceLayer(lgraph, layer.Name, newLayer);
        cnt = cnt + 1;
        if cnt == 10
            break
        end
    end
end

% add dropout to reduce overfitting
drop = dropoutLayer(0.4,'Name','new_drop');
lgraph = addLayers(lgraph,drop);
lgraph = disconnectLayers(lgraph, 'pool5', 'fc1000');
lgraph = connectLayers(lgraph,'pool5','new_drop');
lgraph = connectLayers(lgraph,'new_drop', 'fc1000');

% replace fully connected layer, softmax and prediction
numClasses = 2;
newFc = fullyConnectedLayer(numClasses, 'Name', 'new_fc', 'WeightLearnRateFactor', 10, 'BiasLearnRateFactor', 10);
newFc.Weights = randn([numClasses 512]) * 0.01;  % Init weights (from avgPool size 2048, small std to avoid exploding grads)
newFc.Bias = zeros([numClasses 1]);  % Init bias as zeros
lgraph = replaceLayer(lgraph,'fc1000',newFc);
% lgraph = connectLayers(lgraph, lgraph.Layers(68,1).Name, 'new_fc');

newSoftmax = softmaxLayer('Name', 'new_softmax');
lgraph = replaceLayer(lgraph,'prob',newSoftmax);
% lgraph = connectLayers(lgraph, 'new_fc', 'new_softmax');

classes = categorical([0,1]);  % Or your labels
newClass = classificationLayer('Name', 'new_classoutput', 'Classes', classes);
lgraph = replaceLayer(lgraph,'ClassificationLayer_predictions',newClass);
% lgraph = connectLayers(lgraph, 'new_softmax', 'new_classoutput');

net = assembleNetwork(lgraph);

options = trainingOptions('adam', ...
    'InitialLearnRate', 1e-5, ...  % Low for transfer
    'MaxEpochs', 25, ...
    'ValidationPatience', 3,...
    'L2Regularization', 0.005,...
    'Shuffle', 'every-epoch', ...
    'GradientThresholdMethod', 'l2norm',...
    'GradientThreshold', 1,...
    'MiniBatchSize', 64, ...  % Adjust based on GPU memory
    'ValidationFrequency', 30, ...
    'Verbose', false, ...
    'ExecutionEnvironment','GPU');

% K-fold parameter
K = 10;


%% CV evaluation (leaky) - sanity check only!!
% stratified partition of full set, regardless of subject grouping
cv_leaky = cvpartition(y, 'KFold', K, 'Stratify', true);
acc_leaky = zeros(K,1);
Yp_leaky = cell(K,1);
Yt_leaky = cell(K,1);
Sc_leaky = cell(K,1);

for fold = 1:10
    tic
    trainIdx = training(cv_leaky, fold);
    valIdx = test(cv_leaky, fold);
    imdsTrain = subset(XT_full, trainIdx);
    imdsVal = subset(XT_full, valIdx);
    netFold = trainNetwork(imdsTrain, lgraph, options);  % Use lgraph or updated net
    [YPred, scores] = classify(netFold, imdsVal);
    valAll = readall(imdsVal);
    YTrue = categorical(cell2mat(valAll(:,2)));
    acc_leaky(fold) = mean(YPred == YTrue);
    Yp_leaky{fold} = YPred;
    Yt_leaky{fold} = YTrue;
    Sc_leaky{fold} = scores;
    close all
    disp(['Fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end
meanAcc = mean(acc_leaky);
stdAcc = std(acc_leaky);
fprintf('Leaky 10-Fold CV Mean Accuracy: %.2f%% (SD: %.2f%%)\n', meanAcc*100, stdAcc*100);

%% Exact CV evaluation - leaky vs subject-based
% stratified partition of subjects, then gain all examples
cv_subj = cvpartition(y_subj, 'KFold', K, 'Stratify', true);
acc_subj = zeros(K,1);
Yp_subj = cell(K,1);
Yt_subj = cell(K,1);
Sc_subj = cell(K,1);

for fold = 1:10
    tic
    trainSubj = training(cv_subj, fold);
    valSubj = test(cv_subj, fold);
    trainIdx = logical(ismember(subjIDs,X_subj(trainSubj)));
    valIdx = logical(ismember(subjIDs,X_subj(valSubj)));
    imdsTrain = subset(XT_full, trainIdx);
    imdsVal = subset(XT_full, valIdx);
    netFold = trainNetwork(imdsTrain, lgraph, options);  % Use lgraph or updated net
    [YPred, scores] = classify(netFold, imdsVal);
    valAll = readall(imdsVal);
    YTrue = categorical(cell2mat(valAll(:,2)));
    acc_subj(fold) = mean(YPred == YTrue);
    Yp_subj{fold} = YPred;
    Yt_subj{fold} = YTrue;
    Sc_subj{fold} = scores;
    close all
    disp(['Fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end
meanAcc = mean(acc_subj);
stdAcc = std(acc_subj);
fprintf('Subject-based 10-Fold CV Mean Accuracy: %.2f%% (SD: %.2f%%)\n', meanAcc*100, stdAcc*100);

%% save performance

perf_leaky = struct(...
    'Yt', Yt_leaky,...
    'Yp', Yp_leaky,...
    'Sc', Sc_leaky);

perf_subj = struct(...
    'Yt', Yt_subj,...
    'Yp', Yp_subj,...
    'Sc', Sc_subj);

save('ws_perf_CNN_dataset1.mat','perf_leaky','perf_subj')