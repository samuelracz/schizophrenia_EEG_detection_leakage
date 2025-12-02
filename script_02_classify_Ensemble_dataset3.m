% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

fnames = dir('Dataset3_ML_features/*.mat');
ns = length(fnames);

%% load data
for subj = 1:ns
    load([fnames(subj).folder '/' fnames(subj).name])

    if subj == 1
        subjIDs = table_feat.subjID;
        X_subj = unique([table_feat.subjID]);
        y_subj = unique([table_feat.y]);

        X = removevars(table_feat,{'subjID','y'});
        y = table_feat.y;
    else
        subjIDs = cat(1,subjIDs,table_feat.subjID);
        X_subj = cat(1,X_subj,unique(table_feat.subjID));
        y_subj = cat(1,y_subj,unique(table_feat.y));

        X = cat(1,X,removevars(table_feat,{'subjID','y'}));
        y = cat(1,y,table_feat.y);
    end

end

%% Gradient Boosting Ensemble parameters

feat_full = X.Properties.VariableNames';
Nf = length(feat_full);
Kr = 10; % k parameter for ReliefF
K = 10; % k-fold cross-validation

t = templateTree('MaxNumSplits', 5, 'Surrogate', 'on');
method = 'GentleBoost';
numlearningcycles = 100;
learnrate = 0.1;

cv_leaky = cvpartition(y, 'KFold', K, 'Stratify', true);
cv_subj = cvpartition(y_subj, 'KFold', K, 'Stratify', true);

%% CV evaluation - Leaky train-test, leaky feature selection

acc_leaky_full = zeros(K,1);
Yp_leaky_full = cell(K,1);
Yt_leaky_full = cell(K,1);
Sc_leaky_full = cell(K,1);

% leaky feature selection
% z-score table to normalize features
X_z = table();
for f = 1:Nf
    X_z.(feat_full{f}) = zscore(X.(feat_full{f}));
end

% Perform feature ranking using the ReliefF algorithm
[feat_rank, feat_weight] = relieff(table2array(X_z),y,Kr);

% select only relevant features (exclude negative weights)
if sum(feat_weight>0) > Nf/2
    X_LTTFS = X_z(:,feat_rank(1:floor(Nf/2)));
else
    X_LTTFS = X_z(:,feat_weight>0);
end

for fold = 1:K
    tic
    % train and validation data
    trainIdx = training(cv_leaky, fold);
    valIdx = test(cv_leaky, fold);

    X_train = X_LTTFS(trainIdx,:);
    y_train = y(trainIdx);
    X_val = X_LTTFS(valIdx,:);
    y_true = y(valIdx);

    % train model
    mdlFold = fitensemble(X_train,y_train,method,numlearningcycles,t,'LearnRate',learnrate);
    mdlFold.ScoreTransform = 'logit';

    % predict
    [y_pred, scores] = predict(mdlFold, X_val);

    acc_leaky_full(fold) = mean(y_pred == y_true);
    Yp_leaky_full{fold} = y_pred;
    Yt_leaky_full{fold} = y_true;
    Sc_leaky_full{fold} = scores;

    disp(['Leaky TT-FS: fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end


%% ablation 1: leaky train-test split, leakage-free feature selection

acc_leaky_ablFS = zeros(K,1);
Yp_leaky_ablFS = cell(K,1);
Yt_leaky_ablFS = cell(K,1);
Sc_leaky_ablFS = cell(K,1);

for fold = 1:K
    tic
    % initial train data
    trainIdx = training(cv_leaky, fold);
    valIdx = test(cv_leaky, fold);

    X_train_full = X(trainIdx,:);
    y_train = y(trainIdx);

    % feature selection
    % standardize data (save mean and std for validation set)
    f_mean = zeros(Nf,1);
    f_std = zeros(Nf,1);
    X_train_z = table();
    for f = 1:Nf
        feat = X_train_full.(feat_full{f});
        f_mean(f) = mean(feat);
        f_std(f) = std(feat);
        X_train_z.(feat_full{f}) = (feat - mean(feat))/std(feat);
    end

    % Perform feature ranking using the ReliefF algorithm
    [feat_rank, feat_weight] = relieff(table2array(X_train_z),y_train,K);

    % select only relevant features (exclude negative weights)
    if sum(feat_weight>0) > Nf/2
        X_LTT = X_train_z(:,feat_rank(1:floor(Nf/2)));
    else
        X_LTT = X_train_z(:,feat_weight>0);
    end

    % final train data
    X_train = X_LTT;

    % save top features
    feat_topk = X_train.Properties.VariableNames';
    Nf_topk = length(feat_topk);

    % full validation data
    X_val_full = X(valIdx,:);
    y_true = y(valIdx);

    % select relevant features
    if sum(feat_weight>0) > Nf/2
        X_val_red = X_val_full(:,feat_rank(1:floor(Nf/2)));
    else
        X_val_red = X_val_full(:,feat_weight>0);
    end

    % standardize variables in validation set
    X_val_z = table();
    for f = 1:Nf_topk
        feat_idx = find(strcmp(feat_full,feat_topk{f}));
        feat = X_val_red.(feat_full{feat_idx});
        X_val_z.(feat_full{feat_idx}) = (feat - f_mean(feat_idx))/f_std(feat_idx);
    end

    % final validation data
    X_val = X_val_z;

    mdlFold = fitensemble(X_train,y_train,method,numlearningcycles,t,'LearnRate',learnrate);
    mdlFold.ScoreTransform = 'logit';

    [y_pred, scores] = predict(mdlFold, X_val);

    acc_leaky_ablFS(fold) = mean(y_pred == y_true);
    Yp_leaky_ablFS{fold} = y_pred;
    Yt_leaky_ablFS{fold} = y_true;
    Sc_leaky_ablFS{fold} = scores;

    disp(['Leaky TT: fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end

%% ablation 2: subject-based train-test split, leaky feature selection

acc_leaky_ablTT = zeros(K,1);
Yp_leaky_ablTT = cell(K,1);
Yt_leaky_ablTT = cell(K,1);
Sc_leaky_ablTT = cell(K,1);

% leaky feature selection
% z-score table to normalize features
X_z = table();
for f = 1:Nf
    X_z.(feat_full{f}) = zscore(X.(feat_full{f}));
end

% Perform feature ranking using the ReliefF algorithm
[feat_rank, feat_weight] = relieff(table2array(X_z),y,Kr);

% select only relevant features (exclude negative weights)
if sum(feat_weight>0) > Nf/2
    X_LTTFS = X_z(:,feat_rank(1:floor(Nf/2)));
else
    X_LTTFS = X_z(:,feat_weight>0);
end

for fold = 1:K
    tic
    % train and validation data - subject-based
    trainSubj = training(cv_subj, fold);
    valSubj = test(cv_subj, fold);
    trainIdx = logical(ismember(subjIDs,X_subj(trainSubj)));
    valIdx = logical(ismember(subjIDs,X_subj(valSubj)));

    X_train = X_LTTFS(trainIdx,:);
    y_train = y(trainIdx);
    X_val = X_LTTFS(valIdx,:);
    y_true = y(valIdx);

    % train model
    mdlFold = fitensemble(X_train,y_train,method,numlearningcycles,t,'LearnRate',learnrate);
    mdlFold.ScoreTransform = 'logit';

    % predict
    [y_pred, scores] = predict(mdlFold, X_val);

    acc_leaky_ablTT(fold) = mean(y_pred == y_true);
    Yp_leaky_ablTT{fold} = y_pred;
    Yt_leaky_ablTT{fold} = y_true;
    Sc_leaky_ablTT{fold} = scores;

    disp(['Leaky FS: fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end

%% CV evaluation - Leakage-free
acc_subj = zeros(K,1);
Yp_subj = cell(K,1);
Yt_subj = cell(K,1);
Sc_subj = cell(K,1);
feat_subj = cell(K,1);

for fold = 1:K
    tic
    % train and validation data - subject-based
    trainSubj = training(cv_subj, fold);
    valSubj = test(cv_subj, fold);
    trainIdx = logical(ismember(subjIDs,X_subj(trainSubj)));
    valIdx = logical(ismember(subjIDs,X_subj(valSubj)));

    % initial train data
    X_train_full = X(trainIdx,:);
    y_train = y(trainIdx);

    % feature selection
    % standardize data (save mean and std for validation set)
    f_mean = zeros(Nf,1);
    f_std = zeros(Nf,1);
    X_train_z = table();
    for f = 1:Nf
        feat = X_train_full.(feat_full{f});
        f_mean(f) = mean(feat);
        f_std(f) = std(feat);
        X_train_z.(feat_full{f}) = (feat - mean(feat))/std(feat);
    end

    % Perform feature ranking using the ReliefF algorithm
    [feat_rank, feat_weight] = relieff(table2array(X_train_z),y_train,K);

    % select only relevant features (exclude negative weights)
    if sum(feat_weight>0) > Nf/2
        X_LTT = X_train_z(:,feat_rank(1:floor(Nf/2)));
    else
        X_LTT = X_train_z(:,feat_weight>0);
    end

    % final train data
    X_train = X_LTT;

    % save top features
    feat_topk = X_train.Properties.VariableNames';
    Nf_topk = length(feat_topk);

    % full validation data
    X_val_full = X(valIdx,:);
    y_true = y(valIdx);

    % select relevant features
    if sum(feat_weight>0) > Nf/2
        X_val_red = X_val_full(:,feat_rank(1:floor(Nf/2)));
    else
        X_val_red = X_val_full(:,feat_weight>0);
    end

    % standardize variables in validation set
    X_val_z = table();
    for f = 1:Nf_topk
        feat_idx = find(strcmp(feat_full,feat_topk{f}));
        feat = X_val_red.(feat_full{feat_idx});
        X_val_z.(feat_full{feat_idx}) = (feat - f_mean(feat_idx))/f_std(feat_idx);
    end

    % final validation data
    X_val = X_val_z;

    mdlFold = fitensemble(X_train,y_train,method,numlearningcycles,t,'LearnRate',learnrate);
    mdlFold.ScoreTransform = 'logit';

    [y_pred, scores] = predict(mdlFold, X_val);

    acc_subj(fold) = mean(y_pred == y_true);
    Yp_subj{fold} = y_pred;
    Yt_subj{fold} = y_true;
    Sc_subj{fold} = scores;
    feat_subj{fold} = {feat_topk};

    disp(['Leakage-free CV: fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end

%% saving results

perf_leaky_full = struct(...
    'Yt', Yt_leaky_full,...
    'Yp', Yp_leaky_full,...
    'Sc', Sc_leaky_full);

perf_leaky_ablFS = struct(...
    'Yt', Yt_leaky_ablFS,...
    'Yp', Yp_leaky_ablFS,...
    'Sc', Sc_leaky_ablFS);

perf_leaky_ablTT = struct(...
    'Yt', Yt_leaky_ablTT,...
    'Yp', Yp_leaky_ablTT,...
    'Sc', Sc_leaky_ablTT);

perf_subj = struct(...
    'Yt', Yt_subj,...
    'Yp', Yp_subj,...
    'Sc', Sc_subj);

save('ws_perf_DL_dataset3.mat','perf_leaky_full','perf_leaky_ablFS','perf_leaky_ablTT','perf_subj')