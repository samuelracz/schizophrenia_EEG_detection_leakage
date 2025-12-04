% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

addpath(genpath('functions'))
fnames = dir('Dataset2_ML_features/*.mat');
ns = length(fnames);
rng(42, 'twister') % for reproducibility

%% load data
for subj = 1:ns
    load([fnames(subj).folder '/' fnames(subj).name])

    if subj == 1
        subjIDs = unique([table_feat.subjID]);
        y = unique([table_feat.y]);
        X = mean(removevars(table_feat,{'subjID','y'}),1);
    else
        subjIDs = cat(1,subjIDs,unique(table_feat.subjID));

        y = cat(1,y,unique(table_feat.y));

        X = cat(1,X,mean(removevars(table_feat,{'subjID','y'}),1));
    end
end


%% create extensive LOSO-CV
subjID_hc = subjIDs(y==0);
subjID_sz = subjIDs(y==1);
Nhc = sum(y==0);
Nsz = sum(y==1);
Ncv = Nhc*Nsz;

CVpartition = cell(Ncv,2);
idx = 0;
for g1 = 1:Nhc
    for g2 = 1:Nsz
        idx = idx + 1;
        CVpartition{idx,1} = subjID_hc{g1};
        CVpartition{idx,2} = subjID_sz{g2};
    end
end

feat_full = X.Properties.VariableNames';
Nf = length(feat_full);
K = 10; % k parameter for ReliefF
k_top = 10; % number of top features used in final model

%% CV evaluation - leaky LOSO-CV

% z-score table to normalize features
X_z = X;
for f = 1:Nf
    X_z.(feat_full{f}) = zscore(X.(feat_full{f}));
end

% Perform feature ranking using the ReliefF algorithm
[feat_rank_initial, feat_weight] = relieff(table2array(X_z),y,K);

% select only relevant features (exclude negative weights)
X_rel = X_z(:,feat_weight>0);

% select features using mRMR
feat_rank = fscmrmr(X_rel,y);
X_red = X_rel(:,feat_rank(1:k_top));

% save top features
feat_leaky = X_red.Properties.VariableNames;

% evaluate extensive LOSO-CV performance
Yp_leaky = cell(Ncv,1);
Yt_leaky = cell(Ncv,1);
Sc_leaky = cell(Ncv,1);
acc_leaky = zeros(Ncv,1);

tic
for cv = 1:Ncv
    % select train and validation data
    trainIdx = ~ismember(subjIDs,CVpartition(cv,:));
    valIdx = ismember(subjIDs,CVpartition(cv,:));

    X_train = X_red(trainIdx,:);
    y_train = y(trainIdx);

    X_val = X_red(valIdx,:);
    y_val = y(valIdx);

    % train model
    mdl = fitcsvm(X_train,y_train,'Standardize',false,'KernelFunction','RBF','KernelScale','auto');
    mdl.ScoreTransform = 'logit';

    % predict
    [y_pred, scores] = predict(mdl,X_val);

    % store outcome
    Yp_leaky{cv} = y_pred;
    Yt_leaky{cv} = y_val;
    Sc_leaky{cv} = scores;
    acc_leaky(cv) = mean(y_pred == y_val);

    % display progress
    if rem(cv,20) == 0
        disp(['Leaky CV: Fold ' num2str(cv) ' done in ' num2str(toc,'%.2f') 's.'])
        tic
    elseif cv == Ncv
        disp(['Leaky CV: Fold ' num2str(cv) ' done in ' num2str(toc,'%.2f') 's.'])
    end
end


%% CV evaluation - leake-free LOSO-CV

Yp_subj = cell(Ncv,1);
Yt_subj = cell(Ncv,1);
Sc_subj = cell(Ncv,1);
feat_subj = cell(Ncv,k_top);
acc_subj = zeros(Ncv,1);
tic
for cv = 1:Ncv
    % select train and validation data
    trainIdx = ~ismember(subjIDs,CVpartition(cv,:));
    valIdx = ismember(subjIDs,CVpartition(cv,:));

    % full training data
    X_train_full = X(trainIdx,:);
    y_train = y(trainIdx);

    % standardize data (save mean and std for validation set)
    f_mean = zeros(Nf,1);
    f_std = zeros(Nf,1);
    X_train_z = X_train_full;
    for f = 1:Nf
        feat = X_train_full.(feat_full{f});
        f_mean(f) = mean(feat);
        f_std(f) = std(feat);
        X_train_z.(feat_full{f}) = (feat - mean(feat))/std(feat);
    end

    % Perform feature ranking using the ReliefF algorithm
    [feat_rank_initial, feat_weight] = relieff(table2array(X_train_z),y_train,K);

    % select only relevant features (exclude negative weights)
    X_train_rel = X_train_z(:,feat_weight>0);

    % select features using mRMR - final train data
    feat_rank = fscmrmr(X_train_rel,y_train);
    X_train = X_train_rel(:,feat_rank(1:k_top));

    % save top features
    feat_topk = X_train.Properties.VariableNames';

    % full validation data
    X_val_full = X(valIdx,:);
    y_val = y(valIdx);

    % select relevant features
    X_val_rel = X_val_full(:,feat_weight>0);
    X_val_red = X_val_rel(:,feat_rank(1:k_top));

    % standardize variables in validation set
    X_val_z = X_val_red;
    for f = 1:k_top
        feat_idx = find(strcmp(feat_full,feat_topk{f}));
        feat = X_val_red.(feat_full{feat_idx});
        X_val_z.(feat_full{feat_idx}) = (feat - f_mean(feat_idx))/f_std(feat_idx);
    end

    % final validation data
    X_val = X_val_z;

    % train model
    mdl = fitcsvm(X_train,y_train,'Standardize',false,'KernelFunction','RBF','KernelScale','auto');
    mdl.ScoreTransform = 'logit';

    % predict
    [y_pred, scores] = predict(mdl,X_val);

    % store outcome
    Yp_subj{cv} = y_pred;
    Yt_subj{cv} = y_val;
    Sc_subj{cv} = scores;
    feat_subj(cv,:) = feat_topk';
    acc_subj(cv) = mean(y_pred == y_val);

    % display progress
    if rem(cv,20) == 0
        disp(['Leakage-free CV: Fold ' num2str(cv) ' done in ' num2str(toc,'%.2f') 's.'])
        tic
    elseif cv == Ncv
        disp(['Leakage-free CV: Fold ' num2str(cv) ' done in ' num2str(toc,'%.2f') 's.'])
    end
end

%% save performance

perf_leaky = struct(...
    'Yt', Yt_leaky,...
    'Yp', Yp_leaky,...
    'Sc', Sc_leaky);

perf_subj = struct(...
    'Yt', Yt_subj,...
    'Yp', Yp_subj,...
    'Sc', Sc_subj);

save('ws_perf_ML_dataset2.mat','perf_leaky','perf_subj')
