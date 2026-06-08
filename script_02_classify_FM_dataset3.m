
fnames = dir('Dataset3_REVE_features/*.mat');
ns = length(fnames);
rng(42) % for reproducibility

%% load data
for subj = 1:ns
    load([fnames(subj).folder '/' fnames(subj).name])

    if subj == 1
        subjIDs = repmat({subject_id},n_epochs,1);
        X_subj = {subject_id};
        y_subj = label;

        X = double(embeddings);
        y = repmat(label,n_epochs,1);
    else
        subjIDs = cat(1,subjIDs,repmat({subject_id},n_epochs,1));
        X_subj = cat(1,X_subj,{subject_id});
        y_subj = cat(1,y_subj,label);

        X = cat(1,X,double(embeddings));
        y = cat(1,y,repmat(label,n_epochs,1));
    end

end

%% Cross-validation parameters

K = 10; % k-fold cross-validation

% stratified partition of full set, regardless of subject grouping
cv_leaky = cvpartition(y, 'KFold', K, 'Stratify', true);

% stratified partition of subjects, then gain all examples
cv_subj = cvpartition(y_subj, 'KFold', K, 'Stratify', true);

%% CV evaluation - Leaky train-test
acc_leaky = zeros(K,1);
Yp_leaky = cell(K,1);
Yt_leaky = cell(K,1);
Sc_leaky = cell(K,1);

for fold = 1:K
    tic
    % train and validation data
    trainIdx = training(cv_leaky, fold);
    valIdx = test(cv_leaky, fold);

    X_train = X(trainIdx,:);
    y_train = y(trainIdx);
    X_val = X(valIdx,:);
    y_true = y(valIdx);

    % train model
    mdlFold = fitcsvm(X_train, y_train, 'KernelFunction','linear', 'Standardize',true);

    % predict
    [y_pred, scores] = predict(mdlFold, X_val);

    acc_leaky(fold) = mean(y_pred == y_true);
    Yp_leaky{fold} = y_pred;
    Yt_leaky{fold} = y_true;
    Sc_leaky{fold} = scores;

    disp(['Leaky DP CV: fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
end

meanAcc = mean(acc_leaky);
stdAcc = std(acc_leaky);
fprintf('Leaky 10-Fold CV Mean Accuracy: %.2f%% (SD: %.2f%%)\n', meanAcc*100, stdAcc*100);

%% CV evaluation - Leakage-free
acc_subj = zeros(K,1);
Yp_subj = cell(K,1);
Yt_subj = cell(K,1);
Sc_subj = cell(K,1);

for fold = 1:K
    tic
    % train and validation data - subject-based
    trainSubj = training(cv_subj, fold);
    valSubj = test(cv_subj, fold);
    trainIdx = logical(ismember(subjIDs,X_subj(trainSubj)));
    valIdx = logical(ismember(subjIDs,X_subj(valSubj)));

    % initial train data
    X_train = X(trainIdx,:);
    y_train = y(trainIdx);

    % full validation data
    X_val = X(valIdx,:);
    y_true = y(valIdx);

    mdlFold = fitcsvm(X_train, y_train, 'KernelFunction','linear', 'Standardize',true);

    [y_pred, scores] = predict(mdlFold, X_val);

    acc_subj(fold) = mean(y_pred == y_true);
    Yp_subj{fold} = y_pred;
    Yt_subj{fold} = y_true;
    Sc_subj{fold} = scores;

    disp(['Leakage-free CV: fold ' num2str(fold) ' done in ' num2str(toc,'%.2f') 's.'])
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

save('ws_perf_FM_dataset3.mat','perf_leaky','perf_subj')