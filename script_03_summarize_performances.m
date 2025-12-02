% Author: F. Samuel Racz, 2025, The University of Texas at Austin
% email: fsr324@austin.utexas.edu

addpath(genpath('functions'))

table_tmp1 = table();
table_tmp1.dataset = cell(3,1);
table_tmp1.ACC_Lm = zeros(3,1);
table_tmp1.ACC_Ls = zeros(3,1);
table_tmp1.ACC_Sm = zeros(3,1);
table_tmp1.ACC_Ss = zeros(3,1);
table_tmp1.dACC = zeros(3,1);
table_tmp1.CI_lo = zeros(3,1);
table_tmp1.CI_hi = zeros(3,1);
table_tmp1.p_L = zeros(3,1);
table_tmp1.ch = zeros(3,1);
table_tmp1.p_ch = zeros(3,1);


table_tmp2 = table();
table_tmp2.dataset = cell(3,1);
table_tmp2.ACC_LFm = zeros(3,1);
table_tmp2.ACC_LFs = zeros(3,1);
table_tmp2.ACC_LPm = zeros(3,1);
table_tmp2.ACC_LPs = zeros(3,1);
table_tmp2.ACC_LSm = zeros(3,1);
table_tmp2.ACC_LSs = zeros(3,1);
table_tmp2.ACC_Sm = zeros(3,1);
table_tmp2.ACC_Ss = zeros(3,1);
table_tmp2.dACC = cell(3,1);
table_tmp2.CI = cell(3,1);
table_tmp2.p_L = cell(3,1);
table_tmp2.ch = cell(3,1);
table_tmp2.p_ch = cell(3,1);

disp('==||==||==')

%% summarize CNN performance
load('ws_perf_CNN_dataset1.mat')
CNN_LD1 = perf_leaky;
CNN_LFD1 = perf_subj;

load('ws_perf_CNN_dataset2.mat')
CNN_LD2 = perf_leaky;
CNN_LFD2 = perf_subj;

load('ws_perf_CNN_dataset3.mat')
CNN_LD3 = perf_leaky;
CNN_LFD3 = perf_subj;

table_CNN = table_tmp1;

for i = 1:3
    switch i
        case 1
            dl = CNN_LD1;
            ds = CNN_LFD1;
            table_CNN.dataset{i} = 'MSU';
        case 2 
            dl = CNN_LD2;
            ds = CNN_LFD2;
            table_CNN.dataset{i} = 'RepOD';
        case 3
            dl = CNN_LD3;
            ds = CNN_LFD3;
            table_CNN.dataset{i} = 'SU-SZ';
    end

    K = length(dl);
    acc_leaky = zeros(K,1);
    acc_subj = zeros(K,1);
    for k = 1:K
        acc_leaky(k) = mean(dl(k).Yp==dl(k).Yt);
        acc_subj(k) = mean(ds(k).Yp==ds(k).Yt);
    end
    CI = func_bootstrap_CI(acc_leaky,acc_subj,10000);
    pL = func_permutation_difference(acc_leaky,acc_subj,10000);
    

    table_CNN.ACC_Lm(i) = mean(acc_leaky);
    table_CNN.ACC_Ls(i) = std(acc_leaky);
    table_CNN.ACC_Sm(i) = mean(acc_subj);
    table_CNN.ACC_Ss(i) = std(acc_subj);
    table_CNN.dACC(i) = mean(acc_leaky) - mean(acc_subj);
    table_CNN.CI_lo(i) = min(CI);
    table_CNN.CI_hi(i) = max(CI);
    table_CNN.p_L(i) = pL;

    disp([table_CNN.dataset{i} ' CNN approach: Leaky ' ...
        num2str(table_CNN.ACC_Lm(i)*100,'%.2f') '% (SD: ' num2str(table_CNN.ACC_Ls(i)*100,'%.2f') ...
        '%) vs. Leakage-free ' num2str(table_CNN.ACC_Sm(i)*100,'%.2f') '% (SD:' num2str(table_CNN.ACC_Ss(i)*100,'%.2f') ...
        '%), diff: '  num2str(table_CNN.dACC(i)*100,'%.2f') '% [' ...
        num2str(table_CNN.CI_lo(i)*100,'%.2f') '; ' num2str(table_CNN.CI_hi(i)*100,'%.2f') ...
        '], p=' num2str(table_CNN.p_L(i),'%.4f')])

end

disp('==||==||==')

%% summarize Ensemble performance
load('ws_perf_Ensemble_dataset1.mat')
DL_LD1 = perf_leaky_ablFS;
DL_LFD1 = perf_subj;

load('ws_perf_Ensemble_dataset2.mat')
DL_LD2 = perf_leaky_ablFS;
DL_LFD2 = perf_subj;

load('ws_perf_Ensemble_dataset3.mat')
DL_LD3 = perf_leaky_ablFS;
DL_LFD3 = perf_subj;

table_DL = table_tmp1;

for i = 1:3
    switch i
        case 1
            dl = DL_LD1;
            ds = DL_LFD1;
            table_DL.dataset{i} = 'MSU';
        case 2 
            dl = DL_LD2;
            ds = DL_LFD2;
            table_DL.dataset{i} = 'RepOD';
        case 3
            dl = DL_LD3;
            ds = DL_LFD3;
            table_DL.dataset{i} = 'SU-SZ';
    end

    K = length(dl);
    acc_leaky = zeros(K,1);
    acc_subj = zeros(K,1);
    for k = 1:K
        acc_leaky(k) = mean(dl(k).Yp==dl(k).Yt);
        acc_subj(k) = mean(ds(k).Yp==ds(k).Yt);
    end
    CI = func_bootstrap_CI(acc_leaky,acc_subj,10000);
    pL = func_permutation_difference(acc_leaky,acc_subj,10000);
    

    table_DL.ACC_Lm(i) = mean(acc_leaky);
    table_DL.ACC_Ls(i) = std(acc_leaky);
    table_DL.ACC_Sm(i) = mean(acc_subj);
    table_DL.ACC_Ss(i) = std(acc_subj);
    table_DL.dACC(i) = mean(acc_leaky) - mean(acc_subj);
    table_DL.CI_lo(i) = min(CI);
    table_DL.CI_hi(i) = max(CI);
    table_DL.p_L(i) = pL;

    disp([table_DL.dataset{i} ' DL approach: Leaky ' ...
        num2str(table_DL.ACC_Lm(i)*100,'%.2f') '% (SD: ' num2str(table_DL.ACC_Ls(i)*100,'%.2f') ...
        '%) vs. Leakage-free ' num2str(table_DL.ACC_Sm(i)*100,'%.2f') '% (SD:' num2str(table_DL.ACC_Ss(i)*100,'%.2f') ...
        '%), diff: '  num2str(table_DL.dACC(i)*100,'%.2f') '% [' ...
        num2str(table_DL.CI_lo(i)*100,'%.2f') '; ' num2str(table_DL.CI_hi(i)*100,'%.2f') ...
        '], p=' num2str(table_DL.p_L(i),'%.4f')])

end

disp('==||==||==')

%% summarize ML performance
load('ws_perf_ML_dataset1.mat')
ML_LD1 = perf_leaky;
ML_LFD1 = perf_subj;

load('ws_perf_ML_dataset2.mat')
ML_LD2 = perf_leaky;
ML_LFD2 = perf_subj;

load('ws_perf_ML_dataset3.mat')
ML_LD3 = perf_leaky;
ML_LFD3 = perf_subj;

table_ML = table_tmp1;

for i = 1:3
    switch i
        case 1
            dl = ML_LD1;
            ds = ML_LFD1;
            table_ML.dataset{i} = 'MSU';
        case 2 
            dl = ML_LD2;
            ds = ML_LFD2;
            table_ML.dataset{i} = 'RepOD';
        case 3
            dl = ML_LD3;
            ds = ML_LFD3;
            table_ML.dataset{i} = 'SU-SZ';
    end

    K = length(dl);
    acc_leaky = zeros(K,1);
    acc_subj = zeros(K,1);
    for k = 1:K
        acc_leaky(k) = mean(dl(k).Yp==dl(k).Yt);
        acc_subj(k) = mean(ds(k).Yp==ds(k).Yt);
    end
    CI = func_bootstrap_CI(acc_leaky,acc_subj,10000);
    pL = func_permutation_difference(acc_leaky,acc_subj,10000);
    

    table_ML.ACC_Lm(i) = mean(acc_leaky);
    table_ML.ACC_Ls(i) = std(acc_leaky);
    table_ML.ACC_Sm(i) = mean(acc_subj);
    table_ML.ACC_Ss(i) = std(acc_subj);
    table_ML.dACC(i) = mean(acc_leaky) - mean(acc_subj);
    table_ML.CI_lo(i) = min(CI);
    table_ML.CI_hi(i) = max(CI);
    table_ML.p_L(i) = pL;

    disp([table_ML.dataset{i} ' ML approach: Leaky ' ...
        num2str(table_ML.ACC_Lm(i)*100,'%.2f') '% (SD: ' num2str(table_ML.ACC_Ls(i)*100,'%.2f') ...
        '%) vs. Leakage-free ' num2str(table_ML.ACC_Sm(i)*100,'%.2f') '% (SD:' num2str(table_ML.ACC_Ss(i)*100,'%.2f') ...
        '%), diff: '  num2str(table_ML.dACC(i)*100,'%.2f') '% [' ...
        num2str(table_ML.CI_lo(i)*100,'%.2f') '; ' num2str(table_ML.CI_hi(i)*100,'%.2f') ...
        '], p=' num2str(table_ML.p_L(i),'%.4f')])

end

disp('==||==||==')

    