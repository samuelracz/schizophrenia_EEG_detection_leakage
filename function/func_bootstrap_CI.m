function [ CI ] = func_bootstrap_CI(acc_leaky, acc_subj, K)

bootfun = @(x) mean(x(:,1)) - mean(x(:,2));

CI = bootci(K, {bootfun, [acc_leaky, acc_subj]}, 'Type', 'percentile');


end