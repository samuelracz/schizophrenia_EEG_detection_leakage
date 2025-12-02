function [ p ] = func_permutation_difference(acc_leaky, acc_subj, K)

N = numel(acc_leaky);

% observed difference
d_true = mean(acc_leaky) - mean(acc_subj);

% pooled differences
d_pool = [acc_leaky; acc_subj];

% random permutations
d_perm = zeros(K,1);
for k = 1:K
    perm_idx = randperm(2*N);
    acc1 = d_pool(perm_idx(1:N));
    acc2 = d_pool(perm_idx(N+1:end));
    d_perm(k) = mean(acc1) - mean(acc2);
end

p = (sum(d_perm >= d_true) + 1)/(K+1);

end